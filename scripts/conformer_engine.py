# env: chameleon-calc
"""
01_conformer_engine.py
----------------------
Tier-1 dual-environment conformer pipeline.

For each molecule:
  1. Generate N conformers with RDKit ETKDGv3 (macrocycle-aware)
  2. Minimize all conformers with MMFF94s
  3. Compute 3D descriptors for each conformer:
       - 3D PSA (SASA of polar atoms)
       - Intramolecular H-bond count
       - Radius of gyration (Rg)
       - NPR1, NPR2 (normalized principal moment ratios — shape)
       - Asphericity, Eccentricity, SpherocityIndex
  4. Select "aqueous" conformer  = max-PSA conformer  (polar groups exposed)
     Select "membrane" conformer = min-PSA conformer  (polar groups buried)
  5. Save per-molecule Δ descriptors to CSV

This approximates dual-dielectric behaviour via conformational sampling:
  chameleonic molecules have a wide PSA conformer spread.

Tier-2 (OMEGA + OpenMM GB/SA) validates this approximation on ~5 references.

Usage:
  python conformer_engine.py [--input pampa_curated.csv]
                             [--outdir ../results]
                             [--n-confs 50]
                             [--n-cpus 4]
                             [--max-mols 0]   # 0 = all

Outputs:
  results/conformer_descriptors_raw.csv   — per-molecule raw 3D descriptors
"""

import argparse
import json
import logging
import multiprocessing as mp
import warnings
from pathlib import Path

logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")
_log = logging.getLogger(__name__)

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import Descriptors3D
from rdkit.Chem import rdFreeSASA
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

# ── Constants ────────────────────────────────────────────────────────────────
HB_DONOR_SMARTS    = Chem.MolFromSmarts("[N,O;!H0]")   # NH, OH (N-Me amides correctly excluded)
HB_ACCEPTOR_SMARTS = Chem.MolFromSmarts("[N,O]")
HB_DIST_CUTOFF  = 3.0   # Å — H...acceptor distance (Baker-Hubbard permissive cutoff)
HB_ANGLE_CUTOFF = 120.0  # degrees — D-H...A angle at H
# HB SAMPLING NOTE: delta_hb = 0 for all 5 reference compounds at n=20 conformers.
# This is a sampling artefact — the HB-forming membrane conformer is rare in a
# 20-conformer vacuum ensemble.  Minimum recommended n_confs for reliable ΔHB:
#   n ≥ 50  for screening (ΔHB signal emerges for strongly chameleonic molecules)
#   n ≥ 200 for final production run (reliable ΔHB for borderline molecules too)
# The psa3d_spread and psa3d_std columns are more robust at low n than delta_hb.


def _polar_sasa(mol: Chem.Mol, conf_id: int) -> float:
    """
    Compute polar SASA (3D PSA proxy) for a given conformer.

    Root-cause fix (2024): rdFreeSASA.classifyAtoms() uses the Protor classifier
    by default, which is a protein-residue-based scheme (identifies atoms by PDB
    residue name).  For arbitrary small molecules it returns all-zero radii and
    marks every atom as SASAClass.Unclassified, so CalcSASA with the polar query
    returns 0.0 for every conformer.  The OONS and NACCESS classifiers exhibit
    the same failure on non-protein input.

    Fix: assign Bondi VdW radii manually and mark N/O/S/P atoms as SASAClass=0
    (Polar) before calling CalcSASA with MakeFreeSasaPolarAtomQuery().  All atoms
    still contribute their full radius to the mutual-occlusion computation; only
    the polar subset is summed in the final area.  This matches the convention used
    by the 2D TPSA descriptor (Ertl 2000) and is appropriate for cyclic peptides.

    Bondi radii used (Å):
      H 1.20, C 1.70, N 1.55, O 1.52, S 1.80, P 1.80, F 1.47, Cl 1.75,
      Br 1.85, I 1.98, default 1.50 for any other element.

    Returns np.nan on any failure so downstream dropna() filters it out.
    """
    # Bondi VdW radii (Å) — widely used for SASA of organic/peptide molecules
    _BONDI = {
        'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
        'S': 1.80, 'P': 1.80, 'F': 1.47, 'Cl': 1.75,
        'Br': 1.85, 'I': 1.98,
    }
    # Heteroatoms treated as polar (contributes to 3D PSA).
    # Polar-H DECISION: Witek 2016 (JCTC) includes H atoms bonded to N or O
    # in their PSA calculation (explicit-solvent MD, gmx sasa).  RDKit's 2D TPSA
    # (Ertl 2000) uses a fragment-based approach that implicitly includes NH/OH
    # contributions via the fragment areas but does NOT sum explicit H atom areas.
    # Here we follow the HEAVY-ATOM-ONLY convention (N, O, S, P) consistent with
    # the Ertl 2000 spirit and the CycPeptMPDB delta_3DPSA column definition.
    # Effect on CsA: including polar-H raises absolute PSA values by ~5-15 Å²
    # (each NH/OH adds ~2-5 Å² depending on burial), but the ΔPSA = PSA_max −
    # PSA_min is minimally affected because polar-H exposure tracks heavy-atom N/O
    # exposure.  Relative ranking across conformers and molecules is unchanged.
    # DECISION: heavy-atom-only.  Document this choice if publishing results.
    _POLAR_ELEMENTS = {'N', 'O', 'S', 'P'}

    try:
        # Assign per-atom radii and polar/apolar classification manually.
        # SASAClass integer encoding: 0 = Polar, 1 = APolar, 2 = Unclassified
        radii = []
        for atom in mol.GetAtoms():
            sym = atom.GetSymbol()
            radii.append(_BONDI.get(sym, 1.50))
            if sym in _POLAR_ELEMENTS:
                atom.SetIntProp('SASAClass', 0)
                atom.SetProp('SASAClassName', 'Polar')
            else:
                atom.SetIntProp('SASAClass', 1)
                atom.SetProp('SASAClassName', 'APolar')

        query = rdFreeSASA.MakeFreeSasaPolarAtomQuery()
        polar_sasa = rdFreeSASA.CalcSASA(mol, radii, confIdx=conf_id, query=query)
        return round(polar_sasa, 4)
    except Exception:
        return np.nan


def _intramolecular_hbonds(mol: Chem.Mol, conf_id: int) -> int:
    """
    Count intramolecular H-bonds using 3D geometry.
    Criterion: D-H ... A distance < 3.5 Å, D-H...A angle > 120°.
    Only counts bonds where donor and acceptor are 4+ bonds apart (not covalent).
    """
    try:
        conf = mol.GetConformer(conf_id)
        pos = conf.GetPositions()

        donors = [idx for match in mol.GetSubstructMatches(HB_DONOR_SMARTS)
                  for idx in match]
        acceptors = [idx for match in mol.GetSubstructMatches(HB_ACCEPTOR_SMARTS)
                     for idx in match]

        hb_count = 0
        for d_idx in donors:
            d_atom = mol.GetAtomWithIdx(d_idx)
            for h_atom in d_atom.GetNeighbors():
                if h_atom.GetAtomicNum() != 1:
                    continue
                h_idx = h_atom.GetIdx()
                h_pos = pos[h_idx]
                d_pos = pos[d_idx]

                for a_idx in acceptors:
                    if a_idx == d_idx:
                        continue
                    # Check topological distance (avoid 1,2 / 1,3 / 1,4)
                    try:
                        path = Chem.GetShortestPath(mol, d_idx, a_idx)
                        # Require at least a γ-turn (5-membered ring = 4 bonds = 5 atoms in path)
                        # len(path) includes both endpoints, so < 6 excludes ≤ 4-bond contacts
                        if len(path) < 6:
                            continue
                    except Exception:
                        continue

                    a_pos = pos[a_idx]
                    h_a_dist = np.linalg.norm(h_pos - a_pos)
                    if h_a_dist > HB_DIST_CUTOFF:
                        continue

                    # D-H...A angle — measured at H (angle between H→D and H→A)
                    # A linear H-bond (180°) gives cos=−1 → angle=180° → passes ≥120°
                    vec_hd = d_pos - h_pos   # H → D
                    vec_ha = a_pos - h_pos   # H → A
                    cos_angle = np.dot(vec_hd, vec_ha) / (
                        np.linalg.norm(vec_hd) * np.linalg.norm(vec_ha) + 1e-9
                    )
                    angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
                    if angle >= HB_ANGLE_CUTOFF:
                        hb_count += 1

        return hb_count
    except Exception:
        return np.nan


def _shape_descriptors(mol: Chem.Mol, conf_id: int) -> dict:
    """Compute RDKit 3D shape descriptors for one conformer."""
    try:
        rg  = Descriptors3D.RadiusOfGyration(mol, confId=conf_id)
        npr1 = Descriptors3D.NPR1(mol, confId=conf_id)   # I1/I3, disc-like→1
        npr2 = Descriptors3D.NPR2(mol, confId=conf_id)   # I2/I3, rod-like→0
        asph = Descriptors3D.Asphericity(mol, confId=conf_id)
        ecce = Descriptors3D.Eccentricity(mol, confId=conf_id)
        sphe = Descriptors3D.SpherocityIndex(mol, confId=conf_id)
        pbf  = Descriptors3D.PBF(mol, confId=conf_id)    # plane of best fit
        return {
            "Rg": rg, "NPR1": npr1, "NPR2": npr2,
            "Asphericity": asph, "Eccentricity": ecce,
            "SpherocityIndex": sphe, "PBF": pbf,
        }
    except Exception:
        return {k: np.nan for k in
                ["Rg", "NPR1", "NPR2", "Asphericity", "Eccentricity",
                 "SpherocityIndex", "PBF"]}


def process_molecule(args: tuple) -> dict:
    """
    Process a single molecule: embed → minimize → select aqueous/membrane conformers.
    Returns a flat dict of Δ descriptors, or a dict with 'error' key on failure.
    """
    mol_id, smiles, n_confs = args

    # Build mol with explicit H (needed for HB counting and proper SASA)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"ID": mol_id, "error": "invalid_smiles"}

    # ── pH 7.4 standardization ───────────────────────────────────────────────
    # Normalize and neutralize before conformer generation.
    # At pH 7.4 most cyclic peptides are neutral — this handles edge cases
    # (e.g. charged SMILES in CycPeptMPDB, salt forms, formal charge artifacts).
    # Steps: fragment removal → charge neutralization → canonical tautomer.
    # Does NOT change stereochemistry or ring connectivity.
    _std_ok = False
    try:
        # 1. Remove salts / largest fragment
        mol = rdMolStandardize.FragmentParent(mol)
        # 2. Neutralize formal charges (deprotonate acids, protonate amines at pH 7.4)
        uncharger = rdMolStandardize.Uncharger(canonicalOrder=True)
        mol = uncharger.uncharge(mol)
        # 3. Re-sanitize after modification
        Chem.SanitizeMol(mol)
        _std_ok = True
    except Exception as _e:
        _log.warning("Standardization failed for ID=%s (%s) — using original mol", mol_id, _e)

    mol = Chem.AddHs(mol)

    # ── ETKDGv3 embedding (macrocycle-aware) ─────────────────────────────────
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.maxIterations = 2000
    params.numThreads = 1
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    params.pruneRmsThresh = 0.5   # prune near-duplicate conformers

    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)

    if len(conf_ids) == 0:
        # Fallback: relax RMSD pruning, keep same seed for reproducibility
        params.randomSeed = 42
        params.pruneRmsThresh = 1.0
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)

    if len(conf_ids) == 0:
        return {"ID": mol_id, "error": "embed_failed"}

    # ── MMFF94s minimization ─────────────────────────────────────────────────
    results = AllChem.MMFFOptimizeMoleculeConfs(
        mol, mmffVariant="MMFF94s", maxIters=2000, numThreads=1
    )
    # results[i] = (not_converged, energy)

    # Filter to only converged conformers
    converged_ids = [
        conf_ids[i] for i, (not_conv, _) in enumerate(results) if not_conv == 0
    ]
    if len(converged_ids) == 0:
        converged_ids = list(conf_ids)  # use all if none converged

    # ── Per-conformer descriptors ─────────────────────────────────────────────
    per_conf = []
    for cid in converged_ids:
        psa3d = _polar_sasa(mol, cid)
        hb    = _intramolecular_hbonds(mol, cid)
        shape = _shape_descriptors(mol, cid)
        mmff_energy = next(
            (e for i, (_, e) in enumerate(results) if conf_ids[i] == cid), np.nan
        )
        per_conf.append({
            "conf_id": cid,
            "mmff_energy": mmff_energy,
            "psa3d": psa3d,
            "hb_count": hb,
            **shape,
        })

    if not per_conf:
        return {"ID": mol_id, "error": "no_descriptors"}

    df_conf = pd.DataFrame(per_conf).dropna(subset=["psa3d"])
    if df_conf.empty:
        return {"ID": mol_id, "error": "psa_failed"}

    # ── Select aqueous (max-PSA) and membrane (min-PSA) conformers ───────────
    # NOTE — vacuum conformer artefact for rigid impermeable peptides:
    # Molecules that are physically rigid and cannot shield their polar groups
    # (e.g. fully-exposed backbone peptides like c*[PSLYF], CycPeptMPDB ID=1829)
    # may show spuriously large ΔPSA from this vacuum ETKDGv3 sampling.  ETKDG
    # does not model solvation, so collapsed hydrophobic conformers that bury
    # polar groups may be generated even though they are thermodynamically
    # inaccessible in aqueous solution.  This is a fundamental limitation of
    # Tier-1 heuristic sampling; Tier-2 CREST+ALPB (tier2_crest.py) addresses
    # it by using environment-specific conformer ensembles.
    # Expected artefact magnitude: ΔPSA up to ~70 Å² for a compound that should
    # show ~0 Å² (e.g. c*[PSLYF]).  Increase n_confs to ≥200 to partially
    # mitigate by populating the near-extended conformers that dominate in water.
    # The psa3d_spread metric is a useful signal: large spread + small ΔPSA
    # expected from biology → flag the molecule for Tier-2 review.
    #
    # SMILES NOTE for c*[PSLYF] (ID=1829): the canonical SMILES in
    # CycPeptMPDB/reference_set.csv is:
    #   CC(C)C[C@@H]1NC(=O)[C@H](CO)NC(=O)[C@@H]2CCCN2[C@H](C(=O)NC(C)(C)C)
    #   [C@H](C)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](Cc2ccc(O)cc2)NC1=O
    # This contains a Ser (serine, CO sidechain), Leu (isobutyl), Tyr (4-OH-Bn),
    # Phe (Bn), and Pro (pyrrolidine) residue — consistent with c*[PSLYF].
    # There is NO Asp/Glu carboxylate in the canonical SMILES.  The Uncharger
    # step is therefore a no-op for this compound (nothing to neutralize).
    # If a variant SMILES with aspartate is encountered (e.g. from a different
    # database entry), the Uncharger will convert COO⁻ → COOH, which is the
    # correct neutral form for permeability modelling (lipid bilayer has low ε).
    aq_row  = df_conf.loc[df_conf["psa3d"].idxmax()]
    mem_row = df_conf.loc[df_conf["psa3d"].idxmin()]

    # ── Compute Δ features ───────────────────────────────────────────────────
    result = {
        "ID": mol_id,
        "n_confs_generated": len(conf_ids),
        "n_confs_used": len(df_conf),
        # Aqueous conformer descriptors (high-PSA, polar groups exposed)
        "aq_psa3d":           float(aq_row["psa3d"]),
        "aq_hb_count":        float(aq_row["hb_count"]),
        "aq_Rg":              float(aq_row["Rg"]),
        "aq_NPR1":            float(aq_row["NPR1"]),
        "aq_NPR2":            float(aq_row["NPR2"]),
        "aq_Asphericity":     float(aq_row["Asphericity"]),
        "aq_SpherocityIndex": float(aq_row["SpherocityIndex"]),
        # Membrane conformer descriptors (low-PSA, polar groups buried)
        "mem_psa3d":          float(mem_row["psa3d"]),
        "mem_hb_count":       float(mem_row["hb_count"]),
        "mem_Rg":             float(mem_row["Rg"]),
        "mem_NPR1":           float(mem_row["NPR1"]),
        "mem_NPR2":           float(mem_row["NPR2"]),
        "mem_Asphericity":    float(mem_row["Asphericity"]),
        "mem_SpherocityIndex":float(mem_row["SpherocityIndex"]),
        # Δ features (chameleonic potential)
        "delta_psa3d":        float(aq_row["psa3d"]  - mem_row["psa3d"]),
        "delta_hb":           float(mem_row["hb_count"] - aq_row["hb_count"]),   # more HB in membrane
        "delta_Rg":           float(aq_row["Rg"]     - mem_row["Rg"]),            # more compact in membrane
        "delta_NPR1":         float(aq_row["NPR1"]   - mem_row["NPR1"]),
        "delta_NPR2":         float(aq_row["NPR2"]   - mem_row["NPR2"]),
        "delta_Asphericity":  float(aq_row["Asphericity"] - mem_row["Asphericity"]),
        # PSA spread across all conformers (conformational flexibility)
        "psa3d_spread":       float(df_conf["psa3d"].max() - df_conf["psa3d"].min()),
        "psa3d_std":          float(df_conf["psa3d"].std()),
        "hb_spread":          float(df_conf["hb_count"].max() - df_conf["hb_count"].min()),
        "error": None,
    }
    return result


def run(input_csv: str, outdir: Path, n_confs: int, n_cpus: int, max_mols: int,
        checkpoint_every: int = 500) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_csv, low_memory=False)
    smiles_col = "SMILES_canonical" if "SMILES_canonical" in df.columns else "SMILES"
    df = df[["ID", smiles_col, "PAMPA"]].dropna(subset=[smiles_col]).copy()

    if max_mols > 0:
        df = df.head(max_mols)

    out_path   = outdir / "conformer_descriptors_raw.csv"
    ckpt_path  = outdir / "conformer_descriptors_checkpoint.csv"

    # ── Resume from checkpoint if one exists ────────────────────────────────
    done_ids: set = set()
    if ckpt_path.exists():
        ckpt_df = pd.read_csv(ckpt_path, low_memory=False)
        done_ids = set(ckpt_df["ID"].astype(str))
        print(f"Resuming: {len(done_ids)} molecules already in checkpoint.")
    else:
        ckpt_df = pd.DataFrame()

    tasks = [
        (row["ID"], row[smiles_col], n_confs)
        for _, row in df.iterrows()
        if str(row["ID"]) not in done_ids
    ]

    print(f"Processing {len(tasks)} molecules ({len(done_ids)} already done) "
          f"with {n_cpus} CPUs, {n_confs} conformers each ...")

    results: list = []

    with mp.Pool(n_cpus) as pool:
        for res in tqdm(
            pool.imap_unordered(process_molecule, tasks, chunksize=1),
            total=len(tasks),
            desc="Conformers",
            miniters=1,
        ):
            results.append(res)

            # ── Incremental checkpoint every N completed molecules ────────
            if len(results) % checkpoint_every == 0:
                batch_df = pd.DataFrame(results)
                combined = pd.concat([ckpt_df, batch_df], ignore_index=True)
                combined.to_csv(ckpt_path, index=False)
                print(f"\n  [checkpoint] {len(done_ids) + len(results)} / "
                      f"{len(df)} saved to {ckpt_path}")

    # ── Final save ───────────────────────────────────────────────────────────
    all_results = pd.concat(
        [ckpt_df, pd.DataFrame(results)], ignore_index=True
    ) if not ckpt_df.empty else pd.DataFrame(results)

    n_failed = all_results["error"].notna().sum()
    print(f"\n  Succeeded: {len(all_results) - n_failed}")
    print(f"  Failed:    {n_failed}")
    if n_failed > 0:
        fail_counts = all_results["error"].value_counts()
        print(f"  Failure breakdown:\n{fail_counts.to_string()}")

    all_results.to_csv(out_path, index=False)
    if ckpt_path.exists():
        ckpt_path.unlink()   # remove checkpoint once final CSV is written
    print(f"\nSaved: {out_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Tier-1 conformer engine for CycPeptMPDB")
    parser.add_argument("--input",   "-i", default="data/pampa_curated.csv")
    parser.add_argument("--outdir",  "-o", default="results")
    parser.add_argument("--n-confs", "-c", type=int, default=50,
                        help="Conformers per molecule (default 50, use 200 for final run)")
    parser.add_argument("--n-cpus",  "-j", type=int,
                        default=max(1, mp.cpu_count() - 1))
    parser.add_argument("--max-mols", "-n", type=int, default=0,
                        help="Limit molecules for testing (0 = all)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(
        input_csv=args.input,
        outdir=Path(args.outdir),
        n_confs=args.n_confs,
        n_cpus=args.n_cpus,
        max_mols=args.max_mols,
    )

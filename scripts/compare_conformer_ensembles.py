# env: chameleon-calc
"""
compare_conformer_ensembles.py
------------------------------
Compares conformer ensembles across three sources for each reference compound:

  1. RDKit  — ETKDGv3 + MMFF94 (fast, CPU, seconds)
  2. CREST  — GFN2-xTB iMTD-GC (slow, our pipeline output)
  3. CREMP  — published CREST CHCl3 ensembles (pickle files, 6-mers only)
  4. Crystal / NMR — CIF experimental structures (CsA only, via gemmi)

For each source+molecule:
  - Boltzmann weights from energies (RT = 0.593 kcal/mol at 298.15 K)
  - Trim conformers below 1% cumulative Boltzmann population
  - Pairwise RMSD matrix → hierarchical clustering (Ward) → major macrostates
  - Per-conformer descriptors: PSA (3D SASA), intramolecular H-bonds, Rg
  - Boltzmann-averaged descriptors per ensemble

Outputs:
  results/ensemble_comparison_descriptors.csv   — one row per source × molecule
  results/ensemble_comparison_clusters.csv      — cluster centroids per source × molecule
  results/figures/ensemble_comparison_*.svg     — RMSD heatmaps, descriptor bar charts

Usage:
  python scripts/compare_conformer_ensembles.py
  python scripts/compare_conformer_ensembles.py --outdir results --n-rdkit 500
"""

import argparse
import json
import pickle
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors3D, rdMolDescriptors
from rdkit.Chem import rdFreeSASA
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

RDLogger.DisableLog("rdApp.*")

RT_KCAL = 0.592  # kcal/mol at 298.15 K

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO = Path(__file__).resolve().parents[1]

HEXPEP_SMILES = (
    "CC(C)C[C@@H]1NC(=O)[C@@H](CC(C)C)NC(=O)[C@@H](CC(C)C)NC(=O)"
    "[C@H](Cc2ccc(O)cc2)NC(=O)[C@@H]2CCCN2C(=O)[C@@H](CC(C)C)NC1=O"
)
CSA_SMILES = (
    "C/C=C/C[C@@H](C)[C@@H](O)[C@H]1C(=O)N[C@@H](CC)C(=O)N(C)CC(=O)"
    "N(C)[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N(C)[C@@H](CC(C)C)C(=O)"
    "N[C@@H](C)C(=O)N[C@H](C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)[C@@H](CC(C)C)"
    "C(=O)N(C)[C@@H](C(C)C)C(=O)N1C"
)

MOLECULES = {
    "HexPep": {
        "smiles": HEXPEP_SMILES,
        "crest_water": REPO / "results/crest_runs/run_20260429_141431_0_HexPep/water/crest/crest_conformers.xyz",
        "crest_chcl3": REPO / "results/crest_runs/run_20260429_141431_0_HexPep/chcl3/crest/crest_conformers.xyz",
        "cremp_id": "dL.dL.L.dL.P.Y",
        "cif_water": None,
        "cif_chcl3": None,
    },
    "CsA": {
        "smiles": CSA_SMILES,
        "crest_water": REPO / "data/CREST_CsA_20260512/crest_conformers.xyz",
        "crest_chcl3": None,  # still running on server
        "cremp_id": None,     # 11-mer, not in CREMP
        "cif_water": REPO / "data/experimental_structure_references_CsA/CsA_A1_combined_xray_neutron.cif",
        "cif_chcl3": REPO / "data/experimental_structure_references_CsA/CsA_C1_closed_DEKSAN_CCDC1138505.cif",
    },
}

# ── Descriptor functions (mirrors conformer_engine.py) ────────────────────────
_BONDI = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
    'S': 1.80, 'P': 1.80, 'F': 1.47, 'Cl': 1.75,
    'Br': 1.85, 'I': 1.98,
}
_POLAR = {'N', 'O', 'S', 'P'}
HB_DONOR    = Chem.MolFromSmarts("[N,O;!H0]")
HB_ACCEPTOR = Chem.MolFromSmarts("[N,O]")


def polar_sasa(mol, conf_id):
    try:
        radii = []
        for atom in mol.GetAtoms():
            sym = atom.GetSymbol()
            radii.append(_BONDI.get(sym, 1.50))
            atom.SetIntProp('SASAClass', 0 if sym in _POLAR else 1)
            atom.SetProp('SASAClassName', 'Polar' if sym in _POLAR else 'APolar')
        query = rdFreeSASA.MakeFreeSasaPolarAtomQuery()
        return round(rdFreeSASA.CalcSASA(mol, radii, confIdx=conf_id, query=query), 4)
    except Exception:
        return np.nan


def intramolecular_hbonds(mol, conf_id):
    try:
        conf = mol.GetConformer(conf_id)
        pos  = conf.GetPositions()
        donors    = [i for m in mol.GetSubstructMatches(HB_DONOR)    for i in m]
        acceptors = [i for m in mol.GetSubstructMatches(HB_ACCEPTOR) for i in m]
        count = 0
        for d in donors:
            for h in mol.GetAtomWithIdx(d).GetNeighbors():
                if h.GetAtomicNum() != 1:
                    continue
                hp = pos[h.GetIdx()]
                dp = pos[d]
                for a in acceptors:
                    if a == d:
                        continue
                    if len(Chem.GetShortestPath(mol, d, a)) < 6:
                        continue
                    ap = pos[a]
                    if np.linalg.norm(hp - ap) > 3.5:
                        continue
                    vhd = dp - hp
                    vha = ap - hp
                    cos = np.dot(vhd, vha) / (np.linalg.norm(vhd) * np.linalg.norm(vha) + 1e-9)
                    if np.degrees(np.arccos(np.clip(cos, -1, 1))) >= 120:
                        count += 1
        return count
    except Exception:
        return np.nan


def radius_of_gyration(mol, conf_id):
    try:
        return round(Descriptors3D.RadiusOfGyration(mol, confId=conf_id), 4)
    except Exception:
        return np.nan


def compute_descriptors(mol, conf_id):
    return {
        "psa":   polar_sasa(mol, conf_id),
        "hbonds": intramolecular_hbonds(mol, conf_id),
        "rg":    radius_of_gyration(mol, conf_id),
    }


# ── Boltzmann weights from Hartree energies ───────────────────────────────────
HARTREE_TO_KCAL = 627.509

def boltzmann_weights(energies_hartree):
    energies = np.array(energies_hartree, dtype=float)
    valid = np.isfinite(energies)
    if not valid.any():
        return np.ones(len(energies)) / len(energies)
    e_kcal = (energies - energies[valid].min()) * HARTREE_TO_KCAL
    w = np.zeros(len(energies))
    w[valid] = np.exp(-e_kcal[valid] / RT_KCAL)
    w /= w.sum()
    return w


def trim_by_population(weights, threshold=0.01):
    """Keep conformers until cumulative weight >= (1 - threshold). Returns boolean mask."""
    order = np.argsort(weights)[::-1]
    cumsum = np.cumsum(weights[order])
    keep_n = np.searchsorted(cumsum, 1.0 - threshold) + 1
    mask = np.zeros(len(weights), dtype=bool)
    mask[order[:keep_n]] = True
    return mask


# ── Parse multi-conformer XYZ (CREST format) ─────────────────────────────────
def parse_xyz_ensemble(xyz_path):
    """Returns list of (energy_hartree, coords_array, symbols_list)."""
    lines = Path(xyz_path).read_text().splitlines()
    conformers = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        try:
            n_atoms = int(line)
        except ValueError:
            i += 1
            continue
        i += 1
        energy = np.nan
        if i < len(lines):
            for tok in lines[i].strip().replace('=', ' ').split():
                try:
                    energy = float(tok)
                    break
                except ValueError:
                    continue
        i += 1
        symbols, coords = [], []
        for _ in range(n_atoms):
            if i >= len(lines):
                break
            parts = lines[i].split()
            if len(parts) >= 4:
                symbols.append(parts[0])
                coords.append([float(x) for x in parts[1:4]])
            i += 1
        if coords:
            conformers.append((energy, np.array(coords), symbols))
    return conformers


# ── Build RDKit mol with multiple conformers from XYZ data ────────────────────
def xyz_to_rdkit_mol(smiles, xyz_conformers):
    """
    Embed SMILES into RDKit mol and assign 3D coordinates from XYZ conformers.
    Uses substructure matching to map XYZ atom order to RDKit atom order.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())  # need one conformer for atom order

    ref_conf = mol.GetConformer(0)
    n_heavy = mol.GetNumAtoms()

    mol.RemoveAllConformers()
    added = 0
    for energy, coords, symbols in xyz_conformers:
        if len(coords) != n_heavy:
            continue
        conf = Chem.Conformer(n_heavy)
        for j, (x, y, z) in enumerate(coords):
            conf.SetAtomPosition(j, (x, y, z))
        conf.SetId(added)
        mol.AddConformer(conf, assignId=True)
        added += 1

    return mol, added


# ── RMSD matrix (heavy atoms, no alignment) ───────────────────────────────────
def rmsd_matrix(mol):
    mol_noh = Chem.RemoveHs(mol)
    n = mol_noh.GetNumConformers()
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            r = AllChem.GetBestRMS(mol_noh, mol_noh, i, j)
            mat[i, j] = mat[j, i] = r
    return mat


# ── Cluster conformers via Ward linkage on RMSD ───────────────────────────────
def cluster_conformers(rmsd_mat, weights, cutoff=2.0):
    """
    Returns cluster labels (1-indexed) for each conformer.
    cutoff: RMSD threshold in Angstroms for cluster merging.
    """
    if len(weights) < 2:
        return np.ones(len(weights), dtype=int)
    condensed = squareform(rmsd_mat)
    Z = linkage(condensed, method='ward')
    labels = fcluster(Z, t=cutoff, criterion='distance')
    return labels


def cluster_summary(labels, weights, desc_rows):
    """Returns DataFrame: one row per cluster with population and mean descriptors."""
    df = pd.DataFrame(desc_rows)
    df['weight'] = weights
    df['cluster'] = labels
    rows = []
    for cl in sorted(set(labels)):
        mask = df['cluster'] == cl
        pop  = df.loc[mask, 'weight'].sum()
        row  = {'cluster': cl, 'population': round(pop, 4)}
        for col in ['psa', 'hbonds', 'rg']:
            if col in df.columns:
                vals = df.loc[mask, col].values.astype(float)
                w    = df.loc[mask, 'weight'].values
                w    = w / w.sum()
                row[f'mean_{col}'] = round(float(np.nansum(vals * w)), 3)
        rows.append(row)
    return pd.DataFrame(rows).sort_values('population', ascending=False)


# ── Load CIF via gemmi → RDKit mol with one conformer ─────────────────────────
def load_cif_structure(cif_path, smiles):
    try:
        import gemmi
    except ImportError:
        print("  gemmi not available — skipping CIF")
        return None, None

    structure = gemmi.read_structure(str(cif_path))
    model = structure[0]

    coords_all, symbols_all = [], []
    for chain in model:
        for residue in chain:
            for atom in residue:
                pos = atom.pos
                coords_all.append([pos.x, pos.y, pos.z])
                symbols_all.append(atom.element.name if atom.element.name else 'C')

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

    n_heavy = mol.GetNumAtoms()
    heavy_coords = [c for c, s in zip(coords_all, symbols_all) if s != 'H']

    if len(heavy_coords) != n_heavy:
        print(f"  CIF atom count mismatch: {len(heavy_coords)} vs {n_heavy} — skipping")
        return None, None

    conf = Chem.Conformer(n_heavy)
    for j, (x, y, z) in enumerate(heavy_coords):
        conf.SetAtomPosition(j, (x, y, z))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol, heavy_coords


# ── Load CREMP pickle ─────────────────────────────────────────────────────────
def load_cremp_compound(compound_id, pickle_dir):
    pickle_dir = Path(pickle_dir)
    pkl_path = pickle_dir / f"{compound_id}.pickle"
    if not pkl_path.exists():
        # also try inside tar
        tar_path = pickle_dir.parent / "pickle.tar.gz"
        if tar_path.exists():
            with tarfile.open(tar_path) as tf:
                try:
                    member = tf.getmember(f"pickle/{compound_id}.pickle")
                    f = tf.extractfile(member)
                    return pickle.load(f)
                except KeyError:
                    return None
        return None
    with open(pkl_path, 'rb') as f:
        return pickle.load(f)


# ── Process one ensemble ──────────────────────────────────────────────────────
def process_ensemble(name, source, solvent, mol, weights, desc_label=""):
    """
    Given an RDKit mol (multiple conformers) and Boltzmann weights,
    compute descriptors, RMSD matrix, clusters, and summary stats.
    Returns (summary_row dict, cluster_df).
    """
    print(f"  [{name}] {source} {solvent}: {mol.GetNumConformers()} conformers")

    mask = trim_by_population(weights)
    conf_ids = [c.GetId() for c in mol.GetConformers()]
    kept_ids  = [cid for cid, m in zip(conf_ids, mask) if m]
    kept_w    = weights[mask]
    kept_w   /= kept_w.sum()

    print(f"    After 1% pop trim: {len(kept_ids)} conformers")

    # Build trimmed mol
    mol_trim = Chem.RWMol(mol)
    mol_trim.RemoveAllConformers()
    for cid in kept_ids:
        conf = mol.GetConformer(cid)
        new_conf = Chem.Conformer(mol.GetNumAtoms())
        for j in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(j)
            new_conf.SetAtomPosition(j, p)
        mol_trim.AddConformer(new_conf, assignId=True)
    mol_trim = mol_trim.GetMol()

    # Descriptors
    desc_rows = []
    for i, cid in enumerate(range(mol_trim.GetNumConformers())):
        desc_rows.append(compute_descriptors(mol_trim, cid))

    # RMSD matrix + clusters
    print(f"    Computing RMSD matrix...")
    rmat = rmsd_matrix(mol_trim)
    labels = cluster_conformers(rmat, kept_w, cutoff=2.0)
    n_clusters = len(set(labels))
    print(f"    {n_clusters} clusters (RMSD cutoff=2.0 Å)")

    clust_df = cluster_summary(labels, kept_w, desc_rows)
    clust_df.insert(0, 'molecule', name)
    clust_df.insert(1, 'source', source)
    clust_df.insert(2, 'solvent', solvent)

    # Boltzmann-averaged descriptors
    psa_arr = np.array([d['psa'] for d in desc_rows], dtype=float)
    hb_arr  = np.array([d['hbonds'] for d in desc_rows], dtype=float)
    rg_arr  = np.array([d['rg'] for d in desc_rows], dtype=float)

    summary = {
        'molecule':       name,
        'source':         source,
        'solvent':        solvent,
        'n_confs_raw':    mol.GetNumConformers(),
        'n_confs_kept':   len(kept_ids),
        'n_clusters':     n_clusters,
        'boltz_psa':      round(float(np.nansum(psa_arr * kept_w)), 2),
        'boltz_hbonds':   round(float(np.nansum(hb_arr  * kept_w)), 2),
        'boltz_rg':       round(float(np.nansum(rg_arr  * kept_w)), 3),
        'psa_std':        round(float(np.nanstd(psa_arr)), 2),
        'psa_spread':     round(float(np.nanmax(psa_arr) - np.nanmin(psa_arr)), 2),
        'lowen_psa':      round(float(psa_arr[0]), 2),
        'lowen_hbonds':   int(hb_arr[0]) if np.isfinite(hb_arr[0]) else np.nan,
    }

    return summary, clust_df


# ── Main ──────────────────────────────────────────────────────────────────────
def run(outdir: Path, n_rdkit: int, pickle_dir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    fig_dir = outdir / "figures"
    fig_dir.mkdir(exist_ok=True)

    all_summaries = []
    all_clusters  = []

    for mol_name, cfg in MOLECULES.items():
        smi = cfg["smiles"]
        print(f"\n{'='*60}")
        print(f"Molecule: {mol_name}")
        print(f"{'='*60}")

        # ── 1. RDKit ensemble ─────────────────────────────────────────────
        print("\n[RDKit] Generating conformers...")
        mol_rd = Chem.MolFromSmiles(smi)
        mol_rd = Chem.AddHs(mol_rd)
        params = AllChem.ETKDGv3()
        params.useMacrocycleTorsions = True
        params.numThreads = 0
        AllChem.EmbedMultipleConfs(mol_rd, numConfs=n_rdkit, params=params)
        AllChem.MMFFOptimizeMoleculeConfs(mol_rd, numThreads=0)

        ff_props = AllChem.MMFFGetMoleculeProperties(mol_rd)
        energies_rdkit = []
        for i in range(mol_rd.GetNumConformers()):
            ff = AllChem.MMFFGetMoleculeForceField(mol_rd, ff_props, confId=i)
            energies_rdkit.append(ff.CalcEnergy() if ff else np.nan)

        energies_rdkit = np.array(energies_rdkit)
        e_min = np.nanmin(energies_rdkit)
        w_rdkit = np.exp(-(energies_rdkit - e_min) / RT_KCAL)
        w_rdkit[~np.isfinite(w_rdkit)] = 0
        w_rdkit /= w_rdkit.sum()

        print(f"  {mol_rd.GetNumConformers()} conformers embedded")
        s, c = process_ensemble(mol_name, "RDKit_MMFF", "vacuum", mol_rd, w_rdkit)
        all_summaries.append(s)
        all_clusters.append(c)

        # ── 2. CREST ensembles ────────────────────────────────────────────
        for solvent_label, xyz_path in [("water", cfg["crest_water"]), ("chcl3", cfg["crest_chcl3"])]:
            if xyz_path is None or not Path(xyz_path).exists():
                print(f"\n[CREST/{solvent_label}] Not available — skipping")
                continue
            print(f"\n[CREST/{solvent_label}] Loading {xyz_path.name}...")
            raw = parse_xyz_ensemble(xyz_path)
            energies = np.array([e for e, _, _ in raw])
            w_crest = boltzmann_weights(energies)

            mol_crest, n_added = xyz_to_rdkit_mol(smi, raw)
            if n_added == 0:
                print(f"  No conformers loaded — skipping")
                continue
            # renorm weights to loaded conformers
            w_use = w_crest[:n_added]
            w_use = w_use / w_use.sum()

            s, c = process_ensemble(mol_name, "CREST", solvent_label, mol_crest, w_use)
            all_summaries.append(s)
            all_clusters.append(c)

        # ── 3. CREMP ensemble (CHCl3) ─────────────────────────────────────
        if cfg["cremp_id"] and pickle_dir.exists():
            print(f"\n[CREMP/chcl3] Loading {cfg['cremp_id']}...")
            obj = load_cremp_compound(cfg["cremp_id"], pickle_dir)
            if obj is not None:
                mol_cremp = obj.get('rd_mol')
                confs_meta = obj.get('conformers', [])
                if mol_cremp and mol_cremp.GetNumConformers() > 0:
                    w_cremp = np.array([c.get('boltzmannweight', 0) for c in confs_meta], dtype=float)
                    if len(w_cremp) != mol_cremp.GetNumConformers():
                        w_cremp = np.ones(mol_cremp.GetNumConformers())
                    w_cremp /= w_cremp.sum()
                    s, c = process_ensemble(mol_name, "CREMP", "chcl3", mol_cremp, w_cremp)
                    all_summaries.append(s)
                    all_clusters.append(c)
            else:
                print(f"  CREMP pickle not found for {cfg['cremp_id']}")

        # ── 4. CIF experimental structures ────────────────────────────────
        for env_label, cif_path in [("water_crystal", cfg["cif_water"]), ("chcl3_crystal", cfg["cif_chcl3"])]:
            if cif_path is None or not Path(cif_path).exists():
                continue
            print(f"\n[CIF/{env_label}] Loading {Path(cif_path).name}...")
            mol_cif, _ = load_cif_structure(cif_path, smi)
            if mol_cif is None:
                continue
            # Single structure → weight = 1.0
            w_cif = np.array([1.0])
            desc = compute_descriptors(mol_cif, 0)
            summary = {
                'molecule':     mol_name,
                'source':       'CIF_experimental',
                'solvent':      env_label,
                'n_confs_raw':  1,
                'n_confs_kept': 1,
                'n_clusters':   1,
                'boltz_psa':    desc['psa'],
                'boltz_hbonds': desc['hbonds'],
                'boltz_rg':     desc['rg'],
                'psa_std':      0.0,
                'psa_spread':   0.0,
                'lowen_psa':    desc['psa'],
                'lowen_hbonds': desc['hbonds'],
            }
            all_summaries.append(summary)
            clust = pd.DataFrame([{
                'molecule': mol_name, 'source': 'CIF_experimental', 'solvent': env_label,
                'cluster': 1, 'population': 1.0,
                'mean_psa': desc['psa'], 'mean_hbonds': desc['hbonds'], 'mean_rg': desc['rg'],
            }])
            all_clusters.append(clust)

    # ── Write outputs ─────────────────────────────────────────────────────
    summary_df = pd.DataFrame(all_summaries)
    cluster_df = pd.concat(all_clusters, ignore_index=True)

    summary_path = outdir / "ensemble_comparison_descriptors.csv"
    cluster_path = outdir / "ensemble_comparison_clusters.csv"
    summary_df.to_csv(summary_path, index=False)
    cluster_df.to_csv(cluster_path, index=False)

    print(f"\n{'='*60}")
    print(f"Summary written → {summary_path}")
    print(f"Clusters  written → {cluster_path}")
    print(f"\n{summary_df[['molecule','source','solvent','n_confs_kept','n_clusters','boltz_psa','boltz_hbonds','boltz_rg']].to_string(index=False)}")

    _plot(summary_df, fig_dir)


# ── Plotting ──────────────────────────────────────────────────────────────────
def _plot(df, fig_dir):
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use("Agg")
    except ImportError:
        print("matplotlib not available — skipping plots")
        return

    metrics = [
        ("boltz_psa",    "Boltzmann PSA (Å²)"),
        ("boltz_hbonds", "Boltzmann H-bonds"),
        ("boltz_rg",     "Boltzmann Rg (Å)"),
    ]

    for mol_name, grp in df.groupby("molecule"):
        fig, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 4))
        labels = [f"{r['source']}\n{r['solvent']}" for _, r in grp.iterrows()]

        for ax, (col, ylabel) in zip(axes, metrics):
            vals = grp[col].values
            colors = ["#4c72b0" if "water" in s else "#dd8452" if "chcl3" in s else "#55a868"
                      for s in grp["solvent"].values]
            ax.bar(range(len(vals)), vals, color=colors)
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels, fontsize=7, rotation=30, ha='right')
            ax.set_ylabel(ylabel, fontsize=9)
            ax.set_title(ylabel, fontsize=9)

        fig.suptitle(f"{mol_name} — Ensemble Descriptor Comparison", fontsize=11, y=1.01)
        fig.tight_layout()
        out = fig_dir / f"ensemble_comparison_{mol_name}.svg"
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
        print(f"Figure saved → {out}")


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="Compare RDKit / CREST / CREMP / CIF conformer ensembles")
    p.add_argument("--outdir",     default="results", type=Path)
    p.add_argument("--n-rdkit",    default=500, type=int, help="RDKit conformers to generate")
    p.add_argument("--pickle-dir", default="dependencies/pickle", type=Path,
                   help="Directory containing CREMP pickle files")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(
        outdir=REPO / args.outdir,
        n_rdkit=args.n_rdkit,
        pickle_dir=REPO / args.pickle_dir,
    )

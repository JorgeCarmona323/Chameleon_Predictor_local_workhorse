"""
colab_utils.py
--------------
Utility functions for Tier-2 CREST+ALPB pipeline running on Google Colab.

Upload this file to your Google Drive alongside the compound batch CSVs.
The notebook imports it at runtime.

Key functions:
  - smiles_to_xyz()          : RDKit SMILES → XYZ for CREST input
  - run_crest()              : execute CREST with ALPB solvation
  - parse_crest_best()       : extract lowest-energy conformer from CREST output
  - xyz_block_to_mol()       : apply CREST coordinates back onto RDKit template mol
  - polar_sasa()             : Bondi-radii 3D PSA (identical fix as conformer_engine.py)
  - intramolecular_hbonds()  : geometry-based H-bond counter
  - process_compound_crest() : full per-compound pipeline, returns result dict
"""

import os
import subprocess
import tempfile
import logging
import numpy as np
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors3D, rdFreeSASA
from rdkit.Chem.MolStandardize import rdMolStandardize

logging.basicConfig(level=logging.INFO, format="%(levelname)s [%(funcName)s]: %(message)s")
_log = logging.getLogger(__name__)

# ── Constants (identical to conformer_engine.py) ──────────────────────────────
_BONDI = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
    'S': 1.80, 'P': 1.80, 'F': 1.47, 'Cl': 1.75,
    'Br': 1.85, 'I': 1.98,
}
_POLAR_ELEMENTS = {'N', 'O', 'S', 'P'}

HB_DONOR_SMARTS    = Chem.MolFromSmarts("[N,O;!H0]")
HB_ACCEPTOR_SMARTS = Chem.MolFromSmarts("[N,O]")
HB_DIST_CUTOFF  = 3.0    # Å  H···A distance
HB_ANGLE_CUTOFF = 120.0  # degrees D-H···A at H


# ── SMILES → XYZ ─────────────────────────────────────────────────────────────

def smiles_to_xyz(smiles: str, mol_id, work_dir: Path) -> tuple[Path, Chem.Mol] | tuple[None, None]:
    """
    Generate a 3D XYZ file from SMILES using RDKit ETKDGv3 + MMFF94s.
    Returns (xyz_path, mol_with_H) or (None, None) on failure.

    The returned mol is the template used later to assign CREST coordinates —
    atom ordering in the XYZ matches the RDKit mol exactly.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        _log.warning("ID=%s: invalid SMILES", mol_id)
        return None, None

    # pH 7.4 standardization (same as conformer_engine.py)
    try:
        mol = rdMolStandardize.FragmentParent(mol)
        mol = rdMolStandardize.Uncharger(canonicalOrder=True).uncharge(mol)
        Chem.SanitizeMol(mol)
    except Exception as e:
        _log.warning("ID=%s: standardization failed (%s), using original", mol_id, e)

    mol_h = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useMacrocycleTorsions = True
    params.useSmallRingTorsions = True
    params.maxIterations = 2000

    if AllChem.EmbedMolecule(mol_h, params) != 0:
        _log.warning("ID=%s: embedding failed", mol_id)
        return None, None

    AllChem.MMFFOptimizeMolecule(mol_h, mmffVariant="MMFF94s", maxIters=2000)

    # Write XYZ
    conf = mol_h.GetConformer()
    xyz_lines = [str(mol_h.GetNumAtoms()), f"ID={mol_id}"]
    for atom in mol_h.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        xyz_lines.append(f"{atom.GetSymbol()}  {pos.x:.6f}  {pos.y:.6f}  {pos.z:.6f}")

    xyz_path = work_dir / f"{mol_id}_start.xyz"
    xyz_path.write_text("\n".join(xyz_lines))
    return xyz_path, mol_h


# ── Run CREST ─────────────────────────────────────────────────────────────────

def run_crest(xyz_path: Path, solvent: str, n_threads: int, run_dir: Path) -> Path | None:
    """
    Run CREST with ALPB solvation.

    solvent: 'water' (ε=80) or 'chcl3' (ε=4.8)
    Returns path to CREST output directory, or None on failure.

    Flags used:
      --alpb <solvent>  : analytical linearised Poisson-Boltzmann solvation
      --T <n>           : number of threads
      --quick           : reduced sampling (faster, still thorough for macrocycles)
      --mquick          : even faster; use if --quick times out on large compounds
    """
    solvent_dir = run_dir / solvent
    solvent_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "crest", str(xyz_path),
        "--alpb", solvent,
        "--T", str(n_threads),
        "--quick",
    ]

    log_path = solvent_dir / "crest.log"
    try:
        with open(log_path, "w") as logf:
            result = subprocess.run(
                cmd,
                cwd=str(solvent_dir),
                stdout=logf,
                stderr=subprocess.STDOUT,
                timeout=3600,  # 60 min hard timeout per solvent per compound
            )
        if result.returncode != 0:
            _log.warning("CREST failed (returncode=%d) for %s, solvent=%s",
                         result.returncode, xyz_path.stem, solvent)
            return None
        ensemble = solvent_dir / "crest_conformers.xyz"
        if not ensemble.exists():
            _log.warning("CREST ran but crest_conformers.xyz not found: %s %s",
                         xyz_path.stem, solvent)
            return None
        return solvent_dir
    except subprocess.TimeoutExpired:
        _log.warning("CREST timed out (>60 min): %s solvent=%s — try --mquick",
                     xyz_path.stem, solvent)
        return None
    except FileNotFoundError:
        _log.error("'crest' binary not found. Run the install cell first.")
        return None


# ── Parse CREST output ────────────────────────────────────────────────────────

def parse_crest_best(crest_dir: Path) -> str | None:
    """
    Extract the lowest-energy conformer XYZ block from crest_conformers.xyz.
    CREST sorts conformers by energy ascending — the first block is the best.
    Returns the raw XYZ block string or None.
    """
    ensemble_path = crest_dir / "crest_conformers.xyz"
    if not ensemble_path.exists():
        return None

    text = ensemble_path.read_text()
    blocks = []
    lines = text.strip().split("\n")
    i = 0
    while i < len(lines):
        try:
            n_atoms = int(lines[i].strip())
        except ValueError:
            i += 1
            continue
        block_lines = lines[i: i + n_atoms + 2]
        if len(block_lines) == n_atoms + 2:
            blocks.append("\n".join(block_lines))
        i += n_atoms + 2

    return blocks[0] if blocks else None


# ── XYZ block → RDKit Mol ─────────────────────────────────────────────────────

def xyz_block_to_mol(xyz_block: str, template_mol: Chem.Mol) -> Chem.Mol | None:
    """
    Apply CREST coordinates to the RDKit template mol.

    CREST preserves atom ordering from the input XYZ, so coordinates map
    directly onto the RDKit mol that generated that XYZ. We verify atom counts
    and element symbols match before applying.
    """
    lines = xyz_block.strip().split("\n")
    try:
        n_atoms = int(lines[0].strip())
    except ValueError:
        _log.warning("xyz_block_to_mol: could not parse atom count")
        return None

    if n_atoms != template_mol.GetNumAtoms():
        _log.warning("Atom count mismatch: CREST=%d, template=%d", n_atoms, template_mol.GetNumAtoms())
        return None

    coords = []
    elements = []
    for line in lines[2: 2 + n_atoms]:
        parts = line.split()
        if len(parts) < 4:
            return None
        elements.append(parts[0])
        coords.append((float(parts[1]), float(parts[2]), float(parts[3])))

    # Element-order sanity check
    for i, (atom, elem) in enumerate(zip(template_mol.GetAtoms(), elements)):
        if atom.GetSymbol() != elem:
            _log.warning("Element mismatch at atom %d: template=%s, CREST=%s",
                         i, atom.GetSymbol(), elem)
            return None

    rw = Chem.RWMol(template_mol)
    # Remove existing conformers and add new one
    rw.RemoveAllConformers()
    conf = Chem.Conformer(n_atoms)
    from rdkit.Geometry import rdGeometry
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, rdGeometry.Point3D(x, y, z))
    rw.AddConformer(conf, assignId=True)
    return rw.GetMol()


# ── PSA (identical fix to conformer_engine.py) ───────────────────────────────

def polar_sasa(mol: Chem.Mol, conf_id: int = 0) -> float:
    """
    Compute 3D polar SASA using Bondi radii + manual SASAClass assignment.

    Root cause of rdFreeSASA.classifyAtoms() failure on small molecules:
    all built-in classifiers are protein-residue-based and return zero
    for SMILES-derived molecules without PDB monomerInfo.

    Fix: manually assign SASAClass=Polar for N/O/S/P, APolar for all others.
    All atoms contribute their Bondi radii to mutual occlusion; only Polar
    atoms are summed. Consistent with heavy-atom Ertl (2000) TPSA convention.
    """
    try:
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
        return round(rdFreeSASA.CalcSASA(mol, radii, confIdx=conf_id, query=query), 4)
    except Exception as e:
        _log.debug("polar_sasa failed: %s", e)
        return np.nan


# ── H-bond counter (identical to conformer_engine.py) ────────────────────────

def intramolecular_hbonds(mol: Chem.Mol, conf_id: int = 0) -> int:
    """
    Count intramolecular H-bonds using 3D geometry.
    Criterion: H···A ≤ 3.0 Å, D-H···A angle ≥ 120°, topological path ≥ 6 atoms.
    """
    try:
        conf = mol.GetConformer(conf_id)
        pos = conf.GetPositions()

        donors    = [idx for match in mol.GetSubstructMatches(HB_DONOR_SMARTS) for idx in match]
        acceptors = [idx for match in mol.GetSubstructMatches(HB_ACCEPTOR_SMARTS) for idx in match]

        hb_count = 0
        for d_idx in donors:
            d_atom = mol.GetAtomWithIdx(d_idx)
            for h_atom in d_atom.GetNeighbors():
                if h_atom.GetAtomicNum() != 1:
                    continue
                h_idx = h_atom.GetIdx()
                h_pos, d_pos = pos[h_idx], pos[d_idx]

                for a_idx in acceptors:
                    if a_idx == d_idx:
                        continue
                    try:
                        if len(Chem.GetShortestPath(mol, d_idx, a_idx)) < 6:
                            continue
                    except Exception:
                        continue

                    a_pos = pos[a_idx]
                    if np.linalg.norm(h_pos - a_pos) > HB_DIST_CUTOFF:
                        continue

                    vec_hd = d_pos - h_pos
                    vec_ha = a_pos - h_pos
                    cos_a  = np.dot(vec_hd, vec_ha) / (
                        np.linalg.norm(vec_hd) * np.linalg.norm(vec_ha) + 1e-9
                    )
                    if np.degrees(np.arccos(np.clip(cos_a, -1, 1))) >= HB_ANGLE_CUTOFF:
                        hb_count += 1

        return hb_count
    except Exception as e:
        _log.debug("intramolecular_hbonds failed: %s", e)
        return np.nan


# ── Shape descriptors ─────────────────────────────────────────────────────────

def shape_descriptors(mol: Chem.Mol, conf_id: int = 0) -> dict:
    try:
        return {
            "Rg":             Descriptors3D.RadiusOfGyration(mol, confId=conf_id),
            "NPR1":           Descriptors3D.NPR1(mol, confId=conf_id),
            "NPR2":           Descriptors3D.NPR2(mol, confId=conf_id),
            "Asphericity":    Descriptors3D.Asphericity(mol, confId=conf_id),
            "SpherocityIndex":Descriptors3D.SpherocityIndex(mol, confId=conf_id),
        }
    except Exception:
        return {k: np.nan for k in ["Rg","NPR1","NPR2","Asphericity","SpherocityIndex"]}


# ── Main per-compound pipeline ────────────────────────────────────────────────

def process_compound_crest(
    mol_id,
    smiles: str,
    n_threads: int = 4,
    work_root: str | Path = "/tmp/crest_runs",
) -> dict:
    """
    Full Tier-2 CREST+ALPB pipeline for one compound.

    Steps:
      1. SMILES → 3D XYZ (RDKit ETKDGv3 + MMFF94s)
      2. CREST --alpb water  → lowest-energy aqueous conformer
      3. CREST --alpb chcl3  → lowest-energy membrane-mimetic conformer
      4. Compute PSA, HB, shape for each
      5. Compute Δ features

    Returns a flat dict with all features + 'error' key (None on success).
    """
    work_dir = Path(work_root) / str(mol_id)
    work_dir.mkdir(parents=True, exist_ok=True)

    base = {"ID": mol_id, "SMILES": smiles, "error": None}

    # Step 1: SMILES → XYZ
    xyz_path, template_mol = smiles_to_xyz(smiles, mol_id, work_dir)
    if xyz_path is None:
        return {**base, "error": "embed_failed"}

    # Steps 2 & 3: CREST in each solvent
    results_by_solvent = {}
    for solvent in ("water", "chcl3"):
        crest_dir = run_crest(xyz_path, solvent, n_threads, work_dir)
        if crest_dir is None:
            return {**base, "error": f"crest_failed_{solvent}"}

        xyz_block = parse_crest_best(crest_dir)
        if xyz_block is None:
            return {**base, "error": f"parse_failed_{solvent}"}

        mol_crest = xyz_block_to_mol(xyz_block, template_mol)
        if mol_crest is None:
            return {**base, "error": f"xyz_assign_failed_{solvent}"}

        psa  = polar_sasa(mol_crest, conf_id=0)
        hb   = intramolecular_hbonds(mol_crest, conf_id=0)
        shp  = shape_descriptors(mol_crest, conf_id=0)
        results_by_solvent[solvent] = {"psa": psa, "hb": hb, **shp}

    aq  = results_by_solvent["water"]
    mem = results_by_solvent["chcl3"]

    return {
        **base,
        # Aqueous (water) conformer
        "aq_psa3d":            aq["psa"],
        "aq_hb_count":         aq["hb"],
        "aq_Rg":               aq["Rg"],
        "aq_NPR1":             aq["NPR1"],
        "aq_NPR2":             aq["NPR2"],
        "aq_Asphericity":      aq["Asphericity"],
        "aq_SpherocityIndex":  aq["SpherocityIndex"],
        # Membrane (CHCl3) conformer
        "mem_psa3d":           mem["psa"],
        "mem_hb_count":        mem["hb"],
        "mem_Rg":              mem["Rg"],
        "mem_NPR1":            mem["NPR1"],
        "mem_NPR2":            mem["NPR2"],
        "mem_Asphericity":     mem["Asphericity"],
        "mem_SpherocityIndex": mem["SpherocityIndex"],
        # Δ features
        "delta_psa3d":         float(aq["psa"] - mem["psa"]) if not (np.isnan(aq["psa"]) or np.isnan(mem["psa"])) else np.nan,
        "delta_hb":            float(mem["hb"] - aq["hb"])   if not (np.isnan(aq["hb"])  or np.isnan(mem["hb"]))  else np.nan,
        "delta_Rg":            float(aq["Rg"] - mem["Rg"])   if not (np.isnan(aq["Rg"])  or np.isnan(mem["Rg"]))  else np.nan,
        "delta_NPR1":          float(aq["NPR1"] - mem["NPR1"]) if not (np.isnan(aq["NPR1"]) or np.isnan(mem["NPR1"])) else np.nan,
        "delta_NPR2":          float(aq["NPR2"] - mem["NPR2"]) if not (np.isnan(aq["NPR2"]) or np.isnan(mem["NPR2"])) else np.nan,
        "source":              "CREST+ALPB",
        "solvent_aq":          "water (eps=80)",
        "solvent_mem":         "chcl3 (eps=4.8)",
    }

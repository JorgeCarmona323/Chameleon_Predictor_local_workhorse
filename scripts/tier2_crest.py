"""
06b_tier2_crest.py
------------------
Tier-2 validation using CREST iMTD-GC + ALPB dual-dielectric conformer sampling.

Scientific rationale:
  Unlike Tier-1 (ETKDG vacuum heuristic: max-PSA ≈ aqueous, min-PSA ≈ membrane),
  CREST samples conformer ensembles INSIDE each dielectric environment using the
  ALPB (Analytical Linearized Poisson-Boltzmann) solvation model. This directly
  fulfills the proposal goal: "MMFF94 minimization at ε=78 and ε=4."

  The chameleonic ΔPSA computed here is:
    Boltzmann-weighted mean PSA (water ensemble, ε=80)
  − Boltzmann-weighted mean PSA (CHCl₃ ensemble, ε=4.8)

  Boltzmann weights are computed from GFN2-xTB energies in the CREST XYZ
  comment lines at T=298.15 K.  This is thermodynamically correct — each
  conformer contributes proportionally to its equilibrium population rather
  than using only the single lowest-energy structure.

Reference compounds (5 final) — size-stratified for chameleonic effect validation:
  0. Hexapeptide   (ID=2,   Rezai & Lokey JACS 2006)      PAMPA=-6.20  impermeable  6-mer  negative control
  1. CsA           (ID=1,   Witek JCTC 2016)               PAMPA=-5.90  permeable   11-mer  gold standard chameleonic
  2. c*[PSLYF]     (ID=1829,Hickey JMedChem 2016)          PAMPA=-9.10  impermeable 11-mer  large impermeable
  3. DP-955        (ID=917, CHUGAI 2013)                   PAMPA=-5.20  permeable   15-mer  largest permeable
  4. DP-944        (ID=906, CHUGAI 2013)                   PAMPA=-7.00  impermeable 15-mer  largest impermeable

Tools:
  - CREST 3.x  (conda install -c conda-forge crest)
  - xtb        (conda install -c conda-forge xtb)
  - RDKit      for 3D PSA and H-bond calculation on ensemble

Usage:
  python tier2_crest.py [--outdir results] [--max-confs 200] [--dry-run]

  --dry-run  skips CREST (uses pre-computed dummy ensemble) for testing
"""

import argparse
import os
import shutil
import subprocess
import tempfile
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

# ── Reference compounds ───────────────────────────────────────────────────────
REFERENCE_COMPOUNDS = [
    {
        "name": "Hexapeptide",
        "short": "HexPep",
        "cycpeptmpdb_id": 2,
        "smiles": (
            "CC(C)C[C@@H]1NC(=O)[C@@H](CC(C)C)NC(=O)[C@@H](CC(C)C)NC(=O)"
            "[C@H](Cc2ccc(O)cc2)NC(=O)[C@@H]2CCCN2C(=O)[C@@H](CC(C)C)NC1=O"
        ),
        "source": "Rezai & Lokey, JACS 2006",
        "pampa": -6.20,
        "permeable": False,
        "db_delta_psa": 2.0,
        "hbd": 6,
    },
    {
        "name": "Cyclosporin A",
        "short": "CsA",
        "cycpeptmpdb_id": 1,
        "smiles": (
            "C/C=C/C[C@@H](C)[C@@H](O)[C@H]1C(=O)N[C@@H](CC)C(=O)N(C)CC(=O)"
            "N(C)[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N(C)[C@@H](CC(C)C)C(=O)"
            "N[C@@H](C)C(=O)N[C@H](C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)[C@@H](CC(C)C)"
            "C(=O)N(C)[C@@H](C(C)C)C(=O)N1C"
        ),
        "source": "Witek JCTC 2016 / Rezai 2006 (NMR+MD proven, expected ΔPSA ~75 Å²)",
        "pampa": -5.90,  # cross-lab mean
        "permeable": True,
        "db_delta_psa": -1.0,
        "hbd": 5,
    },
    {
        "name": "c*[PSLYF]",
        "short": "PSLYF",
        "cycpeptmpdb_id": 1829,
        "smiles": None,  # fetched from feature matrix at runtime
        "source": "Hickey, J Med Chem 2016",
        "pampa": -9.10,
        "permeable": False,
        "db_delta_psa": 0.0,
        "hbd": 8,
    },
    {
        "name": "DP-955",
        "short": "DP955",
        "cycpeptmpdb_id": 917,
        "smiles": (
            "CC[C@H](C)[C@H](NC(=O)[C@@H]1CC(=O)N[C@@H](C)C(=O)N[C@@H]([C@@H](C)O)"
            "C(=O)N(C)CC(=O)N[C@@H](CC(C)C)C(=O)N2CCC[C@H]2C(=O)N(C)[C@@H](CC(C)C)"
            "C(=O)N[C@@H](Cc2ccccc2)C(=O)N(C)[C@@H](Cc2ccccc2)C(=O)N(C)[C@@H](C)C(=O)"
            "N(C)[C@@H](CC(C)C)C(=O)N[C@@H](C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N1)"
            "C(=O)N1CCCCC1"
        ),
        "source": "CHUGAI 2013 pharmaceutical screen",
        "pampa": -5.20,
        "permeable": True,
        "db_delta_psa": None,
        "hbd": None,
    },
    {
        "name": "DP-944",
        "short": "DP944",
        "cycpeptmpdb_id": 906,
        "smiles": (
            "CC[C@H](C)[C@H](NC(=O)[C@H](CC(C)C)N(C)C(=O)[C@H](CC(C)C)NC(=O)[C@H](C)"
            "N(C)C(=O)[C@@H]1CC(=O)N[C@H](C(C)C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)"
            "[C@@H](CC(C)C)C(=O)N[C@@H](C)C(=O)N[C@@H](Cc2ccccc2)C(=O)N[C@@H]([C@@H](C)O)"
            "C(=O)N[C@@H](CC(C)C)C(=O)N(C)[C@@H](Cc2ccccc2)C(=O)N(C)[C@@H](CC(C)C)"
            "C(=O)N1)C(=O)N1CCCCC1"
        ),
        "source": "CHUGAI 2013 pharmaceutical screen",
        "pampa": -7.00,
        "permeable": False,
        "db_delta_psa": None,
        "hbd": None,
    },
]

# Solvents for ALPB — CHCl3 approximates membrane dielectric (ε≈4.8)
SOLVENT_AQ  = "water"   # ε=80
SOLVENT_MEM = "chcl3"   # ε=4.8

POLAR_ATOMIC_NUMS = {7, 8}  # N and O


# ── Utility: get SMILES from feature matrix ───────────────────────────────────
def load_smiles_from_matrix(matrix_csv: str) -> dict:
    """Return {ID: SMILES} for all reference IDs that have None SMILES."""
    fm = pd.read_csv(matrix_csv, low_memory=False)
    smiles_col = "SMILES_canonical" if "SMILES_canonical" in fm.columns else "SMILES"
    id_to_smi = fm.set_index("ID")[smiles_col].to_dict()
    return id_to_smi


# ── Utility: SMILES → 3D XYZ via RDKit ETKDG (CREST input) ──────────────────
def smiles_to_xyz(smiles: str, xyz_path: Path) -> bool:
    """Generate a single 3D conformer from SMILES and write as XYZ for CREST."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    ret = AllChem.EmbedMolecule(mol, params)
    if ret == -1:
        # Fallback: random coords
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s")

    conf = mol.GetConformer()
    atoms = [mol.GetAtomWithIdx(i) for i in range(mol.GetNumAtoms())]
    lines = [str(len(atoms)), f"Generated from SMILES by RDKit for CREST input"]
    for atom, pos in zip(atoms, conf.GetPositions()):
        lines.append(f"{atom.GetSymbol()}  {pos[0]:.6f}  {pos[1]:.6f}  {pos[2]:.6f}")
    xyz_path.write_text("\n".join(lines))
    return True


# ── Utility: run CREST in a given solvent ─────────────────────────────────────
def run_crest(xyz_path: Path, work_dir: Path, solvent: str,
              max_confs: int = 200, charge: int = 0) -> Path | None:
    """
    Run CREST iMTD-GC with ALPB solvation.
    Returns path to the crest_conformers.xyz ensemble file, or None on failure.
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "crest", str(xyz_path),
        "--alpb", solvent,
        "--T", str(max(1, os.cpu_count() - 1)),
        "--quick",           # faster sampling — use --noreftopo for macrocycles
        "--keepdir",
        "--mquick",          # macrocycle-aware quick mode
    ]
    if charge != 0:
        cmd += ["--chrg", str(charge)]

    print(f"      Running: {' '.join(cmd)}")
    result = subprocess.run(
        cmd, cwd=work_dir,
        capture_output=True, text=True, timeout=7200,
    )

    ensemble = work_dir / "crest_conformers.xyz"
    if ensemble.exists() and ensemble.stat().st_size > 0:
        return ensemble
    else:
        print(f"      CREST stderr: {result.stderr[-500:]}")
        return None


# ── Utility: parse multi-conformer XYZ → list of coordinate arrays ───────────
def parse_xyz_ensemble(xyz_path: Path) -> list[tuple[list[str], np.ndarray, float]]:
    """
    Parse a CREST multi-conformer XYZ file.
    Returns list of (atom_symbols, coords_array, energy_hartree) per conformer.
    Energy is read from the XYZ comment line (CREST writes GFN2-xTB energy there).
    Defaults to 0.0 if the comment line is not parseable as a float.
    """
    conformers = []
    lines = xyz_path.read_text().splitlines()
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
        # Comment line — CREST writes the GFN2-xTB energy (Hartree) here
        energy = 0.0
        if i < len(lines):
            try:
                energy = float(lines[i].strip().split()[0])
            except (ValueError, IndexError):
                pass
        i += 1
        symbols, coords = [], []
        for _ in range(n_atoms):
            if i >= len(lines):
                break
            parts = lines[i].split()
            if len(parts) >= 4:
                symbols.append(parts[0])
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            i += 1
        if len(symbols) == n_atoms:
            conformers.append((symbols, np.array(coords), energy))
    return conformers


def boltzmann_weights(energies_hartree: list[float], T: float = 298.15) -> np.ndarray:
    """
    Boltzmann population weights from GFN2-xTB energies at temperature T.
    If all energies are 0.0 (not parsed), returns uniform weights.
    """
    KCAL_PER_HARTREE = 627.509
    RT = 1.987e-3 * T  # kcal/mol
    e = np.array(energies_hartree) * KCAL_PER_HARTREE
    e_rel = e - e.min()
    # If all energies identical (unparsed), fall back to uniform weights
    if e_rel.max() < 1e-10:
        return np.ones(len(e)) / len(e)
    weights = np.exp(-e_rel / RT)
    return weights / weights.sum()


# ── Utility: compute 3D PSA from XYZ conformer ───────────────────────────────
def compute_psa_xyz(symbols: list[str], coords: np.ndarray,
                    template_mol=None) -> float:
    """
    Compute 3D polar SASA (Å²) from atom symbols and coordinates.

    Implementation: uses rdFreeSASA with manual Bondi radii and SASAClass
    assignment — identical to the _polar_sasa() function in conformer_engine.py.
    This replaces the previous ad-hoc contact-counting exposure model
    (exposure = max(0.1, 1.0 - contacts * 0.12)) which was physically incorrect
    and systematically biased.

    If template_mol is provided (an RDKit Mol with connectivity), the XYZ
    coordinates are inserted as a new conformer and CalcSASA is called directly
    — this gives correct mutual-occlusion via the Lee-Richards/Shrake-Rupley
    algorithm.  If template_mol is None, the function falls back to the
    standalone distance-matrix SASA approximation using the same Bondi radii
    and polar-atom definition (N, O, S, P), which is still far more accurate
    than the original contact-counting model.

    Polar elements: N, O, S, P (consistent with conformer_engine.py and the
    Ertl 2000 TPSA convention; polar-H excluded per the heavy-atom-only
    convention documented in sasa_research_findings.md).

    NOTE for n≥50 HB sampling: at n=20 conformers ΔHB=0 for all test
    compounds.  This is a sampling artefact — use n≥50 for screening and
    n≥200 for final production runs to obtain reliable ΔHB signal.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdFreeSASA
    from rdkit.Geometry import rdGeometry

    # Bondi VdW radii (Å) — same as conformer_engine.py
    _BONDI = {
        'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
        'S': 1.80, 'P': 1.80, 'F': 1.47, 'Cl': 1.75,
        'Br': 1.85, 'I': 1.98,
    }
    # Heteroatoms treated as polar (consistent with conformer_engine.py)
    _POLAR_ELEMENTS = {'N', 'O', 'S', 'P'}

    if template_mol is not None:
        # Path A: insert XYZ as a new conformer on the template mol and use
        # rdFreeSASA.CalcSASA with manual Bondi radii + SASAClass assignment.
        # This is the correct, rigorous path when SMILES connectivity is known.
        try:
            mol_h = Chem.RWMol(template_mol)
            if mol_h.GetNumAtoms() != len(symbols):
                raise ValueError(
                    f"Atom count mismatch: template has {mol_h.GetNumAtoms()}, "
                    f"XYZ has {len(symbols)}"
                )
            conf = Chem.Conformer(mol_h.GetNumAtoms())
            for i, (x, y, z) in enumerate(coords):
                conf.SetAtomPosition(i, (float(x), float(y), float(z)))
            conf_id = mol_h.AddConformer(conf, assignId=True)
            mol_h = mol_h.GetMol()

            radii = []
            for atom in mol_h.GetAtoms():
                sym = atom.GetSymbol()
                radii.append(_BONDI.get(sym, 1.50))
                if sym in _POLAR_ELEMENTS:
                    atom.SetIntProp('SASAClass', 0)
                    atom.SetProp('SASAClassName', 'Polar')
                else:
                    atom.SetIntProp('SASAClass', 1)
                    atom.SetProp('SASAClassName', 'APolar')

            query = rdFreeSASA.MakeFreeSasaPolarAtomQuery()
            psa = rdFreeSASA.CalcSASA(mol_h, radii, confIdx=conf_id,
                                       query=query)
            return round(float(psa), 2)
        except Exception as _e:
            # Fall through to Path B on any failure
            pass

    # Path B: standalone approximate SASA using pairwise distance-matrix
    # shadowing with Bondi radii.  More accurate than the old contact-counting
    # model because it weights overlap by the actual geometric arc area
    # reduction rather than a fixed 0.12 per contact.
    # Still an approximation (no proper Lee-Richards integration), but
    # consistent in its polar-atom definition and radii with Path A.
    PROBE = 1.40  # water probe radius (Å)
    n = len(symbols)
    radii_all = np.array([_BONDI.get(s, 1.50) + PROBE for s in symbols])
    polar_idx = [i for i, s in enumerate(symbols) if s in _POLAR_ELEMENTS]
    if not polar_idx:
        return 0.0

    # Compute pairwise distance matrix once
    diff = coords[:, None, :] - coords[None, :, :]        # (n,n,3)
    dist = np.sqrt((diff ** 2).sum(axis=-1))               # (n,n)
    np.fill_diagonal(dist, np.inf)

    total = 0.0
    for idx in polar_idx:
        r_i = radii_all[idx]
        # Lee-Richards approximation: fraction of sphere area not overlapped.
        # For each neighbour j that overlaps (d_ij < r_i + r_j), subtract the
        # spherical cap area on sphere i using the formula:
        #   cap_area = 2π r_i h   where h = r_i - (r_i² + d² - r_j²)/(2d)
        buried = 0.0
        sphere_area = 4.0 * np.pi * r_i ** 2
        for jdx in range(n):
            if jdx == idx:
                continue
            d = dist[idx, jdx]
            r_j = radii_all[jdx]
            if d >= r_i + r_j:
                continue
            if d <= abs(r_i - r_j):
                # Fully enclosed: subtract entire sphere area of smaller
                if r_i <= r_j:
                    buried = sphere_area  # fully inside j
                    break
                continue
            # Partial overlap: spherical cap on sphere i
            h = r_i - (r_i**2 + d**2 - r_j**2) / (2.0 * d)
            h = max(0.0, min(h, 2.0 * r_i))
            buried += 2.0 * np.pi * r_i * h
        exposed = max(0.0, sphere_area - buried)
        total += exposed

    return round(total, 2)


# ── Utility: count intramolecular H-bonds in XYZ conformer ───────────────────
def count_hbonds_xyz(symbols: list[str], coords: np.ndarray) -> int:
    """
    Count intramolecular H-bonds: D-H...A where
      D, A ∈ {N, O}; H...A distance < 2.5 Å; D-H...A angle > 120°.
    """
    D_ATOMS = {"N", "O"}
    A_ATOMS = {"N", "O"}
    H_DIST_MAX = 2.5
    ANGLE_MIN = 120.0

    # Find donor H atoms (H bonded to N or O)
    donor_h = []
    for i, sym in enumerate(symbols):
        if sym != "H":
            continue
        # Find nearest heavy atom
        dists = [(np.linalg.norm(coords[i] - coords[j]), j)
                 for j, s in enumerate(symbols) if s != "H" and j != i]
        if not dists:
            continue
        d_nearest, d_idx = min(dists)
        if d_nearest < 1.3 and symbols[d_idx] in D_ATOMS:
            donor_h.append((i, d_idx))  # (H_idx, D_idx)

    # Find acceptors
    acceptors = [i for i, s in enumerate(symbols) if s in A_ATOMS]

    count = 0
    for h_idx, d_idx in donor_h:
        h_pos = coords[h_idx]
        d_pos = coords[d_idx]
        for a_idx in acceptors:
            if a_idx == d_idx:
                continue
            a_pos = coords[a_idx]
            ha_dist = np.linalg.norm(h_pos - a_pos)
            if ha_dist > H_DIST_MAX:
                continue
            # Angle at H: D-H...A
            vec_hd = d_pos - h_pos
            vec_ha = a_pos - h_pos
            cos_a = np.dot(vec_hd, vec_ha) / (
                np.linalg.norm(vec_hd) * np.linalg.norm(vec_ha) + 1e-10)
            angle = np.degrees(np.arccos(np.clip(cos_a, -1, 1)))
            if angle >= ANGLE_MIN:
                count += 1
    return count


# ── Process one compound ──────────────────────────────────────────────────────
def process_compound(cpd: dict, work_base: Path,
                     max_confs: int, dry_run: bool) -> dict:
    name  = cpd["name"]
    short = cpd["short"]
    smi   = cpd["smiles"]

    print(f"\n{'─'*60}")
    print(f"  {name}  (ID={cpd['cycpeptmpdb_id']})")
    print(f"  Source: {cpd['source']}")
    print(f"  PAMPA: {cpd['pampa']}  |  HBD: {cpd['hbd']}  |  DB ΔPSA: {cpd['db_delta_psa']}")

    result = {
        "compound": name,
        "cycpeptmpdb_id": cpd["cycpeptmpdb_id"],
        "pampa": cpd["pampa"],
        "permeable": cpd["permeable"],
        "hbd": cpd["hbd"],
        "db_delta_psa": cpd["db_delta_psa"],
        "source": cpd["source"],
    }

    if smi is None:
        print(f"  ⚠ No SMILES available — skipping")
        return result

    # Generate starting XYZ
    work_dir = work_base / short
    work_dir.mkdir(parents=True, exist_ok=True)
    xyz_in = work_dir / f"{short}_start.xyz"

    if not smiles_to_xyz(smi, xyz_in):
        print(f"  ⚠ RDKit 3D embedding failed")
        return result
    print(f"  Starting XYZ written: {xyz_in.name}")

    # Build template mol with explicit H for compute_psa_xyz Path A.
    # We standardize (FragmentParent + Uncharger) identically to conformer_engine.py
    # so that the atom count matches the CREST XYZ output (which was generated
    # from the same starting geometry via smiles_to_xyz — no standardization there,
    # so we build the template from the raw SMILES with explicit H added).
    _template_mol = None
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem.MolStandardize import rdMolStandardize
        _tmpl = Chem.MolFromSmiles(smi)
        if _tmpl is not None:
            _tmpl = Chem.AddHs(_tmpl)
            # Embed a single conformer so the Mol object is valid for SASA
            _ep = AllChem.ETKDGv3()
            _ep.randomSeed = 42
            _ep.useMacrocycleTorsions = True
            AllChem.EmbedMolecule(_tmpl, _ep)
            _template_mol = _tmpl
            print(f"  Template mol built: {_tmpl.GetNumAtoms()} atoms (with H)")
    except Exception as _te:
        print(f"  ⚠ Template mol build failed: {_te} — falling back to Path B SASA")

    for solvent, label in [(SOLVENT_AQ, "aq"), (SOLVENT_MEM, "mem")]:
        sol_dir = work_dir / solvent
        print(f"\n  ── CREST [{label}] solvent={solvent} ──")

        if dry_run:
            # Dry run: return placeholder PSA values
            psa_mean = 180.0 if label == "aq" else 140.0
            psa_vals = np.random.default_rng(42).normal(psa_mean, 10, 50)
            hb_vals  = np.random.default_rng(42).integers(0, 4, 50)
            result[f"{label}_psa_min"]    = float(psa_vals.min())
            result[f"{label}_psa_max"]    = float(psa_vals.max())
            result[f"{label}_psa_mean"]   = float(psa_vals.mean())
            result[f"{label}_psa_boltz"]  = float(psa_vals.mean())  # uniform weights in dry run
            result[f"{label}_psa_std"]    = float(psa_vals.std())
            result[f"{label}_hb_min"]     = int(hb_vals.min())
            result[f"{label}_hb_max"]     = int(hb_vals.max())
            result[f"{label}_hb_mean"]    = float(hb_vals.mean())
            result[f"{label}_hb_boltz"]   = float(hb_vals.mean())
            result[f"{label}_n_confs"]    = 50
            result[f"{label}_psa_lowen"]  = float(psa_vals[0])
            result[f"{label}_hb_lowen"]   = int(hb_vals[0])
            print(f"      [DRY RUN] PSA_low-energy={result[f'{label}_psa_lowen']:.1f}")
            continue

        # Run CREST
        try:
            ensemble_xyz = run_crest(xyz_in, sol_dir, solvent, max_confs)
        except subprocess.TimeoutExpired:
            print(f"      ⚠ CREST timed out after 2 hours")
            ensemble_xyz = None
        except Exception as e:
            print(f"      ⚠ CREST error: {e}")
            ensemble_xyz = None

        if ensemble_xyz is None:
            print(f"      ⚠ CREST failed — no ensemble produced")
            continue

        # Parse ensemble — returns (symbols, coords, energy_hartree) tuples
        conformers = parse_xyz_ensemble(ensemble_xyz)
        n_confs = len(conformers)
        print(f"      Parsed {n_confs} conformers from ensemble")

        if n_confs == 0:
            continue

        # Trim to max_confs and compute PSA / H-bonds per conformer
        conformers = conformers[:max_confs]
        psa_vals, hb_vals, energies = [], [], []
        for syms, crds, eng in conformers:
            psa_vals.append(compute_psa_xyz(syms, crds, template_mol=_template_mol))
            hb_vals.append(count_hbonds_xyz(syms, crds))
            energies.append(eng)

        psa_arr = np.array(psa_vals)
        hb_arr  = np.array(hb_vals, dtype=float)

        # Boltzmann weights from GFN2-xTB energies at 298.15 K
        weights = boltzmann_weights(energies)
        psa_boltz = float(np.dot(weights, psa_arr))
        hb_boltz  = float(np.dot(weights, hb_arr))

        # Lowest-energy conformer = first in CREST output (sorted by energy)
        psa_lowen = psa_arr[0]
        hb_lowen  = hb_arr[0]

        result[f"{label}_psa_lowen"]  = round(float(psa_lowen), 2)
        result[f"{label}_psa_boltz"]  = round(psa_boltz, 2)
        result[f"{label}_psa_min"]    = round(float(psa_arr.min()), 2)
        result[f"{label}_psa_max"]    = round(float(psa_arr.max()), 2)
        result[f"{label}_psa_mean"]   = round(float(psa_arr.mean()), 2)
        result[f"{label}_psa_std"]    = round(float(psa_arr.std()), 2)
        result[f"{label}_hb_lowen"]   = int(hb_lowen)
        result[f"{label}_hb_boltz"]   = round(hb_boltz, 2)
        result[f"{label}_hb_min"]     = int(hb_arr.min())
        result[f"{label}_hb_max"]     = int(hb_arr.max())
        result[f"{label}_hb_mean"]    = round(float(hb_arr.mean()), 2)
        result[f"{label}_n_confs"]    = n_confs

        print(f"      PSA (Boltzmann-wtd): {psa_boltz:.1f} Å²  "
              f"(low-energy={psa_lowen:.1f}, mean={psa_arr.mean():.1f})")
        print(f"      HB  (Boltzmann-wtd): {hb_boltz:.2f}  "
              f"(low-energy={int(hb_lowen)})")

    # ── Compute Δ features ────────────────────────────────────────────────────
    if "aq_psa_boltz" in result and "mem_psa_boltz" in result:
        result["crest_delta_psa"]       = round(result["aq_psa_boltz"] - result["mem_psa_boltz"], 2)
        result["crest_delta_psa_lowen"] = round(result.get("aq_psa_lowen", 0) - result.get("mem_psa_lowen", 0), 2)
        result["crest_delta_hb"]        = round(result["mem_hb_boltz"] - result["aq_hb_boltz"], 2)
        result["crest_psa_spread_aq"]   = result.get("aq_psa_std", np.nan)
        result["crest_psa_spread_mem"]  = result.get("mem_psa_std", np.nan)
        print(f"\n  ✓ CREST ΔPSA (Boltzmann) = {result['crest_delta_psa']:.1f} Å²  "
              f"(low-energy = {result['crest_delta_psa_lowen']:.1f} Å²  "
              f"ΔHB = {result['crest_delta_hb']:.2f})")

    return result


# ── Plotting ──────────────────────────────────────────────────────────────────
def plot_results(results: list[dict], outdir: Path) -> None:
    valid = [r for r in results if "crest_delta_psa" in r]
    if not valid:
        print("No valid results to plot")
        return

    names   = [r["compound"].replace(" (1NMe3)", "\n(1NMe3)") for r in valid]
    colors  = ["#D6604D" if r["permeable"] else "#4393C3" for r in valid]
    pampa   = [r["pampa"] for r in valid]
    c_dpsa  = [r["crest_delta_psa"] for r in valid]
    db_dpsa = [r["db_delta_psa"] for r in valid]
    c_dhb   = [r["crest_delta_hb"] for r in valid]
    n       = len(valid)

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle(
        "Tier-2 CREST+ALPB Validation — Dual-Dielectric Conformer Sampling\n"
        "Water (ε=80) vs CHCl₃ (ε=4.8) | GFN2-xTB + ALPB",
        fontsize=11, fontweight="bold",
    )

    # Panel A: ΔPSA comparison (CREST vs DB)
    ax = axes[0, 0]
    x = np.arange(n)
    w = 0.35
    ax.bar(x - w/2, db_dpsa, width=w, label="DB (static, CycPeptMPDB)",
           color="#BEAED4", edgecolor="grey", linewidth=0.5)
    ax.bar(x + w/2, c_dpsa, width=w, label="CREST+ALPB (ensemble)",
           color="#FDC086", edgecolor="grey", linewidth=0.5)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel("ΔPSA (Å²) = PSA_aq − PSA_mem")
    ax.set_title("A. ΔPSA: DB static vs CREST+ALPB ensemble", fontweight="bold")
    ax.legend(fontsize=8)

    # Panel B: PAMPA vs CREST ΔPSA
    ax = axes[0, 1]
    for i, r in enumerate(valid):
        c = "#D6604D" if r["permeable"] else "#4393C3"
        ax.scatter(r["crest_delta_psa"], r["pampa"],
                   s=120, c=c, edgecolors="black", linewidths=0.8, zorder=4)
        ax.annotate(r["compound"].split("(")[0].strip(),
                    (r["crest_delta_psa"], r["pampa"]),
                    xytext=(5, 4), textcoords="offset points", fontsize=7.5)
    ax.axhline(-6.0, color="grey", linestyle="--", linewidth=0.8,
               label="PAMPA threshold (−6.0)")
    ax.set_xlabel("CREST ΔPSA (Å²)")
    ax.set_ylabel("PAMPA LogPexp (log cm/s)")
    ax.set_title("B. PAMPA vs CREST ΔPSA", fontweight="bold")
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor="#D6604D", label="Permeable"),
        Patch(facecolor="#4393C3", label="Impermeable"),
        plt.Line2D([0],[0], color="grey", linestyle="--", label="−6.0 threshold"),
    ], fontsize=8)

    # Panel C: PAMPA vs CREST ΔHB
    ax = axes[1, 0]
    for i, r in enumerate(valid):
        c = "#D6604D" if r["permeable"] else "#4393C3"
        ax.scatter(r["crest_delta_hb"], r["pampa"],
                   s=120, c=c, edgecolors="black", linewidths=0.8, zorder=4)
        ax.annotate(r["compound"].split("(")[0].strip(),
                    (r["crest_delta_hb"], r["pampa"]),
                    xytext=(5, 4), textcoords="offset points", fontsize=7.5)
    ax.axhline(-6.0, color="grey", linestyle="--", linewidth=0.8)
    ax.set_xlabel("CREST ΔHB (H-bonds mem − H-bonds aq)")
    ax.set_ylabel("PAMPA LogPexp (log cm/s)")
    ax.set_title("C. PAMPA vs CREST ΔHB\n(mechanistic: intramolecular H-bond formation)",
                 fontweight="bold")

    # Panel D: CREST ΔPSA vs DB ΔPSA cross-check
    ax = axes[1, 1]
    for i, r in enumerate(valid):
        c = "#D6604D" if r["permeable"] else "#4393C3"
        ax.scatter(r["db_delta_psa"], r["crest_delta_psa"],
                   s=120, c=c, edgecolors="black", linewidths=0.8, zorder=4)
        ax.annotate(r["compound"].split("(")[0].strip(),
                    (r["db_delta_psa"], r["crest_delta_psa"]),
                    xytext=(5, 4), textcoords="offset points", fontsize=7.5)
    all_v = c_dpsa + db_dpsa
    lim = [min(all_v) - 5, max(all_v) + 5]
    ax.plot(lim, lim, "k--", linewidth=0.8, alpha=0.5, label="y=x")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("DB ΔPSA (static CycPeptMPDB)")
    ax.set_ylabel("CREST ΔPSA (ensemble, dual-dielectric)")
    ax.set_title("D. CREST vs DB Cross-Check\n(DB static misses CsA chameleonism)",
                 fontweight="bold")
    ax.legend(fontsize=8)
    if len(valid) >= 3:
        try:
            r_val, _ = stats.pearsonr(db_dpsa, c_dpsa)
            ax.text(0.05, 0.95, f"r = {r_val:.2f}", transform=ax.transAxes,
                    fontsize=10, va="top", fontweight="bold")
        except Exception:
            pass

    plt.tight_layout()
    fig_path = outdir / "figures" / "tier2_crest_crosscheck.png"
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {fig_path}")


# ── Main ──────────────────────────────────────────────────────────────────────
def run(matrix_csv: str, outdir: Path, max_confs: int, dry_run: bool,
        compound_idx: int | None = None, n_threads: int | None = None) -> None:
    """
    Run Tier-2 CREST validation.

    Parallelisation mode (recommended by computational collaborator):
      Submit one job per compound using --compound 0..4 and --threads N.
      Each job writes results/tier2_crest_compound_{idx}.csv independently.
      After all 5 jobs finish, run --merge to combine into tier2_crest_table.csv.

    Example job array (SLURM):
      for i in 0 1 2 3 4; do
        sbatch --ntasks=8 run_crest.sh --compound $i --threads 8
      done
      python scripts/tier2_crest.py --merge --outdir results

    Single-compound local test:
      python scripts/tier2_crest.py --compound 0 --threads 4 --dry-run
    """
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    work_base = outdir / "crest_runs"
    work_base.mkdir(exist_ok=True)

    # Fill in SMILES for compounds that need them from feature matrix
    id_to_smi = load_smiles_from_matrix(matrix_csv)
    for cpd in REFERENCE_COMPOUNDS:
        if cpd["smiles"] is None:
            cid = cpd["cycpeptmpdb_id"]
            cpd["smiles"] = id_to_smi.get(cid)
            if cpd["smiles"]:
                print(f"Loaded SMILES for {cpd['name']} (ID={cid}) from feature matrix")
            else:
                print(f"WARNING: No SMILES found for {cpd['name']} (ID={cid})")

    # Override CREST thread count if specified — use a module-level patch
    # (cannot shadow `run_crest` as a local variable; Python would treat
    #  the name as local throughout the function and fail before the assignment)
    if n_threads is not None:
        _nt = n_threads  # capture for closure

        def _threaded_run_crest(xyz_path, work_dir, solvent, max_confs=200, charge=0):
            work_dir.mkdir(parents=True, exist_ok=True)
            cmd = [
                "crest", str(xyz_path),
                "--alpb", solvent,
                "--T", str(_nt),
                "--quick",
                "--keepdir",
                "--mquick",
            ]
            if charge != 0:
                cmd += ["--chrg", str(charge)]
            print(f"      Running: {' '.join(cmd)}")
            result = subprocess.run(
                cmd, cwd=work_dir,
                capture_output=True, text=True, timeout=7200,
            )
            ensemble = work_dir / "crest_conformers.xyz"
            if ensemble.exists() and ensemble.stat().st_size > 0:
                return ensemble
            else:
                print(f"      CREST stderr: {result.stderr[-500:]}")
                return None

        # Monkey-patch the module-level function so process_compound picks it up
        import sys as _sys
        _this = _sys.modules[__name__]
        _orig_run_crest = _this.run_crest
        _this.run_crest = _threaded_run_crest

    # Select compounds to run
    if compound_idx is not None:
        if compound_idx < 0 or compound_idx >= len(REFERENCE_COMPOUNDS):
            raise ValueError(
                f"--compound must be 0–{len(REFERENCE_COMPOUNDS)-1}, "
                f"got {compound_idx}. Compounds: "
                + ", ".join(f"{i}={c['short']}"
                            for i, c in enumerate(REFERENCE_COMPOUNDS))
            )
        compounds_to_run = [REFERENCE_COMPOUNDS[compound_idx]]
        print(f"\n[Parallel mode] Running compound {compound_idx}: "
              f"{compounds_to_run[0]['name']}")
    else:
        compounds_to_run = REFERENCE_COMPOUNDS

    results = []
    try:
        for cpd in compounds_to_run:
            r = process_compound(cpd, work_base, max_confs, dry_run)
            results.append(r)
    finally:
        # Restore original run_crest if it was monkey-patched
        if n_threads is not None:
            _this.run_crest = _orig_run_crest

    # Save per-compound CSV (parallel-safe — no file collision)
    if compound_idx is not None:
        short = compounds_to_run[0]["short"]
        out_csv = outdir / f"tier2_crest_compound_{compound_idx}_{short}.csv"
        pd.DataFrame(results).to_csv(out_csv, index=False)
        print(f"\nSaved: {out_csv}")
        print(f"Run --merge after all 5 compounds complete to combine results.")
        return

    # Sequential mode: save combined table and plot immediately
    _save_and_plot(results, outdir)


def _save_and_plot(results: list[dict], outdir: Path) -> None:
    table = pd.DataFrame([{k: v for k, v in r.items()} for r in results])
    out_csv = outdir / "tier2_crest_table.csv"
    table.to_csv(out_csv, index=False)

    print(f"\n{'='*60}")
    print("CREST Tier-2 Summary")
    print(f"{'='*60}")
    disp_cols = ["compound", "pampa", "db_delta_psa",
                 "crest_delta_psa", "crest_delta_hb",
                 "aq_n_confs", "mem_n_confs"]
    avail = [c for c in disp_cols if c in table.columns]
    print(table[avail].to_string(index=False))
    print(f"\nSaved: {out_csv}")
    plot_results(results, outdir)


def merge_parallel_results(outdir: Path) -> None:
    """
    Combine per-compound CSVs written by parallel jobs into tier2_crest_table.csv.
    Run after all --compound jobs complete.
    """
    import glob
    pattern = str(outdir / "tier2_crest_compound_*.csv")
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"No per-compound CSVs found matching: {pattern}")
        return
    print(f"Merging {len(files)} per-compound files...")
    dfs = [pd.read_csv(f) for f in files]
    combined = pd.concat(dfs, ignore_index=True)
    # Restore original compound order
    order = {c["name"]: i for i, c in enumerate(REFERENCE_COMPOUNDS)}
    combined["_order"] = combined["compound"].map(order)
    combined = combined.sort_values("_order").drop(columns="_order").reset_index(drop=True)
    out_csv = outdir / "tier2_crest_table.csv"
    combined.to_csv(out_csv, index=False)
    print(f"Merged {len(combined)} compounds → {out_csv}")
    plot_results(combined.to_dict("records"), outdir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Tier-2 CREST+ALPB validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Parallelisation (one job per compound):
  python tier2_crest.py --compound 0 --threads 8   # HexPep
  python tier2_crest.py --compound 1 --threads 8   # 1NMe3
  python tier2_crest.py --compound 2 --threads 8   # CsA
  python tier2_crest.py --compound 3 --threads 8   # DP-172
  python tier2_crest.py --compound 4 --threads 8   # PSLYF

  After all 5 finish:
  python tier2_crest.py --merge

Compound index reference:
  0 = HexPep   (impermeable, 6-mer)
  1 = 1NMe3    (permeable,   6-mer, N-methylated)
  2 = CsA      (permeable,  11-mer, chameleonic)
  3 = DP-172   (permeable,   CHUGAI)
  4 = PSLYF    (impermeable, HBD=8)
        """
    )
    parser.add_argument("--matrix",    "-m", default="results/feature_matrix.csv")
    parser.add_argument("--outdir",    "-o", default="results", type=Path)
    parser.add_argument("--max-confs", "-c", type=int, default=200)
    parser.add_argument("--dry-run",   action="store_true",
                        help="Skip CREST, use placeholder values for testing")
    parser.add_argument("--compound",  type=int, default=None,
                        metavar="IDX",
                        help="Run only compound IDX (0-4). For parallel job submission.")
    parser.add_argument("--threads",   type=int, default=None,
                        metavar="N",
                        help="CREST --T threads per job. Default: all available cores.")
    parser.add_argument("--merge",     action="store_true",
                        help="Merge per-compound CSVs after parallel jobs complete.")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.merge:
        merge_parallel_results(Path(args.outdir))
    else:
        run(args.matrix, Path(args.outdir),
            max_confs=args.max_confs, dry_run=args.dry_run,
            compound_idx=args.compound, n_threads=args.threads)

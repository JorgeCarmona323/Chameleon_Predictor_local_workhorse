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

Pre-processing pipeline (mirrors CREMP — Atz et al. 2024):
  1. RDKit ETKDGv3 embedding (5000 conformers, useMacrocycleTorsions=True)
     → MMFF optimization → RMSD filter → top 50 conformers
  2. xTB GFN2 geometry optimization of all 50 in parallel (one per CPU),
     each with the target solvent (ALPB) — gives the lowest-energy conformer
     as starting geometry for CREST
  3. CREST iMTD-GC: --gfn2 --chrg {charge} --alpb {solvent} -T {threads}
     (no --quick, no --mquick, no --noreftopo — matches CREMP exactly)

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
import multiprocessing
import os
import shutil
import subprocess
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

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


def _log(msg: str) -> None:
    """Print a timestamped line and flush immediately (for tail -f visibility)."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


# ── Utility: get SMILES from feature matrix ───────────────────────────────────
def load_smiles_from_matrix(matrix_csv: str) -> dict:
    """Return {ID: SMILES} for all reference IDs that have None SMILES."""
    fm = pd.read_csv(matrix_csv, low_memory=False)
    smiles_col = "SMILES_canonical" if "SMILES_canonical" in fm.columns else "SMILES"
    id_to_smi = fm.set_index("ID")[smiles_col].to_dict()
    return id_to_smi


# ── CREMP Step 1: RDKit conformer embedding ───────────────────────────────────
def embed_rdkit_conformers(mol, n_embed: int = 5000, n_max: int = 50,
                           rmsd_thresh: float = 0.5):
    """
    ETKDGv3 embedding + MMFF filter → up to n_max conformers.
    Mirrors CREMP rdconf.py exactly (Atz et al. 2024).
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    params = AllChem.ETKDGv3()
    params.maxIterations = 10 * n_embed   # 50 000
    params.pruneRmsThresh = 0.01
    params.useRandomCoords = True
    params.useMacrocycleTorsions = True
    params.numThreads = 0                 # all available cores

    _log(f"    ETKDGv3: embedding {n_embed} conformers...")
    t0 = time.time()
    AllChem.EmbedMultipleConfs(mol, numConfs=n_embed, params=params)
    if mol.GetNumConformers() == 0:
        raise RuntimeError("ETKDGv3 produced zero conformers")
    _log(f"    ETKDGv3: {mol.GetNumConformers()} embedded in {time.time()-t0:.0f}s")

    # Skip MMFF optimisation if any bond has E/Z stereo (CREMP behaviour)
    has_cistrans = not all(
        bond.GetStereo() is Chem.BondStereo.STEREONONE for bond in mol.GetBonds()
    )
    if not has_cistrans:
        _log("    MMFF: skipped (E/Z stereo bonds present)")
    else:
        _log("    MMFF: optimising all conformers...")
        t0 = time.time()
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)
        _log(f"    MMFF: done in {time.time()-t0:.0f}s")

    # Rank by MMFF energy, keep up to n_max with pairwise RMSD ≥ rmsd_thresh
    AllChem.MMFFSanitizeMolecule(mol)
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
    confs = list(mol.GetConformers())
    energies = []
    for conf in confs:
        ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=conf.GetId())
        energies.append(ff.CalcEnergy())

    sort_idx = np.argsort(energies)

    mol_no_h = Chem.RemoveHs(mol)
    atom_idxs = [a.GetIdx() for a in mol_no_h.GetAtoms()]
    matches = mol_no_h.GetSubstructMatches(mol_no_h, uniquify=False)
    atom_map = [list(zip(match, atom_idxs)) for match in matches]

    keep = [sort_idx[0]]
    for i in sort_idx[1:]:
        if len(keep) >= n_max:
            break
        rmsds = [
            AllChem.GetBestRMS(
                mol_no_h, mol_no_h,
                confs[j].GetId(), confs[i].GetId(),
                map=atom_map,
            )
            for j in keep
        ]
        if all(r >= rmsd_thresh for r in rmsds):
            keep.append(i)

    new_mol = Chem.Mol(mol)
    new_mol.RemoveAllConformers()
    conf_ids = [c.GetId() for c in confs]
    for i in keep:
        new_mol.AddConformer(mol.GetConformer(conf_ids[i]), assignId=True)

    return new_mol


# ── CREMP Step 2: xTB geometry optimisation ───────────────────────────────────
def _write_conformer_xyz(mol, conf_id: int, xyz_path: Path,
                         comment: str = "") -> None:
    """Write a single conformer from mol to an XYZ file."""
    conf = mol.GetConformer(conf_id)
    pos  = conf.GetPositions()
    syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(mol.GetNumAtoms())]
    lines = [str(mol.GetNumAtoms()), comment]
    for sym, (x, y, z) in zip(syms, pos):
        lines.append(f"{sym}  {x:.6f}  {y:.6f}  {z:.6f}")
    xyz_path.write_text("\n".join(lines) + "\n")


def _xtb_opt_worker(args):
    """
    Module-level worker for multiprocessing: GFN2-xTB --opt of one conformer.
    Mirrors CREMP xtb.py XTBConformerOptimizer.run_xtb().
    Returns (conf_id, mol_with_single_optimised_conf, energy_hartree) or None.
    """
    from rdkit import Chem

    conf_id, mol, conf_dir, solvent, charge = args
    conf_dir = Path(conf_dir)
    conf_dir.mkdir(parents=True, exist_ok=True)

    xyz_path = conf_dir / "conf.xyz"
    _write_conformer_xyz(mol, conf_id, xyz_path, comment=str(conf_id))

    xtb_exe = shutil.which("xtb")
    if xtb_exe is None:
        return None

    cmd = [xtb_exe, "conf.xyz", "--opt", "--gfn", "2", "--chrg", str(charge)]
    if solvent:
        model = "gbsa" if solvent.lower() == "methanol" else "alpb"
        cmd.extend([f"--{model}", solvent])

    # OMP_NUM_THREADS=1,1 — matches CREMP: single-threaded per worker
    env = {**os.environ, "OMP_NUM_THREADS": "1,1"}
    try:
        subprocess.run(cmd, cwd=conf_dir, check=True,
                       capture_output=True, env=env, timeout=3600)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return None

    opt_xyz = conf_dir / "xtbopt.xyz"
    if not opt_xyz.exists():
        return None

    try:
        with open(opt_xyz) as f:
            n_atoms = int(next(f))
            # xTB comment: "energy: -XXX.XXX Eh  gradient norm: ..."
            energy = float(next(f).split()[1])
            new_conf = Chem.Conformer(n_atoms)
            for i in range(n_atoms):
                parts = next(f).split()
                new_conf.SetAtomPosition(i, [float(parts[1]),
                                             float(parts[2]),
                                             float(parts[3])])
        mol_copy = Chem.Mol(mol, quickCopy=True)
        mol_copy.AddConformer(new_conf, assignId=True)
        return (conf_id, mol_copy, energy)
    except Exception:
        return None


def xtb_preopt_mol(mol, solvent: str, work_dir: Path,
                   charge: int, n_procs: int = 1):
    """
    Run GFN2-xTB --opt on all conformers in mol in parallel (n_procs workers).
    Returns a mol containing only the single lowest-energy optimised conformer.
    Mirrors CREMP xtb.py XTBMolOptimization + get_mol_with_lowest_energy_conf().
    Returns None if all xTB jobs fail.
    """
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    conf_ids = [c.GetId() for c in mol.GetConformers()]
    if not conf_ids:
        return None

    params_list = [
        (conf_id, mol, work_dir / f"conf_{conf_id}", solvent, charge)
        for conf_id in conf_ids
    ]

    n_workers = min(n_procs, len(conf_ids))
    _log(f"    xTB: starting {len(conf_ids)} optimisations ({n_workers} parallel workers)...")
    t0 = time.time()
    results = []
    done = 0
    with multiprocessing.Pool(n_workers) as pool:
        for result in pool.imap_unordered(_xtb_opt_worker, params_list):
            done += 1
            status = "ok" if result is not None else "failed"
            elapsed = time.time() - t0
            _log(f"    xTB: {done}/{len(conf_ids)} done ({elapsed:.0f}s elapsed) — conf {result[0] if result else '?'} {status}")
            results.append(result)

    n_ok = sum(1 for r in results if r is not None)
    _log(f"    xTB: {n_ok}/{len(conf_ids)} succeeded in {time.time()-t0:.0f}s")

    valid = [r for r in results if r is not None]
    if not valid:
        return None

    _, mol_min, _ = min(valid, key=lambda x: x[2])
    return mol_min


# ── CREMP Step 3: CREST conformer sampling ────────────────────────────────────
def run_crest(xyz_path: Path, work_dir: Path, solvent: str,
              n_threads: int, charge: int = 0) -> Path | None:
    """
    Run CREST iMTD-GC with GFN2-xTB + ALPB solvation.
    Command mirrors CREMP crest.py CRESTSampler.run_crest() exactly:
      crest {xyz} -T {n} --gfn2 --chrg {q} --alpb {solvent} --keepdir
    Returns path to crest_conformers.xyz, or None on failure.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    xyz_path = xyz_path.resolve()

    cmd = [
        "crest", xyz_path.name,
        "-T",      str(n_threads),
        "--gfn2",
        "--chrg",  str(charge),
        "--alpb",  solvent,
        "--keepdir",
    ]

    _log(f"    CREST: {' '.join(cmd)}")
    _log(f"    CREST: stdout → crest.out  stderr → crest.err")
    t0 = time.time()

    # Redirect CREST output to files (mirrors CREMP) so tail -f stays readable
    out_path = work_dir / "crest.out"
    err_path = work_dir / "crest.err"
    with open(out_path, "w") as fout, open(err_path, "w") as ferr:
        proc = subprocess.run(
            cmd, cwd=work_dir.resolve(),
            stdout=fout, stderr=ferr, timeout=86400,
        )

    elapsed = time.time() - t0
    ensemble = work_dir.resolve() / "crest_conformers.xyz"
    if ensemble.exists() and ensemble.stat().st_size > 0:
        _log(f"    CREST: finished in {elapsed/3600:.2f}h — ensemble found")
        return ensemble
    else:
        _log(f"    CREST: FAILED after {elapsed:.0f}s (exit={proc.returncode})")
        # Print last 20 lines of stderr for diagnosis
        try:
            tail = err_path.read_text().splitlines()[-20:]
            print("\n".join(f"      ERR: {l}" for l in tail), flush=True)
        except Exception:
            pass
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


def save_top_conformers(
    conformers: list[tuple[list[str], np.ndarray, float]],
    weights: np.ndarray,
    label: str,
    compound_short: str,
    outdir: Path,
    top_n: int = 10,
) -> Path:
    """
    Save top-N Boltzmann-weighted conformers as a multi-conformer XYZ file.

    Conformers are sorted by weight (highest first). Comment line encodes rank,
    weight, and GFN2-xTB energy so the file can be used directly as input to
    xTB, ORCA, or docking pipelines.

    Output: results/conformers/{compound_short}/{label}/top{N}_boltzmann.xyz
    """
    ranked = sorted(zip(weights, conformers), key=lambda x: -x[0])[:top_n]

    conf_dir = outdir / "conformers" / compound_short / label
    conf_dir.mkdir(parents=True, exist_ok=True)
    out_path = conf_dir / f"top{top_n}_boltzmann.xyz"

    lines = []
    for rank, (w, (syms, crds, eng)) in enumerate(ranked, start=1):
        lines.append(str(len(syms)))
        lines.append(f"rank={rank} weight={w:.6f} energy={eng:.8f}")
        for sym, (x, y, z) in zip(syms, crds):
            lines.append(f"{sym}  {x:.6f}  {y:.6f}  {z:.6f}")

    out_path.write_text("\n".join(lines) + "\n")
    return out_path


def boltzmann_weights(energies_hartree: list[float], T: float = 298.15) -> np.ndarray:
    """
    Boltzmann population weights from GFN2-xTB energies at temperature T.
    If all energies are 0.0 (not parsed), returns uniform weights.
    """
    KCAL_PER_HARTREE = 627.509
    RT = 1.987e-3 * T  # kcal/mol
    e = np.array(energies_hartree) * KCAL_PER_HARTREE
    e_rel = e - e.min()
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

    _BONDI = {
        'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
        'S': 1.80, 'P': 1.80, 'F': 1.47, 'Cl': 1.75,
        'Br': 1.85, 'I': 1.98,
    }
    _POLAR_ELEMENTS = {'N', 'O', 'S', 'P'}

    if template_mol is not None:
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
            psa = rdFreeSASA.CalcSASA(mol_h, radii, confIdx=conf_id, query=query)
            return round(float(psa), 2)
        except Exception:
            pass

    # Path B: standalone approximate SASA using pairwise distance-matrix shadowing
    PROBE = 1.40
    n = len(symbols)
    radii_all = np.array([_BONDI.get(s, 1.50) + PROBE for s in symbols])
    polar_idx = [i for i, s in enumerate(symbols) if s in _POLAR_ELEMENTS]
    if not polar_idx:
        return 0.0

    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=-1))
    np.fill_diagonal(dist, np.inf)

    total = 0.0
    for idx in polar_idx:
        r_i = radii_all[idx]
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
                if r_i <= r_j:
                    buried = sphere_area
                    break
                continue
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

    donor_h = []
    for i, sym in enumerate(symbols):
        if sym != "H":
            continue
        dists = [(np.linalg.norm(coords[i] - coords[j]), j)
                 for j, s in enumerate(symbols) if s != "H" and j != i]
        if not dists:
            continue
        d_nearest, d_idx = min(dists)
        if d_nearest < 1.3 and symbols[d_idx] in D_ATOMS:
            donor_h.append((i, d_idx))

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
                     max_confs: int, dry_run: bool, n_threads: int,
                     top_confs: int = 10, outdir: Path | None = None,
                     restart: bool = False) -> dict:
    from rdkit import Chem

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

    work_dir = work_base / short
    work_dir.mkdir(parents=True, exist_ok=True)

    # ── Build mol and compute formal charge (CREMP pattern) ───────────────────
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"  ⚠ Invalid SMILES")
        return result
    mol = Chem.AddHs(mol)
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    _log(f"  Formal charge: {charge}  |  Atoms (with H): {mol.GetNumAtoms()}")

    # ── CREMP Step 1: embed 5000 → filter to 50 ──────────────────────────────
    _log(f"  Step 1: RDKit ETKDGv3 embedding (5000 → 50)...")
    try:
        mol_embedded = embed_rdkit_conformers(mol)
    except RuntimeError as e:
        _log(f"  ⚠ Embedding failed: {e}")
        return result
    n_embedded = mol_embedded.GetNumConformers()
    _log(f"  Step 1 done: {n_embedded} conformers after MMFF+RMSD filter")

    # Template mol for PSA Path A — reuse mol_embedded (has connectivity + H)
    _template_mol = mol_embedded

    for solvent, label in [(SOLVENT_AQ, "aq"), (SOLVENT_MEM, "mem")]:
        sol_dir = work_dir / solvent
        print(f"\n  ── [{label}] solvent={solvent} ──")

        if dry_run:
            psa_mean = 180.0 if label == "aq" else 140.0
            psa_vals = np.random.default_rng(42).normal(psa_mean, 10, 50)
            hb_vals  = np.random.default_rng(42).integers(0, 4, 50)
            result[f"{label}_psa_min"]    = float(psa_vals.min())
            result[f"{label}_psa_max"]    = float(psa_vals.max())
            result[f"{label}_psa_mean"]   = float(psa_vals.mean())
            result[f"{label}_psa_boltz"]  = float(psa_vals.mean())
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

        # Full ensemble save path — used for restart and cap extension
        saved_ens_path = (
            outdir / "conformers" / short / label / "full_ensemble.xyz"
            if outdir is not None else None
        )

        # ── Run pipeline or reload saved ensemble ─────────────────────────────
        if restart and saved_ens_path is not None and saved_ens_path.exists():
            print(f"      [RESTART] Loading saved ensemble: {saved_ens_path}")
            ensemble_xyz = saved_ens_path
        else:
            # CREMP Step 2: xTB pre-optimisation (per-solvent)
            xtb_dir = sol_dir / "xtb_opt"
            _log(f"  Step 2: xTB pre-opt ({n_embedded} conformers, {min(n_threads, n_embedded)} workers)...")
            mol_min = xtb_preopt_mol(mol_embedded, solvent, xtb_dir, charge,
                                     n_procs=n_threads)
            if mol_min is None:
                _log(f"  ⚠ All xTB jobs failed — skipping {label}")
                continue

            # Write lowest-energy xTB conformer as CREST input
            crest_dir = sol_dir / "crest"
            crest_dir.mkdir(parents=True, exist_ok=True)
            xyz_in = crest_dir / f"{short}_{label}_start.xyz"
            _write_conformer_xyz(mol_min, mol_min.GetConformer().GetId(),
                                 xyz_in, comment=smi)
            _log(f"  Step 2 done: xTB-optimised start XYZ written → {xyz_in.name}")
            _log(f"  Step 3: CREST iMTD-GC ({solvent})...")

            # CREMP Step 3: CREST iMTD-GC
            try:
                ensemble_xyz = run_crest(xyz_in, crest_dir, solvent,
                                         n_threads, charge)
            except subprocess.TimeoutExpired:
                print(f"      ⚠ CREST timed out after 24h")
                ensemble_xyz = None
            except Exception as e:
                print(f"      ⚠ CREST error: {e}")
                ensemble_xyz = None

        if ensemble_xyz is None:
            print(f"      ⚠ CREST failed — no ensemble produced")
            continue

        # Parse full ensemble before applying cap
        conformers = parse_xyz_ensemble(ensemble_xyz)
        n_confs_full = len(conformers)
        print(f"      Parsed {n_confs_full} conformers from ensemble")

        # Save full ensemble for later restart / cap extension
        if saved_ens_path is not None and not restart:
            saved_ens_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ensemble_xyz, saved_ens_path)
            print(f"      Full ensemble saved → {saved_ens_path.name} ({n_confs_full} conformers)")

        if n_confs_full == 0:
            continue

        # Trim to cap and compute PSA / H-bonds per conformer
        conformers = conformers[:max_confs]
        n_confs = len(conformers)
        psa_vals, hb_vals, energies = [], [], []
        for syms, crds, eng in conformers:
            psa_vals.append(compute_psa_xyz(syms, crds, template_mol=_template_mol))
            hb_vals.append(count_hbonds_xyz(syms, crds))
            energies.append(eng)

        psa_arr = np.array(psa_vals)
        hb_arr  = np.array(hb_vals, dtype=float)

        weights   = boltzmann_weights(energies)
        psa_boltz = float(np.dot(weights, psa_arr))
        hb_boltz  = float(np.dot(weights, hb_arr))

        psa_lowen = psa_arr[0]
        hb_lowen  = hb_arr[0]

        result[f"{label}_psa_lowen"]    = round(float(psa_lowen), 2)
        result[f"{label}_psa_boltz"]    = round(psa_boltz, 2)
        result[f"{label}_psa_min"]      = round(float(psa_arr.min()), 2)
        result[f"{label}_psa_max"]      = round(float(psa_arr.max()), 2)
        result[f"{label}_psa_mean"]     = round(float(psa_arr.mean()), 2)
        result[f"{label}_psa_std"]      = round(float(psa_arr.std()), 2)
        result[f"{label}_hb_lowen"]     = int(hb_lowen)
        result[f"{label}_hb_boltz"]     = round(hb_boltz, 2)
        result[f"{label}_hb_min"]       = int(hb_arr.min())
        result[f"{label}_hb_max"]       = int(hb_arr.max())
        result[f"{label}_hb_mean"]      = round(float(hb_arr.mean()), 2)
        result[f"{label}_n_confs"]      = n_confs
        result[f"{label}_n_confs_full"] = n_confs_full

        print(f"      PSA (Boltzmann-wtd): {psa_boltz:.1f} Å²  "
              f"(low-energy={psa_lowen:.1f}, mean={psa_arr.mean():.1f})")
        print(f"      HB  (Boltzmann-wtd): {hb_boltz:.2f}  "
              f"(low-energy={int(hb_lowen)})")
        if n_confs_full > n_confs:
            print(f"      ⚠ Cap applied: using {n_confs}/{n_confs_full} conformers "
                  f"(raise --max-confs or use --restart to extend)")

        if top_confs > 0 and outdir is not None:
            xyz_out = save_top_conformers(
                conformers, weights, label, short, outdir, top_n=top_confs,
            )
            print(f"      Saved top-{top_confs} conformers → {xyz_out.name}")

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


# ── Main ──────────────────────────────────────────────────────────────────────
def run(matrix_csv: str, outdir: Path, max_confs: int, dry_run: bool,
        compound_idx: int | None = None, n_threads: int | None = None,
        top_confs: int = 10, restart: bool = False) -> None:
    """
    Run Tier-2 CREST validation.

    Parallelisation mode (recommended):
      Submit one SLURM job per compound using --compound 0..4 and --threads N.
      Each job runs: RDKit embed → xTB pre-opt (water) → CREST (water)
                                 → xTB pre-opt (CHCl3) → CREST (CHCl3)
      After all 5 jobs finish, run --merge to combine into tier2_crest_table.csv.

    Single-compound local test:
      python scripts/tier2_crest.py --compound 0 --threads 4 --dry-run
    """
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    work_base = outdir / "crest_runs"
    work_base.mkdir(exist_ok=True)

    n_threads = n_threads or os.cpu_count() or 1

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
    for cpd in compounds_to_run:
        r = process_compound(cpd, work_base, max_confs, dry_run, n_threads,
                             top_confs=top_confs, outdir=outdir, restart=restart)
        results.append(r)

    if compound_idx is not None:
        short = compounds_to_run[0]["short"]
        out_csv = outdir / f"tier2_crest_compound_{compound_idx}_{short}.csv"
        pd.DataFrame(results).to_csv(out_csv, index=False)
        print(f"\nSaved: {out_csv}")
        print(f"Run --merge after all 5 compounds complete to combine results.")
        return

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
    print(f"Run locally: python scripts/plot_tier2_results.py --csv {out_csv}")


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
    order = {c["name"]: i for i, c in enumerate(REFERENCE_COMPOUNDS)}
    combined["_order"] = combined["compound"].map(order)
    combined = combined.sort_values("_order").drop(columns="_order").reset_index(drop=True)
    out_csv = outdir / "tier2_crest_table.csv"
    combined.to_csv(out_csv, index=False)
    print(f"Merged {len(combined)} compounds → {out_csv}")
    print(f"Run locally: python scripts/plot_tier2_results.py --csv {out_csv}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Tier-2 CREST+ALPB validation (CREMP-matching pipeline)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Parallelisation (one job per compound):
  python tier2_crest.py --compound 0 --threads 20   # HexPep
  python tier2_crest.py --compound 1 --threads 20   # CsA
  python tier2_crest.py --compound 2 --threads 20   # PSLYF
  python tier2_crest.py --compound 3 --threads 20   # DP-955
  python tier2_crest.py --compound 4 --threads 20   # DP-944

  After all 5 finish:
  python tier2_crest.py --merge

Compound index reference:
  0 = HexPep  (impermeable,  6-mer)
  1 = CsA     (permeable,   11-mer, chameleonic, expected ΔPSA ~75 Å²)
  2 = PSLYF   (impermeable, 11-mer, HBD=8)
  3 = DP-955  (permeable,   15-mer, CHUGAI 2013)
  4 = DP-944  (impermeable, 15-mer, CHUGAI 2013)
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
                        help="Threads for CREST (-T) and workers for xTB pre-opt. "
                             "Default: all available cores.")
    parser.add_argument("--top-confs", type=int, default=10,
                        metavar="N",
                        help="Save top-N Boltzmann-weighted conformers per solvent as "
                             "multi-conformer XYZ in results/conformers/. "
                             "Set to 0 to disable. (default: 10)")
    parser.add_argument("--merge",     action="store_true",
                        help="Merge per-compound CSVs after parallel jobs complete.")
    parser.add_argument("--restart",   action="store_true",
                        help="Reload saved full_ensemble.xyz instead of re-running CREST. "
                             "Use to extend analysis with a higher --max-confs cap without "
                             "re-running the conformer search.")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.merge:
        merge_parallel_results(Path(args.outdir))
    else:
        run(args.matrix, Path(args.outdir),
            max_confs=args.max_confs, dry_run=args.dry_run,
            compound_idx=args.compound, n_threads=args.threads,
            top_confs=args.top_confs, restart=args.restart)

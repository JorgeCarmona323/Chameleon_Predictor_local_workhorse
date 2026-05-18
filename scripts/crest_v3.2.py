"""
crest_v3.2.py
-------------
Tier-2 validation using CREST iMTD-GC + ALPB dual-dielectric conformer sampling.

Changes from v3.1:
  - Consolidated output structure: results/runs/run_{ts}_{idx}_{short}/water|mem/
    (no separate conformers/ directory — everything lives in the run folder)
  - Parses crest.out → full CREMP-format thermodynamics in ensemble.json
  - Exports ensemble.sdf (all conformers, RDKit connectivity) + ensemble.json per solvent
  - from __future__ import annotations for Python 3.9 compat
  - No max_confs cap by default (use --max-confs to override)
  - Keeps --keepdir for METADYN/NORMMD crash recovery

Output layout per run:
  results/runs/run_{timestamp}_{idx}_{short}/
    water/
      xtb_opt/          ← xTB intermediate files
      crest/            ← CREST output (METADYN*, NORMMD*, crest_conformers.xyz)
      ensemble.xyz      ← full conformer ensemble
      ensemble.sdf      ← CREMP-format: all conformers with RDKit connectivity
      ensemble.json     ← CREMP-format: thermodynamics + PSA/HB per conformer
    mem/
      (same structure)
    {short}_results.csv ← PSA / HB / ΔPSA summary

Usage:
  python crest_v3.2.py --compound 1 --threads 20 --outdir results
"""

from __future__ import annotations

import argparse
import json
import multiprocessing
import os
import re
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
        "source": "Witek JCTC 2016",
        "pampa": -5.90,
        "permeable": True,
        "hbd": 5,
    },
    {
        "name": "Cyclosporin O",
        "short": "CsO",
        "cycpeptmpdb_id": None,
        "smiles": (
            "CCC[C@H]1C(=O)N(CC(=O)N([C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H]"
            "(C(=O)N[C@@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H]"
            "(C(=O)N1)CC(C)C)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)"
            "CC(C)C)C)C"
        ),
        "source": "Horizon-LBA Ono et al. Chem. Sci. 2023; LPE Naylor et al. J. Med. Chem. 2018",
        "pampa": None,
        "permeable": True,
        "hbd": 4,
        "horizon_lba_papp": 3e-6,  # cm/s; similar to CsA per Ono 2023
    },
    {
        "name": "c*[PSLYF]",
        "short": "PSLYF",
        "cycpeptmpdb_id": 1829,
        "smiles": (
            "CC(C)C[C@@H]1NC(=O)[C@H](CO)NC(=O)[C@@H]2CCCN2[C@H](C(=O)NC(C)(C)C)"
            "[C@H](C)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](Cc2ccc(O)cc2)NC1=O"
        ),
        "source": "Hickey, J Med Chem 2016",
        "pampa": -9.10,
        "permeable": False,
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
        "source": "CHUGAI 2013",
        "pampa": -5.20,
        "permeable": True,
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
        "source": "CHUGAI 2013",
        "pampa": -7.00,
        "permeable": False,
        "hbd": None,
    },
    {
        "name": "White_compd3",
        "short": "WhC3",
        "cycpeptmpdb_id": 25,
        "smiles": (
            "CC(C)C[C@@H]1NC(=O)[C@H](Cc2ccc(O)cc2)N(C)C(=O)[C@H]2CCCN2C(=O)"
            "[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)N(C)C(=O)[C@@H](CC(C)C)N(C)C1=O"
        ),
        "source": "White, Nat Chem Biol 2011",
        "pampa": -5.31,  # RRCK assay; no PAMPA measurement available
        "permeable": True,
        "hbd": 3,  # 2 amide NH + Tyr phenol OH
    },
]

SOLVENT_AQ  = "water"
SOLVENT_MEM = "chcl3"
LABEL_AQ    = "water"
LABEL_MEM   = "mem"


def _log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


# ── CREMP Step 1: RDKit conformer embedding ───────────────────────────────────
def embed_rdkit_conformers(mol, n_embed: int = 5000, n_max: int = 50,
                           rmsd_thresh: float = 0.5,
                           mmff_fail_fraction_threshold: float = 0.25):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    params = AllChem.ETKDGv3()
    params.maxIterations = 10 * n_embed
    params.pruneRmsThresh = 0.01
    params.useRandomCoords = True
    params.useMacrocycleTorsions = True
    params.numThreads = 0

    _log(f"    ETKDGv3: embedding {n_embed} conformers...")
    t0 = time.time()
    AllChem.EmbedMultipleConfs(mol, numConfs=n_embed, params=params)
    if mol.GetNumConformers() == 0:
        raise RuntimeError("ETKDGv3 produced zero conformers")
    _log(f"    ETKDGv3: {mol.GetNumConformers()} embedded in {time.time()-t0:.0f}s")

    AllChem.MMFFSanitizeMolecule(mol)
    has_cistrans = any(
        bond.GetStereo() is not Chem.BondStereo.STEREONONE for bond in mol.GetBonds()
    )
    if has_cistrans:
        _log("    MMFF: skipped (E/Z stereo bonds present)")
    else:
        _log("    MMFF: optimising all conformers...")
        t0 = time.time()
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)
        _log(f"    MMFF: done in {time.time()-t0:.0f}s")

    mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
    confs = list(mol.GetConformers())
    energies = []
    mmff_fail_count = 0
    mmff_failed_conf_ids = []
    for conf in confs:
        ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=conf.GetId())
        if ff is None:
            mmff_fail_count += 1
            mmff_failed_conf_ids.append(conf.GetId())
            energies.append(float("inf"))
            continue
        energies.append(ff.CalcEnergy())

    if mmff_fail_count:
        fail_frac = mmff_fail_count / max(1, len(confs))
        _log(f"    MMFF: force-field setup failed for {mmff_fail_count}/{len(confs)} conformers")
        if fail_frac > mmff_fail_fraction_threshold:
            raise RuntimeError(
                f"MMFF force-field setup failed for too many conformers: "
                f"{mmff_fail_count}/{len(confs)} (>{mmff_fail_fraction_threshold:.0%})"
            )

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
    for i in keep:
        new_mol.AddConformer(mol.GetConformer(confs[i].GetId()), assignId=True)
    return new_mol


# ── CREMP Step 2: xTB geometry optimisation ───────────────────────────────────
def _write_conformer_xyz(mol, conf_id: int, xyz_path: Path,
                         comment: str = "") -> None:
    conf = mol.GetConformer(conf_id)
    pos  = conf.GetPositions()
    syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(mol.GetNumAtoms())]
    lines = [str(mol.GetNumAtoms()), comment]
    for sym, (x, y, z) in zip(syms, pos):
        lines.append(f"{sym}  {x:.6f}  {y:.6f}  {z:.6f}")
    xyz_path.write_text("\n".join(lines) + "\n")


def _xtb_opt_worker(args):
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

    env = {**os.environ, "OMP_NUM_THREADS": "1,1"}
    xtb_out = conf_dir / "xtb.out"
    xtb_err = conf_dir / "xtb.err"
    try:
        proc = subprocess.run(
            cmd, cwd=conf_dir, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False,
        )
        xtb_out.write_text(proc.stdout or "")
        xtb_err.write_text(proc.stderr or "")

    if proc.returncode != 0:
        return None

    opt_xyz = conf_dir / "xtbopt.xyz"
    if not opt_xyz.exists():
        return None

    try:
        with open(opt_xyz) as f:
            n_atoms = int(next(f))
            comment = next(f).strip().split()
            energy = None
            for i, tok in enumerate(comment[:-1]):
                if tok.lower().startswith("energy"):
                    try:
                        energy = float(comment[i + 1])
                        break
                    except (ValueError, IndexError):
                        pass
            if energy is None:
                return None
            new_conf = Chem.Conformer(n_atoms)
            for i in range(n_atoms):
                parts = next(f).split()
                new_conf.SetAtomPosition(i, [float(parts[1]),
                                             float(parts[2]),
                                             float(parts[3])])
        mol_copy = Chem.Mol(mol)
        mol_copy.RemoveAllConformers()
        mol_copy.AddConformer(new_conf, assignId=True)
        return (conf_id, mol_copy, energy)
    except Exception:
        return None


def xtb_preopt_mol(mol, solvent: str, work_dir: Path,
                   charge: int, n_procs: int = 1):
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
            _log(f"    xTB: {done}/{len(conf_ids)} done ({elapsed:.0f}s) — conf {result[0] if result else '?'} {status}")
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
    work_dir.mkdir(parents=True, exist_ok=True)
    xyz_path = xyz_path.resolve()

    cmd = [
        "crest", str(xyz_path),
        "-T",      str(n_threads),
        "--gfn2",
        "--chrg",  str(charge),
        "--alpb",  solvent,
        "--keepdir",
    ]

    _log(f"    CREST: {' '.join(cmd)}")
    t0 = time.time()
    out_path = work_dir / "crest.out"
    err_path = work_dir / "crest.err"
    with open(out_path, "w") as fout, open(err_path, "w") as ferr:
        proc = subprocess.run(
            cmd, cwd=work_dir.resolve(),
            stdout=fout, stderr=ferr,
        )

    elapsed = time.time() - t0
    ensemble = work_dir.resolve() / "crest_conformers.xyz"
    if proc.returncode != 0:
        _log(f"    CREST: FAILED after {elapsed:.0f}s (exit={proc.returncode})")
        try:
            tail = err_path.read_text().splitlines()[-20:]
            print("\n".join(f"      ERR: {l}" for l in tail), flush=True)
        except Exception:
            pass
        return None
    if ensemble.exists() and ensemble.stat().st_size > 0:
        _log(f"    CREST: finished in {elapsed/3600:.2f}h — ensemble found")
        return ensemble
    _log(f"    CREST: finished with exit=0 but no ensemble after {elapsed:.0f}s")
    return None


# ── CREST cregen: refinement from existing rotamers, skipping MTDs ────────────
def run_crest_cregen(xyz_path: Path, work_dir: Path, solvent: str,
                     n_threads: int, charge: int = 0) -> Path | None:
    """Run CREST --cregen on existing crest_rotamers_0.xyz, skipping the MTD phase."""
    rotamers = work_dir.resolve() / "crest_rotamers_0.xyz"
    if not rotamers.exists() or rotamers.stat().st_size == 0:
        _log(f"    CREST --cregen: crest_rotamers_0.xyz not found in {work_dir}")
        return None

    cmd = [
        "crest", str(xyz_path.resolve()),
        "--cregen",
        "-T",     str(n_threads),
        "--gfn2",
        "--chrg", str(charge),
        "--alpb", solvent,
    ]

    _log(f"    CREST --cregen: {' '.join(cmd)}")
    _log(f"    CREST --cregen: refining {rotamers.stat().st_size/1e6:.0f} MB rotamers file")
    t0 = time.time()
    out_path = work_dir / "crest_cregen.out"
    err_path = work_dir / "crest_cregen.err"
    with open(out_path, "w") as fout, open(err_path, "w") as ferr:
        proc = subprocess.run(
            cmd, cwd=work_dir.resolve(),
            stdout=fout, stderr=ferr,
        )

    elapsed = time.time() - t0
    ensemble = work_dir.resolve() / "crest_conformers.xyz"
    if proc.returncode != 0:
        _log(f"    CREST --cregen: FAILED after {elapsed:.0f}s (exit={proc.returncode})")
        try:
            tail = err_path.read_text().splitlines()[-20:]
            print("\n".join(f"      ERR: {l}" for l in tail), flush=True)
        except Exception:
            pass
        return None
    if ensemble.exists() and ensemble.stat().st_size > 0:
        _log(f"    CREST --cregen: finished in {elapsed/3600:.2f}h — ensemble found")
        return ensemble
    _log(f"    CREST --cregen: finished with exit=0 but no ensemble after {elapsed:.0f}s")
    return None


# ── Parse crest.out → CREMP-format thermodynamics ─────────────────────────────
def parse_crest_log(crest_out_path: Path) -> dict | None:
    """
    Parse crest.out for ensemble thermodynamics and per-conformer data.
    Adapted from CREMP cremp/utils/postprocess.py parse_ensemble_data().
    """
    ensemble_patterns = {
        "ensembleenergy":     ("ensemble average energy", float),
        "ensembleentropy":    ("ensemble entropy",        float),
        "ensemblefreeenergy": ("ensemble free energy",   float),
        "lowestenergy":       ("E lowest",               float),
        "poplowestpct":       ("population of lowest in %", float),
        "temperature":        ("T /K",                   float),
        "uniqueconfs":        ("number of unique conformers", int),
        "totalconfs":         ("total number unique points", int),
    }
    conf_cols = {
        "relativeenergy":  1,
        "totalenergy":     2,
        "conformerweight": 3,
        "boltzmannweight": 4,
        "set":             5,
        "degeneracy":      6,
    }

    ensemble_data: dict = {}
    conf_data: list[dict] = []
    read = False

    try:
        lines = crest_out_path.read_text(errors="replace").splitlines()
    except Exception:
        return None

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if "Final Geometry Optimization" in line:
            read = True

        if read:
            for key, (pattern, typ) in ensemble_patterns.items():
                if line.startswith(pattern):
                    try:
                        ensemble_data[key] = typ(line.split()[-1])
                    except (ValueError, IndexError):
                        pass

            if line.startswith("Erel/kcal"):
                i += 1
                while i < len(lines):
                    conf_line = lines[i].strip()
                    parts = conf_line.split()
                    if len(parts) < 5:
                        break
                    if len(parts) >= 7:
                        try:
                            single: dict = {
                                "relativeenergy":   float(parts[conf_cols["relativeenergy"]]),
                                "totalenergy":      float(parts[conf_cols["totalenergy"]]),
                                "conformerweights": [float(parts[conf_cols["conformerweight"]])],
                                "boltzmannweight":  float(parts[conf_cols["boltzmannweight"]]),
                                "set":              int(parts[conf_cols["set"]]),
                                "degeneracy":       int(parts[conf_cols["degeneracy"]]),
                            }
                            conf_data.append(single)
                        except (ValueError, IndexError):
                            pass
                    elif conf_data:
                        try:
                            conf_data[-1]["conformerweights"].append(float(parts[3]))
                        except (ValueError, IndexError):
                            pass
                    i += 1
                continue
        i += 1

    if not ensemble_data:
        return None

    ensemble_data["conformers"] = conf_data
    return ensemble_data


# ── Parse multi-conformer XYZ ─────────────────────────────────────────────────
def parse_xyz_ensemble(xyz_path: Path) -> tuple[list, int]:
    conformers = []
    failed_energy_parses = 0
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
        energy = np.nan
        if i < len(lines):
            comment_line = lines[i].strip()
            parsed = False
            for tok in comment_line.replace('=', ' ').replace(':', ' ').split():
                try:
                    energy = float(tok)
                    parsed = True
                    break
                except ValueError:
                    continue
            if not parsed:
                failed_energy_parses += 1
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
    return conformers, failed_energy_parses


def boltzmann_weights(energies_hartree: list[float], T: float = 298.15) -> np.ndarray:
    KCAL_PER_HARTREE = 627.509
    RT = 1.987e-3 * T
    e = np.array(energies_hartree, dtype=float) * KCAL_PER_HARTREE
    valid = np.isfinite(e)
    if not np.any(valid):
        raise RuntimeError("No valid conformer energies for Boltzmann weighting")
    e_rel = e[valid] - np.nanmin(e[valid])
    valid_weights = np.exp(-e_rel / RT)
    valid_weights /= valid_weights.sum()
    weights = np.full(len(e), np.nan, dtype=float)
    weights[valid] = valid_weights
    return weights


# ── 3D PSA ────────────────────────────────────────────────────────────────────
def compute_psa_xyz(symbols: list[str], coords: np.ndarray,
                    template_mol=None) -> float:
    from rdkit import Chem
    from rdkit.Chem import rdFreeSASA

    _BONDI = {
        'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
        'S': 1.80, 'P': 1.80, 'F': 1.47, 'Cl': 1.75,
        'Br': 1.85, 'I': 1.98,
    }
    _POLAR = {'N', 'O', 'S', 'P'}

    if template_mol is not None:
        try:
            mol_h = Chem.RWMol(template_mol)
            if mol_h.GetNumAtoms() != len(symbols):
                raise ValueError("Atom count mismatch")
            conf = Chem.Conformer(mol_h.GetNumAtoms())
            for i, (x, y, z) in enumerate(coords):
                conf.SetAtomPosition(i, (float(x), float(y), float(z)))
            conf_id = mol_h.AddConformer(conf, assignId=True)
            mol_h = mol_h.GetMol()
            radii = []
            for atom in mol_h.GetAtoms():
                sym = atom.GetSymbol()
                radii.append(_BONDI.get(sym, 1.50))
                if sym in _POLAR:
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

    PROBE = 1.40
    n = len(symbols)
    radii_all = np.array([_BONDI.get(s, 1.50) + PROBE for s in symbols])
    polar_idx = [i for i, s in enumerate(symbols) if s in _POLAR]
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


# ── Intramolecular H-bonds ────────────────────────────────────────────────────
def count_hbonds_xyz(symbols: list[str], coords: np.ndarray) -> int:
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


# ── SDF export ────────────────────────────────────────────────────────────────
def export_sdf(conformers: list, smiles: str, out_path: Path,
               template_mol=None) -> bool:
    """Write all conformers to SDF using RDKit connectivity from template mol."""
    from rdkit import Chem

    if template_mol is None:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        if mol is None:
            return False
    else:
        mol = Chem.RWMol(template_mol)
        mol.RemoveAllConformers()

    for syms, crds, _ in conformers:
        if mol.GetNumAtoms() != len(syms):
            continue
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, (x, y, z) in enumerate(crds):
            conf.SetAtomPosition(i, (float(x), float(y), float(z)))
        mol.AddConformer(conf, assignId=True)

    if mol.GetNumConformers() == 0:
        return False

    writer = Chem.SDWriter(str(out_path))
    for cid in range(mol.GetNumConformers()):
        writer.write(mol, confId=cid)
    writer.close()
    return True


# ── JSON export ───────────────────────────────────────────────────────────────
def export_json(conformers: list, psa_vals: list, hb_vals: list,
                weights: np.ndarray, smiles: str, charge: int,
                crest_log: dict | None, out_path: Path) -> None:
    """
    Write CREMP-format JSON extended with per-conformer PSA and HB values.
    If crest_log is available, ensemble thermodynamics are included.
    """
    conf_records = []
    for i, (_, _, energy) in enumerate(conformers):
        rec: dict = {
            "totalenergy":    float(energy) if np.isfinite(energy) else None,
            "boltzmannweight": float(weights[i]) if np.isfinite(weights[i]) else None,
            "psa":            psa_vals[i],
            "hbonds":         int(hb_vals[i]),
        }
        # Merge per-conformer data from crest.out if available and indices align
        if crest_log and i < len(crest_log.get("conformers", [])):
            log_conf = crest_log["conformers"][i]
            rec["relativeenergy"]   = log_conf.get("relativeenergy")
            rec["conformerweights"] = log_conf.get("conformerweights")
            rec["set"]              = log_conf.get("set")
            rec["degeneracy"]       = log_conf.get("degeneracy")
        conf_records.append(rec)

    data: dict = {
        "smiles":    smiles,
        "charge":    charge,
        "n_confs":   len(conformers),
    }

    if crest_log:
        for key in ("ensembleenergy", "ensembleentropy", "ensemblefreeenergy",
                    "lowestenergy", "poplowestpct", "temperature",
                    "uniqueconfs", "totalconfs"):
            if key in crest_log:
                data[key] = crest_log[key]

    valid_w = np.isfinite(weights)
    psa_arr = np.array(psa_vals)
    hb_arr  = np.array(hb_vals, dtype=float)
    data["boltzmann_psa"] = round(float(np.dot(weights[valid_w], psa_arr[valid_w])), 2)
    data["boltzmann_hb"]  = round(float(np.dot(weights[valid_w], hb_arr[valid_w])), 2)
    data["lowen_psa"]     = psa_vals[0]
    data["lowen_hb"]      = int(hb_vals[0])
    data["conformers"]    = conf_records

    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)


# ── Process one compound ──────────────────────────────────────────────────────
def process_compound(cpd: dict, work_base: Path,
                     max_confs: int | None, dry_run: bool,
                     n_threads: int) -> dict:
    from rdkit import Chem

    name  = cpd["name"]
    short = cpd["short"]
    smi   = cpd["smiles"]

    print(f"\n{'─'*60}")
    print(f"  {name}  (ID={cpd['cycpeptmpdb_id']})  PAMPA={cpd['pampa']}")

    result = {
        "compound":       name,
        "cycpeptmpdb_id": cpd["cycpeptmpdb_id"],
        "pampa":          cpd["pampa"],
        "permeable":      cpd["permeable"],
        "run_id":         work_base.name,
    }

    if smi is None:
        print(f"  ⚠ No SMILES — skipping")
        return result

    work_base.mkdir(parents=True, exist_ok=True)
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"  ⚠ Invalid SMILES")
        return result
    mol = Chem.AddHs(mol)
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    _log(f"  Formal charge: {charge}  |  Atoms (with H): {mol.GetNumAtoms()}")

    _log(f"  Step 1: RDKit ETKDGv3 embedding (5000 → 50)...")
    try:
        mol_embedded = embed_rdkit_conformers(mol)
    except RuntimeError as e:
        _log(f"  ⚠ Embedding failed: {e}")
        return result
    _log(f"  Step 1 done: {mol_embedded.GetNumConformers()} conformers")
    _template_mol = mol_embedded

    failed_solvents = []

    for solvent, label in [(SOLVENT_AQ, LABEL_AQ), (SOLVENT_MEM, LABEL_MEM)]:
        sol_dir = work_base / label
        sol_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n  ── [{label}] solvent={solvent} ──")

        if dry_run:
            rng = np.random.default_rng(42 if label == LABEL_AQ else 43)
            psa_mean = 180.0 if label == LABEL_AQ else 140.0
            psa_vals = rng.normal(psa_mean, 10, 50)
            hb_vals  = rng.integers(0, 4, 50)
            result[f"{label}_psa_boltz"] = float(psa_vals.mean())
            result[f"{label}_hb_boltz"]  = float(hb_vals.mean())
            result[f"{label}_n_confs"]   = 50
            result[f"{label}_status"]    = "dry_run"
            continue

        # ── Checkpoint: skip xTB + CREST if ensemble already exists ─────────────
        crest_dir = sol_dir / "crest"
        xyz_in    = crest_dir / f"{short}_{label}_start.xyz"
        existing_ensemble = crest_dir / "crest_conformers.xyz"
        existing_rotamers = crest_dir / "crest_rotamers_0.xyz"

        if existing_ensemble.exists() and existing_ensemble.stat().st_size > 0:
            _log(f"  Steps 2-3: existing CREST ensemble found "
                 f"({existing_ensemble.stat().st_size/1e6:.0f} MB) — skipping xTB + CREST")
            raw_ensemble = existing_ensemble
        elif existing_rotamers.exists() and existing_rotamers.stat().st_size > 0 and xyz_in.exists():
            # MTDs finished and start xyz exists — skip xTB entirely, run cregen only
            _log(f"  Steps 2-3: rotamers found ({existing_rotamers.stat().st_size/1e6:.0f} MB)"
                 f" + start xyz exists — skipping xTB, running --cregen")
            try:
                raw_ensemble = run_crest_cregen(xyz_in, crest_dir, solvent,
                                                n_threads, charge)
            except Exception as e:
                print(f"      ⚠ CREST --cregen error: {e}")
                raw_ensemble = None
        else:
            # Step 2: xTB pre-opt
            xtb_dir = sol_dir / "xtb_opt"
            _log(f"  Step 2: xTB pre-opt ({mol_embedded.GetNumConformers()} conformers)...")
            mol_min = xtb_preopt_mol(mol_embedded, solvent, xtb_dir, charge,
                                     n_procs=n_threads)
            if mol_min is None:
                _log(f"  ⚠ All xTB jobs failed — skipping {label}")
                result[f"{label}_status"] = "failed"
                result[f"{label}_error"]  = "all xTB jobs failed"
                failed_solvents.append(label)
                continue

            crest_dir.mkdir(parents=True, exist_ok=True)

            # Step 3: CREST — cregen from existing rotamers (no start xyz yet), or full iMTD-GC
            if existing_rotamers.exists() and existing_rotamers.stat().st_size > 0:
                _log(f"  Step 3: rotamers found ({existing_rotamers.stat().st_size/1e6:.0f} MB)"
                     f" — running --cregen to skip MTDs")
                _write_conformer_xyz(mol_min, mol_min.GetConformer().GetId(),
                                     xyz_in, comment=smi)
                try:
                    raw_ensemble = run_crest_cregen(xyz_in, crest_dir, solvent,
                                                    n_threads, charge)
                except Exception as e:
                    print(f"      ⚠ CREST --cregen error: {e}")
                    raw_ensemble = None
            else:
                _write_conformer_xyz(mol_min, mol_min.GetConformer().GetId(),
                                     xyz_in, comment=smi)
                _log(f"  Step 3: CREST iMTD-GC ({solvent})...")
                try:
                    raw_ensemble = run_crest(xyz_in, crest_dir, solvent,
                                             n_threads, charge)
                except Exception as e:
                    print(f"      ⚠ CREST error: {e}")
                    raw_ensemble = None

        if raw_ensemble is None:
            result[f"{label}_status"] = "failed"
            result[f"{label}_error"]  = "no ensemble produced"
            failed_solvents.append(label)
            continue

        # Copy ensemble.xyz to sol_dir level for easy access
        ensemble_xyz = sol_dir / "ensemble.xyz"
        shutil.copy2(raw_ensemble, ensemble_xyz)
        _log(f"  Ensemble saved → {ensemble_xyz.name}")

        # ── Parse ensemble ────────────────────────────────────────────────────
        conformers, energy_fail = parse_xyz_ensemble(ensemble_xyz)
        n_confs_full = len(conformers)
        print(f"      Parsed {n_confs_full} conformers")
        if energy_fail:
            _log(f"      Warning: {energy_fail}/{n_confs_full} conformers missing energy")

        if n_confs_full == 0:
            result[f"{label}_status"] = "failed"
            result[f"{label}_error"]  = "empty ensemble"
            failed_solvents.append(label)
            continue

        if max_confs is not None:
            conformers = conformers[:max_confs]
        n_confs = len(conformers)

        psa_vals, hb_vals, energies = [], [], []
        for syms, crds, eng in conformers:
            psa_vals.append(compute_psa_xyz(syms, crds, template_mol=_template_mol))
            hb_vals.append(count_hbonds_xyz(syms, crds))
            energies.append(eng)

        psa_arr = np.array(psa_vals)
        hb_arr  = np.array(hb_vals, dtype=float)
        try:
            weights = boltzmann_weights(energies)
        except RuntimeError as e:
            _log(f"  ⚠ Boltzmann weighting failed ({e}) — skipping {label}")
            result[f"{label}_status"] = "failed"
            result[f"{label}_error"]  = str(e)
            failed_solvents.append(label)
            continue
        valid_w = np.isfinite(weights)

        psa_boltz = float(np.dot(weights[valid_w], psa_arr[valid_w]))
        hb_boltz  = float(np.dot(weights[valid_w], hb_arr[valid_w]))

        result[f"{label}_psa_lowen"]    = round(float(psa_arr[0]), 2)
        result[f"{label}_psa_boltz"]    = round(psa_boltz, 2)
        result[f"{label}_psa_min"]      = round(float(psa_arr.min()), 2)
        result[f"{label}_psa_max"]      = round(float(psa_arr.max()), 2)
        result[f"{label}_psa_mean"]     = round(float(psa_arr.mean()), 2)
        result[f"{label}_psa_std"]      = round(float(psa_arr.std()), 2)
        result[f"{label}_hb_lowen"]     = int(hb_arr[0])
        result[f"{label}_hb_boltz"]     = round(hb_boltz, 2)
        result[f"{label}_hb_min"]       = int(hb_arr.min())
        result[f"{label}_hb_max"]       = int(hb_arr.max())
        result[f"{label}_hb_mean"]      = round(float(hb_arr.mean()), 2)
        result[f"{label}_n_confs"]      = n_confs
        result[f"{label}_n_confs_full"] = n_confs_full
        result[f"{label}_status"]       = "ok"

        print(f"      PSA (Boltzmann): {psa_boltz:.1f} Å²  "
              f"(low-energy={psa_arr[0]:.1f}, mean={psa_arr.mean():.1f})")
        print(f"      HB  (Boltzmann): {hb_boltz:.2f}  (low-energy={int(hb_arr[0])})")

        # ── Parse CREST log for thermodynamics ───────────────────────────────
        # Prefer crest_cregen.out (written by --cregen runs) over crest.out
        crest_log = None
        for log_candidate in [sol_dir / "crest" / "crest_cregen.out",
                               sol_dir / "crest" / "crest.out"]:
            if log_candidate.exists():
                crest_log = parse_crest_log(log_candidate)
                if crest_log:
                    _log(f"      {log_candidate.name} parsed: "
                         f"{crest_log.get('uniqueconfs','?')} unique confs")
                    break
        if not crest_log:
            _log(f"      No parseable CREST log found — JSON will omit thermodynamics")

        # ── Export SDF ────────────────────────────────────────────────────────
        sdf_path = sol_dir / "ensemble.sdf"
        if export_sdf(conformers, smi, sdf_path, template_mol=_template_mol):
            _log(f"      SDF written → {sdf_path.name} ({n_confs} conformers)")
        else:
            _log(f"      ⚠ SDF export failed")

        # ── Export JSON ───────────────────────────────────────────────────────
        json_path = sol_dir / "ensemble.json"
        export_json(conformers, psa_vals, hb_vals, weights,
                    smi, charge, crest_log, json_path)
        _log(f"      JSON written → {json_path.name}")

        # ── Checkpoint: save partial CSV after each solvent ───────────────────
        partial_csv = work_base / f"{short}_results_partial.csv"
        pd.DataFrame([result]).to_csv(partial_csv, index=False)
        _log(f"      Checkpoint saved → {partial_csv.name}")

    if failed_solvents:
        _log(f"  WARNING: {short} failed for solvent(s): {', '.join(sorted(set(failed_solvents)))}")
        result["failed_solvents"] = ",".join(sorted(set(failed_solvents)))

    # ── Δ features ────────────────────────────────────────────────────────────
    aq_key  = f"{LABEL_AQ}_psa_boltz"
    mem_key = f"{LABEL_MEM}_psa_boltz"
    if aq_key in result and mem_key in result:
        result["crest_delta_psa"]      = round(result[aq_key] - result[mem_key], 2)
        aq_lo  = result.get(f"{LABEL_AQ}_psa_lowen",  np.nan)
        mem_lo = result.get(f"{LABEL_MEM}_psa_lowen", np.nan)
        result["crest_delta_psa_lowen"] = round(aq_lo - mem_lo, 2) if pd.notna(aq_lo) and pd.notna(mem_lo) else np.nan
        result["crest_delta_hb"]        = round(result[f"{LABEL_MEM}_hb_boltz"] - result[f"{LABEL_AQ}_hb_boltz"], 2)
        result["crest_psa_spread_aq"]   = result.get(f"{LABEL_AQ}_psa_std",  np.nan)
        result["crest_psa_spread_mem"]  = result.get(f"{LABEL_MEM}_psa_std", np.nan)
        print(f"\n  ✓ CREST ΔPSA (Boltzmann) = {result['crest_delta_psa']:.1f} Å²  "
              f"(low-energy = {result['crest_delta_psa_lowen']:.1f} Å²  "
              f"ΔHB = {result['crest_delta_hb']:.2f})")

    return result


# ── Auto-resume: find most recent partial run for this compound ───────────────
def _find_resume_dir(runs_base: Path, compound_idx: int, short: str) -> Path | None:
    """Return the most recent run dir that has checkpoint data (rotamers or ensemble)."""
    if not runs_base.exists():
        return None
    checkpoint_paths = [
        f"{label}/crest/crest_conformers.xyz"
        for label in ("water", "mem")
    ] + [
        f"{label}/crest/crest_rotamers_0.xyz"
        for label in ("water", "mem")
    ]
    candidates = sorted(
        runs_base.glob(f"run_*_{compound_idx}_{short}"),
        reverse=True,  # most recent timestamp first
    )
    for d in candidates:
        if any((d / p).exists() for p in checkpoint_paths):
            return d
    return None


# ── Main ──────────────────────────────────────────────────────────────────────
def run(outdir: Path, max_confs: int | None, dry_run: bool,
        compound_idx: int | None = None, n_threads: int | None = None) -> None:

    n_threads = n_threads or os.cpu_count() or 1

    if compound_idx is None:
        raise RuntimeError("--compound is required.")
    if compound_idx < 0 or compound_idx >= len(REFERENCE_COMPOUNDS):
        raise ValueError(f"--compound must be 0–{len(REFERENCE_COMPOUNDS)-1}")

    cpd   = REFERENCE_COMPOUNDS[compound_idx]
    short = cpd["short"]

    resume_dir = _find_resume_dir(outdir / "runs", compound_idx, short)
    if resume_dir is not None:
        work_base = resume_dir
        print(f"\n[Compound {compound_idx}] {cpd['name']}  ← resuming")
        print(f"Run directory: {work_base}")
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_id    = f"run_{timestamp}_{compound_idx}_{short}"
        work_base = outdir / "runs" / run_id
        work_base.mkdir(parents=True, exist_ok=True)
        print(f"\n[Compound {compound_idx}] {cpd['name']}")
        print(f"Run directory: {work_base}")

    r = process_compound(cpd, work_base, max_confs, dry_run, n_threads)

    out_csv = work_base / f"{short}_results.csv"
    pd.DataFrame([r]).to_csv(out_csv, index=False)
    print(f"\nSaved: {out_csv}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Tier-2 CREST+ALPB validation (v3.2 — CREMP-format outputs)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Compound index:
  0 = HexPep  (impermeable,  6-mer)
  1 = CsA     (permeable,   11-mer)
  2 = PSLYF   (impermeable, 11-mer)
  3 = DP-955  (permeable,   15-mer)
  4 = DP-944  (impermeable, 15-mer)
  5 = WhC3    (permeable,    6-mer, White 2011 compd.3, RRCK=-5.31)
        """
    )
    parser.add_argument("--outdir",    "-o", default="results", type=Path)
    parser.add_argument("--compound",  type=int, default=None, metavar="IDX")
    parser.add_argument("--threads",   type=int, default=None, metavar="N")
    parser.add_argument("--max-confs", "-c", type=int, default=None,
                        help="Cap conformers used for PSA analysis. Default: all.")
    parser.add_argument("--dry-run",   action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(
        Path(args.outdir),
        max_confs    = args.max_confs,
        dry_run      = args.dry_run,
        compound_idx = args.compound,
        n_threads    = args.threads,
    )

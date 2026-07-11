# env: chameleon-sim
"""
crest_conformers.py
-------------------
CREST iMTD-GC conformer sampling pipeline.
Handles all simulation steps (RDKit embedding, xTB pre-opt, CREST, post-processing)
and ensemble I/O (XYZ parsing, SDF/JSON export).

Imported by crest_v3.2.py — not run directly.
"""

from __future__ import annotations

import json
import multiprocessing
import os
import re
import shutil
import subprocess
import time
import warnings
from datetime import datetime
from pathlib import Path

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import pandas as pd

from phys_descriptors_v2 import boltzmann_weights, compute_psa_xyz, count_hbonds_xyz

warnings.filterwarnings("ignore")

SOLVENT_AQ  = "water"
SOLVENT_MEM = "chcl3"
LABEL_AQ    = "water"
LABEL_MEM   = "mem"


def _log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


# ── RDKit conformer embedding ─────────────────────────────────────────────────
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


# ── xTB geometry optimisation ─────────────────────────────────────────────────
def _write_conformer_xyz(mol, conf_id: int, xyz_path: Path,
                         comment: str = "") -> None:
    conf = mol.GetConformer(conf_id)
    pos  = conf.GetPositions()
    syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(mol.GetNumAtoms())]
    lines = [str(mol.GetNumAtoms()), comment]
    for sym, (x, y, z) in zip(syms, pos):
        lines.append(f"{sym}  {x:.6f}  {y:.6f}  {z:.6f}")
    xyz_path.write_text("\n".join(lines) + "\n")


# ── diazirine N=N constraint (GFN2 stretches it to a spurious ~1.43 Å; see ──────
#    docs/experiments + memory diazirine-review-checklist). True N=N = 1.228 Å (exp,
#    microwave) / 1.230 Å (CCSD(T)). Constrain it so the macrocycle samples freely
#    while the rigid diazirine can't drift into GFN2's single-bond basin.
DIAZIRINE_SMARTS = "[#6]1[#7]=[#7]1"
DIAZIRINE_NN = 1.23   # Å; constrain target (literature 1.228–1.230)


def diazirine_nn_atoms(mol):
    """1-based (N1, N2) indices of the diazirine N=N for xtb/CREST, or None if absent."""
    from rdkit import Chem
    matches = mol.GetSubstructMatches(Chem.MolFromSmarts(DIAZIRINE_SMARTS))
    if not matches:
        return None
    _c, n1, n2 = matches[0]      # SMARTS order: carbon, then the two nitrogens
    return n1 + 1, n2 + 1        # xtb/CREST atom indices are 1-based


def write_constraint_file(path: Path, nn_pair, value: float = DIAZIRINE_NN,
                          fc: float = 0.25) -> Path:
    """Write an xtb/CREST $constrain file pinning the diazirine N=N distance (Fix 1).
    Used for both the xTB pre-opt (--input) and the CREST search (--cinp).

    force constant = 0.25 Eh/Bohr^2 — the CREST-documented value for DISTANCE constraints
    (https://crest-lab.github.io/crest-docs/page/examples/example_4.html#constrained-sampling;
    the docs caution against high values). At r0=1.23 Å this is a ~20 kcal/mol barrier vs the
    1.43 Å drift. Escalate to 0.5 → 1.0 only if the integrity check shows WATCH/FAIL drift.
    No $metadyn block / reference file is needed for a distance constraint (those are only for
    substructure fixing = Fix 2)."""
    n1, n2 = nn_pair
    path.write_text(
        "$constrain\n"
        f"  force constant={fc}\n"
        f"  distance: {n1}, {n2}, {value}\n"
        "$end\n"
    )
    return path


def _xtb_opt_worker(args):
    from rdkit import Chem

    conf_id, mol, conf_dir, solvent, charge, constraint_file = args
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
    if constraint_file:
        cmd.extend(["--input", str(constraint_file)])   # diazirine N=N constraint

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
    except Exception:
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
                   charge: int, n_procs: int = 1, constraint_file: Path | None = None):
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    conf_ids = [c.GetId() for c in mol.GetConformers()]
    if not conf_ids:
        return None

    params_list = [
        (conf_id, mol, work_dir / f"conf_{conf_id}", solvent, charge, constraint_file)
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


# ── CREST conformer sampling ──────────────────────────────────────────────────
def run_crest(xyz_path: Path, work_dir: Path, solvent: str,
              n_threads: int, charge: int = 0, constraint_file: Path | None = None) -> Path | None:
    work_dir.mkdir(parents=True, exist_ok=True)
    xyz_path = xyz_path.resolve()

    cmd = [
        "crest", str(xyz_path),
        "-T",      str(n_threads),
        "--gfn2",
        "--chrg",  str(charge),
        "--alpb",  solvent,
        "--keepdir",
        "--noreftopo",
        "-notopo",
    ]
    if constraint_file:
        cmd.extend(["--cinp", str(constraint_file)])    # diazirine N=N constraint

    _log(f"    CREST: {' '.join(cmd)}")
    t0 = time.time()
    out_path = work_dir / "crest.out"
    err_path = work_dir / "crest.err"
    with open(out_path, "w") as fout, open(err_path, "w") as ferr:
        proc = subprocess.run(cmd, cwd=work_dir.resolve(), stdout=fout, stderr=ferr)

    elapsed = time.time() - t0
    ensemble = work_dir.resolve() / "crest_conformers.xyz"
    if ensemble.exists() and ensemble.stat().st_size > 0:
        if proc.returncode != 0:
            _log(f"    CREST: exit={proc.returncode} but ensemble found — treating as success")
        else:
            _log(f"    CREST: finished in {elapsed/3600:.2f}h — ensemble found")
        return ensemble
    if proc.returncode != 0:
        _log(f"    CREST: FAILED after {elapsed:.0f}s (exit={proc.returncode})")
        try:
            tail = err_path.read_text().splitlines()[-20:]
            print("\n".join(f"      ERR: {l}" for l in tail), flush=True)
        except Exception:
            pass
    else:
        _log(f"    CREST: finished with exit=0 but no ensemble after {elapsed:.0f}s")
    return None


def run_crest_cregen(xyz_path: Path, work_dir: Path, solvent: str,
                     n_threads: int, charge: int = 0, constraint_file: Path | None = None) -> Path | None:
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
    if constraint_file:
        cmd.extend(["--cinp", str(constraint_file)])    # diazirine N=N constraint

    _log(f"    CREST --cregen: {' '.join(cmd)}")
    _log(f"    CREST --cregen: refining {rotamers.stat().st_size/1e6:.0f} MB rotamers file")
    t0 = time.time()
    out_path = work_dir / "crest_cregen.out"
    err_path = work_dir / "crest_cregen.err"
    with open(out_path, "w") as fout, open(err_path, "w") as ferr:
        proc = subprocess.run(cmd, cwd=work_dir.resolve(), stdout=fout, stderr=ferr)

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


# ── Parse crest.out ───────────────────────────────────────────────────────────
def parse_crest_log(crest_out_path: Path) -> dict | None:
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


# ── SDF export ────────────────────────────────────────────────────────────────
def export_sdf(conformers: list, smiles: str, out_path: Path,
               template_mol=None) -> bool:
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
    conf_records = []
    for i, (_, _, energy) in enumerate(conformers):
        rec: dict = {
            "totalenergy":    float(energy) if np.isfinite(energy) else None,
            "boltzmannweight": float(weights[i]) if np.isfinite(weights[i]) else None,
            "psa":            psa_vals[i],
            "hbonds":         int(hb_vals[i]),
        }
        if crest_log and i < len(crest_log.get("conformers", [])):
            log_conf = crest_log["conformers"][i]
            rec["relativeenergy"]   = log_conf.get("relativeenergy")
            rec["conformerweights"] = log_conf.get("conformerweights")
            rec["set"]              = log_conf.get("set")
            rec["degeneracy"]       = log_conf.get("degeneracy")
        conf_records.append(rec)

    data: dict = {"smiles": smiles, "charge": charge, "n_confs": len(conformers)}

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
# ── Registry-free, direct-SMILES entry point (notebook front end) ─────────────
def check_binaries(require=("xtb", "crest")) -> dict:
    """Return {binary: resolved_path}. Raise RuntimeError if any required binary is
    not on PATH, with actionable install guidance. CREST/xTB are external programs;
    the conformer search cannot run without them."""
    found = {b: shutil.which(b) for b in require}
    missing = [b for b, p in found.items() if p is None]
    if missing:
        raise RuntimeError(
            f"Required external binary/binaries not found on PATH: {', '.join(missing)}. "
            f"The CREST/xTB conformer search cannot run without them. Install with "
            f"`conda install -c conda-forge xtb crest` (or load your cluster's modules), "
            f"then restart the kernel/shell so PATH is refreshed.")
    return found


def _safe_short(name: str) -> str:
    """Filesystem-safe short label derived from a molecule name (used in filenames)."""
    s = re.sub(r"[^0-9A-Za-z._-]+", "_", (name or "").strip()).strip("_")
    return (s or "molecule")[:40]


def generate_ensembles(smiles: str, name: str = "molecule",
                       outdir: str | Path = "results/notebook_runs",
                       charge: int | None = None, n_threads: int | None = None,
                       max_confs: int | None = None,
                       check_binaries_first: bool = True) -> dict:
    """Generate water + chloroform (mem) CREST/xTB conformer ensembles for an **arbitrary
    SMILES** — the registry-free front end used by the notebook.

    Wraps the existing `process_compound` engine (no logic duplicated); builds a minimal
    compound record from the SMILES so no entry in `crest_v3.2.REFERENCE_COMPOUNDS` is needed.

    Parameters
    ----------
    smiles   : the molecule to sample (validated with RDKit; ValueError if unparseable)
    name     : optional label (used for the run-dir and file names)
    outdir   : where the run directory is created
    charge   : optional formal-charge override (default: auto-derived from the SMILES)
    n_threads: CPU cores for xTB/CREST (default: all cores)
    max_confs: cap on conformers kept for the chloroform ensemble (default: engine default, 50)
    check_binaries_first : if True (default), fail fast with a clear message when xtb/crest
                           are missing, before doing any work.

    Returns
    -------
    dict with: work_dir, result (the raw engine row), water/mem ensemble paths, and `ok`
    (True iff all four ensemble files were written).
    """
    from rdkit import Chem

    # 1. validate the SMILES up front — clear error rather than a deep failure later
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("`smiles` must be a non-empty string.")
    if Chem.MolFromSmiles(smiles) is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}. Check the string and try again.")

    # 2. validate external binaries early (CREST/xTB are required)
    if check_binaries_first:
        check_binaries()

    # 3. build a registry-free compound record (metadata fields are optional → None)
    short = _safe_short(name)
    cpd = {"name": name, "short": short, "smiles": smiles,
           "cycpeptmpdb_id": None, "pampa": None, "permeable": None}

    # 4. fresh, timestamped run directory
    n_threads = n_threads or os.cpu_count() or 1
    work_base = Path(outdir) / f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{short}"
    work_base.mkdir(parents=True, exist_ok=True)
    _log(f"Direct-SMILES run: {name}  →  {work_base}")

    # 5. run the existing engine (with optional charge override)
    result = process_compound(cpd, work_base, max_confs=max_confs, dry_run=False,
                              n_threads=n_threads, charge_override=charge)

    # 6. collect + verify outputs
    paths = {sol: {"sdf": work_base / sol / "ensemble.sdf",
                   "json": work_base / sol / "ensemble.json"}
             for sol in ("water", "mem")}
    ok = all(p.exists() for sol in paths for p in paths[sol].values())
    return {"work_dir": work_base, "result": result,
            "water": paths["water"], "mem": paths["mem"], "ok": ok}


def process_compound(cpd: dict, work_base: Path,
                     max_confs: int | None, dry_run: bool,
                     n_threads: int, charge_override: int | None = None,
                     solvent_pairs: list[tuple[str, str]] | None = None) -> dict:
    """Run the CREST/xTB pipeline for one compound across one or more solvent legs.

    solvent_pairs : optional list of (xtb_solvent, label) legs overriding the default
        [(water, water), (chcl3, mem)]. The FIRST leg is the polar reference used for
        the ΔPSA/ΔHB deltas; every other leg is treated as apolar (its ensemble is
        capped to max_confs, exactly like the old chloroform "mem" leg). Example for a
        water/cyclohexane partition run:
            [("water", "water"), ("cyclohexane", "cyclohexane")]
    """
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
    charge = (charge_override if charge_override is not None
              else sum(atom.GetFormalCharge() for atom in mol.GetAtoms()))
    _log(f"  {'Charge (override)' if charge_override is not None else 'Formal charge'}: "
         f"{charge}  |  Atoms (with H): {mol.GetNumAtoms()}")

    _log(f"  Step 1: RDKit ETKDGv3 embedding (5000 → 50)...")
    try:
        mol_embedded = embed_rdkit_conformers(mol)
    except RuntimeError as e:
        _log(f"  ⚠ Embedding failed: {e}")
        return result
    _log(f"  Step 1 done: {mol_embedded.GetNumConformers()} conformers")
    _template_mol = mol_embedded

    # Auto-detect a diazirine; if present, build the N=N constraint file (1-based indices,
    # consistent with the xyz atom order) once and reuse it for both solvents.
    nn_atoms = diazirine_nn_atoms(mol_embedded)
    constraint_file = None
    if nn_atoms is not None:
        constraint_file = write_constraint_file(
            (work_base / "diazirine_constrain.inp").resolve(), nn_atoms)
        _log(f"  Diazirine detected → constraining N=N (atoms {nn_atoms[0]},{nn_atoms[1]}) "
             f"to {DIAZIRINE_NN} Å [prevents GFN2's spurious ~1.43 Å]")

    if solvent_pairs is None:
        solvent_pairs = [(SOLVENT_AQ, LABEL_AQ), (SOLVENT_MEM, LABEL_MEM)]
    polar_label = solvent_pairs[0][1]   # first leg = polar reference (deltas + no cap)

    failed_solvents = []

    for solvent, label in solvent_pairs:
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

        crest_dir = sol_dir / "crest"
        xyz_in    = crest_dir / f"{short}_{label}_start.xyz"
        existing_ensemble = crest_dir / "crest_conformers.xyz"
        existing_rotamers = crest_dir / "crest_rotamers_0.xyz"

        if existing_ensemble.exists() and existing_ensemble.stat().st_size > 0:
            _log(f"  Steps 2-3: existing CREST ensemble found "
                 f"({existing_ensemble.stat().st_size/1e6:.0f} MB) — skipping xTB + CREST")
            raw_ensemble = existing_ensemble
        elif existing_rotamers.exists() and existing_rotamers.stat().st_size > 0 and xyz_in.exists():
            _log(f"  Steps 2-3: rotamers found ({existing_rotamers.stat().st_size/1e6:.0f} MB)"
                 f" + start xyz exists — skipping xTB, running --cregen")
            try:
                raw_ensemble = run_crest_cregen(xyz_in, crest_dir, solvent, n_threads, charge,
                                                constraint_file=constraint_file)
            except Exception as e:
                print(f"      ⚠ CREST --cregen error: {e}")
                raw_ensemble = None
        else:
            xtb_dir = sol_dir / "xtb_opt"
            _log(f"  Step 2: xTB pre-opt ({mol_embedded.GetNumConformers()} conformers)...")
            mol_min = xtb_preopt_mol(mol_embedded, solvent, xtb_dir, charge,
                                     n_procs=n_threads, constraint_file=constraint_file)
            if mol_min is None:
                _log(f"  ⚠ All xTB jobs failed — skipping {label}")
                result[f"{label}_status"] = "failed"
                result[f"{label}_error"]  = "all xTB jobs failed"
                failed_solvents.append(label)
                continue

            crest_dir.mkdir(parents=True, exist_ok=True)

            if existing_rotamers.exists() and existing_rotamers.stat().st_size > 0:
                _log(f"  Step 3: rotamers found ({existing_rotamers.stat().st_size/1e6:.0f} MB)"
                     f" — running --cregen to skip MTDs")
                _write_conformer_xyz(mol_min, mol_min.GetConformer().GetId(), xyz_in, comment=smi)
                try:
                    raw_ensemble = run_crest_cregen(xyz_in, crest_dir, solvent, n_threads, charge,
                                                constraint_file=constraint_file)
                except Exception as e:
                    print(f"      ⚠ CREST --cregen error: {e}")
                    raw_ensemble = None
            else:
                _write_conformer_xyz(mol_min, mol_min.GetConformer().GetId(), xyz_in, comment=smi)
                _log(f"  Step 3: CREST iMTD-GC ({solvent})...")
                try:
                    raw_ensemble = run_crest(xyz_in, crest_dir, solvent, n_threads, charge,
                                             constraint_file=constraint_file)
                except Exception as e:
                    print(f"      ⚠ CREST error: {e}")
                    raw_ensemble = None

        if raw_ensemble is None:
            result[f"{label}_status"] = "failed"
            result[f"{label}_error"]  = "no ensemble produced"
            failed_solvents.append(label)
            continue

        ensemble_xyz = sol_dir / "ensemble.xyz"
        shutil.copy2(raw_ensemble, ensemble_xyz)
        _log(f"  Ensemble saved → {ensemble_xyz.name}")

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

        if label != polar_label:   # cap apolar legs (was: solvent == SOLVENT_MEM)
            max_post = max_confs if max_confs is not None else 50
            if len(conformers) > max_post:
                _log(f"  Capping ensemble: {len(conformers)} → {max_post} lowest-energy conformers")
            conformers = conformers[:max_post]
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

        sdf_path = sol_dir / "ensemble.sdf"
        if export_sdf(conformers, smi, sdf_path, template_mol=_template_mol):
            _log(f"      SDF written → {sdf_path.name} ({n_confs} conformers)")
        else:
            _log(f"      ⚠ SDF export failed")

        json_path = sol_dir / "ensemble.json"
        export_json(conformers, psa_vals, hb_vals, weights, smi, charge, crest_log, json_path)
        _log(f"      JSON written → {json_path.name}")

        partial_csv = work_base / f"{short}_results_partial.csv"
        pd.DataFrame([result]).to_csv(partial_csv, index=False)
        _log(f"      Checkpoint saved → {partial_csv.name}")

    if failed_solvents:
        _log(f"  WARNING: {short} failed for solvent(s): {', '.join(sorted(set(failed_solvents)))}")
        result["failed_solvents"] = ",".join(sorted(set(failed_solvents)))

    # Deltas between the polar reference (first leg) and the apolar phase (last leg).
    # Result keys are kept generic ("crest_delta_psa", ..._aq/_mem) so downstream tooling
    # is unchanged whether the apolar phase is chloroform ("mem") or cyclohexane.
    apolar_label = solvent_pairs[-1][1]
    aq_key  = f"{polar_label}_psa_boltz"
    mem_key = f"{apolar_label}_psa_boltz"
    if len(solvent_pairs) >= 2 and aq_key in result and mem_key in result:
        result["crest_delta_psa"]      = round(result[aq_key] - result[mem_key], 2)
        aq_lo  = result.get(f"{polar_label}_psa_lowen",  np.nan)
        mem_lo = result.get(f"{apolar_label}_psa_lowen", np.nan)
        result["crest_delta_psa_lowen"] = round(aq_lo - mem_lo, 2) if pd.notna(aq_lo) and pd.notna(mem_lo) else np.nan
        result["crest_delta_hb"]        = round(result[f"{apolar_label}_hb_boltz"] - result[f"{polar_label}_hb_boltz"], 2)
        result["crest_psa_spread_aq"]   = result.get(f"{polar_label}_psa_std",  np.nan)
        result["crest_psa_spread_mem"]  = result.get(f"{apolar_label}_psa_std", np.nan)
        print(f"\n  ✓ CREST ΔPSA (Boltzmann) = {result['crest_delta_psa']:.1f} Å²  "
              f"(low-energy = {result['crest_delta_psa_lowen']:.1f} Å²  "
              f"ΔHB = {result['crest_delta_hb']:.2f})")

    return result


# ── Resume logic ──────────────────────────────────────────────────────────────
def find_resume_dir(runs_base: Path, compound_idx: int, short: str) -> Path | None:
    """Return the most recent incomplete run dir (has checkpoint data but no final CSV)."""
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
        reverse=True,
    )
    for d in candidates:
        final_csv = d / f"{short}_results.csv"
        if final_csv.exists():
            continue
        if any((d / p).exists() for p in checkpoint_paths):
            return d
    return None

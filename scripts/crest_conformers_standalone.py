# env: chameleon-sim
"""
crest_conformers_standalone.py
------------------------------
Trimmed CREST iMTD-GC conformer sampling pipeline.

This copy keeps only the direct-SMILES path needed for a standalone notebook or
command-line run:

  SMILES -> RDKit embedding -> xTB pre-opt -> CREST -> SDF/JSON export

The retained CREST ensemble is then reduced to a diverse set of unique
conformers so downstream geometry summaries can use min/max/median/range
without overweighting near-duplicate structures.

It intentionally omits the reference-compound registry, resume logic, checkpoint
CSV writing, and other notebook/CLI plumbing that are not required for a fresh
SMILES input.
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

import numpy as np

try:
    from phys_descriptors_v2 import compute_psa_xyz, count_hbonds_xyz
except ImportError as exc:
    raise ImportError(
        "crest_conformers_standalone.py requires `phys_descriptors_v2` to be importable. "
        "Install or expose the Chameleon_Predictor `scripts/` directory on PYTHONPATH, "
        "or copy the descriptor helpers into the same folder as this script."
    ) from exc

warnings.filterwarnings("ignore")

SOLVENT_AQ = "water"
SOLVENT_MEM = "chcl3"
LABEL_AQ = "water"
LABEL_MEM = "mem"
BOLTZMANN_RT_KCAL = 0.592
DEFAULT_UNIQUE_CONFS = 20


def _log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def _aligned_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    """Return Kabsch-aligned RMSD for two coordinate arrays with matching atoms."""
    a = np.asarray(coords_a, dtype=float)
    b = np.asarray(coords_b, dtype=float)
    if a.shape != b.shape or a.ndim != 2 or a.shape[1] != 3:
        return float("inf")
    a0 = a - a.mean(axis=0)
    b0 = b - b.mean(axis=0)
    h = a0.T @ b0
    u, _, vt = np.linalg.svd(h)
    rot = vt.T @ u.T
    if np.linalg.det(rot) < 0:
        vt[-1, :] *= -1
        rot = vt.T @ u.T
    diff = a0 @ rot - b0
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def deduplicate_conformers(conformers: list, max_keep: int = DEFAULT_UNIQUE_CONFS,
                           rmsd_thresh: float = 0.5) -> list:
    """Keep a diverse, energy-ranked subset of unique conformers."""
    ordered = sorted(
        conformers,
        key=lambda item: item[2] if np.isfinite(item[2]) else float("inf"),
    )
    unique = []
    for symbols, coords, energy in ordered:
        if len(unique) >= max_keep:
            break
        if all(_aligned_rmsd(coords, kept_coords) >= rmsd_thresh for _, kept_coords, _ in unique):
            unique.append((symbols, coords, energy))
    return unique


def boltzmann_weights(energies, rt_kcal: float = BOLTZMANN_RT_KCAL) -> np.ndarray:
    """Compute normalized Boltzmann weights from a sequence of energies."""
    e = np.asarray(energies, dtype=float)
    mask = np.isfinite(e)
    weights = np.full_like(e, np.nan, dtype=float)
    if not mask.any():
        return weights
    shifted = e[mask] - np.nanmin(e[mask])
    raw = np.exp(-shifted / rt_kcal)
    weights[mask] = raw / raw.sum()
    return weights


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
    for conf in confs:
        ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=conf.GetId())
        if ff is None:
            mmff_fail_count += 1
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


def _write_conformer_xyz(mol, conf_id: int, xyz_path: Path, comment: str = "") -> None:
    conf = mol.GetConformer(conf_id)
    pos = conf.GetPositions()
    syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(mol.GetNumAtoms())]
    lines = [str(mol.GetNumAtoms()), comment]
    for sym, (x, y, z) in zip(syms, pos):
        lines.append(f"{sym}  {x:.6f}  {y:.6f}  {z:.6f}")
    xyz_path.write_text("\n".join(lines) + "\n")


DIAZIRINE_SMARTS = "[#6]1[#7]=[#7]1"
DIAZIRINE_NN = 1.23


def diazirine_nn_atoms(mol):
    """1-based (N1, N2) indices of a diazirine N=N, or None if absent."""
    from rdkit import Chem

    matches = mol.GetSubstructMatches(Chem.MolFromSmarts(DIAZIRINE_SMARTS))
    if not matches:
        return None
    _c, n1, n2 = matches[0]
    return n1 + 1, n2 + 1


def write_constraint_file(path: Path, nn_pair, value: float = DIAZIRINE_NN,
                          fc: float = 0.25) -> Path:
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
        cmd.extend(["--input", str(constraint_file)])

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
                new_conf.SetAtomPosition(i, [float(parts[1]), float(parts[2]), float(parts[3])])
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
    results = []
    with multiprocessing.Pool(n_workers) as pool:
        for result in pool.imap_unordered(_xtb_opt_worker, params_list):
            results.append(result)

    valid = [r for r in results if r is not None]
    if not valid:
        return None

    _, mol_min, _ = min(valid, key=lambda x: x[2])
    return mol_min


def run_crest(xyz_path: Path, work_dir: Path, solvent: str,
              n_threads: int, charge: int = 0, constraint_file: Path | None = None) -> Path | None:
    work_dir.mkdir(parents=True, exist_ok=True)
    xyz_path = xyz_path.resolve()

    cmd = [
        "crest", str(xyz_path),
        "-T", str(n_threads),
        "--gfn2",
        "--chrg", str(charge),
        "--alpb", solvent,
        "--keepdir",
        "--noreftopo",
        "-notopo",
    ]
    if constraint_file:
        cmd.extend(["--cinp", str(constraint_file)])

    _log(f"    CREST: {' '.join(cmd)}")
    out_path = work_dir / "crest.out"
    err_path = work_dir / "crest.err"
    with open(out_path, "w") as fout, open(err_path, "w") as ferr:
        proc = subprocess.run(cmd, cwd=work_dir.resolve(), stdout=fout, stderr=ferr)

    ensemble = work_dir.resolve() / "crest_conformers.xyz"
    if ensemble.exists() and ensemble.stat().st_size > 0:
        return ensemble
    if proc.returncode != 0:
        _log(f"    CREST: FAILED (exit={proc.returncode})")
    return None


def run_crest_cregen(xyz_path: Path, work_dir: Path, solvent: str,
                     n_threads: int, charge: int = 0, constraint_file: Path | None = None) -> Path | None:
    rotamers = work_dir.resolve() / "crest_rotamers_0.xyz"
    if not rotamers.exists() or rotamers.stat().st_size == 0:
        _log(f"    CREST --cregen: crest_rotamers_0.xyz not found in {work_dir}")
        return None

    cmd = [
        "crest", str(xyz_path.resolve()),
        "--cregen",
        "-T", str(n_threads),
        "--gfn2",
        "--chrg", str(charge),
        "--alpb", solvent,
    ]
    if constraint_file:
        cmd.extend(["--cinp", str(constraint_file)])

    _log(f"    CREST --cregen: {' '.join(cmd)}")
    out_path = work_dir / "crest_cregen.out"
    err_path = work_dir / "crest_cregen.err"
    with open(out_path, "w") as fout, open(err_path, "w") as ferr:
        proc = subprocess.run(cmd, cwd=work_dir.resolve(), stdout=fout, stderr=ferr)

    ensemble = work_dir.resolve() / "crest_conformers.xyz"
    if proc.returncode != 0:
        _log(f"    CREST --cregen: FAILED (exit={proc.returncode})")
        return None
    if ensemble.exists() and ensemble.stat().st_size > 0:
        return ensemble
    return None


def parse_crest_log(crest_out_path: Path) -> dict | None:
    """Parse CREST thermodynamics and per-conformer weights from crest.out."""
    ensemble_patterns = {
        "ensembleenergy": ("ensemble average energy", float),
        "ensembleentropy": ("ensemble entropy", float),
        "ensemblefreeenergy": ("ensemble free energy", float),
        "lowestenergy": ("E lowest", float),
        "poplowestpct": ("population of lowest in %", float),
        "temperature": ("T /K", float),
        "uniqueconfs": ("number of unique conformers", int),
        "totalconfs": ("total number unique points", int),
    }
    conf_cols = {
        "relativeenergy": 1,
        "totalenergy": 2,
        "conformerweight": 3,
        "boltzmannweight": 4,
        "set": 5,
        "degeneracy": 6,
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
                                "relativeenergy": float(parts[conf_cols["relativeenergy"]]),
                                "totalenergy": float(parts[conf_cols["totalenergy"]]),
                                "conformerweights": [float(parts[conf_cols["conformerweight"]])],
                                "boltzmannweight": float(parts[conf_cols["boltzmannweight"]]),
                                "set": int(parts[conf_cols["set"]]),
                                "degeneracy": int(parts[conf_cols["degeneracy"]]),
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


def compare_boltzmann_weights(calculated: np.ndarray, reported: np.ndarray,
                              atol: float = 1e-3, rtol: float = 1e-2) -> dict:
    """Compare calculated and reported weights and return a small diagnostics dict."""
    calc = np.asarray(calculated, dtype=float)
    rep = np.asarray(reported, dtype=float)
    n = min(len(calc), len(rep))
    if n == 0:
        return {"n": 0, "max_abs_diff": np.nan, "mean_abs_diff": np.nan, "ok": False}
    diff = np.abs(calc[:n] - rep[:n])
    ok = bool(np.allclose(calc[:n], rep[:n], atol=atol, rtol=rtol, equal_nan=False))
    return {
        "n": int(n),
        "max_abs_diff": float(np.nanmax(diff)),
        "mean_abs_diff": float(np.nanmean(diff)),
        "ok": ok,
    }


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


def export_sdf(conformers: list, smiles: str, out_path: Path, template_mol=None) -> bool:
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


def export_json(conformers: list, psa_vals: list, hb_vals: list,
                weights: np.ndarray, smiles: str, charge: int, out_path: Path) -> None:
    conf_records = []
    for i, (_, _, energy) in enumerate(conformers):
        conf_records.append({
            "totalenergy": float(energy) if np.isfinite(energy) else None,
            "boltzmannweight": float(weights[i]) if np.isfinite(weights[i]) else None,
            "psa": psa_vals[i],
            "hbonds": int(hb_vals[i]),
        })

    valid_w = np.isfinite(weights)
    psa_arr = np.array(psa_vals)
    hb_arr = np.array(hb_vals, dtype=float)
    energies = np.array([energy for _, _, energy in conformers], dtype=float)
    data = {
        "smiles": smiles,
        "charge": charge,
        "geometry_mode": "unique_deduplicated",
        "n_confs": len(conformers),
        "energy_min": float(np.nanmin(energies)) if np.isfinite(energies).any() else None,
        "energy_median": float(np.nanmedian(energies)) if np.isfinite(energies).any() else None,
        "energy_max": float(np.nanmax(energies)) if np.isfinite(energies).any() else None,
        "boltzmann_psa": round(float(np.dot(weights[valid_w], psa_arr[valid_w])), 2),
        "boltzmann_hb": round(float(np.dot(weights[valid_w], hb_arr[valid_w])), 2),
        "lowen_psa": psa_vals[0],
        "lowen_hb": int(hb_vals[0]),
        "conformers": conf_records,
    }

    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)


def check_binaries(require=("xtb", "crest")) -> dict:
    found = {b: shutil.which(b) for b in require}
    missing = [b for b, p in found.items() if p is None]
    if missing:
        raise RuntimeError(
            f"Required external binary/binaries not found on PATH: {', '.join(missing)}. "
            f"Install with `conda install -c conda-forge xtb crest` or load your cluster's modules."
        )
    return found


def _safe_short(name: str) -> str:
    s = re.sub(r"[^0-9A-Za-z._-]+", "_", (name or "").strip()).strip("_")
    return (s or "molecule")[:40]


def generate_ensembles(smiles: str, name: str = "molecule",
                       outdir: str | Path = "results/notebook_runs",
                       charge: int | None = None, n_threads: int | None = None,
                       max_confs: int | None = None,
                       check_binaries_first: bool = True) -> dict:
    from rdkit import Chem

    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("`smiles` must be a non-empty string.")
    if Chem.MolFromSmiles(smiles) is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")

    if check_binaries_first:
        check_binaries()

    short = _safe_short(name)
    n_threads = n_threads or os.cpu_count() or 1
    work_base = Path(outdir) / f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{short}"
    work_base.mkdir(parents=True, exist_ok=True)
    _log(f"Direct-SMILES run: {name}  →  {work_base}")

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    charge = charge if charge is not None else sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    _log(f"  Formal charge: {charge}  |  Atoms (with H): {mol.GetNumAtoms()}")

    _log("  Step 1: RDKit ETKDGv3 embedding (5000 → 50)...")
    mol_embedded = embed_rdkit_conformers(mol)
    _log(f"  Step 1 done: {mol_embedded.GetNumConformers()} conformers")

    template_mol = mol_embedded
    nn_atoms = diazirine_nn_atoms(mol_embedded)
    constraint_file = None
    if nn_atoms is not None:
        constraint_file = write_constraint_file((work_base / "diazirine_constrain.inp").resolve(), nn_atoms)
        _log(f"  Diazirine detected → constraining N=N (atoms {nn_atoms[0]},{nn_atoms[1]})")

    result = {"compound": name, "run_id": work_base.name, "smiles": smiles, "charge": charge}

    for solvent, label in [(SOLVENT_AQ, LABEL_AQ), (SOLVENT_MEM, LABEL_MEM)]:
        sol_dir = work_base / label
        sol_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n  ── [{label}] solvent={solvent} ──")

        crest_dir = sol_dir / "crest"
        xyz_in = crest_dir / f"{short}_{label}_start.xyz"
        existing_ensemble = crest_dir / "crest_conformers.xyz"
        existing_rotamers = crest_dir / "crest_rotamers_0.xyz"

        if existing_ensemble.exists() and existing_ensemble.stat().st_size > 0:
            raw_ensemble = existing_ensemble
        elif existing_rotamers.exists() and existing_rotamers.stat().st_size > 0 and xyz_in.exists():
            raw_ensemble = run_crest_cregen(xyz_in, crest_dir, solvent, n_threads, charge, constraint_file)
        else:
            xtb_dir = sol_dir / "xtb_opt"
            _log(f"  Step 2: xTB pre-opt ({mol_embedded.GetNumConformers()} conformers)...")
            mol_min = xtb_preopt_mol(mol_embedded, solvent, xtb_dir, charge, n_procs=n_threads,
                                     constraint_file=constraint_file)
            if mol_min is None:
                raise RuntimeError(f"All xTB jobs failed for {label}")

            crest_dir.mkdir(parents=True, exist_ok=True)
            _write_conformer_xyz(mol_min, mol_min.GetConformer().GetId(), xyz_in, comment=smiles)
            _log(f"  Step 3: CREST iMTD-GC ({solvent})...")
            raw_ensemble = run_crest(xyz_in, crest_dir, solvent, n_threads, charge, constraint_file)

        if raw_ensemble is None:
            raise RuntimeError(f"No ensemble produced for {label}")

        ensemble_xyz = sol_dir / "ensemble.xyz"
        shutil.copy2(raw_ensemble, ensemble_xyz)
        conformers, _energy_fail = parse_xyz_ensemble(ensemble_xyz)
        if not conformers:
            raise RuntimeError(f"Empty ensemble for {label}")

        max_post = max_confs if max_confs is not None else DEFAULT_UNIQUE_CONFS
        conformers = deduplicate_conformers(conformers, max_keep=max_post, rmsd_thresh=0.5)
        if not conformers:
            raise RuntimeError(f"No unique conformers retained for {label}")

        psa_vals, hb_vals, energies = [], [], []
        for syms, crds, eng in conformers:
            psa_vals.append(compute_psa_xyz(syms, crds, template_mol=template_mol))
            hb_vals.append(count_hbonds_xyz(syms, crds))
            energies.append(eng)

        weights = boltzmann_weights(energies)
        crest_log = None
        for log_name in ("crest.out", "crest_cregen.out"):
            candidate = crest_dir / log_name
            if candidate.exists():
                crest_log = parse_crest_log(candidate)
                if crest_log is not None:
                    break

        weight_check = None
        if crest_log and crest_log.get("conformers"):
            reported = [conf.get("boltzmannweight", np.nan) for conf in crest_log["conformers"]]
            weight_check = compare_boltzmann_weights(weights, reported)
            if not weight_check["ok"]:
                _log(
                    f"  [{label}] Boltzmann safeguard mismatch: max_abs_diff={weight_check['max_abs_diff']:.3g}, "
                    f"mean_abs_diff={weight_check['mean_abs_diff']:.3g}"
                )

        sdf_path = sol_dir / "ensemble.sdf"
        json_path = sol_dir / "ensemble.json"
        if not export_sdf(conformers, smiles, sdf_path, template_mol=template_mol):
            raise RuntimeError(f"SDF export failed for {label}")
        export_json(conformers, psa_vals, hb_vals, weights, smiles, charge, json_path)

        result[f"{label}_n_confs"] = len(conformers)
        result[f"{label}_status"] = "ok"
        result[f"{label}_sdf"] = str(sdf_path)
        result[f"{label}_json"] = str(json_path)
        if weight_check is not None:
            result[f"{label}_boltzmann_check"] = weight_check

    result["ok"] = all((work_base / sol / "ensemble.sdf").exists() and (work_base / sol / "ensemble.json").exists()
                        for sol in (SOLVENT_AQ, SOLVENT_MEM))
    result["work_dir"] = str(work_base)
    return result


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Standalone CREST/xTB conformer generator")
    parser.add_argument("smiles", help="Molecule SMILES")
    parser.add_argument("--name", default="molecule", help="Label for the run")
    parser.add_argument("--outdir", default="results/notebook_runs", type=Path)
    parser.add_argument("--charge", default=None, type=int, help="Formal charge override")
    parser.add_argument("--threads", default=None, type=int, help="CPU cores to use")
    parser.add_argument("--max-confs", default=None, type=int, help="Cap unique conformers retained per solvent")
    args = parser.parse_args()

    out = generate_ensembles(
        args.smiles,
        name=args.name,
        outdir=args.outdir,
        charge=args.charge,
        n_threads=args.threads,
        max_confs=args.max_confs,
    )
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
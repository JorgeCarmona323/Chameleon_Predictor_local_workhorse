# env: chameleon-sim
"""
crest_engine.py
---------------
Conformer-generation engine.  SMILES → RDKit ETKDG embedding → xTB pre-opt →
CREST iMTD-GC sampling, exported as XYZ + SDF per solvent.

This is a *pure conformer generator*: no descriptors, no PSA/H-bond/Boltzmann
averaging, no ΔPSA/ΔHB, no permeability/PAMPA metadata, no dry-run mockups.
Descriptor calculation lives downstream (see notebooks/pipeline).

Public API
    generate_conformers(smiles, ...)                 registry-free front end
    process_molecule(smiles, name, work_base, ...)   one molecule, N solvent legs
    check_binaries(), _safe_short()                  helpers
    find_resume_dir(...)                             resume a prior incomplete run

Imported by crest_v3.2.py (reference-compound registry + CLI).
"""

from __future__ import annotations

import hashlib
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
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np

warnings.filterwarnings("ignore")

# Solvent legs are (xtb/CREST solvent keyword, output-directory label).
# The label names the sub-folder; the solvent string is the xtb/CREST --alpb keyword.
SOLVENT_PAIRS_DEFAULT: list[tuple[str, str]] = [
    ("water",   "water"),        # folder: water
    ("chcl3",   "chloroform"),   # folder: chloroform, solvent keyword: chcl3
    ("hexane",  "hexane"),       # folder: hexane, solvent keyword: hexane
]

SOURCE_TAG = "CREST GFN2-xTB ALPB"


def _log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


# ── binary / filename helpers ─────────────────────────────────────────────────
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


# ── Cached embedding I/O ──────────────────────────────────────────────────────
# The ETKDG embedding (5000 → 50 diverse conformers) is solvent-independent and the most
# expensive CPU step. Persisting it lets every solvent leg — including separate runs / SLURM
# tasks for the same molecule — start from the SAME seeds (consistency) without repeating the
# embedding (speed). load_or_embed(): if the cache file exists it is reused, else generated.
def _smiles_key(smiles: str) -> str:
    """Short stable hash of the (canonical) SMILES, so a cached embedding is reused only for
    the exact same molecule (guards a short-name collision or an edited SMILES)."""
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    canonical = Chem.MolToSmiles(mol) if mol is not None else smiles
    return hashlib.sha1(canonical.encode("utf-8")).hexdigest()[:10]


def write_embedding(mol, path: Path) -> None:
    """Persist the embedded multi-conformer molecule to SDF (atoms + bonds + all 3D confs).
    Written to a temp file then atomically renamed, so a concurrent reader never sees a
    partial file."""
    from rdkit import Chem
    tmp = path.with_name(path.name + ".tmp")
    writer = Chem.SDWriter(str(tmp))
    for cid in range(mol.GetNumConformers()):
        writer.write(mol, confId=cid)
    writer.close()
    tmp.replace(path)


def read_embedding(path: Path):
    """Load a molecule previously written by write_embedding, rebuilt with all conformers.
    Returns None if the file is unreadable or empty (caller then regenerates)."""
    from rdkit import Chem
    try:
        supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
        mols = [m for m in supplier if m is not None]
    except Exception:
        return None
    if not mols:
        return None
    base = Chem.Mol(mols[0])
    base.RemoveAllConformers()
    for m in mols:
        if m.GetNumConformers() and m.GetNumAtoms() == base.GetNumAtoms():
            base.AddConformer(m.GetConformer(), assignId=True)
    return base if base.GetNumConformers() else None


def load_or_embed(mol, cache_path: Path | None):
    """Reuse a cached embedding if one exists for this molecule, else embed and cache it.
    cache_path=None disables caching (embed every time). May raise RuntimeError on embed
    failure (propagated to the caller)."""
    if cache_path is not None:
        cache_path = Path(cache_path)
        if cache_path.exists() and cache_path.stat().st_size > 0:
            cached = read_embedding(cache_path)
            if cached is not None and cached.GetNumAtoms() == mol.GetNumAtoms():
                _log(f"  Step 1: reusing cached embedding "
                     f"({cached.GetNumConformers()} confs) ← {cache_path.name}")
                return cached
            _log(f"  Step 1: cached embedding {cache_path.name} unusable — regenerating")

    mol_embedded = embed_rdkit_conformers(mol)   # may raise RuntimeError

    if cache_path is not None:
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            write_embedding(mol_embedded, cache_path)
            _log(f"  Step 1: cached embedding → {cache_path}")
        except Exception as e:
            _log(f"  Step 1: could not cache embedding ({e}) — continuing without cache")
    return mol_embedded


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


# ── diazirine N=N constraint (optional; GFN2 stretches it to a spurious ~1.43 Å; ─
#    see docs/experiments + memory diazirine-review-checklist). True N=N = 1.228 Å (exp,
#    microwave) / 1.230 Å (CCSD(T)). Constrain it so the macrocycle samples freely
#    while the rigid diazirine can't drift into GFN2's single-bond basin. Auto-detected
#    per molecule and gated by `use_diazirine_constraint`; a no-op when absent.
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


# ── Metadata export (conformer provenance only — no descriptors) ──────────────
def write_metadata(path: Path, *, smiles: str, name: str, charge: int,
                   solvent: str, label: str, n_conformers: int) -> None:
    """Write a small provenance JSON for one solvent leg. Intentionally carries no
    PSA/H-bond/Boltzmann/permeability values — descriptors are computed downstream."""
    with open(path, "w") as f:
        json.dump({
            "smiles":       smiles,
            "name":         name,
            "charge":       charge,
            "solvent":      solvent,
            "label":        label,
            "n_conformers": n_conformers,
            "source":       SOURCE_TAG,
        }, f, indent=2)


# ── One solvent leg ───────────────────────────────────────────────────────────
def _run_solvent_leg(mol_embedded, template_mol, smiles: str, name: str, short: str,
                     work_base: Path, solvent: str, label: str, charge: int,
                     n_threads: int, max_confs: int | None,
                     constraint_file: Path | None) -> dict:
    """Generate one solvent's ensemble. Directory layout (kept stable for resume):
        <work_base>/<label>/xtb_opt/        xTB pre-opt scratch
        <work_base>/<label>/crest/          CREST run dir (crest.out/err, rotamers, ...)
        <work_base>/<label>/ensemble.xyz    final multi-conformer XYZ
        <work_base>/<label>/ensemble.sdf    final SDF
        <work_base>/<label>/metadata.json   provenance
    """
    sol_dir = work_base / label
    sol_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n  ── [{label}] solvent={solvent} ──")

    crest_dir = sol_dir / "crest"
    xyz_in    = crest_dir / f"{short}_{label}_start.xyz"
    existing_ensemble = crest_dir / "crest_conformers.xyz"
    existing_rotamers = crest_dir / "crest_rotamers_0.xyz"

    raw_ensemble = None
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
            return {"solvent": solvent, "label": label,
                    "status": "failed", "error": "all xTB jobs failed"}

        crest_dir.mkdir(parents=True, exist_ok=True)

        if existing_rotamers.exists() and existing_rotamers.stat().st_size > 0:
            _log(f"  Step 3: rotamers found ({existing_rotamers.stat().st_size/1e6:.0f} MB)"
                 f" — running --cregen to skip MTDs")
            _write_conformer_xyz(mol_min, mol_min.GetConformer().GetId(), xyz_in, comment=smiles)
            try:
                raw_ensemble = run_crest_cregen(xyz_in, crest_dir, solvent, n_threads, charge,
                                                constraint_file=constraint_file)
            except Exception as e:
                print(f"      ⚠ CREST --cregen error: {e}")
                raw_ensemble = None
        else:
            _write_conformer_xyz(mol_min, mol_min.GetConformer().GetId(), xyz_in, comment=smiles)
            _log(f"  Step 3: CREST iMTD-GC ({solvent})...")
            try:
                raw_ensemble = run_crest(xyz_in, crest_dir, solvent, n_threads, charge,
                                         constraint_file=constraint_file)
            except Exception as e:
                print(f"      ⚠ CREST error: {e}")
                raw_ensemble = None

    if raw_ensemble is None:
        return {"solvent": solvent, "label": label,
                "status": "failed", "error": "no ensemble produced"}

    ensemble_xyz = sol_dir / "ensemble.xyz"
    shutil.copy2(raw_ensemble, ensemble_xyz)
    _log(f"  Ensemble saved → {ensemble_xyz.name}")

    conformers, energy_fail = parse_xyz_ensemble(ensemble_xyz)
    n_confs_full = len(conformers)
    print(f"      Parsed {n_confs_full} conformers")
    if energy_fail:
        _log(f"      Warning: {energy_fail}/{n_confs_full} conformers missing energy")
    if n_confs_full == 0:
        return {"solvent": solvent, "label": label,
                "status": "failed", "error": "empty ensemble"}

    # Optional cap to the N lowest-energy conformers (CREST output is energy-ordered).
    if max_confs is not None and len(conformers) > max_confs:
        _log(f"  Capping ensemble: {len(conformers)} → {max_confs} lowest-energy conformers")
        conformers = conformers[:max_confs]
    n_confs = len(conformers)

    sdf_path = sol_dir / "ensemble.sdf"
    if export_sdf(conformers, smiles, sdf_path, template_mol=template_mol):
        _log(f"      SDF written → {sdf_path.name} ({n_confs} conformers)")
    else:
        _log("      ⚠ SDF export failed")

    meta_path = sol_dir / "metadata.json"
    write_metadata(meta_path, smiles=smiles, name=name, charge=charge,
                   solvent=solvent, label=label, n_conformers=n_confs)
    _log(f"      Metadata written → {meta_path.name}")

    return {
        "solvent":           solvent,
        "label":             label,
        "status":            "ok",
        "n_conformers":      n_confs,
        "n_conformers_full": n_confs_full,
        "ensemble_xyz":      str(ensemble_xyz),
        "ensemble_sdf":      str(sdf_path),
        "metadata":          str(meta_path),
    }


# ── Process one molecule across solvent legs ──────────────────────────────────
def process_molecule(smiles: str, name: str, work_base: Path,
                     solvent_pairs: list[tuple[str, str]] | None = None,
                     charge: int | None = None, n_threads: int = 1,
                     max_confs: int | None = None,
                     use_diazirine_constraint: bool = True,
                     embed_cache_dir: str | Path | None = None) -> dict:
    """Generate CREST/xTB conformer ensembles for one molecule across one or more
    solvent legs.

    solvent_pairs : list of (xtb/CREST solvent, output-directory label). Defaults to
        SOLVENT_PAIRS_DEFAULT (water/chloroform/hexane). The label names the
        sub-folder; the solvent string is the xtb/CREST --alpb keyword.
    charge : formal-charge override; default = sum of RDKit formal charges.
    use_diazirine_constraint : auto-detect a diazirine and pin its N=N distance
        (GFN2 otherwise drifts it to ~1.43 Å). A no-op when no diazirine is present.
    embed_cache_dir : if given, the RDKit embedding is cached here as
        <short>_<smiles-hash>.sdf and reused by any later run/solvent leg for the same
        molecule (embed once, reuse everywhere). None disables caching.
    """
    from rdkit import Chem

    if solvent_pairs is None:
        solvent_pairs = list(SOLVENT_PAIRS_DEFAULT)

    short = _safe_short(name)
    print(f"\n{'─'*60}")
    print(f"  {name}")

    result: dict = {
        "name":     name,
        "smiles":   smiles,
        "run_id":   work_base.name,
        "solvents": {},
    }

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("  ⚠ Invalid SMILES")
        result["status"] = "failed"
        result["error"]  = "invalid SMILES"
        return result
    mol = Chem.AddHs(mol)
    charge = (charge if charge is not None
              else sum(atom.GetFormalCharge() for atom in mol.GetAtoms()))
    result["charge"] = charge
    _log(f"  Charge: {charge}  |  Atoms (with H): {mol.GetNumAtoms()}")

    work_base.mkdir(parents=True, exist_ok=True)

    cache_path = None
    if embed_cache_dir is not None:
        cache_path = Path(embed_cache_dir) / f"{short}_{_smiles_key(smiles)}.sdf"

    _log("  Step 1: RDKit ETKDGv3 embedding (5000 → 50)...")
    try:
        mol_embedded = load_or_embed(mol, cache_path)
    except RuntimeError as e:
        _log(f"  ⚠ Embedding failed: {e}")
        result["status"] = "failed"
        result["error"]  = f"embedding failed: {e}"
        return result
    _log(f"  Step 1 done: {mol_embedded.GetNumConformers()} conformers")
    template_mol = mol_embedded

    # Auto-detect a diazirine; if present (and enabled), build the N=N constraint file
    # once (1-based indices, matching the xyz atom order) and reuse it for every solvent.
    constraint_file = None
    if use_diazirine_constraint:
        nn_atoms = diazirine_nn_atoms(mol_embedded)
        if nn_atoms is not None:
            constraint_file = write_constraint_file(
                (work_base / "diazirine_constrain.inp").resolve(), nn_atoms)
            _log(f"  Diazirine detected → constraining N=N (atoms {nn_atoms[0]},{nn_atoms[1]}) "
                 f"to {DIAZIRINE_NN} Å [prevents GFN2's spurious ~1.43 Å]")

    for solvent, label in solvent_pairs:
        result["solvents"][label] = _run_solvent_leg(
            mol_embedded, template_mol, smiles, name, short, work_base,
            solvent, label, charge, n_threads, max_confs, constraint_file)

    failed = [lbl for lbl, info in result["solvents"].items() if info.get("status") != "ok"]
    if failed:
        result["failed_solvents"] = ",".join(failed)
        _log(f"  WARNING: {short} failed for solvent(s): {', '.join(failed)}")
    result["status"] = "ok" if not failed else "partial"
    return result


# ── Registry-free, direct-SMILES front end ────────────────────────────────────
def generate_conformers(smiles: str, name: str = "molecule",
                        outdir: str | Path = "results/notebook_runs",
                        charge: int | None = None, n_threads: int | None = None,
                        max_confs: int | None = None,
                        solvent_pairs: list[tuple[str, str]] | None = None,
                        use_diazirine_constraint: bool = True,
                        embed_cache_dir: str | Path | None = None,
                        check_binaries_first: bool = True) -> dict:
    """Generate CREST/xTB conformer ensembles for an **arbitrary SMILES** across the
    given solvent legs — the registry-free front end (used by the notebook).

    Parameters
    ----------
    smiles        : the molecule to sample (validated with RDKit; ValueError if unparseable)
    name          : optional label (used for the run-dir and file names)
    outdir        : where the run directory is created
    charge        : optional formal-charge override (default: auto-derived from the SMILES)
    n_threads     : CPU cores for xTB/CREST (default: all cores)
    max_confs     : optional cap on conformers kept per solvent (default: keep all)
    solvent_pairs : list of (xtb/CREST solvent, folder label); default SOLVENT_PAIRS_DEFAULT
    use_diazirine_constraint : pin a detected diazirine N=N (no-op if absent)
    embed_cache_dir : where to cache/reuse the RDKit embedding (default: <outdir>/embeddings,
        shared across runs so each molecule is embedded once). Pass a different path to
        relocate it; the engine's process_molecule treats None as "no cache".
    check_binaries_first     : fail fast with a clear message when xtb/crest are missing

    Returns
    -------
    dict with: work_dir, result (the raw engine record), a `solvents` map of
    {label: {xyz, sdf, metadata, n_conformers, status}}, and `ok`
    (True iff every solvent leg produced an ensemble).
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

    if solvent_pairs is None:
        solvent_pairs = list(SOLVENT_PAIRS_DEFAULT)

    # 3. fresh, timestamped run directory; embeddings cached at a stable, shared location
    short = _safe_short(name)
    n_threads = n_threads or os.cpu_count() or 1
    if embed_cache_dir is None:
        embed_cache_dir = Path(outdir) / "embeddings"
    work_base = Path(outdir) / f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{short}"
    work_base.mkdir(parents=True, exist_ok=True)
    _log(f"Conformer run: {name}  →  {work_base}")

    # 4. run the engine
    result = process_molecule(smiles, name, work_base, solvent_pairs=solvent_pairs,
                              charge=charge, n_threads=n_threads, max_confs=max_confs,
                              use_diazirine_constraint=use_diazirine_constraint,
                              embed_cache_dir=embed_cache_dir)

    # 5. collect outputs
    solvents_out = {}
    for label, info in result["solvents"].items():
        solvents_out[label] = {
            "xyz":          Path(info["ensemble_xyz"]) if info.get("ensemble_xyz") else None,
            "sdf":          Path(info["ensemble_sdf"]) if info.get("ensemble_sdf") else None,
            "metadata":     Path(info["metadata"]) if info.get("metadata") else None,
            "n_conformers": info.get("n_conformers"),
            "status":       info.get("status"),
        }
    ok = bool(result["solvents"]) and all(v["status"] == "ok" for v in solvents_out.values())
    return {"work_dir": work_base, "result": result, "solvents": solvents_out, "ok": ok}


# ── Resume logic ──────────────────────────────────────────────────────────────
def find_resume_dir(runs_base: Path, compound_idx: int, short: str,
                    labels: tuple[str, ...] = ("water", "chloroform", "hexane")) -> Path | None:
    """Return the most recent incomplete run dir (has CREST checkpoint data but no
    final manifest)."""
    if not runs_base.exists():
        return None
    checkpoint_paths = (
        [f"{label}/crest/crest_conformers.xyz" for label in labels]
        + [f"{label}/crest/crest_rotamers_0.xyz" for label in labels]
    )
    candidates = sorted(
        runs_base.glob(f"run_*_{compound_idx}_{short}"),
        reverse=True,
    )
    for d in candidates:
        if (d / f"{short}_manifest.json").exists():
            continue
        if any((d / p).exists() for p in checkpoint_paths):
            return d
    return None

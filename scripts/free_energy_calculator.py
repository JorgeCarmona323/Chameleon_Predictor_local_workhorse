#!/usr/bin/env python3
"""free_energy_calculator — per-conformer free energies + ΔG_transfer from CREST ensembles.

Scores an *already-generated* conformer ensemble (no re-search). For each phase you
pass its OWN native ensemble (phase-specific — CREST was run separately per solvent),
runs an xTB single-point per conformer in that solvent, and returns:

  * per-conformer energies + Boltzmann weights (within each phase)
  * per-phase ensemble free energy  G_ens = -RT ln SUM_i exp(-E_i / RT)
  * ΔG_transfer(ref -> phase) = G_ens(phase) - G_ens(ref)   [ref default = water]

ΔG_transfer (= partition free energy, the FlexiSol "log K_a/b" quantity) is the
permeability descriptor; it benefits from cross-solvent error cancellation, so it is
more robust than any single absolute solvation energy. Populations come from the
energies (no geometry cherry-picking); geometry-first distributions live elsewhere.

METHOD (scoring solvation model, single-point only):
  --method cpcmx   (default)  Grimme's CPCM-X, native to xTB (no ORCA/COSMO-RS needed).
                              More accurate for partition ratios (FlexiSol). No gradient
                              -> single-point only, which is exactly what we want here.
  --method alpb               Analytical linearized PB. Faster, less accurate for
                              partition ratios; kept as a comparison / fallback.

Geometry is NOT re-optimized: CREST already relaxed each phase in ALPB, so we score on
those phase-specific geometries. (CPCM-X cannot optimize anyway — no gradient.)

Apolar phase of record is CYCLOHEXANE (Ono 2019: cyclohexane SASA >> chloroform for
permeability). If your xTB build lacks a cyclohexane keyword, use the closest supported
saturated-hydrocarbon surrogate (hexane ε≈1.88 / hexadecane ε≈2.06) and say so.

Runs on the HPC (xTB is Linux-only). Pure stdlib — no rdkit/pandas.

Examples
--------
  # phase-specific CPCM-X scoring: water + cyclohexane, each its own native ensemble
  python free_energy_calculator.py \
      --leg water=path/to/DOPC_3-12-8-12_S/water/ensemble.xyz \
      --leg cyclohexane=path/to/DOPC_3-12-8-12_S/cyclohexane/ensemble.xyz \
      --ref water --out fe_812S.csv

  # same run with ALPB, for the one-time comparison table
  python free_energy_calculator.py --method alpb \
      --leg water=.../water/ensemble.xyz --leg cyclohexane=.../cyclohexane/ensemble.xyz \
      --ref water --out fe_812S_alpb.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# --- constants --------------------------------------------------------------
HARTREE_TO_KCAL = 627.5094740631
R_KCAL = 1.987204258640832e-3   # kcal / mol / K
DEFAULT_T = 298.15              # K
_GSOLV_RE = re.compile(r"(?:CPCM-?X|solvation free energy|->\s*Gsolv)[^\-\d]*(-?\d+\.\d+)",
                       re.IGNORECASE)


# --- multi-frame XYZ --------------------------------------------------------
def parse_xyz_ensemble(path: Path):
    """Yield (index, comment, xyz_text) for each frame of a multi-frame .xyz."""
    lines = Path(path).read_text().splitlines()
    i, idx, frames = 0, 0, []
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            natoms = int(lines[i].split()[0])
        except (ValueError, IndexError):
            i += 1
            continue
        comment = lines[i + 1] if i + 1 < len(lines) else ""
        atom_lines = lines[i + 2 : i + 2 + natoms]
        block = "\n".join([str(natoms), comment, *atom_lines]) + "\n"
        frames.append((idx, comment.strip(), block))
        idx += 1
        i += 2 + natoms
    return frames


# --- one xtb single-point ---------------------------------------------------
def xtb_score(xyz_text, solvent, method, charge, xtb_bin, gfn, uhf, add_gsolv):
    """Run one xtb single-point; return (E_total_Eh, Gsolv_Eh_or_None, status).

    method='alpb'  : total energy from xtbout.json already includes solvation.
    method='cpcmx' : CPCM-X is a post-SCF correction. VERIFY on your build whether the
                     json total already includes it. Default: use json total as-is and
                     ALSO record the parsed Gsolv separately. If your build reports the
                     GAS total under --cpcmx, pass --cpcmx-add-gsolv to add Gsolv back.
    """
    with tempfile.TemporaryDirectory(prefix="xtb_") as td:
        (Path(td) / "conf.xyz").write_text(xyz_text)
        if method == "cpcmx":
            solv_flag = ["--cpcmx", solvent]
        else:
            solv_flag = ["--alpb", solvent]
        cmd = [xtb_bin, "conf.xyz", "--gfn", gfn, *solv_flag,
               "--chrg", str(charge), "--uhf", str(uhf), "--json"]
        try:
            proc = subprocess.run(cmd, cwd=td, capture_output=True, text=True,
                                  timeout=1800)
        except (subprocess.TimeoutExpired, FileNotFoundError) as exc:
            return None, None, type(exc).__name__

        e_total = None
        jpath = Path(td) / "xtbout.json"
        if jpath.exists():
            try:
                data = json.loads(jpath.read_text())
                for key in ("total energy", "total_energy", "energy"):
                    if key in data:
                        e_total = float(data[key])
                        break
            except (json.JSONDecodeError, ValueError, TypeError):
                pass
        if e_total is None:                       # stdout fallback
            for ln in proc.stdout.splitlines():
                if "total energy" in ln.lower():
                    for tok in ln.replace(":", " ").split():
                        try:
                            e_total = float(tok)
                            break
                        except ValueError:
                            continue

        gsolv = None
        m = _GSOLV_RE.search(proc.stdout)
        if m:
            try:
                gsolv = float(m.group(1))
            except ValueError:
                gsolv = None

        if e_total is None:
            return None, gsolv, f"parse-fail(rc={proc.returncode})"
        if method == "cpcmx" and add_gsolv and gsolv is not None:
            e_total = e_total + gsolv
        return e_total, gsolv, "ok"


# --- ensemble statistics ----------------------------------------------------
def ensemble_free_energy(energies_kcal, T=DEFAULT_T):
    e = [x for x in energies_kcal if x is not None]
    if not e:
        return None
    emin = min(e)
    z = sum(math.exp(-(x - emin) / (R_KCAL * T)) for x in e)
    return emin - R_KCAL * T * math.log(z)


def boltzmann_pops(energies_kcal, T=DEFAULT_T):
    e = [x if x is not None else float("inf") for x in energies_kcal]
    emin = min(e)
    w = [math.exp(-(x - emin) / (R_KCAL * T)) for x in e]
    s = sum(w)
    return [x / s for x in w] if s else [0.0] * len(w)


def score_leg(solvent, xyz_path, method, charge, xtb_bin, gfn, uhf, add_gsolv, jobs):
    """Score every conformer of one phase-specific ensemble. Returns (rows, G_ens_kcal)."""
    frames = parse_xyz_ensemble(xyz_path)
    if not frames:
        print(f"  !! no frames in {xyz_path}", file=sys.stderr)
        return [], None
    print(f"  [{solvent}] {xyz_path}: {len(frames)} conformers x {method}",
          file=sys.stderr)

    energies = [None] * len(frames)
    gsolvs = [None] * len(frames)

    def _work(item):
        idx, _c, block = item
        e, g, st = xtb_score(block, solvent, method, charge, xtb_bin, gfn, uhf, add_gsolv)
        return idx, e, g, st

    done = fail = 0
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        for fut in as_completed(ex.submit(_work, f) for f in frames):
            idx, e, g, st = fut.result()
            energies[idx], gsolvs[idx] = e, g
            done += 1
            fail += (e is None)
            if done % 50 == 0 or done == len(frames):
                print(f"    {done}/{len(frames)} ({fail} failed)", file=sys.stderr)

    kcal = [None if e is None else e * HARTREE_TO_KCAL for e in energies]
    pops = boltzmann_pops(kcal)
    emin = min(x for x in kcal if x is not None)
    rows = []
    for (idx, comment, _b) in frames:
        rows.append({
            "solvent": solvent, "conf": idx, "method": method,
            "E_Eh": energies[idx],
            "Gsolv_Eh": gsolvs[idx],
            "relE_kcal": None if kcal[idx] is None else kcal[idx] - emin,
            "pop": pops[idx],
        })
    return rows, ensemble_free_energy(kcal)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--leg", action="append", required=True, metavar="SOLVENT=PATH",
                    help="a phase: ALPB/CPCM-X solvent keyword = path to that phase's "
                         "native ensemble.xyz. Repeat for each phase (e.g. water, "
                         "cyclohexane).")
    ap.add_argument("--ref", default="water",
                    help="reference solvent for ΔG_transfer (default water)")
    ap.add_argument("--method", choices=["cpcmx", "alpb"], default="cpcmx")
    ap.add_argument("--compare", action="store_true",
                    help="score every leg with BOTH cpcmx and alpb and print ΔG_transfer "
                         "side-by-side (the one-time methods check); overrides --method")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--uhf", type=int, default=0)
    ap.add_argument("--gfn", default="2")
    ap.add_argument("--cpcmx-add-gsolv", action="store_true",
                    help="add parsed Gsolv to the json total (only if your xtb build "
                         "reports the GAS total under --cpcmx; verify first)")
    ap.add_argument("--xtb", default="xtb")
    ap.add_argument("--jobs", type=int, default=1, help="parallel xtb workers")
    ap.add_argument("--out", type=Path, default=Path("free_energy.csv"),
                    help="per-conformer CSV (summary -> <out>.summary.csv)")
    args = ap.parse_args(argv)

    if shutil.which(args.xtb) is None:
        sys.exit(f"error: '{args.xtb}' not on PATH — run on the HPC (xtb is Linux-only)")

    legs = {}
    for spec in args.leg:
        if "=" not in spec:
            sys.exit(f"error: --leg must be SOLVENT=PATH, got '{spec}'")
        solv, path = spec.split("=", 1)
        legs[solv.strip()] = Path(path.strip())
    if args.ref not in legs:
        sys.exit(f"error: --ref '{args.ref}' has no matching --leg")

    methods = ["cpcmx", "alpb"] if args.compare else [args.method]

    all_rows, G = [], {}     # G[(method, solvent)] = ensemble free energy (kcal)
    for method in methods:
        for solv, xp in legs.items():
            rows, g_ens = score_leg(solv, xp, method, args.charge, args.xtb,
                                    args.gfn, args.uhf, args.cpcmx_add_gsolv, args.jobs)
            all_rows.extend(rows)
            G[(method, solv)] = g_ens

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()))
        w.writeheader()
        w.writerows(all_rows)
    print(f"wrote {len(all_rows)} conformer rows -> {args.out}", file=sys.stderr)

    summaries = []
    for method in methods:
        s = {"method": method, "ref_solvent": args.ref}
        for solv in legs:
            s[f"Gens_{solv}_kcal"] = G[(method, solv)]
            s[f"n_{solv}"] = sum(1 for r in all_rows
                                 if r["solvent"] == solv and r["method"] == method)
            if solv != args.ref:
                gr, gs = G[(method, args.ref)], G[(method, solv)]
                s[f"dGtransfer_{args.ref}->{solv}_kcal"] = (
                    None if gr is None or gs is None else gs - gr)
        summaries.append(s)

    keys = []
    for s in summaries:
        keys += [k for k in s if k not in keys]
    spath = args.out.with_suffix(".summary.csv")
    with open(spath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(summaries)
    print(f"wrote summary -> {spath}", file=sys.stderr)

    for s in summaries:
        print(f"[{s['method']}]", file=sys.stderr)
        for k, v in s.items():
            if k.startswith("dGtransfer") and v is not None:
                print(f"    {k} = {v:+.2f} kcal/mol", file=sys.stderr)
    if args.compare:
        print("  (compare: run this for R and S — CPCM-X 'wins' if it keeps the same "
              "R-vs-S ΔG_transfer ordering as ALPB but separates them more cleanly)",
              file=sys.stderr)


if __name__ == "__main__":
    main()

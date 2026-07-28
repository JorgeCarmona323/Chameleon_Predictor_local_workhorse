#!/usr/bin/env python3
"""free_energy_orca_smd — ΔG_transfer from CREST ensembles via ORCA DFT + SMD single-points.

The higher-tier DFT counterpart to free_energy_calculator.py (xTB CPCM-X/ALPB). It scores
the SAME phase-specific CREST ensembles, but each single-point is an ORCA DFT calculation
with the SMD continuum model instead of an xTB one — so the resulting ΔG_transfer is directly
comparable to the CPCM-X number and tells us how far off the cheap xTB solvation is.

It deliberately REUSES free_energy_calculator for everything except the single-point engine
(ensemble parsing, --ewin trim, Boltzmann G_ens, populations), so the two methods differ ONLY
in the per-conformer energy — an apples-to-apples methods comparison by construction:

  * per-conformer DFT+SMD energy + Boltzmann weight (within each phase)
  * per-phase ensemble free energy  G_ens = -RT ln SUM_i exp(-E_i / RT)
  * ΔG_transfer(ref -> phase) = G_ens(phase) - G_ens(ref)   [ref default = water]

DFT single-points are ~minutes each (not seconds like xTB), so --top caps each phase to its N
lowest-CREST-energy conformers AFTER the --ewin trim (default 12). Geometry is NOT re-optimized
(CREST already relaxed each phase in ALPB); SMD is a single-point solvation model on top.

Level of theory default is r2SCAN-3c (def2-mTZVPP + gCP + D4) — the efficient composite that is
the standard sweet spot for SMD single-points. ORCA reports the SMD CDS (cavity-dispersion-
solvent-structure) term inside FINAL SINGLE POINT ENERGY, which is the value used; the CDS is
also parsed out and logged per conformer so the first job can be checked (nonzero, solvent-
dependent) to confirm SMD actually engaged.

Runs on the HPC (ORCA is Linux and set up there). Pure stdlib — no rdkit/pandas.

Example
-------
  python scripts/free_energy_orca_smd.py \
      --orca "$HOME/orca_6.1.1/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca" \
      --leg water=.../HexPep/water/ensemble.xyz \
      --leg hexane=.../HexPep/hexane/ensemble.xyz \
      --ref water --ewin 5 --top 12 --jobs 12 \
      --workdir results/free_energy/hexpep_smd_orca \
      --out results/free_energy/hexpep_smd.csv
"""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# reuse the xTB calculator's ensemble/statistics helpers so the two methods share EXACTLY
# the same parsing, energy window, and Boltzmann math (only the single-point engine differs)
sys.path.insert(0, str(Path(__file__).resolve().parent))
import free_energy_calculator as fe  # noqa: E402

# our --leg solvent keys  ->  ORCA SMD solvent names (SMD solvent list in the ORCA manual)
ORCA_SMD = {
    "water": "water", "h2o": "water",
    "hexane": "n-hexane", "n-hexane": "n-hexane",
    "chcl3": "chloroform", "chloroform": "chloroform",
    "octanol": "1-octanol",
}

_FINAL_E_RE = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
# ORCA prints the SMD non-electrostatic term as e.g. "SMD-CDS (Gcds) ... -0.00123 Eh"
_CDS_RE = re.compile(r"SMD[- ]CDS[^\n]*?(-?\d+\.\d+)\s*(?:Eh|a\.u\.|hartree)", re.IGNORECASE)


def orca_smd_score(xyz_text, smd_solvent, charge, mult, level, orca_bin, workdir):
    """Run one ORCA DFT+SMD single-point in `workdir`; return (E_total_Eh, Gcds_kcal, status).

    E_total is ORCA's FINAL SINGLE POINT ENERGY, which in ORCA 5/6 already INCLUDES the SMD CDS
    correction (SCF-with-CPCM/SMD-electrostatics + Gcds). The .inp/.out are kept in `workdir` for
    traceability."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / "conf.xyz").write_text(xyz_text)
    inp = (
        f"! {level}\n"
        f"%maxcore 3000\n"
        f"%pal nprocs 1 end\n"          # serial per job; we parallelize ACROSS conformers
        f"%cpcm\n"
        f"   smd true\n"
        f'   SMDsolvent "{smd_solvent}"\n'
        f"end\n"
        f"* xyzfile {charge} {mult} conf.xyz *\n"
    )
    (workdir / "conf.inp").write_text(inp)
    outpath = workdir / "conf.out"
    try:
        with open(outpath, "w") as fo:
            proc = subprocess.run([str(orca_bin), "conf.inp"], cwd=workdir,
                                  stdout=fo, stderr=subprocess.STDOUT, timeout=7200)
    except (subprocess.TimeoutExpired, FileNotFoundError) as exc:
        return None, None, type(exc).__name__

    text = outpath.read_text(errors="ignore")
    finals = _FINAL_E_RE.findall(text)
    e_total = float(finals[-1]) if finals else None

    cds = None
    mc = _CDS_RE.search(text)
    if mc:
        try:
            cds = float(mc.group(1)) * fe.HARTREE_TO_KCAL
        except ValueError:
            cds = None

    if e_total is None:
        return None, cds, f"parse-fail(rc={proc.returncode})"
    return e_total, cds, "ok"


def select_top(frames, top):
    """Keep the `top` lowest-CREST-energy frames (comment energy = Hartree). None-energy frames
    sort last. If top is falsy or >= len(frames), return frames unchanged."""
    if not top or len(frames) <= top:
        return frames

    def keyfn(fr):
        e = fe._frame_energy_hartree(fr[1])   # fr = (idx, comment, block); comment carries E
        return e if e is not None else float("inf")

    return sorted(frames, key=keyfn)[:top]


def score_leg_smd(solvent, xyz_path, smd_solvent, charge, mult, level, orca_bin, jobs,
                  workroot, ewin, top):
    """DFT+SMD-score one phase-specific ensemble. Returns (rows, G_ens_kcal)."""
    frames = fe.parse_xyz_ensemble(xyz_path)
    if not frames:
        print(f"  !! no frames in {xyz_path}", file=sys.stderr)
        return [], None
    n_total = len(frames)
    if ewin is not None:
        frames = fe.apply_energy_window(frames, ewin)
    frames = select_top(frames, top)
    print(f"  [{solvent}] {n_total} -> {len(frames)} conformers x ORCA-SMD({smd_solvent}) "
          f"@ {level}", file=sys.stderr)
    if not frames:
        return [], None

    n = len(frames)
    energies = [None] * n
    cdss = [None] * n
    stats = [None] * n

    def _work(item):
        pos, (idx, _c, block) = item
        wd = Path(workroot) / f"conf{idx:03d}_{solvent}"
        e, c, st = orca_smd_score(block, smd_solvent, charge, mult, level, orca_bin, wd)
        return pos, e, c, st

    done = fail = 0
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        for fut in as_completed(ex.submit(_work, it) for it in enumerate(frames)):
            pos, e, c, st = fut.result()
            energies[pos], cdss[pos], stats[pos] = e, c, st
            done += 1
            fail += (e is None)
            print(f"    {done}/{n} ({fail} failed)  [{stats[pos]}]", file=sys.stderr)

    kcal = [None if e is None else e * fe.HARTREE_TO_KCAL for e in energies]
    pops = fe.boltzmann_pops(kcal)
    valid = [x for x in kcal if x is not None]
    emin = min(valid) if valid else None
    rows = []
    for pos, (idx, _c, _b) in enumerate(frames):
        rows.append({
            "solvent": solvent, "conf": idx, "method": f"smd/{level}",
            "E_Eh": energies[pos], "Gcds_kcal": cdss[pos], "status": stats[pos],
            "relE_kcal": None if (kcal[pos] is None or emin is None) else kcal[pos] - emin,
            "pop": pops[pos],
        })
    return rows, fe.ensemble_free_energy(kcal)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--leg", action="append", required=True, metavar="SOLVENT=PATH",
                    help="a phase: solvent key (water/hexane/chcl3) = path to that phase's "
                         "native ensemble.xyz. Repeat per phase.")
    ap.add_argument("--ref", default="water", help="reference solvent for ΔG_transfer")
    ap.add_argument("--orca", required=True, help="absolute path to the ORCA binary")
    ap.add_argument("--level", default="r2SCAN-3c",
                    help="ORCA level of theory / simple-input line (default r2SCAN-3c)")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1, help="spin multiplicity (default 1)")
    ap.add_argument("--ewin", type=float, default=None, metavar="KCAL",
                    help="pre-filter to conformers within KCAL of the lowest CREST energy "
                         "before selecting --top (recommended 5)")
    ap.add_argument("--top", type=int, default=12, metavar="N",
                    help="cap each phase to its N lowest-CREST-energy conformers after --ewin "
                         "(DFT is expensive; default 12; 0 = all)")
    ap.add_argument("--jobs", type=int, default=1, help="concurrent ORCA workers (each serial)")
    ap.add_argument("--workdir", default=None,
                    help="keep ORCA .inp/.out here (per-conformer subdirs). Default: a temp dir "
                         "that is deleted on exit.")
    ap.add_argument("--out", type=Path, default=Path("free_energy_smd.csv"),
                    help="per-conformer CSV (summary -> <out>.summary.csv)")
    args = ap.parse_args(argv)

    if not Path(args.orca).exists():
        sys.exit(f"error: --orca '{args.orca}' does not exist")

    legs = {}
    for spec in args.leg:
        if "=" not in spec:
            sys.exit(f"error: --leg must be SOLVENT=PATH, got '{spec}'")
        solv, path = spec.split("=", 1)
        legs[solv.strip()] = Path(path.strip())
    if args.ref not in legs:
        sys.exit(f"error: --ref '{args.ref}' has no matching --leg")

    tmp = None
    if args.workdir:
        workroot = Path(args.workdir)
        workroot.mkdir(parents=True, exist_ok=True)
    else:
        tmp = tempfile.TemporaryDirectory(prefix="orca_smd_")
        workroot = Path(tmp.name)

    all_rows, G = [], {}
    for solv, xp in legs.items():
        smd_name = ORCA_SMD.get(solv.lower())
        if smd_name is None:
            print(f"  !! '{solv}' not in ORCA SMD map — passing it through verbatim; "
                  f"verify it is a valid SMDsolvent name", file=sys.stderr)
            smd_name = solv
        rows, g_ens = score_leg_smd(solv, xp, smd_name, args.charge, args.mult, args.level,
                                    args.orca, args.jobs, workroot, args.ewin, args.top)
        all_rows.extend(rows)
        G[solv] = g_ens

    if not all_rows:
        sys.exit("error: no conformers scored — check the ensemble paths and ORCA setup")

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()))
        w.writeheader()
        w.writerows(all_rows)
    print(f"wrote {len(all_rows)} conformer rows -> {args.out}", file=sys.stderr)

    summary = {"method": f"smd/{args.level}", "ref_solvent": args.ref}
    for solv in legs:
        summary[f"Gens_{solv}_kcal"] = G[solv]
        summary[f"n_{solv}"] = sum(1 for r in all_rows if r["solvent"] == solv)
        if solv != args.ref:
            gr, gs = G[args.ref], G[solv]
            summary[f"dGtransfer_{args.ref}->{solv}_kcal"] = (
                None if gr is None or gs is None else gs - gr)

    spath = args.out.with_suffix(".summary.csv")
    with open(spath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary.keys()))
        w.writeheader()
        w.writerow(summary)
    print(f"wrote summary -> {spath}", file=sys.stderr)
    for k, v in summary.items():
        if k.startswith("dGtransfer") and v is not None:
            print(f"    {k} = {v:+.2f} kcal/mol", file=sys.stderr)

    if tmp is not None:
        tmp.cleanup()


if __name__ == "__main__":
    main()

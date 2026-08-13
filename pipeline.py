#!/usr/bin/env python3
"""pipeline.py -- one command: SMILES -> conformers -> energies -> 3D descriptors.

The single entry point for the 3D-descriptor pipeline. Chains the three stages and hands each
stage's output to the next:

  Stage 1  CREST conformer generation      scripts/crest_engine.generate_conformers()
  Stage 2  CPCM-X dG_transfer energies      scripts/free_energy_calculator.py
  Stage 3  3D physics descriptors           scripts/ensemble_descriptors.py
           (Boltzmann-weighted by the stage-2 SOLVATED populations -- geometry from CREST,
            populations from the solvated single-points)

Usage
-----
  # full run from a SMILES
  python pipeline.py --smiles "<SMILES>" --name Foo

  # full run from the reference registry (scripts/crest_v3.2.py)
  python pipeline.py --compound 24

  # skip generation: score + describe an EXISTING ensemble directory (water/ + apolar/ inside)
  python pipeline.py --run-dir results/runs/run_..._WhC3 --name WhC3

Resumes by default: a stage whose output already exists is skipped (use --force to redo).
Descriptors are computed for every apolar leg present (chloroform and/or hexane) -> both
water->chloroform and water->hexane deltas.

Outputs land in the compound's run directory:
  <run>/water|chloroform|hexane/ensemble.{xyz,sdf}   (stage 1)
  <run>/free_energy_cpcmx.csv[.summary.csv]           (stage 2, with coordinate provenance)
  <run>/descriptors_<apolar>.csv                      (stage 3)
  <run>/pipeline_manifest.json                        (what ran, where the inputs came from)

Environment: needs rdkit + xtb + crest + numpy + pandas on the same PATH. Source scripts/env.sh
first (installs/points at the CPCM-X-enabled xtb). Stages 1-2 are Linux-only (xtb/crest).
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

# possible folder names for the chloroform leg (pipeline writes 'chloroform'; older runs used these)
CHCL3_FOLDERS = ("chloroform", "chcl3", "mem")
BAR = "=" * 72


def load_registry():
    """REFERENCE_COMPOUNDS from crest_v3.2.py (dotted filename -> importlib)."""
    spec = importlib.util.spec_from_file_location("crest_v32", SCRIPTS / "crest_v3.2.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.REFERENCE_COMPOUNDS


def leg_path(base: Path, *names):
    """First existing <base>/<name>/ensemble.xyz among names, else None."""
    for n in names:
        p = base / n / "ensemble.xyz"
        if p.exists():
            return p
    return None


def resolve_charge(charge, smiles, work: Path) -> int:
    """Explicit --charge wins; else formal charge from the SMILES (rdkit); else read a
    solvent metadata.json; else 0."""
    if charge is not None:
        return charge
    if smiles:
        try:
            from rdkit import Chem
            m = Chem.MolFromSmiles(smiles)
            if m is not None:
                return Chem.GetFormalCharge(m)
        except Exception:
            pass
    for s in (*CHCL3_FOLDERS, "water", "hexane"):
        meta = work / s / "metadata.json"
        if meta.exists():
            try:
                d = json.loads(meta.read_text())
                if "charge" in d:
                    return int(d["charge"])
            except Exception:
                pass
    return 0


def run_step(cmd, desc):
    print(f"\n{BAR}\n>> {desc}\n{BAR}\n  $ {' '.join(str(c) for c in cmd)}", flush=True)
    if subprocess.run(cmd, cwd=ROOT).returncode != 0:
        sys.exit(f"ERROR: {desc} failed")
    print(f"OK: {desc}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_argument_group("input (one of)")
    src.add_argument("--smiles", help="molecule SMILES (with --name)")
    src.add_argument("--compound", type=int, help="index into scripts/crest_v3.2.py REFERENCE_COMPOUNDS")
    src.add_argument("--run-dir", help="skip stage 1; use this existing ensemble dir (water/ + apolar/)")
    ap.add_argument("--name", help="compound label (default: registry short name, or 'molecule')")
    ap.add_argument("--charge", type=int, default=None, help="formal charge (default: auto)")
    ap.add_argument("--gfn", default="2", choices=["2", "1", "0", "ff"],
                    help="stage-1 CREST level of theory: 2/1/0 = GFN-n, ff = GFN-FF (fast, for "
                         "large flexible macrocycles). Stage-2 scoring stays GFN2/CPCM-X regardless.")
    ap.add_argument("--method", default="cpcmx", choices=["cpcmx", "alpb"], help="stage-2 solvation model")
    ap.add_argument("--ewin", type=float, default=8.0, help="stage-2 energy-window pre-trim (kcal/mol)")
    ap.add_argument("--threads", type=int, default=None, help="CPU threads for CREST/xtb (default: all)")
    ap.add_argument("--max-confs", type=int, default=None, help="cap conformers kept per solvent")
    ap.add_argument("--outdir", default="results/pipeline", help="where stage-1 creates the run dir")
    ap.add_argument("--force", action="store_true", help="re-run stages even if outputs exist")
    args = ap.parse_args(argv)

    # ---- resolve input -------------------------------------------------------
    smiles, name = args.smiles, args.name
    if args.compound is not None:
        cpd = load_registry()[args.compound]
        smiles = smiles or cpd["smiles"]
        name = name or cpd.get("short") or cpd.get("name")
    name = name or "molecule"

    # ---- stage 1: conformers -------------------------------------------------
    if args.run_dir:
        work = Path(args.run_dir).resolve()
        print(f"[stage 1] skipped -- using existing ensembles at {work}")
    else:
        if not smiles:
            sys.exit("error: provide --smiles (with --name) or --compound, or --run-dir to skip generation")
        import crest_engine as ce
        ce.GFN_METHOD = args.gfn      # 2/1/0/ff — level of theory for CREST + xtb pre-opt
        print(f"[stage 1] CREST generation for {name} (GFN{'-FF' if args.gfn=='ff' else args.gfn}) -> {args.outdir}")
        res = ce.generate_conformers(smiles, name=name, outdir=args.outdir,
                                     charge=args.charge, n_threads=args.threads,
                                     max_confs=args.max_confs)
        work = Path(res["work_dir"]).resolve()
        if not res.get("ok"):
            failed = res.get("result", {}).get("failed_solvents", "?")
            print(f"WARNING: stage 1 -- some solvent legs failed ({failed}); continuing with what generated")

    water = work / "water" / "ensemble.xyz"
    chcl3 = leg_path(work, *CHCL3_FOLDERS)
    hexane = work / "hexane" / "ensemble.xyz"
    if not water.exists():
        sys.exit(f"error: no water ensemble at {water} -- stage 1 incomplete?")
    charge = resolve_charge(args.charge, smiles, work)

    # ---- stage 2: energies ---------------------------------------------------
    fe_csv = work / "free_energy_cpcmx.csv"
    if fe_csv.exists() and not args.force:
        print(f"[stage 2] skipped -- {fe_csv} exists (--force to redo)")
    else:
        legs = ["--leg", f"water={water}"]
        if chcl3:
            legs += ["--leg", f"chloroform={chcl3}"]
        if hexane.exists():
            legs += ["--leg", f"hexane={hexane}"]
        run_step([sys.executable, str(SCRIPTS / "free_energy_calculator.py"),
                  "--method", args.method, "--ewin", str(args.ewin), "--ref", "water",
                  "--charge", str(charge), "--jobs", str(args.threads or 1),
                  *legs, "--out", str(fe_csv)],
                 f"Stage 2 -- {args.method.upper()} dG_transfer (charge {charge})")

    # ---- stage 3: descriptors (one run per apolar leg present) ----------------
    apolar_folders = []
    if chcl3:
        apolar_folders.append(chcl3.parent.name)   # 'chloroform' / 'chcl3' / 'mem'
    if hexane.exists():
        apolar_folders.append("hexane")
    if not apolar_folders:
        sys.exit("error: no apolar leg (chloroform or hexane) -- cannot compute delta-PSA descriptors")

    descriptors = {}
    for folder in apolar_folders:
        out = work / f"descriptors_{folder}.csv"
        if out.exists() and not args.force:
            print(f"[stage 3/{folder}] skipped -- {out} exists")
        else:
            run_step([sys.executable, str(SCRIPTS / "ensemble_descriptors.py"),
                      "--run-dir", str(work), "--apolar", folder, "--name", name,
                      "--energies-csv", str(fe_csv), "-o", str(out)],
                     f"Stage 3 -- 3D descriptors (water vs {folder})")
        descriptors[folder] = str(out)

    # ---- manifest (provenance) ----------------------------------------------
    ensembles = {s: str(work / s / "ensemble.xyz")
                 for s in ("water", "hexane") if (work / s / "ensemble.xyz").exists()}
    if chcl3:
        ensembles[chcl3.parent.name] = str(chcl3)
    manifest = {
        "name": name, "smiles": smiles, "charge": charge,
        "run_dir": str(work), "completed": datetime.now().isoformat(timespec="seconds"),
        "stage1_ensembles": ensembles,
        "stage2_energies": str(fe_csv),
        "stage2_summary": str(fe_csv.with_suffix(".summary.csv")),
        "stage3_descriptors": descriptors,
    }
    (work / "pipeline_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"\n{BAR}\nDONE -- pipeline complete -> {work}\n"
          f"  energies:    {fe_csv.with_suffix('.summary.csv').name}\n"
          f"  descriptors: {', '.join(Path(p).name for p in descriptors.values())}\n"
          f"  manifest:    pipeline_manifest.json\n{BAR}")


if __name__ == "__main__":
    main()

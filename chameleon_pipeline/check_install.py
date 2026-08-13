#!/usr/bin/env python3
"""check_install.py -- verify this machine can run the pipeline.

Checks, in order:
  1. Python packages:  rdkit, numpy, pandas   (import test)
  2. Binaries on PATH: xtb, crest              (the conformer + energy engines)
  3. xtb version
  4. CPCM-X support:   runs a 1-molecule `xtb --gfn 2 --cpcmx water` -- the #1 gotcha, because
                       the conda/CREST xtb is built WITHOUT CPCM-X. This is what actually decides
                       whether the ENERGY stage works.

Exit code 0 = ready to run; nonzero = something's missing (each line says what and how to fix).
Run after `conda activate chameleon_pipeline`, `bash setup_xtb.sh`, and `source env.sh`.
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

OK, WARN, BAD = "OK  ", "WARN", "FAIL"
problems = []


def line(status, msg):
    print(f"  [{status}] {msg}")


print("Chameleon pipeline -- install check\n" + "-" * 44)

# 1. python packages -----------------------------------------------------------
print("Python packages:")
for pkg in ("rdkit", "numpy", "pandas"):
    try:
        __import__(pkg)
        line(OK, f"{pkg} importable")
    except Exception as e:
        line(BAD, f"{pkg} MISSING ({e}) -- `conda env create -f environment.yml`")
        problems.append(pkg)

# 2. binaries ------------------------------------------------------------------
print("\nBinaries on PATH:")
xtb = shutil.which("xtb")
crest = shutil.which("crest")
line(OK if xtb else BAD, f"xtb   -> {xtb or 'NOT FOUND -- run setup_xtb.sh then `source env.sh`'}")
line(OK if crest else BAD, f"crest -> {crest or 'NOT FOUND -- `conda install -c conda-forge crest`'}")
if not xtb:
    problems.append("xtb")
if not crest:
    problems.append("crest")

# 3. xtb version ---------------------------------------------------------------
if xtb:
    try:
        v = subprocess.run([xtb, "--version"], capture_output=True, text=True, timeout=60)
        ver = next((ln.strip() for ln in (v.stdout + v.stderr).splitlines() if "version" in ln.lower()), "?")
        print(f"\nxtb version:\n  [{OK}] {ver}")
    except Exception as e:
        line(WARN, f"could not read xtb version ({e})")

# 4. CPCM-X support (the decisive check for the energy stage) -------------------
print("\nCPCM-X solvation (energy stage):")
if not xtb:
    line(BAD, "skipped -- xtb not on PATH")
    problems.append("cpcmx")
else:
    with tempfile.TemporaryDirectory() as td:
        (Path(td) / "h2o.xyz").write_text(
            "3\n\nO 0.0 0.0 0.0\nH 0.757 0.586 0.0\nH -0.757 0.586 0.0\n")
        try:
            r = subprocess.run([xtb, "h2o.xyz", "--gfn", "2", "--cpcmx", "water", "--chrg", "0"],
                               cwd=td, capture_output=True, text=True, timeout=180)
            out = r.stdout + r.stderr
            if "not included" in out.lower() or "cpcm-x library was not" in out.lower():
                line(BAD, "this xtb was built WITHOUT CPCM-X -- run setup_xtb.sh to install the "
                          "grimme-lab release, then `source env.sh` so it is first on PATH")
                problems.append("cpcmx")
            elif r.returncode == 0:
                line(OK, "'--cpcmx water' ran and finished normally -- CPCM-X is available")
            else:
                line(WARN, f"'--cpcmx water' exited {r.returncode} (check XTBPATH; `source env.sh`)")
        except Exception as e:
            line(BAD, f"CPCM-X test failed to run ({e})")
            problems.append("cpcmx")

# verdict ----------------------------------------------------------------------
print("-" * 44)
if problems:
    print(f"NOT READY -- fix: {', '.join(sorted(set(problems)))}")
    print("Descriptor-only use (stage 3) still works with just rdkit/numpy/pandas if you already "
          "have ensembles; conformer generation + energies need xtb/crest/CPCM-X.")
    sys.exit(1)
print("READY -- all stages can run. Try:  python pipeline.py --smiles \"...\" --name test")
sys.exit(0)

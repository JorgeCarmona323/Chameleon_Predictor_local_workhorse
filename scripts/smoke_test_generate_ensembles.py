# env: chameleon-calc
"""
smoke_test_generate_ensembles.py
--------------------------------
Fast validation of the registry-free direct-SMILES front end added to
`crest_conformers.py` (generate_ensembles / check_binaries / _safe_short).

These checks do NOT run CREST/xTB — they exercise the input-validation and
plumbing layers only, so they pass on any machine with RDKit (no external
binaries needed). The full end-to-end generation is validated separately on a
host that has xtb + crest installed.

Run:  python scripts/smoke_test_generate_ensembles.py
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from crest_conformers import generate_ensembles, check_binaries, _safe_short

SIMPLE = "CCO"  # ethanol
# macrocycle-like: a real DOPC hit (thioether xylene macrocycle)
MACRO = ("C#CCCC(=O)N[C@H]1CSCc2ccccc2CSC[C@@H](C(N)=O)NC(=O)[C@H](c2cccs2)"
         "NC(=O)[C@H](CO)NC(=O)[C@@H]2CCN2C(=O)[C@H](CO)NC1=O")

n_pass = n_fail = 0


def check(desc, cond):
    global n_pass, n_fail
    print(f"  [{'PASS' if cond else 'FAIL'}] {desc}")
    n_pass += cond
    n_fail += (not cond)


print("1. SMILES validation")
try:
    generate_ensembles("not_a_smiles!!!", name="bad", check_binaries_first=False)
    check("bad SMILES raises ValueError", False)
except ValueError:
    check("bad SMILES raises ValueError", True)
except Exception as e:
    check(f"bad SMILES raises ValueError (got {type(e).__name__})", False)

try:
    generate_ensembles("", name="empty", check_binaries_first=False)
    check("empty SMILES raises ValueError", False)
except ValueError:
    check("empty SMILES raises ValueError", True)

print("\n2. RDKit accepts valid inputs (simple + macrocycle)")
from rdkit import Chem
check("simple SMILES parses", Chem.MolFromSmiles(SIMPLE) is not None)
check("macrocycle SMILES parses", Chem.MolFromSmiles(MACRO) is not None)

print("\n3. Filename sanitization")
check("spaces/slashes sanitized", _safe_short("DOPC 3/12 R") == "DOPC_3_12_R")
check("empty name falls back", _safe_short("") == "molecule")

print("\n4. Binary check")
have_bins = shutil.which("xtb") and shutil.which("crest")
if have_bins:
    check("xtb + crest found on PATH", True)
else:
    try:
        check_binaries()
        check("missing binaries raise RuntimeError", False)
    except RuntimeError as e:
        check("missing binaries raise RuntimeError", True)
        check("error names the missing binary", "xtb" in str(e) or "crest" in str(e))
    # and generate_ensembles should fail fast with the same clear error
    try:
        generate_ensembles(SIMPLE, name="nobins")
        check("generate_ensembles fails fast without binaries", False)
    except RuntimeError:
        check("generate_ensembles fails fast without binaries", True)

print(f"\n{'='*40}\n{n_pass} passed, {n_fail} failed")
sys.exit(1 if n_fail else 0)

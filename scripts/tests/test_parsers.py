#!/usr/bin/env python3
"""test_parsers.py — functional tests for the ensemble/energy parsing layer.

Covers the plumbing that connects the stages, on tiny synthetic fixtures (fast, no xtb,
no real ensembles):

  * free_energy_calculator.parse_xyz_ensemble   — multi-frame XYZ -> frames
  * free_energy_calculator.apply_energy_window  — the --ewin pre-filter
  * ensemble_descriptors._energies_from_xyz     — energies out of CREST comment lines
  * ensemble_descriptors._pops_from_energy_csv  — stage-2 populations scattered back by
    conformer index (the --ewin subset case: unscored conformers must get weight 0)

Picked up automatically by scripts/tests/run_all.sh. Exits non-zero on any failure.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SCRIPTS))

_passed = 0
_failed = 0


def check(desc, got, want):
    global _passed, _failed
    if got == want:
        print(f"  PASS  {desc}")
        _passed += 1
    else:
        print(f"  FAIL  {desc}\n          got:  {got!r}\n          want: {want!r}")
        _failed += 1


def approx(desc, got, want, tol=1e-6):
    global _passed, _failed
    try:
        ok = len(got) == len(want) and all(abs(a - b) <= tol for a, b in zip(got, want))
    except TypeError:
        ok = abs(got - want) <= tol
    if ok:
        print(f"  PASS  {desc}")
        _passed += 1
    else:
        print(f"  FAIL  {desc}\n          got:  {got!r}\n          want: {want!r}")
        _failed += 1


# a 2-frame, 2-atom ensemble; CREST writes the energy as the comment line
XYZ = """  2
       -1.00000000
 C   0.0 0.0 0.0
 H   0.0 0.0 1.0
  2
       -0.90000000
 C   0.0 0.0 0.0
 H   0.0 0.0 1.1
"""


def test_free_energy_calculator():
    print("free_energy_calculator")
    import free_energy_calculator as fe

    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "ensemble.xyz"
        p.write_text(XYZ)
        frames = fe.parse_xyz_ensemble(p)
        check("parse_xyz_ensemble finds 2 frames", len(frames), 2)
        check("frame index preserved", [f[0] for f in frames], [0, 1])

        # 0.1 Eh apart == ~62.75 kcal/mol, so a 1 kcal window keeps only the lowest
        kept = fe.apply_energy_window(frames, 1.0)
        check("ewin=1 kcal keeps only the lowest frame", len(kept), 1)
        check("ewin keeps the ORIGINAL conformer index", kept[0][0], 0)

        kept_all = fe.apply_energy_window(frames, 100.0)
        check("ewin=100 kcal keeps both", len(kept_all), 2)

        approx("ensemble_free_energy of one conformer == its energy",
               fe.ensemble_free_energy([0.0]), 0.0)


def test_ensemble_descriptors():
    print("ensemble_descriptors")
    try:
        import ensemble_descriptors as ed
    except ImportError as exc:          # rdkit/pandas missing -> not a parser failure
        print(f"  SKIP  cannot import ensemble_descriptors ({exc})")
        return

    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "ensemble.xyz"
        p.write_text(XYZ)
        approx("_energies_from_xyz reads both comment-line energies",
               ed._energies_from_xyz(p), [-1.0, -0.9])

        # stage-2 CSV: only conformers 0 and 2 survived --ewin, out of 4 total
        csv = Path(td) / "fe.csv"
        csv.write_text(
            "solvent,conf,method,E_Eh,Gsolv_Eh,relE_kcal,pop\n"
            "water,0,cpcmx,-1.0,-0.01,0.0,0.7\n"
            "water,2,cpcmx,-0.99,-0.01,0.6,0.3\n"
            "hexane,0,cpcmx,-1.1,-0.02,0.0,1.0\n"
        )
        w = ed._pops_from_energy_csv(csv, "water", 4)
        approx("pops scattered back by conf index; unscored get 0",
               list(w), [0.7, 0.0, 0.3, 0.0])

        wh = ed._pops_from_energy_csv(csv, "hexane", 4)
        approx("solvent filter selects the right leg", list(wh), [1.0, 0.0, 0.0, 0.0])

        check("unknown solvent -> None", ed._pops_from_energy_csv(csv, "chcl3", 4), None)


if __name__ == "__main__":
    test_free_energy_calculator()
    test_ensemble_descriptors()
    print(f"\n{_passed} passed, {_failed} failed")
    sys.exit(1 if _failed else 0)

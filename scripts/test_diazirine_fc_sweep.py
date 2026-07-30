#!/usr/bin/env python3
"""test_diazirine_fc_sweep.py — find the force constant that pins the diazirine N=N at 1.23 A.

Cheap single-conformer sweep (seconds per point, NO CREST). GFN2-xTB has a spurious minimum that
stretches the diazirine N=N toward ~1.43 A; our harmonic distance restraint (crest_engine, default
force constant 0.25 Eh/Bohr^2) is too weak to win that force balance for strained macrocycles, so it
settles at an intermediate ~1.318 A (FAIL) instead of 1.23. This script embeds a diazirine compound
with RDKit (good N=N ~1.24), then runs GFN2 constrained `--opt` at a ladder of force constants and
reports the resulting N=N — so we pick the force constant ONCE and regenerate ONCE, instead of
burning hours of CREST per guess.

It also measures the terminal alkyne C#C (a SECOND GFN2-distorted group that is currently NOT
constrained), so we can see whether it needs its own constraint for the pendant to be uniformly
correct (uniformity across legs is what makes the diazirine energy cancel in dG_transfer).

Run on the HPC (xtb is Linux-only), e.g. inside the crest env which has both rdkit and xtb:
  conda activate chameleon_crest212
  python scripts/test_diazirine_fc_sweep.py                       # compound 24, water, default ladder
  python scripts/test_diazirine_fc_sweep.py --compound 24 --solvent water --fcs 0 0.25 0.5 1.0 2.0 5.0
  python scripts/test_diazirine_fc_sweep.py --smiles "C#CC...N=N..."   # arbitrary diazirine
"""
from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# reuse the pipeline's constraint writer + diazirine detection so the test matches production
sys.path.insert(0, str(Path(__file__).resolve().parent))
import crest_engine as ce  # noqa: E402

ALKYNE_SMARTS = "[CX2]#[CX2]"
# same verdict bands as verify_diazirine_integrity.py
NN_PASS = (1.22, 1.25)
NN_WATCH_HI = 1.30


def load_reference_compounds():
    """Import REFERENCE_COMPOUNDS from crest_v3.2.py (dotted filename -> importlib)."""
    path = Path(__file__).resolve().parent / "crest_v3.2.py"
    spec = importlib.util.spec_from_file_location("crest_v32", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.REFERENCE_COMPOUNDS


def embed(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        sys.exit(f"error: RDKit could not parse SMILES:\n  {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        # retry with random coords for a stubborn macrocycle
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            sys.exit("error: RDKit embedding failed")
    return mol


def nn_atoms_0based(mol):
    pair = ce.diazirine_nn_atoms(mol)      # 1-based (n1, n2) or None
    if pair is None:
        sys.exit("error: no diazirine (C-N=N ring) found in this molecule")
    return pair, (pair[0] - 1, pair[1] - 1)   # (1-based for xtb, 0-based for measuring)


def alkyne_atoms_0based(mol):
    m = mol.GetSubstructMatches(Chem.MolFromSmarts(ALKYNE_SMARTS))
    return (m[0][0], m[0][1]) if m else None


def read_coords(xyz_text):
    lines = xyz_text.splitlines()
    n = int(lines[0].split()[0])
    return np.array([[float(x) for x in ln.split()[1:4]] for ln in lines[2:2 + n]])


def dist(coords, i, j):
    return float(np.linalg.norm(coords[i] - coords[j]))


def verdict(nn):
    if NN_PASS[0] <= nn <= NN_PASS[1]:
        return "PASS"
    if nn <= NN_WATCH_HI:
        return "WATCH"
    return "FAIL"


def constrained_opt(xyz_text, solvent, charge, fc, nn_1based, xtb_bin, workdir):
    """One GFN2 --opt (with an N=N restraint if fc>0). Returns optimized xyz text, or None."""
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / "conf.xyz").write_text(xyz_text)
    cmd = [xtb_bin, "conf.xyz", "--opt", "--gfn", "2", "--alpb", solvent, "--chrg", str(charge)]
    if fc > 0:
        cf = workdir / "constrain.inp"
        ce.write_constraint_file(cf, nn_1based, value=ce.DIAZIRINE_NN, fc=fc)
        cmd += ["--input", str(cf)]
    try:
        subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=1800)
    except (subprocess.TimeoutExpired, FileNotFoundError) as exc:
        print(f"    xtb failed (fc={fc}): {type(exc).__name__}", file=sys.stderr)
        return None
    opt = workdir / "xtbopt.xyz"
    return opt.read_text() if opt.exists() else None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compound", type=int, default=24,
                    help="index into REFERENCE_COMPOUNDS (default 24 = the worst-failing diazirine)")
    ap.add_argument("--smiles", default=None, help="use this SMILES instead of --compound")
    ap.add_argument("--solvent", default="water", help="ALPB solvent for the opt (default water)")
    ap.add_argument("--charge", type=int, default=None, help="override formal charge")
    ap.add_argument("--fcs", type=float, nargs="+", default=[0, 0.25, 0.5, 1.0, 2.0, 5.0],
                    help="force-constant ladder (Eh/Bohr^2); 0 = unconstrained baseline")
    ap.add_argument("--xtb", default="xtb", help="xtb binary (default: xtb on PATH)")
    args = ap.parse_args(argv)

    if args.smiles:
        smiles, name = args.smiles, "custom"
    else:
        ref = load_reference_compounds()
        if not (0 <= args.compound < len(ref)):
            sys.exit(f"error: --compound {args.compound} out of range (0..{len(ref) - 1})")
        smiles, name = ref[args.compound]["smiles"], ref[args.compound].get("name", "?")

    mol = embed(smiles)
    charge = args.charge if args.charge is not None else Chem.GetFormalCharge(mol)
    nn_1based, nn_0based = nn_atoms_0based(mol)
    alk = alkyne_atoms_0based(mol)

    # starting geometry sanity (RDKit should give ~1.23-1.24)
    start = read_coords(Chem.MolToXYZBlock(mol))
    nn_start = dist(start, *nn_0based)

    print(f"compound {args.compound} ({name}) | solvent {args.solvent} | charge {charge}")
    print(f"N=N atoms (1-based): {nn_1based} | RDKit start N=N = {nn_start:.3f} A "
          f"| target {ce.DIAZIRINE_NN} A (PASS {NN_PASS[0]}-{NN_PASS[1]})")
    if alk:
        print(f"alkyne C#C atoms (0-based): {alk}")
    print(f"\n  {'fc (Eh/Bohr^2)':>16}   {'N=N (A)':>8}   {'verdict':>7}"
          + (f"   {'C#C (A)':>8}" if alk else ""))

    with tempfile.TemporaryDirectory(prefix="fcsweep_") as td:
        xyz0 = Chem.MolToXYZBlock(mol)
        for fc in args.fcs:
            out = constrained_opt(xyz0, args.solvent, charge, fc, nn_1based, args.xtb,
                                  Path(td) / f"fc_{fc}")
            if out is None:
                print(f"  {fc:>16.2f}   {'--':>8}   {'ERR':>7}")
                continue
            c = read_coords(out)
            nn = dist(c, *nn_0based)
            row = f"  {fc:>16.2f}   {nn:>8.3f}   {verdict(nn):>7}"
            if alk:
                row += f"   {dist(c, *alk):>8.3f}"
            if fc == 0:
                row += "   <- unconstrained baseline (expect ~1.43)"
            print(row)

    print("\nPick the SMALLEST fc whose N=N is PASS (1.22-1.25). If C#C stays ~1.30 at all fc, the")
    print("alkyne is unconstrained and distorting too -> add an alkyne constraint before the regen.")


if __name__ == "__main__":
    main()

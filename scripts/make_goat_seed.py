#!/usr/bin/env python3
"""
make_goat_seed.py  --  SMILES -> ETKDG 3D seed.xyz for a GOAT run (FlexiSol-style naive seed).

FlexiSol seeded GOAT from a SMILES-derived, tight-binding-optimized structure -- NOT from a CREST
ensemble. This reproduces that: SMILES (from the SAME REFERENCE_COMPOUNDS registry CREST uses, so
GOAT and CREST start from identical chemistry) -> ETKDGv3 embed (macrocycle-aware) -> MMFF/UFF
clean-up -> seed.xyz. GOAT then does the global search from there. Charge is taken from the SMILES
formal charge (matches what CREST used).

Run in an env with RDKit (e.g. chameleon_crest212). --index pulls SMILES from crest_v3.2.py;
--smiles passes one inline (for molecules not in the registry, e.g. the Fairlie 6-mers).
"""
import argparse, sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

# molecules not in REFERENCE_COMPOUNDS -> inline SMILES (verified against 7L96/7L98)
INLINE = {
    "cmpd4":  "CC(C)C[C@@H]1NC(=O)[C@@H](CC(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc2ccc(O)cc2)N(C)C(=O)[C@H]2CCCN2C(=O)[C@H](CC(C)C)NC1=O",
    "cmpd10": "CC[C@H](C)[C@@H]1NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc2ccccc2)N(C)C(=O)[C@H](C)NC(=O)c2csc(n2)[C@H]([C@@H](C)CC)NC(=O)[C@@H]2CCCN2C1=O",
}
# name -> registry index (SMILES pulled from crest_v3.2.py REFERENCE_COMPOUNDS)
REGISTRY = {
    "CsA": 1, "3-12-8-12_R_xylene": 5, "3-12-8-12_S_xylene": 6, "6-4-4-13_xylene": 7,
    "1-6-4-7_xylene": 20, "2-9-9-8_xylene": 22,
}

def smiles_from_registry(idx):
    import importlib.util
    p = Path(__file__).resolve().parent / "crest_v3.2.py"
    try:
        spec = importlib.util.spec_from_file_location("crest_registry", p)
        mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
        return mod.REFERENCE_COMPOUNDS[idx]["smiles"]
    except Exception as e:
        sys.exit(f"could not read REFERENCE_COMPOUNDS[{idx}] from {p} ({e}); run in the crest env or pass --smiles")

def embed(smiles, out):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        sys.exit(f"invalid SMILES: {smiles}")
    charge = Chem.GetFormalCharge(mol)
    mol = Chem.AddHs(mol)
    p = AllChem.ETKDGv3(); p.randomSeed = 42; p.useMacrocycleTorsions = True
    if AllChem.EmbedMolecule(mol, p) != 0:
        p.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, p) != 0:
            sys.exit("ETKDG embedding failed (macrocycle?) -- seed this one from an existing conformer instead")
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    except Exception:
        try: AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
        except Exception: pass
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), f"seed charge={charge}"]
    for i, a in enumerate(mol.GetAtoms()):
        q = conf.GetAtomPosition(i)
        lines.append(f"{a.GetSymbol():<3s} {q.x:>15.8f} {q.y:>15.8f} {q.z:>15.8f}")
    Path(out).write_text("\n".join(lines) + "\n")
    print(f"seed -> {out}  ({mol.GetNumAtoms()} atoms)")
    print(f"SEED_CHARGE={charge}")   # captured by the wrapper

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", help="cmpd4/cmpd10 or a registry name (CsA, 1-6-4-7_xylene, ...)")
    ap.add_argument("--index", type=int, help="REFERENCE_COMPOUNDS index (overrides --name lookup)")
    ap.add_argument("--smiles", help="explicit SMILES (overrides everything)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.smiles:
        smi = a.smiles
    elif a.index is not None:
        smi = smiles_from_registry(a.index)
    elif a.name in INLINE:
        smi = INLINE[a.name]
    elif a.name in REGISTRY:
        smi = smiles_from_registry(REGISTRY[a.name])
    else:
        sys.exit(f"unknown --name '{a.name}'; known: {list(INLINE)+list(REGISTRY)}")
    embed(smi, a.out)

if __name__ == "__main__":
    main()

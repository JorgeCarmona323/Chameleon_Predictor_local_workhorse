# env: chameleon-calc
import sys
sys.path.insert(0, "scripts")
from importlib import import_module
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
m = import_module("crest_v3.2")
rc = m.REFERENCE_COMPOUNDS
print(f"{len(rc)} compounds loaded")
bad = [(i, c["short"]) for i, c in enumerate(rc) if Chem.MolFromSmiles(c["smiles"]) is None]
print("INVALID:", bad if bad else "none - all SMILES valid")
for i, c in enumerate(rc):
    if i >= 10:
        mol = Chem.MolFromSmiles(c["smiles"])
        chg = sum(a.GetFormalCharge() for a in mol.GetAtoms())
        print(f"  {i:2d} {c['short']:12s} {mol.GetNumAtoms():3d} atoms, charge {chg}")

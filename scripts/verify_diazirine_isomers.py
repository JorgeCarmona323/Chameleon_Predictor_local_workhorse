# env: chameleon-calc
"""Verify the 6 new compounds (3-12-10-12 sarcosine pair + 4 diazirine variants):
validity, CIP at thiophene, stereocenter count, diazirine + CF3 presence, diastereomer pairing."""
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
RDLogger.DisableLog("rdApp.*")

new = {
 "3-12-10-12 R":            "CN(CC(N[C@@H](CO)C(N[C@@H](c1ccc[s]1)C(N[C@@H](CSCc1c(CSC[C@@H](C(N[C@@H]2CO)=O)NC(CCC#C)=O)cccc1)C(N)=O)=O)=O)=O)C2=O",
 "3-12-10-12 S":            "CN(CC(N[C@@H](CO)C(N[C@H](c1ccc[s]1)C(N[C@@H](CSCc1c(CSC[C@@H](C(N[C@@H]2CO)=O)NC(CCC#C)=O)cccc1)C(N)=O)=O)=O)=O)C2=O",
 "3-12-8-12 R Diazirine":   "C#CCCC(N[C@@H](CSCc1cc(C2(C(F)(F)F)N=N2)cc(CSC[C@@H](C(N)=O)NC([C@H](c2ccc[s]2)NC([C@H](CO)NC([C@H](CC2)N2C([C@H](CO)N2)=O)=O)=O)=O)c1)C2=O)=O",
 "3-12-8-12 S Diazirine":   "C#CCCC(N[C@@H](CSCc1cc(C2(C(F)(F)F)N=N2)cc(CSC[C@@H](C(N)=O)NC([C@@H](c2ccc[s]2)NC([C@H](CO)NC([C@H](CC2)N2C([C@H](CO)N2)=O)=O)=O)=O)c1)C2=O)=O",
 "3-12-10-12 R Diazirine":  "CN(CC(N[C@@H](CO)C(N[C@@H](c1ccc[s]1)C(N[C@@H](CSCc1cc(CSC[C@@H](C(N[C@@H]2CO)=O)NC(CCC#C)=O)cc(C3(C(F)(F)F)N=N3)c1)C(N)=O)=O)=O)=O)C2=O",
 "3-12-10-12 S Diazirine":  "CN(CC(N[C@@H](CO)C(N[C@H](c1ccc[s]1)C(N[C@@H](CSCc1cc(CSC[C@@H](C(N[C@@H]2CO)=O)NC(CCC#C)=O)cc(C3(C(F)(F)F)N=N3)c1)C(N)=O)=O)=O)=O)C2=O",
}

DIAZIRINE = Chem.MolFromSmarts("[C]1([CX4](F)(F)F)N=N1")   # 3-CF3-3H-diazirine
CF3 = Chem.MolFromSmarts("[CX4](F)(F)F")

def thio_cip(m, centers):
    for idx, cip in centers:
        for nb in m.GetAtomWithIdx(idx).GetNeighbors():
            if nb.GetIsAromatic():
                for r in m.GetRingInfo().AtomRings():
                    if nb.GetIdx() in r and any(m.GetAtomWithIdx(x).GetSymbol()=="S" for x in r):
                        return cip
    return "?"

print(f"{'compound':26s} {'valid':5s} {'formula':18s} {'nSt':3s} {'thioCIP':7s} {'diazir':6s} {'CF3':3s} {'inchikey':27s}")
rows = {}
for name, smi in new.items():
    m = Chem.MolFromSmiles(smi)
    if m is None:
        print(f"{name:26s} INVALID"); continue
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    centers = Chem.FindMolChiralCenters(m, useLegacyImplementation=False)
    ik = Chem.InchiToInchiKey(Chem.MolToInchi(m))
    rows[name] = (ik, thio_cip(m, centers))
    print(f"{name:26s} {'yes':5s} {rdMolDescriptors.CalcMolFormula(m):18s} "
          f"{len(centers):<3d} {thio_cip(m,centers):7s} "
          f"{'yes' if m.HasSubstructMatch(DIAZIRINE) else 'NO':6s} "
          f"{'yes' if m.HasSubstructMatch(CF3) else 'no':3s} {ik}")

print("\n=== diastereomer pair checks (same skeleton, opposite thiophene CIP) ===")
for base in ["3-12-10-12", "3-12-8-12 R Diazirine"]:
    pass
pairs = [("3-12-10-12 R","3-12-10-12 S"),
         ("3-12-8-12 R Diazirine","3-12-8-12 S Diazirine"),
         ("3-12-10-12 R Diazirine","3-12-10-12 S Diazirine")]
for a,b in pairs:
    ika, ca = rows[a]; ikb, cb = rows[b]
    print(f"  {a:24s} ({ca}) / {b:24s} ({cb}): "
          f"same skeleton={ika[:14]==ikb[:14]}, R=CIP-R & S=CIP-S={ca=='R' and cb=='S'}")

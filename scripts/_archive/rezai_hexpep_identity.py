# env: chameleon-calc
"""
rezai_hexpep_identity.py
------------------------
(1) Confirms HexPep's identity among Rezai 2006's 9 cyclic-hexapeptide
    diastereomers by building each from its Table-1 sequence (with a
    build->CIP self-check) and matching canonical SMILES.
    RESULT: HexPep = Rezai COMPOUND 1 (permeable diastereomer, logP_E -6.2).

(2) Provides `rezai_pattern(mol)`, which returns the backbone Cα residues in
    Rezai ring order [Leu1, Leu2, Leu3, Leu4, Pro5] (Tyr6 = anchor). Reuse this
    to do PER-RESIDUE ³J assignment in validate_hexpep_nmr.py (pair each computed
    Leu_i coupling with its named experimental value instead of a sorted set).

Usage:  python scripts/rezai_hexpep_identity.py
"""
from rdkit import Chem

HEXPEP = ("CC(C)C[C@@H]1NC(=O)[C@@H](CC(C)C)NC(=O)[C@@H](CC(C)C)NC(=O)"
          "[C@H](Cc2ccc(O)cc2)NC(=O)[C@@H]2CCCN2C(=O)[C@@H](CC(C)C)NC1=O")

# Rezai Table 1 (order Leu1,Leu2,Leu3,Leu4,Pro5,Tyr6). D/L per residue.
TABLE1 = {
    1: ["D", "D", "L", "D", "L", "L"], 2: ["D", "D", "D", "D", "L", "L"],
    3: ["L", "L", "L", "D", "L", "L"], 4: ["L", "D", "D", "D", "L", "L"],
    5: ["L", "L", "L", "L", "D", "L"], 6: ["D", "D", "D", "D", "D", "L"],
    7: ["L", "L", "D", "D", "L", "L"], 8: ["L", "D", "L", "D", "D", "L"],
    9: ["L", "D", "L", "L", "D", "L"],
}


def build(seq):
    def c(dl):
        return "@@H" if dl == "L" else "@H"
    r = f"N1[C{c(seq[0])}](CC(C)C)C(=O)"
    for i in (1, 2, 3):
        r += f"N[C{c(seq[i])}](CC(C)C)C(=O)"
    r += f"N2CCC[C{c('D' if seq[4]=='L' else 'L')}]2C(=O)"       # Pro (fragment convention inverted)
    r += f"N[C{c(seq[5])}](Cc3ccc(O)cc3)C1=O"                    # Tyr, closes macro ring
    return r


def rezai_pattern(mol):
    """Ring-ordered D/L for [Leu1,Leu2,Leu3,Leu4,Pro5] (Tyr6 = anchor).
    Returns list of (atom_idx, restype, 'L'/'D') so callers can map couplings."""
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    centers = dict(Chem.FindMolChiralCenters(mol, useLegacyImplementation=False))
    rings = mol.GetRingInfo().AtomRings()
    macro = set(max(rings, key=len))

    def restype(a):
        if any(len(r) == 5 and a.GetIdx() in r for r in rings):
            return "Pro"
        seen, st = set(), [n.GetIdx() for n in a.GetNeighbors() if n.GetIdx() not in macro]
        while st:
            j = st.pop()
            if j in seen:
                continue
            seen.add(j)
            if mol.GetAtomWithIdx(j).GetIsAromatic():
                return "Tyr"
            for nb in mol.GetAtomWithIdx(j).GetNeighbors():
                if nb.GetIdx() not in macro and nb.GetIdx() not in seen:
                    st.append(nb.GetIdx())
        return "Leu"

    ca = {}
    for idx, rs in centers.items():
        a = mol.GetAtomWithIdx(idx)
        if a.GetSymbol() == "C" and idx in macro and any(n.GetSymbol() == "N" for n in a.GetNeighbors()):
            ca[idx] = (restype(a), "L" if rs == "S" else "D")

    def bb_adj(i):
        res, seen, frontier = [], {i}, [n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() in macro]
        while frontier:
            j = frontier.pop(0)
            if j in seen:
                continue
            seen.add(j)
            if j in ca:
                res.append(j)
                continue
            for nb in mol.GetAtomWithIdx(j).GetNeighbors():
                if nb.GetIdx() in macro and nb.GetIdx() not in seen:
                    frontier.append(nb.GetIdx())
        return res

    tyr = [i for i in ca if ca[i][0] == "Tyr"][0]
    leu1 = [i for i in bb_adj(tyr) if ca[i][0] == "Leu"][0]
    order = [tyr, leu1]
    while len(order) < 6:
        nxts = [x for x in bb_adj(order[-1]) if x not in order]
        if not nxts:
            break
        order.append(nxts[0])
    return [(i, ca[i][0], ca[i][1]) for i in order[1:]]  # Leu1..Pro5


def main():
    built = {}
    print("self-check (build -> CIP == Table 1):")
    for n, seq in TABLE1.items():
        m = Chem.MolFromSmiles(build(seq))
        pat = [x[2] for x in rezai_pattern(m)] + ["L"]
        built[n] = Chem.MolToSmiles(m)
        print(f"  cmpd {n}: {'OK' if pat == seq else 'MISMATCH'}")
    mh = Chem.MolFromSmiles(HEXPEP)
    patt = rezai_pattern(mh)
    print("\nHexPep ring order (Leu1..Pro5):", [(x[1], x[2]) for x in patt])
    match = [n for n, s in built.items() if s == Chem.MolToSmiles(mh)]
    print("HexPep canonical-SMILES match:", f"compound {match[0]}" if match else "NONE")


if __name__ == "__main__":
    main()

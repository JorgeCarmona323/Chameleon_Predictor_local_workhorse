# env: chameleon-calc
"""
validate_hexpep_nmr.py
----------------------
Validate the HexPep CREST ensemble against the experimental ³J(HN-Hα) coupling
constants in the Rezai 2006 SI. HexPep is confirmed by canonical SMILES to be
Rezai COMPOUND 1 (their most-permeable diastereomer, logP_E -6.2) — NOT compound 9.
NB: the DB's "permeable: False" is correct (-6.20 < the -6.0 threshold); compound 1 is
just the most-permeable of the Rezai series. Match to cmpd_1 = pass; cmpd_9 = control.

For each backbone amide NH (Leu×4 + Tyr; Pro has no NH) we measure the
H-N-Cα-Hα dihedral per conformer, convert to ³J via a Karplus equation, and
ensemble-average. Compares the predicted coupling SET to the two experimental
diastereomer sets (labelling of the 4 Leu is resolved by sorted-set MAE; Tyr is
anchored by its aromatic side chain). Experimental error is ±2 Hz.

Karplus (Vuister & Bax, JACS 1993): ³J = 6.51 cos²θ − 1.76 cosθ + 1.60,
θ = dihedral(HN, N, Cα, Hα).

Usage:
  python scripts/validate_hexpep_nmr.py \
      --xyz results/conformers/HexPep/aq/full_ensemble.xyz \
      --out results/2026-07-08_hexpep_nmr_validation.txt
"""
import argparse
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolTransforms

RDLogger.DisableLog("rdApp.*")

HEXPEP_SMILES = ("CC(C)C[C@@H]1NC(=O)[C@@H](CC(C)C)NC(=O)[C@@H](CC(C)C)NC(=O)"
                 "[C@H](Cc2ccc(O)cc2)NC(=O)[C@@H]2CCCN2C(=O)[C@@H](CC(C)C)NC1=O")
RT_KCAL = 0.5924096
HARTREE_KCAL = 627.5094740631

# Rezai 2006 SI experimental ³J(HN-Hα), Hz (±2), per diastereomer.
# HexPep is confirmed (canonical SMILES) to be COMPOUND 1 (permeable) — so it
# should match cmpd_1; cmpd_9 (impermeable) is the wrong-diastereomer control.
#   Table 3 (compound 9): Leu 9.0/7.8/8.4/6.6, Tyr 12.0
#   Table 2 (compound 1): Leu 4.0/10.2/9.0/7.2, Tyr 8.4
EXP = {
    "compound_1 (HexPep, permeable)":       {"Leu": [4.0, 10.2, 9.0, 7.2], "Tyr": 8.4},
    "compound_9 (impermeable, control)":    {"Leu": [9.0, 7.8, 8.4, 6.6], "Tyr": 12.0},
}


def karplus(theta_deg):
    t = np.radians(theta_deg)
    return 6.51 * np.cos(t) ** 2 - 1.76 * np.cos(t) + 1.60


def parse_xyz(path):
    lines = Path(path).read_text().splitlines()
    confs, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        comment = lines[i + 1] if i + 1 < len(lines) else ""
        try:
            e = float(comment.split()[0])
        except Exception:
            e = np.nan
        block = lines[i + 2: i + 2 + n]
        els = [b.split()[0] for b in block]
        xyz = np.array([[float(v) for v in b.split()[1:4]] for b in block])
        confs.append((els, xyz, e))
        i += 2 + n
    return confs


def build_template():
    mol = Chem.MolFromSmiles(HEXPEP_SMILES)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    return mol


def find_nh_groups(mol):
    """Return list of dicts: {res, hn, n, ca, ha} for each backbone amide NH."""
    ri = mol.GetRingInfo()
    ring = set(max(ri.AtomRings(), key=len)) if ri.AtomRings() else set()
    groups = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N" or atom.GetIdx() not in ring:
            continue
        h = next((nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetSymbol() == "H"), None)
        if h is None:
            continue  # Pro / N-methyl
        # carbonyl C neighbour (amide) and Cα neighbour (ring C with an H)
        ca = None
        for nb in atom.GetNeighbors():
            if nb.GetSymbol() != "C" or nb.GetIdx() not in ring:
                continue
            is_carbonyl = any(b.GetBondTypeAsDouble() == 2.0 and
                              b.GetOtherAtom(nb).GetSymbol() == "O" for b in nb.GetBonds())
            if is_carbonyl:
                continue
            ca = nb
        if ca is None:
            continue
        ha = next((nb.GetIdx() for nb in ca.GetNeighbors() if nb.GetSymbol() == "H"), None)
        if ha is None:
            continue
        # residue type by Cα side chain: aromatic → Tyr, else Leu
        res = "Leu"
        for nb in ca.GetNeighbors():
            if nb.GetSymbol() == "C" and nb.GetIdx() not in ring:
                if any(a.GetIsAromatic() for a in _reachable(mol, nb.GetIdx(), ring | {ca.GetIdx()})):
                    res = "Tyr"
        groups.append({"res": res, "hn": h, "n": atom.GetIdx(),
                       "ca": ca.GetIdx(), "ha": ha})
    return groups


def _reachable(mol, start, blocked, depth=6):
    seen, frontier = set(), [(start, 0)]
    while frontier:
        idx, d = frontier.pop()
        if idx in seen or d > depth:
            continue
        seen.add(idx)
        for nb in mol.GetAtomWithIdx(idx).GetNeighbors():
            if nb.GetIdx() not in blocked:
                frontier.append((nb.GetIdx(), d + 1))
    return [mol.GetAtomWithIdx(i) for i in seen]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xyz", default="results/conformers/HexPep/aq/full_ensemble.xyz")
    ap.add_argument("--out", default="results/2026-07-08_hexpep_nmr_validation.txt")
    args = ap.parse_args()

    confs = parse_xyz(args.xyz)
    tmpl = build_template()
    tel = [a.GetSymbol() for a in tmpl.GetAtoms()]
    cel = confs[0][0]
    order_ok = (tel == cel)

    lines = []
    lines.append(f"HexPep NMR validation — {len(confs)} conformers, template {len(tel)} atoms")
    lines.append(f"atom-order match template vs xyz: {order_ok}")
    if not order_ok:
        # report first mismatch and abort coupling calc
        diff = next((k for k, (a, b) in enumerate(zip(tel, cel)) if a != b), None)
        lines.append(f"  first diff at atom {diff}: template={tel[diff] if diff is not None else '?'} "
                     f"xyz={cel[diff] if diff is not None else '?'} — cannot map, aborting")
        Path(args.out).write_text("\n".join(lines), encoding="utf-8")
        print("ORDER MISMATCH — see out")
        return

    groups = find_nh_groups(tmpl)
    lines.append(f"backbone NH groups found: {[g['res'] for g in groups]}  (expect 4 Leu + 1 Tyr)")

    E = np.array([c[2] for c in confs], float) * HARTREE_KCAL
    Erel = E - np.nanmin(E)
    w = np.exp(-Erel / RT_KCAL)
    w = w / np.nansum(w)

    # per-group predicted J per conformer
    conf = Chem.Conformer(tmpl.GetNumAtoms())
    J = {id(g): [] for g in groups}
    for els, xyz, e in confs:
        for k, (x, y, z) in enumerate(xyz):
            conf.SetAtomPosition(k, (float(x), float(y), float(z)))
        for g in groups:
            th = rdMolTransforms.GetDihedralDeg(conf, g["hn"], g["n"], g["ca"], g["ha"])
            J[id(g)].append(karplus(th))
    for g in groups:
        arr = np.array(J[id(g)])
        g["J_bw"] = float(np.nansum(w * arr))
        g["J_mean"] = float(np.nanmean(arr))

    tyr = [g for g in groups if g["res"] == "Tyr"]
    leu = [g for g in groups if g["res"] == "Leu"]

    def preds(key):
        pl = sorted(g[key] for g in leu)
        pt = tyr[0][key] if tyr else float("nan")
        return pl, pt

    def errstats(pl, pt, ex):
        e = np.array(pl + [pt]) - np.array(sorted(ex["Leu"]) + [ex["Tyr"]])
        return float(np.mean(np.abs(e))), float(np.sqrt(np.mean(e ** 2)))

    # UNIFORM MEAN is the primary (energy-free) metric; Boltzmann is a secondary
    # sensitivity check. RMSD is primary (matches the Vuister-Bax 0.73 Hz floor).
    for scheme, key in [("UNIFORM MEAN  (primary, energy-free)", "J_mean"),
                        ("Boltzmann     (secondary sensitivity check)", "J_bw")]:
        pl, pt = preds(key)
        lines.append(f"\n=== {scheme} ===")
        lines.append(f"Predicted ³J: Tyr {pt:.2f} | Leu(sorted) {[round(x, 2) for x in pl]}")
        for name, ex in EXP.items():
            mae, rmsd = errstats(pl, pt, ex)
            lines.append(f"  vs {name}: RMSD={rmsd:.2f} Hz  (MAE {mae:.2f})  "
                         f"| exp Tyr {ex['Tyr']} Leu {sorted(ex['Leu'])}")
    lines.append("\n[exp error ±2 Hz; Vuister-Bax Karplus RMSD floor 0.73 Hz]")

    Path(args.out).write_text("\n".join(lines), encoding="utf-8")
    print("done")


if __name__ == "__main__":
    main()

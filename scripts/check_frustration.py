# env: chameleon-calc
"""
check_frustration.py
--------------------
Test the conformational-frustration hypothesis for 3-12-8-12 S (azetidine):
the macrocycle closes around a FIXED azetidine turn in MULTIPLE backbone H-bond
patterns, because the rigid turn conflicts with clean closure.

Three falsifiable predictions of frustration, with controls (8-12 R = open/no
frustration; 10-12 S = sarcosine/clean fold):

  (1) multiple folds sharing the cis turn but DIFFERING in backbone IMHB
  (2) those folds NEAR-DEGENERATE in energy
  (3) per-atom RMSF: azetidine RIGID (low) while backbone is MOBILE (high)
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.ML.Cluster import Butina
import sys; sys.path.insert(0, "scripts")
from phys_descriptors_v3 import imhb_descriptors_mol, backbone_hbond_atoms, macrocycle_atoms

RDLogger.DisableLog("rdApp.*")
HART = 627.509
BASE = Path("results/conformers")
TARGETS = {
    "8-12 S  (Aze, predicted FRUSTRATED)": "DOPC 3-12-8-12 S",
    "8-12 R  (Aze, open = control)":       "DOPC 3-12-8-12 R",
    "10-12 S (Sar, clean fold = control)": "3-12-10-12 S",
}


def load(wd):
    confs = json.load(open(wd / "ensemble.json"))["conformers"]
    e = np.array([c.get("totalenergy", np.nan) for c in confs], float) * HART
    w = np.array([c.get("boltzmannweight", np.nan) for c in confs], float)
    mols = [m for m in Chem.SDMolSupplier(str(wd / "ensemble.sdf"), removeHs=False, sanitize=True) if m]
    n = min(len(mols), len(w))
    return mols[:n], w[:n], e[:n]


def azetidine_atoms(mol):
    for r in mol.GetRingInfo().AtomRings():
        if len(r) == 4 and sum(mol.GetAtomWithIdx(i).GetSymbol() == "N" for i in r) == 1:
            return set(r)
    return set()


def folds(mols, w, e, cap=80, rms=1.0):
    sub = [int(i) for i in np.argsort(w)[::-1][:cap]]
    base = Chem.Mol(mols[sub[0]]); base.RemoveAllConformers()
    for i in sub:
        base.AddConformer(Chem.Conformer(mols[i].GetConformer()), assignId=True)
    heavy = [a.GetIdx() for a in base.GetAtoms() if a.GetAtomicNum() > 1]
    mat = AllChem.GetConformerRMSMatrix(base, atomIds=heavy, prealigned=False)
    cl = Butina.ClusterData(mat, len(sub), rms, isDistData=True, reordering=True)
    bb = backbone_hbond_atoms(mols[0])
    bbimhb = {j: imhb_descriptors_mol(mols[i], -1, bb)["imhb_bb"] for j, i in enumerate(sub)}
    wsub, esub, emin = w[sub], e[sub], np.nanmin(e)
    rows = []
    for c in cl:
        pop = float(np.nansum(wsub[list(c)]) / np.nansum(w))
        if pop < 0.03:
            continue
        bbv = [bbimhb[j] for j in c]
        rows.append((pop, float(np.nanmin(esub[list(c)]) - emin),
                     float(np.mean(bbv)), int(min(bbv)), int(max(bbv))))
    return sorted(rows, reverse=True)


def per_atom_rmsf(mols, w):
    keep = np.where(np.isfinite(w) & (w > 0))[0]
    ref = int(keep[np.argmax(w[keep])])
    order = [ref] + [int(i) for i in keep if i != ref]
    base = Chem.Mol(mols[ref]); base.RemoveAllConformers()
    for i in order:
        base.AddConformer(Chem.Conformer(mols[i].GetConformer()), assignId=True)
    heavy = [a.GetIdx() for a in base.GetAtoms() if a.GetAtomicNum() > 1]
    rdMolAlign.AlignMolConformers(base, atomIds=heavy)
    P = np.array([c.GetPositions() for c in base.GetConformers()])
    ww = w[order] / w[order].sum()
    mean = (ww[:, None, None] * P).sum(0)
    return np.sqrt((ww[:, None] * ((P - mean) ** 2).sum(2)).sum(0))   # per-atom RMSF


def main():
    for label, dname in TARGETS.items():
        mols, w, e = load(BASE / dname / "water")
        print("\n" + "=" * 66)
        print(label)
        print("=" * 66)
        # (1)+(2) folds: population, relative energy, backbone-IMHB spread
        print(" fold  pop%   relE(kcal)  backbone-IMHB(mean[min-max])")
        for k, (pop, relE, mean, lo, hi) in enumerate(folds(mols, w, e), 1):
            print(f"  {k:2d}   {100*pop:5.1f}   {relE:7.2f}      {mean:4.2f} [{lo}-{hi}]")
        # (3) per-atom RMSF localization
        az = azetidine_atoms(mols[0])
        macro = macrocycle_atoms(mols[0])
        rmsf = per_atom_rmsf(mols, w)
        heavy = [a.GetIdx() for a in mols[0].GetAtoms() if a.GetAtomicNum() > 1]
        bb_non_az = [i for i in macro if i not in az and i in heavy]
        if az:
            print(f"  per-atom RMSF (A):  azetidine={np.mean([rmsf[i] for i in az]):.2f}  "
                  f"backbone(non-Aze)={np.mean([rmsf[i] for i in bb_non_az]):.2f}  "
                  f"ratio bb/Aze={np.mean([rmsf[i] for i in bb_non_az])/np.mean([rmsf[i] for i in az]):.1f}x")
        else:
            print(f"  per-atom RMSF (A):  no azetidine; backbone(macrocycle)={np.mean([rmsf[i] for i in macro if i in heavy]):.2f}")


if __name__ == "__main__":
    main()

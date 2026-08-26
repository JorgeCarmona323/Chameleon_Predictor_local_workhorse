#!/usr/bin/env python3
"""
pre_post_compare.py  --  PRE (GFN2) vs POST (r2SCAN-3c/CENSO) descriptor comparison.

For each compound: population-weighted PSA / Rgyr / IMHB / shape on the GFN2 ensemble
(weighted by CPCM-X fe pops) vs the refined ensemble (weighted by r2SCAN-3c pops from
refined/energies.csv). Answers the PI's question: does refinement move the descriptors?

Reuses the exact descriptor definitions from phys_descriptors_v3 (same as every report).
Run locally in rdkit_env after syncing the `refined/` dirs down.

Manifest CSV (--manifest) columns:
    name,pre_sdf,pre_fe_csv,pre_leg,post_dir
      pre_sdf    = <leg>/ensemble.sdf              (GFN2 conformers, bond template + coords)
      pre_fe_csv = results/free_energy/fe_xylene_<label>.csv
      pre_leg    = chloroform | hexane | water     (which leg's CPCM-X pops to weight PRE by)
      post_dir   = <leg>/refined                   (has ensemble.xyz + energies.csv)
"""
import argparse, csv, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from phys_descriptors_v3 import surface_descriptors_mol, count_hbonds_xyz
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Geometry import Point3D
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

def read_xyz(path):
    L = Path(path).read_text().splitlines(); frames = []; i = 0
    while i < len(L):
        s = L[i].strip()
        if not s: i += 1; continue
        try: n = int(s)
        except ValueError: break
        coords = [tuple(map(float, L[i+2+j].split()[1:4])) for j in range(n)]
        frames.append(coords); i += n + 2
    return frames

def desc(mol, cid):
    sd = surface_descriptors_mol(mol, cid)
    conf = mol.GetConformer(cid); co = conf.GetPositions()
    sy = [a.GetSymbol() for a in mol.GetAtoms()]
    hv = np.array([s != "H" for s in sy]); ch = co[hv]
    rg = float(np.sqrt(((ch - ch.mean(0))**2).sum(1).mean()))
    return dict(psa=sd["psa"], rg=rg, imhb=count_hbonds_xyz(sy, co),
                npr1=rdMolDescriptors.CalcNPR1(mol, confId=cid),
                npr2=rdMolDescriptors.CalcNPR2(mol, confId=cid))

def wmean(vals, w):
    vals = np.array(vals, float); w = np.array(w, float); w = w / w.sum()
    return float((w * vals).sum())

def weighted_desc(mols_confs, weights):
    keys = ["psa", "rg", "imhb", "npr1", "npr2"]
    D = {k: [] for k in keys}
    for (mol, cid) in mols_confs:
        d = desc(mol, cid)
        for k in keys: D[k].append(d[k])
    return {k: wmean(D[k], weights) for k in keys}

def load_pre(pre_sdf, fe_csv, leg):
    mols = [m for m in Chem.SDMolSupplier(pre_sdf, removeHs=False) if m]
    pops = {}
    with open(fe_csv) as f:
        for r in csv.DictReader(f):
            if r["solvent"] == leg: pops[int(r["conf"])] = float(r["pop"])
    confs, ws = [], []
    for c in sorted(pops, key=lambda c: -pops[c]):
        if c < len(mols): confs.append((mols[c], mols[c].GetConformer().GetId())); ws.append(pops[c])
    return mols[0], confs, ws

def load_post(template, post_dir):
    coords = read_xyz(Path(post_dir) / "ensemble.xyz")
    pops = {}
    ep = Path(post_dir) / "energies.csv"
    if ep.exists():
        for r in csv.DictReader(open(ep)): pops[int(r["conf"])] = float(r["pop"])
    confs, ws = [], []
    for k, xyz in enumerate(coords):
        if template.GetNumAtoms() != len(xyz):
            raise SystemExit(f"atom-count mismatch: template {template.GetNumAtoms()} vs refined {len(xyz)} (atom order?)")
        m = Chem.Mol(template); conf = m.GetConformer()
        for i, (x, y, z) in enumerate(xyz): conf.SetAtomPosition(i, Point3D(x, y, z))
        confs.append((m, conf.GetId())); ws.append(pops.get(k, 1.0))
    return confs, ws

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", default="results/validation/refine_pre_post")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.manifest)))
    table = []
    for r in rows:
        name = r["name"]
        template, pre_confs, pre_w = load_pre(r["pre_sdf"], r["pre_fe_csv"], r["pre_leg"])
        post_confs, post_w = load_post(template, r["post_dir"])
        pre = weighted_desc(pre_confs, pre_w)
        post = weighted_desc(post_confs, post_w)
        print(f"\n=== {name}  (pre {len(pre_confs)} confs / post {len(post_confs)}) ===")
        for k in ["psa", "rg", "imhb", "npr1", "npr2"]:
            print(f"   {k:5s}  pre {pre[k]:8.2f}   post {post[k]:8.2f}   d {post[k]-pre[k]:+7.2f}")
        table.append(dict(name=name, **{f"pre_{k}": pre[k] for k in pre}, **{f"post_{k}": post[k] for k in post}))
    # write table
    outcsv = Path(a.out + ".csv"); outcsv.parent.mkdir(parents=True, exist_ok=True)
    with open(outcsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(table[0].keys())); w.writeheader(); w.writerows(table)
    # PSA dumbbell (the headline descriptor)
    fig, ax = plt.subplots(figsize=(7, 0.7 * len(table) + 1.5))
    for i, t in enumerate(table):
        ax.plot([t["pre_psa"], t["post_psa"]], [i, i], color="#bbb", lw=2, zorder=1)
        ax.scatter(t["pre_psa"], i, s=90, color="#c7551a", zorder=3, label="pre (GFN2)" if i == 0 else "")
        ax.scatter(t["post_psa"], i, s=90, color="#0d8b96", zorder=3, label="post (r2SCAN-3c)" if i == 0 else "")
    ax.set_yticks(range(len(table))); ax.set_yticklabels([t["name"] for t in table])
    ax.set_xlabel("population-weighted 3D-PSA  (A^2)"); ax.set_title("Refinement effect on PSA (pre GFN2 -> post r2SCAN-3c)")
    ax.legend(loc="best", fontsize=9); plt.tight_layout(); plt.savefig(a.out + ".png", dpi=170)
    print(f"\ntable -> {outcsv}\nfigure -> {a.out}.png")

if __name__ == "__main__":
    main()

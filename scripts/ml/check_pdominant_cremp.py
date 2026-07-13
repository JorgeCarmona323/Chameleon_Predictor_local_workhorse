# env: chameleon-calc
"""
check_pdominant_cremp.py
------------------------
Two sanity checks on the low dominant-conformer population (p_dominant) seen for the
DOPC R/S ensembles, especially 3-12-10-12 (~0.12):

PART A — over-fragmentation check (our ensembles):
  Is a low p_dominant a real diffuse ensemble, or CREST splitting one basin into near-
  duplicates? Cluster the low-energy conformers by heavy-atom RMSD (Butina, 1.0 Å) and
  compare the per-conformer p_dominant to the CLUSTER-level p_dominant (Boltzmann pop of
  the largest cluster). If cluster p_dominant >> conformer p_dominant → fragmentation.
  If they track → genuinely diffuse.

PART B — CREMP 6-mer reference:
  Is a low p_dominant normal for 6-mer macrocycles of this size? CREMP ensembles use the
  same CREST/GFN2 protocol. Read p_dominant / n_eff for all CREMP hexamers (filename =
  dot-separated residues) and report the distribution + where our compounds fall.

Usage:  python scripts/check_pdominant_cremp.py
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina

RDLogger.DisableLog("rdApp.*")

RT = 0.592  # kcal/mol at 298 K
HARTREE_KCAL = 627.509
CONF = Path("results/conformers")
OURS = {  # water-phase ensembles of the non-diazirine pairs
    "3-12-8-12_R":  CONF / "DOPC 3-12-8-12 R" / "water",
    "3-12-8-12_S":  CONF / "DOPC 3-12-8-12 S" / "water",
    "3-12-10-12_R": CONF / "3-12-10-12 R" / "water",
    "3-12-10-12_S": CONF / "3-12-10-12 S" / "water",
}
CREMP_DIR = Path("dependencies/pickle")


def _pdom_neff(weights):
    w = np.asarray(weights, float)
    w = w[np.isfinite(w) & (w > 0)]
    if w.size == 0:
        return np.nan, np.nan
    w = w / w.sum()
    neff = float(np.exp(-(w * np.log(w)).sum()))
    return float(w.max()), neff


# ── PART A: fragmentation check ───────────────────────────────────────────────
def fragmentation(label, wdir, ewin=3.0, cap=80, rms_thresh=1.0):
    jp, sp = wdir / "ensemble.json", wdir / "ensemble.sdf"
    confs = json.load(open(jp)).get("conformers", [])
    e = np.array([c.get("totalenergy", np.nan) for c in confs], float) * HARTREE_KCAL
    w = np.array([c.get("boltzmannweight", np.nan) for c in confs], float)
    mols = [m for m in Chem.SDMolSupplier(str(sp), removeHs=False, sanitize=True) if m]
    n = min(len(mols), len(e))
    e, w, mols = e[:n], w[:n], mols[:n]
    p_conf, neff = _pdom_neff(w)
    erel = e - np.nanmin(e)

    # subset = highest-weight (lowest-energy) conformers, for basin clustering
    sub = [int(i) for i in np.argsort(w)[::-1][:cap]]
    base = Chem.Mol(mols[sub[0]]); base.RemoveAllConformers()
    for i in sub:
        base.AddConformer(Chem.Conformer(mols[i].GetConformer()), assignId=True)
    heavy = [a.GetIdx() for a in base.GetAtoms() if a.GetAtomicNum() > 1]
    rms = AllChem.GetConformerRMSMatrix(base, atomIds=heavy, prealigned=False)
    clusters = Butina.ClusterData(rms, len(sub), rms_thresh, isDistData=True, reordering=True)
    # cluster-level p_dominant over the subset's renormalized weights
    wsub = w[sub] / np.nansum(w[sub])
    cl_pops = sorted((float(np.nansum(wsub[list(c)])) for c in clusters), reverse=True)
    return dict(n=n, p_conf=p_conf, neff=neff,
                within_RT=int((erel <= RT).sum()), within_2=int((erel <= 2).sum()),
                n_sub=len(sub), n_clusters=len(clusters), p_cluster=cl_pops[0])


# ── PART B: CREMP hexamer reference ───────────────────────────────────────────
def cremp_hexamers(sample=600):
    files = sorted(CREMP_DIR.glob("*.pickle"))
    hexes = [f for f in files if f.stem.count(".") + 1 == 6]
    if len(hexes) > sample:
        step = len(hexes) // sample
        hexes = hexes[::step][:sample]
    rows = []
    for f in hexes:
        try:
            obj = pickle.load(open(f, "rb"))
            cm = obj.get("conformers", [])
            w = [c.get("boltzmannweight") for c in cm if c.get("boltzmannweight") is not None]
            if not w:
                continue
            pdom, neff = _pdom_neff(w)
            rows.append((pdom, neff, len(w)))
        except Exception:
            continue
    return np.array(rows), len(hexes)


def cremp_basin(sample=40, cap=80, rms_thresh=1.0):
    """Basin-level p_dominant for CREMP hexamers: same protocol as our ensembles
    (top-cap conformers by weight, Butina @ rms_thresh, largest-cluster Boltzmann pop)."""
    files = sorted(CREMP_DIR.glob("*.pickle"))
    hexes = [f for f in files if f.stem.count(".") + 1 == 6]
    if len(hexes) > sample:
        hexes = hexes[::max(1, len(hexes) // sample)][:sample]
    out = []
    for f in hexes:
        try:
            obj = pickle.load(open(f, "rb"))
            mol, cm = obj.get("rd_mol"), obj.get("conformers", [])
            if mol is None or not cm:
                continue
            w = np.array([c.get("boltzmannweight", np.nan) for c in cm], float)
            confs = list(mol.GetConformers())
            n = min(len(confs), len(w))
            if n < 3:
                continue
            top = [int(i) for i in np.argsort(w[:n])[::-1][:cap]]
            base = Chem.Mol(mol); base.RemoveAllConformers()
            for i in top:
                base.AddConformer(Chem.Conformer(confs[i]), assignId=True)
            heavy = [a.GetIdx() for a in base.GetAtoms() if a.GetAtomicNum() > 1]
            rms = AllChem.GetConformerRMSMatrix(base, atomIds=heavy, prealigned=False)
            cl = Butina.ClusterData(rms, len(top), rms_thresh, isDistData=True, reordering=True)
            wt = w[top] / np.nansum(w[top])
            out.append(max(float(np.nansum(wt[list(c)])) for c in cl))
        except Exception:
            continue
    return np.array(out)


def pct(arr, v):
    a = arr[np.isfinite(arr)]
    return 100.0 * float((a < v).mean())


def main():
    print("=" * 70)
    print("PART A — over-fragmentation check (water ensembles, Butina @ 1.0 A RMSD)")
    print("=" * 70)
    print(f"{'compound':14s} {'n':>4} {'p_conf':>7} {'p_clust':>8} {'n_clust':>8} "
          f"{'n_eff':>6} {'≤RT':>4} {'≤2kcal':>7}")
    ours_pdom, ours_basin = {}, {}
    for label, wdir in OURS.items():
        if not (wdir / "ensemble.json").exists():
            print(f"{label:14s}  (missing)"); continue
        r = fragmentation(label, wdir)
        ours_pdom[label] = (r["p_conf"], r["neff"], r["n"])
        ours_basin[label] = r["p_cluster"]
        print(f"{label:14s} {r['n']:>4} {r['p_conf']:>7.3f} {r['p_cluster']:>8.3f} "
              f"{r['n_clusters']:>8} {r['neff']:>6.1f} {r['within_RT']:>4} {r['within_2']:>7}")
    print("\n  p_conf  = per-conformer dominant population (the reported p_dominant)")
    print("  p_clust = Boltzmann pop of the largest 1.0-A RMSD cluster")
    print("  if p_clust >> p_conf → fragmentation; if similar → genuinely diffuse")

    print("\n" + "=" * 70)
    print("PART B — CREMP hexamer reference (same CREST/GFN2 protocol)")
    print("=" * 70)
    data, n_hex = cremp_hexamers()
    if data.size == 0:
        print("  No CREMP hexamer data parsed."); return
    pdom, neff, nconf = data[:, 0], data[:, 1], data[:, 2]
    q = lambda a, p: float(np.nanpercentile(a, p))
    print(f"  CREMP hexamers parsed: {len(pdom)} (of {n_hex} found)")
    print(f"  p_dominant : median {q(pdom,50):.3f}  IQR [{q(pdom,25):.3f}, {q(pdom,75):.3f}]"
          f"  (10th–90th: {q(pdom,10):.3f}–{q(pdom,90):.3f})")
    print(f"  n_eff      : median {q(neff,50):.1f}   IQR [{q(neff,25):.1f}, {q(neff,75):.1f}]")
    print(f"  n_confs    : median {q(nconf,50):.0f}")
    print(f"  fraction of CREMP hexamers with p_dominant < 0.20 : {100*float((pdom<0.20).mean()):.0f}%")
    print(f"  fraction with p_dominant < 0.15                   : {100*float((pdom<0.15).mean()):.0f}%")
    print("\n  Where our compounds fall in the CREMP hexamer per-conformer p_dominant distribution:")
    for label, (pc, ne, nn) in ours_pdom.items():
        print(f"    {label:14s} p_dominant={pc:.3f}  →  {pct(pdom, pc):4.0f}th percentile  (n_eff={ne:.0f})")
    print("  (CONFOUNDED: CREMP ensembles have ~4x more conformers than ours → see Part C)")

    print("\n" + "=" * 70)
    print("PART C — CREMP hexamer BASIN-level p_dominant (1.0 A RMSD, SAME protocol as ours)")
    print("=" * 70)
    cb = cremp_basin()
    if cb.size == 0:
        print("  No CREMP basin data parsed."); return
    print(f"  CREMP hexamers clustered: {len(cb)}")
    print(f"  basin p_dominant : median {np.median(cb):.3f}  "
          f"IQR [{np.percentile(cb,25):.3f}, {np.percentile(cb,75):.3f}]  "
          f"(10th-90th: {np.percentile(cb,10):.3f}-{np.percentile(cb,90):.3f})")
    print("\n  Our basin p_dominant vs the CREMP hexamer basin distribution:")
    for label, bp in ours_basin.items():
        print(f"    {label:14s} basin p_dom={bp:.3f}  →  {pct(cb, bp):4.0f}th percentile")


if __name__ == "__main__":
    main()

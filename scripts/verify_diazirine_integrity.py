# env: chameleon-calc
"""
verify_diazirine_integrity.py
-----------------------------
Diazirine integrity check (steps 1-4) for CREST/GFN2 conformer ensembles.

The strained 3-membered C-N=N diazirine could be broken or distorted by CREST
metadynamics (400-500 K) + GFN2-xTB. This screens EVERY conformer of each
diazirine ensemble for gross failure BEFORE any dz finding is trusted:

  1. locate diazirine   : SMARTS [#6]1[#7]=[#7]1  (matched once; atom indices reused
                          across all conformers so a broken bond is still measured)
  2. N=N tiered         : PASS 1.22-1.25 A | WATCH <=1.30 | FAIL any >1.35 (stretched/broken)
  3. no N2 extrusion    : both C-N bonds ~1.45 A; FAIL > 1.60 (ring opened / N2 left)
  4. geometry sanity    : N-C-N apex angle ~49 deg; per-conformer consistency (std)
  5. alkyne (monitor)   : terminal/internal C#C ~1.20 A + C-C#C linearity ~180 deg
                          (monitor-only; not constrained unless it flags an artifact)

Does NOT do DFT. It reports each diazirine ensemble's lowest-energy conformer index so a
representative can be picked if a DFT spot-check is ever wanted.

Verdict per ensemble = worst of the per-motif verdicts: PASS / WATCH / FAIL.
Aggregate CSV -> results/diazirine_integrity.csv.

Usage:
  python verify_diazirine_integrity.py                 # auto-scan results/conformers/*Diazirine*
  python verify_diazirine_integrity.py --run-dir "results/conformers/3-12-8-12 Diazirine R" --name 812dz_R
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

DIAZIRINE_SMARTS = "[#6]1[#7]=[#7]1"   # C-N=N three-membered ring
ALKYNE_SMARTS = "[CX2]#[CX2]"          # C#C (terminal or internal) — monitor-only motif
REF_NN, REF_CN, REF_NCN = 1.23, 1.47, 49.0   # reference 3,3-disubstituted diazirine
REF_CC = 1.20                          # reference alkyne C#C
# Diazirine N=N tiered acceptance (reviewer-set, 2026-06-14):
NN_PASS = (1.22, 1.25)   # PASS band
NN_WATCH_HI = 1.30       # PASS_hi..WATCH_hi = WATCH (holding but drifting)
NN_FAIL = 1.35           # any conformer N=N above this = FAIL
CN_FAIL = 1.60           # C-N above this = ring opened (toward N2 extrusion)
STD_WARN = 0.04          # per-ensemble N=N std (A) above this = inconsistent
# Alkyne (monitor only — not constrained unless this flags a reproducible artifact):
CC_FAIL = 1.30           # C#C above this = distorted
CC_ANGLE_FAIL = 155.0    # C-C#C linearity below this (deg) = bent


def _measure(coords, c, n1, n2):
    d = lambda a, b: float(np.linalg.norm(coords[a] - coords[b]))
    nn, cn1, cn2 = d(n1, n2), d(c, n1), d(c, n2)
    v1, v2 = coords[n1] - coords[c], coords[n2] - coords[c]
    cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9)
    ang = float(np.degrees(np.arccos(np.clip(cos, -1, 1))))
    return nn, cn1, cn2, ang


def _angle(coords, a, b, c) -> float:
    v1, v2 = coords[a] - coords[b], coords[c] - coords[b]
    cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9)
    return float(np.degrees(np.arccos(np.clip(cos, -1, 1))))


def _find_diazirine(mols):
    patt = Chem.MolFromSmarts(DIAZIRINE_SMARTS)
    for m in mols:
        ms = m.GetSubstructMatches(patt)
        if ms:
            return ms[0]                      # (C, N1, N2)
    return None


def _find_alkyne(mols):
    """Return (C1, C2, neighbor1, neighbor2) for an alkyne, or None.
    Neighbors (the non-alkyne atom bonded to each sp carbon) give the C-C#C linearity angle."""
    patt = Chem.MolFromSmarts(ALKYNE_SMARTS)
    for m in mols:
        ms = m.GetSubstructMatches(patt)
        if ms:
            c1, c2 = ms[0]
            nb1 = next((a.GetIdx() for a in m.GetAtomWithIdx(c1).GetNeighbors() if a.GetIdx() != c2), None)
            nb2 = next((a.GetIdx() for a in m.GetAtomWithIdx(c2).GetNeighbors() if a.GetIdx() != c1), None)
            return c1, c2, nb1, nb2
    return None


def _load(sdf: Path, jsonp: Path):
    energies = weights = None
    if jsonp.exists():
        with open(jsonp) as f:
            confs = json.load(f).get("conformers", [])
        energies = np.array([c.get("totalenergy", np.nan) for c in confs], float)
        weights = np.array([c.get("boltzmannweight", np.nan) for c in confs], float)
    mols = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=True) if m is not None]
    return mols, energies, weights


def check_ensemble(sdf: Path, jsonp: Path, label: str) -> dict | None:
    mols, energies, weights = _load(sdf, jsonp)
    if not mols:
        return {"label": label, "verdict": "NO_DATA", "n_conf": 0}

    dz = _find_diazirine(mols)        # (C, N1, N2) or None
    alk = _find_alkyne(mols)          # (C1, C2, nb1, nb2) or None
    if dz is None and alk is None:
        return {"label": label, "verdict": "NO_MOTIF", "n_conf": len(mols)}

    nn_l, cn_l, ncn_l, nn_broken = [], [], [], 0
    cc_l, cca_l = [], []
    for m in mols:
        coords = m.GetConformer().GetPositions()
        if dz is not None:
            c, n1, n2 = dz
            nn, cn1, cn2, ang = _measure(coords, c, n1, n2)
            nn_l.append(nn); cn_l.append(max(cn1, cn2)); ncn_l.append(ang)
            if nn > NN_FAIL or cn1 > CN_FAIL or cn2 > CN_FAIL:
                nn_broken += 1
        if alk is not None:
            c1, c2, nb1, nb2 = alk
            cc_l.append(float(np.linalg.norm(coords[c1] - coords[c2])))
            angs = [a for a in (
                _angle(coords, nb1, c1, c2) if nb1 is not None else None,
                _angle(coords, c1, c2, nb2) if nb2 is not None else None) if a is not None]
            cca_l.append(min(angs) if angs else float("nan"))

    row = {"label": label, "n_conf": len(mols)}
    verdicts = []

    if dz is not None:
        nn_a, cn_a, ncn_a = np.array(nn_l), np.array(cn_l), np.array(ncn_l)
        nn_mean, nn_std = float(nn_a.mean()), float(nn_a.std())
        if nn_broken > 0:
            v = "FAIL"
        elif NN_PASS[0] <= nn_mean <= NN_PASS[1] and nn_std <= STD_WARN and cn_a.max() <= 1.55:
            v = "PASS"
        elif nn_mean <= NN_WATCH_HI:
            v = "WATCH"
        else:
            v = "FAIL"
        verdicts.append(v)
        lo = (int(np.nanargmin(energies[:len(mols)]))
              if energies is not None and np.isfinite(energies).any() else 0)
        row.update({
            "nn_mean": round(nn_mean, 3), "nn_min": round(float(nn_a.min()), 3),
            "nn_max": round(float(nn_a.max()), 3), "nn_std": round(nn_std, 3),
            "cn_max": round(float(cn_a.max()), 3), "ncn_ang_mean": round(float(ncn_a.mean()), 1),
            "nn_broken": nn_broken, "lowE_idx": lo, "lowE_nn": round(nn_l[lo], 3),
        })

    if alk is not None:
        cc_a, cca_a = np.array(cc_l), np.array(cca_l)
        cc_max, cc_ang_min = float(np.nanmax(cc_a)), float(np.nanmin(cca_a))
        verdicts.append("FAIL" if (cc_max > CC_FAIL or cc_ang_min < CC_ANGLE_FAIL) else "PASS")
        row.update({
            "cc_mean": round(float(np.nanmean(cc_a)), 3), "cc_max": round(cc_max, 3),
            "cc_ang_min": round(cc_ang_min, 1),
        })

    order = {"PASS": 0, "WATCH": 1, "FAIL": 2}
    row["verdict"] = max(verdicts, key=lambda x: order.get(x, 0))
    return row


def _jobs(args):
    if args.run_dir:
        return [(args.name[i] if i < len(args.name) else Path(rd).name, Path(rd))
                for i, rd in enumerate(args.run_dir)]
    base = Path("results/conformers")
    return [(p.name, p) for p in sorted(base.glob("*Diazirine*"))
            if (p / "water").exists() or (p / "mem").exists()]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", action="append", default=[])
    ap.add_argument("--name", action="append", default=[])
    ap.add_argument("-o", "--out", default="results/diazirine_integrity.csv", type=Path)
    args = ap.parse_args()

    jobs = _jobs(args)
    if not jobs:
        print("No diazirine ensembles found under results/conformers/*Diazirine*")
        return

    rows = []
    for name, rdir in jobs:
        for solvent in ("water", "mem"):
            sdf, jsonp = rdir / solvent / "ensemble.sdf", rdir / solvent / "ensemble.json"
            if not sdf.exists():
                continue
            res = check_ensemble(sdf, jsonp, f"{name}/{solvent}")
            if res:
                rows.append(res)

    if not rows:
        print("No ensemble.sdf files found in the diazirine run dirs.")
        return

    df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)

    print(f"Reference: diazirine N=N {REF_NN} A / C-N {REF_CN} A / N-C-N {REF_NCN} deg; "
          f"alkyne C#C {REF_CC} A")
    print(f"Diazirine N=N: PASS {NN_PASS[0]}-{NN_PASS[1]} | WATCH <= {NN_WATCH_HI} | "
          f"FAIL any > {NN_FAIL} or C-N > {CN_FAIL}")
    print(f"Alkyne (monitor): FAIL if C#C > {CC_FAIL} or C-C#C angle < {CC_ANGLE_FAIL} deg\n")
    cols = [c for c in ["label", "verdict", "n_conf", "nn_mean", "nn_min", "nn_max", "nn_std",
                        "cn_max", "ncn_ang_mean", "nn_broken", "cc_mean", "cc_ang_min",
                        "lowE_idx", "lowE_nn"] if c in df.columns]
    print(df[cols].to_string(index=False))

    counts = df["verdict"].value_counts().to_dict()
    print(f"\nSUMMARY: {len(df)} ensembles | "
          f"{counts.get('PASS',0)} PASS, {counts.get('WATCH',0)} WATCH, {counts.get('FAIL',0)} FAIL"
          f"  -> {args.out}")
    flagged = df[~df["verdict"].isin(["PASS"])]["label"].tolist()
    if flagged:
        print(f"Review (not PASS): {flagged}")
    if "lowE_idx" in df.columns:
        print("\nIf a DFT spot-check is ever wanted, use each ensemble's lowE_idx conformer.")


if __name__ == "__main__":
    main()

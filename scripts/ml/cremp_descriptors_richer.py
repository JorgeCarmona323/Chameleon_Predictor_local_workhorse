# env: chameleon-calc
"""
cremp_descriptors_richer.py
---------------------------
Compute the VALIDATED-CORE 3D descriptors over the CREMP chloroform ensembles —
the set our descriptor review actually endorsed (apolar 3D-PSA + radius of
gyration + backbone-transannular IMHB + shape), rather than the PSA-only feature
set in the old cremp_deltapsa.csv.

ENERGY-FREE BY DESIGN. Every descriptor is an UNWEIGHTED ensemble statistic
(mean / min / max / spread over sampled conformers). No Boltzmann weighting, no
ensemble energy, no lowest-population — those are on hold until the energy reruns,
so nothing here depends on the (untrusted) per-conformer energies.

Reuses the exact definitions from phys_descriptors_v3.py (same PSA + IMHB the
DOPC R/S reports use); Rgyr / NPR / asphericity come straight from RDKit.

NOTE: CREMP is chloroform-only, so there is still no true ΔG_transfer / dual-solvent
ΔPSA here — that gap is structural to CREMP and flagged in the writeup.

Usage:
  python scripts/cremp_descriptors_richer.py \
      --pickledir dependencies/pickle \
      --out results/2026-07-07_cremp_descriptors_richer.csv [--max-confs 100]
"""
import argparse
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # scripts/ (for phys_descriptors_v3)
from phys_descriptors_v3 import (
    surface_descriptors_mol,
    imhb_descriptors_mol,
    backbone_hbond_atoms,
)

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")


def _sample_conf_ids(mol, cap):
    ids = [c.GetId() for c in mol.GetConformers()]
    if len(ids) <= cap:
        return ids
    step = max(1, len(ids) // cap)
    return ids[::step][:cap]


def _stats(a, prefix):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {f"{prefix}_mean": np.nan, f"{prefix}_min": np.nan,
                f"{prefix}_max": np.nan, f"{prefix}_spread": np.nan}
    return {f"{prefix}_mean": round(float(a.mean()), 4),
            f"{prefix}_min": round(float(a.min()), 4),
            f"{prefix}_max": round(float(a.max()), 4),
            f"{prefix}_spread": round(float(a.max() - a.min()), 4)}


def process(name, obj, cap):
    cid = Path(name).stem
    mol = obj.get("rd_mol")
    smiles = obj.get("smiles", "")
    if mol is None or mol.GetNumConformers() == 0:
        return {"compound_id": cid, "smiles": smiles, "error": "no_rdmol"}

    n_confs = mol.GetNumConformers()
    conf_ids = _sample_conf_ids(mol, cap)
    bb = backbone_hbond_atoms(mol)  # once per compound

    psa, rg, imhb_tot, imhb_bb, npr1, npr2, asph = [], [], [], [], [], [], []
    for i in conf_ids:
        s = surface_descriptors_mol(mol, i)
        psa.append(s["psa"])
        h = imhb_descriptors_mol(mol, i, bb)
        imhb_tot.append(h["imhb"])
        imhb_bb.append(h["imhb_bb"])
        try:
            rg.append(rdMolDescriptors.CalcRadiusOfGyration(mol, confId=i))
            npr1.append(rdMolDescriptors.CalcNPR1(mol, confId=i))
            npr2.append(rdMolDescriptors.CalcNPR2(mol, confId=i))
            asph.append(rdMolDescriptors.CalcAsphericity(mol, confId=i))
        except Exception:
            rg.append(np.nan); npr1.append(np.nan)
            npr2.append(np.nan); asph.append(np.nan)

    row = {"compound_id": cid, "smiles": smiles, "n_confs": n_confs,
           "n_sampled": len(conf_ids), "error": None}
    row.update(_stats(psa, "psa"))       # apolar 3D-PSA (Ono/Begnini def)
    row.update(_stats(rg, "rg"))         # radius of gyration
    # IMHB + shape: unweighted means (counts/ratios — spread less meaningful)
    def _m(a):
        a = np.asarray(a, float); a = a[np.isfinite(a)]
        return round(float(a.mean()), 4) if a.size else np.nan
    row["imhb_total_mean"] = _m(imhb_tot)
    row["imhb_bb_mean"] = _m(imhb_bb)          # backbone transannular
    row["imhb_total_max"] = (int(np.nanmax(imhb_tot)) if np.isfinite(np.nanmax(imhb_tot)) else np.nan)
    row["npr1_mean"] = _m(npr1)
    row["npr2_mean"] = _m(npr2)
    row["asphericity_mean"] = _m(asph)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pickledir", default="dependencies/pickle")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-confs", type=int, default=100)
    args = ap.parse_args()

    files = sorted(Path(args.pickledir).glob("*.pickle"))
    print(f"Found {len(files)} CREMP pickles in {args.pickledir}")

    rows, failed = [], 0
    for f in tqdm(files, desc="CREMP richer descriptors"):
        try:
            with open(f, "rb") as fh:
                obj = pickle.load(fh)
            rows.append(process(f.name, obj, args.max_confs))
        except Exception as e:
            rows.append({"compound_id": f.stem, "smiles": "", "error": str(e)})
            failed += 1

    df = pd.DataFrame(rows)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    n_ok = df["error"].isna().sum()
    print(f"\nSaved {out}  |  ok={n_ok}/{len(df)}  failed={failed}")
    if n_ok:
        ok = df[df["error"].isna()]
        print("psa_mean:", ok["psa_mean"].describe().round(1).to_dict())
        print("rg_mean :", ok["rg_mean"].describe().round(2).to_dict())
        print("imhb_bb_mean:", ok["imhb_bb_mean"].describe().round(2).to_dict())


if __name__ == "__main__":
    main()

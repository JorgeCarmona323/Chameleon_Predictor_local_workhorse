# env: chameleon-ml
"""
dump_feature_matrices.py
------------------------
Stage 1 of the decoupled TabPFN benchmark. Compute the F1-F7 feature matrices +
labels EXACTLY as feature_benchmark.py does (reuses its generators + merge), then
dump them to a single .npz.

Why separate: TabPFN needs torch, and on Windows torch intermittently fails to
load its DLLs (WinError 127) when co-imported with the mkl-heavy rdkit/mordred
stack inside one long-running process. So we generate features here (rdkit/mordred,
no torch) and score TabPFN in a torch-clean process (tabpfn_score.py) that reads
this .npz. RF/LightGBM already ran inside feature_benchmark.py — this is only to
feed the offline TabPFN lane, so results stay identical (same seed, same splits).

Usage:
  python scripts/dump_feature_matrices.py \
      --cremp results/archive/cremp_deltapsa.csv \
      --matrix results/archive/feature_matrix.csv \
      --out <scratch>/feature_matrices.npz
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parent))
import feature_benchmark as fb  # reuse generators, constants, F7 columns


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cremp", default="results/archive/cremp_deltapsa.csv")
    ap.add_argument("--matrix", default="results/archive/feature_matrix.csv")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # ── Load + merge: identical to feature_benchmark.main() ──────────────────
    cremp = pd.read_csv(args.cremp)
    cremp = cremp[cremp["error"].isna()].copy()
    fm = pd.read_csv(args.matrix)

    perm_col = next((c for c in fm.columns if c.lower() == "permeability"), None)
    if perm_col is None:
        perm_col = next((c for c in fm.columns if c.lower() == "pampa"), None)
    smiles_col = next((c for c in fm.columns if c.lower() == "smiles"), None)

    fm_labeled = fm[[smiles_col, perm_col, "Source"]].dropna(subset=[perm_col]).copy()
    fm_labeled["permeable"] = (fm_labeled[perm_col] >= fb.PAMPA_THRESHOLD).astype(int)

    def _canonical(smi):
        try:
            mol = Chem.MolFromSmiles(str(smi))
            return Chem.MolToSmiles(mol) if mol else None
        except Exception:
            return None

    cremp["canon_smiles"] = cremp["smiles"].apply(_canonical)
    fm_labeled["canon_smiles"] = fm_labeled[smiles_col].apply(_canonical)
    cremp = cremp.dropna(subset=["canon_smiles"])
    fm_labeled = fm_labeled.dropna(subset=["canon_smiles"])

    merged = cremp.merge(fm_labeled, on="canon_smiles", how="inner")
    merged = merged.drop_duplicates(subset="canon_smiles").copy()
    print(f"  Merged (F7 + labels): {len(merged)} compounds")
    print(f"  Permeable: {merged['permeable'].sum()} ({merged['permeable'].mean()*100:.1f}%)")

    smiles = merged["canon_smiles"].tolist()
    y = merged["permeable"].values.astype(np.int64)
    sources = merged["Source"].values.astype(str)

    # ── Build the same feature sets ─────────────────────────────────────────
    feats = {}
    print("  F1: Morgan bit 2048");    feats["F1_morgan_bit_2048"]   = fb.morgan_bits(smiles, fb.MORGAN_BIT_DIM)
    print("  F2: Morgan count 500");   feats["F2_morgan_count_500"]  = fb.morgan_counts(smiles, fb.MORGAN_COUNT_DIM_SMALL)
    print("  F3: Morgan count 2048");  feats["F3_morgan_count_2048"] = fb.morgan_counts(smiles, fb.MORGAN_COUNT_DIM_LARGE)
    print("  F4: Atom-pair 2048");     feats["F4_mapc_2048"]         = fb.mapc_features(smiles, fb.MAPC_DIM)
    print("  F5: Mordred 2D (slow)");  feats["F5_mordred_2d"]        = fb.mordred_2d(smiles)[0]
    print("  F6: Mordred 2D+3D (slow)"); feats["F6_mordred_2d3d"]    = fb.mordred_2d3d(smiles)[0]
    f7 = [c for c in fb.F7_COLS if c in merged.columns]
    feats["F7_crest_chcl3"] = merged[f7].fillna(0).values.astype(np.float32)
    print(f"  F7: CREST CHCl3 ({len(f7)} cols)")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, y=y, sources=sources, **feats)
    print(f"Saved {out}  ({', '.join(f'{k}{v.shape}' for k, v in feats.items())})")


if __name__ == "__main__":
    main()

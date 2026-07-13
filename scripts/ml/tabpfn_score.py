# env: chameleon-ml  (TORCH-CLEAN process: imports torch first, no rdkit/mordred)
"""
tabpfn_score.py
---------------
Stage 2 of the decoupled TabPFN benchmark. Runs ONLY TabPFN v2, on the feature
matrices dumped by dump_feature_matrices.py, in a process that imports torch
FIRST and never imports rdkit/mordred — sidestepping the Windows shm.dll
(WinError 127) conflict that killed TabPFN inside the monolithic benchmark.

Same subset, same CV as feature_benchmark.py so the AUCs are directly
comparable to its RF/LightGBM rows:
  - random 5-fold stratified CV (seed 42)
  - leave-source-out CV holding out 2020_Townsend
  - wide feature sets (>100 dims) are PCA-reduced to 100 (as build_tabpfn did)

Output schema matches feature_benchmark_results.csv, so the rows can be
concatenated directly.

Usage:
  python scripts/tabpfn_score.py --npz <scratch>/feature_matrices.npz \
      --out results/2026-07-07_tabpfn_results.csv
"""
import torch  # FIRST — initialize torch DLLs before anything else
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from tabpfn import TabPFNClassifier

CV_FOLDS = 5
CV_SEED = 42
HOLDOUT_SOURCE = "2020_Townsend"
MAX_FEAT = 100  # TabPFN feature cap used by feature_benchmark.build_tabpfn


def build_pipe(n_feat: int) -> Pipeline:
    steps = [("imp", SimpleImputer(strategy="median")), ("scl", StandardScaler())]
    if n_feat > MAX_FEAT:
        steps.append(("pca", PCA(n_components=MAX_FEAT, random_state=CV_SEED)))
    steps.append(("clf", TabPFNClassifier(device="cpu")))
    return Pipeline(steps)


def random_cv(X, y):
    skf = StratifiedKFold(n_splits=CV_FOLDS, shuffle=True, random_state=CV_SEED)
    aucs = []
    for tr, te in skf.split(X, y):
        pipe = build_pipe(X.shape[1])
        pipe.fit(X[tr], y[tr])
        aucs.append(roc_auc_score(y[te], pipe.predict_proba(X[te])[:, 1]))
    return float(np.mean(aucs)), float(np.std(aucs))


def source_cv(X, y, sources):
    te = sources == HOLDOUT_SOURCE
    tr = ~te
    if tr.sum() == 0 or te.sum() == 0:
        return float("nan")
    pipe = build_pipe(X.shape[1])
    pipe.fit(X[tr], y[tr])
    return float(roc_auc_score(y[te], pipe.predict_proba(X[te])[:, 1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    print(f"torch {torch.__version__}  cuda={torch.cuda.is_available()}")
    data = np.load(args.npz, allow_pickle=True)
    y = data["y"].astype(int)
    sources = data["sources"].astype(str)
    feat_names = [k for k in data.files if k not in ("y", "sources")]
    print(f"Loaded {len(y)} compounds, {len(feat_names)} feature sets")

    rows = []
    for name in feat_names:
        X = data[name].astype(np.float32)
        print(f"\n{name}  shape={X.shape}")
        print("  random CV...", flush=True)
        m, s = random_cv(X, y)
        print(f"    AUC={m:.3f} ± {s:.3f}")
        rows.append(dict(feature_set=name, model="TabPFN", cv_type="random_5fold",
                         auc_mean=m, auc_std=s, n_features=X.shape[1], n_compounds=len(y)))
        print("  source-stratified (leave-Townsend-out)...", flush=True)
        a = source_cv(X, y, sources)
        print(f"    AUC={a:.3f}")
        rows.append(dict(feature_set=name, model="TabPFN",
                         cv_type="source_stratified_townsend",
                         auc_mean=a, auc_std=np.nan, n_features=X.shape[1],
                         n_compounds=len(y)))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()

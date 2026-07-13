# env: chameleon-ml
"""
feature_benchmark.py
--------------------
Compare F1-F7 feature sets x 3 models on CREMP permeability subset (n=3,258).

Feature sets:
  F1 — Morgan bit-based  (r=2, 2048-dim)
  F2 — Morgan count-based (r=2, 500-dim)
  F3 — Morgan count-based (r=2, 2048-dim)
  F4 — Atom-pair fingerprint (2048-dim, RDKit hashed bit vector; chirality-aware)
  F5 — Mordred 2D only
  F6 — Mordred 2D+3D (single ETKDG conformer)
  F7 — CREST CHCl3 ensemble descriptors (from cremp_deltapsa.csv)

Models: TabPFN, LightGBM, RandomForest

Usage:
  python scripts/feature_benchmark.py
  python scripts/feature_benchmark.py --cremp results/cremp_deltapsa.csv
                                       --matrix results/feature_matrix.csv
                                       --outdir results
"""

import os
import argparse
import warnings
import logging
from pathlib import Path

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"  # Windows OpenMP conflict between torch and rdkit

import numpy as np
import pandas as pd
from rdkit import Chem  # module-level: main()'s _canonical() + F7 merge need it
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from tqdm import tqdm

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.WARNING)

# ── Constants ─────────────────────────────────────────────────────────────────

PAMPA_THRESHOLD  = -6.0          # log cm/s — permeable if >= threshold
CV_FOLDS         = 5
CV_SEED          = 42
MORGAN_RADIUS    = 2
MAPC_DIM         = 2048
MORGAN_BIT_DIM   = 2048
MORGAN_COUNT_DIM_SMALL = 500
MORGAN_COUNT_DIM_LARGE = 2048

# Richer validated-core F7 (2026-07-07): apolar 3D-PSA + radius of gyration +
# backbone-transannular IMHB + shape, from cremp_descriptors_richer.py. All
# ENERGY-FREE unweighted ensemble stats — the 3 energy-weighted columns of the
# old PSA-only set (bw_psa3d, ensemble_energy, pop_lowest_pct) are intentionally
# dropped until the energy reruns land. Old set kept in git history / the
# 2026-07-07 experiment writeup.
F7_COLS = [
    "psa_mean", "psa_min", "psa_max", "psa_spread",
    "rg_mean", "rg_min", "rg_max", "rg_spread",
    "imhb_total_mean", "imhb_bb_mean", "imhb_total_max",
    "npr1_mean", "npr2_mean", "asphericity_mean",
]

# ── Feature generators ────────────────────────────────────────────────────────

def morgan_bits(smiles_list: list[str], n_bits: int = MORGAN_BIT_DIM) -> np.ndarray:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    rows = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rows.append(np.zeros(n_bits))
        else:
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, nBits=n_bits)
            rows.append(np.array(fp))
    return np.array(rows, dtype=np.float32)


def morgan_counts(smiles_list: list[str], n_bits: int = MORGAN_COUNT_DIM_LARGE) -> np.ndarray:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    from rdkit import DataStructs
    rows = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rows.append(np.zeros(n_bits))
        else:
            fp = rdMolDescriptors.GetHashedMorganFingerprint(mol, MORGAN_RADIUS, nBits=n_bits)
            arr = np.zeros(n_bits, dtype=np.float32)
            DataStructs.ConvertToNumpyArray(fp, arr)
            rows.append(arr)
    return np.array(rows, dtype=np.float32)


def mapc_features(smiles_list: list[str], dim: int = MAPC_DIM) -> np.ndarray:
    # Atom-pair fingerprint — chirality-aware, captures long-range pairwise
    # relationships. Conceptually equivalent to MAPC/MAP4 but via RDKit's
    # GetHashedAtomPairFingerprintAsBitVect (avoids mhfp hash encoding issues).
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    rows = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rows.append(np.zeros(dim, dtype=np.float32))
        else:
            fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=dim)
            rows.append(np.array(fp, dtype=np.float32))
    return np.array(rows)


def mordred_2d(smiles_list: list[str]) -> tuple[np.ndarray, list[str]]:
    from rdkit import Chem
    from mordred import Calculator, descriptors as desc
    calc = Calculator(desc, ignore_3D=True)
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    df = calc.pandas(mols, nproc=1, quiet=True)
    df = df.select_dtypes(include=[np.number]).fillna(0).replace([np.inf, -np.inf], 0)
    return df.values.astype(np.float32), list(df.columns)


def mordred_2d3d(smiles_list: list[str]) -> tuple[np.ndarray, list[str]]:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from mordred import Calculator, descriptors as desc
    calc = Calculator(desc, ignore_3D=False)
    mols = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            mol = Chem.AddHs(mol)
            try:
                AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                AllChem.MMFFOptimizeMolecule(mol)
                mol = Chem.RemoveHs(mol)
            except Exception:
                mol = Chem.RemoveHs(mol)
        mols.append(mol)
    df = calc.pandas(mols, nproc=1, quiet=True)
    df = df.select_dtypes(include=[np.number]).fillna(0).replace([np.inf, -np.inf], 0)
    return df.values.astype(np.float32), list(df.columns)


# ── Model builders ────────────────────────────────────────────────────────────

def build_rf():
    return Pipeline([
        ("imp", SimpleImputer(strategy="median")),
        ("scl", StandardScaler()),
        ("clf", RandomForestClassifier(n_estimators=300, random_state=CV_SEED, n_jobs=-1)),
    ])


def build_lgbm():
    from lightgbm import LGBMClassifier
    return Pipeline([
        ("imp", SimpleImputer(strategy="median")),
        ("scl", StandardScaler()),
        ("clf", LGBMClassifier(n_estimators=300, learning_rate=0.05,
                               num_leaves=31, random_state=CV_SEED,
                               verbosity=-1, n_jobs=-1)),
    ])


def build_tabpfn(X_train):
    from tabpfn import TabPFNClassifier
    n_feat = X_train.shape[1]
    max_feat = 100
    if n_feat > max_feat:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=max_feat, random_state=CV_SEED)
        return Pipeline([
            ("imp", SimpleImputer(strategy="median")),
            ("scl", StandardScaler()),
            ("pca", pca),
            ("clf", TabPFNClassifier(device="cpu")),
        ])
    return Pipeline([
        ("imp", SimpleImputer(strategy="median")),
        ("scl", StandardScaler()),
        ("clf", TabPFNClassifier(device="cpu")),
    ])


# ── Cross-validation ──────────────────────────────────────────────────────────

def _build_pipe(model_name: str, X_tr: np.ndarray):
    if model_name == "TabPFN":
        return build_tabpfn(X_tr)
    elif model_name == "LightGBM":
        return build_lgbm()
    return build_rf()


def cross_val_auc(X: np.ndarray, y: np.ndarray, model_name: str) -> dict:
    skf = StratifiedKFold(n_splits=CV_FOLDS, shuffle=True, random_state=CV_SEED)
    aucs = []
    for train_idx, test_idx in skf.split(X, y):
        X_tr, X_te = X[train_idx], X[test_idx]
        y_tr, y_te = y[train_idx], y[test_idx]
        pipe = _build_pipe(model_name, X_tr)
        pipe.fit(X_tr, y_tr)
        prob = pipe.predict_proba(X_te)[:, 1]
        aucs.append(roc_auc_score(y_te, prob))
    return {"mean": float(np.mean(aucs)), "std": float(np.std(aucs)), "folds": aucs}


def source_stratified_auc(X: np.ndarray, y: np.ndarray, sources: np.ndarray,
                           model_name: str, holdout_source: str = "2020_Townsend") -> dict:
    """Train on all non-holdout sources, test on holdout source."""
    test_mask = sources == holdout_source
    train_mask = ~test_mask
    if train_mask.sum() == 0 or test_mask.sum() == 0:
        return {"mean": float("nan"), "std": float("nan")}
    X_tr, X_te = X[train_mask], X[test_mask]
    y_tr, y_te = y[train_mask], y[test_mask]
    pipe = _build_pipe(model_name, X_tr)
    pipe.fit(X_tr, y_tr)
    prob = pipe.predict_proba(X_te)[:, 1]
    auc = roc_auc_score(y_te, prob)
    return {"mean": float(auc), "std": float("nan"), "folds": [auc]}


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cremp",  default="results/cremp_deltapsa.csv")
    parser.add_argument("--matrix", default="results/feature_matrix.csv")
    parser.add_argument("--outdir", default="results")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Load CREMP CREST data (F7) ─────────────────────────────────────────
    print("Loading CREMP delta-PSA data...")
    cremp = pd.read_csv(args.cremp)
    cremp = cremp[cremp["error"].isna()].copy()
    print(f"  CREMP usable: {len(cremp)} compounds")

    # ── Load feature_matrix for permeability labels ────────────────────────
    print("Loading feature matrix for permeability labels...")
    fm = pd.read_csv(args.matrix)

    # find permeability column — prefer exact match first
    perm_col = next((c for c in fm.columns if c.lower() == "permeability"), None)
    if perm_col is None:
        perm_col = next((c for c in fm.columns if c.lower() == "pampa"), None)
    if perm_col is None:
        perm_col = next((c for c in fm.columns if "pampa" in c.lower()), None)
    smiles_col = next((c for c in fm.columns if c.lower() == "smiles"), None)
    if smiles_col is None:
        smiles_col = next((c for c in fm.columns if "smiles" in c.lower()), None)
    if perm_col is None or smiles_col is None:
        raise ValueError(f"Cannot find SMILES or permeability column. Cols: {list(fm.columns[:20])}")
    print(f"  Using permeability col: '{perm_col}', SMILES col: '{smiles_col}'")

    fm_labeled = fm[[smiles_col, perm_col, "Source"]].dropna(subset=[perm_col]).copy()
    fm_labeled["permeable"] = (fm_labeled[perm_col] >= PAMPA_THRESHOLD).astype(int)

    # ── Canonicalize SMILES on both sides before merging ──────────────────
    def _canonical(smi):
        try:
            mol = Chem.MolFromSmiles(str(smi))
            return Chem.MolToSmiles(mol) if mol else None
        except Exception:
            return None

    cremp["canon_smiles"] = cremp["smiles"].apply(_canonical)
    fm_labeled["canon_smiles"] = fm_labeled[smiles_col].apply(_canonical)

    n_cremp_before = len(cremp)
    cremp = cremp.dropna(subset=["canon_smiles"])
    fm_labeled = fm_labeled.dropna(subset=["canon_smiles"])
    print(f"  CREMP after canonicalization: {len(cremp)} (dropped {n_cremp_before - len(cremp)} unparseable)")

    # ── Merge to get compounds with both F7 and labels ─────────────────────
    merged = cremp.merge(fm_labeled, on="canon_smiles", how="inner")

    # Deduplicate on canonical SMILES — keep first occurrence
    n_before = len(merged)
    merged = merged.drop_duplicates(subset="canon_smiles").copy()
    print(f"  Merged (F7 + labels): {len(merged)} compounds ({n_before - len(merged)} duplicates removed)")
    print(f"  Permeable: {merged['permeable'].sum()} ({merged['permeable'].mean()*100:.1f}%)")
    print(f"  Source breakdown: {merged['Source'].value_counts().to_dict()}")
    print(f"\n  WARNING: 2020_Townsend = {(merged['Source']=='2020_Townsend').sum()} / {len(merged)} compounds ({(merged['Source']=='2020_Townsend').mean()*100:.1f}%)")
    print(f"  Random CV may be optimistic — also running source-stratified CV (leave-source-out on Townsend)\n")

    smiles = merged["canon_smiles"].tolist()
    y = merged["permeable"].values
    sources = merged["Source"].values

    # ── Build feature sets ─────────────────────────────────────────────────
    print("\nGenerating features...")

    feature_sets = {}

    print("  F1: Morgan bit-based 2048-dim")
    feature_sets["F1_morgan_bit_2048"] = morgan_bits(smiles, MORGAN_BIT_DIM)

    print("  F2: Morgan count-based 500-dim")
    feature_sets["F2_morgan_count_500"] = morgan_counts(smiles, MORGAN_COUNT_DIM_SMALL)

    print("  F3: Morgan count-based 2048-dim")
    feature_sets["F3_morgan_count_2048"] = morgan_counts(smiles, MORGAN_COUNT_DIM_LARGE)

    print("  F4: Atom-pair fingerprint 2048-dim (chirality-aware)")
    feature_sets["F4_mapc_2048"] = mapc_features(smiles, MAPC_DIM)

    print("  F5: Mordred 2D only (this may take a few minutes)...")
    X_f5, _ = mordred_2d(smiles)
    feature_sets["F5_mordred_2d"] = X_f5

    print("  F6: Mordred 2D+3D single ETKDG conformer (this may take several minutes)...")
    X_f6, _ = mordred_2d3d(smiles)
    feature_sets["F6_mordred_2d3d"] = X_f6

    print("  F7: CREST CHCl3 ensemble descriptors")
    f7_available = [c for c in F7_COLS if c in merged.columns]
    feature_sets["F7_crest_chcl3"] = merged[f7_available].fillna(0).values.astype(np.float32)
    print(f"    Using {len(f7_available)} CREST features: {f7_available}")

    # ── Run benchmark ──────────────────────────────────────────────────────
    models = ["RandomForest", "LightGBM", "TabPFN"]
    rows = []

    for feat_name, X in feature_sets.items():
        print(f"\n{'─'*60}")
        print(f"  {feat_name}  shape={X.shape}")
        for model_name in models:
            # Random stratified CV
            print(f"    {model_name} [random CV]...", end=" ", flush=True)
            try:
                result = cross_val_auc(X, y, model_name)
                print(f"AUC={result['mean']:.3f} ± {result['std']:.3f}")
                rows.append({
                    "feature_set": feat_name,
                    "model": model_name,
                    "cv_type": "random_5fold",
                    "auc_mean": result["mean"],
                    "auc_std": result["std"],
                    "n_features": X.shape[1],
                    "n_compounds": len(y),
                })
            except Exception as e:
                import traceback
                print(f"FAILED: {e}")
                traceback.print_exc()
                rows.append({
                    "feature_set": feat_name,
                    "model": model_name,
                    "cv_type": "random_5fold",
                    "auc_mean": np.nan,
                    "auc_std": np.nan,
                    "n_features": X.shape[1],
                    "n_compounds": len(y),
                })

            # Source-stratified CV (train on non-Townsend, test on Townsend)
            print(f"    {model_name} [source-stratified]...", end=" ", flush=True)
            try:
                result_ss = source_stratified_auc(X, y, sources, model_name)
                print(f"AUC={result_ss['mean']:.3f}")
                rows.append({
                    "feature_set": feat_name,
                    "model": model_name,
                    "cv_type": "source_stratified_townsend",
                    "auc_mean": result_ss["mean"],
                    "auc_std": result_ss["std"],
                    "n_features": X.shape[1],
                    "n_compounds": len(y),
                })
            except Exception as e:
                import traceback
                print(f"FAILED: {e}")
                traceback.print_exc()
                rows.append({
                    "feature_set": feat_name,
                    "model": model_name,
                    "cv_type": "source_stratified_townsend",
                    "auc_mean": np.nan,
                    "auc_std": np.nan,
                    "n_features": X.shape[1],
                    "n_compounds": len(y),
                })

    # ── Save results ───────────────────────────────────────────────────────
    results_df = pd.DataFrame(rows)
    out_path = outdir / "feature_benchmark_results.csv"
    results_df.to_csv(out_path, index=False)
    print(f"\n{'='*60}")
    print(f"Results saved to {out_path}")
    for cv_type in results_df["cv_type"].unique():
        subset = results_df[results_df["cv_type"] == cv_type]
        print(f"\n--- {cv_type} ---")
        try:
            print(subset.pivot(index="feature_set", columns="model", values="auc_mean").round(3).to_string())
        except Exception:
            print(subset[["feature_set","model","auc_mean"]].to_string())


if __name__ == "__main__":
    main()

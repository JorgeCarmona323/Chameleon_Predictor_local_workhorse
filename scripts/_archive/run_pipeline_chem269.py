"""
run_pipeline.py — Single-entry-point for the full Tier-1 pipeline.

Steps:
  1. curate_data          — filter PAMPA subset, curate reference set
  2. conformer_engine     — ETKDG + MMFF94s Δ descriptors (Tier-1)
  3. build_feature_matrix — merge all features
  4. correlation_analysis — Pearson/Spearman r, AUC-ROC, LR importance
  5. umap_visualization   — UMAP panels A/B/C with Leiden + LogPexp coloring
  6. tier2_validation     — cross-check on reference set

Usage:
  python run_pipeline.py                          # full run
  python run_pipeline.py --max-mols 100           # quick test (100 molecules)
  python run_pipeline.py --skip-conformers        # use DB 3DPSA only (fastest)

All outputs go to results/ and figures/.
"""

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent
DATA_CSV = ROOT / "CycPeptMPDB_Peptide_All (2).csv"
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
SCRIPTS = ROOT / "scripts"


def run_step(name: str, cmd: list[str]) -> None:
    print(f"\n{'='*70}")
    print(f"STEP: {name}")
    print(f"{'='*70}")
    result = subprocess.run(cmd, cwd=ROOT)
    if result.returncode != 0:
        print(f"\n❌ Step '{name}' failed (exit code {result.returncode})")
        sys.exit(result.returncode)
    print(f"✓ {name} complete")


def main(args: argparse.Namespace) -> None:
    print("=" * 70)
    print("CHEM 269 Final Project — Cyclic Peptide Dual-Dielectric Pipeline")
    print("Jorge Carmona | March 2026")
    print("=" * 70)

    py = sys.executable

    # Step 1: Data curation
    run_step("Data curation", [
        py, str(SCRIPTS / "curate_data.py"),
        "--input", str(DATA_CSV),
        "--outdir", str(DATA_DIR),
    ])

    # Step 2: Conformer engine (Tier-1)
    if not args.skip_conformers:
        n_confs = str(args.n_confs)
        max_mols = str(args.max_mols)
        n_cpus = str(args.n_cpus)
        run_step("Conformer generation (Tier-1)", [
            py, str(SCRIPTS / "conformer_engine.py"),
            "--input", str(DATA_DIR / "pampa_curated.csv"),
            "--outdir", str(RESULTS_DIR),
            "--n-confs", n_confs,
            "--n-cpus", n_cpus,
            "--max-mols", max_mols,
        ])
    else:
        print("\nSkipping conformer generation (--skip-conformers flag set)")
        print("Using DB 3DPSA values (delta_3DPSA_db) as primary 3D feature")

    # Step 3: Feature matrix
    # Only pass --conformers if the file exists (i.e., Step 2 was not skipped)
    conformer_csv = RESULTS_DIR / "conformer_descriptors_raw.csv"
    conformer_args = ["--conformers", str(conformer_csv)] if conformer_csv.exists() else []
    if not conformer_args:
        print("  Note: Tier-1 conformer file not found — feature matrix will use DB 3DPSA + 2D descriptors only.")
        print("  To generate Tier-1 Δ features, re-run without --skip-conformers.")
    run_step("Build feature matrix", [
        py, str(SCRIPTS / "build_feature_matrix.py"),
        "--pampa", str(DATA_DIR / "pampa_curated.csv"),
        "--outdir", str(RESULTS_DIR),
        *conformer_args,
    ])

    # Step 4: Correlation analysis
    run_step("Correlation analysis", [
        py, str(SCRIPTS / "correlation_analysis.py"),
        "--matrix", str(RESULTS_DIR / "feature_matrix.csv"),
        "--outdir", str(RESULTS_DIR),
    ])

    # Step 5: UMAP visualization
    run_step("UMAP visualization", [
        py, str(SCRIPTS / "umap_visualization.py"),
        "--matrix", str(RESULTS_DIR / "feature_matrix.csv"),
        "--outdir", str(RESULTS_DIR),
    ])

    # Step 6: Tier-2 validation
    ref_csv = DATA_DIR / "reference_set.csv"
    if ref_csv.exists():
        run_step("Tier-2 validation (reference set cross-check)", [
            py, str(SCRIPTS / "tier2_validation.py"),
            "--matrix", str(RESULTS_DIR / "feature_matrix.csv"),
            "--refset", str(ref_csv),
            "--outdir", str(RESULTS_DIR),
        ])

    print("\n" + "=" * 70)
    print("✅ Pipeline complete!")
    print(f"   Results: {RESULTS_DIR}/")
    print(f"   Figures: {RESULTS_DIR}/figures/")
    print("=" * 70)


def parse_args() -> argparse.Namespace:
    import multiprocessing
    parser = argparse.ArgumentParser(description="Full CHEM269 cyclic peptide pipeline")
    parser.add_argument("--max-mols", "-n", type=int, default=0,
                        help="Limit to N molecules for testing (0 = all)")
    parser.add_argument("--n-confs", "-c", type=int, default=50,
                        help="Conformers per molecule (default=50; use 200 for final)")
    parser.add_argument("--n-cpus", "-j", type=int,
                        default=max(1, multiprocessing.cpu_count() - 1))
    parser.add_argument("--skip-conformers", action="store_true",
                        help="Skip Tier-1 conformer generation, use DB 3DPSA only")
    return parser.parse_args()


if __name__ == "__main__":
    main(parse_args())

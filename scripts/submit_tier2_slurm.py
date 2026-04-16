"""
submit_tier2_slurm.py
---------------------
Submits 5 independent SLURM jobs for tier2_crest.py (one per reference compound),
then queues a merge job that runs automatically after all 5 complete.

Reference compounds (indexed 0–4):
  0 — Hexapeptide  (impermeable,  6-mer, HBD=6)
  1 — 1NMe3        (permeable,    6-mer, N-Me, HBD=3)
  2 — CsA          (permeable,    11-mer, chameleonic, expected ΔPSA ~75 Å²)
  3 — DP-172       (permeable,    CHUGAI 2013)
  4 — c*[PSLYF]    (impermeable,  HBD=8)

Usage:
  # Dry run — preview the SLURM scripts without submitting
  python scripts/submit_tier2_slurm.py --dry-run

  # Submit all 5 jobs + merge
  python scripts/submit_tier2_slurm.py

  # Submit only specific compounds (e.g. re-run failed jobs)
  python scripts/submit_tier2_slurm.py --compounds 2 4

  # Custom resource settings
  python scripts/submit_tier2_slurm.py --cpus 16 --mem 32G --time 12:00:00

Requirements on the cluster:
  - conda environment with crest, xtb, rdkit installed (set CONDA_ENV below)
  - This repo cloned to the same relative path on the cluster
  - feature_matrix.csv in results/ (needed for DP-172 and PSLYF SMILES)

Author: Jorge Carmona, Chameleon_Predictor project
"""

import argparse
import subprocess
import sys
import textwrap
from pathlib import Path

# ── Configuration — edit these for the cluster ────────────────────────────────
PARTITION   = "all"               # SLURM partition (run `sinfo` to list available)
CONDA_ENV   = "chameleon"         # conda environment with crest + xtb + rdkit
CONDA_SH    = "~/miniconda3/etc/profile.d/conda.sh"  # adjust if miniforge3
CPUS        = 8                   # CPUs per compound job (= --threads passed to CREST)
MEM         = "16G"               # RAM per job
TIME        = "08:00:00"          # Wall-time limit per compound (HH:MM:SS)
# CsA (11-mer) is the most complex — if 8h is not enough, increase only for job 2

# Paths — relative to repo root on the cluster
SCRIPT_PATH = "scripts/tier2_crest.py"
OUTDIR      = "results"
MATRIX_CSV  = "results/feature_matrix.csv"
LOGS_DIR    = "results/slurm_logs"

# ── Compound metadata (mirrors REFERENCE_COMPOUNDS in tier2_crest.py) ─────────
COMPOUNDS = [
    {"idx": 0, "short": "HexPep",  "name": "Hexapeptide",       "permeable": False},
    {"idx": 1, "short": "1NMe3",   "name": "N-Me Hexapeptide",  "permeable": True},
    {"idx": 2, "short": "CsA",     "name": "Cyclosporin A",     "permeable": True},
    {"idx": 3, "short": "DP172",   "name": "DP-172",            "permeable": True},
    {"idx": 4, "short": "PSLYF",   "name": "c*[PSLYF]",        "permeable": False},
]


def build_crest_script(cpd: dict, cpus: int, mem: str, time_limit: str,
                       repo_root: str) -> str:
    """
    Return a SLURM batch script string for one compound.
    repo_root: absolute path to the Chameleon_Predictor directory on the cluster.
    """
    idx   = cpd["idx"]
    short = cpd["short"]
    name  = cpd["name"]

    return textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name=crest_{short}
        #SBATCH --partition={PARTITION}
        #SBATCH --cpus-per-task={cpus}
        #SBATCH --mem={mem}
        #SBATCH --time={time_limit}
        #SBATCH --output={LOGS_DIR}/crest_{idx}_{short}_%j.out
        #SBATCH --error={LOGS_DIR}/crest_{idx}_{short}_%j.err

        # ── Environment ──────────────────────────────────────────────────────
        source "${{HOME}}/{CONDA_SH.lstrip('~/')}"
        conda activate {CONDA_ENV}

        cd {repo_root}

        echo "======================================================="
        echo "  Compound {idx}: {name}"
        echo "  Node    : $(hostname)"
        echo "  CPUs    : {cpus}"
        echo "  Started : $(date)"
        echo "======================================================="

        python {SCRIPT_PATH} \\
            --compound {idx} \\
            --threads  {cpus} \\
            --matrix   {MATRIX_CSV} \\
            --outdir   {OUTDIR}

        EXIT_CODE=$?
        echo "Finished: $(date)  |  exit=$EXIT_CODE"
        exit $EXIT_CODE
    """)


def build_merge_script(job_ids: list[int], repo_root: str) -> str:
    """
    Return a SLURM batch script that merges per-compound CSVs.
    Depends on all 5 compound jobs via afterok: only runs if all succeed.
    """
    dep = ":".join(str(j) for j in job_ids)
    return textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name=crest_merge
        #SBATCH --partition={PARTITION}
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=4G
        #SBATCH --time=00:10:00
        #SBATCH --dependency=afterok:{dep}
        #SBATCH --output={LOGS_DIR}/crest_merge_%j.out
        #SBATCH --error={LOGS_DIR}/crest_merge_%j.err

        source "${{HOME}}/{CONDA_SH.lstrip('~/')}"
        conda activate {CONDA_ENV}

        cd {repo_root}

        echo "Merging per-compound results..."
        python {SCRIPT_PATH} --merge --outdir {OUTDIR}

        echo "Done: $(date)"
    """)


def submit_script(script_text: str, dry_run: bool) -> int | None:
    """
    Write script to a temp file and submit via sbatch.
    Returns SLURM job ID, or None if dry_run.
    """
    import tempfile, os

    with tempfile.NamedTemporaryFile(mode="w", suffix=".sh",
                                     delete=False, prefix="slurm_") as f:
        f.write(script_text)
        tmp = f.name

    if dry_run:
        print(script_text)
        os.unlink(tmp)
        return None

    result = subprocess.run(
        ["sbatch", tmp],
        capture_output=True, text=True,
    )
    os.unlink(tmp)

    if result.returncode != 0:
        print(f"ERROR: sbatch failed\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    # sbatch prints: "Submitted batch job 123456"
    job_id = int(result.stdout.strip().split()[-1])
    return job_id


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Submit tier2_crest.py CREST jobs to SLURM (one job per compound)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples:
              python scripts/submit_tier2_slurm.py --dry-run
              python scripts/submit_tier2_slurm.py
              python scripts/submit_tier2_slurm.py --compounds 2 4 --cpus 16
        """),
    )
    p.add_argument(
        "--compounds", nargs="+", type=int, default=None,
        metavar="IDX",
        help="Compound indices to submit (default: 0 1 2 3 4). "
             "Use to re-run failed jobs without resubmitting all.",
    )
    p.add_argument(
        "--cpus", type=int, default=CPUS,
        help=f"CPUs per compound job — passed as --threads to CREST (default: {CPUS})",
    )
    p.add_argument(
        "--mem", default=MEM,
        help=f"Memory per job, e.g. 16G (default: {MEM})",
    )
    p.add_argument(
        "--time", default=TIME, dest="time_limit",
        help=f"Wall-time per job, HH:MM:SS (default: {TIME}). "
             "CsA (idx=2) is slowest — 8h is conservative.",
    )
    p.add_argument(
        "--repo-root", default=None,
        help="Absolute path to Chameleon_Predictor on the cluster. "
             "Defaults to the parent of this script's directory.",
    )
    p.add_argument(
        "--no-merge", action="store_true",
        help="Skip submitting the merge job (submit compounds only).",
    )
    p.add_argument(
        "--dry-run", action="store_true",
        help="Print SLURM scripts without submitting.",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # Determine repo root
    if args.repo_root:
        repo_root = args.repo_root
    else:
        # Resolve relative to this script
        repo_root = str(Path(__file__).resolve().parent.parent)

    # Create log directory (if running locally before pushing; safe to re-run)
    logs = Path(repo_root) / LOGS_DIR
    logs.mkdir(parents=True, exist_ok=True)

    # Which compounds to run
    indices = args.compounds if args.compounds is not None else [c["idx"] for c in COMPOUNDS]
    selected = [c for c in COMPOUNDS if c["idx"] in indices]

    if not selected:
        print("No matching compound indices.", file=sys.stderr)
        sys.exit(1)

    print(f"{'[DRY RUN] ' if args.dry_run else ''}Submitting {len(selected)} compound job(s):")
    for c in selected:
        label = "permeable" if c["permeable"] else "impermeable"
        print(f"  [{c['idx']}] {c['name']:30s}  ({label})")
    print()

    job_ids = []
    for cpd in selected:
        script = build_crest_script(
            cpd, args.cpus, args.mem, args.time_limit, repo_root
        )
        if args.dry_run:
            print(f"{'─'*60}")
            print(f"  SCRIPT for compound {cpd['idx']} ({cpd['short']})")
            print(f"{'─'*60}")
        job_id = submit_script(script, args.dry_run)
        if job_id is not None:
            print(f"  Submitted compound {cpd['idx']} ({cpd['short']}): job {job_id}")
            job_ids.append(job_id)

    # Merge job — only if all 5 compounds were submitted
    all_five_submitted = (not args.no_merge and
                          not args.dry_run and
                          set(indices) == {0, 1, 2, 3, 4})

    if args.dry_run and not args.no_merge and set(indices) == {0, 1, 2, 3, 4}:
        print(f"\n{'─'*60}")
        print("  SCRIPT for merge job (dependency: afterok:JOB1:JOB2:JOB3:JOB4:JOB5)")
        print(f"{'─'*60}")
        submit_script(build_merge_script([0, 0, 0, 0, 0], repo_root), dry_run=True)

    elif all_five_submitted and len(job_ids) == 5:
        merge_script = build_merge_script(job_ids, repo_root)
        merge_id = submit_script(merge_script, dry_run=False)
        if merge_id is not None:
            dep_str = ":".join(str(j) for j in job_ids)
            print(f"\n  Submitted merge job: {merge_id}  (depends on {dep_str})")

    elif not args.no_merge and not args.dry_run and set(indices) != {0, 1, 2, 3, 4}:
        print(
            "\nNote: merge job not submitted — only a subset of compounds was run.\n"
            "Run manually after all 5 complete:\n"
            f"  python {SCRIPT_PATH} --merge --outdir {OUTDIR}"
        )

    if not args.dry_run:
        print(f"\nMonitor with:  squeue -u $USER")
        print(f"Check logs in: {LOGS_DIR}/")
        print(f"After merge:   cat {OUTDIR}/tier2_crest_table.csv")


if __name__ == "__main__":
    main()

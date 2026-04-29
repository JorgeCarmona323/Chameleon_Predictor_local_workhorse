"""
submit_tier2_slurm.py
---------------------
Submit independent SLURM jobs for crest_v3.1.py (one per reference compound).
Each job runs a single compound and writes its own run folder under results/crest_runs/.

Examples:
  python scripts/submit_tier2_slurm.py --dry-run
  python scripts/submit_tier2_slurm.py
  python scripts/submit_tier2_slurm.py --compounds 2 4
  python scripts/submit_tier2_slurm.py --cpus 8 --mem 24G --time 24:00:00
"""

import argparse
import subprocess
import sys
import textwrap
from datetime import datetime
from pathlib import Path

PARTITION = "all"
CONDA_ENV = "chameleon"
CONDA_SH = "~/miniconda3/etc/profile.d/conda.sh"
CPUS = 4
MEM = "16G"
TIME = None

SCRIPT_PATH = "scripts/crest_v3.1.py"
OUTDIR = "results"
LOGS_BASE = "results/slurm_logs"

COMPOUNDS = [
    {"idx": 0, "short": "HexPep", "name": "Hexapeptide",   "permeable": False},
    {"idx": 1, "short": "CsA",    "name": "Cyclosporin A", "permeable": True},
    {"idx": 2, "short": "PSLYF",  "name": "c*[PSLYF]",     "permeable": False},
    {"idx": 3, "short": "DP955",  "name": "DP-955",        "permeable": True},
    {"idx": 4, "short": "DP944",  "name": "DP-944",        "permeable": False},
]


def build_crest_script(cpd: dict, cpus: int, mem: str, time_limit: str | None,
                       repo_root: str, logs_dir: str, dry_run: bool) -> str:
    idx = cpd["idx"]
    short = cpd["short"]
    name = cpd["name"]
    dry_flag = " --dry-run" if dry_run else ""
    time_line = f"#SBATCH --time={time_limit}\n" if time_limit else ""

    return textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name=crest_{short}
        #SBATCH --partition={PARTITION}
        #SBATCH --ntasks=1
        #SBATCH --cpus-per-task={cpus}
        #SBATCH --mem={mem}
        {time_line}#SBATCH --output={logs_dir}/crest_{idx}_{short}_%j.out
        #SBATCH --error={logs_dir}/crest_{idx}_{short}_%j.err

        source "${{HOME}}/{CONDA_SH.lstrip('~/')}"
        conda activate {CONDA_ENV}
        export OMP_NUM_THREADS={cpus}
        export OPENBLAS_NUM_THREADS=1
        export MKL_NUM_THREADS=1

        cd {repo_root}

        echo "======================================================="
        echo "  Compound {idx}: {name}"
        echo "  Node    : $(hostname)"
        echo "  CPUs    : {cpus}"
        echo "  Started : $(date)"
        echo "======================================================="

        python {SCRIPT_PATH} \\
            --compound {idx} \\
            --threads {cpus} \\
            --outdir {OUTDIR}{dry_flag}

        EXIT_CODE=$?
        echo "Finished: $(date)  |  exit=$EXIT_CODE"
        exit $EXIT_CODE
    """)


def submit_script(script_text: str, preview_only: bool) -> int | None:
    import os
    import tempfile

    with tempfile.NamedTemporaryFile(mode="w", suffix=".sh", delete=False, prefix="slurm_") as f:
        f.write(script_text)
        tmp = f.name

    if preview_only:
        print(script_text)
        os.unlink(tmp)
        return None

    result = subprocess.run(["sbatch", tmp], capture_output=True, text=True)
    os.unlink(tmp)

    if result.returncode != 0:
        print(f"ERROR: sbatch failed\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    return int(result.stdout.strip().split()[-1])


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Submit crest_v3.1.py jobs to SLURM (one job per compound)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples:
              python scripts/submit_tier2_slurm.py --dry-run
              python scripts/submit_tier2_slurm.py
              python scripts/submit_tier2_slurm.py --compounds 2 4 --cpus 8
        """),
    )
    p.add_argument(
        "--compounds", nargs="+", type=int, default=None, metavar="IDX",
        help="Compound indices to submit (default: 0 1 2 3 4).",
    )
    p.add_argument("--cpus", type=int, default=CPUS,
                   help=f"CPUs per job and --threads value passed to CREST (default: {CPUS})")
    p.add_argument("--mem", default=MEM,
                   help=f"Memory per job, e.g. 16G (default: {MEM})")
    p.add_argument("--time", default=TIME, dest="time_limit",
                   help="Wall-time per job, HH:MM:SS (default: partition default)")
    p.add_argument(
        "--repo-root", default=None,
        help="Absolute path to Chameleon_Predictor on the cluster. Defaults to the repo root.",
    )
    p.add_argument("--dry-run", action="store_true",
                   help="Print SLURM scripts without submitting them.")
    p.add_argument("--submit-dry-run-jobs", action="store_true",
                   help="Submit real SLURM jobs that run crest_v3.1.py with its --dry-run flag.")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    repo_root = args.repo_root or str(Path(__file__).resolve().parent.parent)

    run_ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    logs_dir = f"{LOGS_BASE}/run_{run_ts}"
    Path(repo_root, logs_dir).mkdir(parents=True, exist_ok=True)

    indices = args.compounds if args.compounds is not None else [c["idx"] for c in COMPOUNDS]
    selected = [c for c in COMPOUNDS if c["idx"] in indices]

    if not selected:
        print("No matching compound indices.", file=sys.stderr)
        sys.exit(1)

    preview_only = args.dry_run and not args.submit_dry_run_jobs
    job_mode = "DRY-RUN job scripts" if preview_only else (
        "dry-run jobs" if args.submit_dry_run_jobs else "real jobs"
    )

    print(f"Submitting {len(selected)} compound {job_mode}:")
    for c in selected:
        label = "permeable" if c["permeable"] else "impermeable"
        print(f"  [{c['idx']}] {c['name']:30s} ({label})")
    print(f"  Logs -> {logs_dir}/\n")

    for cpd in selected:
        script = build_crest_script(
            cpd=cpd,
            cpus=args.cpus,
            mem=args.mem,
            time_limit=args.time_limit,
            repo_root=repo_root,
            logs_dir=logs_dir,
            dry_run=args.submit_dry_run_jobs,
        )
        if preview_only:
            print("-" * 60)
            print(f"SCRIPT for compound {cpd['idx']} ({cpd['short']})")
            print("-" * 60)
        job_id = submit_script(script, preview_only=preview_only)
        if job_id is not None:
            mode_label = "dry-run" if args.submit_dry_run_jobs else "real"
            print(f"  Submitted {mode_label} job for compound {cpd['idx']} ({cpd['short']}): job {job_id}")

    if not preview_only:
        print("\nMonitor with: squeue -u $USER")
        print(f"Check logs in: {logs_dir}/")
        if args.submit_dry_run_jobs:
            print("Each job will run crest_v3.1.py with --dry-run and still create its own run folder under results/crest_runs/.")
        else:
            print("Each job will write its own run folder under results/crest_runs/.")


if __name__ == "__main__":
    main()

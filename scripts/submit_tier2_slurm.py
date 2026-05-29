"""
submit_tier2_slurm.py
---------------------
Submits 5 independent SLURM jobs for tier2_crest.py (one per reference compound).
Each job runs water then CHCl3 CREST sequentially and writes its own CSV.
After all 5 complete, run the merge step locally:
  python scripts/tier2_crest.py --merge --outdir results

Reference compounds (indexed 0–4):
  0 — Hexapeptide  (impermeable,  6-mer,  HBD=6)
  1 — CsA          (permeable,   11-mer,  chameleonic, expected ΔPSA ~75 Å²)
  2 — c*[PSLYF]   (impermeable, 11-mer,  HBD=8)
  3 — DP-955       (permeable,   15-mer,  CHUGAI 2013)
  4 — DP-944       (impermeable, 15-mer,  CHUGAI 2013)

Usage:
  # Dry run — preview the SLURM scripts without submitting
  python scripts/submit_tier2_slurm.py --dry-run

  # Submit all 5 jobs
  python scripts/submit_tier2_slurm.py

  # Submit only specific compounds (e.g. re-run failed jobs)
  python scripts/submit_tier2_slurm.py --compounds 2 4

  # Custom resource settings
  python scripts/submit_tier2_slurm.py --cpus 20 --mem 32G

Requirements on the cluster:
  - conda environment with crest, xtb, rdkit installed (set CONDA_ENV below)
  - This repo cloned to the same relative path on the cluster
  - feature_matrix.csv in results/ (needed for PSLYF SMILES)

Author: Jorge Carmona, Chameleon_Predictor project
"""

import argparse
import subprocess
import sys
import textwrap
from datetime import datetime
from pathlib import Path

# ── Configuration — edit these for the cluster ────────────────────────────────
PARTITION   = "all"               # SLURM partition (run `sinfo` to list available)
CONDA_ENV   = "chameleon"         # conda environment with crest + xtb + rdkit
CONDA_SH    = "~/miniconda3/etc/profile.d/conda.sh"
CPUS        = 20                  # CPUs per compound job (= --threads passed to CREST)
MEM         = "32G"               # RAM per job
TIME        = None                # Wall-time limit per compound — None = no limit (cluster default)

# Paths — relative to repo root on the cluster
SCRIPT_PATH = "scripts/crest_v3.2.py"
OUTDIR      = "results"
LOGS_BASE   = "results/slurm_logs"

# ── Compound metadata (mirrors REFERENCE_COMPOUNDS in crest_v3.2.py) ──────────
COMPOUNDS = [
    {"idx": 0, "short": "HexPep",  "name": "Hexapeptide",         "permeable": False},
    {"idx": 1, "short": "CsA",     "name": "Cyclosporin A",       "permeable": True},
    {"idx": 2, "short": "CsO",     "name": "Cyclosporin O",       "permeable": True},
    {"idx": 3, "short": "PSLYF",   "name": "c*[PSLYF]",          "permeable": False},
    {"idx": 4, "short": "WhC3",    "name": "White_compd3",        "permeable": True},
    {"idx": 5, "short": "DOPC_R",  "name": "DOPC_3-12-8-12_R",   "permeable": True},
    {"idx": 6, "short": "DOPC_S",  "name": "DOPC_3-12-8-12_S",   "permeable": True},
    {"idx": 7, "short": "Brain1",  "name": "Brain_6-4-4-13",      "permeable": True},
    {"idx": 8, "short": "DOPC2",   "name": "DOPC_6-5-8-12",       "permeable": True},
    {"idx": 9, "short": "CsA_v2", "name": "Cyclosporin A (rerun)", "permeable": True},
]


def build_crest_script(cpd: dict, cpus: int, mem: str, time_limit: str | None,
                       repo_root: str, logs_dir: str) -> str:
    """
    Return a SLURM batch script string for one compound.
    repo_root: absolute path to the Chameleon_Predictor directory on the cluster.
    time_limit: HH:MM:SS string or None (omit --time, use partition default).
    logs_dir: timestamped subdirectory for this submission's logs.
    """
    idx   = cpd["idx"]
    short = cpd["short"]
    name  = cpd["name"]
    time_line = f"#SBATCH --time={time_limit}\n        " if time_limit else ""

    return textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name=crest_{short}
        #SBATCH --partition={PARTITION}
        #SBATCH --cpus-per-task={cpus}
        #SBATCH --mem={mem}
        {time_line}#SBATCH --output={logs_dir}/crest_{idx}_{short}_%j.out
        #SBATCH --error={logs_dir}/crest_{idx}_{short}_%j.err

        # ── Environment ──────────────────────────────────────────────────────
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
            --threads  {cpus} \\
            --outdir   {OUTDIR}

        EXIT_CODE=$?
        echo "Finished: $(date)  |  exit=$EXIT_CODE"
        exit $EXIT_CODE
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
        help="Wall-time per job, HH:MM:SS (default: no limit — uses partition default). "
             "Pass e.g. --time 12:00:00 to impose a cap.",
    )
    p.add_argument(
        "--repo-root", default=None,
        help="Absolute path to Chameleon_Predictor on the cluster. "
             "Defaults to the parent of this script's directory.",
    )
    p.add_argument(
        "--dry-run", action="store_true",
        help="Print SLURM scripts without submitting.",
    )
    p.add_argument(
        "--test-slurm", action="store_true",
        help="Submit a minimal SLURM job that runs CREST directly (no Python) "
             "to verify CREST works in the SLURM environment.",
    )
    return p.parse_args()


def build_test_script(repo_root: str) -> str:
    xyz = f"{repo_root}/results/crest_runs/HexPep/HexPep_start.xyz"
    return (
        "#!/bin/bash\n"
        f"#SBATCH --job-name=test_crest\n"
        f"#SBATCH --partition={PARTITION}\n"
        "#SBATCH --cpus-per-task=1\n"
        "#SBATCH --mem=4G\n"
        f"#SBATCH --output={repo_root}/test_crest_%j.out\n"
        f"#SBATCH --error={repo_root}/test_crest_%j.err\n"
        f'source "${{HOME}}/{CONDA_SH.lstrip("~/")}" \n'
        f"conda activate {CONDA_ENV}\n"
        "export OMP_NUM_THREADS=1\n"
        "export OPENBLAS_NUM_THREADS=1\n"
        "export MKL_NUM_THREADS=1\n"
        "echo \"=== ENV ===\"\n"
        "which crest\n"
        "which xtb\n"
        "echo \"PWD=$PWD\"\n"
        f"cd {repo_root}\n"
        f"crest {xyz} --alpb water --T 1 --noreftopo --keepdir\n"
    )


def main() -> None:
    args = parse_args()

    if args.repo_root:
        repo_root = args.repo_root
    else:
        repo_root = str(Path(__file__).resolve().parent.parent)

    if args.test_slurm:
        script = build_test_script(repo_root)
        job_id = submit_script(script, dry_run=False)
        print(f"Submitted test job: {job_id}")
        print(f"Watch: tail -f {repo_root}/test_crest_{job_id}.out")
        return

    # Timestamped log directory — one folder per submission
    run_ts  = datetime.now().strftime("%Y%m%d_%H%M%S")
    logs_dir = f"{LOGS_BASE}/run_{run_ts}"
    Path(repo_root, logs_dir).mkdir(parents=True, exist_ok=True)

    indices = args.compounds if args.compounds is not None else [c["idx"] for c in COMPOUNDS]
    selected = [c for c in COMPOUNDS if c["idx"] in indices]

    if not selected:
        print("No matching compound indices.", file=sys.stderr)
        sys.exit(1)

    print(f"{'[DRY RUN] ' if args.dry_run else ''}Submitting {len(selected)} compound job(s):")
    for c in selected:
        label = "permeable" if c["permeable"] else "impermeable"
        print(f"  [{c['idx']}] {c['name']:30s}  ({label})")
    print(f"  Logs → {logs_dir}/\n")

    for cpd in selected:
        script = build_crest_script(
            cpd, args.cpus, args.mem, args.time_limit, repo_root, logs_dir
        )
        if args.dry_run:
            print(f"{'─'*60}")
            print(f"  SCRIPT for compound {cpd['idx']} ({cpd['short']})")
            print(f"{'─'*60}")
        job_id = submit_script(script, args.dry_run)
        if job_id is not None:
            print(f"  Submitted compound {cpd['idx']} ({cpd['short']}): job {job_id}")

    if not args.dry_run:
        print(f"\nMonitor with:  squeue -u $USER")
        print(f"Check logs in: {logs_dir}/")
        print(f"\nAfter jobs complete, scp results back to local and analyse.")


if __name__ == "__main__":
    main()

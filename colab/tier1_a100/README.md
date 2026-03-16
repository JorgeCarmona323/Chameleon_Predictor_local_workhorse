# Tier-1 Colab — A100 High RAM Setup

## What's in this folder

| File | Purpose |
|---|---|
| `tier1_a100.ipynb` | The notebook — open directly in Colab |
| `README.md` | This file |

No other local dependencies. All processing code is self-contained in the notebook.

---

## Step 1 — Set the runtime

In Colab: **Runtime → Change runtime type**
- Hardware accelerator: **A100 GPU**
- Runtime shape: **High RAM**

> The GPU itself is not used for computation. A100 + High RAM gives you **12 CPU cores** and **83 GB RAM**, which is what speeds this up.

---

## Step 2 — Upload files to Google Drive

Create this folder structure in your Drive **before running**:

```
MyDrive/
└── chem269_tier1/
    ├── pampa_curated.csv       ← upload from data/pampa_curated.csv
    └── results/                ← create this empty folder
```

**To resume from the local run's checkpoint** (recommended):
Also upload `results/conformer_descriptors_checkpoint.csv` into `MyDrive/chem269_tier1/results/`.
The notebook will detect it and skip already-completed molecules.

---

## Step 3 — Open and run the notebook

1. Go to [colab.research.google.com](https://colab.research.google.com)
2. **File → Upload notebook** → select `tier1_a100.ipynb`
3. Run cells **top to bottom**

---

## Output files (written to `MyDrive/chem269_tier1/results/`)

| File | When written | Description |
|---|---|---|
| `conformer_descriptors_checkpoint.csv` | Every 200 molecules | Incremental save — safe to interrupt |
| `run_log.csv` | Every 200 molecules | Per-molecule wall time + errors |
| `conformer_descriptors_raw.csv` | On completion | **Final output** — copy to `results/` locally |
| `run_summary.txt` | On completion | Descriptor stats + failure breakdown |

**After the run:** download `conformer_descriptors_raw.csv` and place it in your local `results/` folder so the rest of the pipeline can use it.

---

## Resuming after a disconnect

Colab Pro sessions can disconnect after ~12 hours.

1. Re-open the notebook in Colab (same session or new)
2. Re-run **all cells top to bottom** — Cell 5 reads the checkpoint and skips completed molecules
3. Cell 6 picks up where it left off

---

## Estimated runtime (A100, 12 CPUs, n_confs=50)

| Molecules | Estimated time |
|---|---|
| 7,298 (full dataset) | ~3–6 hours |
| Remaining after local checkpoint (~6,800) | ~3–5 hours |

Actual speed depends on molecule size. Large macrocycles (>70 heavy atoms) take longer per molecule.

---

## Troubleshooting

**`FileNotFoundError: pampa_curated.csv`**
Drive isn't mounted or the file path is wrong. Check that the file is at exactly `MyDrive/chem269_tier1/pampa_curated.csv`.

**`WARNING: RAM < 50 GB`**
You're not on the High RAM runtime. Go to Runtime → Change runtime type and select High RAM.

**Session disconnects mid-run**
Normal for long jobs. Resume as described above — all completed molecules are already saved.

**`embed_failed` errors in the log**
Expected for a small fraction of molecules with unusual ring systems. These are skipped and flagged in `run_log.csv`.

# Tier-2 CREST+ALPB — Google Colab Setup

**Purpose:** Run dual-dielectric CREST conformer sampling on ~400 stratified CycPeptMPDB compounds using 2 Google Colab Pro accounts (2 notebooks each, 4 total).

---

## Step 0 — Prerequisites (Sunday night / Monday AM)

1. Verify Tier-1 has finished:
   ```
   ls -lh results/conformer_descriptors_raw.csv
   ```
   If only a checkpoint exists (`results/conformer_descriptors_checkpoint.csv`), wait for Tier-1 to complete.

2. Generate the 4 batch CSVs (run locally):
   ```
   python colab/generate_batches.py
   ```
   Outputs: `colab/batches/batch_1.csv` through `batch_4.csv` (~100 compounds each)

---

## Step 1 — Google Drive folder structure

Create this exact folder in **each** Google Drive account:

```
MyDrive/
└── chem269_tier2/
    ├── colab_utils.py       ← upload from colab/colab_utils.py
    ├── batch_1.csv          ← upload to Account 1 only
    ├── batch_2.csv          ← upload to Account 1 only
    ├── batch_3.csv          ← upload to Account 2 only
    ├── batch_4.csv          ← upload to Account 2 only
    └── results/             ← create empty folder; results saved here automatically
```

### Files to upload

| File | Where to get it | Upload to |
|---|---|---|
| `colab_utils.py` | `colab/colab_utils.py` | Both accounts |
| `batch_1.csv` | `colab/batches/batch_1.csv` | Account 1 |
| `batch_2.csv` | `colab/batches/batch_2.csv` | Account 1 |
| `batch_3.csv` | `colab/batches/batch_3.csv` | Account 2 |
| `batch_4.csv` | `colab/batches/batch_4.csv` | Account 2 |

**Do not upload the notebook itself** — open it directly in Colab (see Step 2).

---

## Step 2 — Open the notebook in Colab

1. Go to [colab.research.google.com](https://colab.research.google.com)
2. File → Upload notebook → select `colab/tier2_colab.ipynb`
3. Make sure you are signed into the correct Google account before uploading

Repeat for each of the 4 notebook instances (2 per account).

---

## Step 3 — Configure each notebook (Cell 4)

Before running, edit the config cell to match the batch:

| Notebook | Account | `BATCH_CSV` | `BATCH_NAME` |
|---|---|---|---|
| Notebook 1 | Account 1 | `.../chem269_tier2/batch_1.csv` | `batch_1` |
| Notebook 2 | Account 1 | `.../chem269_tier2/batch_2.csv` | `batch_2` |
| Notebook 3 | Account 2 | `.../chem269_tier2/batch_3.csv` | `batch_3` |
| Notebook 4 | Account 2 | `.../chem269_tier2/batch_4.csv` | `batch_4` |

Also verify:
- `UTILS_PATH` = `/content/drive/MyDrive/chem269_tier2/colab_utils.py`
- `RESULTS_DIR` = `/content/drive/MyDrive/chem269_tier2/results/`
- `N_THREADS` = `4` (Colab Pro; use `8` on Pro+)

---

## Step 4 — Run the notebook

Run cells **in order, top to bottom**:

| Cell | Action | Notes |
|---|---|---|
| Cell 1 | Install condacolab | **Runtime restarts automatically** — normal |
| Cell 2 | Install CREST + xtb + RDKit | Takes ~5-8 min; only once per session |
| Cell 3 | Mount Google Drive | Sign in when prompted |
| Cell 4 | Set configuration | Edit batch name/path here |
| Cell 5 | Import colab_utils | Confirms functions loaded |
| Cell 6 | Load batch | Shows compound count and PAMPA range |
| Cell 7 | **Main processing loop** | Runs CREST on every compound |
| Cell 8 | Results summary | Run after loop completes |
| Cell 9 | Download results | Optional — results also on Drive |

**After the Cell 1 restart:** Cells 1 and 2 are already complete — skip to Cell 3.

---

## Step 5 — Run all 4 notebooks simultaneously

Start all 4 notebooks at roughly the same time:
- Account 1: open two browser tabs, one per notebook
- Account 2: same on a different browser or incognito window

Each notebook is independent. They write to separate result files:
- `results/tier2_crest_batch_1.csv`
- `results/tier2_crest_batch_2.csv`
- `results/tier2_crest_batch_3.csv`
- `results/tier2_crest_batch_4.csv`

Results are written to Drive **after every compound** — safe to interrupt and resume.

---

## Step 6 — Resume after interruption

Colab Pro sessions disconnect after ~12 hours. To resume:

1. Re-run Cells 1–6 (Cell 1 may not need reinstall if session was just paused)
2. Cell 6 will detect the existing results file and skip completed compounds
3. Re-run Cell 7 — only remaining compounds will be processed

---

## Step 7 — Download and merge results (Monday evening)

Once all 4 notebooks finish, download the results from Drive or use Cell 9 in each notebook.

Then merge locally:
```python
import pandas as pd
from pathlib import Path

dfs = [pd.read_csv(f) for f in Path("results/").glob("tier2_crest_batch_*.csv")]
merged = pd.concat(dfs, ignore_index=True)
merged.to_csv("results/tier2_crest_all.csv", index=False)
print(f"Total compounds: {len(merged)}")
print(f"Successful: {merged['error'].isna().sum()}")
```

---

## Estimated runtime

| Scenario | Time per compound | 100 compounds |
|---|---|---|
| Small peptide (~10 heavy atoms) | 8-12 min | ~16-20 h |
| Medium cyclic peptide (~50 heavy atoms) | 15-25 min | ~25-40 h |
| Large cyclic peptide (>70 heavy atoms) | 25-45 min | ~40-75 h |

With 4 notebooks running in parallel: **expect 10-20 hours total** for 400 compounds.

---

## Troubleshooting

**"crest: command not found" after Cell 2**
- Re-run Cell 2; mamba install sometimes silently fails on first attempt.

**"No module named colab_utils"**
- Check that `colab_utils.py` is at exactly `/content/drive/MyDrive/chem269_tier2/colab_utils.py`
- Re-run Cell 3 (mount Drive) then Cell 5.

**CREST exits with non-zero code / timeout**
- Compounds with >100 heavy atoms or very complex ring systems may timeout (60 min limit).
- These appear in results with `error='crest_timeout'` or `error='crest_failed'` — they are skipped.

**Runtime disconnects mid-loop**
- Resume procedure in Step 6 above. All completed compounds are already saved to Drive.

**Drive quota warning**
- Each result CSV is ~5-10 KB. 400 compounds = ~2-4 MB total. Not a concern.

# Chameleon Predictor — Reproducible Pipeline (3 notebooks)

A user-friendly, transferable walkthrough of the core pipeline as **three notebooks**:

```
   SMILES ──▶ 1. GENERATE ──▶ unique conformers ──▶ 2. CALCULATE ──▶ per-conformer + ──▶ 3. PLOT ──▶ distribution
              (CREST/xTB,       water/ mem/           (geometry stats)   summary CSVs        figures
               RMSD dedup)
```

**Geometry-first.** The CREST ensemble is reduced by RMSD to a diverse set of **unique conformers** (default 20 per solvent), and the analysis summarizes the **distribution** of geometric descriptors over that set — medians, IQRs, ranges — *not* Boltzmann-weighted means. Energies and Boltzmann weights are retained only as **metadata / QC**.

| # | Notebook | Role | Standalone? |
|---|---|---|---|
| 1 | `01_ensemble_generation.ipynb` | SMILES → deduplicated unique-conformer ensembles (`ensemble.sdf`+`ensemble.json`) | Engine-backed: imports `scripts/crest_conformers_standalone.py`, needs `xtb`+`crest` |
| 2 | `02_descriptor_calculation.ipynb` | ensembles → per-conformer table + geometric summary (median/min/max/range/IQR) | ✅ fully inlined, no `scripts/` import |
| 3 | `03_report_figures.ipynb` | per-conformer table → 4 distribution figures (boxplots, medians, points, PMI scatter) | ✅ fully inlined, no `scripts/` import |

**Standalone means:** Notebooks 2 and 3 carry their own descriptor/plotting code inlined in a cell — they do **not** import from `scripts/`, so you can run them from a bare checkout (or upload just the notebook + its input file). The inlined cells are *copies*; if the repo scripts change, re-inline. Notebook 1 is the exception: it drives the external `xtb`/`crest` binaries through the repo engine, so it needs `scripts/` + the binaries.

---

## Environments — what you need

| Stage | Needs to run | Python packages | External binaries |
|---|---|---|---|
| **1 · generate** (sampling) | Linux / macOS / WSL (or Colab) | `rdkit`, `numpy` | **`xtb`**, **`crest`** on `PATH` |
| **1 · inspect** (notebook cells) | any machine | `rdkit`, `numpy` | — |
| **2 · calculate** | any machine | `rdkit`, `numpy`, `pandas` | — |
| **3 · plot** | any machine | `pandas`, `matplotlib`, `scipy` (+ `rdkit` for Fig 1) | — |
| *all notebooks (kernel)* | any machine | `ipykernel` / Jupyter | — |

**Local kernel:** on this repository the ready kernel is the **`base`** conda env (Jupyter + rdkit + matplotlib + scipy). CREST has no native Windows build, so the *generation* step (Notebook 1) needs Linux/macOS/WSL — Google Colab works (there's a Colab setup cell in NB1).

Recreate the analysis env anywhere:
```bash
conda create -n chameleon -c conda-forge python=3.11 rdkit numpy pandas matplotlib scipy jupyter ipykernel
# to also generate ensembles (Notebook 1), on a Linux host add:
conda install -c conda-forge xtb crest
```

---

## How to run

```bash
conda activate base            # or your recreated 'chameleon' env
jupyter lab notebooks/pipeline/
```

Run in order (1 → 2 → 3): each stage reads the previous stage's **output file**. You set input/output paths in each notebook's Step 1. Notebook 1 can be skipped for a demo — Notebooks 2 and 3 run on the example ensembles already in `results/conformers/` (note: those predate deduplication, so the example figures are dense; a fresh ~20-conformer run gives clean plots).

---

## Note on the 3D-PSA definition

The 3D-PSA of record (Notebooks 2–3) follows the **Ono 2019 / Begnini 2021** convention: SASA over N/O + polar hydrogens, with the Ertl sulfur rule (reduced sulfur contributes 0). The raw `psa` field inside `ensemble.json` uses an older definition and is **not** used downstream.

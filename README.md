# Chameleon Predictor

**3D conformational descriptors + dual-solvent free energies for cyclic-peptide membrane permeability.**

Cyclic peptides can cross membranes despite high polarity by folding into compact, H-bond-shielded
conformations in lipophilic environments and re-exposing polar groups in water — *chameleonic*
behavior. 2D descriptors (TPSA, MolLogP) are conformationally blind to this. This repo samples real
conformer **ensembles per solvent** with CREST/GFN2-xTB, scores them with solvated single-points, and
derives Boltzmann-weighted 3D descriptors (ΔPSA, ΔG_transfer, IMHB, Rg, shape) from the result.

```
   SMILES ──▶ 1. GENERATE ──▶ 2. SCORE ENERGIES ──▶ 3. DESCRIBE ──▶ 4. PLOT
              CREST/GFN2-xTB    xtb single-points    Boltzmann-wtd    figures
              + ALPB, per       (CPCM-X / ALPB)      3D descriptors
              solvent           → ΔG_transfer        → ΔPSA, IMHB, Rg
```

## The pipeline

| # | Stage | Script | Reads | Writes |
|---|-------|--------|-------|--------|
| 1 | **Generate** | `scripts/crest_v3.2.py` (CLI + compound registry) → `scripts/crest_engine.py` (engine) | SMILES | `<solvent>/ensemble.xyz`, `ensemble.sdf`, `metadata.json` |
| 2 | **Score energies** | `scripts/free_energy_calculator.py` | `ensemble.xyz` per solvent | per-conformer CSV (energies, pops) + `.summary.csv` (**ΔG_transfer**) |
| 3 | **Describe** | `scripts/ensemble_descriptors.py` (+ `phys_descriptors_v3.py` lib) | `ensemble.sdf` + energies (or stage-2 CSV) | descriptor CSV |
| 4 | **Plot** | *being consolidated — see Status* | descriptor CSV | figures |

**Stages 1 and 2 need Linux** (xtb/CREST have no native Windows build). Stages 3 and 4 are pure
Python + RDKit and run anywhere, including Windows.

## Quickstart

### 1. Clone and create the environment

```bash
git clone <this-repo> && cd Chameleon_Predictor
conda env create -f environment.yml      # creates env: chameleon
conda activate chameleon
```

> **If `conda env create` fails with a `libmamba` / solver error** — e.g.
> `module 'libmambapy' has no attribute 'QueryFormat'` or *"solver backend (libmamba) was not
> recognized; choose one of: classic"* — your base conda's libmamba solver is broken (a common
> conda issue, unrelated to this repo). Switch to the classic solver and retry:
> `conda config --set solver classic` then re-run `conda env create -f environment.yml`.
> (Or repair it: `conda update -n base conda conda-libmamba-solver libmamba libmambapy`.)

> **CREST version pin:** `environment.yml` is pinned to **`crest=2.12` + `xtb=6.7.1`** — the
> confirmed-working combo (mirrors the cluster env `chameleon_crest212`). **CREST 3.x crashes
> reproducibly** on these macrocycles during iMTD-GC, so don't bump it without re-validating.

Split envs are also available if you prefer one env per role (`envs/sim.yml` = sampling,
`envs/calc.yml` = descriptors, `envs/ml.yml` = ML); see [`docs/environments.md`](docs/environments.md).

### 2. Install xtb **with CPCM-X** — required for stage 2

⚠️ **This is the one non-obvious step.** The conda-forge xtb is compiled **without** the CPCM-X
solvation library (`"CPCM-X library was not included in this version of xTB"`), so stage 2 fails on
every conformer if you rely on it. This script installs the official grimme-lab release binary,
which has CPCM-X built in, and verifies it:

```bash
bash scripts/setup_xtb.sh
```

Idempotent — safe to re-run. It downloads to `~/xtb-dist`, then proves CPCM-X works with a real
single-point before declaring success.

### 3. Verify your setup

```bash
bash scripts/tests/run_all.sh
```

Checks python, rdkit, xtb on PATH, `XTBPATH`, and **that xtb actually has CPCM-X** (ORCA is optional
and only warns). Green light here means you're ready.

### 4. Run

`scripts/env.sh` is the single source of truth for the runtime — it activates conda (if present) and
puts the CPCM-X xtb first on `PATH`. Every SLURM wrapper sources it, so generation and scoring always
use the same binary.

#### On an HPC (SLURM)

```bash
sbatch scripts/crest_hexane_hexpep_slurm.sh              # stage 1: one compound, hexane
sbatch scripts/crest_hexane_array_slurm.sh               # stage 1: 14-compound array (throttled %3)
sbatch scripts/slurm_free_energy.sh                      # stage 2: ΔG_transfer over results/conformers/
sbatch scripts/solvent_model_comparison_hexpep_slurm.sh  # stage 2: ALPB vs CPCM-X comparison
```

#### On Google Colab / any Linux box (no scheduler)

The SLURM wrappers are just thin shells around plain Python — call the stages directly:

```bash
source scripts/env.sh          # or skip conda entirely and pip install rdkit numpy pandas

# stage 1 — generate ensembles (needs xtb + crest on PATH)
python scripts/crest_v3.2.py --compound 0 --threads 4 --outdir results --solvents water,chcl3,hexane

# stage 2 — score them (needs the CPCM-X xtb from setup_xtb.sh)
python scripts/free_energy_calculator.py --method cpcmx --ewin 8 --ref water --jobs 4 \
    --leg water=results/conformers/HexPep/water/ensemble.xyz \
    --leg hexane=results/conformers/HexPep/hexane/ensemble.xyz \
    --out results/free_energy/fe_HexPep.csv

# stage 3 — descriptors, weighted by the stage-2 populations
python scripts/ensemble_descriptors.py --run-dir results/conformers/HexPep --apolar hexane \
    --energies-csv results/free_energy/fe_HexPep.csv --name HexPep -o results/descriptors.csv
```

On Colab, `scripts/setup_xtb.sh` works as-is (Linux x86_64). **CREST is separate** — install it with
`conda install -c conda-forge crest=2.12` (via `condacolab`) or drop a CREST release binary on `PATH`.
You only need CREST for stage 1; if you upload existing ensembles, stages 2–4 run without it.

There are also notebooks in [`notebooks/pipeline/`](notebooks/pipeline/) that walk the same stages
interactively.

## Key conventions

- **Solvents.** `--solvents` takes `LABEL=SOLVENT` or a bare name (`water,chloroform=chcl3,hexane`).
  LABEL is the output folder; SOLVENT is the xtb/CREST keyword. The apolar transfer phase of record
  is **hexane** — parameterized in both ALPB (generation) and CPCM-X (scoring), so no surrogate is
  needed. (ALPB has no cyclohexane.)
- **ΔG_transfer** = G_ens(apolar) − G_ens(water), the partition free energy. It benefits from
  cross-solvent error cancellation, so it's more robust than any single absolute solvation energy.
- **Phase-specific geometry.** CREST is run separately per solvent; each phase is scored on *its own*
  ensemble. Geometry is never re-optimized during scoring.
- **`--ewin 8`** pre-trims each ensemble to within 8 kcal/mol of its lowest CREST energy before the
  expensive single-points. Conformers beyond that carry negligible Boltzmann weight.
- **CPCM-X energies already include Gsolv.** xtb writes `res%e_total = dG_solv + E_gas` before the
  JSON is emitted (`main.F90:996` → `cpx.F90:108` → `json.F90:160`), so **never** pass
  `--cpcmx-add-gsolv` — it double-counts.

## Repository layout

```
scripts/
  env.sh                    # single source of truth: conda + CPCM-X xtb  (source this)
  setup_xtb.sh              # one-time: install + verify the CPCM-X xtb
  tests/                    # preflight suite — run_all.sh auto-discovers check_*.sh / test_*.py
  crest_engine.py           # stage 1 engine (RDKit embed → xtb pre-opt → CREST → SDF/XYZ)
  crest_v3.2.py             # stage 1 CLI + REFERENCE_COMPOUNDS registry
  free_energy_calculator.py # stage 2
  ensemble_descriptors.py   # stage 3
  phys_descriptors_v3.py    # shared descriptor library (PSA, IMHB, SASA)
  *_slurm.sh                # HPC wrappers (all source env.sh)
  ml/                       # separate research layer: CREMP mining, TabPFN, MACE benchmarks
  _archive/                 # superseded scripts (kept, not deleted)
notebooks/pipeline/         # interactive walkthrough of the stages
docs/                       # methodology, experiments, reports
results/
  conformers/<name>/<solvent>/   # ensembles
  free_energy/                   # stage-2 CSVs
  embeddings/                    # cached RDKit embeddings (reused across solvent legs)
```

## Status

| Stage | State |
|-------|-------|
| 1 · Generate | ✅ Working — 14 hexane ensembles + water/chcl3 legs generated |
| 2 · Score energies | ✅ Working — ALPB validated; CPCM-X unlocked via `setup_xtb.sh` |
| 3 · Describe | ⚠️ Rewired for the current ensemble format + stage-2 populations — **needs a verification run** |
| 4 · Plot | ❌ Not yet consolidated — figure code currently lives as one-offs in `scripts/_archive/` |

## History

This repo grew out of a class project whose write-up (CycPeptMPDB, ETKDG Tier-1 descriptors, UMAP
analysis, AUC results) is preserved at
[`docs/chem269_final_project.md`](docs/chem269_final_project.md). Note that its conclusion
*"Tier-2 CREST: Failed"* reflects Colab resource limits at the time and is **no longer true** — CREST
sampling is the production method here.

## References

- Pracht et al. (2020). Automated exploration of the low-energy chemical space with CREST. *PCCP*
- Grimme et al. (2025). FlexiSol: solvation/partition benchmark on flexible molecules. *Chem. Sci.*
- Stahn et al. CPCM-X — open-source COSMO-RS-family solvation. *J. Phys. Chem.*
- Rezai et al. (2006). Conformational flexibility, IMHB, and passive membrane permeability. *JACS*
- Ono et al. (2019) / Begnini et al. (2021). 3D-PSA conventions for macrocycles.

# Chameleon Pipeline

A physics-based pipeline for **passive membrane permeability of macrocycles / cyclic peptides**
(the beyond-rule-of-5, "bRo5", space). Give it a SMILES; it returns solvent-dependent 3D
conformer ensembles, a water→apolar **transfer free energy (ΔG_transfer)**, and **3D descriptors**
(3D PSA, radius of gyration, intramolecular H-bonds, ΔPSA) that capture *molecular chameleonicity*
— the ability to shield polar surface in a low-dielectric (membrane-like) environment.

```
SMILES ──► [1] CREST conformers ──► [2] CPCM-X ΔG_transfer ──► [3] 3D descriptors
           (per solvent, GFN2/GFN-FF)  (solvated single-points)   (Boltzmann-weighted)
```

One command runs all three: **`python pipeline.py --smiles "…" --name Foo`**.

---

## What's in here

| file | role |
|---|---|
| `pipeline.py` | the wrapper — chains the three stages, one command |
| `scripts/crest_engine.py` | **stage 1** — SMILES → RDKit embed → xTB pre-opt → CREST ensembles per solvent |
| `scripts/free_energy_calculator.py` | **stage 2** — CPCM-X (or ALPB) single-points → ΔG_transfer + per-conformer energies |
| `scripts/descriptors_calculator.py` | **stage 3** — 3D descriptors, Boltzmann-weighted by the stage-2 populations |
| `scripts/descriptor_equations.py` | the "equation sheet" — the PSA/Rg/IMHB/SASA math `descriptors_calculator.py` imports (not run directly) |
| `environment.yml` · `setup_xtb.sh` · `env.sh` | the runtime (see Install) |
| `check_install.py` | verifies everything is present, incl. the CPCM-X build |

---

## Install (local)

Needs **Linux** (xtb + CREST are Linux binaries). Descriptor-only use (stage 3) works anywhere with rdkit.

```bash
# 1. python deps + CREST
conda env create -f environment.yml
conda activate chameleon_pipeline

# 2. the CPCM-X-enabled xtb (conda-forge xtb lacks CPCM-X — the energy stage needs it)
bash setup_xtb.sh

# 3. put xtb + the env on PATH (source this every session, or from a job script)
source env.sh

# 4. confirm it's all ready — including that CPCM-X actually works
python check_install.py
```

`check_install.py` prints a per-item PASS/FAIL and exits nonzero if anything is missing.

---

## Quick start

```bash
# full run from a SMILES (water + chloroform + hexane)
python pipeline.py --smiles "O=C1N[C@@H](C(OC)=O)CSCC(C=CC=C2)=C2C(OC[C@H]1NC([C@@H]3CCCN3C(C)=O)=O)=O" --name Begnini1

# large / very flexible macrocycle → use GFN-FF so CREST finishes in reasonable time
python pipeline.py --smiles "<big SMILES>" --name CsA --gfn ff

# already have ensembles? skip generation, just score + describe
python pipeline.py --run-dir results/pipeline/<run>_Foo --name Foo
```

---

## CLI toggles

| flag | stage | default | what it does |
|---|---|---|---|
| `--smiles "<S>"` | input | — | the molecule (with `--name`) |
| `--name <str>` | input | `molecule` | label for the run/output |
| `--run-dir <dir>` | input | — | **skip stage 1**; score + describe an existing ensemble dir (`water/` + apolar/) |
| `--charge <int>` | — | auto | formal charge (auto from SMILES, else 0) |
| `--gfn {2,1,0,ff}` | **1** | `2` | CREST level of theory. `2`=GFN2-xTB (accurate). `ff`=**GFN-FF** (≈100–1000× faster; use for large flexible macrocycles). *Stage-2 scoring stays GFN2/CPCM-X regardless.* |
| `--max-confs <int>` | **1** | keep all | cap conformers kept per solvent |
| `--method {cpcmx,alpb}` | **2** | `cpcmx` | solvation model for the energy single-points. CPCM-X is more accurate for partition ratios; ALPB is a faster fallback |
| `--ewin <kcal>` | **2** | `8` | pre-trim conformers to within this window of the lowest energy before scoring (throughput; conformers beyond carry negligible Boltzmann weight) |
| `--threads <int>` | 1&2 | all cores | CPU threads for CREST/xTB |
| `--outdir <dir>` | 1 | `results/pipeline` | where stage 1 creates the run directory |
| `--force` | all | off | re-run stages even if their outputs already exist (default = resume/skip) |

Run `python pipeline.py --help` for the same list.

**Solvents:** the reference is always **water**; the apolar legs are **chloroform** (the
literature-standard membrane mimic) and **hexane**. Descriptors and ΔG_transfer are produced for
every apolar leg present → both `water→chloroform` and `water→hexane`.

---

## Outputs

Everything lands in the run directory (`--outdir/<timestamp>_<name>/`):

```
water/ chloroform/ hexane/ensemble.{xyz,sdf}   # stage 1 — the conformer ensembles
free_energy_cpcmx.csv                            # stage 2 — per-conformer energies + weights
free_energy_cpcmx.summary.csv                    #          ΔG_transfer + coordinate provenance
descriptors_chloroform.csv                       # stage 3 — 3D descriptors (water vs chloroform)
descriptors_hexane.csv                           #          3D descriptors (water vs hexane)
pipeline_manifest.json                           # what ran, and the exact input paths
```

- **`*.summary.csv`** — the headline: `dGtransfer_water->chloroform_kcal`, `dGtransfer_water->hexane_kcal`,
  plus `source_dir`/`src_<solvent>` recording exactly which ensembles were scored.
- **`descriptors_*.csv`** — Boltzmann-weighted PSA3D, Rg, IMHB, SASA breakdown, and the cross-solvent
  `delta_*` terms (the chameleonicity signal).

---

## Notes & limits

- **xtb/CREST are Linux-only.** Stage 3 (descriptors) needs only rdkit/numpy/pandas and runs anywhere,
  so you can generate ensembles on Linux and analyze anywhere.
- **CREST is the slow step.** A small 6-mer on GFN2 is quick; large/floppy macrocycles can take many
  hours — reach for `--gfn ff` and/or `--max-confs`.
- **CPCM-X is the install gotcha.** If the energy stage errors with "CPCM-X library was not included",
  your `xtb` is the conda build — run `setup_xtb.sh` and `source env.sh` so the grimme-lab release is
  first on PATH. `check_install.py` catches this up front.
- `--compound <N>` (run a built-in reference compound by index) needs the optional `scripts/crest_v3.2.py`
  registry, which is not in this core install — use `--smiles` instead.

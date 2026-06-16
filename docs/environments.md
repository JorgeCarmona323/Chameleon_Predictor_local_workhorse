# Project environments — one env per role

We split conda envs by **role**, not by script, so dependency stacks never collide
(rdkit's numpy/boost pins vs. pytorch vs. xtb). Specs live in `envs/`; create each with
`conda env create -f envs/<role>.yml`.

```
                 ┌─────────────────┐   ensembles    ┌──────────────────┐  feature CSV  ┌────────────────┐
   compounds ──► │ chameleon-sim   │ ─────────────► │ chameleon-calc   │ ────────────► │ chameleon-ml   │
                 │ (CREST / xtb)   │  water/ mem/   │ (rdkit descrip.) │ results/*.csv │ (sklearn/xgb)  │
                 │  HPC cluster    │                │  local + cluster │               │  models, SHAP  │
                 └─────────────────┘                └──────────────────┘               └────────────────┘
```

The seam that makes this clean: **calc produces feature CSVs, ml consumes them.** The ML
env therefore needs *no* rdkit, and the two never fight over numpy pins.

## Roles

| Env | Spec | Built where | Purpose |
|---|---|---|---|
| **chameleon-sim** | `envs/sim.yml` | HPC cluster (Linux) — mirrors existing `chameleon_crest212` | GFN2-xTB / CREST 2.12 conformer search; also the regime-2 `xtb --ohess`/CENSO Hessian job |
| **chameleon-calc** | `envs/calc.yml` | local + cluster — mirrors existing `rdkit_env` | 2D/3D descriptors, surface areas, shape, figures → **writes feature CSVs** |
| **chameleon-ml** | `envs/ml.yml` | local | benchmark, feature importance, UMAP → **reads feature CSVs** |
| *MACE* (separate) | `MACE` env (pip/GPU) | local GPU | MACE-OFF / OpenMM — a *distinct* simulation stack, intentionally not merged |

## Which env runs which script

**chameleon-sim** (CREST / xtb — cluster):
`crest_conformers.py`, `crest_v3.*.py`, `submit_tier2_slurm*.py`

**chameleon-calc** (rdkit descriptors / figures — local + cluster):
`phys_descriptors_v2.py`, `phys_descriptors_v3.py`, `ensemble_descriptors.py`,
`compute_2d_descriptors.py`, `build_feature_matrix.py`, `compare_conformer_ensembles.py`,
`compare_methods_csa.py`, `verify_diazirine_isomers.py`, `cremp_deltapsa.py`,
`curate_data.py`, `conformer_engine.py`, `analyze_dataset_coverage.py`,
`reference_compounds.py`, `_validate_compounds.py`,
`make_isomer_figures.py`, `plot_isomer_comparison.py`, `plot_csa_validation.py`,
`pipeline_overview_slides.py`, `export_docx.py`, `export_slides.py`, `validate_csa_water.py`

**chameleon-ml** (modeling — local):
`feature_benchmark.py`, `correlation_analysis.py`, `umap_visualization.py`,
`plot_tier2_results.py`

**MACE** (separate GPU sim): `benchmark_mace_vs_xtb.py`, `mace_tip3p_feasibility.py`

## Convention

Each active script declares its env on a comment at the very top of the file (line 1,
above the docstring — a comment isn't a statement, so the module docstring is preserved):

```python
# env: chameleon-calc
"""... module docstring ..."""
```

This is greppable (`grep -rl "env: chameleon-ml" scripts/`) and self-documenting — no
runner framework. Rule of thumb if a new script is unclear: **needs rdkit → calc;
runs crest/xtb → sim; trains/evaluates a model → ml.**

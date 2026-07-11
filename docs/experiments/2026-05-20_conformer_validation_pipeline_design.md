# Experiment Design — Conformer Validation & Descriptor Pipeline
**Date:** 2026-05-20

---

## Motivation

CREST/xTB runs take 10+ days for 11-mer cyclopeptides on 20 CPU threads (no GPU path). Before committing to the full compound library we need to know: (1) how well CREST ensembles capture the real conformational landscape vs cheaper methods, and (2) which 3D descriptors actually predict permeability. This session designs the validation experiment and feature selection workflow.

---

## Experiment 1 — RDKit vs CREST Conformer Quality

**Goal:** Establish whether fast RDKit conformer generation is sufficient for descriptor extraction, or whether CREST-level sampling is required.

### Steps

1. **RDKit fast conformers** — ETKDGv3 + MMFF94 energy minimization (seconds per molecule)
   - Generate ~500–5000 conformers per molecule
   - Compute MMFF94 energies → Boltzmann weights at 300 K
   - Prune low-population conformers (< 1% cumulative Boltzmann weight threshold)

2. **RMSD ensemble check** — build pairwise RMSD matrix over the pruned ensemble
   - Identifies spread and redundancy
   - Flag if top-2 conformers RMSD < 0.5 Å (likely over-pruned)

3. **Cluster major conformations** — hierarchical clustering (Ward linkage) on RMSD matrix
   - Target: 3–6 macrostates per molecule
   - Label each cluster by dominant cis-amide pattern and H-bond motif post-hoc

4. **Compare RDKit clusters vs CREST clusters vs experimental geometry**
   - For CsA: compare against A1 crystal (water, CCDC2149649) and DEKSAN (CHCl3, CCDC1138505)
   - For HexPep (when data transferred): compare against Rezai 2006 NMR structure if available
   - Metric: RMSD of cluster centroid to experimental; Boltzmann population of the experimental-matching cluster

5. **Compare descriptors across methods** — extract the same descriptor set from RDKit ensemble, CREST ensemble, and experimental CIF geometry; report deviation table

### Molecules
- CsA (water ensemble ready: 23 conformers; CHCl3 pending job 259118)
- HexPep (user reports both solvents complete on server — need `scp` transfer)
- CsO (water + CHCl3 queued, ~10-day ETA)

---

## Experiment 2 — NMR-Derived Descriptor Comparison

**Goal:** Use published solution NMR data as a ground-truth conformational reference independent of crystal packing.

### NMR observables to extract (where published data exists)

| Observable | Source | Descriptor proxy |
|---|---|---|
| NOE inter-proton distances | Bhatt JACS 2022 (CsA A1) | Boltzmann-averaged H–H distance from ensemble |
| Scalar couplings (³J_HNα) | Karplus equation → φ dihedral | Boltzmann-averaged φ per residue |
| Chemical shift deviation | Published δ vs random coil | Proxy for H-bond / shielding environment |

**For CsA:** Bhatt JACS 2022 SI has NOE data for the A1 (water) conformer. Compare Boltzmann-averaged CREST H–H distances vs experimental NOEs.

**For HexPep:** Rezai & Lokey JACS 2006 may have NMR J-coupling or NOE data — needs literature check once paper is accessible.

### Comparison metric
Mean absolute error (MAE) between predicted (Boltzmann-averaged) and experimental distance/dihedral. Lower MAE = better ensemble quality.

---

## Experiment 3 — Descriptor Correlation & Feature Importance

**Goal:** Identify which 3D descriptors carry independent predictive signal for permeability.

### Steps

1. **Descriptor extraction** — run full descriptor set on all available CREST ensembles:
   - Boltzmann PSA (3D SASA-based, not 2D TPSA)
   - Boltzmann nHBD_exposed (SASA per HBD)
   - ΔPSA(water − CHCl3)
   - ΔnHBD(water − CHCl3)
   - Cis-amide propensity per bond (all backbone amide bonds)
   - Rg, anisotropy (inertia tensor)
   - ΔΔG bias index (blocked until CHCl3 ensembles complete)
   - Congruent conformer population (blocked until cross-ensemble RMSD clustering implemented)

2. **Correlation matrix** — Pearson + Spearman between all descriptor pairs
   - Identify collinear groups (|r| > 0.9 → keep one representative)
   - Visualize as heatmap

3. **Tree model feature importance** — XGBoost or Random Forest on descriptor matrix vs binary permeability label (PAMPA threshold)
   - **Dataset constraint: more molecules than descriptors** — target ≥ 2× ratio; with ~50 compounds from CycPeptMPDB + our 7 reference compounds, limit descriptor count to ≤ 20 for initial model
   - SHAP values for per-feature contribution
   - Leave-one-out CV or 5-fold CV for generalization estimate

4. **Report** — which descriptors survive collinearity pruning AND appear in top-5 SHAP importance

---

## Implementation Gaps to Close First

| Gap | Needed for | Priority |
|---|---|---|
| FreeSASA / rdFreeSASA wired to ensemble | 3D PSA_exposed, HBD SASA | High |
| Joint RMSD clustering (MDAnalysis) | Macrostate definition, congruent pop | High |
| `compare_crest_vs_experimental.py` (gemmi → CIF load) | Exp. 1 & 2 | High |
| MACE-OFF CUDA fix (upgrade PyTorch to cu128) | Speed up future CREST replacement | Medium |
| CCS prediction tool | Shape descriptor | Low |

---

## Cluster Job Status

| Compound | Solvent | Status | Notes |
|---|---|---|---|
| CsA | water | Done | 23 conformers, `data/CREST_CsA_20260512/` |
| CsA | CHCl3 | Running | Job 259118, day 10+ |
| HexPep | water | Done (server) | Need `scp` to local |
| HexPep | CHCl3 | Done (server) | Need `scp` to local |
| CsO | water | Queued | `sbatch scripts/cso_crest_slurm.sh` (--compound 2 runs both solvents) |
| CsO | CHCl3 | Queued | Same job as above |
| PSLYF | water + CHCl3 | **Not submitted** | `--compound 3` |
| DP-955 | water + CHCl3 | **Not submitted** | `--compound 4` |
| DP-944 | water + CHCl3 | **Not submitted** | `--compound 5` |
| WhC3 | water + CHCl3 | **Not submitted** | `--compound 6` |

### Commands to queue remaining jobs (run on server after `git pull`):
```bash
sbatch scripts/cso_crest_slurm.sh   # already done — CsO compound 2
# Create equivalent scripts for compounds 3–6 or pass --compound directly:
python scripts/crest_v3.2.py --compound 3 --threads 20 --outdir results  # PSLYF
python scripts/crest_v3.2.py --compound 4 --threads 20 --outdir results  # DP-955
python scripts/crest_v3.2.py --compound 5 --threads 20 --outdir results  # DP-944
python scripts/crest_v3.2.py --compound 6 --threads 20 --outdir results  # WhC3
```
Each needs its own SLURM script (no `--time` limit, 20 CPUs, 32 GB).

---

## Open Questions

- Does Rezai 2006 JACS have solution NMR NOE or J-coupling data for the hexapeptide? (needs literature check)
- CsO NMR data (CHCl3) is image-based PDF — user will translate manually; water NMR not yet located
- Is 50 compounds from CycPeptMPDB sufficient for the tree model, or do we need to expand to the full database?

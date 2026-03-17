# Project Roadmap — Cyclic Peptide Permeability via 3D Conformational Descriptors

**Started:** 2026-03-17
**Status:** Active — CHEM 269 final project submitted, continuing as independent research

---

## Current State (2026-03-17)

| What | Status |
|---|---|
| 7,298-compound PAMPA dataset (CycPeptMPDB) | Curated |
| Tier-1 conformers (ETKDGv3+MMFF94s, n=20) | 1,502 / 7,298 — full run in progress on Colab A100 |
| Feature matrix, correlation analysis, UMAP | Complete on 1,502-compound subset |
| Best predictor: psa3d_std AUC=0.749, δPSA AUC=0.744 | +11 AUC pts over best 2D (MolLogP=0.631) |
| CREST Tier-2 | Failed — compute/environment limitations |
| xtb+GBSA single-structure | Completed but insufficient (δPSA ≈ 0 for all) |
| UMAP cluster stability | Unstable (ARI 0.07–0.38) — continuous permeability + selection bias |

**Key validated insight:** Ensemble-derived δPSA captures chameleonic behavior; single-structure methods cannot. This is the paper's core finding.

---

## Phase 0 — Immediate (by 2026-03-17, tonight)

- [ ] Submit CHEM 269 final report
- [ ] Download full Tier-1 CSV (7,298 rows) from Colab when run completes
- [ ] Re-run pipeline scripts with full dataset
- [ ] Archive current results

---

## Phase 1 — Short Term (weeks 1–4 post-submission)

### 1.1 Finish the full Tier-1 pipeline
- Run `build_feature_matrix.py`, `correlation_analysis.py`, `umap_visualization.py` on 7,298-compound Tier-1
- Expected: AUC values stabilize, UMAP clustering may improve slightly
- Check whether psa3d_std AUC holds above 0.74 at full scale

### 1.2 Expand the dataset
- **Primary target:** Search literature for non-CycPeptMPDB cyclic peptide permeability data
  - Lokey lab publications (UCSF) — Caco-2 and PAMPA series
  - Bhardwaj et al. 2022 (Science) — de novo designed cyclic peptides with experimental permeability
  - ChEMBL: search `cyclic peptide` + `membrane permeability` or `Papp`
  - PubChem BioAssay: PAMPA assays for cyclic scaffolds
- **Why this matters:** CycPeptMPDB is 66.4% permeable (selection bias). An unbiased dataset would sharpen AUC and reveal true decision boundaries.
- **Goal:** Compile 500–2,000 new compounds with experimental Caco-2/PAMPA, balance permeable/non-permeable ~50/50

### 1.3 N-methylation systematic study
- The database has N-methylation labels — run subgroup analysis
- Hypothesis: N-Me degree correlates with δPSA independently of sequence
- This is a quick analysis on the existing feature matrix — low effort, potentially publishable sub-finding

### 1.4 Sequence-based features
- Encode amino acid sequence as one-hot or property vectors (hydrophobicity, charge, size)
- Combine with 3D δ features in a simple random forest or XGBoost model
- Question: do sequence features + δPSA outperform δPSA alone?

---

## Phase 2 — Medium Term (months 1–3)

### 2.1 Better conformer sampling — CREST on HPC
- **What failed:** CREST 2.12 on Colab — memory/timeout constraints, no persistent compute
- **Fix:** Use an HPC cluster (university allocation) or AWS/GCP spot instances
  - CREST with `--alpb water` and `--alpb chcl3` for dual-dielectric sampling
  - Target: 50–100 reference compounds first, then scale
- **Expected payoff:** CREST δPSA values for CsA literature-validated at ~75 Å² — would allow direct comparison to Tier-1 ETKDGv3 values and quantify how much force-field-level sampling matters
- **Resources to look into:**
  - XSEDE/ACCESS allocations (free for academic users)
  - University HPC — check if CHEM dept has allocation
  - Grimme group CREST documentation for memory optimization flags

### 2.2 OpenMM explicit solvent MD
- Replace CREST with short (10–50 ns) OpenMM MD in explicit water + implicit membrane (GBSA or Membrane Builder)
- More physically rigorous than GBSA implicit solvent
- Use the existing RDKit conformer as starting geometry
- **Packages:** OpenMM 8.x, OpenFF force field (SMIRNOFF), ParmEd for setup
- Start with the 5 reference compounds, validate δPSA against Tier-1 and CREST
- OpenMM runs well on Colab A100 for small cyclic peptides (<30 heavy atoms)

### 2.3 Graph neural network / MPNN
- Move beyond hand-crafted features entirely
- Use a pretrained molecular GNN (e.g., ChemProp, DimeNet++, or Uni-Mol) fine-tuned on permeability
- Input: SMILES → implicit 3D geometry via GNN encoder
- Compare to Tier-1 δPSA: does end-to-end learning recover the same signal?
- **Interesting question:** Does a GNN trained on permeability implicitly learn δPSA?

### 2.4 Interpretability analysis
- SHAP values on a trained XGBoost model using the full feature matrix
- Which residues / substructures drive high δPSA?
- Map SHAP contributions back onto molecular structures (scaffold decomposition)
- This is the "rules governing permeability" angle — data-driven SAR

---

## Phase 3 — Long Term (months 3–12)

### 3.1 Permeability-guided de novo design
- Use δPSA as an objective function in a generative model
- Reinforce cyclic peptide sequences that maximize predicted δPSA while maintaining drug-likeness
- **Approach:** REINVENT (Astra Zeneca, open source) or RDKit + genetic algorithm
- Produce a set of novel in silico designed cyclic peptides predicted to be highly permeable

### 3.2 Explicit lipid bilayer simulations
- CHARMM-GUI Membrane Builder + OpenMM or GROMACS
- Simulate permeation free energy via umbrella sampling or metadynamics
- Gold standard but expensive (~days/compound on HPC)
- Target: 3–5 high-δPSA compounds from your dataset for proof of concept

### 3.3 Wet lab validation (if collaborators available)
- Synthesize 2–3 designed peptides from Phase 3.1
- Run PAMPA assay in-house or through collaborator (Lokey lab publishes protocols)
- Close the loop: design → predict → synthesize → measure

### 3.4 Manuscript
- Target journal: *Journal of Chemical Information and Modeling* or *J. Med. Chem.*
- Core story: ensemble δPSA > single-structure > 2D; CREST vs ETKDGv3 comparison (Phase 2.1); dataset bias analysis
- Co-authors: advisor + anyone who contributes HPC/wet lab

---

## Data Sources to Monitor

| Source | What to look for | URL/DOI |
|---|---|---|
| CycPeptMPDB updates | New PAMPA entries | cycpeptmpdb.com |
| ChEMBL | Cyclic peptide Papp assays | chembl.ebi.ac.uk |
| Lokey lab (UCSF) | Caco-2 series publications | pubmed search: Lokey cyclic peptide permeability |
| Bhardwaj et al. | De novo designed peptides | DOI: 10.1126/science.abo1940 |
| BindingDB | Cyclic peptide bioactivity | bindingdb.org |

---

## Tools to Get Working

| Tool | Purpose | Blocker | Fix |
|---|---|---|---|
| CREST 2.12 | High-quality conformer ensemble | Colab memory/timeout | HPC allocation |
| OpenMM | Explicit solvent MD | Setup complexity | OpenFF + SMIRNOFF tutorial |
| ChemProp | MPNN for permeability | None — pip install | Fine-tune on CycPeptMPDB |
| REINVENT | Generative design | Needs trained prior | Use pretrained ChEMBL prior |
| GROMACS | Lipid bilayer PMF | HPC required | After CREST works |

---

## Open Questions

1. Does the δPSA–permeability relationship hold for **non-PAMPA** assays (Caco-2, RRCK)? PAMPA has no active transport — Caco-2 does.
2. Is there a **δPSA threshold** above which permeability is near-guaranteed? The current data suggests ~40–50 Å² as a soft cutoff but needs more data.
3. How much does **ring size** modulate the δPSA–permeability relationship? CycPeptMPDB spans 4–20 residue rings.
4. Can **N-methylation pattern** (not just count) be encoded as a structural feature that improves model performance?
5. Do GNN models trained on PAMPA **transfer** to RRCK or Caco-2 endpoints without retraining?

---

## Log

| Date | Milestone |
|---|---|
| 2026-03-17 | CHEM 269 final submission. Tier-1 1,502 compounds. δPSA AUC=0.744 > MolLogP AUC=0.631. CREST failed. xtb single-structure failed. UMAP unstable. Full 7,298-compound Tier-1 run in progress. |
| — | *(update this as you go)* |

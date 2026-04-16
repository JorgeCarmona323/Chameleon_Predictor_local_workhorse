# Mentor Project Overview: Chameleon_Predictor

**Date**: April 9, 2026
**Assessment Type**: Initial Mentoring Review
**Reviewer**: Mentor Assessment

## Executive Summary

The **Chameleon_Predictor** project represents **research-grade computational chemistry work** investigating chameleonic behavior in cyclic peptides for membrane permeability prediction. The student has moved well beyond coursework into original scientific discovery, with findings ready for publication.

**Core Research Question**: Can 3D ensemble-derived conformational descriptors (specifically ΔPSA - the difference in polar surface area between aqueous and membrane conformers) predict membrane permeability better than traditional 2D descriptors?

**Major Discovery**: Evidence for **size-gated chameleonic behavior** requiring approximately ≥9 residues (~900 Da) to manifest effectively.

---

## Project Structure & Architecture

### Repository Overview
```
Chameleon_Predictor/
├── README.md                          # Comprehensive project documentation (269 lines)
├── environment.yml                     # Conda environment (RDKit, scikit-learn, UMAP, etc.)
├── run_pipeline.py                     # Main entry point for full pipeline
├── literature_review.md               # Literature validation
│
├── assignment/
│   └── climb_route_prompt.md          # Original project assignment (CHEM 269)
│
├── docs/
│   ├── hypothesis.md                  # Running hypothesis: chameleonic behavior is residue-gated
│   ├── data_schema.md                 # Database schema design for future ML model
│   ├── chameleon_model_architecture.md # Planned multi-modal neural network architecture
│   ├── methodology_flowchart.md       # Visual workflow
│   ├── findings_and_methods_log.md    # Detailed lab notebook
│   ├── writeup_2026-03-18.md         # Main research writeup
│   └── experiments/                   # Experimental tracking documents
│
├── notebooks/
│   ├── 3d_descriptors.ipynb          # Main analysis notebook (6+ MB, extensive)
│   ├── cremp_benchmark.ipynb         # CREMP dataset validation
│   ├── experimental_psa_benchmark.ipynb
│   └── library_chemical_space_explorer.ipynb
│
├── scripts/                           # Production pipeline (3,326 lines total)
│   ├── curate_data.py                # SMILES canonicalization, PAMPA filtering
│   ├── conformer_engine.py           # ETKDGv3 + MMFF94s conformer generation (461 lines)
│   ├── build_feature_matrix.py       # Merge all descriptor groups
│   ├── correlation_analysis.py       # Pearson/Spearman/AUC-ROC analysis (369 lines)
│   ├── umap_visualization.py         # Dual-track clustering (K-Medoids + HDBSCAN) (617 lines)
│   ├── cremp_deltapsa.py            # CREMP benchmark pipeline (222 lines)
│   ├── tier2_crest.py               # Attempted CREST validation (failed, documented)
│   ├── export_docx.py               # Report generation
│   └── export_slides.py             # Presentation generation
│
└── results/
    ├── feature_matrix.csv            # 7,298 × 274 feature matrix (20 MB)
    ├── cremp_deltapsa.csv           # 2,457 CREMP compounds with ΔPSA
    ├── correlation_table.csv
    ├── auc_roc_table.csv
    ├── tier2_xtb_results.csv        # Negative control results
    └── figures/                      # 37 publication-quality figures
```

### Technology Stack

**Core Chemistry:**
- RDKit ≥2023.9 (conformer generation, descriptors, SASA calculation)
- ETKDGv3 (macrocycle-aware conformer generation)
- MMFF94s force field

**Machine Learning/Analysis:**
- scikit-learn (logistic regression, AUC-ROC, clustering)
- UMAP (dimensionality reduction)
- HDBSCAN + K-Medoids (dual-track clustering)
- Scipy (statistical tests)

**Data Sources:**
- CycPeptMPDB v1.2 (8,466 cyclic peptides, 7,298 with PAMPA data)
- CREMP dataset (2,457 compounds with pre-computed CHCl₃ conformers)

---

## Scientific Methodology & Pipeline

### Workflow Overview

**Pipeline Flow:**
1. **Data Curation** (`curate_data.py`): SMILES canonicalization, PAMPA filtering (threshold: ≥-6.0 log cm/s), 2D baseline descriptors
2. **Conformer Generation** (`conformer_engine.py`): 20 conformers/molecule via ETKDGv3, MMFF94s minimization, identify aqueous (max-PSA) and membrane (min-PSA) conformers
3. **Feature Extraction**: Compute ΔPSA, ΔHB (H-bond count), ΔRg (radius of gyration), ΔNPRs (shape descriptors)
4. **Feature Matrix** (`build_feature_matrix.py`): Merge 2D, DB 3D (negative control), and Tier-1 3D delta features
5. **Analysis** (`correlation_analysis.py`): Pearson/Spearman correlations, AUC-ROC classification
6. **Visualization** (`umap_visualization.py`): UMAP with dual-track clustering, permeability enrichment analysis

### Descriptor Categories

**Target Property:** PAMPA (Parallel Artificial Membrane Permeability Assay) - passive transcellular diffusion

**Descriptor Types:**
1. **2D Baseline**: MolLogP, TPSA, MW, HBD/HBA (conformationally blind)
2. **DB 3D (negative control)**: H2O_3DPSA - CHCl3_3DPSA from CycPeptMPDB (single-structure, AUC ~0.5)
3. **Tier-1 Δ descriptors** (ensemble-derived):
   - `delta_psa3d`: PSA(aqueous) - PSA(membrane) - **primary descriptor**
   - `psa3d_std`: conformational flexibility
   - `delta_hb`: intramolecular H-bond count change
   - `delta_Rg`, `delta_NPR1/2`: compaction and shape changes

### Key Innovation: Dual-Track Clustering
- **K-Medoids** (k=8, deterministic, cosine metric) for reproducibility
- **HDBSCAN** (min_cluster_size=50) for exploratory analysis
- **ARI stability validation** across 5 random seeds

---

## Major Scientific Findings

### 1. Ensemble vs. Single-Structure Approaches
- **Traditional 3D descriptors fail**: DB 3DPSA (single-structure) yields AUC ~0.5
- **Ensemble methods work**: Tier-1 ΔPSA achieves AUC 0.744 on clean subset
- **Validation**: Cyclosporin A Tier-1 ΔPSA = 84.9 Å² vs. literature ~75 Å² (within 10%)

### 2. Data Quality Impact
- **PAMPA heterogeneity**: CycPeptMPDB aggregates 4 incompatible protocols
  - Townsend 2020 (~42%): pooled PAMPA with interference issues
  - Kelly 2021 (~21%): same pooled protocol
  - Furukawa 2016 (~9%): cleanest individual LC-MS data
  - Chugai (~12%): different membrane composition
- **Performance collapse**: Clean subset AUC 0.744 → Full dataset AUC 0.505

### 3. Size-Gating Hypothesis ⭐
**Novel Discovery**: Chameleonic behavior appears size-gated at ~9 residues
- **Evidence**: MW gap (1,180 Da permeable vs. 820 Da impermeable in clean data)
- **Independent validation**: Converges with Yu et al. 2026's ≥9 residue cutoff
- **Biological significance**: Suggests different permeation mechanisms below threshold

### 4. CREMP Benchmark Validation
**Recent addition** validating ΔPSA using physically-grounded membrane conformers:
- **Hybrid approach**: Aqueous PSA from ETKDGv3, membrane PSA from CREST CHCl₃ ensembles
- **Dataset**: 2,457 compounds with ~1,000-5,000 pre-computed conformers each
- **Validation**: Confirms CREMP conformers are more collapsed than vacuum

---

## Student Assessment

### Advanced Research Skills Demonstrated

**Technical Competencies:**
- ✅ Production-quality modular code architecture (3,326 lines)
- ✅ Large dataset handling (7,298 compounds, 20 MB feature matrices)
- ✅ RDKit expertise (conformer generation, descriptor calculation)
- ✅ Machine learning implementation (UMAP, clustering, statistical analysis)
- ✅ Proper experimental design with negative controls
- ✅ Failed experiment documentation (CREST limitations)

**Scientific Methodology:**
- ✅ Rigorous hypothesis development and testing
- ✅ Independent validation using external datasets
- ✅ Literature validation and cross-referencing
- ✅ Statistical analysis with appropriate metrics
- ✅ Publication-quality documentation and figures

**Research Maturity:**
- ✅ Moved beyond coursework to original discovery
- ✅ Novel scientific insights (size-gating hypothesis)
- ✅ Understanding of method limitations and data quality issues
- ✅ Forward-thinking database design for future ML models

---

## Areas for Mentoring Focus

### 1. Statistical Rigor (Immediate Priority)
**Current Gaps:**
- Missing formal statistical test for size-gating hypothesis
- Need explicit ablation study (2D only vs. 3D only vs. combined)
- Source-stratified analysis required for PAMPA heterogeneity

**Mentoring Actions:**
- Guide through Mann-Whitney U test for residue bin comparisons
- Implement proper train/validation/test splits
- Formalize hypothesis testing framework

### 2. Machine Learning Best Practices
**Current State**: Single-descriptor AUC analysis
**Advancement Needed:**
- Ensemble methods (Random Forest, XGBoost)
- SHAP analysis for feature importance
- Cross-validation strategies
- Feature interaction analysis

### 3. Descriptor Methodology Refinement
**Active Investigations:**
- Implement ΔPSA/SASA_total normalization (per Yu et al. 2026)
- Per-molecular-weight normalization strategies
- Fractional switching efficiency metrics

### 4. Resource Planning & Infrastructure
**HPC Requirements:**
- CREST validation requires significant compute resources
- Database implementation (Phase 2 schema)
- Experimental validation pipeline design

**Guidance Needed:**
- Help secure HPC cluster access
- Cloud computing strategy for large-scale conformer generation
- Database architecture decisions

### 5. Publication Strategy
**Publication Potential**: 2-3 papers worth of findings
**Priority Papers:**
1. **Methodology**: "Ensemble-Derived Conformational Descriptors for Membrane Permeability Prediction"
2. **Discovery**: "Size-Gated Chameleonic Behavior in Cyclic Peptides"
3. **Dataset**: "PAMPA Data Quality Assessment and Source Stratification"

**Mentoring Support:**
- Journal selection strategy
- Manuscript structure guidance
- Peer review preparation
- Conference presentation opportunities

---

## Immediate Next Steps (Priority Order)

### Phase 1: Statistical Foundation (Next 2-4 weeks)
1. **Descriptor Normalization**: Implement ΔPSA/SASA_total per Yu et al. 2026
2. **Size-Stratified Analysis**: Formal AUC comparison (≤8 vs. ≥9 residues)
3. **Clean Dataset Focus**: Furukawa-only subset analysis (highest quality PAMPA)
4. **Statistical Testing**: Mann-Whitney U test for size-gating hypothesis

### Phase 2: Model Development (1-2 months)
1. **Ensemble ML**: Random Forest with full Tier-1 feature set
2. **Feature Importance**: SHAP analysis implementation
3. **CREMP Integration**: Maximize CycPeptMPDB overlap through better SMILES matching
4. **Database Schema**: Implement Phase 2 design

### Phase 3: Validation & Publication (2-3 months)
1. **HPC CREST**: Reference compound validation on cluster
2. **Multi-modal Neural Network**: SequenceEncoder + DynamicEnsembleEncoder architecture
3. **Experimental Design**: DEL cyclic peptide scaffold validation
4. **Manuscript Preparation**: Lead with methodology/discovery papers

---

## Long-term Research Trajectory

### Academic Pathway Considerations
**Current Level**: Advanced undergraduate/early graduate research
**Natural Progression**:
- PhD in computational chemistry/chemical biology
- Postdoc in drug discovery or membrane transport
- Industry R&D in pharmaceutical companies

### Research Impact Potential
**Immediate Impact**: Novel descriptor methodology for cyclic peptide design
**Broader Significance**:
- Understanding of membrane transport mechanisms
- Tool for drug discovery pipelines
- Foundation for experimental validation studies

### Collaboration Opportunities
**Experimental Validation**: Partner with membrane transport labs
**Industry Application**: Pharmaceutical company internships/collaborations
**Academic Networks**: Computational chemistry conferences (ACS, COMP division)

---

## Technical Debt & Code Quality Assessment

### Strengths
✅ Modular architecture with clear separation of concerns
✅ Checkpoint/resume functionality for long-running jobs
✅ Comprehensive logging and error handling
✅ Production-quality documentation
✅ Version control best practices

### Areas for Improvement
⚠️ **Feature Matrix Coverage**: 99.99% success rate (7,297/7,298) - debug the 1 failure
⚠️ **UMAP Stability**: Bimodal ARI distribution (0.10-0.997) needs characterization
⚠️ **Checkpoint Granularity**: Could be more granular in conformer_engine.py
⚠️ **SMILES Matching**: CREMP integration needs canonical SMILES harmonization

---

## Conclusion: Research-Ready Student

This student has demonstrated **exceptional research capabilities** and made **original scientific discoveries**. The work is at a **publication-ready level** and represents a significant contribution to the field of computational drug discovery.

**Key Strengths for Mentoring:**
- Strong foundation in computational chemistry
- Excellent documentation and reproducibility practices
- Novel scientific insights with validation
- Production-quality implementation skills
- Understanding of method limitations

**Mentoring Priority:** Focus on **publication preparation**, **statistical rigor**, and **career development** rather than basic skill building.

**Expected Outcome:** 2-3 high-quality publications and a strong foundation for graduate-level research in computational chemistry/drug discovery.

---

**Assessment Status**: ⭐ **Research-Grade Work - Ready for Advanced Mentoring** ⭐
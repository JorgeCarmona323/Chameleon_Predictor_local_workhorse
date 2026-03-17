# Route CX-001: The Chameleon Traverse

**RouteID:** CX-001
**Wall:** Cheminformatics & Molecular Modeling (W14)
**Grade:** 5.12a
**Routesetter:** Jorge Carmona
**Time:** ~8–12 hours (compute-heavy; Colab A100 recommended)
**You'll need:** Python 3.10+, RDKit, conda, Google Colab (A100), CycPeptMPDB dataset

---

## Why this route exists

Some molecules break the rules. Cyclosporin A weighs 1,203 Da, has 11 H-bond donors and acceptors, and yet it passively crosses cell membranes like a small drug. It does this by shape-shifting — burying its polar groups in a lipid bilayer and exposing them again in water. This is called *chameleonic behavior*, and standard 2D descriptors like TPSA and LogP are blind to it.

CycPeptMPDB gives us experimental permeability data for >8,000 cyclic peptides. The database even provides 3D PSA values for a subset. But do those values actually encode the chameleonic switch? And can we do better with our own conformer ensembles?

In this route, you'll build a pipeline that generates 3D conformer ensembles for cyclic peptides, extracts delta descriptors that quantify the polar surface area switch between aqueous and membrane environments, and tests whether those descriptors predict experimental PAMPA permeability better than any 2D feature — and better than the database's own 3D values.

---

## 🦶 First hold

### Exercise 0: Environment and Data Setup

**Goal:** Get the data and build the environment before touching any chemistry.

**Step 1: Download the dataset**
Download the full CycPeptMPDB CSV from [cycpeptmpdb.com](https://cycpeptmpdb.com). Place it in your project root. It should be named something like `CycPeptMPDB_Peptide_All.csv`.

**Step 2: Set up the environment**

```bash
conda env create -f environment.yml
conda activate chem269_cycpep
```

Key packages: `rdkit`, `numpy`, `pandas`, `scikit-learn`, `umap-learn`, `matplotlib`, `seaborn`

**Step 3: Explore the dataset**

Load it in a notebook and answer:
- How many compounds are in the full database?
- How many have PAMPA measurements?
- What is `H2O_3DPSA`? What is `CHCl3_3DPSA`? What would `delta_3DPSA = H2O_3DPSA - CHCl3_3DPSA` represent physically?
- What fraction of PAMPA compounds are "permeable" (LogPexp ≥ −6.0)?

**Success check:**
- CSV loads, PAMPA subset identified (~7,298 compounds)
- You can articulate what delta 3DPSA is supposed to measure and why it might predict permeability

---

## 🤚 Undercling up

### Exercise 1: Curate the PAMPA Subset

**Goal:** Clean the dataset and establish your 2D baseline.

**Step 1: Canonicalize SMILES**

Use RDKit to parse every SMILES, extract the largest fragment (removes counterions), and re-canonicalize. Drop any compounds that fail to parse.

**Step 2: Filter to PAMPA subset**

Keep only compounds with a valid PAMPA LogPexp value. Apply the binary label:
```python
df['permeable'] = (df['PAMPA_LogPexp'] >= -6.0).astype(int)
```

**Step 3: Compute 2D baseline descriptors**

For each molecule compute: MolWt, MolLogP, TPSA, HBA, HBD, RotatableBonds, FractionCSP3, RingCount

**Step 4: Evaluate the 2D baseline**

Compute AUC-ROC for each 2D descriptor against the permeable label. Record the best one.

*Hint: What does MolLogP encode? Why might it predict PAMPA permeability for cyclic peptides?*

**Success check:**
- Curated PAMPA CSV with 2D descriptors
- AUC-ROC table for all 2D descriptors
- Best 2D AUC recorded (you'll beat it in Exercise 3)

---

## 🧱 Stem up

### Exercise 2: Expose the Database 3DPSA Baseline

**Goal:** Test whether CycPeptMPDB's existing 3D PSA values are useful.

**Step 1: Compute delta_3DPSA_db**

```python
df['delta_3DPSA_db'] = df['H2O_3DPSA'] - df['CHCl3_3DPSA']
```

**Step 2: Evaluate**

Compute AUC-ROC and Spearman ρ for `delta_3DPSA_db` vs. the permeable label.

**Step 3: Interpret**

What AUC did you get? Is it better than random (0.5)?

Now look at the actual values for a few well-known permeable cyclic peptides (e.g., Cyclosporin A, if it's in the database). What is their `delta_3DPSA_db`?

*Expected finding: delta_3DPSA_db ≈ AUC 0.50 — essentially random. The database values come from single optimized structures and cannot capture the chameleonic conformational switch.*

**Key question to answer:** Why does a single-structure approach fail here, even though the structures are 3D? What physical reality is it missing?

**Success check:**
- AUC-ROC for delta_3DPSA_db computed
- You can articulate why single-structure 3D PSA is insufficient for chameleonic molecules

---

## 🦵 Heel hook up

### Exercise 3: Build the Tier-1 Conformer Engine

**Goal:** Generate conformer ensembles and extract Δ descriptors that actually capture chameleonic behavior.

**Step 1: ETKDGv3 conformer generation**

For each molecule:
1. Generate 20 conformers with ETKDGv3 (enable macrocycle torsion library: `ETKDGv3.useSmallRingTorsions=True, useBasicKnowledge=True`)
2. Minimize all conformers with MMFF94s force field
3. Compute 3D polar SASA for each conformer using Bondi radii (`rdFreeSASA`)

**Step 2: Extract delta descriptors**

- `delta_psa3d` = PSA(max-PSA conformer) − PSA(min-PSA conformer)
- `psa3d_std` = standard deviation of PSA across all conformers
- `delta_hb` = H-bond count difference between min-PSA and max-PSA conformers
- `delta_Rg` = radius of gyration difference

**Step 3: Reference compound validation**

Before running at scale, validate on Cyclosporin A. Literature reports a ΔPSA of ~75 Å². What do you get? This is your sanity check before investing compute.

*If your CsA ΔPSA is near 0, something is wrong with your conformer generation or PSA calculation. Debug here before scaling.*

**Step 4: Scale up (Colab A100)**

The full dataset takes 13–15 hours on an A100. Use a checkpoint/resume system — save progress every N molecules so you can resume if the session drops.

**Success check:**
- CsA ΔPSA ≈ 70–90 Å² (literature ~75 Å²)
- All 5 reference compounds computed
- Checkpoint system working

---

## 🤸 Dyno up

### Exercise 4: Feature Matrix and Correlation Analysis

**Goal:** Merge all feature groups and evaluate predictive power.

**Step 1: Build feature matrix**

Merge:
- 2D baseline descriptors (Exercise 1)
- DB delta_3DPSA (Exercise 2)
- Tier-1 Δ descriptors (Exercise 3)

**Step 2: Correlation analysis**

For every feature compute:
- Pearson r and Spearman ρ vs. PAMPA LogPexp
- AUC-ROC vs. permeable label

**Step 3: Feature group comparison**

Create a bar chart comparing AUC-ROC across feature groups. Which group wins? By how much?

**Step 4: Interpret**

- Does Tier-1 ΔPSA beat the best 2D descriptor?
- Does DB delta_3DPSA beat the best 2D descriptor?
- What does the comparison between Tier-1 and DB 3D tell you about conformer sampling methodology?

**Success check:**
- Feature matrix CSV with all descriptor groups
- AUC-ROC bar chart showing group-level comparison
- You can articulate the hierarchy: Tier-1 > 2D >> DB 3D

---

## 🦵 Gaston to next move

### Exercise 5: Chemical Space Visualization

**Goal:** See where permeable and non-permeable compounds live in conformational descriptor space.

**Step 1: UMAP on 3 feature sets**
Run UMAP separately on:
- Panel A: 2D descriptors only
- Panel B: Tier-1 Δ descriptors only
- Panel C: All features combined

**Step 2: Dual-track clustering**

For each panel, apply two clustering approaches:
- K-Medoids (k=8, deterministic, cosine metric) — primary
- HDBSCAN (min_cluster_size=50) — exploratory

**Step 3: Stability validation**

Run UMAP 5 times with different random seeds. Compute ARI between clusterings. Report the stability.

**Step 4: Overlay reference compounds**

Mark CsA and the other reference compounds on all 3 panels. Where do they land?

**Key question:** Does Panel B (3D Δ features) show better permeable/non-permeable separation than Panel A (2D)? If clusters are unstable across seeds, what does that tell you about the nature of the permeability–structure relationship?

**Success check:**
- 3 UMAP panels (PNG figures)
- Stability ARI reported per panel
- Reference compounds visible on plots

---

## 🤸 Mantle to crux

### Exercise 6 (Stretch): Tier-2 Validation

**Goal:** Attempt higher-level conformer sampling to validate Tier-1.

This is the crux. Tier-1 ETKDGv3+MMFF94s is fast but uses a molecular mechanics force field that doesn't explicitly model dielectric environments. Higher-level approaches include:

**Option A: CREST + ALPB solvation**
```bash
crest molecule.xyz --T 4 --alpb water --mquick
crest molecule.xyz --T 4 --alpb chcl3 --mquick
```
CREST samples the conformational landscape at the GFN2-xTB semiempirical level with implicit ALPB solvation — far more physically rigorous than MMFF94s.

**Option B: xtb + GBSA (single-structure)**
```bash
xtb molecule.xyz --opt --gbsa water --gfn 2
xtb molecule.xyz --opt --gbsa chcl3 --gfn 2
```
Single-structure optimization at the GFN2 level. Faster than CREST but misses ensemble sampling.

**What to do:**
1. Attempt CREST on your 5 reference compounds. Document what happens.
2. If CREST fails, run xtb+GBSA as a fallback. Report the ΔPSA values.
3. Compare Tier-2 results against Tier-1. Are they consistent?

**Expected learning regardless of outcome:** Even if CREST fails, the comparison of xtb single-structure ΔPSA vs. Tier-1 ΔPSA demonstrates whether single-structure solvation is sufficient. The xtb result is its own finding.

*Note: CREST requires significant compute (4+ hours for large cyclic peptides) and may fail in Colab due to memory constraints. This is expected — document it.*

**Success check:**
- Attempt documented with method and outcome
- xtb+GBSA ΔPSA values for reference compounds reported
- Comparison to Tier-1 values discussed

---

## 🧗 Send it!

### Submission

**Submit:**
1. This climb route prompt (`assignment/climb_route_prompt.md`)
2. Your complete route solution:
   - `notebooks/3d_descriptors.ipynb` — main analysis notebook
   - `results/` — feature matrix, correlation tables, AUC tables
   - `results/figures/` — all plots
   - `docs/findings_and_methods_log.md` — findings and methods log
3. GitHub repository link

**Your notebook must include:**
- [ ] Dataset exploration (compound counts, permeability distribution, DB 3DPSA inspection)
- [ ] 2D baseline AUC-ROC table
- [ ] DB delta_3DPSA evaluation with interpretation
- [ ] Conformer engine validation (CsA ΔPSA vs literature)
- [ ] Full AUC-ROC comparison: Tier-1 vs DB 3D vs 2D
- [ ] UMAP panels (A, B, C) with clustering and reference overlays
- [ ] UMAP stability ARI results with interpretation
- [ ] Reflection: Do 3D ensemble descriptors outperform 2D? Does single-structure 3D PSA work? Why or why not?
- [ ] Tier-2 attempt documentation (outcome + interpretation)

**Reflection questions (answer in markdown cells):**

1. CycPeptMPDB provides 3D PSA values. Why do they fail to predict permeability while your Tier-1 ΔPSA succeeds? What physical reality is being missed?

2. Your UMAP clusters were unstable across random seeds. Does this mean permeability is unpredictable? Or does it tell you something else about the structure of this chemical space?

3. Tier-1 ΔPSA achieves AUC ~0.75. What are the remaining 25% of variance that it misses? Name at least two factors.

4. If you had unlimited compute, what would you do next to improve either the conformer sampling or the predictive model?

---

## 🎉 Summit!

You made it to the top.

You built a pipeline that:
- Exposed a methodological flaw in the database's own 3D features
- Validated chameleonic behavior computationally against experimental literature
- Demonstrated that ensemble sampling is a requirement, not a refinement
- Characterized the limits of force-field-level dual-dielectric modeling

This is real research. The findings are defensible, the negative results are documented, and the pipeline is reproducible. That's the summit.

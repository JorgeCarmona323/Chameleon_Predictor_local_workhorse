# Scientific Framing & Project Logic

**Project:** CHEM 269 Final — 3D Conformational Descriptors and Dual-Dielectric Solvent Modeling
**Author:** Jorge Carmona | March 2026

---

## The Core Scientific Question

> Can 3D conformational descriptors computed from a conformer ensemble reveal the **physical determinants** of membrane permeability in cyclic peptides — and do they outperform 2D descriptors in doing so?

This is a **feature discovery question**, not a single correlation claim. The goal is not to prove that ΔPSA predicts PAMPA. The goal is to determine which 3D features — if any — best separate permeable from non-permeable cyclic peptides, and what that reveals about the physics of membrane crossing.

---

## Why 2D Descriptors Fail

Standard 2D descriptors (MW, cLogP, TPSA, HBD count) miss the chameleonic mechanism entirely:

- They describe the molecule's **average or maximum** polar surface, not what it actually exposes in a membrane environment
- They cannot capture intramolecular H-bond formation — the NMR-proven mechanism by which cyclic peptides shield polar groups in low-dielectric environments (Rezai 2006, White & Lokey 2011)
- The existing 2D UMAP of the CycPeptMPDB dataset fails to separate permeable from non-permeable compounds — the direct motivation for this project

---

## The Chameleonic Mechanism (What We Are Measuring)

Cyclic peptides that cross membranes do so by **conformational switching**:

1. In **aqueous solution** (ε=78): the molecule adopts an extended conformation, exposing backbone NH donors to water for H-bonding with solvent
2. In the **membrane** (ε≈4): the molecule collapses into a compact conformation, burying NH donors in intramolecular H-bonds — shielding polar surface from the hydrophobic lipid interior

This switching is "chameleonic" — the molecule presents a different face to each environment. It is **experimentally proven** by NMR for:
- Hexapeptide → N-Me Hexapeptide pair (Rezai 2006, White & Lokey 2011): NMR NOESY in CDCl₃ directly shows intramolecular H-bond formation upon N-methylation
- Cyclosporin A (Kessler 1990, Wenger 1994, Witek 2016): NMR + MD showing ~75 Å² ΔPSA between open (aqueous) and closed (membrane) conformations

---

## Why ΔPSA Is A Proxy, Not The Mechanism

ΔPSA (polar surface area in aqueous conformer minus polar surface area in membrane conformer) is the structural **consequence** of chameleonic switching, not the cause. The actual mechanism is intramolecular H-bond formation (ΔHB).

**Important distinction in how ΔPSA is computed:**

| Method | What it actually measures |
|---|---|
| Tier-1 ΔPSA (max_PSA − min_PSA, vacuum ensemble) | Conformational flexibility — range of PSA accessible in vacuum. A proxy, not the true chameleonic switch. |
| Tier-2 ΔPSA (PSA_aq_min − PSA_mem_min, solvated ensembles) | True chameleonic switch — difference between environment-specific minimum energy conformers. This is what the proposal intends. |

Tier-1 approximates Tier-2 using a heuristic (max-PSA conformer ≈ aqueous form, min-PSA conformer ≈ membrane form). This is reasonable but is explicitly an approximation. Tier-2 CREST with ALPB solvation computes the correct quantity.

---

## The Full Feature Set Being Evaluated

We compute multiple Δ features to avoid reducing a complex physical phenomenon to a single number:

| Feature | Physical meaning | Mechanistic basis |
|---|---|---|
| **PSA_mem** | 3D polar surface in membrane conformer | Direct measure of polar exposure in the membrane — Witek 2016 showed this is the best single predictor |
| **ΔPSA** | PSA_aq − PSA_mem | Magnitude of chameleonic switch |
| **ΔHB** | H-bonds_aq − H-bonds_mem | Intramolecular H-bond formation — the actual mechanism proven by NMR |
| **ΔRg** | Radius of gyration change | Global compaction in membrane environment |
| **ΔNPR1, ΔNPR2** | Shape change (rod/sphere/disk axes) | 3D shape switching |
| **PSA_spread** | Std dev of PSA across ensemble | Conformational flexibility — ability to switch |
| **2D baseline** | MW, cLogP, TPSA, HBD, RotBonds | Reference — what 3D features must outperform |

---

## Expected Findings and Their Interpretation

| Finding | Physical interpretation |
|---|---|
| ΔHB correlates better than HBD count alone | It's not just HOW MANY donors you have — it's whether they get shielded via intramolecular H-bonds in the membrane |
| PSA_mem outperforms 2D TPSA | The membrane-form conformation matters more than the static average structure |
| ΔPSA adds predictive value over 2D alone | Chameleonic switching is a real, measurable determinant of permeability |
| ΔRg/ΔNPR support ΔPSA | Permeability requires global molecular compaction, not just local H-bond burial |
| 2D cLogP and HBD have limited discriminating power | Lipophilicity and donor count alone are insufficient for bRo5 macrolides |

---

## The Narrative Arc

```
Observation:
  2D UMAP of CycPeptMPDB fails to separate permeable from non-permeable
  cyclic peptides — 2D descriptors are insufficient
        ↓
Hypothesis:
  Permeability is a 3D conformational problem driven by chameleonic
  switching between aqueous and membrane environments
        ↓
Approach:
  Compute 3D Δ features across two dielectric environments
  (ε=78 aqueous, ε=4 membrane-mimetic) for ~7,298 compounds
        ↓
Validation (Tier-2):
  On 5 reference compounds where chameleonic behavior is NMR/MD-proven,
  CREST+ALPB (and OMEGA+GB/SA) Δ features recapitulate the known
  conformational switching — validating the computational method
        ↓
Scale (Tier-1):
  Apply validated 3D feature pipeline to full PAMPA subset
  Correlate features vs. LogPexp; AUC-ROC; UMAP colored by LogPexp
        ↓
Answer:
  Which 3D features are the strongest determinants of permeability?
  Do they outperform 2D descriptors?
  Is the chameleonic trend visible at scale across 7,298 compounds?
```

---

## What Success Looks Like

**Minimum viable result (backup plan from proposal):**
- Tier-1 Δ features computed for PAMPA subset
- Correlation table showing 3D features outperform 2D (even modestly)
- UMAP with reference compound overlay
- Honest limitations section

**Full result:**
- Tier-2 CREST validates Tier-1 directional agreement on reference compounds
- AUC-ROC: 3D combined features > 2D baseline
- UMAP: 3D features separate permeable/impermeable clusters better than 2D
- PSA_mem and ΔHB emerge as strongest individual predictors

**If correlations are weak:**
That is still a publishable result. It means: "3D Δ features capture chameleonic potential but permeability has additional determinants (efflux, flip-flop kinetics, assay heterogeneity) that a purely conformational model cannot explain." The Tier-2 mechanistic validation on reference compounds still stands independently.

---

## Key References

| Paper | Relevance |
|---|---|
| Rezai & Lokey, *JACS* **2006** | NMR proof of chameleonic mechanism; Hexapeptide reference |
| White & Lokey, *Nat Chem Biol* **2011** | NMR proof of N-methylation enabling chameleonic switch; 1NMe3 reference |
| Witek et al., *JCTC* **2016** | MD-quantified ΔPSA for CsA (~75 Å²); PSA_mem as best predictor |
| Jiang et al., *J Chem Inf Model* **2023** | CycPeptMPDB; PAMPA threshold −6.0 log cm/s |
| Pracht et al., *PCCP* **2020** | CREST conformer sampling; macrocycle benchmarks |
| Ehlert et al., *JCTC* **2021** | ALPB solvation model; ε-tunable implicit solvent |
| Naylor et al., *J Med Chem* **2018** | Systematic cyclic peptide permeability; multi-feature analysis |
| Hawkins et al., *J Chem Inf Model* **2017** | OMEGA macrocycle conformer benchmarks |

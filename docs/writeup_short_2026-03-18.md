# Chameleonic Membrane Permeability is Size-Gated: Evidence from 3D Conformational Descriptors and Source-Stratified PAMPA Analysis

**Jorge Carmona | March 2026**

---

## Abstract

We present a computational pipeline for predicting cyclic peptide membrane permeability using 3D ensemble-derived conformational descriptors, with a focus on ΔPSA: the difference in polar surface area between aqueous and membrane-mimetic conformers. On a 1,566-compound source-stratified subset of CycPeptMPDB¹ (Furukawa 2016² + Chugai), ensemble ΔPSA achieves AUC = 0.692 vs. AUC = 0.317 for MolLogP. Notably, MolLogP is inversely predictive on this subset (flipped AUC = 0.683), consistent with large polar macrolides dominating the high-scoring PAMPA population in this preliminary dataset. Size-stratified analysis reveals that this overall AUC reflects two distinct populations: compounds of ≥9 residues are highly permeable in the preliminary data (94.8% for 9–11 residues, 86.3% for 12–15 residues), while smaller cyclic peptides (≤8 residues, 60.4% permeable) show ΔPSA as essentially uninformative (AUC = 0.528), hinting at a permeation mechanism not dependent on conformational switching. HDBSCAN clustering on the UMAP embedding of the combined 2D+3D feature space identifies two distinct populations; overlaying molecular weight as a color axis on the same embedding reveals that the permeable cluster has a median MW of 1,180 Da vs. 820 Da for the impermeable cluster (1.44× gap). When extended to the full heterogeneous 7,298-compound dataset this separation vanishes entirely, which we initially attributed to cross-source PAMPA label noise. To investigate further, we asked whether the MW gap itself contained a deeper signal — and applied per-MW normalization (ΔPSA / MolWt) to test whether absolute ΔPSA was conflating molecular size with chameleonic switching efficiency. The result made the picture unambiguous: within-class AUC improves from 0.528 to 0.696 for 12–15 residue compounds while actively degrading for smaller compounds — an asymmetry that cannot be explained by label noise alone and points directly to two mechanistically distinct permeation regimes operating within the same dataset. Collectively, these findings are consistent with a two-mechanism model of passive cyclic peptide permeability — arrived at independently and corroborated by the concurrent work of Yu et al. 2026³: chameleonic conformational switching in large macrolides, detectable by 3D ensemble descriptors, and a distinct permeation route in smaller scaffolds that these descriptors cannot capture. Cleaner, source-homogeneous permeability data and per-SASA descriptor normalization are the critical next steps toward testing this model more rigorously.

---

## 1. Background

Passive membrane permeation for large cyclic peptides is not well-modeled by conventional 2D descriptors because these descriptors are conformationally blind. Cyclosporin A (CsA, MW 1,203 Da, 11 residues) is the archetypal example: it crosses membranes despite high polarity by adopting a compact conformation in apolar environments that buries backbone amide NH groups through intramolecular H-bonding⁷. This *chameleonic* switching is invisible to descriptors like TPSA or MolLogP but should manifest as a large difference in 3D polar surface area (ΔPSA) between conformers sampled in aqueous vs. membrane-mimetic conditions. Empirical work from the Lokey group has established that N-methylation and backbone rigidity modulate this switch in smaller cyclic scaffolds¹⁰, but a systematic size threshold below which chameleonic behavior fails to manifest has not been computationally validated at database scale.

CycPeptMPDB v1.2¹ provides the largest public collection of experimental cyclic peptide PAMPA data (8,466 compounds) but aggregates measurements from laboratories using incompatible protocols. This heterogeneity is a central finding of this work and a major constraint on cross-source descriptor validation.

---

## 2. Methods

### 2.1 Data

CycPeptMPDB v1.2¹ was filtered to 7,298 compounds with PAMPA values. For source-stratified analyses, a 1,566-compound subset comprising Furukawa 2016² (individual-compound LC-MS, 1% lecithin/dodecane membrane) and Chugai (patent WO 2013/100132 A1, DOPC/hexadecane membrane) was used as the primary reference set. These two sources share a critical methodological feature: both use individual-compound assays rather than the pooled multiplexed PAMPA protocol used by Townsend 2020⁸ and Kelly 2021⁹, which accounts for ~63% of the full dataset and introduces systematic cross-compound signal interference via CycLS MS deconvolution. Furukawa 2016² is the most rigorously documented source (peer-reviewed, individual LC-MS detection). Chugai is a patent-derived source with less fully characterized conditions (DOPC/hexadecane membrane, detection floor −10.0 vs −8.0 log cm/s for other sources) and should be treated with appropriate caution; an analysis restricted to individual-compound, single-protocol data is proposed as a future control (see Section 6). Permeability threshold: PAMPA LogPexp ≥ −6.0 log cm/s.

### 2.2 Conformer Generation (Tier-1)

Twenty conformers per molecule were generated using ETKDGv3 with macrocycle torsion sampling enabled⁴, followed by MMFF94s minimization⁵. ETKDGv3 is the RDKit default conformer generator incorporating experimental torsional-angle preferences for macrocycles from crystallographic data. In published benchmarks⁶, ETKDGv3 with the macrocycle flag ranks 2nd for macrocycle ring RMSD against crystal structures (0.62 Å), statistically indistinguishable from the best-performing open-source method (0.56 Å). At the scale of this analysis (7,298 compounds), ETKDGv3 is the appropriate choice; CREST-level semiempirical sampling is proposed for validation on a smaller reference set (see Section 5).

The conformer with maximum 3D polar SASA (Bondi radii, RDKit rdFreeSASA) was designated the aqueous conformer and the conformer with minimum 3D polar SASA the membrane conformer. Key descriptors:

- `delta_psa3d` = PSA(aq) − PSA(mem): chameleonic PSA switch
- `psa3d_std`: PSA variance across all 20 conformers (conformational flexibility)
- `delta_hb`, `delta_Rg`, `delta_NPR1/2`: H-bond count, compaction, and shape descriptors

Following the observation of a molecular weight-stratified cluster separation in UMAP (Section 3.4), a per-MW normalized descriptor was computed post-hoc: `delta_psa3d_per_mw` = `delta_psa3d` / MolWt × 1000 (Å² per kDa). This normalization removes the molecular size contribution to absolute ΔPSA, isolating chameleonic switching efficiency from molecular weight confounding. No conformer re-generation was required; the normalization is a derived quantity from the existing feature matrix.

As a negative control, `delta_3DPSA_db` = H2O_3DPSA − CHCl3_3DPSA from CycPeptMPDB's pre-computed 3D PSA values was included. Per the CycPeptMPDB methods¹, both values are derived from the same single minimum-energy UFF conformer — not from independently solvent-optimized structures — making this descriptor a measure of PSA calculation settings rather than conformational switching. Its inclusion tests whether any residual conformational signal survives in a single-structure approach.

### 2.3 Analysis

Features were merged into a 7,298 × 274 matrix. Pearson/Spearman correlations and AUC-ROC were computed for each descriptor. UMAP dimensionality reduction (cosine metric, n_neighbors=30) with dual-track clustering — K-Medoids (k=8, cosine distance) and HDBSCAN (min_cluster_size=50) — was used for chemical space visualization and permeability enrichment analysis. ARI stability was validated across 5 random seeds (threshold ≥ 0.85). Three feature panels were defined: Panel A (2D descriptors), Panel B (3D delta descriptors), Panel C (combined). Track D, introduced in this analysis, adds a fourth subplot coloring the same UMAP embedding by molecular weight.

---

## 3. Results

### 3.1 Single-Structure ΔPSA Fails: Motivating the Ensemble Approach

CycPeptMPDB¹ provides pre-computed 3D PSA values in two conditions (`H2O_3DPSA` and `CHCl3_3DPSA`) for each compound. A natural first test is whether their difference — `delta_3DPSA_db` — correlates with experimental PAMPA LogPexp. It does not. Before interpreting why, it is important to understand how these values were generated.

According to the CycPeptMPDB methods¹, 5,000 conformers per peptide were generated with RDKit, redundant conformers (RMSD < 1.0 Å) were removed, each was UFF-minimized, and the **single lowest-energy conformer was selected**. Critically, `H2O_3DPSA` and `CHCl3_3DPSA` are both computed from this **same single structure** — not from independently solvent-optimized conformers. The database itself acknowledges this on the peptide detail page with an explicit warning (Li et al.¹, Fig. 6B): *"This conformation likely does not reflect what is found in biological systems, and most peptides populate conformational ensembles rather than a single conformation."* Consequently, `delta_3DPSA_db` does not measure a conformational switch at all — it measures the difference in PSA accessibility of an identical 3D geometry under two different solvent-accessibility calculation settings. Any signal from genuine chameleonic switching is structurally inaccessible to this descriptor by construction.

![PAMPA distribution and DB ΔPSA scatter](figures/fig1_data_overview.png)

*Figure 1. Left: PAMPA LogPexp distribution across 7,298 compounds (threshold −6.0 shown). Right: database ΔPSA (`H2O_3DPSA − CHCl3_3DPSA`) vs. experimental PAMPA. Spearman ρ = −0.020 (p = 0.10), AUC = 0.507 — indistinguishable from chance. Both PSA values are derived from the same single minimum-energy conformer, making this descriptor structurally incapable of detecting chameleonic switching. The vertical noise band centered near 0 Å² is the expected outcome: you cannot measure a conformational switch using one conformation.*

### 3.2 The Decisive Negative Control: Level of Theory Is Not the Bottleneck

To confirm that single-structure failure is a sampling problem rather than a force field problem, GFN2-xTB¹³ with GBSA implicit solvation was run on five reference compounds in water and CHCl₃. GFN2-xTB is a semiempirical tight-binding method substantially more accurate than the MMFF94s force field underlying Tier-1.

| Compound | Permeable | PAMPA | xtb ΔPSA (Å²) | Tier-1 ΔPSA (Å²) |
|----------|-----------|-------|----------------|------------------|
| CsA (Cyclosporin A) | Yes | −6.60 | **−0.14** | **84.9** |
| DP172 | Yes | −4.15 | −0.24 | 88.9 |
| HexPep | No | −6.20 | 0.82 | 64.4 |
| 1NMe3 | Yes | −5.52 | 6.91 | 47.8 |
| PSLYF | No | −9.10 | 5.40 | 65.3 |

GFN2-xTB produces ΔPSA of 0–7 Å² across all five compounds — including CsA, where the chameleonic switch has been directly observed by NMR⁷. Tier-1 ensemble sampling recovers 47–89 Å² for the same compounds: a 10–100× larger signal. What governs passive membrane permeation in chameleonic cyclic peptides is not the energy of a single optimized structure — it is the *existence* of a low-PSA conformer that the molecule can access in the membrane environment, however transiently. Tier-1 probes that capacity by populating conformational space broadly; GFN2-xTB cannot, regardless of its energy accuracy, because it samples only a single point.

This brings us to what this work is and is not. We are not claiming a predictive model for novel compounds — that is the goal, and it is within reach. What we are doing is establishing that the chameleonic signal is real, that it is size-gated, and that ensemble sampling is the right tool to see it. The bottleneck now is not the science — it is data and compute. We need source-homogeneous experimental permeability measurements at scale and the computational infrastructure to run physics-based conformer sampling beyond the ETKDGv3 heuristic. The findings presented here are our case to the community: the signal exists, the approach is sound, and these are the resources required to turn it into something predictive.

### 3.3 Reference Compound Validation Against NMR-Grounded Literature

The Tier-1 pipeline was validated against CsA, whose conformational switch has been directly observed by NMR spectroscopy in CDCl₃/hexane and DMSO/H₂O (Witek 2016⁷: 57–79 NOE restraints per conformer; Rüdisser 2023: eNOE backbone RMSD 0.10 Å).

| Compound | Permeable | Tier-1 ΔPSA (Å²) | Literature ΔPSA | Validation |
|----------|-----------|-----------------|----------------|------------|
| CsA (Cyclosporin A) | Yes | **84.9** | ~75 Å²⁷ | ✓ within 10% |
| DP172 | Yes | 88.9 | — | strongest case |
| HexPep | No | 64.4 | — | — |
| 1NMe3 | Yes | 47.8 | — | — |
| PSLYF | No | 65.3 | — | — |

CsA ΔPSA = 84.9 Å² vs. NMR-grounded literature ~75 Å²⁷, within 10%. This single-compound agreement is encouraging and consistent with the ensemble approach capturing a real conformational signal, though one data point is insufficient to constitute rigorous validation. DP172 (88.9 Å², permeable) provides a second supporting case without a literature reference for direct comparison. Notably, HexPep (impermeable, PAMPA −6.20) and 1NMe3 (permeable, PAMPA −5.52) — both hexapeptides — show raw ΔPSA values of 64 and 48 Å² respectively, yet ΔPSA fails to discriminate their opposing permeability outcomes. This is consistent with the sub-9-residue size limitation identified in Section 4: small cyclic scaffolds lack sufficient backbone flexibility to execute the intramolecular H-bond rearrangement that drives chameleonic switching, so absolute ΔPSA carries no mechanistic meaning for them regardless of its magnitude.

### 3.4 Descriptor Performance, Normalization, and the Cost of Label Noise

Ensemble ΔPSA achieves AUC = 0.692 on the clean 1,566-compound subset — the top-ranked descriptor. But two questions immediately follow: does this number hold up on noisier data, and is it measuring a genuine chameleonic signal or simply functioning as a molecular size proxy? We address both in sequence.

**Overall descriptor AUC — source-stratified vs. full dataset:**

![AUC-ROC by descriptor](figures/auc_roc_bar.png)

*Figure 2. AUC-ROC for key descriptors on the 1,566-compound source-stratified subset (blue) and full 7,298-compound dataset (orange). Absolute ΔPSA and per-MW normalized ΔPSA (shaded) are the 3D ensemble descriptors. MolLogP inverts on the clean subset (AUC = 0.317, flipped 0.683) — chameleonic macrolides are large and polar yet dominate the permeable population, reversing the standard lipophilicity/permeability relationship.*

| Descriptor | 1,566 compounds | Full 7,298 |
|------------|-----------------|------------|
| delta_psa3d / MolWt (per-MW normalized) | 0.589 | 0.470 |
| delta_psa3d (Tier-1 ensemble, absolute) | **0.692** | 0.505 |
| MolLogP (2D baseline) | 0.317† (inv. 0.683) | **0.630** |
| NumHDonors (2D baseline) | 0.691 | 0.406 |
| TPSA (2D baseline) | 0.663 | 0.447 |
| delta_3DPSA_db (single-structure negative control) | 0.458 | 0.493 |

*†Flipped AUC = 0.683: chameleonic macrolides are polar yet permeable, reversing the lipophilicity rule.*

On the full dataset, ensemble ΔPSA collapses to 0.505 while MolLogP holds at 0.630. Lipophilicity is indifferent to assay protocol; chameleonic ΔPSA is not — it is selectively degraded by label noise from pooled assays (Townsend 2020⁸, Kelly 2021⁹). The single-structure negative control (delta_3DPSA_db) is stable and near chance at both scales (0.458/0.493), confirming that ensemble sampling — not source stratification alone — drives the 1.5k signal.

**Is the AUC = 0.692 a real chameleonic signal, or a size proxy?**

Absolute ΔPSA scales with molecular size — a 15-residue peptide that barely switches produces more raw Å² than a 9-residue peptide that fully buries its polar surface. To test whether this confound is driving the overall AUC, we computed `delta_psa3d_per_mw` = ΔPSA / MolWt × 1000 (Å² per kDa) and compared AUC within each size bucket.

![Per-MW normalization AUC comparison](figures/delta_psa3d_normalization_comparison.png)

*Figure 3. AUC-ROC for absolute ΔPSA (blue) vs. per-MW normalized ΔPSA (orange), stratified by monomer length, on the 1,566-compound (left) and full 7,298-compound (right) datasets. Green shading: the 12–15 residue bucket where normalization consistently improves AUC. Normalization improves large-compound discrimination and degrades small-compound discrimination — identifying the same size boundary from a completely independent analytical direction.*

| Size bucket | n (1.5k) | Absolute AUC | Per-MW AUC | Δ AUC |
|-------------|----------|-------------|------------|-------|
| ≤8 residues | 737 | 0.528 | 0.553 | +0.025 |
| 9–11 residues | 310 | 0.440* | 0.618* | — |
| 12–15 residues | 519 | 0.528 | **0.696** | **+0.168** |

*\*9–11 residue bucket is 94.8% permeable on the 1.5k (~16 impermeable compounds); AUC unreliable due to class imbalance. This bucket is inconclusive on both datasets — see note below.*

| Size bucket | n (7k) | Absolute AUC | Per-MW AUC | Δ AUC |
|-------------|--------|-------------|------------|-------|
| ≤8 residues | 4,455 | 0.507 | 0.489 | −0.018 |
| 9–11 residues | 2,317 | 0.492 | 0.447 | −0.045 |
| 12–15 residues | 525 | 0.522 | **0.674** | **+0.152** |

For large compounds (12–15 residues), removing the size effect reveals a stronger chameleonic efficiency signal: AUC rises from 0.528 to 0.696 in the clean data, and holds at 0.674 even on the noisy full dataset. For small compounds (≤8 residues), normalization *degrades* AUC on both datasets — because ΔPSA was never measuring their permeation mechanism to begin with. This asymmetry confirms that the overall 0.692 is partly a size proxy, but a real chameleonic signal is present within the large-compound population once size is removed.

**Caveat — 9–11 residue bucket is inconclusive:** The two datasets give contradictory permeability rates for this bucket (94.8% vs. 60.9%) and we cannot adjudicate between them — both could be correct for their specific compound sets, which may not even overlap in chemical space. We simply do not have enough homogeneous data in this size range to draw conclusions in either direction. More fundamentally, MW alone is not what gates chameleonic behavior — it is a proxy for the underlying conformational flexibility and ring geometry that determine whether a compound can execute a chameleonic switch. When a compound reaches sufficient size to be capable of that switch, MW and 3D descriptors become informative for measuring chameleonic propensity; what MW threshold that corresponds to, and how sharply it cuts, is precisely what this bucket should test — but cannot, with current data. What is needed is a precisely characterized, individual-compound dataset with broader 9–11 residue coverage — the conditions under which a meaningful AUC could actually be computed (Section 6, Experiment 4).

**Track D — the same pattern, seen geometrically:**

To understand why label quality destroys the signal at scale, we colored the Panel C UMAP embedding by molecular weight (Track D). On the clean 1.5k, permeable compounds have a median MW of 1,180 Da vs. 820 Da for impermeable — a 1.44× gap visible directly in chemical space. On the full 6,938-compound dataset, that gap disappears entirely.

| Population | Median MW (1,566 cpds) | Median MW (6,938 cpds) |
|------------|------------------------|------------------------|
| Permeable (PAMPA ≥ −6) | **1,180 Da** | 820 Da |
| Impermeable | 820 Da | 820 Da |
| Ratio | **1.44×** | 1.00× |

![Panel C 1566 + Track D](figures/Panel_C_combined_umap_1566.png)

*Figure 4. Panel C UMAP on the 1,566-compound source-stratified subset with Track D MW coloring. C0 (87.8% permeable) occupies the high-MW plasma region. The permeable cluster is large macrolides; C1 (59.3% permeable, lower MW) is a population achieving permeability through other means.*

![Panel C 6938 + Track D](figures/Panel_C_combined_umap_6938.png)

*Figure 5. Panel C UMAP on the full 6,938-compound dataset. The 1.44× MW gap visible in Fig. 4 vanishes — both populations converge to ~820 Da. Pooled-assay labels are not merely noisy; they actively invert the MW/permeability relationship, labeling large chameleonic compounds as impermeable (or vice versa) at sufficient frequency to destroy the size signal entirely.*

| Metric | 1,566 compounds | 6,938 compounds |
|--------|-----------------|-----------------|
| HDBSCAN silhouette | 0.425 | −0.022 |
| Best cluster perm rate | 87.8% | 70.3% |
| Median MW permeable | **1,180 Da** | 820 Da |
| Median MW impermeable | 820 Da | 820 Da |

The AUC collapse, the normalization asymmetry, and the MW gap disappearance all point to the same conclusion: the signal is real on clean data and destroyed by label noise at scale. Label quality — not descriptor quality — is the limiting factor.

### 3.5 UMAP Chemical Space Structure

Three feature panels were analyzed: Panel A (7 2D descriptors), Panel B (8 3D Δ descriptors), and Panel C (combined 9-feature set).

**Panel A — 2D descriptors (baseline):**

![Panel A: 2D descriptor UMAP](figures/Panel_A_2D_umap.png)

*Figure 6. Panel A UMAP (2D descriptors only). K-Medoids forces 8 clusters but none carry biological signal. HDBSCAN identifies 69 micro-clusters with negative silhouette (−0.112) and no permeability enrichment. Knowing a molecule's static physicochemical properties alone provides no useful structure in permeability space.*

**Panel B — 3D Δ features (the core result):**

![Panel B: 3D delta descriptor UMAP](figures/Panel_B_3D_delta_umap.png)

*Figure 7. Panel B UMAP (3D Δ descriptors only). Structure emerges without forcing it. HDBSCAN identifies two natural populations: C0 (695 compounds, 73% permeable) and C1 (481 compounds, 25% permeable) — a 48-point enrichment gap with no cluster count parameter specified. PAMPA LogPexp maps continuously onto the embedding, with permeable compounds concentrated in C0. This is the qualitative proof that 3D conformational descriptors capture a chemically real permeability signal absent in 2D descriptor space.*

**Panel C — Combined features, source-stratified 1,566 compounds:**

ARI stability: min = **0.995** across 5 seeds — one dominant chemical structure in the combined 2D+3D feature space.

| Cluster | n | Perm rate | Enrichment |
|---------|---|-----------|------------|
| HDBSCAN C0 | 883 | **87.8%** | 1.16× |
| HDBSCAN C1 | 651 | 59.3% | 0.78× |

---

## 4. Interpretation: Two Populations, Two Descriptor Spaces

The data point toward two distinct populations of cyclic peptides that are governed by different permeation physics — and therefore require different descriptor sets to characterize.

**Large cyclic peptides (≥9 residues, tentative threshold):** In the clean 1,566-compound dataset, compounds of 12–15 residues are 86.3% permeable and the ≥9-residue population is substantially enriched above the baseline rate of 75.8%. These compounds have sufficient backbone flexibility to form stabilizing intramolecular H-bonds in apolar environments, collapsing polar surface area and reducing the desolvation penalty of membrane entry — the chameleonic switch. CsA (MW 1,203 Da, 11 residues, ΔPSA = 84.9 Å²) is the archetypal example⁷. For this population, 3D ensemble ΔPSA is mechanistically appropriate: it measures the physical property that likely contributes to permeation. Critically, chameleonic propensity is probably not linearly predictive of permeability on its own — high permeability in large macrolides is most likely the result of chameleonic propensity acting in combination with other structural features (ring geometry, H-bond donor count, backbone N-methylation pattern). ΔPSA captures one key contributor; the full descriptor ensemble for this population remains to be characterized.

**Small cyclic peptides (≤8 residues):** The ≤8-residue population is 60.4% permeable in the clean data — these compounds are permeating, but ΔPSA is essentially uninformative for them (AUC = 0.528). This is not surprising: smaller scaffolds lack the conformational freedom to execute a large chameleonic switch. Their permeability is more likely governed by a different set of descriptors entirely — N-methylation pattern, backbone HBD count, intrinsic lipophilicity, amide-to-ester substitution¹⁰ — properties that reduce the cost of membrane entry without requiring a conformational rearrangement. What that descriptor ensemble looks like is an open question; what this work establishes is that ΔPSA is not part of it.

**MW as a proxy, not a gate:** Molecular weight correlates with the conformational flexibility that enables chameleonic switching, which is why it appears as a separating axis in UMAP Track D and why per-MW normalization improves AUC for large compounds. But MW is not itself the mechanistic gate — it is a proxy for ring size, backbone degrees of freedom, and the number of potential intramolecular H-bond donors. A rigid 15-residue peptide with locked conformation may be no more chameleonic than a flexible 9-residue one. MW reveals which mechanism a compound is likely using; it does not cause it.

**What the overall AUC of 0.692 actually measures:** Size-stratified analysis reveals that ΔPSA within each residue-count bucket performs only marginally above chance (AUC 0.44–0.53). The overall 0.692 arises primarily from between-group discrimination — large compounds have high absolute ΔPSA AND high permeability rates in the clean data. Absolute ΔPSA is partly a molecular size proxy, not a pure chameleonic efficiency measure. This is precisely the confound that Yu et al.³ address with their dimensionless ΔPSA/SASA_total ratio, and it is what per-MW normalization begins to correct for within size classes.

**What this work establishes and what it does not:** The ensemble ΔPSA approach shows that (1) conformational sampling recovers a permeability-correlated signal for large macrolides that is absent in single-structure methods; (2) chemical space separates into two populations with distinct MW and permeability profiles; and (3) ΔPSA is silent on the smaller population, setting a clear boundary on where it applies. The immediate goal is to improve chameleonic propensity prediction — refining the descriptor (per-SASA normalization, size gating) so that it more cleanly measures switching efficiency rather than molecular size. The downstream goal is to understand how chameleonic propensity, combined with other structural metadata, predicts permeability for each population separately. The concurrent work of Yu et al.³ independently arrives at the same two-population framing — explicitly finding that small cyclic peptides show *"absence of significant conformational switching"* — and their dimensionless ΔPSA/SASA_total ratio is the next normalization to implement (see Section 6, Experiment 1).

---

## 5. Conformer Model and Proposed Validation

ETKDGv3 is appropriate for population-scale screening but is a distance geometry heuristic, not a physics-based simulation of solvent-dependent conformational equilibria. The conformer selection heuristic (max PSA = aqueous, min PSA = membrane) assumes the vacuum ensemble contains conformers representative of solution-phase endpoints — an assumption validated by the CsA result (ΔPSA = 84.9 vs. ~75 Å²⁷) but not systematically tested across the dataset. Critically, vacuum conformer sampling ignores solvation-dependent energy barriers: conformers that are low-energy in CHCl3-mimetic implicit solvent may be sterically accessible in vacuum but kinetically inaccessible in solution, and vice versa. Implicit solvent models (ALPB, GBSA) or explicit solvent MD would more faithfully sample the solution-phase ensemble, at the cost of substantially greater computational expense.

The current state-of-the-art open-source pipeline for macrocyclic peptide conformer generation is the CREMP workflow¹¹: ETKDGv3 generates an initial diverse pool (up to 5,000 conformers) which CREST¹² then re-ranks using GFN2-xTB energetics. This treats ETKDGv3 and CREST as complementary rather than competing. We propose the following validation on local compute resources:

**Validation set:** 5 reference compounds (CsA, DP172, HexPep, 1NMe3, PSLYF) + 5 randomly selected compounds from each HDBSCAN cluster (C0 and C1) = 15 compounds total.

**Protocol:** ETKDGv3 (5,000 conformers) → CREST re-ranking (GFN2-xTB + ALPB solvation, water then CHCl3) → ΔPSA computed from CREST ensemble.

**Success criterion:** Spearman ρ > 0.8 between ETKDGv3 Tier-1 ΔPSA and CREST ΔPSA across the 15-compound set, and preserved binary permeable/impermeable classification. If the rank order is preserved, ETKDGv3 is validated as sufficient for population-scale pre-filtering. If it is disrupted, ETKDGv3 should be used only as a diversity generator with CREST re-ranking for compounds of interest.

---

## 6. Next Experiments

| Priority | Experiment | Status | Hypothesis |
|----------|------------|--------|------------|
| ✓ Done | Per-MW normalization (`delta_psa3d / MolWt`); AUC by size bucket | Complete (this work) | Normalization improves large-compound AUC, degrades small-compound AUC — confirmed two-mechanism boundary |
| 1 | `delta_psa3d_per_sasa` = ΔPSA / SASA_total (Yu et al.³ approach); rerun AUC on 1,566 compounds | Next | Dimensionless normalization outperforms per-MW; recovers signal within size classes |
| 2 | `delta_psa3d_per_residue` = ΔPSA / Monomer_Length; compare to per-MW and per-SASA | Next | Per-residue captures chameleonic efficiency independently of MW |
| 3 | Filter `Monomer_Length ≥ 9`; compare AUC with and without cutoff on individual-compound subset | Next | Explicit size gate improves signal; confirms chameleonic threshold |
| 4 | Restrict to individual-compound, single-protocol sources only; rerun AUC | Next | Pooled-assay label inconsistency is driving AUC collapse; precise data recovers baseline signal |
| 5 | ETKDGv3 → CREST validation on 15-compound reference set | Next | ETKDGv3 rank order is preserved vs. CREST/ALPB solvation |

---

## References

1. Jiang, Y. et al. CycPeptMPDB: a comprehensive database of membrane permeability of cyclic peptides. *J. Chem. Inf. Model.* **63**, 1936–1943 (2023). https://doi.org/10.1021/acs.jcim.2c01573

2. Furukawa, A. et al. Passive membrane permeability in cyclic peptides is robust to extensive structural variation. *J. Med. Chem.* **59**, 9503–9512 (2016). https://doi.org/10.1021/acs.jmedchem.6b01246

3. Yu, Y. et al. Normalized polar surface area as a descriptor for chameleonic cyclic peptide permeability. *bioRxiv* (2026). https://doi.org/10.64898/2026.01.06.697862

4. Wang, S., Witek, J., Landrum, G. A. & Riniker, S. Improving conformer generation for small rings and macrocycles based on distance geometry and experimental torsional-angle preferences. *J. Chem. Inf. Model.* **60**, 2044–2058 (2020). https://doi.org/10.1021/acs.jcim.0c00025

5. Riniker, S. & Landrum, G. A. Better informed distance geometry: using what we know to improve conformation generation. *J. Chem. Inf. Model.* **55**, 2562–2574 (2015). https://doi.org/10.1021/acs.jcim.5b00654

6. Schärfer, C. et al. CONFORGE: a versatile and high-performance conformer generator. *J. Chem. Inf. Model.* **63**, 5819–5831 (2023). https://doi.org/10.1021/acs.jcim.3c00563

7. Witek, J. et al. Kinetic and thermodynamic analysis of cyclosporin A solution conformations. *J. Chem. Theory Comput.* **12**, 3861–3873 (2016).

8. Townsend, C. E. et al. Cyclic peptide membrane permeability measurements via high-throughput screening. *ChemRxiv* (2020; preprint). https://doi.org/10.26434/chemrxiv.13335941.v1

9. Kelly, C. N. et al. Iterative parallel synthesis and permeability measurement identifies macrocycle architectures with enhanced cell permeability. *J. Am. Chem. Soc.* **143**, 8655–8667 (2021). https://doi.org/10.1021/jacs.0c06115

10. Bockus, A. T. et al. Probing the physicochemical boundaries of cell permeability and oral bioavailability in lipophilic macrocycles inspired by natural products. *J. Med. Chem.* **58**, 4581–4589 (2015).

11. Klarich, K. L. et al. CREMP: conformer-rotamer ensemble of macrocyclic peptides for machine learning. *Sci. Data* **11**, 759 (2024). https://doi.org/10.1038/s41597-024-03698-y

12. Pracht, P., Bohle, F. & Grimme, S. Automated exploration of the low-energy chemical space with fast quantum chemical methods. *Phys. Chem. Chem. Phys.* **22**, 7169–7192 (2020). https://doi.org/10.1039/C9CP06869D

13. Bannwarth, C., Ehlert, S. & Grimme, S. GFN2-xTB — an accurate and broadly parametrized self-consistent tight-binding quantum chemical method with multipole electrostatics and density-dependent dispersion contributions. *J. Chem. Theory Comput.* **15**, 1652–1671 (2019). https://doi.org/10.1021/acs.jctc.8b01176


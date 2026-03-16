# Reference Compound Selection — Tier-2 Validation Set

**Project:** CHEM 269 Final Project — Cyclic Peptide Dual-Dielectric Pipeline
**Author:** Jorge Carmona
**Date:** March 2026
**Purpose:** Justify the 5-compound reference set used for Tier-2 OMEGA validation and Tier-1 cross-check.

---

## Scientific Rationale

The central hypothesis of this project is that **3D conformational Δ features** — specifically ΔPSA computed across two dielectric environments (aqueous ε=78, membrane ε=4) — can quantify chameleonic potential and predict PAMPA permeability better than 2D descriptors alone.

To test this with Tier-2 OMEGA, we need reference compounds where the chameleonic mechanism has been **directly proven experimentally**, not merely inferred from PAMPA numbers. The 5 compounds below were selected to satisfy three criteria:

1. **Experimental proof of conformational behavior** (NMR, MD, or X-ray)
2. **Unambiguous PAMPA signal** (reproduced across independent labs where possible)
3. **Coverage of the full permeability range** and both chameleonic / non-chameleonic structural archetypes

---

## Final 5 Reference Compounds

### 1. Hexapeptide — compound.1

| Property | Value |
|---|---|
| CycPeptMPDB ID | 2 |
| Primary Source | Rezai, T.; Yu, B.; Millhauser, G.L.; Jacobson, M.P.; **Lokey, R.S.** *J. Am. Chem. Soc.* **2006**, *128*, 2510–2511. |
| Cross-validated by | Wang 2015 (1 lab) |
| PAMPA | −6.20 log cm/s |
| Classification | **Impermeable** (below −6.0 threshold) |
| MW / HBD / LogP | 713 Da / 6 / 2.55 |
| Ring size | Cyclic hexapeptide (6 residues) |
| DB ΔPSA (H₂O − CHCl₃) | +2.0 Å² |

**Experimental evidence for inclusion:**
Rezai et al. 2006 used ¹H NMR in CDCl₃ to show that compound.1 (the fully unprotected hexapeptide) displays **exposed backbone amides with no evidence of intramolecular H-bonding**. The molecule cannot adopt a compact conformation in a membrane-like dielectric environment — it is conformationally rigid in the sense that all NH donors remain solvent-exposed. PAMPA = −6.20 confirms non-permeability.

**Role in the set:** Impermeable, non-chameleonic baseline. Forms a matched pair with 1NMe3 — same scaffold, same ring size, directly comparable.

---

### 2. N-Me Hexapeptide — 1NMe3

| Property | Value |
|---|---|
| CycPeptMPDB ID | 980 |
| Primary Source | White, T.R.; Renzelman, C.M.; Rand, A.C.; ... **Lokey, R.S.** *Nature Chemical Biology* **2011**, *7*, 810–817. DOI: 10.1038/nchembio.664 |
| Cross-validated by | 9 independent labs (Rand 2012, Hewitt 2015, Lewis 2015, Marelli 2015, Nielsen 2015, Wang 2015, Hickey 2016, Naylor 2018, Bockus 2015) |
| PAMPA | −5.52 (primary); range −4.40 to −6.40 across labs |
| Classification | **Permeable** |
| MW / HBD / LogP | 755 Da / 3 / 3.57 |
| Ring size | Cyclic hexapeptide (6 residues), 3× N-methylated |
| DB ΔPSA (H₂O − CHCl₃) | +1.0 Å² |

**Experimental evidence for inclusion:**
White & Lokey 2011 performed **NMR NOESY** analysis in CDCl₃ showing formation of intramolecular H-bonds in the N-methylated analogs. N-methylation at positions 1, 2, and 3 (1NMe3) eliminates 3 NH donors, enabling the remaining backbone to fold into a compact conformation that shields polar surface from the membrane environment. This is the direct experimental proof of the chameleonic mechanism: the same hexapeptide scaffold becomes chameleonic upon N-methylation.

**Role in the set:** Permeable, chameleonic counterpart to compound.1. The most experimentally reproduced cyclic peptide permeability reference in the literature (9 labs). The compound.1 → 1NMe3 pair is the **canonical demonstration that N-methylation drives chameleonic switching** and is the core scientific story of this reference set.

---

### 3. Cyclosporin A (CsA)

| Property | Value |
|---|---|
| CycPeptMPDB ID | 1 (primary); also 932, 1822, 2356, 7188, 7353 |
| Primary Sources | Rezai & Lokey, *JACS* 2006; Witek et al., *J. Chem. Theory Comput.* **2016**, *12*, 4025–4039 |
| Cross-validated by | 8 independent labs across 6 CycPeptMPDB entries |
| PAMPA | −6.60 (Rezai 2006 primary entry); mean = **−5.90** across all 6 measurements |
| Classification | Borderline — **use mean PAMPA = −5.90 (permeable)** |
| MW / HBD / LogP | 1203 Da / 5 / 3.27 |
| Ring size | Cyclic undecapeptide (11 residues), 7× N-methylated |
| DB ΔPSA (H₂O − CHCl₃) | −1.0 Å² (entry ID=1); range −26 to +1 across entries |
| Oral bioavailability (human) | ~29% |

**Experimental evidence for inclusion:**
CsA has **three independent levels of experimental proof** for chameleonic behavior:

1. **NMR (1990–1994):** Kessler et al. (1990, *JACS*) and Wenger (1994, *Angew. Chem.*) demonstrated that CsA adopts two distinct conformations — a compact conformation in CDCl₃ with intramolecular H-bonds shielding 5 backbone NHs, and an extended conformation in DMSO/H₂O. These are among the first direct NMR proofs of chameleonic behavior in any macrolide.
2. **MD simulations (2016):** Witek et al. *JCTC* 2016 quantified the conformational ensemble using replica exchange MD, showing a ΔPSA of **~75 Å²** between the membrane-form (closed) and aqueous-form (open) conformations. This provides the ground truth ΔPSA magnitude that Tier-1 (ETKDG) and Tier-2 (OMEGA) should recapitulate.
3. **Clinical validation:** 29% oral bioavailability in humans despite MW=1203 and HBD=5 (both violate Lipinski Ro5) — the definitive real-world proof that chameleonic behavior enables oral absorption of bRo5 compounds.

**Important caveat:** The Rezai 2006 primary CycPeptMPDB entry (ID=1) records PAMPA = −6.60, classifying CsA as *impermeable* by the −6.0 threshold. The mean across 6 measurements is −5.90 (permeable). For analysis, **use the cross-lab mean** and note the measurement variability. The −6.60 value reflects PAMPA assay sensitivity at the detection limit for this large compound.

**Role in the set:** Gold-standard chameleonic benchmark. The most studied macrocycle in the permeability literature. Provides the only quantitative MD-based ΔPSA ground truth (~75 Å²) against which both Tier-1 and Tier-2 can be directly validated.

---

### 4. DP-172

| Property | Value |
|---|---|
| CycPeptMPDB ID | 183 |
| Primary Source | Chugai Pharmaceutical Co. internal screen; reported in Fouché, M. et al. *J. Med. Chem.* **2016** (CHUGAI 2013 dataset in CycPeptMPDB) |
| Cross-validated by | Ohta 2023 (1 additional lab) |
| PAMPA | −4.15 log cm/s |
| Classification | **Strongly permeable** |
| MW / HBD / LogP | 1243 Da / 5 / 3.00 |
| Ring size | Cyclic undecapeptide (11 residues) |
| DB ΔPSA (H₂O − CHCl₃) | **−47.0 Å²** (H₂O=155, CHCl₃=202) |

**Experimental evidence for inclusion:**
DP-172 represents the pharmaceutical industry's proof-of-concept for chameleonic bRo5 macrolides. With MW=1243 and HBD=5, it should be completely impermeable by Lipinski/bRo5 rules — yet PAMPA = −4.15 (strongly permeable). Its DB ΔPSA = −47 Å² is the largest absolute chameleonic signal of any compound in the reference set and is available directly from the static CycPeptMPDB 3DPSA data, making it the **strongest single-compound argument for why 3D descriptors outperform 2D** in predicting bRo5 permeability.

**Role in the set:** High-permeability anchor (PAMPA = −4.15); only compound with an unambiguously large DB ΔPSA; demonstrates generalizability beyond the Lokey hexapeptide scaffold; provides pharmaceutical industry validation that chameleonic compounds reach drug development pipelines. No published NMR conformational data (limitation noted below).

---

### 5. c\*[PSLYF] — compound 9

| Property | Value |
|---|---|
| CycPeptMPDB ID | 1829 |
| Primary Source | Hickey, J.L. et al. *J. Med. Chem.* **2016**, *59*, 5342–5358. |
| Cross-validated by | No independent measurements |
| PAMPA | −9.10 log cm/s |
| Classification | **Strongly impermeable** |
| MW / HBD / LogP | 764 Da / 8 / 0.42 |
| Ring size | Cyclic pentapeptide (5 residues) |
| DB ΔPSA (H₂O − CHCl₃) | **0.0 Å²** (H₂O=203, CHCl₃=203) |

**Experimental evidence for inclusion:**
Hickey et al. 2016 conducted a systematic structure-permeability study of cyclic pentapeptides in which they specifically mapped backbone HBD count, N-methylation pattern, and ring rigidity to PAMPA. c\*[PSLYF] represents the maximal-HBD, no-N-methylation extreme: 8 H-bond donors, LogP=0.42 (highly hydrophilic), and a rigid pentapeptide ring. DB ΔPSA = 0 Å² confirms it has no chameleonic potential. PAMPA = −9.10 is the most negative value in the reference set (by 0.95 log units from the next compound), placing it unambiguously at the non-permeable floor.

**Role in the set:** Strongly impermeable, non-chameleonic negative control. Maximum contrast to CsA and 1NMe3. Anchors the bottom of the PAMPA axis. The zero DB ΔPSA paired with PAMPA = −9.10 provides the clearest possible negative data point for the chameleonic hypothesis.

---

## Summary Table

| # | ID | Compound | Source | PAMPA | Permeable | DB ΔPSA | HBD | MW | NMR/MD Proof |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 2 | Hexapeptide (compound.1) | Rezai & Lokey, *JACS* 2006 | −6.20 | No | +2 Å² | 6 | 713 | NMR (CDCl₃) |
| 2 | 980 | N-Me Hexapeptide (1NMe3) | White & Lokey, *Nat Chem Biol* 2011 | −5.52 | Yes | +1 Å² | 3 | 755 | NMR NOESY (CDCl₃) |
| 3 | 1 | Cyclosporin A | Rezai 2006 / Witek, *JCTC* 2016 | −5.90 (mean) | Yes | −1 Å² | 5 | 1203 | NMR + MD (~75 Å²) |
| 4 | 183 | DP-172 | CHUGAI 2013 | −4.15 | Yes | −47 Å² | 5 | 1243 | None published |
| 5 | 1829 | c\*[PSLYF] (9) | Hickey, *J Med Chem* 2016 | −9.10 | No | 0 Å² | 8 | 764 | Systematic SAR |

---

## Limitations of This Reference Set

### 1. Small n for statistical analysis
Five compounds cannot support a statistically meaningful Pearson or Spearman correlation. Any r or ρ reported from this set alone should be treated as **directional evidence only**, not a significance test. The Tier-2 validation is explicitly a mechanistic cross-check, not a regression.

### 2. DP-172 lacks published conformational data
DP-172 is from a pharmaceutical internal screen (CHUGAI). Unlike compounds 1–3, there is **no published NMR or MD study** confirming its chameleonic mechanism. The large DB ΔPSA (−47 Å²) is compelling but is a static single-structure computation. Its inclusion in Tier-2 OMEGA is therefore a **prediction to be tested**, not a validation against known ground truth.

### 3. CsA PAMPA measurement heterogeneity
CsA PAMPA ranges from −6.60 to −5.01 across 6 labs (1.6 log units). This is the largest inter-lab spread in the set and reflects genuine sensitivity issues at the PAMPA detection limit for large macrolides. Using any single measurement introduces bias. **Use the cross-lab mean (−5.90) and report the range as uncertainty.**

### 4. c\*[PSLYF] has no independent replication
Compound 5 has a single PAMPA measurement (Hickey 2016). At PAMPA = −9.10, it is at or below many PAMPA assay detection limits (~−9 log cm/s), meaning the true value may be even more negative. This is unlikely to change its classification but the absolute value has higher uncertainty than compounds 2–4.

### 5. Scaffold bias toward peptidic macrolides
All 5 compounds are cyclic peptides (amide/N-methylamide backbones). The hypothesis may not generalize to depsipeptides, macrolactones, or heterocyclic macrocycles without additional references. Conclusions should be scoped to **cyclic peptide permeability**, not macrocycles in general.

### 6. No compound in the PAMPA range −5.0 to −5.5 with large ΔPSA
There is a gap in the reference set between 1NMe3 (−5.52, ΔPSA ≈ +1) and DP-172 (−4.15, ΔPSA = −47). A compound at PAMPA ≈ −5.0 to −5.2 with intermediate ΔPSA would strengthen the dose-response relationship between chameleonic potential and permeability. This gap is acceptable for a course project but would need to be addressed in a full publication.

### 7. DB ΔPSA does not capture ensemble chameleonism for most compounds
Four of five compounds show DB ΔPSA ≈ 0 to +2 Å². CsA — the canonical chameleonic reference — shows DB ΔPSA = −1 Å², despite a known ~75 Å² conformational switch from MD. This is expected (the DB value is a single static structure), but it means the **DB feature alone cannot validate the chameleonic hypothesis** for this reference set. The Tier-1 and Tier-2 ensemble methods are essential, not optional.

---

## What This Set Proves vs. What It Cannot Prove

| Can prove | Cannot prove |
|---|---|
| N-methylation drives chameleonic switching (compounds 1→2, NMR-confirmed) | Generalizability to non-peptidic macrocycles |
| Large ΔPSA correlates with high permeability (DP-172) | Statistical significance (n=5 is insufficient) |
| Zero ΔPSA corresponds to non-permeability (c\*[PSLYF]) | Causality (correlation ≠ mechanism without the NMR/MD data) |
| Tier-1 ETKDG directionally captures known chameleonic trends | Quantitative accuracy of Tier-1 absolute ΔPSA values |
| DB static 3DPSA misses ensemble chameleonism (CsA) | That OMEGA is superior to all other conformer methods |

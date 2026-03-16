# CSO Scientific Review — CHEM 269 Final Project
## Dual-Dielectric 3D Conformational Descriptors for Cyclic Peptide Permeability

**Reviewer:** Chief Scientific Officer (Claude Sonnet 4.6)
**Date:** March 15, 2026
**Status:** Pre-production — overnight Tier-1 run pending

---

## Executive Summary

The pipeline operationalizes a legitimate and well-grounded scientific hypothesis with a technically sound implementation. The core chameleonic framing is experimentally defensible and the Tier-1 / Tier-2 two-tier design is a reasonable scope for a course project. However, three issues require direct acknowledgment in the final writeup: (1) the Tier-1 vacuum heuristic has a known artifact for rigid impermeable peptides that inflates ΔPSA in a direction that partially undermines the hypothesis; (2) the ΔPSA metric is used as both a primary feature and a validation target against Witek 2016, but the methodological gap between vacuum max-min sampling and explicit-solvent Boltzmann-weighted MD is large enough that quantitative agreement should not be claimed; and (3) the compound.1 / 1NMe3 separation at n=20 conformers (3.5 Å²) is scientifically fragile — this specific pair is the narrative centerpiece of the reference set but is the weakest Tier-1 result. All three issues are manageable within the project scope if framed honestly.

---

## 1. Scientific Validity

### 1.1 Chameleonic Operationalization in Tier-1

The heuristic of treating the max-PSA conformer as the aqueous form and the min-PSA conformer as the membrane form is the standard approximation used in the literature when explicit-solvent ensembles are unavailable (e.g., Naylor et al. 2018 use a closely analogous procedure). It is defensible in a course context provided it is presented as an approximation rather than as equivalent to the Witek protocol. The heuristic implicitly assumes that the conformer with the largest polar exposure is thermodynamically stabilized by aqueous solvation and that the conformer with the smallest polar exposure is the one stabilized by the membrane. Both assumptions have physical backing but are not guaranteed — particularly for rigid molecules where both the max-PSA and min-PSA conformers are vacuum artifacts rather than thermodynamically meaningful states.

**Verdict:** Acceptable approximation, but must be explicitly labeled as such in every figure caption and table header that uses ΔPSA from Tier-1.

### 1.2 ΔPSA vs PSA_mem as the Primary Metric

The scientific framing correctly identifies that Witek et al. 2016 found PSA_mem to be the strongest single predictor of passive permeability, superior to ΔPSA alone. The pipeline computes mem_psa3d and exports it; there is no obstacle to making PSA_mem (labeled as `mem_psa3d`) the lead metric in the correlation analysis, with ΔPSA reported as a secondary measure of conformational flexibility. Reporting PSA_mem as the headline predictor would more accurately reflect the Witek finding and would be more defensible in an oral presentation. The framing document does acknowledge this in the feature table but the narrative arc foregrounds ΔPSA. This ordering should be reconsidered before the final writeup.

**Verdict:** Reframe PSA_mem as co-primary metric alongside ΔPSA. If only one number is cited in the conclusion, it should be PSA_mem.

### 1.3 Bondi Heavy-Atom PSA vs Ertl 2000 TPSA Convention

The heavy-atom-only convention used in `_polar_sasa()` (N, O, S, P, no polar-H) is internally consistent and is correctly documented in both `sasa_research_findings.md` and the code docstring. The key implication is that all absolute PSA values will be systematically lower than Witek 2016 values (which include polar-H) by an estimated 5–15 Å² per NH/OH group. For CsA with 5 NH donors, this could represent a 25–75 Å² systematic offset in absolute PSA, which means the absolute `aq_psa3d = 130.7 Å²` should not be directly compared to Witek's water-ensemble PSA without applying a correction or noting the methodological difference.

Critically, the relative ΔPSA is unaffected by this offset as long as the convention is applied consistently across all conformers — which it is. The Ertl 2000 TPSA convention is fragment-based and does not have an exact analog in 3D SASA; using Bondi heavy atoms is the community-standard workaround and is correct here.

**Verdict:** Scientifically sound. One sentence of disclosure is required in the methods section.

### 1.4 Vacuum Ensemble and Conformational Space Coverage

The vacuum ETKDGv3 + MMFF94s ensemble does not model solvation, which means it samples the intrinsic torsional and ring-flip preferences of the molecule without the thermodynamic bias imposed by the environment. This is the root cause of the c*[PSLYF] artifact (see Section 3.2). For chameleonic molecules, the vacuum ensemble is an imperfect but workable proxy because the chameleonic conformations span a wide PSA range and the max and min conformers bracket the aqueous and membrane forms reasonably well. For conformationally rigid molecules that are non-chameleonic, the vacuum ensemble will produce spurious PSA variation because ETKDGv3 will still generate diverse ring pucker conformers even for a molecule that cannot access low-PSA states in solution.

This is the most fundamental limitation of the Tier-1 approach and the primary reason Tier-2 CREST+ALPB exists.

### 1.5 The c*[PSLYF] Artifact: What It Means for the Hypothesis

c*[PSLYF] (ID=1829) produced ΔPSA = 67 Å² at n=20, against an expected value of ~0 Å². This is a false positive for chameleonic potential. In the context of the full 7,298-compound dataset, this type of artifact will affect all molecules that are conformationally rigid and non-chameleonic — the exact class the hypothesis predicts should have low ΔPSA. If the Tier-1 ΔPSA metric assigns high values to both truly chameleonic compounds (correct) and rigid impermeable peptides (false positive), the discriminating power of the metric will be diluted.

This does not invalidate the project. It means the analysis section must carefully distinguish between psa3d_spread (which can flag molecules with large vacuum conformational diversity regardless of mechanism) and the hypothesis-driven prediction (that permeable compounds should have higher ΔPSA than impermeable ones). Increasing n_confs to 200 will partially mitigate the artifact by better populating the extended, high-PSA conformers that dominate in aqueous solution, which should push the min-PSA conformer of c*[PSLYF] upward toward its physical value.

---

## 2. Reference Set Assessment

### 2.1 Is the Reference Set Well Chosen?

The five compounds provide good coverage of the permeability space: one strongly impermeable rigid control (c*[PSLYF], PAMPA = −9.1), one borderline impermeable non-chameleonic peptide (compound.1, PAMPA = −6.2), one borderline permeable canonical chameleonic compound (1NMe3, PAMPA = −5.5), one well-characterized large chameleonic macrolide (CsA, PAMPA = −5.9), and one strongly permeable pharmaceutical compound (DP-172, PAMPA = −4.15). The inclusion of compounds with NMR or MD-proven chameleonic behavior (compounds 1–3) as the mechanistic anchors is correct and scientifically sound.

The limitations are documented appropriately: n=5 is insufficient for statistical analysis, DP-172 lacks published conformational data, and CsA has inter-lab PAMPA variability of 1.6 log units. These are honest limitations and do not compromise the project's viability.

### 2.2 The compound.1 / 1NMe3 Pair at n=20

The 3.5 Å² separation (39.3 vs 35.8 Å² ΔPSA) at n=20 conformers is scientifically insufficient to make a clean claim about chameleonic discrimination between these two compounds. The framing document correctly identifies this pair as "the canonical demonstration that N-methylation drives chameleonic switching," and the N-methylation story is experimentally proven by NMR (White & Lokey 2011). However, the Tier-1 numerical result at n=20 does not recapitulate this story clearly. At n=200 conformers, the separation should widen, and if it does not, the result should be reported as: "Tier-1 vacuum sampling does not clearly discriminate this pair; Tier-2 CREST+ALPB provides the mechanistic validation."

This pair's story cannot rest on Tier-1 alone. It depends on Tier-2 delivering a clean result.

### 2.3 CsA ΔPSA = 59.6 Å² vs Witek's ~75 Å²

A 15 Å² shortfall (80% of the literature value) from 20 vacuum conformers is within the expected range for a Tier-1 approximation. The direction is correct, the magnitude is plausible, and the deficit is fully explained by: (a) the polar-H exclusion (~5–15 Å² for 5 NH groups), (b) fewer conformers than the Witek MD trajectory, and (c) the absence of explicit solvation bias. This result should be presented as "consistent with" the Witek value, not as "recapitulating" it. At n=200 conformers the value should approach or exceed 65 Å², which would be a reasonable Tier-1 approximation.

---

## 3. Methodological Soundness

### 3.1 ETKDGv3 + MMFF94s for Cyclic Peptides

ETKDGv3 is the current standard for 3D conformer generation in macrocycles within RDKit. The `useMacrocycleTorsions = True` flag enables the macrocycle-aware torsion potential introduced in ETKDGv3 specifically for ring systems of this size. MMFF94s (the "static" variant of MMFF94) is appropriate for energy minimization of peptide-like structures. The pruneRmsThresh = 0.5 Å setting is appropriate for small cyclic peptides but may undersample large macrolides (DP-172, CsA); the fallback to 1.0 Å is a reasonable safeguard. For a course project, this is a fully defensible conformer generation strategy. OMEGA would be preferred in an industrial setting, but the CREST+ALPB Tier-2 run for the 5 reference compounds provides a methodologically superior comparison point.

### 3.2 The Vacuum Artifact: When to Escalate to Tier-2

The c*[PSLYF] artifact is inherent to any vacuum-only approach. It should be flagged explicitly in the methods section rather than treated as a data quality issue. The appropriate response in the analysis is: (a) note that psa3d_spread (the full range) conflates chameleonic potential with vacuum conformational sampling noise; (b) show that among the reference compounds, the artifact is concentrated in the non-chameleonic, high-HBD class; and (c) use Tier-2 CREST+ALPB for c*[PSLYF] as the direct comparison showing the vacuum artifact disappears when environment-specific sampling is used.

### 3.3 H-Bond Counting Geometry

The topological path length cutoff of ≥6 atoms (excluding 1,2 through 1,4 contacts) is appropriate for macrocyclic peptides where short-range ring closure contacts are not physiologically meaningful H-bonds. The H...A distance cutoff of 3.0 Å is standard (Baker-Hubbard criterion uses 2.5 Å H...A, so 3.0 Å is slightly permissive, which is appropriate for computational sampling where conformer geometries may not be optimally hydrogen-bond-shaped). The D-H...A angle cutoff of 120° is standard. The geometry is defensible.

The zero ΔHB result across all 5 reference compounds at n=20 is a sampling artifact, not a methodological failure. At n=200, the signal should emerge. If it does not emerge for CsA at n=200, the angle or distance cutoff may need relaxation.

Note: the Tier-2 `count_hbonds_xyz()` function uses H...A < 2.5 Å rather than 3.0 Å, creating an inconsistency with the Tier-1 criterion. This should be harmonized before comparing Tier-1 and Tier-2 ΔHB values in the same table.

### 3.4 Tier-2 CREST+ALPB vs OMEGA+GB/SA

CREST+ALPB is clearly the more scientifically sound choice for this validation. The original OMEGA+GB/SA proposal used implicit solvent energy minimization rather than conformer sampling under implicit solvent, which is a weaker implementation of the dual-dielectric concept. CREST iMTD-GC with ALPB directly samples the conformational space in each dielectric environment using a GFN2-xTB Hamiltonian, producing environment-specific Boltzmann populations. This is much closer to the Witek 2016 protocol (which used explicit-solvent replica exchange MD) than OMEGA+GB/SA would have been.

The `--quick` and `--mquick` flags reduce sampling thoroughness in exchange for runtime feasibility. This is acceptable for a 48-hour deadline. The result will be less comprehensive than a full iMTD-GC run, but directionally valid for the reference compounds.

The `compute_psa_xyz()` function has been correctly updated from the ad-hoc contact-counting model to a proper FreeSASA-based approach (Path A) with a geometric fallback (Path B). Path B is still an approximation (spherical cap model without full Lee-Richards integration) but is far more accurate than the original `exposure = max(0.1, 1.0 - contacts * 0.12)` formula. This fix is adequate for the Tier-2 validation.

---

## 4. Deliverable Assessment

### Timeline (48 hours to Tuesday 11:59 PM)

**Tonight (Sunday → Monday morning):** Tier-1 full dataset run at n_confs=100–200. At approximately 2–5 seconds per molecule with 4–8 CPUs, 7,298 molecules at n_confs=100 represents roughly 4–10 hours of compute, which fits overnight.

**Monday:** Tier-2 CREST runs (5 compounds × 2 solvents = 10 CREST jobs). With `--quick` mode, each job should complete in 30–90 minutes for these molecular sizes. Running all 10 in parallel (2 per compound, staggered) is feasible in 2–4 hours on a multi-core machine.

**Monday–Tuesday:** Analysis, correlation tables, UMAP, final notebook. This is the most time-constrained component. The minimum viable analysis (correlation table + UMAP + reference compound overlay) should be achievable in 8–12 hours if the Tier-1 run succeeds.

**Minimum viable result if Tier-2 CREST fails:** The Tier-1 correlation analysis on 7,298 compounds, with the 5 reference compounds highlighted, is sufficient for the minimum viable deliverable. The Tier-2 mechanistic validation strengthens the story significantly but is not required for a passing grade given the backup plan documented in the scientific framing.

**Risk:** If the Tier-1 overnight run encounters widespread MMFF94s convergence failures or embedding failures (possible for very large or unusual macrocycles in the dataset), the fallback is to run at n_confs=50, which reduces compute time proportionally but retains directional signal in ΔPSA (ΔHB signal will be weaker).

---

## 5. Recommendations (Ordered by Scientific Impact)

**1. Run Tier-1 at n_confs=200, not 50.** At n=20, the c*[PSLYF] artifact is 67 Å² (should be ~0), and the compound.1/1NMe3 separation is only 3.5 Å². Both issues are expected to improve substantially at n=200 based on the ETKDGv3 conformational space sampling literature. The extra compute time (approximately 4× longer) is the single highest-impact investment before the deadline. If wall-clock time is the binding constraint, use n_confs=100 as the floor.

**2. Report PSA_mem (mem_psa3d) as a co-primary metric.** The Witek 2016 finding that PSA_mem is the best single predictor of permeability is the most direct experimental benchmark for this pipeline. If the correlation analysis shows that mem_psa3d outperforms delta_psa3d in separating permeable from impermeable compounds, that is a stronger and more defensible conclusion than a claim based on ΔPSA alone. The correlation table should include both.

**3. Harmonize the H-bond distance cutoff between Tier-1 and Tier-2.** The Tier-1 conformer engine uses H...A ≤ 3.0 Å; the Tier-2 `count_hbonds_xyz()` uses H...A < 2.5 Å. Any joint table comparing Tier-1 ΔHB to Tier-2 ΔHB will be methodologically inconsistent. Standardize on 2.5 Å (the Baker-Hubbard criterion) for both tiers before the final run, or at minimum add a footnote to every table that reports ΔHB values from both tiers simultaneously.

**4. Pre-flag the c*[PSLYF] artifact as a named limitation in the analysis notebook.** Rather than treating the 67 Å² ΔPSA result as a data quality issue, present it as a transparent methodological finding: "Tier-1 vacuum sampling assigns high ΔPSA to conformationally rigid impermeable peptides (c*[PSLYF]: 67 Å², expected ~0). This false-positive arises from ETKDGv3 sampling collapsed ring-pucker conformers without aqueous solvation. Tier-2 CREST+ALPB resolves this artifact." This framing converts a weakness into a result that strengthens the scientific story about why Tier-2 is necessary.

**5. Use the cross-lab mean PAMPA for CsA (−5.90) and report the measurement range (−6.60 to −5.01) as an error bar in all figures.** The single-lab Rezai 2006 value of −6.60 would classify CsA as impermeable, which contradicts its oral bioavailability data and the clinical record. Any figure that plots PAMPA vs ΔPSA for the reference compounds should use the mean and show the range, or the CsA data point will occupy the wrong side of the −6.0 threshold in a misleading way.

---

## Flags for Oral Presentation Defense

The following specific claims, if made in a 10-minute presentation, would require immediate qualification:

- "ΔPSA = 59.6 Å² for CsA matches Witek 2016": **Incorrect as stated.** Should be "our Tier-1 vacuum approximation yields 59.6 Å², which is directionally consistent with the 75 Å² value from Witek's explicit-solvent MD, representing 80% of the literature value."
- "ΔPSA separates compound.1 from 1NMe3": **Only at n≥100.** At n=20 the separation is 3.5 Å² which is within sampling noise. Verify this holds at n=200 before claiming it.
- "c*[PSLYF] shows zero chameleonic potential": **Tier-1 says 67 Å²; the true value is expected to be ~0.** Only claimable after Tier-2 CREST result.
- "Our 3D PSA method is consistent with the Ertl 2000 TPSA convention": **Partially correct.** Bondi heavy-atom SASA is a 3D analog, not a direct implementation of the fragment-based TPSA. State this carefully.

---

## Overall Assessment

**Scientific soundness:** Sound. The hypothesis is experimentally grounded, the pipeline correctly operationalizes the chameleonic mechanism, and the known limitations are documented with appropriate candor throughout the project documentation.

**Minimum bar for course credit:** The pipeline, as implemented and documented, will produce a defensible result. The PSA bug fix is verified, CsA and DP-172 produce plausible Tier-1 values, and the framing is scientifically literate.

**Strongest single result available right now:** CsA ΔPSA = 59.6 Å² (Tier-1, n=20) is directionally correct and within 80% of the literature benchmark. With n=200 conformers this will likely improve, and the Tier-2 CREST run will provide the methodologically superior comparison point. If CsA's CREST ΔPSA lands between 65–85 Å², the project will have a strong mechanistic validation story.

**Primary risk:** The compound.1 / 1NMe3 story at Tier-1. If n=200 does not produce a separation larger than ~10 Å² between these two compounds, the narrative claim that "N-methylation is detectable in Tier-1 ΔPSA" will be difficult to defend. In that case, lead with PSA_mem rather than ΔPSA as the discriminating feature, and rely on Tier-2 for the mechanistic N-methylation story.

---

## 6. Strategic Alignment with Consultant rNPV Analysis

### 6.1 CSO Position on Option B

I concur with the Option B recommendation, with one important qualification. The consultant correctly identifies that live CREST on an unverified WSL2 install, with a 50-hour deadline, carries asymmetric downside risk: if CREST fails at 2 PM Monday, the student arrives at Option B anyway but with 18 fewer hours for analysis and writing. The rNPV math is sound. However, the recommendation to skip CREST was conditioned on the `compute_psa_xyz` function being an ad-hoc contact-counting approximation — a scientifically disqualifying flaw. That condition has since changed (see 6.3 below), which modestly revises the calculus but does not overturn the Option B recommendation given the timeline.

### 6.2 The DB ΔPSA Cross-Check as a Strengthening Result

Yes, I agree scientifically. The consultant correctly identifies this as a genuine positive finding rather than a fallback. The comparison is stark and unambiguous: CsA DB static ΔPSA = −1 Å² versus ensemble Tier-1 ΔPSA = 60 Å². This directly demonstrates, on the gold-standard chameleonic benchmark, that single-structure 3D descriptors are blind to the conformational switch that defines the chameleonic mechanism. The same pattern holds for 1NMe3 (DB = +1 Å² vs. ensemble = 39 Å²). A table showing all five reference compounds with DB versus ensemble ΔPSA next to experimental PAMPA makes the core argument of this project more concisely than any regression result could. This is not a methodological consolation prize — it is a primary result.

### 6.3 PSA Fix Changes the Risk Calculus for Option A

The rewritten `compute_psa_xyz` using proper rdFreeSASA with Bondi radii and manual SASAClass assignment (matching `_polar_sasa()` in `conformer_engine.py`) eliminates the only scientific objection to Option A's Tier-2 output. This raises Option A's expected quality ceiling from "questionable" to "valid." However, it does not change the installation and runtime risk for CREST on WSL2, which the consultant estimated at 45% failure probability. A scientifically valid function that never runs produces the same result as a broken one. The PSA fix improves Option A's upside; it does not reduce its dominant risk. Option B remains the recommended path, but if CREST is already installed and has been smoke-tested, the rNPV gap between A and B narrows enough (~55 vs. 67) that attempting Tier-2 on the two cleanest reference compounds (CsA and 1NMe3 only, not all five) is a reasonable opportunistic addition.

### 6.4 The c*[PSLYF] SMILES Mismatch

If the SMILES used in all Tier-1 tests corresponds to a different molecule than CycPeptMPDB ID=1829, every result reported for c*[PSLYF] — including the 67 Å² ΔPSA artifact — is potentially invalid as a characterization of that compound. The immediate consequence is that the SMILES must be verified against the CycPeptMPDB entry before the overnight run. For prioritization of the Tier-2 reference set: deprioritize c*[PSLYF] until its SMILES is confirmed. The three compounds with hardcoded literature SMILES (Hexapeptide, 1NMe3, CsA) are the highest-confidence targets for any Tier-2 run. If only two CREST jobs are feasible, run CsA (the quantitative benchmark against Witek ~75 Å²) and 1NMe3 (the N-methylation story). These two compounds anchor the entire scientific narrative and their SMILES provenance is directly cited to peer-reviewed sources.

### 6.5 Final Integrated Recommendation: Next 50 Hours

1. **Immediately (30 min):** Verify the c*[PSLYF] SMILES against CycPeptMPDB ID=1829 before any production run. Correct if wrong. Run the 10-molecule early-warning diagnostic with verbose logging.
2. **Tonight:** Launch Tier-1 at n_confs=100 (not 50 — the CSO recommends the higher floor given the 3.5 Å² separation risk at n=50 and the c*[PSLYF] artifact; the extra compute costs 2–3 hours and is within the overnight window).
3. **Monday AM:** Inspect diagnostics. If psa_failed rate is below 20% and CsA ΔPSA ≥ 55 Å², proceed to correlation analysis. If CREST is confirmed installed and functional, run it on CsA and 1NMe3 only as an opportunistic Tier-2 supplement — not as a dependency.
4. **Monday:** Correlation table (PSA_mem and ΔPSA as co-primary), AUC-ROC, UMAP, reference compound DB-vs-ensemble table.
5. **Tuesday:** Notebook writeup and submission. Lead with the DB-versus-ensemble contrast as the validation section — it is the strongest result in hand regardless of whether CREST runs.

The science is defensible. Execute cleanly and submit on time.

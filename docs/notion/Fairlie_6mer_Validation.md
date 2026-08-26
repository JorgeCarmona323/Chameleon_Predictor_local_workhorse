# Structural Validation of the Chameleon Pipeline — Fairlie 6-mers (+ CsA context)

**Date:** 2026-08-26 · **Author:** Jorge · **Status:** diagnostic phase complete (CPCM-X-scored) — see **§9** for the roadmap (refinement closes out the tool; ML experiments decide whether training even needs it)

> **TL;DR.** Against deposited NMR ensembles for two Fairlie cyclic hexapeptides (7L96, 7L98) plus CsA (11-mer), the pipeline **reproduces the backbone fold** (RMSD covers every NMR model < 2 Å; Rgyr spot-on) but **systematically under-exposes polar surface** (3D-PSA low by 17–27 Å², +1 intramolecular H-bond). The error is **systematic, directional, and grows with ring size**. A controlled re-weighting test proves it lives in **conformer *generation* (implicit solvent), not scoring** — CPCM-X populations barely move it. Crucially, the bias **does not break permeability ranking** (Test A: ρ = −0.33 pooled, −0.44 within 6-mers, N = 3,256), so for *ranking* the pipeline is usable as-is and the offset is **calibratable** — the expensive explicit-solvent fix is only needed for absolute PSA / ≥9-mers.

---

## 1. Question

Does the CREST → CPCM-X → 3D-descriptor pipeline reproduce **experimentally-determined** cyclic-peptide conformations, and if not, **what is the structure of the error** and does it matter for building a permeability model?

## 2. Validation set

| compound | PDB | size | NMR solvent | why chosen |
|---|---|---|---|---|
| Fairlie **cmpd 4** | 7L96 | 6-mer | CDCl₃ | crystal-quality NMR ensemble, membrane-mimic solvent |
| Fairlie **cmpd 10** | 7L98 | 6-mer | DMSO | most permeable analog (PAMPA 11.0), thiazole scaffold |
| **CsA** | CCDC (DEKSAN) | 11-mer | — (X-ray) | large-flexible reference; GFN-FF regime |

Experimental references stored at `data/experimental_structure_references_Fairlie/` (full 4–5 model NMR ensembles).

## 3. Method

- **Stage 1 — geometry:** CREST iMTD-GC, **GFN2-xTB + ALPB** implicit solvent, per NMR/assay solvent (water / chloroform / hexane; DMSO for cmpd 10). Thousands of conformers/leg.
- **Stage 2 — energy:** **CPCM-X** single-point solvation free energy per conformer (fixed geometry, no re-opt) → ΔG_transfer + Boltzmann populations.
- **Descriptors:** 3D-PSA (Ono/Ertl definition), radius of gyration, geometric intramolecular H-bonds (IMHB), PMI shape (NPR1/NPR2).
- **Ensemble-to-ensemble comparison** (our thermodynamic ensemble vs the experimental restraint-satisfying ensemble), three views:
  - **A — descriptor distributions:** population-weighted mean ± spread vs the range over all NMR models.
  - **B — RMSD coverage:** best-match RMSD of our populated conformers to *every* deposited NMR model.
  - **D — PMI shape:** NPR1/NPR2 cloud vs the NMR models.
- **Controlled re-weighting test:** identical conformer sets, **CREST-ALPB weights vs CPCM-X weights**, isolating whether the error is in *weighting* or *geometry*.

## 4. Results

### 4.1 The fold is reproduced (B — RMSD coverage)

Our populated ensemble contains a conformer within ~1.5–2 Å of **every** deposited NMR model:

| | cmpd 4 (5 models) | cmpd 10 (4 models) | CsA |
|---|---|---|---|
| best / worst match | 1.48 – 2.12 Å | 1.77 – 1.93 Å | 1.51 Å |

→ Backbone conformation and overall size are **correct**.

### 4.2 …but polar surface is systematically under-exposed (A — descriptors)

| descriptor | cmpd 4 ours / exp | cmpd 10 ours / exp | CsA ours / exp |
|---|---|---|---|
| **3D-PSA (Å²)** | 117 / **135** (−18) | 88 / **102** (−15) | 85 / **112** (−27) |
| **Rgyr (Å)** | 4.70 / 5.01 | 4.72 / **4.73 ✓** | 5.86 / 5.92 ✓ |
| **IMHB** | 3.0 / 2.0 (+1) | 1.9 / 1.5 (+1) | 5 / 4 (+1) |

**Same signature every time:** PSA biased *low*, IMHB *+1*, Rgyr *correct*. The molecule is the right **size** but over-shields its polar groups via one spurious transannular H-bond. The bias **grows with ring size** (6-mers −15/−18 Å²; 11-mer −27 Å²).

### 4.3 The error is in the geometry, not the weighting (decisive test)

Same conformers, two weightings:

| | CREST-ALPB wt | **CPCM-X wt** | exp | corrected? |
|---|---|---|---|---|
| cmpd 4 PSA | 117.4 | **116.6** | 134.8 | ✗ (−0.8) |
| cmpd 4 IMHB | 3.0 | **3.0** | 2.0 | ✗ |
| cmpd 10 PSA | 84.8 | **87.7** | 102.3 | partial (+2.9) |
| cmpd 10 IMHB | 2.83 | **1.91** | 1.5 | ✓ |

Proper CPCM-X populations **barely touch** the PSA gap and don't pull the shape cloud onto the NMR points (below). The *only* thing re-weighting fixed was cmpd 10's H-bond count. **Conclusion: the conformers are born collapsed under implicit solvent — a better *score* on a pre-shielded *pool* stays shielded.** The error originates in **stage-1 sampling**, and no stage-2 scoring rescues it.

![PMI shape, CPCM-X-weighted](fairlie_pmi_cpcmx.png)

*Our CPCM-X-weighted ensemble (teal, point size ∝ population) vs the NMR models (orange). cmpd 4 stays parked toward the sphere corner (over-globularized); cmpd 10 keeps heavy weight on a collapsed disc state the NMR ensemble never visits.*

### 4.4 ΔG_transfer (CPCM-X)

| compound | water→membrane-mimic | water→hexane |
|---|---|---|
| cmpd 4 | −14.03 (→CHCl₃) | **−5.70** |
| cmpd 10 | −4.16 (→DMSO, *not* a permeability proxy) | **−5.96** |

Both partition favorably into the membrane-core mimic; **cmpd 10 is more membrane-favorable (−5.96 vs −5.70)**, consistent with it being the more permeable analog (PAMPA 11.0). DMSO scored cleanly under CPCM-X (it's a supported solvent), but water→DMSO is *not* a permeability signal — it exists only to weight the DMSO-ensemble descriptors that match the 7L98 structure.

### 4.5 Does the bias matter for a model? (Test A — PSA vs permeability)

Joined our **biased** 3D-PSA on 3,258 CREMP peptides to **CycPeptMPDB permeability** by canonical SMILES (N = 3,256 matched):

| descriptor (ours, biased) | Spearman ρ vs permeability |
|---|---|
| psa_mean | **−0.334** (p = 1.6e-85) |
| psa_max (most-exposed conf) | **−0.396** |
| psa_spread (ΔPSA proxy) | −0.231 |
| imhb_bb_mean | −0.255 |
| rg_mean | −0.046 (≈ 0 — not just "small = permeable") |
| — *within 6-mers* | **−0.439** |
| — *within 7-mers* | −0.233 |
| *baseline: 2D TPSA* | −0.398 |

![Test A](test_a_psa_vs_permeability.png)

**The biased PSA still ranks permeability**, strongly and significantly — the systematic offset does **not** destroy the ranking (as expected for a monotonic bias). *Honest caveat:* on this **pooled** set our 3D psa_mean (−0.33) does **not** beat cheap 2D TPSA (−0.40). The 3D pipeline's value therefore is **not** pooled PSA ranking — it's (i) **within-size / chameleon** cases, (ii) **ΔG_transfer** (a signal 2D descriptors don't have), and (iii) **robustness under distribution shift** (the earlier CREMP leave-source-out result, where 2D collapses and 3D holds). Notably **psa_max** — the most-exposed conformer's polar surface — is the strongest single 3D predictor, which is mechanistically sensible (the open-state desolvation penalty governs membrane entry).

## 5. Error model (what to trust)

**Root cause.** Continuum solvent has no discrete solute–solvent H-bonds, so nothing "pays" to expose an NH/carbonyl. GFN2+ALPB therefore over-forms internal H-bonds during sampling → over-shields polar surface → over-globularizes, **without** collapsing the radius or breaking the fold. (A secondary suspect is GFN2 itself over-stabilizing H-bonds; only an explicit-solvent test separates them — see §6.)

| quantity | trust? | why |
|---|---|---|
| backbone fold / RMSD | ✅ | reproduced < 2 Å |
| Rgyr, molecular size | ✅ | within 0.3 Å |
| **absolute 3D-PSA** | ⚠️ biased low 15–27 Å² | continuum over-shielding |
| **ΔPSA / chameleonicity** | ⚠️ compressed | dynamic range of exposure shrunk (CsA ΔPSA ~0 vs ~48) |
| IMHB count | ⚠️ +1 | same error, counted |
| **permeability *ranking*** | ✅ | survives the bias (Test A) |

**Applicability domain.** Small rings (6–7-mers): usable on **calibration**. ≥9-mers / CsA-class: implicit solvent may miss the open state entirely (a *population* failure, not just a shift) → needs explicit solvation.

## 6. How to improve (ranked) — and whether it's worth it

Scoring is ruled out (§4.3); the fix must change the **generated geometry**:

1. **Calibrate the bias (cheapest, ML-native).** Treat the size-dependent −15→−27 Å² offset as a known correction learned from this ladder. Test A says ranking already survives → **sufficient for a ranking model.**
2. **QCG micro-solvation** (`crest --qcg`, explicit water shell) — the decisive, affordable diagnostic on **one** compound. Adds the missing explicit H-bonds *and* separates "solvent model" from "GFN2 force field": if PSA opens → solvent model; if not → GFN2.
3. **Explicit-solvent MD** (MACE-OFF + OpenMM) — physically correct, GPU-scalable; the real answer if 1–2 fall short, and the *only* route to CsA-class population failures.

> **Note on the solvent model:** ALPB is already the best implicit model usable for sampling (GBSA is worse; CPCM-X has no gradients and cannot drive CREST). A "better continuum" cannot help — the missing physics is explicit H-bonds, so improvement requires leaving the continuum entirely.

## 7. Implications for model development

- **Ranking model (near-term):** proceed on the biased descriptors + calibration. The bias is systematic and rank-preserving; explicit-solvent compute is **not** required to start.
- **The 3D pipeline must justify itself** where 2D can't: ΔG_transfer, within-size chameleon discrimination, and distribution-shift robustness — *not* pooled PSA correlation.
- **Feature to prioritize:** `psa_max` / open-state exposure (strongest single 3D predictor) and ΔG_transfer.
- **Applicability-domain gate:** flag ≥9-mer / high-flexibility inputs as outside the calibrated domain until explicit-solvent sampling is validated there.

## 8. Next steps

- [ ] **`collect_fairlie.py`** — lock this table (RMSD coverage + CREST-vs-CPCM-X descriptors + ΔG) as a one-command reproducible deliverable.
- [ ] **QCG probe** on one 6-mer (§6.2) — does explicit water recover PSA? separates solvent vs GFN2.
- [ ] **Test A, extended** — within-size Spearman on a leave-source-out split; does calibrated 3D beat 2D on chameleons specifically?
- [ ] Run WhC3 (3-solvent ensemble in hand) through the same pipeline for the size ladder.

## 9. Direction & roadmap

This report closes the **diagnostic** phase. The path forward is two parallel tracks toward two large milestones.

**Track 1 — close out the pipeline as a computational tool.** Adding the **r²SCAN-3c + CPCM (CENSO) refinement** stage completes a validated, end-to-end toolkit for *solid computational exploration of cyclic peptides*: **sample (CREST/GFN2) → refine (r²SCAN-3c/CPCM) → score (CPCM-X ΔG) → 3D descriptors**, benchmarked against experimental NMR ensembles. This is the version anyone can pick up and use — with the honest label that **with refinement it is medium-throughput** (reference-tier, ~days/molecule; see §1 runtime), built for careful case studies, not screening thousands. The QCG + r²SCAN-3c experiments then let us **decompose the residual error (GFN2 vs continuum)** and lock the final **conclusions & limitations / applicability domain**.

**Track 2 — decide whether the *model* even needs refinement, then build it.** In parallel we **schedule the ML experiments** (extended Test A: does calibrated cheap-GFN2 descriptor ranking suffice, within-size / leave-source-out?). If the model doesn't need refined descriptors — likely, given Test A already survives the bias — then **production training stays on the cheap, high-throughput GFN2 pipeline** and refinement remains only the *validation anchor*. That clears the way to **focus on actually developing the model** (the layered predict → explain → design architecture).

> **The through-line:** refinement *closes out the exploration tool* (medium-throughput, publishable); the ML experiments tell us whether we can *keep training cheap*; then we *build the model*. These are big milestones, but that is the intended direction.

---

*Reproduce:* ensembles at `results/conformers/fairlie_6mer_cmpd{4,10}/`, energies at `results/free_energy/fairlie_6mer_cmpd{4,10}/`, figures at `results/validation/`.

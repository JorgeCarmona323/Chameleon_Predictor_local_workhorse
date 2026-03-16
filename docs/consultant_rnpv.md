# Risk-Adjusted NPV Analysis — CHEM 269 Final Project
## Strategic Decision Brief

**Prepared:** Sunday, March 15, 2026, ~9 PM
**Deadline:** Tuesday, March 17, 2026, 11:59 PM
**Time remaining:** ~50 hours
**Analyst:** Strategic Risk Assessment (Claude Sonnet 4.6)

---

## Situation Summary

The pipeline is partially operational. The critical `_polar_sasa` bug is fixed and produces physically plausible ΔPSA values (CsA: 59.6 Å², directionally correct vs. 75 Å² MD reference). The conformer engine has not run at scale. Tier-2 CREST has dry-run-passed but not executed for real. Two known failure modes remain: (1) c\*[PSLYF] shows ΔPSA=67 Å² at n=20 when the expected value is ~0 — a vacuum sampling artefact; (2) ΔHB=0 for all five test molecules at n=20, likely a sampling depth issue. The only clean results in existence are five rows from the n=20 smoke test. All prior CSV outputs must be discarded.

The scientific framing is strong: the chameleonic mechanism is NMR-proven, the reference compounds are well-chosen, and the feature set is physically motivated. The project has a coherent narrative from observation to hypothesis to answer. The science is defensible. The execution risk is the question.

---

## Framework Definition

Since the currency is not money, define:

- **Value (V)** = P(submit quality result) × Quality Multiplier (QM), where QM: 1.0 = minimal viable, 1.5 = solid, 2.0 = strong, 3.0 = exceptional
- **Cost (C)** = Expected hours × P(wasted, run fails or produces uninterpretable output)
- **rNPV** = (V × 100) − (C × 10), normalized to 0–100 scale
- **Expected Value (EV)** = P(success) × QM × 100, the gross value before cost discount

**Grade mapping:** A = EV ≥ 65 after cost; A- = EV 55–64; B+ = EV 45–54; B = EV 35–44

---

## Option Analysis

---

### Option A — Current Path (Tier-1 overnight n=50, Tier-2 CREST Monday, analysis Tuesday)

**Timeline:**
- Sunday night: Launch Tier-1 at n\_confs=50, ~7,298 compounds, 3 CPUs → estimated 4–10 hours (wide uncertainty)
- Monday AM: Inspect results, verify ΔPSA distributions are sensible, fix any issues
- Monday afternoon: Run Tier-2 CREST on 5 compounds × 2 solvents → 2–6 hours
- Monday evening: Correlation analysis, AUC-ROC, UMAP, figures
- Tuesday: Notebook writeup, polish, submit

**Estimated total active hours:** 18–30 hours
**Estimated elapsed real time needed:** 42–50 hours (tight against the 50-hour window)

**Success probability decomposition:**

| Risk event | P(occurs) | Recovery available? |
|---|---|---|
| Tier-1 run completes on time (≤8h) | 0.70 | Use n\_confs=20 fallback |
| Tier-1 ΔPSA distributions are physically interpretable | 0.75 | Filter c\*[PSLYF]-class outliers |
| c\*[PSLYF] vacuum artefact persists at n=50 | 0.65 | Flag as known artefact, exclude or annotate |
| ΔHB signal nonzero in ≥1 compound at n=50 | 0.55 | Deemphasize ΔHB, lead with ΔPSA |
| CREST installs and runs cleanly on WSL2 | 0.55 | Skip Tier-2 (drops to Option B) |
| CREST produces interpretable ΔPSA for 3/5 refs | 0.50 | Use Tier-1 only |
| Full analysis completes Monday evening | 0.60 | Compress to Tuesday morning |
| Notebook polished and submitted by 11:59 PM | 0.80 | Submit partial notebook |

**Compound success probability (all major steps succeed):**
P(A works end-to-end) ≈ 0.70 × 0.75 × 0.55 × 0.80 ≈ **0.23**

**But partial success is also valuable.** If CREST fails and you fall back to Tier-1 only (Option B territory), P(that path succeeds) is much higher. The real question is: what is the expected quality across all outcomes?

| Outcome | P | QM | Contribution |
|---|---|---|---|
| Full pipeline works, strong results | 0.20 | 2.5 | 50 |
| Tier-1 works, CREST fails, good Tier-1 results | 0.35 | 1.8 | 63 |
| Tier-1 works, results weak/noisy, honest write-up | 0.25 | 1.2 | 30 |
| Tier-1 fails or produces garbage | 0.20 | 0.5 | 10 |

**EV(A) = 0.20×250 + 0.35×180 + 0.25×120 + 0.20×50 = 50 + 63 + 30 + 10 = 153 raw units**
**Normalized EV = 153/3.0 = 51 (out of 100)**

**Cost:** High. If Tier-1 takes 10 hours and produces garbage (P=0.20), that is 10 hours lost overnight with no recovery time. If CREST fails after 6 hours on Monday, another day is lost. Worst-case total wasted hours: 16. P(≥8h wasted) ≈ 0.25.

**Risk-adjusted cost penalty:** 16 × 0.25 × 10 = 40 points → net rNPV = 51 − (40/10) = **47**

**rNPV(A) ≈ 47 / 100**

**Key risk:** The pipeline has never run at scale. At n=20 on 5 molecules, it ran cleanly. Whether it runs cleanly on 7,298 molecules for 4–10 hours without memory issues, crashes, or silent errors (wrong CSV schema, SMILES parse failures causing mass `psa_failed`) is unknown. This is the single largest binary risk event in the entire decision tree.

---

### Option B — Tier-1 Only (Skip CREST, use DB 3DPSA as Tier-2 proxy)

**What this means:** Run Tier-1 at n\_confs=50 overnight. Monday: run correlation analysis and UMAP using Tier-1 ΔPSA features. For Tier-2 "validation," use the pre-computed 3DPSA values already in CycPeptMPDB (H₂O PSA and CHCl₃ PSA columns) as a single-structure proxy — explicitly acknowledged as a static approximation rather than an ensemble method. Reference compound section compares Tier-1 ΔPSA vs. DB static ΔPSA and explains why the ensemble method recovers the CsA signal that the static method misses (CsA: DB ΔPSA = −1 Å², Tier-1 ΔPSA = 60 Å²; this contrast is itself a result).

**Timeline:**
- Sunday night: Launch Tier-1 at n\_confs=50 → 4–10 hours
- Monday AM: Inspect results, validate reference compounds against DB ΔPSA
- Monday AM–PM: Build correlation table, AUC-ROC, UMAP, figures (3–6 hours)
- Monday PM: Notebook writeup (4–6 hours)
- Tuesday: Polish and submit

**Estimated total active hours:** 14–22 hours. Comfortable margin.

**Success probability decomposition:**

| Risk event | P(occurs) | Consequence |
|---|---|---|
| Tier-1 run completes in time (≤8h) | 0.70 | If longer, use n\_confs=20 results |
| Tier-1 ΔPSA signal is meaningful (CsA/1NMe3 > Hexapeptide) | 0.80 | More conformers help |
| Correlation analysis shows 3D ≥ 2D at AUC-ROC | 0.55 | Negative result still publishable |
| UMAP separates permeable/impermeable with 3D features | 0.60 | — |
| Notebook completed and submitted | 0.90 | Strong time margin helps |

**P(B works end-to-end)** ≈ 0.70 × 0.80 × 0.90 ≈ **0.50**

| Outcome | P | QM | Contribution |
|---|---|---|---|
| Tier-1 works, strong 3D > 2D result, well-written | 0.30 | 2.0 | 60 |
| Tier-1 works, modest signal, honest analysis | 0.35 | 1.5 | 52.5 |
| Tier-1 works, weak/null signal, honest write-up | 0.20 | 1.0 | 20 |
| Tier-1 fails completely | 0.15 | 0.4 | 6 |

**EV(B) = 0.30×200 + 0.35×150 + 0.20×100 + 0.15×40 = 60 + 52.5 + 20 + 6 = 138.5 raw units**
**Normalized EV = 138.5/2.0 = 69 (out of 100)**

**Cost:** Moderate. Tier-1 only. If it fails, you lose one overnight slot (4–10 hours) and pivot to Option C with 36+ hours remaining. Risk is manageable.

**Risk-adjusted cost penalty:** P(Tier-1 fails and Option C pivot required) = 0.15 × 10h × 10 / 10 = 1.5 points.

**rNPV(B) ≈ 69 − 2 = 67 / 100**

**Scientific defensibility note:** Explicitly acknowledging that DB ΔPSA is a single-structure proxy and contrasting it with ensemble ΔPSA (where CsA goes from −1 Å² to 60 Å²) is itself a strong result — it directly demonstrates that static 3D values miss the chameleonic signal. This makes the "Tier-2 as DB cross-check" framing a legitimate scientific contribution, not a compromise.

---

### Option C — 2D + DB 3DPSA Only (No conformer generation)

**What this means:** Do not run conformer_engine.py at all. Use: (1) 2D descriptors from CycPeptMPDB (MW, cLogP, TPSA, HBD, RotBonds); (2) the pre-computed DB H₂O-PSA and CHCl₃-PSA static values already in the database; (3) compute DB ΔPSA = H₂O-PSA − CHCl₃-PSA as the sole 3D-proxy feature. Run correlation analysis and UMAP immediately. Submit a notebook that frames this as "static 3D descriptors vs. 2D" with honest limitations.

**Timeline:**
- Sunday night / Monday AM: Build feature matrix from existing CSV, run analysis (2–4 hours)
- Monday: Writeup and notebook (4–8 hours)
- Tuesday: Polish and submit

**Estimated total active hours:** 8–14 hours. Very comfortable margin.

**Success probability:** Very high for completion. P(submit something) ≈ 0.95. But scientific quality ceiling is low.

| Outcome | P | QM | Contribution |
|---|---|---|---|
| DB ΔPSA shows meaningful signal over 2D | 0.40 | 1.3 | 52 |
| DB ΔPSA weak, honest analysis + limitations | 0.45 | 0.9 | 40.5 |
| Analysis runs, no interpretable findings | 0.15 | 0.6 | 9 |

**EV(C) = 0.40×130 + 0.45×90 + 0.15×60 = 52 + 40.5 + 9 = 101.5 raw units**
**Normalized EV = 101.5/1.3 = 78 (raw), but quality ceiling is 1.3 = B-range**

**The fundamental problem with Option C:** The proposal explicitly motivates ensemble 3D descriptors over static structures. The scientific framing document notes that DB ΔPSA = −1 Å² for CsA, despite the compound having a known ~75 Å² conformational switch. Using DB ΔPSA as the sole "3D" feature directly contradicts the project's core motivation. A grader who reads the proposal and the notebook will notice this gap. The ceiling is B+, not A.

**rNPV(C) ≈ 42 / 100** (high completion certainty, low quality ceiling, proposal-to-delivery mismatch)

---

### Option D — Hybrid: Tier-1 at n=20, supplement with DB 3DPSA

**What this means:** Use the n=20 smoke test results already in hand (partially — need to re-run on the full dataset but at n=20 instead of n=50, cutting runtime roughly in half). Accept weaker ΔPSA signal (1NMe3 vs. Hexapeptide gap is only 3.5 Å² — marginal). Supplement with DB ΔPSA as a validation column. Analyze Monday AM after a 2–4 hour run.

**Timeline:**
- Sunday night: Re-run Tier-1 at n\_confs=20, full dataset → ~2–4 hours
- Monday AM: Results in hand, begin analysis
- Monday: Correlation, UMAP, figures, notebook (8–12 hours)
- Tuesday: Polish and submit

**Key problem with n=20:** The 1NMe3 vs. Hexapeptide ΔPSA gap is 3.5 Å² (39.3 vs. 35.8 Å²). This is the canonical demonstration that the chameleonic mechanism is detectable. At 3.5 Å², it is within the noise of a 20-conformer vacuum ensemble. The c\*[PSLYF] artefact (67 Å² when expected ~0) is also unresolved at n=20. The entire discriminatory power of the feature depends on the signal being real. At n=20, it is questionable.

**Success probability:**

| Outcome | P | QM | Contribution |
|---|---|---|---|
| n=20 provides enough signal, clean results | 0.30 | 1.6 | 48 |
| n=20 signal marginal, honest analysis | 0.40 | 1.1 | 44 |
| n=20 produces noise/artefacts, weak paper | 0.30 | 0.7 | 21 |

**EV(D) = 0.30×160 + 0.40×110 + 0.30×70 = 48 + 44 + 21 = 113 raw units**
**Normalized EV = 113/1.6 = 71 (raw)**

**However:** Option D gives up the main methodological advantage (sufficient conformer sampling) without eliminating the main risk (pipeline failure on the full dataset). You still run the pipeline at scale, you just get weaker results. It is a dominated strategy relative to Option B (same risk, lower reward at n=50). The only case where D beats B is if n=50 takes longer than 8 hours and timeline collapses — and in that case, D's n=20 run still takes 2–4 hours, providing a partial result.

**rNPV(D) ≈ 52 / 100**

---

## Scoring Summary Table

| Option | Description | Est. Hours | P(Success) | Quality Ceiling | EV (raw) | Risk Penalty | rNPV |
|---|---|---|---|---|---|---|---|
| **A** | Full pipeline: Tier-1 n=50 + CREST Tier-2 | 18–30h | 0.50 end-to-end | Exceptional (A) | 51 | −4 | **47** |
| **B** | Tier-1 n=50 + DB ΔPSA as Tier-2 proxy | 14–22h | 0.65 end-to-end | Strong (A-/A) | 69 | −2 | **67** |
| **C** | 2D + DB static only, no conformers | 8–14h | 0.90 completion | B+/B | 42 | −0.5 | **42** |
| **D** | Tier-1 n=20 + DB supplement | 10–16h | 0.55 end-to-end | B+/A- | 52 | −3 | **49** |

**Winner: Option B (rNPV = 67)**

---

## RECOMMENDATION

**Execute Option B: Run Tier-1 at n\_confs=50 tonight, skip live CREST, use DB static ΔPSA as Tier-2 cross-check.**

This is not a compromise. It is the correct strategic choice given the risk-time profile. Here is the reasoning:

### Why not Option A (full pipeline)?

CREST on WSL2 is an unverified dependency with a 45% chance of failing. If CREST fails at 2 PM Monday after you have spent Sunday night and Monday morning on Tier-1, you have ~32 hours left and no Tier-2 results. You are forced into Option B anyway — but now with less time for analysis and writing. Option A's upside (rNPV=47 at best) does not justify the downside of arriving at Option B on Tuesday morning instead of Monday morning.

Additionally, the `compute_psa_xyz` function in `tier2_crest.py` uses an ad-hoc contact-counting approximation (not FreeSASA) that will produce systematically biased PSA values even if CREST runs. The Tier-2 output is scientifically questionable without fixing that function first. Attempting to fix it and run CREST in the same overnight window is too much scope.

**CREST is not worth the risk this week.** It can be a future direction.

### Why Option B is scientifically defensible (not just strategically safe)

The DB ΔPSA cross-check is not a fallback — it is a positive finding:

- CsA DB ΔPSA = −1 Å². CsA Tier-1 ensemble ΔPSA = 60 Å². The static value fails completely; the ensemble value recovers the known ~75 Å² chameleonic switch. This comparison IS the validation.
- c\*[PSLYF] DB ΔPSA = 0 Å². If Tier-1 shows 67 Å² (artefact), the contrast with the DB value isolates the vacuum sampling limitation. That is a methodologically honest, scientifically interesting finding: "vacuum ETKDGv3 overestimates chameleonicity for rigid impermeable peptides; solvated sampling would correct this."
- 1NMe3 DB ΔPSA = +1 Å². Tier-1 ensemble = 39 Å². Again, ensemble recovers the chameleonic signal that the static structure completely misses.

You have a table that writes itself: **static DB ΔPSA vs. ensemble Tier-1 ΔPSA for 5 reference compounds, where the ensemble method agrees with experimental NMR/MD evidence and the static method does not.** This directly answers the project's core question (why 2D/static descriptors fail) and validates the pipeline approach. It is Option B masquerading as Option A's validation story.

### The execution plan for Option B

**Sunday night (now):**

1. Set `n_confs=50`, 3 CPUs. Launch Tier-1 on the full CycPeptMPDB dataset. Pipe stdout to a log file. Set a 10-hour timeout alarm for Monday morning.
2. Before launching, verify: (a) delete/archive the old `conformer_descriptors_raw.csv`; (b) confirm the `_polar_sasa` fix is in place; (c) run a 3-molecule smoke test at n=50 to confirm no new errors before committing to 7,298 molecules.

**Monday morning (8–10 AM):**

3. Inspect results. Key diagnostic checks: (a) CsA ΔPSA should be 60–80 Å²; (b) 1NMe3 ΔPSA > Hexapeptide ΔPSA by >5 Å²; (c) c\*[PSLYF] ΔPSA should ideally be lower at n=50 but may still be inflated — document this explicitly; (d) distribution of ΔPSA across all 7,298 molecules should be roughly log-normal, not bimodal or zero-inflated.
4. If results pass diagnostic: proceed to correlation analysis, AUC-ROC, UMAP.
5. If results fail (mass errors, zero PSA, garbage distribution): pivot to Option C immediately. With 36+ hours remaining, Option C can still produce a solid B+ submission.

**Monday (12 PM – 8 PM):**

6. Build correlation table: 3D features vs. 2D features vs. PAMPA (Spearman ρ). Confirm PSA_mem and ΔPSA outperform TPSA and HBD.
7. Run AUC-ROC: permeable vs. impermeable (threshold −6.0 log cm/s). Target: AUC\_3D > AUC\_2D.
8. UMAP: first with 2D features only (expected: weak cluster separation), then with 3D features. Reference compound overlay.
9. Reference compound table: DB static ΔPSA vs. Tier-1 ensemble ΔPSA vs. PAMPA. This is your Tier-2 substitute and it tells a compelling story.

**Monday evening – Tuesday:**

10. Notebook writeup. Key sections: background (chameleonic mechanism), methods (ETKDGv3 + MMFF94s, Bondi radii PSA), results (correlation table, AUC-ROC, UMAP, reference compound validation table), discussion (limitations: vacuum sampling, ΔHB signal weakness, n=50 vs. explicit solvent), conclusions.
11. Submit.

---

## Quick Wins That Dramatically Improve Option A's Expected Value

Since Option B is recommended, these quick wins apply equally — they raise Option B's ceiling to near-Option-A quality:

**QW1 (high impact, 30 min): Fix `compute_psa_xyz` in `tier2_crest.py` before running CREST.**
If you decide to attempt CREST anyway, the contact-counting PSA approximation must be replaced with FreeSASA before results are scientifically valid. Without this fix, CREST results are unusable regardless of whether CREST runs successfully.

**QW2 (high impact, 15 min): Add molecule-size-dependent n\_confs.**
DP-172 is an 11-residue macrolide (~1500 Da). With `pruneRmsThresh=0.5` and n=50, it may still be severely undersampled. A simple rule (n\_confs = max(50, MW/20)) would allocate ~60–75 conformers to large macrolides automatically. This costs almost no time and meaningfully improves the reference compound benchmarks.

**QW3 (high impact, 20 min): Audit and report c\*[PSLYF] separately.**
Rather than hoping the n=50 artefact resolves, explicitly flag in the pipeline output that molecules with ≥7 HBD and DB ΔPSA=0 may show inflated Tier-1 ΔPSA due to vacuum sampling. Compute a "vacuum artefact risk score" (HBD count + rigidity proxy) and report it alongside ΔPSA. This turns a liability into a finding.

**QW4 (medium impact, 10 min): Relax HB distance cutoff to 3.5 Å H...A.**
The current 3.0 Å cutoff may be missing borderline H-bonds that are real. Relaxing to 3.5 Å at n=50 might recover the ΔHB signal for CsA and DP-172. Low-cost, high-upside change if ΔHB=0 persists at n=50.

**QW5 (medium impact, done already): The PSA bug fix + reference compound table.**
You already have the most important quick win in hand: the `_polar_sasa` fix with Bondi radii and manual SASAClass. The n=20 results (CsA=60 Å², 1NMe3=39 Å²) are scientifically plausible and ready to be included in the reference compound validation table regardless of what the full-scale run shows. These five data points are a floor — you will not go below them.

---

## The Single Biggest Risk to Option B and Its Contingency

**Biggest risk: Tier-1 full-scale run fails silently — produces a CSV where 40–60% of molecules show `psa_failed`, leaving too few valid rows for a meaningful statistical analysis.**

This is not hypothetical. The previous pipeline run produced only 3 rows before the PSA bug was discovered. The fix is verified on 5 molecules at n=20, not on 7,298 molecules. Edge cases that could cause mass failures include: SMILES parse failures for unusual macrolide scaffolds; MMFF94s force field failures on molecules with unrecognized atom types (common in N-methylated peptides with nonstandard connectivity); memory errors from large macrolides (DP-172-class compounds can have 100+ heavy atoms); or the `Uncharger` step producing impossible valence states.

**Contingency protocol:**

1. Check `psa_failed` rate after the first 100 molecules (set a checkpoint log). If >30% fail in the first 100, kill the run and diagnose before proceeding to 7,298.
2. If the PAMPA subset alone is smaller (check how many rows in CycPeptMPDB have a PAMPA value — likely 1,000–2,000), prioritize running the PAMPA subset first. A 1,500-molecule run at n=50 takes roughly 1–2 hours and gives you the entire analysis-ready dataset even if the full 7,298-compound run fails.
3. If failure rate is high but the PAMPA subset (your analysis target) succeeds, ignore the full-dataset failure entirely. The correlation analysis only needs the PAMPA subset. The UMAP can be run on the PAMPA subset. The science does not require all 7,298 compounds.
4. Nuclear contingency (Option C pivot): if 60%+ of the PAMPA subset fails, stop all conformer work at 8 AM Monday, pivot to Option C with 40 hours remaining. Option C will produce a B+/A- submission. This is acceptable.

**Early warning system:** Before the overnight run, run the first 10 molecules with verbose logging. If all 10 pass in under 5 minutes, scale confidence is high. If any of the first 10 fail with new error types, investigate before committing to the overnight run.

---

## Is the Chameleonic Story Scientifically Compelling at Option D (n=20)?

**Yes, barely — but only if you are rigorous about what you claim.**

At n=20, the reference compound results are:
- CsA ΔPSA = 59.6 Å² (literature: ~75 Å²) — borderline pass, directionally correct
- 1NMe3 ΔPSA = 39.3 Å² > Hexapeptide 35.8 Å² — correct ordering, gap only 3.5 Å²
- c\*[PSLYF] ΔPSA = 67.3 Å² (expected: ~0 Å²) — artefact, fails

The fundamental problem is not the small gap between 1NMe3 and Hexapeptide (3.5 Å²) — that could be real physics at n=20. The problem is c\*[PSLYF] at 67 Å². If your negative control shows a ΔPSA larger than your positive controls (CsA=60, 1NMe3=39, c\*[PSLYF]=67), the feature has no discriminatory validity. A grader who reads the reference compound table will see that the compound with PAMPA=−9.1 and zero chameleonic potential shows the largest ΔPSA in the validation set. That is a broken compass, not a signal.

**However:** The chameleonic story is scientifically compelling as a narrative even if the quantitative results are imperfect, because:

1. The mechanism is NMR-proven and well-cited. The story does not depend on your numbers being perfect.
2. The reference compound table with DB ΔPSA vs. ensemble ΔPSA makes the point about static vs. ensemble methods regardless of the n=20 artefact.
3. If you explicitly diagnose the c\*[PSLYF] failure as "vacuum ETKDGv3 cannot correctly sample rigid impermeable peptides — this is a known limitation of the Tier-1 heuristic and is why Tier-2 solvated sampling is necessary," you turn the artefact into a discussion point rather than a refutation.

**The chameleonic story survives n=20 IF and only IF you use it as a mechanistic framing, not a quantitative claim.** Lead with the mechanism (NMR, MD, literature), use the n=20 numbers as directional evidence for CsA and 1NMe3, flag the c\*[PSLYF] artefact honestly as a pipeline limitation, and focus the quantitative analysis on the full-dataset correlation (which has 1,000+ compounds and can produce statistically meaningful AUC-ROC).

At n=50 (Option B), the story gets stronger because: (a) the 1NMe3/Hexapeptide gap should widen; (b) c\*[PSLYF] may partially resolve (more extended conformers sampled); (c) ΔHB may become nonzero for at least CsA. The incremental scientific value of n=50 over n=20 is worth the 2–4 extra hours of compute, which is why Option B dominates Option D.

---

## Final Decision Matrix

```
Current time: Sunday 9 PM, 50 hours to deadline

IF (willing to accept 35% chance of major pivot Monday AM):
    Run Tier-1 at n_confs=50, skip CREST → Option B

    Expected outcome: A- to A range, scientifically defensible
    Worst case: pivot to Option C Monday 8 AM, still submit B+

ELIF (want guaranteed submission with minimum risk):
    Use Option C (2D + DB static), start writing tonight
    Expected outcome: B+ to A- range

DO NOT:
    Run CREST without first fixing compute_psa_xyz (results unusable)
    Run n_confs=20 at scale without fixing c*[PSLYF] artefact (broken compass)
    Attempt CREST + Tier-1 simultaneously overnight (overcommitted)
```

---

## Bottom Line

**Run Option B.** Launch Tier-1 at n\_confs=50 tonight with early-warning diagnostics on the first 10 molecules. Do not touch CREST this week. Use the DB static ΔPSA as your Tier-2 reference comparison — it is genuinely informative and directly supports the project's core argument that static methods miss ensemble chameleonicity. If Tier-1 fails, Option C is a dignified fallback with 40 hours of cushion.

The science is strong enough to carry a B+ submission even with imperfect results. Your job between now and Tuesday night is to execute cleanly, report honestly, and not run out of time trying to make Option A work.

Do not let perfect be the enemy of submitted.

---

*rNPV analysis prepared March 15, 2026. All probability estimates are qualitative judgments based on the described pipeline state and should not be interpreted as precise frequentist probabilities.*

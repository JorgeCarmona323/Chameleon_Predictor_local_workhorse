# Pipeline Validation Status & Handoff — 2026-07-10

**Jorge Carmona (Hu Lab, SDSU)** · session handoff for a fresh agent to run a proper, publication-grade validation with the data already in hand.

> **Goal:** show that the conformer pipeline produces *physically faithful* ensembles, anchored to experimental NMR, and map **where that holds vs. breaks by macrocycle size** (the "validity envelope"). Publication likely needs only **coupling RMSD (+ maybe RMSF) + a visual overlay** — not a large new campaign.

---

## Where we are

| System | Size | Result | Evidence |
|---|---|---|---|
| **HexPep** (= Rezai **compound 1**) | 6-mer | **PASS** | ³J(HN-Hα) **RMSD 0.71 Hz** (CDCl₃, uniform-mean) ≈ Karplus floor 0.73 Hz; discriminates its diastereomer 3.4× (0.71 vs 2.40 Hz vs compound 9); energy-independent (Boltzmann 0.69 ≈ uniform 0.71) |
| **CsA** | 11-mer chameleon | **FAIL** | 0% cis MeVal11–MeBmt1 (A1 fold never sampled), over-collapsed water PSA — implicit-solvent limit; needs explicit water |
| CREMP ML benchmark | 6–7-mers | context | CREST 3D descriptors most robust under leave-source-out CV (fingerprints collapse) |

**Interpretation:** the cheap implicit pipeline is quantitatively trustworthy at the small/rigid end (6-mer, at the Karplus floor) and breaks for large chameleons (11-mer). That bracketing is the validation story.

## Locked validation method

- **Aggregation = uniform mean** over the ensemble (energy-free — consistent with the pipeline-wide convention; do NOT Boltzmann-weight until energies are trusted, and even then only pipeline-wide). Boltzmann is a secondary sensitivity check only.
- **Metric = RMSD** (matches the Vuister-Bax Karplus floor; report MAE only as secondary).
- **Karplus:** Vuister & Bax, *JACS* 1993, 115, 7772 — ³J = 6.51 cos²θ − 1.76 cosθ + 1.60, θ = H–N–Cα–Hα dihedral (measured directly). **Floor = 0.73 Hz RMSD** (their Karplus-vs-ubiquitin-crystal fit — the irreducible residual).
- **Solvent-matched** (Rezai NMR = CDCl₃ → use the `mem` ensemble). For chameleons, validate in both phases.
- **Diastereomer discrimination** as a second axis (match own > alternate).

## Proper-validation tasks (for the new agent)

### 1. Per-residue Leu assignment (upgrade from sorted-set)
Currently the 4 Leu couplings are compared as a *sorted set* — lenient. Instead, order the NHs by ring position (Rezai numbering) and pair each with its named experimental value:
- Anchor at **Tyr6** (unique aromatic side chain); its two backbone neighbors are **Leu1** and **Pro5**. From Leu1, traverse *away from Tyr6*: Leu1→Leu2→Leu3→Leu4→Pro5.
- Validated ring-traversal + identity code: `scripts/rezai_hexpep_identity.py` (`rezai_pattern`/`bb_adj`) — also confirms HexPep = compound 1 (canonical-SMILES match; all 9 diastereomers build→CIP self-checked).
- **Experimental (compound 1, CDCl₃, ±2 Hz):** Leu1 = 4.0, Leu2 = 10.2, Leu3 = 9.0, Leu4 = 7.2, Tyr6 = 8.4 Hz (SI Table 2). Compute RMSD with fixed pairing (expect ≥ the 0.71 sorted-set number — the honest value).

### 2. NOE / D2O-exchange as an independent second observable
- **NOEs:** counts only in the SI (compound 1: 23 intra / 12 sequential / 4 interresidue) — **no assigned peak list**, so not quantitatively usable. Skip or use qualitatively.
- **D2O amide exchange (usable — it's in the Rezai *main-paper* prose):** for compound 1, **only D-Leu2 NH exchanges** over 19 h; **Leu1, Leu3, Leu4, Tyr6 are protected/H-bonded**. Validation: compute per-NH intramolecular-H-bond occupancy (uniform-mean over the ensemble) and check the ensemble H-bonds predominantly on Leu1/Leu3/Leu4/Tyr6 with D-Leu2 exposed. A clean, discrete H-bond-pattern check independent of the couplings.

### 3. Publication deliverables
- **Coupling RMSD** (per-residue) vs the 0.73 Hz Karplus floor — the headline number.
- **RMSF** (per-atom ensemble fluctuation, energy-free) to quantify conformational spread; optionally relate to the NMR CYANA-family spread.
- **Visual overlay:** for CsA we have the A1 crystal (CCDC 2149649) to overlay against the CREST ensemble. For HexPep there is *no deposited structure* (Rezai's 20-member CYANA family coordinates are not in the PDF), so the overlay would be the CREST representative/ensemble illustrating the turn + H-bond pattern consistent with the exchange data — or a predicted-vs-experimental ³J parity plot. Decide per figure.

## Out of scope (do NOT do)
- **Generating compound 9** — it's the *worse-permeable* diastereomer and adds nothing we need; HexPep (= compound 1) already anchors the 6-mer point.
- **GFN-FF vs GFN2-xTB comparison** — Jorge is handling this separately.

## Pitfalls corrected this session (don't reintroduce)
- HexPep is **compound 1**, not compound 9 (SMILES-confirmed).
- Coupling sets: **compound 1 = Leu 4.0/10.2/9.0/7.2, Tyr 8.4 (SI Table 2)**; compound 9 = Leu 9.0/7.8/8.4/6.6, Tyr 12.0 (SI Table 3).
- Use **RMSD**, not MAE. Floor is **0.73 Hz** (not "~1").
- The `reference_compounds` entry is **correct** (not a bug): compound 1, PAMPA −6.20, `permeable: False` is right because −6.20 < the −6.0 threshold, even though it's Rezai's most-permeable diastereomer.

## Files, data, literature
- **Scripts:** `validate_hexpep_nmr.py` (couplings; uniform-mean + RMSD), `build_ensemble_from_crest.py` (raw crest xyz → sdf/json), `rezai_hexpep_identity.py` (identity + ring order).
- **Ensembles:** `results/conformers/HexPep/{aq,mem}/full_ensemble.xyz` (mem = CDCl₃, matched).
- **Outputs:** `results/2026-07-08_hexpep_nmr_validation{,_chcl3}.txt`.
- **Reports:** `docs/experiments/2026-07-08_hexpep_nmr_validation.md`, `2026-07-07_csa_v2_validation.md`, `2026-07-08_nmr_validation_design_and_ff_pivot.md`.
- **Literature (`docs/literature/`):** Rezai main + SI (HexPep), Vuister-Bax (Karplus), Danelius 2020 + the bRo5 sampling paper (`NMR_Validation_Literature/`), Pegasus + SI (`ML literature/`).

## Extending the envelope later (candidate targets, need their SIs)
- **Danelius 2020** — roxithromycin, telithromycin, spiramycin, rifampicin; **dual-solvent (water + CHCl₃) NAMFIS** ensembles (NOE + J + temperature coefficients). Mid–large ring sizes.
- **bRo5 sampling paper** — 10 drugs; NMR/NAMFIS solution conformers + X-ray → **conformer-recovery RMSD** template (does our pool contain the NMR conformers, per solvent).

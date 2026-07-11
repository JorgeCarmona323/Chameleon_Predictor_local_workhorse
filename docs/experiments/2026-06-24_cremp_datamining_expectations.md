# CREMP Datamining — Expectations & Hypotheses (2026-06-24)

**Purpose.** Record the priors and hypotheses to test when we datamine CREMP (Grambow 2024:
3,258 permeability-annotated macrocyclic peptides, generated with our exact CREST / GFN2-xTB
protocol). These come from the in-house **N = 2 R/S case study** (DOPC 3-12-8-12 and 3-12-10-12,
water vs chloroform ensembles). They are **hypotheses to test at scale, not conclusions** — every
"always" below means "always in our 4 species." This is Track B of
[descriptor path forward](2026-06-23_descriptor_review_path_forward.md).

## Source data (the 4 species, Δ = water − chloroform; + = drops in chloroform)

| | 3D-PSA Δ | Rgyr Δ | IMHB(total) Δ | mem-PSA (abs) | mem-Rg | NPR1, NPR2 (mem) |
|---|---|---|---|---|---|---|
| 8-12 R  | +28.9 | −0.22 | −2.52 | 198 | 4.66 | 0.58, 0.84 |
| 8-12 S  | +11.4 | −0.08 | −0.14 | **154** | **4.59** | 0.50, 0.91 |
| 10-12 R | +42.6 | −0.36 | −2.22 | 201 | 4.65 | 0.56, 0.87 |
| 10-12 S | +43.8 | −0.21 | −0.92 | **154** | **4.55** | 0.52, 0.96 |

S = predicted more permeable in both pairs (lower absolute mem-PSA, lower mem-Rg).

---

## 1. Supported in our set — carry as priors, confirm at scale

1. **3D-PSA drops in the apolar/membrane ensemble (4/4).** Primary, most reliable signal.
2. **Total IMHB rises in the apolar ensemble (4/4).** Mechanism behind the PSA burial (less
   solvent competition → more internal H-bonds). *Total* is the robust signal; the
   backbone-vs-sidechain split is **not** consistent (e.g. 10-12 R shifts backbone→sidechain).
3. **Absolute apolar PSA discriminates permeability better than ΔPSA.** The permeable epimer
   *reaches* the lowest apolar PSA (~154) even when it buried **less** than its partner
   (8-12 S buries only 11.4 but ends at 154; R buries 28.9 but ends at 198). ΔPSA is confounded
   by the water starting point → **lead with absolute apolar PSA; treat ΔPSA as secondary.**
4. **Rgyr rises into a tight apolar band (~4.55–4.66 Å); permeability does NOT require compaction
   on transfer.** Within the band the permeable epimer sits at the low edge → it's an
   **absolute Rg window + relative position**, not "smaller is better."

## 2. Speculative — test, do not assume (N = 2 and/or weak descriptors)

5. **Optimal Rg window, size-dependent (not "smaller = better").** Needs size-diverse data to
   define edges. See §3.
6. **"Pre-organized beats chameleonic."** The already-shielded epimer (low water PSA, high IMHB,
   small cross-solvent Δ's) outperformed the bigger chameleon (large ΔPSA/ΔRg/amphi). Caution
   against using chameleonicity *magnitude* as the permeability proxy. Could be scaffold-specific.
7. **Permeable epimer resists flattening — stays near the rod–sphere edge (high NPR2); the
   less-permeable one drifts disc-ward (lower NPR2).** Softest claim (shape = demoted descriptor).

## 3. Shape vs Rgyr — the size-normalization test (the key shape question)

Hypothesis: **shape (NPR) becomes useful iff, in large data, a size-normalized Rg window and a
size-dependent preferred shape both emerge and shape adds signal beyond PSA + normalized Rg.**
Three sub-claims:

- **(a) Normalize Rg by size; NPR is already normalized.** Rg ≈ N^ν grows with size, so a fixed
  Rg band can't transfer across different-size peptides — condition on size (per-residue, per-MW,
  or Rg/N^⅓). NPR1/NPR2 are inertia-moment *ratios* → largely size-independent by construction,
  which is the main argument for keeping shape as a *transferable* descriptor.
- **(b) Preferred shape may shift with size.** Small macrocycles are geometrically forced toward
  rod/elongated; larger ones can access disc and globular. The "permeable" PMI region may move
  with size — only visible with a size-diverse library.
- **(c) Incremental-value test (the crux).** In our 4 points higher Rg tracks with more disc-like
  (R: high Rg, low NPR2; S: low Rg, high NPR2) → Rg and shape may be **partly redundant**. Shape
  earns promotion **only if** NPR1/NPR2 improve a permeability model that already contains
  apolar-PSA + size-normalized Rg. If not, shape stays diagnostic.

## 4. Concrete CREMP analyses / acceptance criteria

- Compute the validated core (apolar 3D-PSA, Rg, total IMHB) + NPR1/NPR2 + asphericity over the
  CREMP permeability-annotated subset (their chloroform ensembles).
- **Absolute apolar-PSA threshold/score** vs measured permeability (calibrate the cutoff; compare
  to ΔPSA's added value — expect ΔPSA to add little once absolute PSA is in).
- **Size-conditioned Rg window:** bin by residue count / MW, fit permeable Rg/size band per bin.
- **Shape:** does a preferred NPR region exist, does it shift with size, and does it improve a
  PSA + normalized-Rg model (incremental R² / AUC)? Decision rule: promote shape only on a
  positive incremental test.
- Baseline per Liu 2025: physics-core regression vs a DMPNN graph model; our 3D descriptors are
  the value-add / physics-informed baseline.

---

*N = 2 caveat applies throughout. Related: [descriptor path forward](2026-06-23_descriptor_review_path_forward.md),
ML-descriptor implications (Limbach 2025 / Witek 2016).*

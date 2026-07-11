# HexPep NMR Validation — CREST vs. Rezai 2006 ³J couplings

**2026-07-08 · Jorge Carmona (Hu Lab, SDSU) · method validation**

> **HexPep is Rezai compound 1 (the *permeable* diastereomer, logP_E −6.2), confirmed by canonical SMILES.** Its CREST ensemble reproduces compound 1's experimental ³J(HN–Hα) couplings at **RMSD 0.71 Hz (uniform-mean, energy-free)** — essentially *at* the Vuister-Bax Karplus floor (0.73 Hz) — and matches its own diastereomer **3.4× better** than compound 9. So CREST is quantitatively faithful at the 6-mer scale (contrast: the 11-mer chameleon CsA, which it fails).

## Result

**Aggregation = uniform mean (energy-free);** metric = **RMSD** (to match the Karplus floor). Predicted ³J = per-conformer H–N–Cα–Hα dihedral → Karplus (Vuister & Bax 1993) → averaged over the CREST ensemble. Experimental couplings from the Rezai 2006 SI (±2 Hz, CDCl₃).

**Chloroform ensemble (solvent-matched, 3,295 conformers):**

| ³J(HN–Hα), Hz | Tyr | Leu (sorted) |
|---|---|---|
| CREST predicted (uniform mean) | 8.5 | 2.9, 6.3, 8.5, 9.8 |
| Experimental — compound 1 (= HexPep) | 8.4 | 4.0, 7.2, 9.0, 10.2 |

| CREST vs. | RMSD | MAE |
|---|---|---|
| **compound 1 (HexPep's diastereomer)** | **0.71 Hz** | 0.61 |
| compound 9 (impermeable control) | 2.40 Hz | 1.93 |

*Boltzmann-weighted (sensitivity check): 0.69 Hz vs compound 1 — ≈ uniform mean, so the result is energy-independent. Water ensemble gives the same picture. Vuister-Bax Karplus RMSD floor = 0.73 Hz.*

## What it means

- **Quantitative pass at the Karplus limit.** 0.71 Hz RMSD ≈ the 0.73 Hz floor (Vuister & Bax's own Karplus-vs-crystal-structure RMSD). The residual is the equation, not the ensemble — CREST captures HexPep's real solution conformation.
- **Resolves stereochemistry.** Matches its own diastereomer (compound 1) 3.4× better than the alternate (compound 9). Not a generic backbone — the correct fold for a fixed stereoisomer, which is the premise behind the R/S-epimer work.
- **Energy-free and robust.** Uniform mean ≈ Boltzmann (0.71 vs 0.69), so no dependence on the (still-pending) energies. Consistent with the pipeline-wide uniform-mean convention.
- **Anchors the validity envelope.** 6-mer → at the Karplus floor (good enough); 11-mer chameleon CsA → fails (0% cis, over-collapsed, needs explicit water).

## Identity & data-integrity notes

- **HexPep = Rezai compound 1**, confirmed by canonical-SMILES match after building all 9 Table-1 diastereomers (build→CIP self-check passed for all 9). Pattern: cyclo[D-Leu-D-Leu-L-Leu-D-Leu-L-Pro-L-Tyr].
- **The database is *not* wrong** (I initially thought it was): the entry's SMILES (compound 1), PAMPA (−6.20), and `permeable: False` are all consistent — −6.20 is below the project's −6.0 binary threshold, so "impermeable by threshold" is correct, even though compound 1 is Rezai's *most-permeable diastereomer*. The real error was mine: assuming "impermeable" ⇒ "compound 9."
- **Coupling-set assignment (corrected):** compound 1 = Leu 4.0/10.2/9.0/7.2, Tyr 8.4 (SI Table 2); compound 9 = Leu 9.0/7.8/8.4/6.6, Tyr 12.0 (SI Table 3). An earlier draft crossed these *and* mis-identified HexPep as compound 9; the two errors cancelled numerically, and both are now corrected.

## Method notes & remaining refinements

- **Karplus:** Vuister & Bax, *JACS* 1993, 115, 7772 — A = 6.51, B = −1.76, C = 1.60; 0.73 Hz RMSD vs. ubiquitin crystal φ (the floor). Parameterized on proteins; a standard choice.
- **Sorted-set Leu comparison** is still used (labeling-robust). Now that HexPep = compound 1 with a known fold, a **per-residue assignment** (Leu1–4 by ring position) is the next refinement.
- **NOE / D2O-exchange:** the SI's NOEs are counts-only; exchange is described in the main-paper prose (compound 1: only D-Leu2 exposed, Leu1/3/4/Tyr6 shielded) — a discrete H-bond-pattern check we can add as a second observable.
- **To complete the permeable-vs-impermeable pair:** generate **compound 9** (impermeable) — HexPep already covers the permeable side.

## Files & references

- `scripts/validate_hexpep_nmr.py` · `scripts/build_ensemble_from_crest.py`
- Ensembles: `results/conformers/HexPep/{aq,mem}/full_ensemble.xyz` · outputs: `results/2026-07-08_hexpep_nmr_validation{,_chcl3}.txt`
- Rezai, Yu, Millhauser, Jacobson, Lokey, *JACS* **2006**, 128, 2510 (+ SI). Vuister & Bax, *JACS* **1993**, 115, 7772 (Karplus). Companion: [CsA v2 validation](2026-07-07_csa_v2_validation.md).

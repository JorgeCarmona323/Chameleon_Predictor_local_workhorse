# DOPC 3-12-10-12 R/S — stereochemistry-driven conformational divergence resolved only by solvent-resolved 3D ensemble descriptors

**Report · 2026-06-16 · DOPC 3-12-10-12 R vs S (sarcosine backbone, xylylene linker)**

> **Bottom line.** The R and S stereoisomers are **identical on every 2D and lipophilicity
> descriptor**, yet their solvent-resolved 3D conformational ensembles diverge by up to ~70% —
> concentrated in water. The discriminating signal is the **amphipathic moment, shape anisotropy
> (asphericity), and intramolecular H-bonding (backbone, transannular)** — descriptors a 2D/TPSA
> model cannot see (the *solvent-accessible 3D*-PSA, however, does separate them; §2.2). **R is
> solvent-responsive** (opens in water, closes in membrane); **S is
> pre-organized** (closed in both). Notably, the dominant-conformer population (`p_dominant`) does
> *not* distinguish these isomers — the surface/H-bond/shape descriptors do. Hypothesis-generating;
> the experimental R-vs-S permeability is the decisive next datum.

---

## 1. Methods

**Compounds.** DOPC 3-12-10-12 R and S differ *only* at one backbone stereocenter (CIP R vs S);
the 6-residue macrocycle (sarcosine substituted for azetidine), the *m*-xylylene bis-thioether
linker, and every side chain are identical. Stereochemistry verified by RDKit CIP assignment; the
isomers share an identical 2D graph and differ only in the stereo layer.

> **[FIGURE 1 — insert: R/S 2D structure]** *Place the 2D structure panel with the stereocenter
> highlighted — the single point of difference between R and S. (Draw in ChemDraw/Illustrator.)*

**Conformer generation (CREST).** Independent ensembles per isomer in **water** and **chloroform**
(CHCl₃, ε≈4.8, membrane-interior mimic). Pipeline: RDKit ETKDGv3 (5000→50, MMFF94 prune) →
**GFN2-xTB/ALPB** pre-optimization → **CREST 2.12 iMTD-GC** (GFN2-xTB, ALPB, `--noreftopo --notopo`;
energy window 6 kcal/mol, cregen RMSD 0.125 Å), MTD 400–500 K → Boltzmann weighting at 298 K.
Protocol mirrors the CREMP dataset. Ensembles: R 567 / S 565 conformers (pooled over both solvents).

**Descriptors** (per conformer, Boltzmann-weighted, whole-molecule; `phys_descriptors_v3`):
2D/lipophilicity (MolWt, TPSA, Crippen LogP, HBD, HBA, rot. bonds, FractionCSP3); 3D shape (Rg,
asphericity, spherocity, NPR1/2); 3D surface (3D-PSA; **SA_HD** = donor-H surface, **SA_HA** =
acceptor surface; hydrophobic SASA; **amphipathic moment** = Å separation of polar vs apolar SASA
centroids); **weighted RMSF** (threshold-free ensemble flexibility); **IMHB** (geometric
H···A<2.5 Å, ∠>120°) split into **backbone (transannular)** vs
**side-chain** and into **donors/acceptors** engaged; cross-solvent Δ(water−CHCl₃). All values are
GFN2-energy-weighted (CREMP-consistent).

---

## 2. Results

### 2.1 2D and lipophilicity descriptors cannot distinguish the isomers

Every 2D/lipophilicity descriptor is identical between R and S. Any model on 2D/lipophilicity
features predicts identical permeability for both isomers.

### 2.2 The 3D ensembles diverge — and only in water

> **[FIGURE 2 — insert: `results/figures/isomers/reldiff_3-12-10-12.svg`]** *Relative |R−S|
> difference per descriptor (water phase). The 2D/lipophilicity descriptors (incl. TPSA) sit at ~0;
> the 3D ensemble descriptors extend out. Top discriminators: **asphericity (~70%), amphipathic
> moment (~59%), IMHB (~46%)**, with the **solvent-accessible 3D-PSA also separating them (~21%)**.
> Note the **dominant-conformer population separates them by only ~2%** — see §2.3.*

Membrane-phase descriptors barely differ; the signal is water-localized and multi-axis.

### 2.3 Mechanism — R solvent-responsive, S pre-organized; and a descriptor caveat

> **[FIGURE 3 — insert: `results/figures/isomers/hbonds_3-12-10-12.svg`]** *Per-conformer IMHB
> distributions (R/S × water/membrane). Same R-open/S-closed motif: in water R is open (low IMHB)
> and S is closed; the **backbone** panel shows S ≈ 3.0 transannular backbone H-bonds in water vs
> ≈ 2.0 for R.*

> **[FIGURE 4 — insert: `results/figures/isomers/overlap3d_3-12-10-12.svg`]** *Robust-scaled
> (median-centred, IQR units) distributions of the continuous 3D descriptors in water, R vs S — the
> boxes pull apart across asphericity, the surface descriptors, and the amphipathic moment.*

Key water-phase values:

| (water) | R | S | rel. diff |
|---|---|---|---|
| total IMHB | 3.58 | 5.70 | 46% |
| **backbone IMHB** | 1.98 | 2.98 | ~40% |
| side-chain IMHB | 1.61 | 2.72 | 53% |
| SA_HD (Å²) | 87.7 | 67.9 | 25% |
| amphipathic moment (Å) | 3.54 | 1.92 | 59% |
| asphericity | 0.043 | 0.089 | 70% |
| ensemble RMSF (Å, flexibility) | 0.55 | 0.44 | 22% |

**Flexibility / descriptor caveat (this scaffold).** Per-conformer `p_dominant` is **blind here
(≈2%)** — but that is a *discretization artifact*: it counts individual conformers, not folds, and
both isomers' single most-populated conformer holds only ~12% simply because CREST splits one fold
into many near-identical sub-states. The robust, threshold-free **weighted RMSF shows S is more
rigid than R in water (0.44 vs 0.55 Å)**, and 1-Å basin clustering confirms S sits in **one
dominant fold (98% of the population)** vs R's more multi-fold ensemble — i.e. S *is* the more
pre-organized isomer, consistent with its higher H-bonding. The discriminating descriptors are the
robust ones (amphipathic moment, asphericity, IMHB, SA_HD, RMSF); `p_dominant`/`n_eff` are not used.

> **[FIGURE 5 — insert: `results/figures/isomers/3d/openclosed_3-12-10-12_R.png`]** *PyMOL render of
> the R **water (open, marine)** vs **membrane (closed, salmon)** dominant-fold representatives.
> Optionally pair with `overlay_3-12-10-12_R_water.png` (top-20 conformational fan). Note S is the
> rigid one here (one fold, 98%) — `overlay_3-12-10-12_S_water.png` visually shows that tight fold.*

---

## 3. Interpretation

The stereocenter — invisible to every 2D/lipophilicity descriptor — drives a large, multi-axis,
water-localized difference in conformational behavior. The backbone/side-chain IMHB split localizes
S's pre-organization to the **transannular backbone** H-bond network — the canonical closed,
membrane-ready cyclic-peptide fold.

**Frustration hypothesis.** The stereocenter sets the pseudo-axial/equatorial disposition of the
side chain. In **S** the closed (high backbone-IMHB) fold and the side chain's preferred solvent
exposure are compatible → closed state strongly favored → pre-organized, solvent-insensitive. In
**R** the two cannot be satisfied at once → the ring opens in water and re-closes only in low
dielectric → solvent-responsive. A competing explanation (direct backbone φ/ψ stereoelectronic
preference) is distinguishable by a per-conformer side-chain-SASA vs backbone-IMHB correlation
(planned).

---

## 4. Reliability and next step

The **relative** R-vs-S differences are reproducible and come from the uncapped water ensembles
(moderate confidence). Absolute values are less certain — single-start CREST, implicit solvent,
sub-threshold 6-mer; no experimental structure exists. This is **hypothesis-generating**, not proven.

The decisive datum is the **experimental R-vs-S permeability**: if S (pre-organized) is more
permeable, it supports the pre-organization mechanism; if R (solvent-responsive) is, it points to
chameleonic behavior even at this size.

---

*Figures: `scripts/plot_isomer_figures.py` → `results/figures/isomers/reldiff_3-12-10-12.svg`,
`hbonds_3-12-10-12.svg`, `overlap3d_3-12-10-12.svg` (SVG = Illustrator-editable). Descriptor library:
`scripts/phys_descriptors_v3.py`; v3 rationale: `docs/experiments/2026-06-13_descriptor_literature_review.md`.*

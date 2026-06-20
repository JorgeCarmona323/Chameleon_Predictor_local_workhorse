# DOPC 3-12-8-12 R/S — stereochemistry-driven conformational divergence resolved only by solvent-resolved 3D ensemble descriptors

**Report · 2026-06-16 · DOPC 3-12-8-12 R vs S (azetidine backbone, xylylene linker)**

> **Bottom line.** The R and S stereoisomers are **provably identical on every 2D and lipophilicity
> descriptor, including TPSA**, yet their solvent-resolved 3D conformational ensembles diverge by up
> to ~60% — concentrated entirely in water. The discriminating signal is **intramolecular H-bonding
> (specifically the backbone, transannular network), the amphipathic moment, and — once donor-H
> exposure is counted (the Ono/Begnini solvent-accessible 3D-PSA definition) — the 3D-PSA itself**;
> all are invisible to a 2D/TPSA model. **R is solvent-responsive** (opens in water, closes in
> membrane); **S is pre-organized** (closed in both). Hypothesis-generating; the experimental
> R-vs-S permeability is the decisive next datum.

---

## 1. Methods

**Compounds.** DOPC 3-12-8-12 R and S differ *only* at one backbone stereocenter (CIP R vs S);
the 6-residue macrocycle, the *m*-xylylene bis-thioether linker, and every side chain are identical.
Stereochemistry verified by RDKit CIP assignment; the isomers share an identical 2D graph and differ
only in the stereo layer.

> **[FIGURE 1 — insert: R/S 2D structure]** *Place the 2D structure panel with the stereocenter
> highlighted — the single point of difference between R and S. (Draw in ChemDraw/Illustrator.)*

**Conformer generation (CREST).** Independent ensembles per isomer in **water** and **chloroform**
(CHCl₃, ε≈4.8, membrane-interior mimic). Pipeline: RDKit ETKDGv3 (5000→50, MMFF94 prune) →
**GFN2-xTB/ALPB** pre-optimization → **CREST 2.12 iMTD-GC** (GFN2-xTB, ALPB, `--noreftopo --notopo`;
energy window 6 kcal/mol, cregen RMSD 0.125 Å), MTD 400–500 K → Boltzmann weighting at 298 K.
Protocol mirrors the CREMP dataset. Ensembles: R 529 / S 481 conformers (pooled over both solvents).

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

Every 2D/lipophilicity descriptor is identical between R and S (MolWt 801.97, TPSA 249.36, Crippen
LogP −1.49, HBD 8, HBA 12, …) — including the topological PSA. Any model on 2D/TPSA features
predicts identical permeability for both isomers. (The *solvent-accessible 3D*-PSA is a different
quantity and *does* separate them — see §2.3 — because it counts donor-H exposure, which the
stereocenter modulates.)

### 2.2 The 3D ensembles diverge — and only in water

> **[FIGURE 2 — insert: `results/figures/isomers/reldiff_3-12-8-12.svg`]** *Relative |R−S|
> difference per descriptor (water phase). The 2D/lipophilicity descriptors (incl. TPSA) sit at ~0;
> the 3D ensemble descriptors extend out. Top discriminators: **backbone IMHB (~59%), amphipathic
> moment (~57%), IMHB acceptors (~56%), SA_HD (~54%)**, with the **solvent-accessible 3D-PSA also
> separating the epimers (~31%)** — the surface and H-bond-resolved descriptors, not the 2D terms.*

Membrane-phase descriptors barely differ; the signal is water-localized and multi-axis.

### 2.3 Mechanism — R solvent-responsive, S pre-organized (backbone H-bonding)

> **[FIGURE 3 — insert: `results/figures/isomers/hbonds_3-12-8-12.svg`]** *Per-conformer IMHB
> distributions (R/S × water/membrane). In water R is open (low IMHB, diffuse) and S is closed (high
> IMHB); both converge to the closed form in membrane. The **backbone** panel localizes S's
> pre-organization: S ≈ 2.9 transannular backbone H-bonds in water vs ≈ 1.6 for R.*

> **[FIGURE 4 — insert: `results/figures/isomers/overlap3d_3-12-8-12.svg`]** *Robust-scaled
> (median-centred, IQR units) distributions of the continuous 3D descriptors in water, R vs S — the
> boxes pull apart across surface and shape axes simultaneously.*

Key water-phase values:

| (water) | R | S | rel. diff |
|---|---|---|---|
| total IMHB | 4.34 | 6.75 | 44% |
| **backbone IMHB** | 1.60 | 2.94 | ~59% |
| side-chain IMHB | 2.74 | 3.82 | 33% |
| SA_HD (Å²) | 78.8 | 45.0 | 54% |
| amphipathic moment (Å) | 3.42 | 1.91 | 57% |
| ensemble RMSF (Å, flexibility) | 0.62 | 0.96 | 43% |
| 3D-PSA (Å², Ono/Begnini def) | 227 | 166 | 31% |
| TPSA (2D) | 249.4 | 249.4 | 0% |

*Flexibility note:* by **weighted RMSF** (threshold-free), S samples a *broader* geometric ensemble
than R in water (0.96 vs 0.62 Å) **despite being more H-bonded** — i.e. S explores many closed-ish
sub-folds rather than one rigid fold. (Per-conformer `p_dominant` is discretization-sensitive — it
counts individual conformers, not folds — so it is not used here; weighted RMSF and 1-Å basin
clustering are the robust flexibility metrics. See `docs/experiments/2026-06-16_*` analysis.)

> **[FIGURE 5 — insert: `results/figures/isomers/3d/openclosed_3-12-8-12_R.png`]** *PyMOL render of
> the R **water (open, marine)** vs **membrane (closed, salmon)** dominant-fold representatives — the
> solvent response in 3D. Optionally pair with `overlay_3-12-8-12_R_water.png` (the conformational
> "fan" of the top-20 water conformers, dominant fold = 86% of population) to show ensemble spread.*

---

## 3. Interpretation

The stereocenter — invisible to every 2D/lipophilicity descriptor — drives a large, multi-axis,
water-localized difference in conformational behavior. A useful subtlety distinguishes the two
notions of "polar surface": the **topological PSA (TPSA) is identical** between R and S (same 2D
graph), but the **solvent-accessible 3D-PSA, which counts the donor-H surface, separates them**
(R 227 vs S 166 Å²). S's extra intramolecular H-bonds bury its donor *hydrogens*, lowering its
3D-PSA — exactly the pattern Begnini 2021 reported for their more-permeable epimer (lower 3D-PSA).
This is why a 2D/TPSA model is blind here while the 3D, donor-resolved surface (3D-PSA, SA_HD, IMHB)
sees the stereocenter.

The backbone/side-chain IMHB split localizes S's pre-organization to the **transannular backbone**
H-bond network — the canonical closed, membrane-ready cyclic-peptide fold.

**Frustration hypothesis.** The stereocenter sets the pseudo-axial/equatorial disposition of the
side chain. In **S** the closed (high backbone-IMHB) fold and the side chain's preferred solvent
exposure are compatible → closed state strongly favored → pre-organized, solvent-insensitive. In
**R** the two cannot be satisfied at once → the ring opens in water to relieve the conflict and
re-closes only in low dielectric → solvent-responsive. A competing explanation (direct backbone
φ/ψ stereoelectronic preference) is distinguishable by a per-conformer side-chain-SASA vs
backbone-IMHB correlation (planned).

---

## 4. Reliability and next step

The **relative** R-vs-S differences are reproducible and come from the uncapped water ensembles
(moderate confidence). Absolute values (especially ΔPSA) are less certain — single-start CREST,
implicit solvent, sub-threshold 6-mer; no experimental structure exists. This is
**hypothesis-generating**, not proven.

The decisive datum is the **experimental R-vs-S permeability**: if S (pre-organized) is more
permeable, it supports the pre-organization mechanism; if R (solvent-responsive) is, it points to
chameleonic behavior even at this size.

---

*Figures: `scripts/plot_isomer_figures.py` → `results/figures/isomers/reldiff_3-12-8-12.svg`,
`hbonds_3-12-8-12.svg`, `overlap3d_3-12-8-12.svg` (SVG = Illustrator-editable). Descriptor library:
`scripts/phys_descriptors_v3.py`; v3 rationale: `docs/experiments/2026-06-13_descriptor_literature_review.md`.*

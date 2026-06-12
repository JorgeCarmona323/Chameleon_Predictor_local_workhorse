# DOPC 3-12-8-12 R/S Isomers — Method & Results

**For PI review · 2026-06-10 · keep updated**

> **Bottom line:** Two stereoisomers that are *provably identical* on every 2D and lipophilicity descriptor differ by up to **98%** in their 3D conformational ensembles — concentrated entirely in the aqueous environment. The stereocenter makes **R solvent-responsive** (opens in water) and **S pre-organized** (closed in both water and membrane). Only solvent-resolved 3D ensemble descriptors capture this; no 2D or PSA-based model can. This is hypothesis-generating evidence; the experimental R-vs-S permeability is the decisive next data point.

---

## Methods

**Compounds.** DOPC 3-12-8-12 R and S are stereoisomers differing only at the thiophene-bearing backbone stereocenter (CIP-R vs CIP-S; assignments verified by RDKit CIP, identical InChIKey skeleton `KAPGDKUZKGUOTI`, differing only in the stereo layer). All other atoms — the 6-residue macrocycle, the dibenzyl-thioether bridge, every side chain — are identical.

> **[FIGURE 1 — molecule structures]** *Place the R/S 2D structure panel here (thiophene stereocenter highlighted). Shows the single point of difference between the isomers.*

**Conformer generation (CREST conditions).** For each isomer, conformational ensembles were generated independently in two solvents — **water** (high dielectric, aqueous phase) and **chloroform / CHCl₃** (low dielectric, ε ≈ 4.8, membrane-interior mimic). The per-molecule pipeline:

1. **Initial geometries — RDKit ETKDGv3:** embed up to 5,000 conformers (`numConfs=5000`, `useRandomCoords=True`, macrocycle torsion preferences), MMFF94 minimize, sort by energy, prune by heavy-atom RMSD (0.5 Å) to ≤ 50 representatives.
2. **xTB pre-optimization:** each retained conformer optimized with **GFN2-xTB** + **ALPB** implicit solvent (the run's solvent); the lowest-energy optimized structure seeds CREST.
3. **CREST conformer search — version 2.12:** iterative metadynamics + genetic structure crossing (**iMTD-GC**), **GFN2-xTB**, **ALPB** solvent, flags **`--noreftopo --notopo`** (required so the macrocycle can flip/sample cis-trans amide states without false "broken-ring" termination). CREST defaults retained: energy window **6.0 kcal/mol**, cregen RMSD threshold **0.125 Å**, inter-conformer energy threshold **0.05 kcal/mol**. Metadynamics/MD runs at elevated temperatures (**400–500 K**) for enhanced barrier crossing; final conformer energies are GFN2-xTB values.
4. **Ensembles:** Boltzmann-weighted at **298 K** (RT = 0.592 kcal/mol) from the GFN2-xTB energies. The CHCl₃ ensemble post-processing was capped at the 50 lowest-energy conformers (negligible effect on Boltzmann-weighted means). Output per solvent: `ensemble.sdf` (3D coordinates) + `ensemble.json` (energies, weights, per-conformer PSA/HB).

CREST 2.12 was used deliberately — CREST 3.x produced reproducible `crest_xtbsp` crashes on these macrocycles. The overall protocol mirrors the published CREMP dataset workflow.

**3D descriptor generation (`ensemble_descriptors.py`).** Computed per solvent, **Boltzmann-weighted over the full ensemble**, on the **whole molecule** (backbone + all side chains):
- **3D polar surface area (PSA):** solvent-accessible surface of polar heavy atoms (N, O, S, P) via RDKit `rdFreeSASA` with Bondi radii.
- **Intramolecular H-bonds:** geometric count (donor H···acceptor < 2.5 Å, D–H···A angle > 120°) over all N/O donors and acceptors.
- **Shape:** radius of gyration, NPR1/NPR2, asphericity, spherocity (RDKit `Descriptors3D`).
- **Conformational concentration:** Boltzmann population of the dominant (highest-weight) conformer.
- **cis-amide propensity:** ω dihedral per backbone ring amide bond (cis if |ω| < 30°).
- **Cross-solvent:** Δ(water − CHCl₃) for each, plus `norm_delta_psa` (ΔPSA / total SASA).

**2D / lipophilicity descriptors (`compute_2d_descriptors.py`).** Computed from the 2D graph only (conformer-independent): MolWt, TPSA, Crippen LogP, Crippen MolMR, HBD, HBA, rotatable bonds, FractionCSP3, LabuteASA, QED.

---

## Results

### 1. 2D and lipophilicity descriptors cannot distinguish the isomers

All 13 canonical 2D / lipophilicity descriptors are identical between R and S (e.g. MolWt 801.97, TPSA 249.36, Crippen LogP −1.49, MolMR 203.24, HBD 8, HBA 12 — all the same). By construction, any model built on 2D or lipophilicity features predicts identical permeability for both isomers.

> **[FIGURE 2 — 2D overlap]** *Place `overlap_2d.svg` here. Bars for R and S coincide exactly across all descriptors — visual proof of 2D blindness to the stereocenter.*

### 2. 3D conformational ensembles differ markedly — and only in water

Relative differences between the isomers are large for several independent 3D descriptors, all concentrated in the aqueous ensemble: dominant-conformer population (98%), intramolecular H-bonds (44%), shape anisotropy (asphericity 18%, spherocity 15%). Membrane-phase descriptors barely differ, and even the 3D water PSA is identical. The signal is multi-dimensional (H-bonding, shape, and conformational concentration) and solvent-specific.

> **[FIGURE 3 — 2D vs 3D relative difference]** *Place `rel_diff_2d_vs_3d.svg` here. All 2D descriptors at 0%; 3D ensemble descriptors extend to 98%.*

### 3. R is solvent-responsive; S is pre-organized

Per-conformer H-bond distributions reveal the mechanism. In water, R's dominant conformer is open (4 intramolecular H-bonds) with a diffuse ensemble (no dominant state, p = 0.08); S's is closed (7 H-bonds) and concentrated (p = 0.23). In membrane, both isomers converge to the closed, 7-H-bond form.

| | water | membrane |
|---|---|---|
| **R** | open, 4 H-bonds, diffuse | closed, 7 H-bonds |
| **S** | closed, 7 H-bonds, concentrated | closed, 7 H-bonds |

→ **R changes conformation with its environment (opens in water, closes in membrane); S is locked in the closed, membrane-ready form in both.**

> **[FIGURE 4 — H-bond box plots]** *Place `box_hbonds.svg` here. R-water is the lone low/spread group; S-water and both membrane forms cluster high (closed).*

> **[FIGURE 5 — optional: representative 3D conformers]** *If desired, place the PyMOL render of R-water (open) vs R-membrane (closed) here to illustrate the solvent response in 3D. (Raster image; trace in AI if vector editing is needed.)*

---

## Interpretation

The stereocenter — invisible to every 2D and lipophilicity descriptor — drives a large, multi-axis, water-localized difference in conformational behavior. This fits the size-gated permeability hypothesis: as 6-mers below the ~9-residue chameleonicity threshold, the discriminating feature is **pre-organization / intramolecular H-bonding**, not chameleonic ΔPSA switching. It also demonstrates that solvent-resolved 3D ensemble descriptors carry signal orthogonal to (and undetectable by) the standard 2D feature set.

**A useful subtlety:** S forms more intramolecular H-bonds in water yet has identical water PSA to R. PSA measures exposure of polar *heavy atoms*; an intramolecular H-bond sequesters the *hydrogen* while the heavy atoms stay solvent-exposed. The two descriptors are therefore complementary, not redundant.

---

## Reliability & next step

The R-vs-S *relative* difference is reproducible and comes from the uncapped water ensembles (moderate confidence). Absolute values (especially ΔPSA) are less certain — single-start CREST, implicit solvent, sub-threshold 6-mer. There is no experimental structure for these isomers, so this is **hypothesis-generating evidence**, not a proven result.

The decisive next data point is the **experimental R-vs-S permeability**: if S (pre-organized) is more permeable, it supports the pre-organization mechanism; if R (solvent-responsive) is, it would point to chameleonic behavior even at this size.

---

*Figure source files (`results/figures/isomers/`, SVG = Illustrator-editable): `overlap_2d.svg`, `rel_diff_2d_vs_3d.svg`, `box_hbonds.svg`, plus `population_profile.svg` / `psa_distribution.svg` as supplements. Detailed log: `docs/experiments/2026-06-05_dopc_rs_3d_vs_2d_descriptors.md`.*

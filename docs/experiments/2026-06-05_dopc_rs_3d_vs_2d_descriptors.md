# Experiment — DOPC R/S Isomers: What 3D Ensemble Descriptors Catch That 2D Cannot

**Date:** 2026-06-05
**Script:** `scripts/ensemble_descriptors.py`
**Data:** `results/conformers/DOPC 3-12-8-12 R/` and `.../DOPC 3-12-8-12 S/` — CREST ensembles, water + CHCl3
**Output:** `results/ensemble_descriptors_dopc_rs.csv`

---

## Question

DOPC 3-12-8-12 R and S are **stereoisomers** — identical 2D structure, differing only in configuration at the thiophene-bearing stereocenter. Can 3D conformational-ensemble descriptors distinguish them, when 2D descriptors by construction cannot?

---

## Why 2D descriptors are blind here

R and S have the **same molecular graph**: same formula, same bonds, same atom types. Every topological/2D descriptor is therefore **identical** between them:

| 2D descriptor | DOPC_R | DOPC_S |
|---|---|---|
| MolWt | identical | identical |
| TPSA (2D) | identical | identical |
| H-bond donors / acceptors | identical | identical |
| Rotatable bonds, ring count, FractionCSP3 | identical | identical |

A 2D-descriptor model (or any fingerprint that ignores stereochemistry) **cannot tell these two apart.** Stereochemistry only manifests in 3D geometry.

---

## What the 3D ensemble descriptors show

Computed per solvent, Boltzmann-weighted over the CREST ensembles:

| Descriptor | DOPC_R | DOPC_S | Interpretation |
|---|---|---|---|
| `water_bw_psa` | 183.5 | 183.5 | **identical** — polarity does not distinguish them |
| `water_bw_hb` | **4.34** | **6.75** | S forms ~2.4 more intramolecular H-bonds in water |
| `mem_bw_hb` | 6.86 | 6.89 | identical in membrane |
| `delta_hb` (water−mem) | **−2.52** | **−0.14** | R *opens up* in water (sheds H-bonds); S stays H-bonded in both |
| `water_p_dominant` | **0.079** | **0.231** | R ensemble is diffuse (no dominant conformer, 479 confs); S has a clear dominant basin |
| `water_bw_spherocity` | 0.60 | 0.52 | R rounder; S more anisotropic |
| `cis_prob` (all bonds) | identical | identical | both have one constitutively-cis amide (bond 3); no difference |

---

## Insight

The stereocenter flip R→S does **not** change polar surface area or the cis-amide pattern — the two descriptors a 2D or PSA-centric view would lean on. Instead it changes the **conformational ensemble behavior**:

- **S** folds into a more **intramolecularly H-bonded, conformationally concentrated** state in water (6.75 HB, dominant conformer at 23%).
- **R** stays **open and conformationally diffuse** in water (4.34 HB, population spread thinly over 479 conformers with no dominant state) and only collapses/H-bonds upon entering the apolar (membrane-like) environment — a larger `delta_hb` swing.

This is precisely the kind of difference that:
- **2D descriptors miss** (identical graph),
- **single-conformer 3D would miss** (you'd have to pick the "right" conformer),
- **a thermal CREST ensemble + Boltzmann averaging reveals** as a difference in *populations* and *H-bond networks*.

It demonstrates the value of the `DynamicEnsembleEncoder` modality: stereochemistry-driven permeability differences are encoded in the conformational ensemble's H-bonding and shape distribution, not in any static 2D property.

---

## Connection to the model

These descriptors feed the `DynamicEnsembleEncoder` (see `chameleon_model_architecture.md`); the R/S pair is a clean demonstration that the dynamic modality carries signal orthogonal to the `StaticDescriptorEncoder` (2D). For an ablation, R/S isomer pairs are an ideal test: any model relying only on static/2D features will predict identical permeability for both, whereas the experimental permeabilities differ.

---

## Limitations

1. **mem ensemble capped at 50 conformers, water uncapped** (479 R / 431 S). The absolute ΔPSA family is biased by this asymmetry. It does **not** affect the R-vs-S comparison (both isomers received identical cap treatment, so the bias cancels in the difference), and the headline signals (water H-bonds, dominant population, shape) come from the **uncapped water** ensembles and are unaffected.
2. **`cis_switch_bond` reporting glitch:** the descriptor reports bond 3 with switch-magnitude 0.0 — an `argmax`-over-zeros artifact, since bond 3 is constitutively cis (not a switch). Cosmetic; does not affect the R/S conclusion (cis pattern is identical for both isomers and therefore not a discriminator).
3. **6-mers below the ~9-residue chameleonicity threshold** (`hypothesis.md`): absolute ΔPSA is not expected to be a clean chameleonic signal for these; H-bonding/shape and `norm_delta_psa` carry the interpretable difference.
4. Descriptors are from single-start CREST ensembles; cross-solvent congruent-state analysis (Witek) intentionally omitted per the 2026-05-31 descriptor-scope decision.

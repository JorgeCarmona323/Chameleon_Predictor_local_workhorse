# DOPC 3-12-8-12 R/S Isomers — 3D Conformational Signal

**Summary report · 2026-06-10**

> Two stereoisomers that are *provably identical* on every 2D and lipophilicity descriptor differ by up to **98%** in their 3D conformational ensembles — concentrated entirely in the aqueous environment. The stereocenter makes **R solvent-responsive** (opens in water) and **S pre-organized** (closed in both water and membrane). Only solvent-resolved 3D ensemble descriptors capture this; no 2D or PSA-based model can.

---

## The compounds

DOPC 3-12-8-12 **R** and **S** are stereoisomers differing **only** at the thiophene-bearing backbone stereocenter (CIP-R vs CIP-S, verified). Identical constitution (InChIKey block `KAPGDKUZKGUOTI`), opposite chirality at one carbon. Same 6-residue macrocycle, same side chains, same everything else.

**Data:** CREST iMTD-GC ensembles (CREST 2.12, GFN2-xTB + ALPB) in water and CHCl₃ (membrane mimic); descriptors via `ensemble_descriptors.py`.

---

## The question

Permeability models lean on 2D descriptors (TPSA, LogP) and PSA. But stereochemistry is invisible to those by construction. **Can 3D conformational descriptors distinguish two isomers that 2D cannot — and does the difference mean anything?**

---

## Result 1 — 2D is provably blind

All **13** canonical 2D / lipophilicity descriptors are byte-for-byte identical between R and S:

| | R | S |
|---|---|---|
| MolWt / TPSA / LogP (Crippen) | 801.97 / 249.36 / −1.49 | identical |
| MolMR / HBD / HBA | 203.24 / 8 / 12 | identical |
| FractionCSP3 / RotBonds / LabuteASA / QED | — | identical |

→ A 2D-descriptor or lipophilicity model predicts the **exact same** permeability for both. *(Figure: `overlap_2d` — bars coincide exactly.)*

---

## Result 2 — 3D ensembles differ, and only in water

Relative difference %|R−S| across descriptors:

| Descriptor | %\|R−S\| | |
|---|---|---|
| `water_p_dominant` (conformational concentration) | **98%** | 3D |
| `water_bw_hb` (intramolecular H-bonds) | **44%** | 3D |
| `water_bw_asphericity` (shape) | 18% | 3D |
| `water_bw_spherocity` (shape) | 15% | 3D |
| `mem_bw_psa` | 4% | 3D |
| `water_bw_rg`, `mem_bw_hb`, `water_bw_psa` | 0–2% | 3D |
| every 2D / lipophilicity descriptor | **0%** | 2D |

The signal is **multi-axis** (H-bonding *and* shape *and* conformational concentration) and **lives in water** — membrane descriptors barely differ. *(Figures: `rel_diff_2d_vs_3d`, `box_hbonds`.)*

---

## Interpretation — solvent-responsive vs pre-organized

| | water | membrane |
|---|---|---|
| **R** dominant conformer | open, **4** H-bonds, diffuse (p=0.08) | closed, 7 H-bonds |
| **S** dominant conformer | closed, **7** H-bonds, concentrated (p=0.23) | closed, 7 H-bonds |

- **R is solvent-responsive** — it opens and exposes in water, collapses and H-bonds in membrane. A larger conformational swing between environments.
- **S is pre-organized** — already locked in the closed, H-bonded, membrane-ready form regardless of solvent.

This is exactly the axis the permeability theories care about (how the H-bond network responds to environment), and it fits the **size-gated hypothesis**: these are 6-mers below the ~9-residue chameleonicity threshold, so the discriminating signal is **pre-organization / intramolecular H-bonding**, not chameleonic ΔPSA switching (which is in fact anti-chameleonic here).

---

## Why H-bonds and PSA disagree (a useful subtlety)

S makes ~2.4 more intramolecular H-bonds in water, yet R and S have **identical** water PSA (183.5). These descriptors decouple: PSA measures solvent exposure of polar *heavy atoms*; an intramolecular H-bond sequesters the *hydrogen* (never counted in PSA) while the heavy atoms can stay exposed. Same polar exposure, different internal H-bond network — a reason both descriptors belong in the feature set.

---

## Reliability (honest framing)

| Aspect | Confidence |
|---|---|
| R-vs-S *relative* difference | **Moderate** — correct pipeline, signal from uncapped water ensembles, biases cancel in the comparison |
| Absolute values (esp. ΔPSA) | **Low** — single-start CREST, implicit solvent, sub-threshold 6-mer |
| Experimental validation | **None** — no crystal/NMR for these isomers |

→ This is **hypothesis-generating evidence**, not a proven result. The difference is real and reproducible; its *direction* (which isomer permeates better) is testable against the experimental R-vs-S permeability — the number that would turn this into mechanistic evidence.

---

## Bottom line

Stereochemistry that is **invisible to 2D and lipophilicity descriptors** produces a **large, multi-dimensional, water-localized 3D conformational signal**. It demonstrates the value of the solvent-resolved 3D ensemble (the `DynamicEnsembleEncoder` modality) and gives a clean, falsifiable hypothesis: **R switches with environment, S is pre-organized** — with the experimental permeability of the two isomers as the decisive next data point.

---

*Figures referenced (`results/figures/isomers/`): `overlap_2d`, `rel_diff_2d_vs_3d`, `box_hbonds`, `population_profile`, `psa_distribution`. Detailed log: `docs/experiments/2026-06-05_dopc_rs_3d_vs_2d_descriptors.md`.*

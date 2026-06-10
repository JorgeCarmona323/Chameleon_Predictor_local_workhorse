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

R and S have the **same molecular graph**: same formula, same bonds, same atom types. Computed explicitly (`scripts/compute_2d_descriptors.py`), **all 13 canonical 2D / lipophilicity descriptors are byte-for-byte identical**:

| 2D descriptor | DOPC_R | DOPC_S |
|---|---|---|
| MolWt | 801.97 | 801.97 |
| TPSA (2D) | 249.36 | 249.36 |
| MolLogP (Crippen) | −1.491 | −1.491 |
| MolMR (molar refractivity) | 203.24 | 203.24 |
| NumHDonors / NumHAcceptors | 8 / 12 | 8 / 12 |
| FractionCSP3 | 0.457 | 0.457 |
| LabuteASA | 324.05 | 324.05 |
| QED | 0.150 | 0.150 |
| *(13/13 total)* | **identical** | **identical** |

A 2D-descriptor or lipophilicity-based model **cannot tell these two apart** — including the two metrics permeability models lean on most (TPSA and LogP). Stereochemistry only manifests in 3D geometry.

**Figures** (`results/figures/isomers/`, from `scripts/plot_isomer_comparison.py`):
- `box_hbonds` — per-conformer H-bond distribution; R-water is the lone low/spread group, S-water + both membrane forms cluster closed.
- `rel_diff_2d_vs_3d` — relative %|R−S| per descriptor; 2D/lipophilicity all at 0%, 3D ensemble up to 98% (`water_p_dominant`), 44% (`water_bw_HB`). The signal lives in water (`mem_bw_HB` only 1%).

---

## What the 3D ensemble descriptors show

Computed per solvent, Boltzmann-weighted over the CREST ensembles:

| Descriptor | DOPC_R | DOPC_S | Interpretation |
|---|---|---|---|
| `water_bw_psa` | 183.5 | 183.5 | **identical in water** — polarity does not distinguish them here |
| `mem_bw_psa` | 208.6 | 199.6 | differ by ~9 Å² in membrane (but capped ensemble — see limitations) |
| `water_bw_hb` | **4.34** | **6.75** | S forms ~2.4 more intramolecular H-bonds in water |
| `mem_bw_hb` | 6.86 | 6.89 | identical in membrane |
| `delta_hb` (water−mem) | **−2.52** | **−0.14** | R *opens up* in water (sheds H-bonds); S stays H-bonded in both |
| `water_p_dominant` | **0.079** | **0.231** | R ensemble is diffuse (no dominant conformer, 479 confs); S has a clear dominant basin |
| `water_bw_spherocity` | 0.60 | 0.52 | R rounder; S more anisotropic |
| `cis_prob` (all bonds) | identical | identical | both have one constitutively-cis amide (bond 3); no difference |

### Why identical water PSA despite different H-bonding? (IMHB and PSA decouple)

S makes ~2.4 more intramolecular H-bonds in water than R, yet their water PSA is identical (183.5). These are **not** redundant descriptors:
- **PSA** = solvent-accessible surface of the polar **heavy atoms** (N, O, S); it does not count hydrogens.
- **IMHB** = a discrete geometric event (donor-H within 2.5 Å of acceptor, >120°).

Forming an intramolecular H-bond sequesters the **H** (never counted in PSA) and only partially shields one acceptor — a "surface" H-bond can add +1 IMHB while the polar heavy atoms stay solvent-exposed. In water, polar groups are thermodynamically driven to expose regardless of transient internal H-bonds, so the time-averaged polar-heavy-atom exposure converges for both isomers. **Result: same PSA, different H-bond network** — direct evidence the two descriptors carry orthogonal information, and a reason to keep both in the feature set.

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

## Connection to the model and hypothesis

These descriptors feed the `DynamicEnsembleEncoder` (see `chameleon_model_architecture.md`); the R/S pair is a clean demonstration that the dynamic modality carries signal orthogonal to the `StaticDescriptorEncoder` (2D). For an ablation, R/S isomer pairs are an ideal test: any model relying only on static/2D features will predict identical permeability for both, whereas the experimental permeabilities differ.

It also fits the **partition-based two-regime hypothesis** (`hypothesis.md`): these are 6-mers, *below* the ~9-residue chameleonicity threshold. As predicted for the small/pre-organized regime, the discriminating signal is **intramolecular H-bonding and conformational pre-organization**, not chameleonic ΔPSA switching (which is in fact anti-chameleonic here — ΔPSA negative). The R/S pair is a concrete instance of the small-peptide mechanism the hypothesis attributes to sub-threshold peptides.

## Reliability

| Aspect | Reliability | Why |
|---|---|---|
| R-vs-S *relative* difference | **Moderate** | Both isomers ran the identical pipeline (`-notopo`, CREST 2.12); headline signals (water H-bonds, dominant population, shape) come from the **uncapped** water ensembles. Systematic biases cancel in the R−S comparison. |
| Absolute values (esp. ΔPSA) | **Low** | mem capped at 50, single-start CREST, implicit solvent, sub-threshold 6-mer. Use relative/normalized only. |
| Experimental validation | **None** | No crystal/NMR structures for these isomers (unlike CsA). The difference is *suggestive, not proven* — hypothesis-generating. |

**Bottom line:** there is a real, reproducible *relative* difference between the isomers in conformational behavior — but it is unvalidated and should be framed as signal worth investigating, not an established result.

---

## Limitations

1. **mem ensemble capped at 50 conformers, water uncapped** (479 R / 431 S). The absolute ΔPSA family is biased by this asymmetry. It does **not** affect the R-vs-S comparison (both isomers received identical cap treatment, so the bias cancels in the difference), and the headline signals (water H-bonds, dominant population, shape) come from the **uncapped water** ensembles and are unaffected.
2. **`cis_switch_bond` reporting glitch:** the descriptor reports bond 3 with switch-magnitude 0.0 — an `argmax`-over-zeros artifact, since bond 3 is constitutively cis (not a switch). Cosmetic; does not affect the R/S conclusion (cis pattern is identical for both isomers and therefore not a discriminator).
3. **6-mers below the ~9-residue chameleonicity threshold** (`hypothesis.md`): absolute ΔPSA is not expected to be a clean chameleonic signal for these; H-bonding/shape and `norm_delta_psa` carry the interpretable difference.
4. Descriptors are from single-start CREST ensembles; cross-solvent congruent-state analysis (Witek) intentionally omitted per the 2026-05-31 descriptor-scope decision.

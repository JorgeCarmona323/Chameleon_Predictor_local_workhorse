# Descriptor Review & Finalized Path Forward (2026-06-23)

**Trigger:** Edison ranked, evidence-graded literature review
(`docs/literature/3d_descriptors literature/edison_3d_descriptors_review.pdf`),
which scores every conformation-dependent 3D descriptor by *experimental*
validation strength for cyclic-peptide / bRo5 passive permeability and maps them
onto our pipeline. This doc records the decisions that review drove — the goal is
to **cut speculative depth, lead with validated descriptors, and finalize the ML
track**.

---

## 1. Decisions (locked)

| Question | Decision |
|---|---|
| Descriptor scope | **Demote, don't delete.** Reports lead with the validated core; weak/exploratory descriptors move to an optional *diagnostics* table with no interpretive narrative. Code stays intact. |
| ML target | **Pivot to CREMP.** Our compounds = mechanistic/validation case study (small N, never ML training). The predictive ML model is built on CREMP's 3,258 permeability-annotated macrocycles. |
| Solvent | **Superseded 2026-07-02 — see §3 update.** Now: native per-phase CREST in **water + cyclohexane** (Ono's permeability solvent); chloroform demoted to NMR/IMHB anchor. Score with CPCM-X. |

---

## 2. The validated core (what reports lead with)

The review's minimal non-redundant set, all computed from the **apolar
(chloroform / hydrocarbon) Boltzmann-weighted ensemble**:

1. **SA 3D-PSA (apolar)** — polarity burden. Best-reproduced single descriptor
   (Begnini: correctly ranked 7 isomers; reproduced in PROTAC series). *We compute this.*
2. **Radius of gyration (apolar)** — compactness/folding; validated companion to PSA. *We compute this.*
3. **IMHB count, backbone-transannular (apolar)** — donor-shielding mechanism
   (Rezai ~100-fold). Qualitative; weak standalone but mechanistically core. *We compute this.*
4. **Cross-solvent Δ / ΔG_transfer (water→apolar)** — chameleonicity; the single
   *strongest* literature correlation (Kamenik r=0.92, ρ=1.0; collapses to r=0.50
   with single structures). **The #1 descriptor we are currently missing.** See §3.

Supporting (keep, secondary): total SASA (apolar permeability / water solubility),
hydrophobic SASA.

### Demoted to "diagnostics" — compute, tabulate, do NOT build narrative on
Per the review, these have weak/absent primary validation as permeability predictors:

- **Amphipathic / integy moment** — *no primary experimental validation found* for macrocycle permeability.
- **Weighted RMSF** — ensemble *diagnostic*, not a predictor.
- **Kier Φ** — *applicability-domain filter* (Φ<10 = "is this method reliable here"), not a predictor.
- **Asphericity / NPR shape** — scaffold-contingent.
- **SA_HD / SA_HA** — mechanistically sound but rarely quantified as standalone predictors.

> These stay in `phys_descriptors_v3.py` and the descriptor CSV. They move out of
> the report's interpretive sections into a single compact diagnostics table.

---

## 3. Two-state apolar re-evaluation — a CREST **feature matrix**, not one ΔG

> **UPDATE 2026-07-02 — un-deferred, and simplified. This supersedes the deferred status
> and the "chloroform + re-opt" recipe below.** Decisions this session:
> - **ΔG_transfer is the active permeability descriptor**, computed in-house at the cheap
>   tier — the FlexiSol benchmark (Grimme, *Chem. Sci.* 2025, [[reference_flexisol_grimme_2025]])
>   shows partition ratios are *more* robust than absolute solvation energies (cross-solvent
>   error cancellation), and that lowest-E ≈ full-Boltzmann while a *random/cherry-picked*
>   single conformer fails — so weight by energy, never by a hand-picked geometry (that
>   retires the earlier min-PSA idea).
> - **Apolar phase = cyclohexane, not hexadecane.** Verified from Ono 2019 (our SASA source
>   paper): explicit cyclohexane SASA correlated *excellently* with permeability, chloroform
>   *weaker*. Hexadecane appears nowhere in this literature; use hexane/hexadecane only as a
>   *documented surrogate* if xTB lacks a cyclohexane keyword. Chloroform → NMR/IMHB anchor only.
> - **Native per-phase, never retrofit.** Run CREST separately per phase (water + cyclohexane);
>   do NOT re-score/re-opt a chloroform pool into cyclohexane. FlexiSol: phase-specific
>   geometries matter most for flexible IMHB molecules (= our chameleons).
> - **Score with CPCM-X (native to xTB — no ORCA/COSMO-RS).** `xtb --cpcmx <solvent>`, no
>   gradient → single-point only (which is what we want; CREST/ALPB already did the geometry).
>   ALPB kept as a one-time comparison. GFN2 electronic energy is sufficient (FlexiSol: level
>   of theory barely moves partition ratios).
> - **No CENSO/COSMO-RS ceiling.** Goal is qualitative sorting + relative trends; the ceiling
>   is experimental data / the collaborator hand-off.
> - **Implementation:** `scripts/free_energy_calculator.py` (phase-specific legs, CPCM-X
>   default, Boltzmann weights + ΔG_transfer) + `scripts/slurm_free_energy.sh`. Runs on the
>   **HPC** (xtb is Linux-only). Open check: confirm `cyclohexane` + CPCM-X support in the
>   cluster's xtb build (`xtb --version` ≥ 6.6; `xtb --help`).
>
> _Historical (2026-06-23) design retained below for the CREMP ML feature-matrix track._

> **Status (2026-06-23): deferred / not built for the current hits.** A collaborating
> computational chemist will run higher-fidelity (explicit-solvent) simulations for our
> hits, so we will **not** force an implicit-ALPB ΔG_transfer that we'd only have to
> caveat. The in-house per-molecule reports therefore do **not** include a ΔG_transfer /
> cross-solvent feature. The design below is retained as the recipe for the *CREMP ML
> track* (§4), where a learned feature matrix is the point — not for these select hits.

**Reframe (important).** Implicit-continuum CREST/xTB energies are *not* a real-world
ΔG_transfer — the physics gap (no explicit solvent, approximate entropy) makes the raw
ΔE unreliable as a calibrated free energy. So we do **not** treat raw ΔE as ΔG_transfer.
Instead we treat CREST's structural + electronic outputs as a **feature matrix** and let
the ML regressor (Track B) learn the correction to experimental permeability. This still
implements the strongest missing *signal* without a second full conformer search.

1. **Foundational ensemble = chloroform** (already generated by Tier-2 CREST). Gives
   the diverse library of accessible intermediate shapes.
2. **Apolar re-evaluation:** take the top-N lowest-energy chloroform conformers
   (~top 20), and for each run xTB in water, chloroform, and an apolar core solvent
   (`--alpb {cyclohexane|hexadecane}`). **Quick `--opt` (not bare `--sp`)** in the
   apolar solvent so the geometry is allowed to collapse — otherwise ΔSASA between
   chloroform and "cyclohexane" is degenerate (same coordinates → same surface).

**Feature matrix to extract per molecule (Boltzmann-weighted over top-N):**
- **⟨ΔE_solv / Gsolv⟩** across water vs chloroform vs apolar (the chameleon energy axis;
  a *feature*, not a calibrated ΔG). Parse `Gsolv` from the xTB ALPB summary.
- **ΔSASA** between the exposed water conformer and the collapsed apolar conformer
  (RDKit `rdFreeSASA` on each geometry). Also Δ(3D-PSA), ΔRg.
- **Dipole-moment fluctuation** across environments (xTB prints the molecular dipole;
  chameleonic molecules suppress their dipole in the apolar core). Use spread across solvents.
- **HOMO–LUMO gap** per solvent (parse from xTB output) — electronic / polarizability info.
- Optional RDKit shape deltas, no new dependency: Δ`CalcWHIM`, Δ`CalcMORSE`.

**Caveats (state in methods, do not overclaim):**
- Implicit ALPB, *not* Kamenik-grade GIST. These are **ML features**, not ΔG predictions.
- **Verify cyclohexane is in the ALPB set** on the cluster (`xtb --alpb`). If not, use
  **hexadecane (ε≈2.06)** — a standard, supported bilayer-core surrogate the review lists
  alongside cyclohexane.

---

## 4. ML track — finalized

**The constraint the review makes explicit:** the best literature correlations are
N=8 (Ono), N=6 (Kamenik), N=7 (Begnini). Our own compound set is a handful of R/S
pairs — **too small to train an ML model.** Continuing to add descriptors to our
own data polishes a set that will never be ML-trainable.

**Two tracks, previously conflated:**

- **Track A — our compounds (mechanistic case study).** Validate the
  permeability–solubility trade-off across R/S epimers using the validated core.
  Small N, descriptor *validation* only, drives the per-molecule reports.
- **Track B — the ML project (CREMP).** Train on CREMP's **3,258
  permeability-annotated macrocycles** (31.3M conformers, generated with our exact
  CREST/GFN2-xTB protocol, in chloroform). Compute the validated core (+ ΔG_transfer)
  on *their* ensembles. Liu 2025 benchmark: 2D representations fail (no 3D burial);
  graph models (DMPNN) win → our 3D descriptors are the value-add as features or a
  physics-informed baseline.

**Next ML actions (separate work items):**
- Pull CREMP permeability-annotated subset; confirm conformer/format compatibility.
- Run `phys_descriptors_v3.py` (validated core + ΔG_transfer) over the CREMP subset.
- Baseline: regression of permeability on the 4-descriptor core; compare to a DMPNN
  graph baseline per Liu 2025.

---

## 4b. Publication section template — validate-then-extend (Begnini 2021)

**Standing precedent for the descriptor section of the paper** (Begnini et al., *ACS Med.
Chem. Lett.* **2021**, *12* (6), 983–990). How that paper argues it:

1. **Validate on a small set with NMR/NAMFIS** (their compounds 1 vs 2). The CREST-sampled
   ensemble does *not* reproduce every full-molecule conformer, but it reproduces the
   **macrocyclic core** well, and the two **validated descriptors (Rgyr, SA 3D-PSA)** agree
   between the NMR-derived and the computational ensembles — and that descriptor difference
   matches the measured permeability difference (their Fig 4).
2. **Extend the purely computational ranking prospectively** to a test set (their 5–7) that
   has no NMR; the sampled descriptors again track permeability.

**Mapping to us:** CREST is the engine; NMR/NAMFIS is only affordable for a few hits. For a
representative subset, show (a) small RMSD of the macrocyclic core (CREST vs NMR) and (b)
agreement on Rgyr / SA 3D-PSA. That agreement **licenses extending** the computational R/S
ranking to the remaining analogs/linkers without NMR. This is the structure to build that
section around as the publication develops.

**Current report status:** deliberately **omitted from the per-molecule reports for now**
(decision 2026-06-23) — we do not yet have subset NMR for the 3-12-x-12 hits (the in-progress
NMR effort is CsA/reference validation). The 3-12-8-12 / 3-12-10-12 reports therefore stay
purely computational with no forward reference to NMR; upgrade to past-tense validation only
once subset NMR lands.

---

## 5. Action items

- [x] `make_molecule_report.py`: keep all descriptors in the tables, but mark the
      weak/unvalidated ones (SA_HD/HA, amphipathic moment, asphericity, RMSF) with † and
      a note that their permeability correlation is currently low — retained for ML, not
      driving the call. (Done 2026-06-23; reports need regeneration to pick it up.)
- [x] Generated publication 3D figures (Fig 4 S-ensemble, Fig 5 R-ensemble, Fig 6 R/S
      min-energy overlay; water | chloroform panels) for both molecules via
      `scripts/make_isomer_3d_views.py` (PyMOL ray-traced, 2400 px, per-solvent intra-fit
      fans, PNG + SVG). Done 2026-06-23.
- [x] Regenerated both per-molecule reports (3-12-8-12, 3-12-10-12) + docx — Figs 4/5/6 wired
      in, † weak-descriptor marks + low-correlation note present, 6 images embedded per docx.
- [~] ~~`compute_dg_transfer.py`~~ — **deferred to collaborator** (explicit-solvent sims). Not built for current hits.
- [ ] (Track B) CREMP subset already in repo (`results/archive/cremp_deltapsa.csv`,
      `scripts/cremp_deltapsa.py`, `notebooks/cremp_benchmark.ipynb`) — descriptor run is the ML work.

---

*Sources mapped: Ono 2019 (JCIM, SASA cyclohexane R²=0.872), Begnini 2021 (SA 3D-PSA,
Rg, isomer ranking), Kamenik 2020 (ΔG_transfer r=0.92), Rezai 2006 (IMHB ~100-fold),
Whitty 2016 / Caron 2024 (ΔPSA chameleonicity), Liu 2025 (2D-ML failure, DMPNN),
Grambow 2024 (CREMP). Full DOIs in the review PDF, References §F.*

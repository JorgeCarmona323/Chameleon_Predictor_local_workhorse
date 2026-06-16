# 3D Descriptor Literature Review — what to add, what to keep, what to benchmark

**2026-06-13 · informs `scripts/phys_descriptors_v3.py`**

> **Decision (this doc records it):** We are **not** trimming descriptors preemptively on
> literature grounds. We **add** the literature-validated descriptors we are currently
> missing, **keep** everything we already compute (including the collinear shape pair),
> and let **in-house ML benchmarking on our own compound set** decide what to drop. The
> literature tells us what is *worth measuring*; our data tells us what is *predictive for
> our chemistry*. Those are different questions and we have only answered the first one.

---

## Why this review

We have never benchmarked descriptors against permeability on our own set. Before we do,
we wanted to confirm (a) which of our current descriptors are literature-validated, and
(b) which validated descriptors we are *missing* and could compute from the 3D ensembles
we already generate (no new simulation needed). Seven papers from the
`3d_descriptors literature/` folder were reviewed.

## Papers reviewed

| # | Source | Scope | Takeaway for us |
|---|--------|-------|-----------------|
| 1 | **Rzepiela 2022** (*J. Med. Chem.*) | ~3,600 macrocycles, MD-derived 3D descriptors vs PAMPA | **SA_HBD** (solvent-accessible surface of H-bond donors) was the **single most important** descriptor. Also surfaces ASA_H (hydrophobic surface area). |
| 2 | **Kim 2025** | Conformational descriptors for membrane permeability | **Radius of gyration ranked #1**; 3D-PSA and IMHB strong. All three we already compute. |
| 3 | **García Jiménez 2024** | Integy / amphipathic moments | **Spatial segregation** of polar vs nonpolar surface (the "integy moment") tracks permeability — a *distribution* descriptor, not a *count* descriptor. We have nothing like it. |
| 4 | **Poongavanam 2020** | Chameleonicity / solvent-dependent 3D-PSA | Validates our solvent-resolved ΔPSA / Δ(descriptor) cross-solvent design. |
| 5 | **Wicker & Cooper** (nConf20) | Conformational flexibility from RDKit conformer counts | **Flexibility** (how many accessible conformers) is itself predictive. We capture ensemble shape via `p_dominant` but have no explicit flexibility count. |
| 6 | **Severoglu 2025** | 3D descriptor survey | Corroborates shape + surface families; no new must-have. |
| 7 | **Sugita 2025** | Enhanced-sampling permeability modeling | Corroborates solvent-resolved ensemble approach; no new must-have. |

## What we already have (literature-validated — keep)

| Our descriptor | Validated by | Status |
|---|---|---|
| 3D-PSA (`bw_psa`) | Kim 2025, Poongavanam 2020 | keep |
| Intramolecular H-bonds (`bw_hb`) | Kim 2025 (IMHB) | keep |
| Radius of gyration (`bw_rg`) | **Kim 2025 (#1)** | keep |
| NPR1/NPR2, asphericity, spherocity | Severoglu 2025, Kim 2025 | keep all (see note) |
| cis-amide propensity (`cis_prob_*`) | scaffold-specific, our own | keep |
| Dominant-conformer population (`p_dominant`) | flexibility proxy | keep |
| Cross-solvent Δ (`delta_*`, `norm_delta_psa`) | Poongavanam 2020 | keep |

**Note on the shape pair.** Spherocity and asphericity (and NPR1/NPR2, eccentricity) are
all derived from the **same principal moments of inertia**, so they are mathematically
collinear — confirmed. The original instinct to drop one is *defensible*, but we are
**keeping both**: collinearity hurts coefficient interpretability in a linear model, not
predictive accuracy in a tree/ensemble model, and dropping a feature before we have a
single benchmark number is premature. The benchmark's feature-importance / correlation
matrix is the right place to make that call, on our data.

## What we are missing (literature-validated — ADD in v3)

| New descriptor | Definition | Why | Source |
|---|---|---|---|
| **SA_HBD** (`*_bw_hbd_sasa`) | Solvent-accessible surface area restricted to H-bond **donor** atoms (polar H on N/O, plus their heavy atom) | "Single most important" over 3,600 macrocycles. Distinct from PSA: PSA counts *all* polar heavy-atom exposure; SA_HBD isolates the donors that actually pay the desolvation penalty crossing the membrane. | Rzepiela 2022 |
| **Hydrophobic SASA** (`*_bw_hydrophobic_sasa`) | SASA of apolar atoms (C and H-on-C) | The nonpolar surface that *favors* membrane partitioning — the complement of PSA, and not recoverable from PSA alone once a molecule folds. | Rzepiela 2022 (ASA_H) |
| **Amphipathic moment** (`*_bw_amphi_moment`) | Distance (Å) between the SASA-weighted centroid of polar surface and that of apolar surface | Captures *spatial segregation* of polar vs nonpolar surface — a chameleon can have fixed PSA but reorganize *where* the polar surface points. A pure count/area descriptor is blind to this. | García Jiménez 2024 (integy moment) |
| **Effective conformer count** (`*_n_eff`) | exp(Shannon entropy of Boltzmann weights) = effective # of populated conformers | In-ensemble analog of nConf20 flexibility. `p_dominant` says how peaked the top state is; `n_eff` says how many states are really populated. | Wicker & Cooper (nConf20, analog) |

All four are computable from the **ensembles we already have** — no new CREST runs.

## Explicitly *not* adding (yet)

- **Witek congruent population / kinetic barrier** — prior 2026-05-31 decision; needs
  clustering + barrier estimation we have deferred. Unchanged.
- **3D-CCS / collision cross section** — needs a CCS calculator; revisit if we get IM-MS data.
- **Literal nConf20** — would require a *separate* RDKit conformer-generation count with a
  fixed energy window; our `n_eff` reads the CREST ensemble we already trust instead.

## How this feeds the benchmark

`phys_descriptors_v3.py` is an **additive superset** of v2: every v2 function is preserved
and re-exported, and the four new per-conformer surface/shape descriptors are added.
`ensemble_descriptors.py` Boltzmann-weights them into new CSV columns alongside the
existing ones. The benchmark then trains on the **full** column set and reports
feature importance + a correlation matrix; *that* output — not this literature review —
decides the final descriptor set for our system.

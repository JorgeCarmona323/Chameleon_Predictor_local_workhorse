# Experiment — CsA 3D Descriptors: Experimental Structures vs CREST V1 Ensemble

**Date:** 2026-06-05
**Script:** `scripts/ensemble_descriptors.py` (CREST side) + ad-hoc gemmi loader (experimental side)
**Data:** `data/CREST_CsA_20260512/` (CREST V1 water, 23 conformers) vs `data/experimental_structure_references_CsA/*.cif`

---

## Question

Does the CREST V1 water ensemble of cyclosporin A reproduce the 3D descriptors of the experimentally determined aqueous conformer (A1), and where does it diverge?

The A1 conformer (Bhatt et al., *JACS* 2022; solution NMR ROE + X-ray/neutron crystal) is the aqueous-relevant form, defined by a **cis amide at MeVal11–MeBmt1**. CREST V1 was run *without* the `-notopo` flag, which we suspected suppressed cis/trans interconversion during sampling.

---

## Method

- **CREST V1:** Boltzmann-weighted descriptors over the 23-conformer water ensemble, plus the dominant conformer (46.3% population). Descriptors via `ensemble_descriptors.py`.
- **Experimental:** three CSD crystal structures loaded with `gemmi.read_small_structure` (fractional → Cartesian):
  - A1 X-ray (CCDC 2149649) — aqueous conformer
  - A1 neutron (CCDC 2149650) — aqueous conformer, all H resolved
  - C1 DEKSAN (CCDC 1138505, 1985) — closed/apolar conformer
- Coordinate-only descriptors (Rg, intramolecular H-bonds) computed directly from element + coordinates.
- Crystal waters removed by keeping the largest covalently-connected fragment (= CsA, 196 atoms) before H-bond counting.

---

## Results

### Atom counts revealed explicit cavity waters
CsA = C62H111N11O12 = **196 atoms**. CREST V1 and the closed C1 crystal both have exactly 196. The aqueous A1 structures have **more**:
- A1 X-ray: 199 (+3 ≈ 1 water)
- A1 neutron: 206 (+10 ≈ 3 waters, neutron resolves all H)

These extras are **water molecules in/around the macrocycle cavity** — part of what stabilizes the open aqueous conformer. CREST uses **implicit ALPB solvent**, so it has no explicit waters and cannot form these specific water bridges.

### Descriptor comparison (crystal waters stripped; 3D PSA via rdFreeSASA)

| Structure | 3D PSA (Å²) | intramolecular HB | Rg (Å) | cis MeVal11–MeBmt1 |
|---|---|---|---|---|
| A1 aqueous (X-ray) | **137.5** | ~2 | 6.15 | **cis** |
| C1 closed (DEKSAN) | **95.9** | 4 | 6.42 | trans |
| CREST V1 dominant (46%) | 67.9 | 3 | 6.27 | **trans** |
| CREST V1 Boltzmann | 84.1 | 2.31 | 6.34 | **trans (cis_prob = 0, all 11 bonds)** |
| CREST V1 min-PSA conf | 51.2 | — | — | trans |
| CREST V1 max-PSA conf | 146.2 | — | — | trans |

- **PSA (now computed properly):** A1 aqueous exposes the most polar surface (137.5) — open, water-loving. C1 closed buries it (95.9) — membrane form. **CREST V1 is over-collapsed:** its dominant conformer (67.9) is *more polar-buried than even the closed crystal* (95.9), and the ensemble average (84.1) sits far below the aqueous A1 (137.5). In water, CREST V1 behaves like a low-dielectric environment.
  - *(Experimental PSA computed via a minimal element-only RDKit mol + rdFreeSASA polar query — same method as the pipeline; bonds not needed since SASA uses element radii + coordinates.)*
- **H-bonds:** raw experimental counts (4–5) were inflated by CsA–water bridges. Intramolecular-only, A1 has ~2 — and CREST V1's 2.31 matches it.
- **Rg:** CREST V1 (6.27–6.34) sits between aqueous A1 (6.15) and closed C1 (6.42), skewing toward the closed fold.
- **cis-amide (the headline):** CREST V1 is **100% trans on all 11 backbone amides** (`cis_prob = 0`). A1 has a **cis MeVal11–MeBmt1**. CREST V1 categorically misses it.

### Does anything close to A1 exist in the CREST V1 ensemble?

**Partly by openness, not at all by the defining geometry:**
- **By PSA/openness:** yes — the max-PSA conformer (146.2) is *more* exposed than A1 (137.5). Open conformers do exist in the ensemble.
- **By cis-amide:** **no** — every open conformer is all-trans. They are "open trans" forms, not the cis A1.
- **By population:** the open conformers are rare; the dominant/populated states are hyper-collapsed (PSA 68). The ensemble is dominated by closed-like trans structures.

So CREST V1 samples the right *openness range* but lands in the wrong *basins*: it never adopts A1's cis backbone, and its populated states are collapsed. The dominant conformer is **not** A1 — by Rg and PSA it most resembles a hyper-closed version of the C1 conformer.

---

## Conclusion

Two failures, both pointing to the same fixes:

1. **Wrong backbone geometry:** CREST V1 is 100% trans on all 11 amides; the aqueous A1 is defined by cis MeVal11–MeBmt1. Without `-notopo`, CREST could not cross the cis/trans barrier. No conformer in the ensemble adopts A1's cis geometry.

2. **Over-collapsed in water:** CREST V1's populated conformers bury *more* polar surface (dominant PSA 67.9, ensemble 84.1) than even the closed crystal form (95.9), and far more than the open aqueous A1 (137.5). The ALPB implicit solvent — with no explicit cavity waters to scaffold the open conformer — lets the molecule collapse in "water" as if it were in a low-dielectric medium.

The dominant CREST V1 conformer is **not** A1; by Rg and PSA it most resembles a hyper-closed C1. Open conformers (PSA up to 146) exist but are rare and still all-trans.

This is the quantitative justification for **CsA_v2 (`--noreftopo -notopo`)** — we expect non-zero `cis_prob` at MeVal11–MeBmt1 and a higher, A1-like water PSA — and flags **implicit solvation** as a second-order limitation that may need explicit-water (Tier-2 OpenMM) sampling to fully reproduce the open aqueous conformer.

**Takeaway for the descriptor pipeline:** `cis_prob` catches the backbone-geometry failure, and 3D PSA catches the over-collapse — two complementary descriptors, both earning their place in the feature set. Bulk HB/Rg looked fine; PSA + cis-amide exposed the real problems.

---

## Limitations

1. **count_hbonds_xyz** uses a fixed geometric cutoff (H···A < 2.5 Å, D–H···A > 120°), which differs from crystallographic H-bond assignment — absolute counts are method-dependent (comparison is internally consistent).
2. **Water stripping** via covalent connected-components is imperfect: A1 X-ray stripped to 195/196 (one atom short), A1 neutron to 200/196 (over-kept, HB unreliable → neutron HB discounted). X-ray strip trusted.
3. **PSA method note** — the analytic SASA fallback returns 0 for 200-atom molecules; experimental PSA therefore computed via rdFreeSASA on a minimal element-only mol (radii + coordinates, no bond perception needed). This matches the pipeline's PSA method, so the comparison is apples-to-apples.
4. **Single structure vs thermal ensemble** — comparing one crystal snapshot to a Boltzmann-averaged ensemble.
5. **Implicit vs explicit solvent** — CREST ALPB cannot reproduce the explicit cavity waters present in the A1 crystal, a known source of divergence for the aqueous conformer.

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

### Descriptor comparison (crystal waters stripped)

| Structure | intramolecular HB | Rg (Å) | cis MeVal11–MeBmt1 |
|---|---|---|---|
| A1 aqueous (X-ray) | ~2 | 6.15 | **cis** |
| C1 closed (DEKSAN) | 4 | 6.42 | trans |
| CREST V1 dominant | 3 | 6.27 | **trans** |
| CREST V1 Boltzmann | 2.31 | 6.34 | **trans (cis_prob = 0, all 11 bonds)** |

- **H-bonds:** the raw experimental counts (4–5) were inflated by CsA–water bridges. Intramolecular-only, the aqueous A1 has ~2 H-bonds — and **CREST V1's 2.31 matches it well.** CREST V1 is *not* under-H-bonded relative to the aqueous conformer.
- **Rg:** CREST V1 (6.27–6.34) sits between aqueous A1 (6.15) and closed C1 (6.42), skewing slightly toward the closed fold.
- **PSA:** CREST V1 `bw_psa = 84.07 Å²` (pipeline rdFreeSASA, reliable). Experimental PSA not computed — the analytic SASA fallback returns 0 for large molecules; would require rdFreeSASA with perceived connectivity (deferred).
- **cis-amide (the headline):** CREST V1 is **100% trans on all 11 backbone amides** (`cis_prob = 0` everywhere). The experimental A1 has a **cis MeVal11–MeBmt1**. CREST V1 categorically misses it.

---

## Conclusion

CREST V1's **bulk** 3D descriptors (size, intramolecular H-bond count) are actually a reasonable match to the experimental aqueous A1 conformer. The single, decisive divergence is the **cis MeVal11–MeBmt1 amide** — the geometric feature that *defines* the A1 state. Without `-notopo`, CREST could not cross the cis/trans barrier and sampled an all-trans, slightly-too-closed ensemble.

This is the quantitative justification for the **CsA_v2 rerun with `--noreftopo -notopo`**: we expect a non-zero `cis_prob` at MeVal11–MeBmt1 in water, and a `cis_switch_bond` between water and CHCl3.

**Takeaway for the descriptor pipeline:** bulk descriptors can look correct while the mechanistically-critical backbone geometry is wrong. The `cis_prob` descriptor is the one that catches it — confirming it belongs in the feature set.

---

## Limitations

1. **count_hbonds_xyz** uses a fixed geometric cutoff (H···A < 2.5 Å, D–H···A > 120°), which differs from crystallographic H-bond assignment — absolute counts are method-dependent (comparison is internally consistent).
2. **Water stripping** via covalent connected-components is imperfect: A1 X-ray stripped to 195/196 (one atom short), A1 neutron to 200/196 (over-kept, HB unreliable → neutron HB discounted). X-ray strip trusted.
3. **PSA not compared** — analytic SASA fallback unusable for 200-atom molecules; needs rdFreeSASA + connectivity.
4. **Single structure vs thermal ensemble** — comparing one crystal snapshot to a Boltzmann-averaged ensemble.
5. **Implicit vs explicit solvent** — CREST ALPB cannot reproduce the explicit cavity waters present in the A1 crystal, a known source of divergence for the aqueous conformer.

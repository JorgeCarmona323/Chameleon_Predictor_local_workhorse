# ML Descriptor Implications from Limbach 2025 + Witek 2016

## Limbach et al. J. Med. Chem. 2025 — Biased Equilibrium / Goldilocks

Permeability requires a *biased* equilibrium between aqueous conformer (A, cis MeVal11-MeBmt1) and membrane conformer (C, cis MeLeu9-MeLeu10).

- CycH is impermeable despite having both conformers because ΔG(O1↔C1) ≈ 2.35 kcal/mol (Keq ≈ 1, balanced equilibrium)
- CycA is permeable because ΔG(A1↔C1) ≈ 13.4 kcal/mol (strongly biased toward A in water, toward C in membrane)
- Goldilocks principle: barrier must be high enough to maintain long-lived states but not trap the molecule in the membrane
- β-hydroxyl of MeBmt1 is structurally required for the A conformer; analogs without it lose permeability

**Key new descriptor:** ΔE(dominant_water − dominant_CHCl3) = proxy for equilibrium bias. Large → biased (permeable). Small → balanced (impermeable even if both conformers exist).

## Witek et al. J. Chem. Inf. Model. 2016 — Congruent Conformations

"Congruent conformations" = metastable states significantly populated in BOTH solvents. C1/W4 (closed pair) and C2/W1 (half-open pair) are the congruent states for CsA.

- Permeability correlates with population of congruent states in water
- Single-start MD cannot cross cis-trans barrier (confirmed at 100 ns at 300K, even 400K insufficient)
- Must seed from both CRYSTC (closed) and CRYSTO (open) structures — exact justification for multi-start CREST (v3.3)

**Key new descriptor:** `congruent_pop_water` = Boltzmann fraction of water-ensemble conformers that have an RMSD match (< threshold) in the CHCl3 ensemble.

## ML Implications

- With <100 compounds: classical ML (XGBoost, RF) with 3D descriptors is the right approach, not large model training
- Must add: ΔE between dominant conformers across solvents, cis-amide fraction in water, cis-amide count/type per conformer
- PSA/HB remain useful but cannot distinguish CycH-type false positives (permeable-shaped but balanced equilibrium)

## Descriptor Scope Decision (2026-05-31) — ML focus, not SAR/theory

Goal is **predictive features for the ML model**, not mechanistic completeness. Decisions:

- **Limbach (ΔΔG biased equilibrium) > Witek (congruent states).** Limbach subsumes Witek: CycH has congruent-like states but balanced equilibrium (ΔG≈2.35 vs CycA's 13.4 kcal/mol) → impermeable. Witek's `congruent_frac` would flag CycH as a false positive. **→ Drop the Witek congruent descriptor.**
- **Kinetic barrier (Goldilocks) deferred.** CREST gives thermodynamics only — no transition states, no real barriers. The Eyring-from-ΔG proxy isn't a true barrier. **→ Only revisit if model performance demonstrably needs it.** Not in the first feature set.
- **No RMSD clustering needed for the baseline.** We ran two separate CREST simulations (water-dielectric, CHCl3-dielectric); the solvent label *is* the state assignment. Descriptors = Boltzmann-average each solvent ensemble, take water−CHCl3 differences. Clustering was only required for congruent population (now dropped) and within-solvent multi-basin splitting (defer).

**Locked feature set (computable now from `ensemble.json` + `ensemble.sdf`, both solvents):**
- Per solvent: `bw_psa`, `bw_hb`, `bw_rg`, `bw_npr1/2`, `bw_asphericity`, `bw_spherocity`, `psa_std`, `psa_spread`, `sasa_total`, `p_dominant`, `cis_prob_i` (per amide bond)
- Cross-solvent: `delta_psa`, `norm_delta_psa` (Yu 2026), `delta_hb`, `delta_rg`, `delta_npr1/2`, `ddG` (from dominant populations), `cis_switch_bond`, `delta_cis_prob_i`, `cis_entropy`

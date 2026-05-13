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
- Must add: ΔE between dominant conformers across solvents, cis-amide fraction in water, congruent conformer population (water), cis-amide count/type per conformer
- PSA/HB remain useful but cannot distinguish CycH-type false positives (permeable-shaped but balanced equilibrium)
- Congruent state computation requires cross-ensemble RMSD clustering — only possible once multi-start CREST generates correct basin coverage for both solvents

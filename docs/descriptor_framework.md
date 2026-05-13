# Generalized 3D Descriptor Framework for Permeability ML Model

Based on Limbach 2025 (J. Med. Chem.) but generalized away from CsA-specific positions. Core concept: find the two-state biased equilibrium (State W in water, State M in CHCl3) wherever it exists in any cyclic peptide scaffold — no hardcoded residue positions.

## Key Design Decisions

- Compute cis-amide propensity for **all** backbone amide bonds (N-Me, regular, Pro) — not just N-Me. Regular amide propensity will be ~0 unless ring strain forces it; no prefiltering needed.
- Cis amides beyond N-methylation: proline (~1–2 kcal/mol cis-trans gap), and ring-strained macrocycles can force non-NMe bonds toward cis. Pipeline must handle these.
- Do NOT hardcode MeBmt1 β-OH — search all protonated HBDs (OH, NH) for the one most exposed in State W/water and buried in State M/CHCl3. The β-OH analog may be serine, threonine, or another side chain in different scaffolds.
- Macrostates defined data-driven via joint RMSD clustering of water + CHCl3 ensembles, not by fixed residue positions.

## Descriptor Groups

### 1. Macrostate Definition
Joint RMSD clustering of water + CHCl3 ensembles. State W = dominant water cluster, State M = dominant CHCl3 cluster. Characterized post-hoc by which cis bonds and H-bond motifs co-occur. Requires MDAnalysis or MDTraj + scipy Ward linkage.

### 2. Cis-Amide Propensity (all amide bonds)
- `cis_prob_i(solvent)` = Boltzmann-weighted fraction of conformers where bond i has ω < 30°
- `Δcis_prob_i` = cis_prob_i(water) − cis_prob_i(CHCl3)
- `cis_switch_bond` = argmax(|Δcis_prob_i|) — which bond most drives the two-state switch
- `cis_entropy(solvent)` = −Σ p_i log(p_i) over all amide bonds; low = one dominant switching bond

### 3. Macrostate Bias
- `p(W, solvent)`, `p(M, solvent)`, `ΔG(W−M, solvent)` per solvent
- `ΔΔG` = ΔG(W−M, water) − ΔG(W−M, CHCl3) — **core bias index**; large → permeable (CycA-like), ~0 → impermeable (CycH-like)
- `Has_bias` = 1 if |ΔΔG| > ~2 kcal/mol

### 4. Exterior HBD Search (Generalized β-OH)
For each protonated HBD k (OH, NH — not NMe):
- `SASA_k` per conformer
- `internal_HB_k` = 1 if forming intramolecular H-bond (distance + angle cutoff)

Aggregated:
- `max_ΔSASA_HBD` = max over k of [⟨SASA_k⟩_W,water − ⟨SASA_k⟩_M,CHCl3]
- Finds the HBD most exposed in State W/water and buried in State M/CHCl3 regardless of which residue it's on

### 5. Boltzmann Polarity per Macrostate
- `⟨PSA_exposed⟩_W`, `⟨PSA_exposed⟩_M` per solvent — 3D SASA-based on polar atoms (not 2D TPSA)
- `⟨nHBD_exposed⟩_W`, `⟨nHBD_exposed⟩_M` per solvent
- `ΔPSA(W−M, water)`

### 6. Shape
- `⟨Rg⟩_W`, `⟨Rg⟩_M`, `ΔRg(W−M)` — pure numpy from coordinates
- `Anisotropy_W`, `Anisotropy_M` — inertia tensor diagonalization (principal moment ratio)
- `⟨CCS⟩_W`, `⟨CCS⟩_M`, `ΔCCS` — optional, requires MOBCAL or ML-CCS tool

### 7. Kinetic Proxies (optional, approximate)
- Eyring τ(W↔M) per solvent using ΔG as barrier proxy (no actual transition state — rough approximation)
- `log τ_ratio` = log[τ(water)/τ(CHCl3)]

## Implementation Status

| Descriptor | Status | Notes |
|---|---|---|
| Cis-amide propensity | Ready | 3D coords + RDKit → omega dihedral per conformer |
| Boltzmann populations | Ready | ensemble.json has weights |
| ΔG per solvent | Ready (per solvent) | ΔΔG blocked until CHCl3 ensemble complete |
| Internal H-bond detection | Ready | Distance + angle cutoff from 3D coords |
| Rg | Ready | Numpy from coordinates |
| Inertia tensor / anisotropy | Ready | Linear algebra on 3D coords |
| 3D PSA_exposed (SASA-based) | Gap | Need FreeSASA or MSMS |
| Exterior HBD SASA | Gap | Same SASA requirement |
| Joint RMSD clustering | Gap | Need MDAnalysis + scipy |
| CCS prediction | Missing | Completely separate tool |
| Kinetic proxies | Approximation | No actual barrier heights |
| Congruent conformer population | Gap | Cross-ensemble RMSD matching not implemented |

## Key Limitations

- **Ensemble completeness**: single-start CREST misses cis-containing basins; v3.3 multi-start is the fix
- **xTB accuracy**: 1 kcal/mol error at 300K (RT ≈ 0.6 kcal/mol) can shift a 20% population conformer to 5% or 50%
- **ALPB solvation**: implicit solvent misses specific solute-solvent interactions; SASA of HBDs is a proxy for explicit H-bonding
- **Kinetic descriptors**: Eyring from ΔG is useful for relative rankings only, not absolute rates
- **Dataset size**: <100 compounds → use regularized models (RF with depth limits, XGBoost with early stopping); feature selection matters more than model choice

## Implementation Order
1. Add SASA calculation (FreeSASA) — unblocks PSA_exposed and exterior HBD features
2. Implement joint RMSD clustering — unblocks macrostate assignment
3. Wait for CHCl3 ensembles — unblocks ΔΔG
4. CCS tool — optional, add last

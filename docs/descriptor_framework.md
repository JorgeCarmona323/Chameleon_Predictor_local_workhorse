# Generalized 3D Descriptor Framework for Permeability ML Model

Based on Limbach 2025 (J. Med. Chem.) but generalized away from CsA-specific positions. Core concept: find the two-state biased equilibrium (State W in water, State M in CHCl3) wherever it exists in any cyclic peptide scaffold — no hardcoded residue positions.

*Last updated: 2026-05-31 — SASA gap closed (rdFreeSASA wired in `phys_descriptors_v2.py`); pipeline integration added.*

---

## Where This Fits in the Permeability Pipeline

These descriptors are the **featurization for the `DynamicEnsembleEncoder`** of `CyclicPermeabilityModel` (see `chameleon_model_architecture.md`). They are computed per-compound from the solvent-derived CREST ensembles (water + CHCl3) and written to the `dynamic_features` table (see `data_schema.md`).

**Training-data constraint (locked):** only solvent-derived ensembles feed the model. Our CREST/GFN2-xTB + ALPB ensembles (water → `environment_id=3`, CHCl3 → `environment_id=2`) qualify; vacuum ETKDGv3 conformers do **not**. This framework operates strictly on CREST output.

**Two implementation tiers:**
- **Tier 1 (Phase 1, no clustering):** `aq` = max-PSA conformer, `mem` = min-PSA conformer; descriptors are whole-ensemble Boltzmann averages and spreads. Maps directly to the existing `dynamic_features` schema columns. Buildable today.
- **Tier 2 (Phase 2, macrostate-based):** State W / State M defined by joint RMSD clustering; ΔΔG and congruent population computed across macrostates. Requires the clustering gap (group 1) to be closed.

**Size-gating caveat (hypothesis.md + Yu 2026):** chameleonic switching is gated at ~9 residues. For sub-threshold scaffolds (4–6mers, e.g. the DOPC/Brain hits), absolute ΔPSA is a poor descriptor — `norm_delta_psa` (ΔPSA / SASA_total) and shape/cis-amide descriptors are expected to carry more signal. Always compute the normalized form alongside the absolute.

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
- `norm_delta_psa` = ΔPSA / SASA_total — dimensionless fractional switching ratio (Yu 2026); the size-robust form that should outperform absolute ΔPSA, especially below the 9-residue threshold

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
| Cis-amide propensity | Ready | 3D coords + RDKit → omega dihedral per conformer (not yet coded) |
| Boltzmann populations | Ready | `ensemble.json` has weights; `boltzmann_weights()` in `phys_descriptors_v2.py` |
| ΔG per solvent | Ready | ΔΔG now unblocked for compounds with both water + CHCl3 ensembles (6-mers done) |
| Internal H-bond detection | **Done** | `count_hbonds_xyz()` in `phys_descriptors_v2.py` |
| Rg / shape | Ready | `Descriptors3D` in `conformer_engine.py`; not yet wired to CREST ensembles |
| Inertia tensor / anisotropy | Ready | Linear algebra on 3D coords |
| 3D PSA_exposed (SASA-based) | **Done** | `compute_psa_xyz()` (rdFreeSASA, Bondi radii) in `phys_descriptors_v2.py` |
| Exterior HBD SASA | Partial | rdFreeSASA available; per-HBD atom loop not yet written |
| Joint RMSD clustering | Gap | Need scipy Ward linkage on heavy-atom RMSD matrix |
| norm_delta_psa (Yu 2026) | Gap | ΔPSA / SASA_total — trivial once total SASA added (rdFreeSASA no query) |
| CCS prediction | Missing | Completely separate tool (MOBCAL / ML-CCS) — optional |
| Kinetic proxies | Approximation | No actual barrier heights |
| Congruent conformer population | Gap | Cross-ensemble RMSD matching; needs both ensembles + clustering |

## Key Limitations

- **Ensemble completeness**: single-start CREST misses cis-containing basins; v3.3 multi-start is the fix
- **xTB accuracy**: 1 kcal/mol error at 300K (RT ≈ 0.6 kcal/mol) can shift a 20% population conformer to 5% or 50%
- **ALPB solvation**: implicit solvent misses specific solute-solvent interactions; SASA of HBDs is a proxy for explicit H-bonding
- **Kinetic descriptors**: Eyring from ΔG is useful for relative rankings only, not absolute rates
- **Dataset size**: <100 compounds → use regularized models (RF with depth limits, XGBoost with early stopping); feature selection matters more than model choice

## Implementation Order (revised 2026-05-31)
1. ~~Add SASA calculation~~ — **done** (`compute_psa_xyz` rdFreeSASA)
2. ~~Wait for CHCl3 ensembles~~ — **done** for 6-mers (WhC3, DOPC R/S, Brain1, DOPC2); CsA_v2 11-mer pending
3. **Now buildable (Tier 1):** wire shape descriptors + cis-amide omega + Boltzmann PSA/HB into a single ensemble-descriptor script reading `ensemble.json` + `ensemble.sdf`; emit `dynamic_features`-schema rows for both solvents incl. `norm_delta_psa`
4. Implement joint RMSD clustering (scipy Ward) — unblocks Tier 2 macrostate ΔΔG + congruent population
5. CCS tool — optional, add last

## Target Descriptor Output (Tier 1 — maps to `dynamic_features` schema)

Per compound, per solvent (water, CHCl3), Boltzmann-weighted over the CREST ensemble:

| Column | Source group | Notes |
|---|---|---|
| `bw_psa3d` | 5 | Boltzmann-mean polar SASA |
| `aq_psa3d` / `mem_psa3d` | 5 | max-/min-PSA conformer (Tier 1 proxy for State W/M) |
| `psa3d_std`, `psa3d_spread` | 5 | ensemble polarity flexibility |
| `delta_psa3d`, `norm_delta_psa` | 5 | absolute + size-normalized switching |
| `bw_hb`, `aq_hb_count`, `mem_hb_count`, `delta_hb` | 4 | intramolecular H-bonds |
| `bw_rg`, `aq_rg`, `mem_rg`, `delta_rg` | 6 | radius of gyration |
| `aq_npr1/2`, `mem_npr1/2`, `asphericity`, `spherocity` | 6 | shape anisotropy |
| `cis_prob_i` (per bond), `cis_switch_bond`, `cis_entropy` | 2 | cis-amide propensity |
| `dG_W_M` (per solvent), `ddG`, `has_bias` | 3 | macrostate bias (Tier 1: dominant-conformer proxy) |

Cross-solvent (computed once both ensembles present): `ΔΔG`, `Δcis_prob_i`, `congruent_pop_water` (Tier 2).

# Math and Methods Background for Descriptor Development

Reference material for implementing and QC-ing the 3D descriptor pipeline.

## Priority Math Topics

1. **Statistical mechanics / Boltzmann statistics** — partition functions, ΔG = −RT ln(p1/p2), sensitivity to energy errors at 300K (RT ≈ 0.6 kcal/mol). A 1 kcal/mol xTB error can shift a 20% population conformer to 5% or 50%.

2. **Linear algebra** — inertia tensor construction and diagonalization (for Rg and anisotropy descriptors), Kabsch algorithm (RMSD-optimal structural alignment, needed for clustering), PCA for ensemble dimensionality reduction.

3. **Geometric calculations** — dihedral angle formula using atan2 and cross products (4-atom omega angle), RMSD between structures, H-bond distance and angle criteria.

4. **Clustering** — hierarchical (Ward linkage) and k-means, silhouette score for choosing cutoff, RMSD-based clustering specifically (linkage between geometry and cluster assignment).

5. **Information theory** — Shannon entropy (−Σ p log p) for the cis-bond diversity descriptor. Low entropy = one dominant switching bond; high = disordered.

6. **Thermodynamic cycles** — how ΔΔG relates to relative permeability across environments, Eyring equation and its assumptions (when it's valid as a barrier approximation and when it isn't).

## Reference Book

**Frenkel & Smit — "Understanding Molecular Simulation"** covers all of the above at the right level for computational chemistry implementation.

## People to Consult

| Person | Institution | Expertise | Accessibility |
|---|---|---|---|
| PI (Hu lab) | SDSU | Scope and priorities | First stop |
| Mike Gilson | UCSD Skaggs School of Pharmacy | Free energy methods, implicit solvation | On campus, directly relevant |
| Rommie Amaro | UCSD | MD simulations, enhanced sampling, drug discovery | On campus |
| Stefan Riniker | ETH Zürich | CsA conformational sampling, cyclic peptides | Cold email viable with validation data |
| Vittorio Limbach | First author Limbach 2025 | Biased equilibrium, CsA analogs | Reachable if CsA validation results are in hand |

## Practical Implementation Resources

- **MDAnalysis** — RMSD clustering, trajectory analysis, SASA
- **MDTraj** — fast trajectory analysis, dihedral calculation
- **FreeSASA** — standalone SASA calculation (fills the gap in current pipeline)
- MDAnalysis and MDTraj community forums for implementation questions without needing a full collaboration

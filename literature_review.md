# Literature Review and Workflow Validation Report

**Project:** 3D Conformational Descriptors and Dual-Dielectric Solvent Modeling to Decode Cyclic Peptide Membrane Permeation
**Course:** CHEM 269 Final Project | Jorge Carmona | March 2026

---

## 1. Chameleonic Behavior — Well-supported

The chameleonic behavior hypothesis is robustly supported. Rezai et al. (2006, *J. Am. Chem. Soc.* 128, 14073) established the mechanistic picture: CycloA and N-methylated cyclic hexapeptides shield polar NH groups through intramolecular H-bonds in apolar environments, substantially reducing effective PSA. Naylor et al. (2018, *J. Med. Chem.* 61, 11169) showed H-bond donor count in the apolar conformer (NMR in CDCl3) strongly correlates with passive permeability — more so than 2D TPSA. Witek et al. (2016, *J. Chem. Theory Comput.* 12, 4398) benchmarked MD-derived ensembles in explicit water vs. chloroform for CycloA analogs, validating the dual-solvent ensemble approach.

**Assessment:** No serious flaws. Max-PSA/min-PSA conformer selection is a computationally tractable simplification of a well-validated concept.

---

## 2. CycPeptMPDB + PAMPA — Appropriate

CycPeptMPDB (Jiang et al., 2023, *J. Chem. Inf. Model.* 63, 2240) is the largest publicly available cyclic peptide permeability dataset (~8,466 entries, ~7,298 PAMPA). PAMPA measures passive transcellular diffusion and is appropriate for chameleonic behavior studies. The database's precomputed CHCl3_3DPSA and H2O_3DPSA make it ideal for validating Tier-1 approximations.

**PAMPA limitations to acknowledge:**
- No active transport (P-gp efflux not captured)
- Assay heterogeneity across contributing labs (different lipid compositions, pH, time points)
- Solubility artifacts at high concentrations for lipophilic peptides

---

## 3. ETKDGv3 for Macrocycles — Appropriate, with caveats

ETKDGv3 (Wang, Riniker, Landrum, *J. Chem. Inf. Model.* 2020, 60, 2044) includes ring-template macrocycle sampling and outperforms ETKDGv2 for 8–15-atom rings by coverage of crystallographic conformations.

**Limitations:**
- 50 conformers may undersample 10–12-residue cyclic peptides (200–500 recommended for rigorous work; 50 is acceptable for Tier-1 extremes)
- Vacuum geometry — does not model solution ensemble
- N-methylated residue MMFF94s parameter coverage can be incomplete

---

## 4. MMFF94s Dual-Dielectric Approximation — Valid approximation; must be caveated

MMFF94s minimization is in vacuum. Max-PSA/min-PSA selection identifies structural extremes, not thermodynamic populations in each dielectric. The correspondence holds approximately (high-PSA ≈ water-favorable) but can fail for molecules with competing intramolecular H-bond networks. Witek et al. showed GB/SA-reweighted ensembles correlate better with experiment.

**Assessment:** Valid as Tier-1; prominently caveat that Δ descriptors are from vacuum-ensemble extremes, not equilibrium populations. Tier-2 (GB/SA cross-check vs. DB values) appropriately addresses this.

---

## 5. UMAP + Leiden (kNN graph) — Methodologically correct

UMAP (McInnes et al., *arXiv:1802.03426*, 2018) with cosine metric is established in cheminformatics (Probst & Reymond, *J. Cheminform.* 2020, 12, 12). Leiden clustering (Traag et al., *Sci. Rep.* 2019, 9, 5233) on the kNN graph — **not** the 2D embedding — is the correct approach. This pipeline gets this right and is a notable methodological strength.

---

## 6. 3D-PSA vs. 2D-TPSA — Strongly supported

3D-PSA on apolar conformers correlates better with passive permeability than 2D TPSA for cyclic peptides (Ertl et al. 2D-TPSA: *J. Med. Chem.* 2000, 43, 3714 was developed for linear small molecules). The CycPeptMPDB paper (Jiang et al. 2023) explicitly benchmarks CHCl3_3DPSA vs. simpler descriptors.

---

## 7. Workflow Limitations to Flag in Write-up

1. **Thermodynamic vs. structural extremes:** Max/min-PSA conformer ≠ Boltzmann-weighted ensemble mean
2. **50-conformer undersampling:** Flag as Tier-1 limitation for large macrocycles
3. **PAMPA assay heterogeneity:** Use assay-source as a covariate or restrict to PAMPA-only entries
4. **Spearman r preferred over Pearson r** for permeability data spanning log units (non-linear)
5. **ΔNPR (shape descriptors):** Weaker direct literature support — treat as hypothesis-generating

---

## Key References

| Citation | Relevance |
|----------|-----------|
| Rezai et al., *JACS* 2006, 128, 14073 | Chameleonic behavior, CsA conformational flexibility |
| Naylor et al., *J. Med. Chem.* 2018, 61, 11169 | HB donor count + permeability |
| Witek et al., *JCTC* 2016, 12, 4398 | MD ensemble in water vs. CHCl3 |
| Bockus et al., *J. Med. Chem.* 2015, 58, 4581 | Cyclic peptide permeability SAR |
| Jiang et al., *JCIM* 2023, 63, 2240 | CycPeptMPDB |
| Wang, Riniker, Landrum, *JCIM* 2020, 60, 2044 | ETKDGv3 macrocycle benchmarking |
| Riniker, Landrum, *JCIM* 2015, 55, 2562 | ETKDG original |
| Ertl et al., *J. Med. Chem.* 2000, 43, 3714 | TPSA |
| Veber et al., *J. Med. Chem.* 2002, 45, 2615 | Oral bioavailability rules |
| McInnes et al., *arXiv:1802.03426* 2018 | UMAP |
| Traag et al., *Sci. Rep.* 2019, 9, 5233 | Leiden clustering |
| Probst & Reymond, *J. Cheminform.* 2020, 12, 12 | UMAP for chemical space |

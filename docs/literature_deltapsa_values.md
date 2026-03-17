# Literature ΔPSA Values — Chameleonic Cyclic Peptides

**Compiled:** 2026-03-17
**Purpose:** Reference values for validating Tier-1 ETKDGv3+MMFF94s ΔPSA against published literature

---

## Critical Note: Three Incompatible ΔPSA Definitions

| Definition | Formula | Notes |
|------------|---------|-------|
| **A — Whitty/Doak (static)** | TPSA − MPSAnp | 2D TPSA minus 3D PSA of single lowest-energy nonpolar conformer |
| **B — Ensemble (your Tier-1)** | PSA(water ensemble) − PSA(nonpolar ensemble) | Boltzmann-averaged or min/max over conformational ensemble. Physically meaningful. |
| **C — Yu 2026 (normalized ratio)** | PSA(H₂O)/SASA(H₂O) − PSA(CHCl₃)/SASA(CHCl₃) | Normalized fraction, not absolute Å². Not directly comparable. |

**Your Tier-1 uses Definition B.** CsA ~75–80 Å² (Witek 2016) is the correct Definition B comparator. CsA 174 Å² (Doak 2016) is Definition A — do not mix.

---

## Tier 1: Experimentally Proven (NMR + PAMPA/Caco-2 + PSA documented)

### 1. Cyclosporin A (CsA)
**NMR evidence:** Kessler et al. (1985) established the closed CDCl₃ conformer with 4 intramolecular H-bonds and antiparallel β-sheet at residues 7–11. Witek et al. (2016) used ROESY in CHX-d₁₂ and HEX-d₁₄ with 57–79 NOEs per solvent confirming closed conformer in all apolar media and a distinct open conformer in aqueous/DMSO. Rüdisser et al. (2023) achieved backbone RMSD 0.10 Å using exact NOE (eNOE) restraints in CDCl₃/hexadecane vs. MeOH/H₂O — the most precise experimental structures to date. PDB deposited.

**ΔPSA:** ~75–80 Å² (Definition B, Witek 2016); 174 Å² (Definition A, Doak 2016)
**Permeability:** PAMPA logPe = **−5.01** (Ahlbach et al. 2015)

**Citations:**
- Kessler H et al. *Helv. Chim. Acta* 1985, 68:661. DOI: 10.1002/hlca.19850680318
- Witek J et al. *J. Chem. Inf. Model.* 2016, 56:1547. DOI: 10.1021/acs.jcim.6b00251
- Rüdisser SH et al. *J. Am. Chem. Soc.* 2023, 145:27601. DOI: 10.1021/jacs.3c09367
- Doak BC et al. *Drug Discov. Today* 2016, 21:1389. DOI: 10.1016/j.drudis.2016.02.014
- Ahlbach CL et al. *Future Med. Chem.* 2015, 7:2121. DOI: 10.4155/fmc.15.78

---

### 2. Roxithromycin, Telithromycin, Spiramycin (Macrolide Antibiotics)
**NMR evidence:** Danelius et al. (2020) obtained NOE distances and vicinal scalar couplings in CDCl₃, DMSO-d₆, and DMSO-d₆/D₂O (10:1). All three macrolides populated significantly less polar, more compact conformational ensembles in CDCl₃ vs. DMSO/water — experimentally proven chameleonic switching. Rifampicin showed less pronounced switching (useful negative control). Rossi Sebastiano et al. (2018) corroborated via crystal structure analysis — minimum 3D PSA in nonpolar environments strongly correlated with efflux-inhibited Caco-2 permeability.

**ΔPSA:** Exact Å² values in paywalled tables; minimum nonpolar 3D PSA ≤140 Å² documented
**Permeability:** Passive efflux-inhibited Caco-2 Papp (Rossi Sebastiano 2018); roxithromycin oral bioavailability ~72–85%

**Citations:**
- Danelius E et al. *Chem. Eur. J.* 2020, 26:5231. DOI: 10.1002/chem.201905599
- Rossi Sebastiano M et al. *J. Med. Chem.* 2018, 61:4189. DOI: 10.1021/acs.jmedchem.8b00347

---

### 3. Cyclo-Leu Diastereomers (Rezai series)
**NMR evidence:** Rezai et al. (2006) characterized 11 cyclic peptides by NMR H/D exchange in CDCl₃. Permeable diastereomers showed 4/5 backbone amides H-bonded in CDCl₃; impermeable ones showed only 2/5. Directly experimental proof of differential intramolecular H-bonding driving permeability.

**ΔPSA:** Computed from low/high-ε conformers anchored to NMR structures; R² = 0.96 vs. PAMPA for 11 compounds
**Permeability:** PAMPA logPe spanning ~−4 to −8 for the series

**Citations:**
- Rezai T et al. *J. Am. Chem. Soc.* 2006, 128:2510. DOI: 10.1021/ja0563455
- Rezai T et al. *J. Am. Chem. Soc.* 2006, 128:14073. DOI: 10.1021/ja063076p

---

### 4. Cyclic Hexapeptide Diastereomers (Ono et al. 2019)
**NMR evidence:** Eight diastereomers (MW 712.93, TPSA 186 Å²) studied by Δδ/ΔT NMR in CDCl₃. Compound 3: 1 exposed NH; compounds 7 and 8: 0 exposed NH — directly observable intramolecular H-bonding by temperature-coefficient NMR.

**ΔPSA:** ~**32 Å²** (cyclohexane ~140 Å² vs. water ~172 Å²; from 1.40 vs. 1.72 nm² in paper)
**Permeability:** PAMPA 1.66–7.74 × 10⁻⁶ cm/s; MDCK 0.4–19.3 × 10⁻⁶ cm/s

**Citation:**
- Ono S et al. *J. Chem. Inf. Model.* 2019, 59:2952. DOI: 10.1021/acs.jcim.9b00217. PMC: PMC7751304

---

## Tier 2: Strong Evidence (NMR structural characterization; permeability discussed)

### 5. Omphalotin A (Rüdisser et al. 2023) — ⭐ Compute PSA yourself
**NMR evidence:** Most technically rigorous modern example. Exact NOEs (eNOE) in CDCl₃/hexadecane-d₃₄ → two slow-exchange conformers C1 and C2 with defined H-bonds. In CD₃OH/H₂O → "indole-in" and "indole-out" states with no stable H-bonds. Backbone RMSD 0.49 Å. **PDB structures deposited: 8QAQ (apolar C1), 8QAS (polar)** — you can compute PSA_apolar and PSA_polar directly in RDKit.

**ΔPSA:** Not reported in paper but computable from PDB 8QAQ vs 8QAS
**Permeability:** Not measured in this paper (mechanistic study)

**Citation:**
- Rüdisser SH et al. *J. Am. Chem. Soc.* 2023, 145:27601. DOI: 10.1021/jacs.3c09367. PDB: 8QAQ, 8QAS

---

### 6. Sanguinamide A Analogs — Lokey Lab
**NMR evidence:** VT-¹H NMR in DMSO-d₆ and CDCl₃ for 17 analogs. Parent SA-A1: NOT chameleonic (no backbone change). SA-B1, SA-B3: backbone H-bond network changes between solvents. N-methylated analogs showed solvent-dependent flexibility.

**ΔPSA:** ~8% PSA reduction for SA-B3 in CDCl₃; ~0 for parent SA-A1 (negative control)
**Permeability:** PAMPA + Caco-2 for all 17 analogs (paywalled values)

**Citation:**
- Bockus AT et al. *J. Med. Chem.* 2015, 58:7409. DOI: 10.1021/acs.jmedchem.5b00919

---

### 7. CsA, Alisporivir, Cyclosporin H (Taming Conformational Heterogeneity)
**NMR evidence:** 2D ROESY in CD₃OD/H₂O and CD₃CN. CsA and Alisporivir (both permeable) shift from open aqueous to closed apolar conformer. CsH (impermeable) does NOT access the closed conformer — direct NMR proof that chameleonic switching is mechanistically required for permeability in the cyclosporin family.

**Citation:**
- *JACS Au* 2024, 4:1458. DOI: 10.1021/jacsau.4c00011. PMC: PMC11040698

---

## Summary Table

| Molecule | NMR Proven | ΔPSA (Å²) | Permeable | Best Reference |
|----------|-----------|-----------|-----------|---------------|
| CsA | Yes (Kessler 1985, Witek 2016, Rüdisser 2023) | ~75–80 (B) / 174 (A) | Yes (PAMPA −5.01) | Witek 2016, Doak 2016 |
| Roxithromycin | Yes (Danelius 2020) | paywalled | Yes | Danelius 2020 |
| Telithromycin | Yes (Danelius 2020) | paywalled | Yes | Danelius 2020 |
| Spiramycin | Yes (Danelius 2020) | paywalled | Yes | Danelius 2020 |
| Rifampicin | Yes (Danelius 2020) | small (negative ctrl) | Moderate | Danelius 2020 |
| Cyclo-Leu diastereomers | Yes (Rezai 2006) | computed ~32–80 | Mixed | Rezai 2006 JACS |
| Cyclic hexapeptide 1–8 | Yes (Ono 2019) | ~32 | Mixed | Ono 2019 JCIM |
| Omphalotin A | Yes (Rüdisser 2023) | computable (PDB 8QAQ/8QAS) | Not measured | Rüdisser 2023 JACS |
| Sanguinamide A (parent) | Yes — negative ctrl | ~0 | Moderate | Bockus 2015 |
| Alisporivir | Yes (2024 JACS Au) | similar to CsA | Yes | JACS Au 2024 |
| Cyclosporin H | Yes — negative ctrl | low (no switch) | No | JACS Au 2024 |

---

## Your Tier-1 vs Literature

| Compound | Tier-1 ΔPSA (Å²) | Literature ΔPSA | Notes |
|----------|-----------------|-----------------|-------|
| CsA | 84.9 | ~75–80 (Witek 2016) | ✓ Within ~10% |
| DP172 | 88.9 | None found | Strongest individual case |
| HexPep | 64.4 | None found | Non-permeable despite high ΔPSA — likely <9-residue effect (Yu 2026) |
| 1NMe3 | 47.8 | None found | Permeable with lower ΔPSA — <9-residue effect |
| PSLYF | 65.3 | None found | Non-permeable |

**Important:** Yu et al. 2026 (*bioRxiv* DOI: 10.64898/2026.01.06.697862) finds ΔPSA is unreliable for <9-residue peptides. HexPep, 1NMe3, PSLYF are hexapeptides (6 residues) — the counterintuitive ΔPSA values may reflect this limitation.

---

## Action Item: Compute Omphalotin A PSA from PDB

```python
# Download 8QAQ (apolar) and 8QAS (polar) from RCSB
# Load with RDKit and compute 3D SASA
# ΔPSA = PSA(8QAS) - PSA(8QAQ)
# This would give you an experimentally grounded ΔPSA for a second molecule
```

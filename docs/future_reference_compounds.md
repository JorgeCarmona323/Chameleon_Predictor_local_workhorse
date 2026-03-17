# Future Reference Compound Validation

**Priority:** High — strengthens the experimental grounding of the ΔPSA validation story
**Status:** Pending computation / journal access

---

## What We Have Now

| Molecule | Tier-1 ΔPSA | Literature ΔPSA | Status |
|----------|------------|-----------------|--------|
| CsA | 84.9 Å² | ~75–80 Å² (Witek 2016) | ✓ Validated |

---

## Priority Calculations Needed

### 1. Omphalotin A — from PDB 8QAQ / 8QAS
- PDB already downloaded to `/tmp/8QAQ.pdb` (apolar) and `/tmp/8QAS.pdb` (polar)
- Need to fix RDKit PSA computation from PDB block (polar atom radii issue)
- Once computed: run Tier-1 ETKDGv3 on Omphalotin A SMILES and compare
- Reference: Rüdisser et al. *JACS* 2023, 145:27601. DOI: 10.1021/jacs.3c09367
- **This would give a second NMR-grounded ΔPSA benchmark**

### 2. Cyclic Hexapeptide Diastereomers (Ono 2019)
- ΔPSA ~32 Å² already published (cyclohexane ~140 Å² vs water ~172 Å²)
- Get SMILES for compounds 1–8 from Ono et al. *J. Chem. Inf. Model.* 2019, 59:2952
- Run Tier-1 on each and compare to NMR-derived values
- PAMPA values available in paper for direct validation
- **Best multi-compound validation set available**

### 3. Cyclosporin H (CsH) — Negative Control
- NMR proven does NOT undergo chameleonic switch (JACS Au 2024, DOI: 10.1021/jacsau.4c00011)
- Should show LOW Tier-1 ΔPSA — prediction: <30 Å²
- Check if in CycPeptMPDB; if not, compute from SMILES
- **Key negative control: same scaffold as CsA, no switch, not permeable**

### 4. Roxithromycin, Telithromycin, Spiramycin
- NMR chameleonic switching proven (Danelius 2020, DOI: 10.1002/chem.201905599)
- Exact PSA values paywalled — need journal access or compute from SMILES
- Passive Caco-2 Papp available in Rossi Sebastiano 2018
- **Three-compound positive control set**

### 5. Rifampicin — Negative Control
- Less chameleonic than macrolides in Danelius 2020
- Lower Caco-2 permeability
- Run Tier-1 and expect lower ΔPSA than roxithromycin
- **Negative control for macrolide comparison**

### 6. Sanguinamide A (parent) — Negative Control
- VT-NMR shows ~0 PSA change between solvents (Bockus 2015, DOI: 10.1021/acs.jmedchem.5b00919)
- PAMPA measured in paper
- Run Tier-1 and expect near-zero ΔPSA
- **Negative control within the cyclic peptide class**

---

## Proposed Final Reference Set

Once computations are done, the reference set becomes:

| Molecule | Evidence | Role |
|----------|----------|------|
| CsA | NMR + PAMPA + Tier-1 validated | Positive control anchor |
| Omphalotin A | eNOE PDB + Tier-1 | Second positive control |
| CsH | NMR proven no switch | Negative control (same scaffold) |
| Cyclic hexapeptide 1–8 (Ono) | NMR + PAMPA, R²=0.96 | Multi-compound validation |
| Sanguinamide A | VT-NMR ~0 ΔPSA | Negative control (cyclic peptide) |
| Roxithromycin | NMR + Caco-2 | Positive control (macrolide) |

---

## Technical Notes

- PDB → PSA issue: `rdFreeSASA.CalcSASA` returns total float, not per-atom. Need to zero non-polar radii approach or use `MakeFreeSasaPolarAtomQuery` with correct atom typing from PDB
- All SMILES-based calculations can use the existing `conformer_engine.py` directly
- CsH SMILES: available from ChemDraw/PubChem — differs from CsA at position 2 (L-Abu → D-Ala) and lacks some N-methylation

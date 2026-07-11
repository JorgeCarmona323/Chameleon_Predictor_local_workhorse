# Cyclosporin A — Aqueous Conformer Validation

**Project:** Chameleon Predictor  
**Compound:** Cyclosporin A (CsA) — CycPeptMPDB ID 1  
**Date started:** 2026-05-08  
**Reference:** Limbach et al., *J. Am. Chem. Soc.* 2022, 144, 12602–12607  

---

## Objective

Validate that the GFN2-xTB/CREST water-phase conformer ensemble captures the experimentally known aqueous conformer (A1) of CsA, as characterized by X-ray/neutron diffraction and NMR by Limbach et al. 2022.

---

## Background

### Why CsA?

CsA is an 11-residue *N*-methylated macrocyclic peptide and a prototypical "chameleon" — it switches conformations between polar and apolar environments, directly enabling passive membrane permeability despite high molecular weight (~1202 Da).

### The A1 conformer (Limbach et al. 2022)

The paper reports the first atomic-resolution aqueous conformer (A1) via single-crystal neutron diffraction + NMR titration in methanol/water mixtures.

| Feature | Description |
|---|---|
| Cis-amide | MeVal11−MeBmt1 (residues 11→1) |
| H-bond 1 | Abu2 (NH) ··· MeLeu10 (C=O) |
| H-bond 2 | Val5 (NH) ··· Ala7 (C=O) |
| Cavity waters | Two water ligands inside macrocycle |
| Crystal system | Orthorhombic P2₁2₁2₁ (CCDC 2149649) |
| NMR signature | Hα(MeVal11)−Hα(MeBmt1) ROE cross-peak |

> **Key finding:** A1 is the dominant conformer at 90% water (v/v). The previously assumed "closed" (CHCl₃) and "open" (CypA-bound) forms do not explain aqueous behavior — A1 does.

---

## Computational Setup

### Method

- **Conformer sampling:** CREST iMTD-GC, GFN2-xTB, ALPB solvation (water)
- **Pre-optimization:** GFN2-xTB/ALPB(water)
- **Refinement:** `--cregen` on existing rotamers (69,986 rotamers from prior MTD run)
- **Cluster:** SDSU GPU cluster, 20 CPUs, SLURM job 259118

### CsA sequence (IUPAC)

```
cyclo(MeBmt1−Abu2−Sar3−MeLeu4−Val5−MeLeu6−Ala7−D-Ala8−MeLeu9−MeLeu10−MeVal11)
```

Residue numbering follows standard CsA convention (MeBmt1 = (4R)-4-[(E)-2-butenyl]-4,N-dimethyl-L-threonine).

---

## Results — Water Phase

### Cregen run summary

| Metric | Value |
|---|---|
| Input rotamers | 69,986 |
| Unique conformers | **23** |
| Runtime | 1 min 24 s |
| Dominant conformer population | **46.3%** (Boltzmann) |
| Job status | ✅ Complete |

> 23 unique conformers is consistent with published literature for CsA in polar solvents (Ono et al. *J. Chem. Inf. Model.* 2021 reported similar ensemble size).

### Ensemble file

```
results/crest_runs/run_20260503_150449_1_CsA/water/crest/crest_conformers.xyz
```

---

## Validation Checklist

### Structural features to verify in ensemble

- [ ] **Cis-amide at MeVal11−MeBmt1** — omega dihedral ≈ 0° (±30°) in dominant conformer(s)
- [ ] **H-bond Abu2···MeLeu10** — N···O distance ≤ 3.5 Å, N−H···O angle ≥ 120°
- [ ] **H-bond Val5···Ala7** — N···O distance ≤ 3.5 Å, N−H···O angle ≥ 120°
- [ ] **All other amide bonds are trans** — omega ≈ 180° for residues 2–10
- [ ] **Boltzmann population of A1-like conformers** — expect dominant fraction ≥ 30–50%

### NMR-derivable checks (qualitative)

- [ ] Short Hα(MeVal11)−Hα(MeBmt1) distance (< 4 Å) in cis-amide conformers → would give ROE cross-peak
- [ ] Abu2 NH and Val5 NH **not** solvent-exposed in H-bonded conformers (intramolecular HB shields from water)

---

## Validation Script Plan

Script to write: `scripts/validate_csa_water.py`

**Inputs:**
- `crest_conformers.xyz` (23 conformers, multi-model XYZ)
- `crest_cregen.out` (energies for Boltzmann weighting)
- CsA SMILES / RDKit mol (for atom indexing)

**Outputs:**
- Table: conformer | energy | Boltzmann pop | cis/trans at 11→1 | HB1 dist | HB2 dist
- Summary: fraction of population with A1 features
- PNG: energy vs. omega scatter plot (color-coded cis/trans)

**VMD visualization:**
- Multi-frame DCD or multi-model PDB from ensemble
- TCL script to color by omega(MeVal11−MeBmt1): red = cis, blue = trans
- Highlight H-bond donor/acceptor pairs

---

## CHCl₃ Phase (Pending)

| Status | Detail |
|---|---|
| Job | SLURM 259118 |
| Expected runtime | ~15–20 h from session start |
| Expected result | "Closed" conformer (4 intramolecular HBs, cis MeLeu9−MeLeu10) |

> The CHCl₃ ensemble should capture the **closed** conformer (C1) — the opposite chameleon state. Comparing water vs. CHCl₃ ensembles side-by-side is the core of the chameleon story for the PI.

---

## Next Steps

- [ ] Check job 259118 status (`squeue -u $USER`)
- [ ] Transfer water `crest_conformers.xyz` to local (`scp`)
- [ ] Run `validate_csa_water.py` — check cis-amide + H-bonds
- [ ] Once mem job finishes: transfer CHCl₃ ensemble
- [ ] Compare water vs CHCl₃: PSA, HB count, SASA
- [ ] Build VMD visualization for PI presentation
- [ ] Resubmit DP955 (compound 3) and DP944 (compound 4) — no checkpoints, need full reruns

---

## References

1. Limbach, M. N. et al. Atomic View of Aqueous Cyclosporine A: Unpacking a Decades-Old Mystery. *J. Am. Chem. Soc.* **2022**, *144*, 12602–12607. https://doi.org/10.1021/jacs.2c01743  
2. Ono, S. et al. Cyclosporin A: Conformational Complexity and Chameleonicity. *J. Chem. Inf. Model.* **2021**, *61*, 5601–5613.  
3. Witek, J. et al. Kinetic Models of Cyclosporin A in Polar and Apolar Environments. *J. Chem. Inf. Model.* **2016**, *56*, 1547–1562.  

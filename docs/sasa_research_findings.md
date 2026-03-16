# SASA Computation: State-of-the-Art Research Findings

**Project:** CHEM 269 Final — Cyclic Peptide Dual-Dielectric Pipeline
**Author:** Jorge Carmona
**Date:** March 2026
**Purpose:** Literature + GitHub survey of how polar SASA / 3D PSA is computed in practice. Used to inform the `_polar_sasa` fix in `scripts/conformer_engine.py`.

---

## 1. Standard Algorithm

Both **Shrake-Rupley** (1973) and **Lee-Richards** (1971) are standard and give nearly identical results for small molecules. Both are implemented in FreeSASA / rdFreeSASA via `SASAAlgorithm.ShrakeRupley` and `SASAAlgorithm.LeeRichards`.

**Probe radius: 1.4 Å (0.14 nm) is universal.** This represents a water molecule and is used without variation across the entire drug discovery literature.

---

## 2. Polar Atom Definition

The consensus in the cyclic peptide / macrocycle literature (Witek 2016, Riniker group, Ertl 2000) is:

> **Polar atoms = N, O, and H atoms covalently bonded to N or O**

This is identical to the definition used by RDKit's 2D TPSA descriptor (`rdMolDescriptors.CalcTPSA`). Some implementations additionally include S and P (the OONS classifier, RDKit TPSA with `includeSandP=True`).

| Convention | Polar atoms |
|---|---|
| Strict N/O (Witek 2016, Ertl 2000 TPSA default) | N, O, polar-H |
| Broad (OONS classifier) | All non-carbon atoms: N, O, S, P, halogens |
| Protein-residue-based (Protor, NACCESS) | Residue/atom-name dependent — **fails on small molecules** |

---

## 3. Critical rdFreeSASA Bug (Confirmed, GitHub #1827 and #5197)

**`rdFreeSASA.classifyAtoms()` silently fails for small molecules built from SMILES.**

Root cause: The Protor, NACCESS, and OONS classifiers inside FreeSASA identify atoms by **PDB residue name and atom name** (e.g., "ALA", "CA"). For arbitrary small molecules assembled from SMILES, no `monomerInfo` (residue/chain/atom name metadata) is present. The function returns either all-zero radii or `None`, and sets every atom to `SASAClass.Unclassified`.

Consequence: `MakeFreeSasaPolarAtomQuery()` builds a query on `SASAClass = Polar`. If no atoms are classified as Polar (all Unclassified), the query matches nothing → `CalcSASA(..., query=polar_query)` returns **0.0 for every conformer of every molecule**.

This is a **silent failure** — no exception is raised, no warning is printed.

### Confirmed workaround (community-established)

Manually assign `SASAClass` and `SASAClassName` atom properties before calling `CalcSASA`, bypassing `classifyAtoms()` entirely:

```python
for atom in mol.GetAtoms():
    if atom.GetSymbol() in {'N', 'O', 'S', 'P'}:
        atom.SetIntProp('SASAClass', 0)       # 0 = Polar
        atom.SetProp('SASAClassName', 'Polar')
    else:
        atom.SetIntProp('SASAClass', 1)       # 1 = APolar
        atom.SetProp('SASAClassName', 'APolar')
```

Also: use **Bondi VdW radii** manually (not from `classifyAtoms()`), as `classifyAtoms()` may return zeros:

| Element | Bondi radius (Å) |
|---|---|
| H | 1.20 |
| C | 1.70 |
| N | 1.55 |
| O | 1.52 |
| S | 1.80 |
| P | 1.80 |
| F | 1.47 |
| Cl | 1.75 |
| Br | 1.85 |
| I | 1.98 |
| Default | 1.50 |

---

## 4. Witek et al. 2016 Method (ΔPSA ~75 Å² for CsA)

The landmark Riniker group papers (Witek 2016 *JCTC*, related PMC7388155, PMC7751304) use:

- **Software**: GROMACS `gmx sasa`
- **Algorithm**: Shrake-Rupley / rolling probe
- **Probe radius**: 1.4 Å
- **Polar atom selection**: N, O, and H atoms bonded to N or O (explicit hydrogen model)
- **Protocol**: Full MD simulations in **explicit** water and **explicit** chloroform; 5,000+ frames from equilibrated trajectories; PSA computed per frame and averaged (Boltzmann-weighted)
- **ΔPSA**: PSA(water ensemble mean) − PSA(chloroform ensemble mean)

**The ~75 Å² ΔPSA for CsA is an ensemble-averaged result from explicit-solvent MD, not a single conformer.** Our Tier-1 (ETKDG + MMFF) approximates this by using max-PSA vs min-PSA across a vacuum conformer ensemble. Tier-2 CREST+ALPB is closer to the Witek approach (environment-specific ensembles) but still uses implicit solvent.

### Reproducibility in Python

Exact reproduction requires GROMACS + force field parameterization. Reasonable approximation:
1. ETKDG conformer ensemble (Tier-1) or CREST+ALPB (Tier-2)
2. Manual Bondi radii + SASAClass assignment
3. `rdFreeSASA.CalcSASA(..., query=MakeFreeSasaPolarAtomQuery())`
4. Select max-PSA conformer (≈ water form) and min-PSA conformer (≈ membrane form)
5. ΔPSA = PSA_max − PSA_min

---

## 5. Alternative Tools

| Tool | Polar SASA? | Small molecule support | Notes |
|---|---|---|---|
| `rdFreeSASA` (RDKit) | Yes (with manual SASAClass fix) | Yes (after fix) | Fastest for in-pipeline use |
| `freesasa` Python package | Yes (custom Classifier subclass) | Yes | Cleaner API; requires structure building from coordinates |
| MDTraj `shrake_rupley()` | Manual post-filter by element | Yes (from trajectory) | Best for MD trajectories; used by Riniker group |
| `gmx sasa` (GROMACS) | Yes (selection by atom name) | With parameterization | Ground truth for Witek-style ΔPSA |
| OpenMM / AMBER GB/SA | Total SASA only | Yes | No built-in polar decomposition |

---

## 6. Applied Fix in This Pipeline

**File:** `scripts/conformer_engine.py`, function `_polar_sasa()`

**Fix applied:** Bypass `classifyAtoms()` entirely. Manually assign Bondi radii and `SASAClass=Polar` for N, O, S, P atoms. Use `MakeFreeSasaPolarAtomQuery()` as intended. All atoms retain their radii for mutual occlusion — only polar atoms are summed. Polar H excluded (only heavy atom contributions; consistent with Ertl 2000 heavy-atom TPSA convention).

**Note on polar-H inclusion:** Witek 2016 includes polar H in their PSA calculation. Ertl 2000 TPSA uses a fragment-based approach that implicitly accounts for NH/OH groups. Our implementation excludes explicit polar-H from the PSA sum (heavy atoms only). This means our absolute PSA values will be slightly lower than Witek's MD values, but the **relative ΔPSA ranking across conformers** is unaffected. This is acceptable for Tier-1 comparative analysis.

---

## 7. Key References

| Source | Key finding |
|---|---|
| Witek et al., *JCTC* **2016** | ΔPSA ~75 Å² for CsA from explicit-solvent MD; N/O/polar-H definition |
| Ertl et al., *J Med Chem* **2000** | 2D TPSA definition: N and O contributions; polar-H included in fragment sums |
| rdkit/rdkit GitHub issue #1827 | `classifyAtoms()` silently fails for non-protein molecules; manual SASAClass fix documented |
| rdkit/rdkit GitHub issue #5197 | ArgumentError on list-of-float radii (fixed in RDKit ≥ 2022.03.4) |
| FreeSASA library docs | Protor/NACCESS/OONS classifiers all residue-name-based; OONS broadest but still fails without monomerInfo |
| bioRxiv 2026 "Delta PSA" preprint | Validates ΔPSA (MDAnalysis + GROMACS) correlates with experimental permeability across 500 macrocycles |

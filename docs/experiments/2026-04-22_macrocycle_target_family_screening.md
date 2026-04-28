# 2026-04-22 — Macrocycle Target-Family Screening Workflow

## Question

Given a macrocyclic or peptide-like hit compound with no known target, which protein families are plausible binders — and can a first-pass computational workflow prioritize them before committing to expensive selectivity assays?

This experiment designs a reproducible pipeline: CREST conformer ensemble → optional CENSO refinement → ensemble docking against a curated protein family panel → rank by cross-conformer docking consistency → advance top combinations to MD and SEEKR binding-rate estimation.

---

## Rationale

Macrocyclic peptides occupy an unusual chemical space: large enough to form multiple specific contacts with protein surfaces, flexible enough to undergo chameleonic conformational switching, and orally bioavailable (in favorable cases) without a defined target class. When a hit emerges from a permeability-first screen (e.g., high ΔPSA compound from CREMP analysis), there is no direct readout of what it binds.

The conformational ensemble is the key asset. Unlike small molecules where a single docked pose often suffices, macrocycles can present different binding faces depending on conformation. A docking workflow that ignores ensemble diversity will systematically miss binding modes accessible only from the minor (e.g., water) conformer. The approach here docks the full ensemble and ranks targets by how *consistently* the compound scores well — not by the single best pose.

SEEKR (Simulation-Enabled Estimation of Kinetic Rates) is added downstream because kon/koff ratios, not just ΔG, determine whether a hit is therapeutically relevant. A compound that binds fast and dissociates fast may not score well in an assay with a short incubation window.

---

## Literature Basis

### Ketzel et al. 2025 (*JACS* or similar) — Heterophyllin B Conformational Analysis

Cyclic octapeptide Heterophyllin B analyzed using CREST (GFN-FF pre-optimization + GFN2-xTB refinement, 6 kcal/mol energy window) → CENSO (PW6B95-D4/def2-TZVP/SMD solvation, single-point DFT correction) → final 2-conformer ensemble validated against RDC, NOE, and J-coupling NMR data.

**Relevance to this experiment:**
- Validates that CREST at the GFN2-xTB level correctly identifies the dominant solution conformers of a cyclic peptide without DFT geometry optimization
- Establishes 6 kcal/mol as a reasonable energy window for macrocyclic conformer sampling (captures NMR-visible minor conformer)
- CENSO DFT refinement improved energy ordering but did not change which conformers were dominant — suggests GFN2-xTB Boltzmann weights are reasonable for ensemble docking input
- Yang Hu (PI) is co-author — this workflow has implicit group endorsement

**Key methodological parameter:** In Ketzel 2025 the final ensemble was 2 conformers out of ~200+ sampled. This supports aggressive pruning before docking: cluster by RMSD, take one representative per cluster.

### CREST (Grimme group, *JCTC* 2019 + updates)

iMTD-GC algorithm for conformational sampling via GFN-FF + xTB. Already used in CREMP (CHCl3 ALPB) and our tier2 pipeline (CHCl3 + aqueous ALPB). For target-screening purposes:
- Membrane conformer (CHCl3) → likely binding-competent conformation for intracellular targets
- Aqueous conformer → relevant for extracellular/allosteric sites exposed to solvent

### Veber et al. 2002 (PSA + rotatable bonds)

Oral bioavailability cutoffs (PSA ≤140 Å², rotbonds ≤10) define which hits are worth pursuing orally. Any compound entering this screening workflow should first pass Veber filter — no point running SEEKR on a compound that won't reach its target.

---

## Proposed Workflow

### Step 1 — Conformer Ensemble Generation (CREST)

Input: SMILES of hit compound(s).  
Protocol: GFN-FF pre-optimization → GFN2-xTB iMTD-GC in CHCl3 ALPB (membrane-mimetic) and H2O ALPB (aqueous) as two separate runs. Energy window: **6 kcal/mol** (motivated by Ketzel 2025).  
Output: Two raw ensembles per compound (membrane + aqueous).

**Open question:** Should we use the 6 kcal/mol window from Ketzel (NMR-validated for solution conformers) or the CREMP default (which was tuned for computational efficiency over NMR coverage)?

### Step 2 — Ensemble Pruning and Clustering

Cluster each ensemble by backbone RMSD (suggested cutoff: 1.0 Å). Take Boltzmann-weighted representative per cluster (lowest-energy member as representative, cluster weight = sum of member Boltzmann weights).  
Target: ≤20 representative conformers per solvent condition per compound, weighted.

**Open question:** RMSD cutoff for macrocycles is ill-defined when ring pucker changes without backbone shift. Consider Torsional Fingerprint Deviation (TFD) as an alternative or supplement.

### Step 3 — Optional CENSO Refinement

For the top N conformers by Boltzmann weight (suggested: top 5 per solvent), apply CENSO single-point DFT correction (PW6B95-D4/def2-TZVP/SMD) to re-rank energies.  
This step is **optional** because:
- Requires TURBOMOLE or ORCA (not available on Jinich cluster without additional setup)
- Ketzel 2025 shows GFN2-xTB Boltzmann weights are qualitatively correct even without CENSO
- For target-family screening (not structure determination), GFN2-xTB ranking is likely sufficient

**Decision rule:** Skip CENSO unless a compound reaches MD stage (Step 6). Apply CENSO then for higher-accuracy input geometries.

### Step 4 — Protein Family Panel Curation

Curate a panel of ≥10 representative protein structures spanning pharmacologically relevant macrocycle target families:

| Family | Representative PDB | Rationale |
|---|---|---|
| Serine protease | e.g., thrombin, 1ppb | Classic macrocyclic drug targets |
| Metalloprotease (MMP) | e.g., MMP-3, 1g4k | Peptide-based inhibitor precedent |
| GPCR (extracellular loop) | e.g., CCR5, 4mbs | Cyclic peptide antagonists known |
| Kinase (allosteric) | e.g., CDK2, 1aq1 | Macrocycle kinase inhibitors in clinic |
| PPI surface (MDM2/p53) | e.g., MDM2, 1rv1 | Stapled peptide precedent; large flat surface |
| Nuclear receptor (LBD) | e.g., AR LBD, 2am9 | Peptidic modulators known |
| Epigenetic reader (BRD) | e.g., BRD4, 3mxf | Well-validated macrocycle target class |
| Ion channel (peptide-gated) | TBD | Cyclic peptide toxins bind here |

**Open question:** Which families are most relevant given the compound's origin (CREMP library = diverse cyclic peptides)? May need to restrict panel based on MW/PSA compatibility with each target class.

### Step 5 — Ensemble Docking

For each conformer representative (from Step 2 or 3), dock rigidly into each protein binding site. Rigid docking of pre-generated conformers is preferred over flexible docking because:
- Macrocycle ring flexibility during docking is poorly handled by most engines
- Pre-computed CREST ensemble already captures conformational diversity
- Docking speed is acceptable for ≤20 conformers × ≤10 proteins

Suggested engines (in order of preference for macrocycles): Glide SP/XP (Schrödinger), Gnina (open-source, uses CNN scoring), AutoDock-GPU (free, fast).

Score aggregation: **weighted mean docking score across ensemble**, using Boltzmann weights from Step 2. Also record the best-scoring conformer pose for visualization.

### Step 6 — Ranking and Advancement Criteria

Rank protein targets by:
1. **Weighted mean docking score** (primary)
2. **Pose consistency** — fraction of conformers that dock in a geometrically similar pose (cluster top poses by RMSD ≤2 Å; consistent binding = majority of conformers in same pose cluster)
3. **Chameleonic match** — does the best-docking conformer come from the membrane or aqueous ensemble? A membrane-conformer binder suggests the compound is pre-organized for binding upon membrane passage.

Advancement to MD: top 1–2 protein targets by combined score + consistency.

### Step 7 — MD and SEEKR

Standard MD (AMBER or GROMACS) for 100–500 ns on top compound–protein combinations. Assess binding stability; extract kon/koff via SEEKR milestoning protocol.

**SEEKR** (Simulation-Enabled Estimation of Kinetic Rates): milestoning-based approach that estimates kinetic rates from short MD trajectories anchored to milestones along a reaction coordinate. Provides kon, koff, and KD without requiring a full dissociation event. Available on SDSC resources.

---

## Key Hypothesis

A macrocyclic compound with high ΔPSA (chameleonic) will show better docking consistency across its conformer ensemble than a non-chameleonic compound of similar MW, because the chameleonic compound's membrane conformer is pre-organized for binding (buried polarity = exposed hydrophobic contacts). The aqueous conformer's exposed polar groups serve a different function (membrane crossing), not target binding.

---

## Open Questions

1. How many conformers should enter docking — all representatives, or only membrane-ensemble conformers for intracellular targets?
2. What energy window (kcal/mol) maximizes NMR-relevant conformer coverage without explosion of conformer count? (Ketzel: 6 kcal/mol → 2 relevant conformers; CREMP default may be tighter)
3. Which docking engine handles macrocycles best at our size range (8–13-mer cyclic peptides, MW 800–2000)?
4. How should docking scores be aggregated across a Boltzmann-weighted ensemble — weighted mean, best score, or a consistency-penalized score?
5. What criteria advance a compound to MD? Minimum weighted docking score threshold, or just top-N from panel?
6. Is CENSO refinement necessary before docking, or does GFN2-xTB geometry quality suffice for the scoring functions we use?

---

## Suggested Outputs

- Per-compound: conformer count per solvent, energy range, top-5 Boltzmann weights
- Docking summary table: compound × protein × weighted score × pose consistency × best-conformer identity (membrane vs. aqueous)
- Ranked protein family list with visualization of top-scoring pose per target
- SEEKR output: estimated kon, koff, KD for top 1–2 combinations
- Comparison of weighted-ensemble docking vs. single-lowest-energy-conformer docking (to validate whether ensemble weighting adds information)

---

## Notes

- This workflow is designed as a **first-pass triage**, not a hit-to-lead protocol. The goal is to rule out implausible target families and identify 1–2 worth pursuing experimentally.
- CENSO is deliberately deferred until the MD stage. DFT-level geometry refinement is appropriate for structure determination (Ketzel 2025 use case) and for input to high-accuracy scoring functions, but not for initial target-family screening.
- Protein panel selection should be guided by the compound's chemical matter. A compound with many aromatic side chains might prioritize BRD4/MDM2-type PPI surfaces; a compound with metal-chelating residues might prioritize metalloproteases.
- If the compound is from the CREMP PAMPA dataset, its permeability class is already known — a high-permeability ΔPSA-positive compound is the ideal input for this workflow because passive membrane transit is confirmed.

---

## References

- Ketzel et al. 2025 — Heterophyllin B conformational analysis (CREST + CENSO + NMR)
- Grimme et al. *JCTC* 2019 — CREST iMTD-GC algorithm
- Veber et al. *JMCA* 2002 — PSA + rotbonds oral bioavailability rules
- Navia & Chaturvedi *Drug Discovery Today* 1996 — hydrophilic collapse hypothesis
- SEEKR2: Votapka et al. *JPCB* 2022

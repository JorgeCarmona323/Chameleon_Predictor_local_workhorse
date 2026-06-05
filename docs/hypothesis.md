# Running Hypothesis: Chameleonic Behavior is Residue-Gated

**Status: Supported by preliminary evidence — under active investigation**
*Last updated: 2026-03-18*

---

## The Hypothesis

Chameleonic membrane permeability in cyclic peptides is not a continuous property of all cyclic scaffolds — it is **gated by molecular size**, specifically by the number of residues available to form stabilizing intramolecular hydrogen bonds in apolar environments. Below a critical residue threshold, cyclic peptides may achieve PAMPA permeability through other mechanisms (N-methylation, intrinsic lipophilicity) but do not undergo the conformational switching that defines true chameleonic behavior. Above that threshold, ΔPSA becomes a meaningful predictor of passive permeability precisely because the molecule has enough backbone flexibility to collapse and bury its polar surface area.

**The threshold is approximately 9 residues (~900 Da).**

---

## Evidence

### 1. MW gap in the 1,566-compound source-stratified dataset (this work, 2026-03-18)

Running UMAP Panel C on the clean Furukawa + Chugai subset and coloring by molecular weight reveals two well-separated populations:

| Population | Median MW | HDBSCAN perm rate |
|------------|-----------|-------------------|
| Permeable cluster (C0) | **1,180 Da** | 87.8% |
| Impermeable cluster (C1) | 820 Da | 59.3% |

The 1.44× MW gap between permeable and impermeable compounds means that in a clean, homogeneous dataset, **large cyclic peptides are disproportionately permeable**. The permeable cluster is not just "more permeable compounds" — it is physically larger molecules. This is consistent with size-gated chameleonicity: compounds large enough to form an intramolecular H-bond network in apolar solvent can collapse and cross membranes; smaller ones cannot.

Critically, this signal **completely disappears on the full 7k dataset** (median MW permeable = impermeable = 820 Da), confirming it is a real biological signal being masked by cross-source PAMPA label noise — not a statistical artifact.

See: [`docs/experiments/2026-03-18_panel_c_1566_track_d_mw.md`](experiments/2026-03-18_panel_c_1566_track_d_mw.md)

### 2. Yu et al. 2026 — independent computational evidence

Yu et al. (*bioRxiv*, DOI: [10.64898/2026.01.06.697862](https://doi.org/10.64898/2026.01.06.697862)) independently arrive at the same threshold through a different route:

- They compute **ΔPSA/SASA_total** (a dimensionless fractional switching ratio) rather than absolute ΔPSA
- They restrict analysis to compounds with **≥9 residues**, explicitly excluding smaller cyclic peptides on the grounds that chameleonic behavior does not reliably manifest below this size
- Their normalized descriptor is predictive where absolute ΔPSA fails — precisely because it removes the confounding size scaling that inflates raw ΔPSA for large peptides regardless of how chameleonic they actually are

The convergence is notable: Yu et al. arrive at the ≥9-residue cutoff from the descriptor side (normalization performance), while our MW gap result arrives at the same conclusion from the clustering/permeability side (empirical enrichment). Two independent lines of evidence pointing at the same threshold strengthens the hypothesis considerably.

### 3. CsA as the archetype

Cyclosporin A — the canonical chameleonic molecule — is 11 residues, MW 1,203 Da. It sits squarely above the threshold. Our Tier-1 ΔPSA for CsA = **84.9 Å²** vs. literature ~75 Å² (within 10%), and CsA falls in the high-MW permeable cluster in Track D. It is the proof of concept that large cyclic peptides can and do switch conformational state between aqueous and apolar environments.

---

## What the hypothesis does NOT claim

- Small cyclic peptides (hexapeptides, <9 residues) are impermeable. They can achieve PAMPA permeability — the C1 cluster is 59.3% permeable. They just do it through a different mechanism: N-methylation reducing HBD count, or intrinsic lipophilicity. Absolute ΔPSA is the wrong descriptor for them.
- MW alone predicts permeability. MW is a proxy for residue count, which is a proxy for the conformational freedom needed to collapse. The actual mechanistic variable is the ability to form intramolecular H-bonds in apolar solvent — ΔPSA/SASA_total is the better descriptor for that.
- The threshold is exactly 9 residues. That is Yu et al.'s cutoff on their dataset. Our data spans 6–15 residues (606–1,778 Da); the transition appears around 900 Da, which corresponds roughly to 8–9 residues at ~100–110 Da/residue for typical cyclic peptide monomers. The exact value needs to be tested empirically on clean data.

---

## Predictions (testable next experiments)

| Prediction | Test | Status |
|------------|------|--------|
| Normalized ΔPSA (per residue or per SASA) outperforms absolute ΔPSA | Add `delta_psa3d_per_residue`, `delta_psa3d_per_sasa` to conformer_engine.py; rerun AUC on 1,566-compound subset | Not yet run |
| AUC improves when analysis is restricted to ≥9-residue compounds | Filter `Monomer_Length >= 9`; compare AUC with full 1,566 | Not yet run |
| The MW gap is sharper when Chugai (unverified protocol) is excluded | Rerun on Furukawa only | Not yet run |
| The permeable cluster in Track D is enriched for known chameleonic scaffolds (CsA analogs, Choi macrolides) | Cross-reference C0 compound IDs against literature | Not yet run |

---

---

## Refinement (2026-06-05): Partition-Based Two-Regime Model

The original "threshold at 9 residues" framing assumes a sharp cutoff. A more realistic and testable version **partitions the dataset by residue count (6, 7, 8, 9, 10, 11, 12-mer)** and asks, per partition, *which descriptor family predicts permeability*. This lets the data reveal *where* and *how gradually* the mechanism transitions, rather than assuming a hard line at 9.

### The two-regime claim

Permeability operates via (at least) two size-dependent mechanisms:

- **Chameleonic regime (larger peptides):** solvent-driven switching between an open aqueous conformer and a collapsed membrane conformer. Governed by **ΔPSA, ΔΔG, cis-amide switching**. CsA (11-mer) is the archetype — experimental A1 (open, PSA 137) ↔ C1 (closed, PSA 96).
- **Pre-organized regime (smaller peptides):** permeability via *pre-organized* intramolecular H-bonding and intrinsic lipophilicity, **without** solvent-driven switching. Governed by **IMHB stability, exposed-HBD count, lipophilicity**.

### Prediction across partitions

For each residue-count bin, fit permeability against both descriptor families and record which dominates:

| Partition | Expected dominant signal |
|---|---|
| 6, 7, 8-mer | pre-organized: IMHB, lipophilicity, exposed HBD (chameleonic descriptors ~flat) |
| 9, 10-mer | **transition zone** — chameleonic descriptors begin to gain weight |
| 11, 12-mer | chameleonic: ΔPSA, ΔΔG, cis-switch dominate |

The **transition partition** (where chameleonic descriptors start carrying predictive weight) is the empirical "threshold" — likely a gradient, not a step.

### Mechanistic evidence from reference compounds (existence proof, not statistics)

- **CsA (11-mer):** experimental structures show the chameleonic two-state switch (A1 open ↔ C1 closed). *Note: our CREST V1 failed to reproduce this — see reliability caveats below.*
- **DOPC R/S, Brain1, DOPC2 (6-mers, permeable hits):** **anti-chameleonic** — ΔPSA *negative* (more polar surface in membrane than water); they do NOT switch, yet permeate. The R/S difference was in **intramolecular H-bonding and conformational pre-organization**, not ΔPSA. Direct support for a distinct small-peptide mechanism.
  → See `docs/experiments/2026-06-05_dopc_rs_3d_vs_2d_descriptors.md`

### Division of labor

- **Reference compounds (CREST):** validate the mechanism, justify the descriptor set, anchor partitions where available. *Cannot* establish the threshold statistically (too few; sampling/solvent caveats).
- **Full CycPeptMPDB (`feature_matrix.csv`, 1,566–7,000 compounds):** where the partition analysis and predictive claim are tested. Reference compounds currently anchor only the 6-mer and 11-mer bins; 7/8/9/10/12-mer partitions are populated by the database.

### Implementation note

The model carries **both** descriptor families for every compound and either (a) includes a residue-count × descriptor interaction term, or (b) is fit per-partition. A single chameleonic descriptor will fail on small peptides; a single lipophilicity descriptor will fail on large ones.

---

## Reliability of Current Evidence (2026-06-05)

| Evidence | Reliability | Why |
|---|---|---|
| CsA experimental A1/C1 (PSA, cis, Rg) | **High** | Crystallography (X-ray/neutron, CCDC) |
| CsA CREST V1 ensemble | **Low** | Missed cis MeVal11–MeBmt1 (no `-notopo`), over-collapsed (ALPB implicit), single-start. Wrong ensemble; superseded by CsA_v2. |
| Exp-vs-CREST V1 comparison (as a diagnostic) | **High** | Correctly identified the V1 failures |
| DOPC R/S relative comparison | **Moderate** | Correct pipeline (`-notopo`, CREST 2.12); R-vs-S differences from uncapped water ensembles. But single-start, implicit solvent, no experimental validation. |
| DOPC R/S absolute ΔPSA | **Low** | mem capped at 50, over-collapse, sub-threshold 6-mer; use relative/normalized only |

---

## References

- Yu et al. (2026). *bioRxiv*. DOI: 10.64898/2026.01.06.697862
- Witek et al. (2016). *J. Chem. Theory Comput.* — CsA conformational ensemble in polar/apolar solvents
- Limbach et al. (2025). *J. Med. Chem.* — biased equilibrium / Goldilocks barriers
- Bockus et al. (2015). *J. Med. Chem.* — decoding chameleonic properties of macrocycles
- Rezai et al. (2006). *J. Am. Chem. Soc.* — conformational flexibility and passive permeability

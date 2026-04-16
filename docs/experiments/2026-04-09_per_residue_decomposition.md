# 2026-04-09 — Per-Residue ΔPSA Decomposition: Framework and Linker Experiment

## The Core Question

Is per-residue ΔPSA decomposition necessary for the thesis, or does whole-molecule ΔPSA
capture sufficient signal? The linker swap experiment tests this empirically before committing
to a complex implementation.

---

## Linker Swap Experiment (Test First)

**Hypothesis:** If swapping the cyclization linker while holding all other building blocks
constant produces a large change in permeability, then:
- The linker is an independent conformational determinant
- Whole-molecule ΔPSA conflates linker and residue contributions
- Per-residue decomposition is essential for design rules

**If linker swap → small Δpermeability:**
- Whole-molecule ΔPSA is sufficient for prediction
- Per-residue decomposition is scientifically interesting but not practically blocking

**How to run it now:** The Macrocycle DEL library UMAP (8,507 compounds, MAPchiral +
Mordred already computed) likely contains scaffold-matched pairs — same residues at X1–X4,
different cyclization chemistry. Check if such pairs exist in different UMAP regions or show
divergent predicted ΔPSA values.

**This is a thesis sub-question:** "Does linker chemistry independently modulate chameleonic
behavior, or does it merely constrain the scaffold's intrinsic conformational landscape?"

---

## Per-Residue ΔPSA: Implementation Path

Computationally straightforward — rdFreeSASA already gives per-atom SASA. Just need
atom → residue mapping.

### Layer 1 — Per-Residue ΔPSA (what changes)

For each molecule and its aqueous/membrane conformer pair:
1. Run `rdFreeSASA.CalcSASA()` per atom on both conformers
2. Map atoms to residues via RDKit atom properties or SMARTS-based residue detection
3. Sum polar atom SASA within each residue → residue_i_PSA(aq), residue_i_PSA(mem)
4. ΔPSA_i = residue_i_PSA(aq) - residue_i_PSA(mem)
5. Linker treated as its own "position" with ΔPSA_linker

Note: ΔPSA is additive per conformer (sum of per-atom SASA = total SASA), so
sum(ΔPSA_i) = ΔPSA_total. Decomposition is internally consistent.

**Interpretation of ΔPSA_i:**
- High |ΔPSA_i|: residue i is conformationally active — its polar exposure changes strongly
- ΔPSA_i ≈ 0: passive — exposure unchanged between conformers
- ΔPSA_i < 0: anti-chameleonic — more exposed in membrane conformer than aqueous

**Limitation:** A residue can have high |ΔPSA_i| because:
  a) It directly forms/breaks IMHBs (mechanistically active)
  b) Neighboring residues pack around it (passive passenger)
  These cannot be distinguished from ΔPSA alone.

### Layer 2 — IMHB Inventory (why it changes)

For each conformer, catalog intramolecular hydrogen bonds using RDKit's
`rdMolDescriptors.CalcNumHBD/HBA` or a custom distance/angle criterion.

- Compare aqueous vs membrane conformer IMHB patterns
- Residues that **gain** IMHBs in the membrane conformer are active drivers
- Residues that are buried but don't gain IMHBs are passive passengers
- Linker atoms participating in IMHBs in the membrane conformer → linker is a chameleonic element

This is the distinction between prediction and design — Layer 1 tells you what, Layer 2 tells you why.

### Layer 3 — Torsion Tracking (how it transitions)

Backbone torsion angles (φ, ψ) for each residue in each conformer.
- Track which torsions differ most between aq and mem conformer
- Residues with large Δtorsion are mechanical pivot points
- Concerted torsion changes across multiple positions indicate cooperative switching

ICoN-v1 (Hung et al., bioRxiv 2026, in literature/) does this with deep learning on MD
trajectories. For CREST ensembles, simpler: extract φ/ψ per residue per conformer and
compute circular standard deviation across the ensemble.

---

## The AGDIFF Constraint and Why Aqueous CREST Is Non-Negotiable

AGDIFF (Wu & Zou, JCIM 2026) was trained exclusively on CREST CHCl3 conformers.
It cannot generate aqueous conformers — polar group exposure in water is absent from its
learned distribution.

**Consequence for the ΔPSA framework:**

| Method | Aqueous conformer | Membrane conformer | ΔPSA reliability |
|---|---|---|---|
| CREST CHCl3 only | Absent | Good | Underestimates aq PSA |
| CREST aqueous only | Good | Absent | Underestimates mem PSA |
| CREST CHCl3 + aqueous | Good | Good | Full signal |
| AGDIFF | None | Good | Severely compressed — not valid |
| ETKDG | Vacuum extremes | Vacuum extremes | Weakest |

**AGDIFF role:** Potential fast membrane-conformer engine for 800K DEL compounds, but
cannot replace dual-CREST for ΔPSA. Could contribute F9 (membrane-only signal) in the
feature benchmark as an upper-bound test for single-solvent utility.

**Aqueous CREST is not optional.** It is the only path to valid aqueous conformers and
therefore the only way to compute physically meaningful ΔPSA at scale.

**Thesis argument:** Existing models (including AGDIFF) fail at full chameleonic ΔPSA
because they are single-solvent. Dual-CREST is a genuine methodological contribution —
not an incremental improvement, but the minimal requirement for the physics to be correct.

---

## Decision Tree

```
Run linker swap experiment
        │
        ├── Linker strongly modulates permeability?
        │         YES → per-residue ΔPSA is essential
        │               implement Layer 1 + 2 (ΔPSA + IMHB)
        │               linker = independent design variable
        │
        └── Linker weakly modulates permeability?
                  NO → whole-molecule ΔPSA sufficient for prediction
                       per-residue adds interpretability but not prediction
                       implement Layer 1 only, use for design rules post-hoc
```

---

## References

- Wu & Zou, *JCIM* 2026. DOI: 10.1021/acs.jcim.5c03236 (AGDIFF)
- Hung, Venkatesan & Chang, bioRxiv 2026. DOI: 10.64898/2026.03.12.711417 (ICoN-v1)
- Yu et al., bioRxiv 2026. DOI: 10.64898/2026.01.06.697862 (Delta PSA normalization)
- Grambow et al., *Scientific Data* 2024. DOI: 10.1038/s41597-024-03698-y (CREMP)

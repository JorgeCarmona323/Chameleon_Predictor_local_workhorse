# Permeability Design-Loop Architecture — Concept

**Date:** 2026-07-28
**Status:** Concept / design proposal (not yet implemented)
**Author:** Jorge (with Claude)

---

## 1. Motivation

Three near-simultaneous works define the current cyclic-peptide permeability landscape, and each
leaves the same gap:

- **EnsembleCycPerm** (Wen 2026, *JCIM*): solvent-dependent conformational ensembles → permeability
  ML. Independently found **ΔPSA3D(CHCl₃−H₂O)** is the dominant descriptor. Ensembles from
  RDKit ETKDGv3 + GFN2/ALPB (cheap search, 5 reps/solvent). Interpretability = gradients on a
  black-box neural net. **Teacher is a neural net.**
- **MycoPermeNet** (mycobacterial outer-membrane permeation): a strong message-passing NN
  **teacher** whose predictions label a large space, then an interpretable **Random Forest student**
  on 2D features + SHAP elucidates which chemical features drive permeation. **Teacher is a neural net.**
- **PEGASUS**: large-scale AI permeability model on massively parallel assays. Data + model; no
  physics, no glass-box distillation.

**None of them has a physics teacher, and none designs/makes molecules.** Our pipeline (CREST
iMTD-GC + CPCM-X ΔG_transfer, NMR-validated) is a *mechanistic* label generator, and the Hu lab
closes the loop by synthesizing and assaying. That is the opening this architecture exploits.

## 2. The architecture in one sentence

A **multi-fidelity, experiment-calibrated 3D→permeability decision-maker**, made interpretable by a
**2D student** that doubles as a **fast forward-surrogate**, wrapped in a **two-tier 2D-mutation
design loop** that proposes the next molecule to synthesize.

```
        (thousands of peptides: CycPeptMPDB + PEGASUS + our hits)
                               │
   TEACHER   ┌─────────────────┴──────────────────┐
  (physics)  │  CREST + CPCM-X  →  ΔG_transfer +   │   dense LOW-fidelity labels
             │  3D descriptors (PSA3D, Rg, IMHB…)  │
             └─────────────────┬──────────────────┘
                               │   + sparse HIGH-fidelity EXPERIMENTAL PAMPA labels
   DECISION  ┌─────────────────┴──────────────────┐
   MAKER  f  │  3D descriptors → permeable?        │   fused: physics shapes manifold,
             │  (multi-fidelity: see §3)           │   experiment calibrates truth
             └─────────────────┬──────────────────┘
                               │
   STUDENT   ┌─────────────────┴──────────────────┐
  (2D, RF)   │  2D features (Mordred + monomer id  │   (a) SHAP interpretability
             │  + STEREO + linker) → 3D signal /   │   (b) FAST 2D→3D surrogate for the loop
             │  permeability                       │
             └─────────────────┬──────────────────┘
                               │
   DESIGN    ┌─────────────────┴──────────────────┐
   LOOP      │  mutate 2D → predict 3D signal      │   Tier 1 (fast, surrogate): 100s of mutants
             │  → score f → keep improvements       │   Tier 2 (slow, real physics): top ~5 only
             └─────────────────┬──────────────────┘
                               │
                        → propose for synthesis + assay → new experimental labels (loop closes)
```

## 3. Component 1 — the multi-fidelity decision-maker (the crux)

Experimental PAMPA (high-fidelity, sparse) and physics ΔG_transfer (low-fidelity, dense) are
**different quantities on different scales** — they must NOT be pooled as one label. Three valid
fusion strategies, in increasing sophistication:

1. **Physics-as-feature (start here):** ΔG_transfer + 3D descriptors are *inputs*; train `f` only on
   experimental labels. Physics informs the features; experiment is the sole truth. Simplest, hardest
   to get wrong.
2. **Pretrain → finetune:** pretrain `f` on the dense physics labels (learns the *shape* of the
   descriptor→permeability manifold — the "in-between points fill the gaps" intuition), then fine-tune
   on the sparse experimental points (calibrate to truth).
3. **Δ-learning (most data-efficient):** model `perm_exp = ΔG_transfer + δ(features)`; learn only the
   small correction δ on experimental data. Physics does the heavy lifting; the model learns only
   where physics is wrong. Aligns with `multifidelity_vacuum_implicit`.

## 4. Component 2 — the interpretable 2D student (dual role)

Following MycoPermeNet's distill-to-glass-box pattern, but with **our physics as the teacher**: label
a large peptide space with the physics/decision-maker, fit an interpretable **RF/GBM student** on
**readable 2D features**, and SHAP it → *which structural features produce the rewarded 3D-descriptor
states*. Training on the teacher's **dense** label surface (not the sparse experimental set) is what
makes SHAP importances **stable despite small experimental data** ("label amplification").

**Key insight — the student is also the fast forward-model.** The same 2D→3D-descriptor map that
gives interpretability is the millisecond surrogate that makes the design loop (§5) computationally
possible. Interpretability and speed come from one model.

**Features MUST include explicit stereo (R/S) and linker (xylene/diazirine) descriptors** — otherwise
the 2D layer is blind to our actual design axes.

## 5. Component 3 — the two-tier design loop

CREST + CPCM-X is **hours per molecule** — the vacuum/RDKit descriptors are cheap, the propensity
(CPCM-X ΔG) is not. So the loop cannot run full physics on every candidate. It runs in two tiers
(cascade, per `multifidelity_vacuum_implicit`):

- **Tier 1 (fast, inside the loop):** mutate 2D structure → predict 3D descriptors from 2D via the
  student surrogate → score with `f` → keep improvements. Hundreds of mutants, milliseconds each.
- **Tier 2 (slow, winners only):** run real CREST + CPCM-X on the top ~5 to *confirm* the 3D signal
  moved as the surrogate predicted, then propose for synthesis.

## 6. Hard problems (eyes open)

1. **Permeability is not a deterministic function of a few 3D descriptors.** EnsembleCycPerm calls it
   "context-dependent"; the goldilocks paper (Limbach 2025) frames it as a **kinetic/barrier**
   phenomenon (membrane-flip rate), not a single descriptor. The decision-maker has a **ceiling** —
   treat it as "enrich for likely-permeable," a strong prior/filter, not an oracle.
2. **Stereochemistry — our main design axis — is the hardest thing for the 2D surrogate.** 2D
   fingerprints barely separate R/S, yet our epimers differ in 3D descriptors AND permeability. So
   the 2D→3D map for *stereo mutations* is where the surrogate is weakest. Mitigation: give the
   surrogate explicit stereo descriptors, and treat **stereo flips as a physics-only move** (always
   Tier-2 validated). Linker/side-chain swaps are more 2D-visible and safe for Tier 1.
3. **Optimizing against our own model risks chasing model artifacts.** The loop must close
   out-of-distribution — Tier-2 physics, then synthesis + assay — not just self-consistency.
4. **Two-label-scale fusion** (see §3) — the most common way to get this wrong is naive pooling.

## 7. Staged de-risking plan

1. **Fusion sanity check** — physics-as-feature RF/GBM on experimental labels; confirm ΔG_transfer +
   descriptors predict PAMPA on CycPeptMPDB (we already know ΔPSA3D does; EnsembleCycPerm confirms).
2. **Build the fast 2D→3D-descriptor surrogate**; measure error, **especially on epimer pairs** — this
   tells us whether stereo moves can ever be Tier 1.
3. **SHAP the student** → interpretable "which 2D features → good 3D signal" readout.
4. **Closed loop on linker/side-chain moves only** (stereo excluded from Tier 1); validate top 3–5
   designs with full physics → propose for synthesis.

## 8. Why this is ours to publish

Only this design has (a) a **physics teacher** — a mechanistic ΔG_transfer, not a learned
approximation — and (b) a **wet lab that closes the loop** by making the molecule. The output is
simultaneously an interpretable ML result, grounded in thermodynamics not a black box, and an
**actionable design rule** (next probe: R or S? xylene or diazirine?) that neither EnsembleCycPerm,
MycoPermeNet, nor PEGASUS produces.

## References

- Memories: `project_model_strengthening_three_repos`, `reference_ensemblecycperm_2026`,
  `project_multifidelity_vacuum_implicit`, `project_model_vision_layered`,
  `project_dgtransfer_literature_comparison`, `project_dpsa_mechanism_scope`.
- Repos: `reference repos/EnsembleCycPerm`, `reference repos/Mycomembrane-permeability-project`,
  `reference repos/PEGASUS1`.

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
**2D student** that doubles as a **fast forward-surrogate**, wrapped in a **cost-tiered, two-branch
2D-mutation design loop** — fronted by a **structural-change dial** — that rescues a user's
non-permeable lead by proposing the minimal edits that flip its 3D chameleon signal.

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
  (2D, RF)   │  2D DESCRIPTORS (Mordred + monomer  │   (a) SHAP interpretability
             │  id + stereo tags) → 3D signal      │   (b) FAST 2D→3D surrogate for the loop
             └─────────────────┬──────────────────┘
                               │
   DESIGN    ┌─────────────────┴──────────────────┐
   LOOP      │  USER INPUT: a lead they want made  │   Tier 0  structural-distance dial (2D, free)
             │  permeable                          │   Tier 1  surrogate scores survivors (fast)
             │  → 2D-notation mutations (2 branches)│   Tier 2  real CREST+CPCM-X, top ~5 (slow)
             │  → rank by GAIN-PER-EDIT            │
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
**readable 2D descriptors**, and SHAP it → *which structural features produce the rewarded
3D-descriptor states*. Training on the teacher's **dense** label surface (not the sparse experimental
set) is what makes SHAP importances **stable despite small experimental data** ("label amplification").

**Key insight — the student is also the fast forward-model.** The same 2D→3D-descriptor map that
gives interpretability is the millisecond surrogate that makes the design loop (§5) tractable.
Interpretability and speed come from one model.

**Descriptors vs notation — do not conflate them:**
- **2D descriptors** (Mordred/fingerprints) are largely **stereo-blind** → the fast *scoring
  surrogate* cannot feel an R→S flip. Mitigate with explicit stereo/CIP-tagged features, but the
  effect is inherently 3D (see §5 branch B).
- **2D notation** (SMILES/HELM) is fully **chirality-aware** → the *move generator* enumerates stereo
  mutations with no loss. The stereo problem is a *scoring* problem, never a *generation* problem.

## 5. Component 3 — the design loop (two branches, three cost tiers)

CREST + CPCM-X is **hours per molecule** — the vacuum/RDKit descriptors are cheap, the propensity
(CPCM-X ΔG) is not. So we never spend a 3D permeability judgment on every mutant. The loop is a
**cost cascade**:

- **Tier 0 — structural-distance dial (2D only, free):** enumerate mutations in 2D notation, compute
  structural distance (Tanimoto / monomer edit-distance) to the input, and keep only those inside the
  user-chosen band (see §6). No model, no physics — this collapses the combinatorics first.
- **Tier 1 — surrogate (fast):** the 2D→3D student predicts the 3D signal / permeability for survivors.
  Hundreds of candidates, milliseconds each.
- **Tier 2 — physics (slow, winners only):** real CREST + CPCM-X on the top ~5 to *confirm* the 3D
  signal moved as the surrogate predicted, then propose for synthesis.

**Two-branch move-set.** At each position the decision-maker flags as high-leverage (monomer-level
attribution / SHAP), spawn two branches of variation with *different cost profiles*:

- **Branch A — constitutional** (swap monomer / side chain / linker): moves 2D descriptors, so the
  **surrogate can pre-screen** many candidates in Tier 1.
- **Branch B — stereochemical** (invert the center at that site): the descriptor surrogate is blind,
  but there are **very few** stereo variants per site (usually 1), so route them **straight to Tier 2
  physics**. Bounded → affordable. *Notation generates it; physics scores it.*

Then score in 3D and **keep whichever branch wins at that site**. Note: attribution tells you *where*
a change matters, not *which direction* helps — both branches must be *scored* to find the gain.

## 6. Product / UX layer — the structural-change dial

**Service thesis (the input contract):** a user submits a molecule they believe is **not permeable**
and want **made permeable**. Nobody submits an already-permeable lead — so the objective is always
"**rescue this lead**", a sharper, more sellable value prop than generic permeability prediction.

**The dial.** Rather than pay a 3D permeability judgment for every possible edit, the user first sets
**how far they are willing to change the structure** — a slider (single value or range, continuous
similarity threshold *or* discrete number of monomer edits) that drives the **Tier-0 gate**. This is
the medicinal-chemistry "**how far can I stray from my lead**" control: small edits preserve target
activity/IP but may not move the 3D signal enough; large edits can swing permeability more but risk
losing binding and wandering where the model is less reliable. The dial navigates that
**conservatism-vs-impact tradeoff** explicitly and cheaply.

**Ranking = gain-per-edit.** Because a *tiny* 2D edit (a stereo flip, one N-methylation) can produce a
*large* 3D/permeability change while a big edit may barely move it, the objective inside the chosen
budget is **permeability-gain per unit structural distance** — surfacing high-leverage, minimally
invasive edits. The dial sets the budget; the 3D layer picks the winners within it.

## 7. Hard problems (eyes open)

1. **Permeability is not a deterministic function of a few 3D descriptors.** EnsembleCycPerm calls it
   "context-dependent"; the goldilocks paper (Limbach 2025) frames it as a **kinetic/barrier**
   phenomenon (membrane-flip rate), not a single descriptor. The decision-maker has a **ceiling** —
   treat it as "enrich for likely-permeable", a strong prior/filter, not an oracle.
2. **Structural distance is a search-control, not the objective.** "Most different" ≠ "most permeable";
   it only bounds the neighborhood. The permeability judgment always comes from the 3D layer.
3. **Optimizing against our own model risks chasing model artifacts.** The loop must close
   out-of-distribution — Tier-2 physics, then synthesis + assay — not just self-consistency.
4. **Two-label-scale fusion** (§3) — the most common way to get this wrong is naive pooling.

## 8. Staged de-risking plan

1. **Fusion sanity check** — physics-as-feature RF/GBM on experimental labels; confirm ΔG_transfer +
   descriptors predict PAMPA on CycPeptMPDB (we already know ΔPSA3D does; EnsembleCycPerm confirms).
2. **Build the fast 2D→3D-descriptor surrogate**; measure error, **especially on epimer pairs** — this
   tells us whether stereo moves can ever be Tier 1 or must stay Tier 2.
3. **SHAP the student** → interpretable "which 2D features → good 3D signal" readout.
4. **Tier-0 gate + dial prototype** — structural-distance enumeration + slider on a toy scaffold.
5. **Closed loop** — two-branch moves, rank by gain-per-edit; validate top 3–5 with full physics →
   propose for synthesis.

## 9. Why this is ours to publish

Only this design has (a) a **physics teacher** — a mechanistic ΔG_transfer, not a learned
approximation — and (b) a **wet lab that closes the loop** by making the molecule. The output is
simultaneously an interpretable ML result, grounded in thermodynamics not a black box, and an
**actionable design service**: give it a non-permeable lead and a structural-change budget, get back
the minimal edits predicted to make it permeable — which neither EnsembleCycPerm, MycoPermeNet, nor
PEGASUS produces.

## References

- Memories: `project_design_loop_architecture`, `project_model_strengthening_three_repos`,
  `reference_ensemblecycperm_2026`, `project_multifidelity_vacuum_implicit`,
  `project_model_vision_layered`, `project_dgtransfer_literature_comparison`,
  `project_dpsa_mechanism_scope`.
- Repos: `reference repos/EnsembleCycPerm`, `reference repos/Mycomembrane-permeability-project`,
  `reference repos/PEGASUS1`.

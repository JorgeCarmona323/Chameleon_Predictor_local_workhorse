# Diazirine N=N constraint — design overview (for review)

**2026-06-14 · reviewer-facing.** Implements a CREST/xTB constraint to stop GFN2-xTB from
corrupting the diazirine N=N during conformer sampling. Please sanity-check the design,
especially the "open questions" at the end.

---

## 1. Problem

GFN2-xTB/CREST relaxes the diazirine **N=N to ~1.43 Å** (an N–N *single*-bond length, with the
3-ring going equilateral, ∠N–C–N ≈ 60° instead of ~49°) in **5 of 8** ensembles. It's *systematic*,
not random thermal breakage — `nn_std` ≈ 0.001 Å across hundreds of conformers, i.e. the bad
geometry was baked into the CREST seed and frozen through the whole search. Knock-on effect: the
distorted ring perturbs the backbone sampling and produced an **artifactual descriptor result**
(a spurious backbone-IMHB collapse in `3-12-10-12_dz_R`).

## 2. Root cause (isolated by elimination)

| stage | diazirine N=N | conclusion |
|---|---|---|
| SMILES (2D graph) | bond order = 2.0 (all 8) | topology correct — not an encoding bug |
| RDKit ETKDG / MMFF | 1.24–1.27 Å | starting geometry correct — not embedding |
| **GFN2 pre-opt / CREST** | **1.43 Å (5/8), 1.23 Å (3/8)** | **GFN2 corrupts it**, path-dependent |

Literature truth: parent 3H-diazirine N=N = **1.228 Å** (microwave, rs) / **1.230 Å**
(CCSD(T)/cc-pVQZ); the CF₃-phenyl diazirine (TPD) sits the same (X-ray + B3LYP). So 1.43 Å is
unambiguously a GFN2 artifact — a spurious elongated-N=N basin reachable on some optimization
paths (all four `3-12-8-12 dz` runs fell in; only `10-12 dz R/water` did for that compound).
**Why GFN2 struggles here:** a 3-membered ring = extreme angle strain, two adjacent N heteroatoms,
azo-like electronic structure — the regime where semi-empirical tight-binding parameters are weakest.

## 3. Fix (implemented — "Fix 1", distance constraint)

**Auto-detection** (`crest_conformers.py::diazirine_nn_atoms`): SMARTS `[#6]1[#7]=[#7]1` on the
H-added RDKit mol → the two ring N atoms → `+1` to convert to the **1-based** indexing xtb/CREST use.
Returns `None` (no-op) if no diazirine, so all other compounds are byte-for-byte unchanged.
Indices are consistent because the same RDKit atom order is written to every `.xyz` the tools read.

**Constraint file** (`write_constraint_file`), built once per compound, reused for both solvents:
```
$constrain
  force constant=0.25
  distance: <N1>, <N2>, 1.23
$end
```
(verified on the real 103-atom macrocycle: detects N=N at indices 19,20.)
**force constant = 0.25** Eₕ/Bohr² — the CREST-documented value for *distance* constraints
([CREST docs, example 4](https://crest-lab.github.io/crest-docs/page/examples/example_4.html#constrained-sampling);
the docs caution against high values). No `$metadyn` block or reference file is needed for a
distance constraint (those are only for substructure fixing = Fix 2). At r₀=1.23 Å, 0.25 already
puts a ~20 kcal/mol barrier in front of the 1.43 Å drift; escalate to 0.5→1.0 only on WATCH/FAIL.

**Applied at every GFN2/ALPB stage** of the pipeline:
- xTB pre-optimization → `xtb ... --input <file>`
- CREST iMTD-GC search → `crest ... --cinp <file>`
- CREST `--cregen` refinement → `crest --cregen ... --cinp <file>`

So the rigid diazirine is pinned to its physical N=N from the very first optimization, while the
macrocycle backbone and side chains sample freely.

## 4. Why the constraint does not bias the result

1. **Orthogonal to sampling.** N=N is a stiff, local degree of freedom of a rigid pendant.
   Pinning the N=N distance removes only an unphysical bond-stretching degree of freedom; it does
   not restrict the macrocycle torsions, side-chain rotamers, or pendant orientation that define the
   conformational ensemble.
2. **Near-constant in the energetics (not "exact cancellation").** The same 1.23 Å target is enforced
   in every conformer, so any restraint contribution is expected to be near-constant and should not
   drive Boltzmann reweighting among macrocycle conformers. (Stated carefully: the local diazirine
   environment could couple weakly to nearby sterics/electrostatics, so the restraint removes the
   unphysical GFN2 basin and contributes *at most* a near-constant local energy term — it does not
   "cancel exactly.")
3. **Physical target.** 1.23 Å is the literature-correct value, so the constraint holds the molecule
   at its true geometry — no artificial strain. It only walls GFN2 out of its spurious basin.

## 5. Verification

`scripts/verify_diazirine_integrity.py` (no sims; reads ensemble SDFs): per-conformer N=N, both C–N
distances, ∠N–C–N, and per-ensemble consistency (std). Pre-fix run: 5 FAIL / 3 PASS.
**Expected post-fix: all PASS, N=N ≈ 1.23 Å.**

Post-fix acceptance window (reviewer-set, 2026-06-14):
- **PASS:** N=N ≈ 1.22–1.25 Å
- **WATCH:** 1.25–1.30 Å (constraint holding but drifting — note, don't necessarily re-run)
- **FAIL:** any conformer > 1.35 Å, or C–N > 1.60 Å (ring opened)
- plus ∠N–C–N and C–N monitored as before.

The integrity check now also monitors the **terminal alkyne** (C≡C ≈ 1.20 Å, C–C≡C linearity ≈ 180°)
as a *monitor-only* motif — not constrained unless it flags a reproducible artifact.

## 6. Escalation path (if Fix 1 is insufficient)

If a post-fix integrity re-check shows the N=N held but the **C–N bonds distorted**, escalate to
**Fix 2 (rigid ring)** — freeze the 3 ring atoms' positions against a clean reference:
```
$constrain
  atoms: <C>,<N1>,<N2>
  force constant=1.0
  reference=coord.ref
$end
```
This preserves the exact ring geometry while still letting the C–C bonds to the macrocycle flex.
~2 lines to add; not enabled by default (minimal-constraint principle).

## 7. Open questions for the reviewer

1. **Force constant** — set to **0.25** Eₕ/Bohr² (the CREST-documented value for distance
   constraints; see §9). Must hold against CREST's metadynamics "kicks" without over-biasing;
   escalate to 0.5→1.0 only if the integrity check shows drift.
2. **`$metadyn atoms:` exclusion?** Should we also exclude the two N's from the RMSD bias (list all
   *other* atoms in a `$metadyn` block)? The MTD kicks are the breaking mechanism, so excluding the
   diazirine from the bias is belt-and-suspenders; the fc=1.0 distance constraint may already suffice.
   Currently NOT included (matches the minimal `$constrain`-only recipe).
3. **Other unusual motifs — constrain or just monitor?** The molecule also has a **terminal alkyne**
   (pent-4-ynoic acid handle, C≡C) and a **xylylene linker** (aromatic + benzylic thioether). Current
   stance: *monitor, don't constrain* — the xylene is validated (the non-diazirine analogs share it and
   sampled cleanly); the alkyne is low-risk (linear sp, unstrained) but unverified. Plan is to add C≡C
   length (~1.20 Å) + linearity (~180°) to the integrity check and only constrain if it flags. Agree?
4. **`--cregen` constraint** — we pass `--cinp` to the cregen step too (it may re-optimize). Harmless?

## 8. Files

- `scripts/crest_conformers.py` — `diazirine_nn_atoms`, `write_constraint_file`, threaded through
  `xtb_preopt_mol` / `run_crest` / `run_crest_cregen` / `process_compound`.
- `scripts/verify_diazirine_integrity.py` — the integrity gate.
- `memory/project_diazirine_review_checklist` — running record of the finding + fix.

## 9. Review outcome (2026-06-14, second-agent review — APPROVED)

Proceed with **Fix 1 as implemented.** Specific dispositions of the open questions:
1. **force constant** — reviewer suggested 1.0, but the CREST docs (example 4) recommend **0.25** for
   *distance* constraints and caution against high values, so we set **0.25** (corrected post-review).
   This keeps the reviewer's "minimal constraint, escalate if needed" logic with the docs-correct
   start; stiffen to 0.5→1.0 only if N=N drifts. Acceptance window adopted (PASS 1.22–1.25 /
   WATCH 1.25–1.30 / FAIL >1.35 Å), now in the integrity checker.
2. **`$metadyn atoms:` exclusion** — do NOT add yet (keeps the sampling definition clean / one fewer methodological difference). Escalate only if N=N stays unstable despite the distance constraint, or C–N distorts while N=N is pinned.
3. **Other motifs** — monitor, don't constrain. Added alkyne C≡C distance + linearity to the checker; aromatic ring sanity only if failures appear.
4. **`--cregen` + `--cinp`** — yes, keep (protects the diazirine through the final refinement; harmless).
5. **Escalation** — rigid 3-atom ring (Fix 2) only if C–N / ∠N–C–N still fail post-fix. Not enabled.

The "cancels exactly" wording in §4 was softened to "near-constant local term" per reviewer (more
reviewer-proof; the local diazirine environment can couple weakly to nearby sterics/electrostatics).

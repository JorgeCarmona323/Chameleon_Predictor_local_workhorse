# DOPC isomer analysis — per-molecule logic, building blocks, and the cis-amide signature

**2026-06-17 · analysis note for the 3-12-8-12 and 3-12-10-12 R/S pairs (normal xylene linker)**

This note records the *reasoning workflow* behind the two isomer reports, the verified building-block
assignment, and a cis-amide consistency check. It is the analytical backbone the reports rest on.

---

## 0. Building blocks — VERIFIED from the structures (not assumed)

The two scaffolds differ at one backbone residue. Verified directly from the SMILES the descriptors
were computed on (`scripts` substructure check):

| compound | residue | structural test | character |
|---|---|---|---|
| **3-12-8-12** | **L-Azetidine-2-carboxylic acid** | 4-membered N-ring present; no N-methyl | rigid ring, proline-homolog, **strong cis-amide inducer**, no backbone N–H |
| **3-12-10-12** | **Sarcosine** (N-methylglycine) | N-methyl amide present; no 4-ring | flexible hinge, weak cis tendency, no backbone N–H |

(Both are N-substituted → **both lack a backbone N–H donor**; the difference is conformational, not
donor count. An earlier draft had these swapped; corrected here against the data.)

---

## 1. The analysis workflow

To separate signal from fluff and avoid cross-talk, each molecule is analyzed **independently first**:
1. State the observation from each descriptor family (H-bonds, surface, shape, flexibility, cross-solvent).
2. Check **internal consistency** — do the families agree on one physical picture? (Agreement across
   independent descriptor families is the validity check; disagreement flags an unreliable descriptor.)
3. Draw a per-molecule conclusion from the consistent majority.
4. *Only then* compare/contrast epimers, and across scaffolds.

This makes differences **emerge** from the data rather than being imposed, and exposes which descriptors
are real signal vs. discretization/sampling artifacts.

---

## 2. Per-molecule conclusions (water phase)

**3-12-8-12 R (azetidine):** backbone IMHB low (1.6); SA_HD high (79, donors exposed); amphi high (3.4);
one fold (basin 0.86, RMSF 0.62); closes entering membrane. *Consistent:* polar groups point outward.
→ **single open fold that displays donors/polar face; solvent-responsive (chameleon).**

**3-12-8-12 S (azetidine):** backbone IMHB high (2.9); SA_HD low (45, buried); amphi low (1.9); **floppy**
(RMSF 0.96, many folds, basin 0.47). *Consistent:* donors sequestered internally — but reached via many
geometries. → **closed but conformationally diverse (floppy-closed).**

**3-12-10-12 R (sarcosine):** backbone IMHB low–mid (2.0); SA_HD high (88); amphi high (3.5); one fold
(basin 0.85, RMSF 0.55); chameleonic. → **open, donor-exposed, single fold — same character as R-8-12.**

**3-12-10-12 S (sarcosine):** backbone IMHB high (3.0); SA_HD lower (68); amphi low (1.9); **rigid**
(RMSF 0.44, single fold 0.98); high asphericity. *Strongest cross-family agreement of all four.*
→ **a single rigid, closed, pre-organized fold ("locked").**

---

## 3. Epimer contrast — what holds across both scaffolds

**Robust epimer trend (both pairs):** **R = open / donors-exposed / polar-face-out / chameleonic;
S = closed / donors-buried.** Because the H-bond, surface, *and* cross-solvent families independently
agree in both pairs, this is real signal — the stereocenter sets open(R) vs closed(S). 2D descriptors
and TPSA are identical between epimers, but the **solvent-accessible 3D-PSA** (donor-H-inclusive,
Ono/Begnini definition) *does* separate them (R-exposed vs S-buried, ~21–31%), consistent with
Begnini 2021. `p_dominant`/`n_eff`/`psa_spread` are discretization/sampling artifacts (fail the
internal-consistency test) and are not used.

---

## 4. Building blocks: local vs. global rigidity (the corrected insight)

The scaffolds differ specifically in the **S** isomer's global flexibility, and it is **counter** to the
naive "rings rigidify" intuition:

- **Azetidine (8-12):** a **local** rigid element — it **locks a cis-amide turn** (§5) — yet the macrocycle
  is **globally floppy** (S RMSF 0.96, many folds). The forced local turn appears to drive global
  conformational diversity (frustration).
- **Sarcosine (10-12):** a **local** flexible hinge that lets the whole macrocycle settle into **one global
  fold** (S RMSF 0.44, single fold 0.98) → global rigidity.

→ **The flexible residue (sarcosine) produced the more pre-organized macrocycle; the rigid ring
(azetidine) produced the floppier one.** Local constraint ≠ global pre-organization.

---

## 5. Cis-amide check — independent validation of the building blocks

Boltzmann-weighted per-backbone-amide cis probability (`cis_prob`, ω<30°):

- **3-12-8-12 (azetidine): one amide locked at 100% cis** in all states (R, S, water, mem) — exactly the
  proline-homolog cis-induction azetidine is known for. The cis lock is **invariant to the stereocenter**
  (same in R and S) → a fixed structural backdrop while the stereocenter modulates the rest.
- **3-12-10-12 (sarcosine): no locked cis** — consistent with sarcosine's weaker tendency.

So the cis descriptor independently distinguishes Aze vs Sar correctly — a consistency check passed.
(Caveat: direct ω measurement on a fused ring-junction N is fiddly; the locked-cis read for the
azetidine compound is chemistry-consistent, but confirming exactly which residue is bond #3 is a
small TODO before quoting it in a paper.)

---

## 6. Permeability / solubility hypothesis (direction, not magnitude)

| axis | physical driver | descriptor | read |
|---|---|---|---|
| permeability | desolvation in membrane state | SA_HD, IMHB (membrane) | both residues remove a donor → similar |
| permeability | low entropy cost to insert (pre-organization) | RMSF, basin | **sarcosine 10-12 globally rigid → favored; esp. S-10-12** |
| permeability | chameleonicity (closes for membrane) | Δ(water−mem) | R is the bigger switcher |
| solubility | polar/donor exposure in water | SA_HD, amphi (water) | **R > S** (R displays donors/polar face) |

**Predictions (falsifiable):**
- Epimer-wise, **S more permeable, R more soluble** within each pair (pre-organization vs solubility trade-off).
- Scaffold-wise, **3-12-10-12 (sarcosine) leans more permeable** (globally pre-organized), **3-12-8-12
  (azetidine) more soluble** (globally floppy) — with the azetidine's locked cis turn a wildcard (it
  pre-sets one turn that *could* aid permeation if it is the membrane-relevant geometry).
- The single most diagnostic measurement remains the **experimental R-vs-S permeability + kinetic solubility.**

---

## 7. Honest caveats
- Hypothesis-generating: single-start CREST, implicit solvent, sub-threshold 6-mers, no experimental structure.
- Descriptors are not yet calibrated to a permeability scale (that is the ML benchmark's job); this note
  gives mechanism + direction.
- The flexibility metric is **weighted RMSF** (threshold-free) + **1-Å basin clustering**, *not* raw
  `p_dominant`/`n_eff` (discretization-sensitive). Reasoning recorded in the 2026-06-16 p_dominant/CREMP analysis.

*Source data: `results/descriptors_8_factorial.csv`; figures: `results/figures/isomers/`; reports:
`docs/experiments/2026-06-16_dopc_3-12-8-12_isomer_rs_report.md` and `..._3-12-10-12_...md`.*

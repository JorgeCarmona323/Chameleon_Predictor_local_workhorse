# NMR Validation Design + GFN-FF Pivot (Pegasus-informed)

**2026-07-08 · Jorge Carmona (Hu Lab, SDSU) · path-forward / experiment design**

> **Thesis.** Adopt Pegasus's scalable engine (GFN-FF metadynamics → GFN2-xTB/CPCM-X
> single-points) but add the rigor Pegasus lacks: **anchor the ensembles to experimental
> NMR.** Build the validation chain **NMR ≈ xTB ≈ GFN-FF ≈ (implicit)** and map *where it
> holds* (small/rigid) vs *breaks* (large chameleons → explicit water). Then run everything
> on GFN-FF for scale, scoring energies with xTB CPCM-X.

---

## 1. Why (Pegasus, Baker *et al.* JMC 2026, 69, 5175)

Pegasus's physics pipeline: RDKit → GFN-FF opt → **GFN-FF metadynamics** (10×100 ps/solvent,
ALPB **water + hexane**) → uniform-sample 100 conformers → **GFN2-xTB single-points** for
features; headline descriptor **ΔG(hexane)/ΔG(water)** predicts ChromlogD at R = 0.91.
Speedups: ~20× (GFN-FF vs GFN2 MTD), ~110× (implicit vs explicit).

**Validation gap:** Pegasus validates only *computation-vs-computation* (SI Table S2: implicit
vs explicit — CsA IMHB 0.95 vs 0.95, CsE 1.14 vs 0.50, SASA/Rg/ΔGsolv within 2–10% / 1–4
kcal/mol; Tables S3–S5: GFN-FF vs GFN2 speed) and *descriptor-vs-bulk-experiment* (ΔG ratio
vs ChromlogD). **No NMR, no crystal, no conformer-level experimental validation anywhere.**
That is our opening.

## 2. Validation chain

`NMR (experiment) ≈ xTB (GFN2 CREST) ≈ GFN-FF ≈ implicit-at-scale`

- Pegasus proved the middle links (FF≈xTB, implicit≈explicit).
- **We add the bottom anchor (xTB ≈ NMR)** — the credibility Pegasus never established.
- Known break point: passes at 6-mer (HexPep, ³J MAE 0.59 Hz) but fails at 11-mer chameleon
  CsA (A1 cis fold, 0% — needs explicit water). Deliverable = a **map of where the chain
  holds**, not a blanket claim. See [HexPep](2026-07-08_hexpep_nmr_validation.md),
  [CsA v2](2026-07-07_csa_v2_validation.md).

## 3. NMR reference data (from `docs/literature/NMR_Validation_Literature`)

| Source | Molecules | Solvents | Observables → method | Status |
|---|---|---|---|---|
| Rezai 2006 | HexPep cmpd 1 (perm.) + 9 (imperm.) | CDCl₃ | ³J(HN-Hα) tables; NOE counts only; D2O exchange (spectra) | J usable |
| Danelius 2020 | roxithromycin, telithromycin, spiramycin, rifampicin | **water (pH 7) + chloroform** | δ, J, temperature coeff., T1, NOE → **NAMFIS** ensembles (800 MHz) | needs SI (J/NOE tables or conformer coords) |
| Sampling paper (bRo5, Kihlberg grp) | 10 bRo5 drugs | polar + apolar | NMR/NAMFIS solution conformers + X-ray | recovery-RMSD template; needs conformer coords |
| Ono 2019 | cyclic hexapeptides | water/CHCl₃/cyclohexane | MD + SASA (no NMR) | **not** an NMR source; keep as cyclohexane-permeability ref |

Key techniques to emulate: **NAMFIS** (deconvolutes NOE + J + temperature coefficients + T1
into a solution conformer *population*); dual-solvent measurement (polar + apolar) to probe
chameleonicity; conformer recovery vs. NMR/crystal by RMSD.

## 4. Validation experiment (two metrics)

1. **Observable back-calculation** (what we did for HexPep): from the Boltzmann ensemble,
   back-calculate ³J (Karplus), and where restraint lists exist, NOE distances (⟨r⁻⁶⟩) and
   the H-bond pattern (temperature coefficients / exchange); MAE vs experiment. Do it **in the
   matched solvent** and, for chameleons, **in both phases**.
2. **Conformer recovery (RMSD)** — does our ensemble *contain* the published NMR/NAMFIS
   solution conformers (heavy-atom RMSD ≲ 2 Å), **per solvent**? (The sampling-paper metric;
   requires the reference conformer coordinates.)

**Target ladder (size → placement on the validity envelope):**
HexPep 6-mer (Rezai) → Danelius drug macrocycles 15–25-atom rings (dual-solvent NAMFIS) →
CsA 11-mer chameleon. Fill 7–10-mer as data allows.

## 5. HexPep redo (immediate, peptide scaffold)

- **Both diastereomers, CDCl₃-matched** (Rezai NMR is CDCl₃): cmpd 9 (= our HexPep, set_B,
  done 0.59 Hz) **and cmpd 1 (permeable, set_A: Tyr 12.0, Leu 6.6/7.8/8.4/9.0)**.
- The permeable-vs-impermeable **pair** is the discrimination test: does the pipeline separate
  1 from 9 as experiment does?
- **Prerequisites for cmpd 1:** (a) its SMILES — pull the exact D/L pattern from Rezai (the
  paper defines the 1/9 diastereomer series); (b) generate its CREST ensemble (we have cmpd 9
  only). Then back-calculate ³J vs set_A.
- Optional strengthening: NOE/exchange for HexPep is not quantitatively extractable from the
  Rezai SI (counts + spectra only), so ³J carries it unless we obtain the deposited restraints.

## 6. GFN-FF pivot — minimal change to the current pipeline

Our `crest_v3.2` already *is* metadynamics (iMTD-GC) with genetic crossing + cregen dedup —
richer than Pegasus's raw 10×100 ps xTB-MTD. **The only change to sample at FF speed is the
level-of-theory flag:**

- **Sampling:** `crest <seed>.xyz --gfnff --alpb <solvent> ...` (was `--gfn2`). ~20× faster;
  no SCF → no `crest_xtbsp` underflow crashes (the reason we downgraded to CREST 2.12).
- **Scoring (new step):** GFN2-xTB single-points with **CPCM-X** on the sampled conformers →
  ΔGsolv per phase → **ΔG(cyclohexane or hexane)/ΔG(water)** transfer descriptor + the
  validated-core geometry descriptors. (Upgrade on Pegasus, which used ALPB ΔGsolv.)
- **Optional geometry refinement:** a quick GFN2 `--opt` on the ~100 sampled conformers before
  scoring recovers geometry quality lost at FF level (hybrid), for PSA/IMHB-sensitive features.
- **Internal check first (Pegasus Tables S3–S5 analog):** on a few compounds, compare GFN-FF vs
  GFN2 ensembles on the descriptor set (IMHB, SASA, Rg, ΔGsolv, PSA) + wall-clock, to confirm
  FF≈xTB before switching wholesale.

## 7. Next actions

1. **Redo HexPep** in CDCl₃ for cmpd 9 (have) — then obtain cmpd 1 SMILES + generate its
   ensemble → ³J validation of the permeable/impermeable pair.
2. **Pull Danelius + sampling-paper SIs** for NAMFIS conformer coordinates / J-NOE tables →
   enables the dual-solvent recovery-RMSD and back-calculation on drug-sized macrocycles.
3. **GFN-FF vs GFN2 internal comparison** on 3–5 compounds (descriptors + timing).
4. **CPCM-X single-point** scoring script (also try SMD / COSMO) → ΔG_transfer; the
   FlexiSol replication for our system.
5. Colab-compatibility pass on `notebooks/pipeline` (clone repo → input SMILES → pick solvents
   → run) — deferred, tracked separately.

## References

Baker et al., *J. Med. Chem.* 2026, 69, 5175 (+ SI). Rezai & Lokey, *JACS* 2006, 128, 2510.
Danelius et al., *Chem. Eur. J.* 2020 (chameleonic drugs, NAMFIS). Poongavanam/Kihlberg,
conformational-sampling-of-macrocyclic-drugs (bRo5 sampling benchmark). Ono et al., *JCIM*
2019 (cyclohexane SASA). Vuister & Bax, *JACS* 1993 (Karplus). Grimme, GFN-FF, *Angew. Chem.*
2020.

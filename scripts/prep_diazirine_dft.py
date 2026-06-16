# env: chameleon-calc
"""
prep_diazirine_dft.py
---------------------
Generate ORCA DFT inputs for the diazirine N=N calibration (step 5 of the integrity
check). GFN2-xTB relaxes the diazirine N=N to a spurious ~1.43 A in 5/8 ensembles;
DFT on the local diazirine motif tells us the true value (~1.23 A expected) and whether
1.43 A is even a stationary point.

Model: 3-(trifluoromethyl)-3-phenyl-3H-diazirine (TPD) — the canonical photoaffinity
diazirine. The N=N geometry is local, so this ~21-atom fragment answers the question at
a fraction of the cost of the full ~150-atom macrocycle. Aryl/alkyne handle truncated to
phenyl (far from N=N, no effect on the ring geometry).

Writes to dft/diazirine_calibration/:
  orca_diazirine_opt.inp   geometry optimization + frequencies (true-minimum check)
  orca_diazirine_scan.inp  relaxed N=N scan 1.15->1.50 A (is 1.43 a real basin?)

Run:  python scripts/prep_diazirine_dft.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

FRAG_SMILES = "FC(F)(F)C1(c2ccccc2)N=N1"      # 3-(trifluoromethyl)-3-phenyl-3H-diazirine
METHOD = "wB97X-D3 def2-TZVP TightSCF"         # robust geometry; swap if colleague prefers
NPROCS, MAXCORE = 8, 3000
OUTDIR = Path("dft/diazirine_calibration")
# Match our two ensemble solvents (ALPB water / ALPB chloroform) + gas reference.
# N=N is a local bond, so all three should agree (~1.23 A) — this rules out a solvent artifact.
SOLVENTS = [("gas", None), ("water", "water"), ("chcl3", "chloroform")]


def build():
    m = Chem.AddHs(Chem.MolFromSmiles(FRAG_SMILES))
    AllChem.EmbedMolecule(m, randomSeed=1)
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    return m


def xyz_block(m, coords=None) -> str:
    conf = m.GetConformer()
    out = []
    for a in m.GetAtoms():
        if coords is None:
            p = conf.GetAtomPosition(a.GetIdx())
            x, y, z = p.x, p.y, p.z
        else:
            x, y, z = coords[a.GetIdx()]
        out.append(f"{a.GetSymbol():2s} {x:14.8f} {y:14.8f} {z:14.8f}")
    return "\n".join(out)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    m = build()
    c, n1, n2 = m.GetSubstructMatches(Chem.MolFromSmarts("[#6]1[#7]=[#7]1"))[0]
    xyz = xyz_block(m)
    p = m.GetConformer().GetPositions()
    nn0 = float(np.linalg.norm(p[n1] - p[n2]))

    header = (f"# {FRAG_SMILES}  (3-CF3-3-phenyl-3H-diazirine, TPD model)\n"
              f"# diazirine N atoms = 0-based indices {n1} and {n2} in the geometry below\n"
              f"# GFN2-xTB gives this N=N ~1.43 A; DFT expected ~1.23 A. Start (MMFF) = {nn0:.3f} A\n")

    written = []
    for label, solv in SOLVENTS:
        cpcm = f" CPCM({solv})" if solv else ""
        freq = " Freq" if solv is None else ""   # one gas-phase freq is enough for the minimum check
        note = f"# solvent: {solv or 'gas phase'}\n"
        opt = (f"! {METHOD} Opt{freq}{cpcm}\n"
               f"%maxcore {MAXCORE}\n%pal nprocs {NPROCS} end\n{header}{note}"
               f"* xyz 0 1\n{xyz}\n*\n")
        fn = f"orca_diazirine_opt_{label}.inp"
        (OUTDIR / fn).write_text(opt, encoding="utf-8")
        written.append(fn)

    # opt STARTING from the GFN2-distorted 1.43 A geometry: if DFT relaxes it back to ~1.23 A,
    # the optimized XYZ alone proves 1.43 A was a spurious basin (no energies / scan needed).
    pert = m.GetConformer().GetPositions().copy()
    v = pert[n2] - pert[n1]
    pert[n2] = pert[n1] + (v / np.linalg.norm(v)) * 1.43
    opt143 = (f"! {METHOD} Opt Freq\n"
              f"%maxcore {MAXCORE}\n%pal nprocs {NPROCS} end\n{header}"
              f"# START geometry has N=N forced to 1.43 A (GFN2's value). Does DFT relax it to ~1.23?\n"
              f"* xyz 0 1\n{xyz_block(m, pert)}\n*\n")
    (OUTDIR / "orca_diazirine_opt_from143_gas.inp").write_text(opt143, encoding="utf-8")
    written.append("orca_diazirine_opt_from143_gas.inp")

    # OPTIONAL relaxed N=N scan (gas) — only if the full energy curve is wanted; needs energies, not XYZ
    scan = (f"! {METHOD} Opt\n"
            f"%maxcore {MAXCORE}\n%pal nprocs {NPROCS} end\n{header}"
            f"# relaxed scan of N=N 1.15 -> 1.50 A (15 pts): is 1.43 A a real minimum?\n"
            f"%geom Scan\n  B {n1} {n2} = 1.15, 1.50, 15\n  end\nend\n"
            f"* xyz 0 1\n{xyz}\n*\n")
    (OUTDIR / "orca_diazirine_scan.inp").write_text(scan, encoding="utf-8")
    written.append("orca_diazirine_scan.inp")

    print(f"Wrote ORCA inputs -> {OUTDIR}")
    print(f"  start N=N (MMFF) = {nn0:.3f} A | N indices {n1},{n2} | atoms={m.GetNumAtoms()}")
    for fn in written:
        print(f"  {fn}")


if __name__ == "__main__":
    main()

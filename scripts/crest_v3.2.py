# env: chameleon-sim
"""
crest_v3.2.py
-------------
Entry point for the CREST iMTD-GC conformer sampling pipeline.
Defines reference compounds and dispatches to crest_conformers.py.

Usage:
  python crest_v3.2.py --compound 1 --threads 20 --outdir results
  python crest_v3.2.py --compound 4 --threads 20 --outdir results --resume

Compound index:
  0 = HexPep   (impermeable,  6-mer)
  1 = CsA      (permeable,   11-mer)
  2 = CsO      (permeable,   11-mer)
  3 = PSLYF    (impermeable, 11-mer)
  4 = WhC3     (permeable,    6-mer, White 2011 compd.3, RRCK=-5.31)
  5 = DOPC_R   (permeable,    6-mer, DOPC 3-12-8-12 R isomer)
  6 = DOPC_S   (permeable,    6-mer, DOPC 3-12-8-12 S isomer)
  7 = Brain1   (permeable,    6-mer, Brain 6-4-4-13, naphthalene)
  8 = DOPC2    (permeable,    6-mer, DOPC 6-5-8-12, naphthalene)
  9 = CsA_v2  (permeable,   11-mer, CsA rerun with --noreftopo -notopo)
 10 = DOPCsar_R   (DOPC 3-12-10-12 R, sarcosine for azetidine)
 11 = DOPCsar_S   (DOPC 3-12-10-12 S, sarcosine for azetidine)
 12 = DOPCdz_R    (DOPC 3-12-8-12 R + CF3-diazirine)
 13 = DOPCdz_S    (DOPC 3-12-8-12 S + CF3-diazirine)
 14 = DOPCsardz_R (DOPC 3-12-10-12 R sarcosine + CF3-diazirine)
 15 = DOPCsardz_S (DOPC 3-12-10-12 S sarcosine + CF3-diazirine)
"""

from __future__ import annotations

import argparse
import os
from datetime import datetime
from pathlib import Path

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import pandas as pd

from crest_conformers import process_compound, find_resume_dir

# ── Reference compounds ───────────────────────────────────────────────────────
REFERENCE_COMPOUNDS = [
    {
        "name": "Hexapeptide",
        "short": "HexPep",
        "cycpeptmpdb_id": 2,
        "smiles": (
            "CC(C)C[C@@H]1NC(=O)[C@@H](CC(C)C)NC(=O)[C@@H](CC(C)C)NC(=O)"
            "[C@H](Cc2ccc(O)cc2)NC(=O)[C@@H]2CCCN2C(=O)[C@@H](CC(C)C)NC1=O"
        ),
        "source": "Rezai & Lokey, JACS 2006",
        "pampa": -6.20,
        "permeable": False,
        "hbd": 6,
    },
    {
        "name": "Cyclosporin A",
        "short": "CsA",
        "cycpeptmpdb_id": 1,
        "smiles": (
            "C/C=C/C[C@@H](C)[C@@H](O)[C@H]1C(=O)N[C@@H](CC)C(=O)N(C)CC(=O)"
            "N(C)[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N(C)[C@@H](CC(C)C)C(=O)"
            "N[C@@H](C)C(=O)N[C@H](C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)[C@@H](CC(C)C)"
            "C(=O)N(C)[C@@H](C(C)C)C(=O)N1C"
        ),
        "source": "Witek JCTC 2016",
        "pampa": -5.90,
        "permeable": True,
        "hbd": 5,
    },
    {
        "name": "Cyclosporin O",
        "short": "CsO",
        "cycpeptmpdb_id": None,
        "smiles": (
            "CCC[C@H]1C(=O)N(CC(=O)N([C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H]"
            "(C(=O)N[C@@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H]"
            "(C(=O)N1)CC(C)C)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)"
            "CC(C)C)C)C"
        ),
        "source": "Horizon-LBA Ono et al. Chem. Sci. 2023; LPE Naylor et al. J. Med. Chem. 2018",
        "pampa": None,
        "permeable": True,
        "hbd": 4,
        "horizon_lba_papp": 3e-6,
    },
    {
        "name": "c*[PSLYF]",
        "short": "PSLYF",
        "cycpeptmpdb_id": 1829,
        "smiles": (
            "CC(C)C[C@@H]1NC(=O)[C@H](CO)NC(=O)[C@@H]2CCCN2[C@H](C(=O)NC(C)(C)C)"
            "[C@H](C)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](Cc2ccc(O)cc2)NC1=O"
        ),
        "source": "Hickey, J Med Chem 2016",
        "pampa": -9.10,
        "permeable": False,
        "hbd": 8,
    },
    {
        "name": "White_compd3",
        "short": "WhC3",
        "cycpeptmpdb_id": 25,
        "smiles": (
            "CC(C)C[C@@H]1NC(=O)[C@H](Cc2ccc(O)cc2)N(C)C(=O)[C@H]2CCCN2C(=O)"
            "[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)N(C)C(=O)[C@@H](CC(C)C)N(C)C1=O"
        ),
        "source": "White, Nat Chem Biol 2011",
        "pampa": -5.31,
        "permeable": True,
        "hbd": 3,
    },
    {
        "name": "DOPC_3-12-8-12_R",
        "short": "DOPC_R",
        "cycpeptmpdb_id": None,
        "smiles": (
            "C#CCCC(=O)N[C@H]1CSCc2ccccc2CSC[C@@H](C(N)=O)NC(=O)[C@H](c2cccs2)"
            "NC(=O)[C@H](CO)NC(=O)[C@@H]2CCN2C(=O)[C@H](CO)NC1=O"
        ),
        "source": "Hu lab — DOPC 3-12-8-12 R isomer",
        "pampa": None,
        "permeable": True,
        "hbd": None,
    },
    {
        "name": "DOPC_3-12-8-12_S",
        "short": "DOPC_S",
        "cycpeptmpdb_id": None,
        "smiles": (
            "C#CCCC(=O)N[C@H]1CSCc2ccccc2CSC[C@@H](C(N)=O)NC(=O)[C@@H](c2cccs2)"
            "NC(=O)[C@H](CO)NC(=O)[C@@H]2CCN2C(=O)[C@H](CO)NC1=O"
        ),
        "source": "Hu lab — DOPC 3-12-8-12 S isomer",
        "pampa": None,
        "permeable": True,
        "hbd": None,
    },
    {
        "name": "Brain_6-4-4-13",
        "short": "Brain1",
        "cycpeptmpdb_id": None,
        "smiles": (
            "C#CCCC(=O)N[C@H]1CSCc2ccccc2CSC[C@@H](C(N)=O)NC(=O)C[C@@H](c2cccc3ccccc23)"
            "NC(=O)[C@H](C)NC(=O)[C@@H]2CCCCN2C(=O)[C@@H](CO)NC1=O"
        ),
        "source": "Hu lab — Brain 6-4-4-13 (naphthalene)",
        "pampa": None,
        "permeable": True,
        "hbd": None,
    },
    {
        "name": "DOPC_6-5-8-12",
        "short": "DOPC2",
        "cycpeptmpdb_id": None,
        "smiles": (
            "C#CCCC(=O)N[C@H]1CSCc2ccccc2CSC[C@@H](C(N)=O)NC(=O)C[C@@H](c2cccc3ccccc23)"
            "NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@@H]2CCN2C(=O)[C@H](CO)NC1=O"
        ),
        "source": "Hu lab — DOPC 6-5-8-12 (naphthalene)",
        "pampa": None,
        "permeable": True,
        "hbd": None,
    },
    {
        "name": "Cyclosporin A",
        "short": "CsA_v2",
        "cycpeptmpdb_id": 1,
        "smiles": (
            "C/C=C/C[C@@H](C)[C@@H](O)[C@H]1C(=O)N[C@@H](CC)C(=O)N(C)CC(=O)"
            "N(C)[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N(C)[C@@H](CC(C)C)C(=O)"
            "N[C@@H](C)C(=O)N[C@H](C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)[C@@H](CC(C)C)"
            "C(=O)N(C)[C@@H](C(C)C)C(=O)N1C"
        ),
        "source": "Witek JCTC 2016 — rerun with --noreftopo -notopo",
        "pampa": -5.90,
        "permeable": True,
        "hbd": 5,
    },
    {
        "name": "DOPC_3-12-10-12_R",
        "short": "DOPCsar_R",
        "cycpeptmpdb_id": None,
        "smiles": "CN(CC(N[C@@H](CO)C(N[C@@H](c1ccc[s]1)C(N[C@@H](CSCc1c(CSC[C@@H](C(N[C@@H]2CO)=O)NC(CCC#C)=O)cccc1)C(N)=O)=O)=O)=O)C2=O",
        "source": "Hu lab — DOPC 3-12-10-12 R (sarcosine for azetidine)",
        "pampa": None, "permeable": True, "hbd": None,
    },
    {
        "name": "DOPC_3-12-10-12_S",
        "short": "DOPCsar_S",
        "cycpeptmpdb_id": None,
        "smiles": "CN(CC(N[C@@H](CO)C(N[C@H](c1ccc[s]1)C(N[C@@H](CSCc1c(CSC[C@@H](C(N[C@@H]2CO)=O)NC(CCC#C)=O)cccc1)C(N)=O)=O)=O)=O)C2=O",
        "source": "Hu lab — DOPC 3-12-10-12 S (sarcosine for azetidine)",
        "pampa": None, "permeable": True, "hbd": None,
    },
    {
        "name": "DOPC_3-12-8-12_R_diazirine",
        "short": "DOPCdz_R",
        "cycpeptmpdb_id": None,
        "smiles": "C#CCCC(N[C@@H](CSCc1cc(C2(C(F)(F)F)N=N2)cc(CSC[C@@H](C(N)=O)NC([C@H](c2ccc[s]2)NC([C@H](CO)NC([C@H](CC2)N2C([C@H](CO)N2)=O)=O)=O)=O)c1)C2=O)=O",
        "source": "Hu lab — DOPC 3-12-8-12 R + CF3-diazirine photocrosslinker",
        "pampa": None, "permeable": True, "hbd": None,
    },
    {
        "name": "DOPC_3-12-8-12_S_diazirine",
        "short": "DOPCdz_S",
        "cycpeptmpdb_id": None,
        "smiles": "C#CCCC(N[C@@H](CSCc1cc(C2(C(F)(F)F)N=N2)cc(CSC[C@@H](C(N)=O)NC([C@@H](c2ccc[s]2)NC([C@H](CO)NC([C@H](CC2)N2C([C@H](CO)N2)=O)=O)=O)=O)c1)C2=O)=O",
        "source": "Hu lab — DOPC 3-12-8-12 S + CF3-diazirine photocrosslinker",
        "pampa": None, "permeable": True, "hbd": None,
    },
    {
        "name": "DOPC_3-12-10-12_R_diazirine",
        "short": "DOPCsardz_R",
        "cycpeptmpdb_id": None,
        "smiles": "CN(CC(N[C@@H](CO)C(N[C@@H](c1ccc[s]1)C(N[C@@H](CSCc1cc(CSC[C@@H](C(N[C@@H]2CO)=O)NC(CCC#C)=O)cc(C3(C(F)(F)F)N=N3)c1)C(N)=O)=O)=O)=O)C2=O",
        "source": "Hu lab — DOPC 3-12-10-12 R (sarcosine) + CF3-diazirine",
        "pampa": None, "permeable": True, "hbd": None,
    },
    {
        "name": "DOPC_3-12-10-12_S_diazirine",
        "short": "DOPCsardz_S",
        "cycpeptmpdb_id": None,
        "smiles": "CN(CC(N[C@@H](CO)C(N[C@H](c1ccc[s]1)C(N[C@@H](CSCc1cc(CSC[C@@H](C(N[C@@H]2CO)=O)NC(CCC#C)=O)cc(C3(C(F)(F)F)N=N3)c1)C(N)=O)=O)=O)=O)C2=O",
        "source": "Hu lab — DOPC 3-12-10-12 S (sarcosine) + CF3-diazirine",
        "pampa": None, "permeable": True, "hbd": None,
    },
]


def run(outdir: Path, max_confs: int | None, dry_run: bool,
        compound_idx: int | None = None, n_threads: int | None = None,
        resume: bool = False) -> None:

    n_threads = n_threads or os.cpu_count() or 1

    if compound_idx is None:
        raise RuntimeError("--compound is required.")
    if compound_idx < 0 or compound_idx >= len(REFERENCE_COMPOUNDS):
        raise ValueError(f"--compound must be 0–{len(REFERENCE_COMPOUNDS)-1}")

    cpd   = REFERENCE_COMPOUNDS[compound_idx]
    short = cpd["short"]

    resume_dir = find_resume_dir(outdir / "runs", compound_idx, short) if resume else None
    if resume_dir is not None:
        work_base = resume_dir
        print(f"\n[Compound {compound_idx}] {cpd['name']}  ← resuming")
        print(f"Run directory: {work_base}")
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_id    = f"run_{timestamp}_{compound_idx}_{short}"
        work_base = outdir / "runs" / run_id
        work_base.mkdir(parents=True, exist_ok=True)
        print(f"\n[Compound {compound_idx}] {cpd['name']}")
        print(f"Run directory: {work_base}")

    r = process_compound(cpd, work_base, max_confs, dry_run, n_threads)

    out_csv = work_base / f"{short}_results.csv"
    pd.DataFrame([r]).to_csv(out_csv, index=False)
    print(f"\nSaved: {out_csv}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Tier-2 CREST+ALPB validation (v3.2 — CREMP-format outputs)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--outdir",    "-o", default="results", type=Path)
    parser.add_argument("--compound",  type=int, default=None, metavar="IDX")
    parser.add_argument("--threads",   type=int, default=None, metavar="N")
    parser.add_argument("--max-confs", "-c", type=int, default=None,
                        help="Cap conformers used for PSA analysis. Default: 50.")
    parser.add_argument("--dry-run",   action="store_true")
    parser.add_argument("--resume",    action="store_true",
                        help="Resume a previous incomplete run instead of starting fresh.")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(
        Path(args.outdir),
        max_confs    = args.max_confs,
        dry_run      = args.dry_run,
        compound_idx = args.compound,
        n_threads    = args.threads,
        resume       = args.resume,
    )

# env: chameleon-calc
"""
reference_compounds.py
----------------------
Single source of truth for the reference compound set. Both crest_v3.2.py and
submit_tier2_slurm.py load from data/reference_set.csv via load() here, so adding
a compound means appending one CSV row -- no more editing three files.

The CSV path data/reference_set.csv is git-tracked (allowlisted in .gitignore
despite data/ being ignored), so it ships to the HPC with every git pull.

  load()              -> list[dict] matching the old REFERENCE_COMPOUNDS structure
  export(compounds)   -> write a compound list to the CSV (bootstrap / regenerate)

Run `python scripts/reference_compounds.py` to validate the CSV round-trips.
"""
from __future__ import annotations

import csv
from pathlib import Path

_DEFAULT_CSV = Path(__file__).resolve().parent.parent / "data" / "reference_set.csv"
COLUMNS = ["short", "name", "cycpeptmpdb_id", "smiles", "source",
           "pampa", "permeable", "hbd", "horizon_lba_papp"]


def _coerce(val, typ):
    if val is None or str(val).strip() == "":
        return None
    if typ is bool:
        return str(val).strip().lower() in ("true", "1", "yes")
    return typ(val)


def load(csv_path: str | Path | None = None) -> list[dict]:
    """Read the master CSV -> list of compound dicts (same shape as the old hardcoded list)."""
    path = Path(csv_path) if csv_path else _DEFAULT_CSV
    out: list[dict] = []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            cpd = {
                "name":           row["name"],
                "short":          row["short"],
                "cycpeptmpdb_id": _coerce(row.get("cycpeptmpdb_id"), int),
                "smiles":         row["smiles"],
                "source":         row.get("source") or "",
                "pampa":          _coerce(row.get("pampa"), float),
                "permeable":      _coerce(row.get("permeable"), bool),
                "hbd":            _coerce(row.get("hbd"), int),
            }
            hl = _coerce(row.get("horizon_lba_papp"), float)
            if hl is not None:
                cpd["horizon_lba_papp"] = hl
            out.append(cpd)
    return out


def export(compounds: list[dict], csv_path: str | Path | None = None) -> Path:
    """Write a compound list to the master CSV (used to bootstrap / regenerate)."""
    path = Path(csv_path) if csv_path else _DEFAULT_CSV
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(COLUMNS)
        for c in compounds:
            row = []
            for k in COLUMNS:
                v = c.get(k)
                row.append("" if v is None else v)
            w.writerow(row)
    return path


if __name__ == "__main__":
    # Validate the CSV round-trips and report the set.
    rc = load()
    print(f"{len(rc)} compounds in {_DEFAULT_CSV.name}:")
    for i, c in enumerate(rc):
        print(f"  {i:2d}  {c['short']:12s}  {c['name']}")

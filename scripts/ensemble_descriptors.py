"""
ensemble_descriptors.py
-----------------------
Compute ML-ready 3D conformational descriptors from completed CREST ensembles.
Featurizer for the DynamicEnsembleEncoder (see docs/chameleon_model_architecture.md);
output columns map to the dynamic_features schema (docs/data_schema.md).

Operates directly on the two solvent ensembles — water/ and mem/ — using the
solvent label as the state assignment (no clustering, no max/min-PSA proxy).
Per-solvent descriptors are Boltzmann-weighted over the real CREST ensemble;
cross-solvent descriptors are water − CHCl3 differences.

Descriptor groups (see docs/descriptor_framework.md): 2 (cis-amide), 3 (ddG),
5 (Boltzmann polarity), 6 (shape). Witek congruent + kinetic barrier intentionally
omitted (see docs/ml_descriptor_implications.md, 2026-05-31 decision).

Usage:
  python ensemble_descriptors.py --run-dir results/runs/run_..._5_DOPC_R --name DOPC_R
  python ensemble_descriptors.py --water-dir data/CREST_CsA_20260512 --name CsA_v1
  python ensemble_descriptors.py --run-dir <dir> --name X --run-dir <dir2> --name Y -o out.csv
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors3D, rdFreeSASA

from phys_descriptors_v2 import compute_psa_xyz, count_hbonds_xyz, boltzmann_weights

RDLogger.DisableLog("rdApp.*")

RT_KCAL = 0.592  # kcal/mol at 298.15 K
CIS_OMEGA_CUTOFF = 30.0  # |omega| < 30 deg → cis


# ── geometry helpers ──────────────────────────────────────────────────────────
def _dihedral(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / (np.linalg.norm(b1) + 1e-12)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def amide_ring_bonds(mol) -> list[tuple]:
    """Return (Ca_acyl, C, N, Ca_amino) atom-index tuples for backbone amide bonds
    that lie in the macrocycle ring. Omega dihedral over these = cis/trans state."""
    patt = Chem.MolFromSmarts("[C;X3](=O)[N;X3]")
    out = []
    for match in mol.GetSubstructMatches(patt):
        c, o, n = match
        bond = mol.GetBondBetweenAtoms(c, n)
        if bond is None or not bond.IsInRing():
            continue
        ca_acyl = None
        for nb in mol.GetAtomWithIdx(c).GetNeighbors():
            j = nb.GetIdx()
            if j in (o, n):
                continue
            b = mol.GetBondBetweenAtoms(c, j)
            if nb.GetSymbol() == "C" and b is not None and b.IsInRing():
                ca_acyl = j
                break
        ca_amino = None
        for nb in mol.GetAtomWithIdx(n).GetNeighbors():
            j = nb.GetIdx()
            if j == c:
                continue
            b = mol.GetBondBetweenAtoms(n, j)
            if nb.GetSymbol() == "C" and b is not None and b.IsInRing():
                ca_amino = j
                break
        if ca_acyl is not None and ca_amino is not None:
            out.append((ca_acyl, c, n, ca_amino))
    return out


def total_sasa(mol, conf_id) -> float:
    try:
        radii = rdFreeSASA.classifyAtoms(mol)
        return float(rdFreeSASA.CalcSASA(mol, radii, confIdx=conf_id))
    except Exception:
        return float("nan")


# ── ensemble loading ──────────────────────────────────────────────────────────
def load_ensemble(sdf_path: Path, json_path: Path):
    """Return dict with mol-per-conformer list, weights, energies, and json psa/hb."""
    with open(json_path) as f:
        data = json.load(f)
    confs = data.get("conformers", [])
    energies = np.array([c.get("totalenergy", np.nan) for c in confs], dtype=float)
    weights = np.array([c.get("boltzmannweight", np.nan) for c in confs], dtype=float)
    if not np.all(np.isfinite(weights)) or weights.sum() <= 0:
        weights = boltzmann_weights(energies.tolist())
    psa_json = np.array([c.get("psa", np.nan) for c in confs], dtype=float)
    hb_json = np.array([c.get("hbonds", np.nan) for c in confs], dtype=float)

    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
    mols = [m for m in supplier if m is not None]

    n = min(len(mols), len(confs))
    return {
        "mols": mols[:n],
        "weights": weights[:n] / weights[:n].sum(),
        "energies": energies[:n],
        "psa_json": psa_json[:n],
        "hb_json": hb_json[:n],
        "smiles": data.get("smiles"),
        "n": n,
    }


# ── per-solvent descriptors ───────────────────────────────────────────────────
def solvent_descriptors(sdf_path: Path, json_path: Path, prefix: str) -> dict:
    ens = load_ensemble(sdf_path, json_path)
    mols, w = ens["mols"], ens["weights"]
    if len(mols) == 0:
        return {}

    amide_tuples = amide_ring_bonds(mols[0])
    n_amide = len(amide_tuples)

    rg, npr1, npr2, asph, sphe, tsasa = [], [], [], [], [], []
    cis_flags = np.zeros((len(mols), n_amide), dtype=float)

    for i, m in enumerate(mols):
        cid = m.GetConformer().GetId()
        try:
            rg.append(Descriptors3D.RadiusOfGyration(m, confId=cid))
            npr1.append(Descriptors3D.NPR1(m, confId=cid))
            npr2.append(Descriptors3D.NPR2(m, confId=cid))
            asph.append(Descriptors3D.Asphericity(m, confId=cid))
            sphe.append(Descriptors3D.SpherocityIndex(m, confId=cid))
        except Exception:
            rg.append(np.nan); npr1.append(np.nan); npr2.append(np.nan)
            asph.append(np.nan); sphe.append(np.nan)
        tsasa.append(total_sasa(m, cid))
        coords = m.GetConformer().GetPositions()
        for j, tup in enumerate(amide_tuples):
            omega = _dihedral(*[coords[k] for k in tup])
            cis_flags[i, j] = 1.0 if abs(omega) < CIS_OMEGA_CUTOFF else 0.0

    def bw(arr):
        arr = np.asarray(arr, dtype=float)
        mask = np.isfinite(arr)
        if not mask.any():
            return float("nan")
        ww = w[mask] / w[mask].sum()
        return float(np.dot(ww, arr[mask]))

    psa = ens["psa_json"]
    hb = ens["hb_json"]
    cis_prob = (w[:, None] * cis_flags).sum(axis=0)  # per amide bond

    # dominant conformer = highest Boltzmann weight
    dom = int(np.argmax(w))

    out = {
        f"{prefix}_n_confs": ens["n"],
        f"{prefix}_bw_psa": round(bw(psa), 2),
        f"{prefix}_bw_hb": round(bw(hb), 3),
        f"{prefix}_bw_rg": round(bw(rg), 3),
        f"{prefix}_bw_npr1": round(bw(npr1), 4),
        f"{prefix}_bw_npr2": round(bw(npr2), 4),
        f"{prefix}_bw_asphericity": round(bw(asph), 4),
        f"{prefix}_bw_spherocity": round(bw(sphe), 4),
        f"{prefix}_psa_std": round(float(np.nanstd(psa)), 2),
        f"{prefix}_psa_spread": round(float(np.nanmax(psa) - np.nanmin(psa)), 2),
        f"{prefix}_sasa_total": round(bw(tsasa), 1),
        f"{prefix}_p_dominant": round(float(w[dom]), 3),
        f"{prefix}_n_amide": n_amide,
    }
    for j in range(n_amide):
        out[f"{prefix}_cis_prob_{j}"] = round(float(cis_prob[j]), 3)
    # store internal for cross-solvent
    out[f"_{prefix}_cis_prob_vec"] = cis_prob.tolist()
    out[f"_{prefix}_dom_energy"] = float(ens["energies"][dom]) if np.isfinite(ens["energies"][dom]) else None
    return out


def cross_solvent(water: dict, mem: dict) -> dict:
    out = {}
    pairs = [
        ("delta_psa", "bw_psa"), ("delta_hb", "bw_hb"), ("delta_rg", "bw_rg"),
        ("delta_npr1", "bw_npr1"), ("delta_npr2", "bw_npr2"),
        ("delta_asphericity", "bw_asphericity"),
    ]
    for out_key, feat in pairs:
        wv, mv = water.get(f"water_{feat}"), mem.get(f"mem_{feat}")
        if wv is not None and mv is not None:
            out[out_key] = round(wv - mv, 3)

    # normalized delta PSA (Yu 2026): ΔPSA / mean total SASA
    if "delta_psa" in out:
        sasa = np.nanmean([water.get("water_sasa_total", np.nan),
                           mem.get("mem_sasa_total", np.nan)])
        if np.isfinite(sasa) and sasa > 0:
            out["norm_delta_psa"] = round(out["delta_psa"] / sasa, 5)

    # ddG between dominant conformers across solvents (Hartree → kcal)
    ew, em = water.get("_water_dom_energy"), mem.get("_mem_dom_energy")
    if ew is not None and em is not None:
        out["ddG_dom_kcal"] = round((ew - em) * 627.509, 2)

    # cis-amide switch: which bond changes most between solvents
    cw = water.get("_water_cis_prob_vec")
    cm = mem.get("_mem_cis_prob_vec")
    if cw and cm and len(cw) == len(cm):
        d = np.array(cw) - np.array(cm)
        mag = float(np.max(np.abs(d)))
        out["cis_switch_mag"] = round(mag, 3)
        # only report a switching bond if a bond actually changes cis-state between solvents
        out["cis_switch_bond"] = int(np.argmax(np.abs(d))) if mag > 1e-6 else None
        # cis entropy in water (Shannon over bonds with any cis population)
        p = np.array(cw)
        p = p[p > 1e-6]
        if p.size:
            p = p / p.sum()
            out["cis_entropy_water"] = round(float(-(p * np.log(p)).sum()), 3)
    return out


def process_compound(name: str, water_dir: Path | None, mem_dir: Path | None) -> dict:
    row = {"compound": name}
    wd = {}
    md = {}
    if water_dir and (water_dir / "ensemble.json").exists():
        wd = solvent_descriptors(water_dir / "ensemble.sdf", water_dir / "ensemble.json", "water")
        row.update({k: v for k, v in wd.items() if not k.startswith("_")})
    if mem_dir and (mem_dir / "ensemble.json").exists():
        md = solvent_descriptors(mem_dir / "ensemble.sdf", mem_dir / "ensemble.json", "mem")
        row.update({k: v for k, v in md.items() if not k.startswith("_")})
    if wd and md:
        row.update(cross_solvent(wd, md))
        row["has_both_solvents"] = 1
    else:
        row["has_both_solvents"] = 0
    return row


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", action="append", default=[],
                   help="CREST run dir containing water/ and mem/ subdirs")
    p.add_argument("--water-dir", action="append", default=[],
                   help="Directory holding a water-only ensemble.{json,sdf}")
    p.add_argument("--name", action="append", default=[], help="Compound name (one per --run-dir/--water-dir)")
    p.add_argument("-o", "--out", default="results/ensemble_descriptors.csv", type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    jobs = []
    for i, rd in enumerate(args.run_dir):
        nm = args.name[i] if i < len(args.name) else Path(rd).name
        jobs.append((nm, Path(rd) / "water", Path(rd) / "mem"))
    offset = len(args.run_dir)
    for i, wdir in enumerate(args.water_dir):
        nm = args.name[offset + i] if offset + i < len(args.name) else Path(wdir).name
        jobs.append((nm, Path(wdir), None))

    rows = []
    for nm, wdir, mdir in jobs:
        print(f"[{nm}] water={wdir} mem={mdir}")
        rows.append(process_compound(nm, wdir, mdir))

    df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"\nSaved {len(df)} compound(s) → {args.out}")
    # console summary
    cols = [c for c in df.columns if not c.startswith("_")]
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(df[cols].T)


if __name__ == "__main__":
    main()

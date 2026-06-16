# env: chameleon-calc
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

v3 additions (2026-06-14; populations stay GFN2-xTB ENERGY-weighted, i.e. CREMP-consistent).
Per docs/experiments/2026-06-13_descriptor_literature_review.md:
  H-bonds (recomputed from geometry; renamed from bw_hb): bw_IMHB (total), bw_IMHBD/bw_IMHBA
    (distinct donors/acceptors engaged), bw_IMHB_bb/bw_IMHB_res (backbone-transannular vs
    side-chain). bb + res == IMHB.
  surface (phys_descriptors_v3): bw_SA_HD (donor-H surface, was bw_hbd_sasa), bw_SA_HA
    (acceptor-atom surface), bw_hydrophobic_sasa, bw_amphi_moment
  regime-1 ensemble fluctuation: std_rg, std_asphericity, std_amphi_moment (Boltzmann-
    weighted std = configurational fluctuation amplitude), omega_circvar (backbone
    cis/trans circular variance), n_eff (effective conformer count, nConf20 analog)
Cross-solvent deltas mirror these (delta_IMHB, delta_IMHB_bb/res, delta_SA_HD, delta_SA_HA, ...).
NOTE: bw_hb -> bw_IMHB and bw_hbd_sasa -> bw_SA_HD are renames; update any consumer that
read the old names (done: plot_isomer_comparison.py).
Regime 2 (Hessian/CENSO free-energy reweighting) is deliberately NOT applied here — it
would change populations off the CREMP footing; see the literature-review doc.

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

from phys_descriptors_v3 import (
    compute_psa_xyz, count_hbonds_xyz, boltzmann_weights,
    surface_descriptors_mol, effective_nconf,
    imhb_descriptors_mol, backbone_hbond_atoms,
)

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
    backbone_atoms = backbone_hbond_atoms(mols[0])  # ring + carbonyl O's = backbone for IMHB

    rg, npr1, npr2, asph, sphe, tsasa = [], [], [], [], [], []
    hbd_l, hba_l, hydro_l, amphi_l = [], [], [], []        # v3 surface descriptors
    imhb_l, imhbd_l, imhba_l, imhb_bb_l, imhb_res_l = [], [], [], [], []  # IMHB breakdown
    cis_flags = np.zeros((len(mols), n_amide), dtype=float)
    omega_vals = np.full((len(mols), n_amide), np.nan, dtype=float)  # regime-1 circular variance

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
        sd = surface_descriptors_mol(m, cid)    # v3: SA_HD, SA_HA, hydrophobic SASA, amphi moment
        hbd_l.append(sd["hbd_sasa"]); hba_l.append(sd["hba_sasa"])
        hydro_l.append(sd["hydrophobic_sasa"])
        amphi_l.append(sd["amphi_moment"])
        ih = imhb_descriptors_mol(m, cid, backbone_atoms)   # IMHB total + donor/acceptor + bb/res
        imhb_l.append(ih["imhb"]); imhbd_l.append(ih["imhbd"]); imhba_l.append(ih["imhba"])
        imhb_bb_l.append(ih["imhb_bb"]); imhb_res_l.append(ih["imhb_res"])
        coords = m.GetConformer().GetPositions()
        for j, tup in enumerate(amide_tuples):
            omega = _dihedral(*[coords[k] for k in tup])
            omega_vals[i, j] = omega
            cis_flags[i, j] = 1.0 if abs(omega) < CIS_OMEGA_CUTOFF else 0.0

    def bw(arr):
        arr = np.asarray(arr, dtype=float)
        mask = np.isfinite(arr)
        if not mask.any():
            return float("nan")
        ww = w[mask] / w[mask].sum()
        return float(np.dot(ww, arr[mask]))

    def bw_std(arr):
        """Boltzmann-weighted std = regime-1 configurational fluctuation amplitude."""
        arr = np.asarray(arr, dtype=float)
        mask = np.isfinite(arr)
        if mask.sum() < 2:
            return float("nan")
        ww = w[mask] / w[mask].sum()
        mean = np.dot(ww, arr[mask])
        var = np.dot(ww, (arr[mask] - mean) ** 2)
        return float(np.sqrt(max(var, 0.0)))

    def omega_circvar():
        """Mean Boltzmann-weighted circular variance of backbone omega dihedrals.
        0 = rigid (all conformers same cis/trans geometry), →1 = freely rotating."""
        if n_amide == 0:
            return float("nan")
        cvs = []
        for j in range(n_amide):
            ang = np.deg2rad(omega_vals[:, j])
            mask = np.isfinite(ang)
            if mask.sum() < 2:
                continue
            ww = w[mask] / w[mask].sum()
            R = np.hypot(np.dot(ww, np.cos(ang[mask])), np.dot(ww, np.sin(ang[mask])))
            cvs.append(1.0 - R)
        return float(np.mean(cvs)) if cvs else float("nan")

    psa = ens["psa_json"]
    cis_prob = (w[:, None] * cis_flags).sum(axis=0)  # per amide bond

    # dominant conformer = highest Boltzmann weight
    dom = int(np.argmax(w))

    out = {
        f"{prefix}_n_confs": ens["n"],
        f"{prefix}_bw_psa": round(bw(psa), 2),
        # ── intramolecular H-bonds (recomputed from geometry; bb+res == IMHB) ──
        f"{prefix}_bw_IMHB": round(bw(imhb_l), 3),       # total (was bw_hb)
        f"{prefix}_bw_IMHBD": round(bw(imhbd_l), 3),     # distinct donors engaged
        f"{prefix}_bw_IMHBA": round(bw(imhba_l), 3),     # distinct acceptors engaged
        f"{prefix}_bw_IMHB_bb": round(bw(imhb_bb_l), 3), # backbone (transannular)
        f"{prefix}_bw_IMHB_res": round(bw(imhb_res_l), 3),  # side-chain (residue)
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
        # ── v3 surface descriptors (Boltzmann-weighted means) ──
        f"{prefix}_bw_SA_HD": round(bw(hbd_l), 2),              # donor-H surface (Rzepiela 2022)
        f"{prefix}_bw_SA_HA": round(bw(hba_l), 2),              # acceptor-atom surface
        f"{prefix}_bw_hydrophobic_sasa": round(bw(hydro_l), 2),  # ASA_H
        f"{prefix}_bw_amphi_moment": round(bw(amphi_l), 3),      # integy/amphipathic moment
        # ── regime-1 ensemble fluctuation (same E-weighted ensemble, variance not mean) ──
        f"{prefix}_std_rg": round(bw_std(rg), 3),
        f"{prefix}_std_asphericity": round(bw_std(asph), 4),
        f"{prefix}_std_amphi_moment": round(bw_std(amphi_l), 3),
        f"{prefix}_omega_circvar": round(omega_circvar(), 4),
        f"{prefix}_n_eff": effective_nconf(w),                   # nConf20 analog
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
        ("delta_psa", "bw_psa"), ("delta_IMHB", "bw_IMHB"), ("delta_rg", "bw_rg"),
        ("delta_IMHB_bb", "bw_IMHB_bb"), ("delta_IMHB_res", "bw_IMHB_res"),
        ("delta_npr1", "bw_npr1"), ("delta_npr2", "bw_npr2"),
        ("delta_asphericity", "bw_asphericity"),
        # v3 surface deltas (cross-solvent change in donor/acceptor exposure, hydrophobicity, amphipathicity)
        ("delta_SA_HD", "bw_SA_HD"),
        ("delta_SA_HA", "bw_SA_HA"),
        ("delta_hydrophobic_sasa", "bw_hydrophobic_sasa"),
        ("delta_amphi_moment", "bw_amphi_moment"),
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

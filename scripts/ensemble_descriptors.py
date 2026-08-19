# env: chameleon-calc
"""
ensemble_descriptors.py
-----------------------
Compute ML-ready 3D conformational descriptors from completed CREST ensembles.
Featurizer for the DynamicEnsembleEncoder (see docs/chameleon_model_architecture.md);
output columns map to the dynamic_features schema (docs/data_schema.md).

Operates directly on two solvent ensembles — water/ plus one apolar leg (chcl3/ or
hexane/, see --apolar) — using the solvent label as the state assignment (no clustering,
no max/min-PSA proxy). Per-solvent descriptors are Boltzmann-weighted over the real CREST
ensemble; cross-solvent descriptors are water − apolar differences (delta_* column names
stay prefix-free whichever apolar phase is used).

WEIGHTS — --energies-csv is REQUIRED. Descriptors are Boltzmann-weighted by the SOLVATED
CPCM-X single-point populations from free_energy_calculator.py (stage 2). This is the
Role-2 -> Role-3 hand-off: geometry from CREST, populations from the solvated scoring.
There is NO ALPB/GFN2 fallback — those energies are the wrong footing for a solvated
Boltzmann ensemble, so a missing/broken populations file is an error, not a downgrade.

Reads BOTH ensemble layouts: current crest_engine.py (ensemble.sdf + ensemble.xyz +
metadata.json) and the legacy one (ensemble.sdf + ensemble.json).

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

Usage (--energies-csv is REQUIRED for every run):
  # water vs chloroform (default apolar leg), weighted by CPCM-X populations
  python ensemble_descriptors.py --run-dir results/conformers/HexPep --name HexPep \
      --energies-csv results/free_energy/fe_HexPep.csv

  # water vs hexane, weighted by the CPCM-X populations from the scoring run
  python ensemble_descriptors.py --run-dir results/conformers/HexPep --apolar hexane \
      --energies-csv results/free_energy/fe_HexPep.csv --name HexPep

  # several compounds into one table
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
    surface_descriptors_mol, effective_nconf, weighted_rmsf, kier_flexibility,
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
def _energies_from_xyz(xyz_path: Path) -> list[float]:
    """Per-conformer energies (Hartree) from a multi-frame CREST ensemble.xyz: the first
    float token on each frame's comment line, as CREST writes it. The current
    crest_engine.py emits ensemble.xyz (+ metadata.json) and no longer writes the legacy
    ensemble.json, so this is where energies come from now."""
    energies: list[float] = []
    lines = xyz_path.read_text().splitlines()
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            natoms = int(lines[i].split()[0])
        except (ValueError, IndexError):
            i += 1
            continue
        comment = lines[i + 1] if i + 1 < len(lines) else ""
        e = float("nan")
        for tok in comment.replace("=", " ").replace(":", " ").split():
            try:
                e = float(tok)
                break
            except ValueError:
                continue
        energies.append(e)
        i += 2 + natoms
    return energies


def _pops_from_energy_csv(csv_path: Path, solvent: str, n_conf: int):
    """Boltzmann populations for ONE solvent leg from free_energy_calculator.py's
    per-conformer CSV (CPCM-X / ALPB single-points) — the Role-2 -> Role-3 hand-off.

    The CSV's `conf` column holds the ORIGINAL conformer index, and with --ewin only the
    low-energy subset was scored. We therefore scatter pops back by index over the full
    ensemble; conformers outside the window keep weight 0, which is what their Boltzmann
    weight is anyway. Returns None if the CSV can't supply usable weights."""
    try:
        df = pd.read_csv(csv_path)
    except Exception:
        return None
    if "solvent" in df.columns:
        df = df[df["solvent"] == solvent]
    if df.empty or "pop" not in df.columns or "conf" not in df.columns:
        return None
    w = np.zeros(n_conf, dtype=float)
    idx = df["conf"].to_numpy(dtype=int)
    pops = df["pop"].to_numpy(dtype=float)
    ok = (idx >= 0) & (idx < n_conf) & np.isfinite(pops)
    if not ok.any():
        return None
    w[idx[ok]] = pops[ok]
    return w if w.sum() > 0 else None


def load_ensemble(solv_dir: Path, energies_csv: Path | None = None,
                  solvent_key: str | None = None):
    """Load one solvent leg: geometry + per-conformer energies + Boltzmann weights.

    Handles BOTH ensemble layouts:
      * current  (crest_engine.py): ensemble.sdf + ensemble.xyz + metadata.json
      * legacy:                     ensemble.sdf + ensemble.json (totalenergy/boltzmannweight)

    Weights come ONLY from the CPCM-X solvated single-point populations in `energies_csv`
    (free_energy_calculator.py `pop`, keyed by solvent). There is NO fallback to ALPB/GFN2
    ensemble energies -- those are the wrong footing for a solvated Boltzmann ensemble. If the
    populations can't be loaded for a leg, this raises (stage 2 not run, or broken path).
    """
    sdf_path = solv_dir / "ensemble.sdf"
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
    mols = [m for m in supplier if m is not None]
    n = len(mols)
    if n == 0:
        return {"mols": [], "weights": np.array([]), "energies": np.array([]),
                "psa_json": np.array([]), "hb_json": np.array([]), "smiles": None, "n": 0}

    energies = np.full(n, np.nan)
    weights = None
    psa_json = np.full(n, np.nan)
    hb_json = np.full(n, np.nan)
    smiles = None

    json_path = solv_dir / "ensemble.json"
    if json_path.exists():                                   # legacy layout
        with open(json_path) as f:
            data = json.load(f)
        confs = data.get("conformers", [])
        m = min(n, len(confs))
        energies[:m] = [confs[i].get("totalenergy", np.nan) for i in range(m)]
        w = np.array([confs[i].get("boltzmannweight", np.nan) for i in range(m)], dtype=float)
        if m == n and np.all(np.isfinite(w)) and w.sum() > 0:
            weights = w
        psa_json[:m] = [confs[i].get("psa", np.nan) for i in range(m)]
        hb_json[:m] = [confs[i].get("hbonds", np.nan) for i in range(m)]
        smiles = data.get("smiles")
    else:                                                    # current layout
        xyz_path = solv_dir / "ensemble.xyz"
        if xyz_path.exists():
            e = np.asarray(_energies_from_xyz(xyz_path), dtype=float)
            m = min(n, e.size)
            energies[:m] = e[:m]
        meta_path = solv_dir / "metadata.json"
        if meta_path.exists():
            try:
                smiles = json.loads(meta_path.read_text()).get("smiles")
            except Exception:
                pass

    # Boltzmann weights MUST come from the CPCM-X solvated single-point populations
    # (free_energy_calculator.py, stage 2). ALPB/GFN2 ensemble energies are the wrong
    # footing for a SOLVATED Boltzmann ensemble, so there is deliberately NO fallback:
    # if the CPCM-X populations can't be loaded, either stage 2 was not run for this leg
    # or the energies path is broken -- both are bugs to fix, not to mask with worse weights.
    if energies_csv is None or not solvent_key:
        raise ValueError(
            f"{solv_dir}: CPCM-X populations are required for weighting, but no energies CSV "
            f"/ solvent key was provided. Run stage 2 (free_energy_calculator.py) and pass "
            f"--energies-csv.")
    weights = _pops_from_energy_csv(Path(energies_csv), solvent_key, n)
    if weights is None:
        raise ValueError(
            f"{solv_dir}: could not load CPCM-X populations for solvent '{solvent_key}' from "
            f"{energies_csv} (n_conf={n}). Either that leg was not scored in stage 2, or the "
            f"energies path is broken (missing pop/conf/solvent columns, or a conf-index "
            f"mismatch with ensemble.sdf). NOT falling back to ALPB weights.")
    weights = np.nan_to_num(np.asarray(weights, dtype=float), nan=0.0)
    total = weights.sum()
    if total <= 0:
        raise ValueError(f"{solv_dir}: CPCM-X populations for '{solvent_key}' sum to 0.")
    weights = weights / total

    return {"mols": mols, "weights": weights, "energies": energies,
            "psa_json": psa_json, "hb_json": hb_json, "smiles": smiles, "n": n}


# ── per-solvent descriptors ───────────────────────────────────────────────────
def solvent_descriptors(solv_dir: Path, prefix: str, energies_csv: Path | None = None,
                        solvent_key: str | None = None) -> dict:
    ens = load_ensemble(solv_dir, energies_csv, solvent_key or prefix)
    mols, w = ens["mols"], ens["weights"]
    if len(mols) == 0:
        return {}

    amide_tuples = amide_ring_bonds(mols[0])
    n_amide = len(amide_tuples)
    backbone_atoms = backbone_hbond_atoms(mols[0])  # ring + carbonyl O's = backbone for IMHB

    rg, npr1, npr2, asph, sphe, tsasa = [], [], [], [], [], []
    hbd_l, hba_l, hydro_l, amphi_l, psa_l = [], [], [], [], []   # v3 surface descriptors
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
        sd = surface_descriptors_mol(m, cid)    # literature 3D-PSA, SA_HD/HA, hydrophobic, amphi
        hbd_l.append(sd["hbd_sasa"]); hba_l.append(sd["hba_sasa"])
        hydro_l.append(sd["hydrophobic_sasa"])
        amphi_l.append(sd["amphi_moment"])
        psa_l.append(sd["psa"]); tsasa.append(sd["total_sasa"])   # consistent Bondi-radii SASA
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

    cis_prob = (w[:, None] * cis_flags).sum(axis=0)  # per amide bond

    # dominant conformer = highest Boltzmann weight
    dom = int(np.argmax(w))

    out = {
        f"{prefix}_n_confs": ens["n"],
        f"{prefix}_bw_psa": round(bw(psa_l), 2),    # literature 3D-PSA (N/O + polar H + oxidized S)
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
        f"{prefix}_psa_std": round(bw_std(psa_l), 2),     # Boltzmann-weighted (was unweighted)
        f"{prefix}_psa_spread": round(float(np.nanmax(psa_l) - np.nanmin(psa_l)), 2),
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
        f"{prefix}_rmsf": weighted_rmsf(mols, w),               # threshold-free flexibility (preferred)
        # NOTE: p_dominant + n_eff are discretization-sensitive (conformer- not fold-level);
        # superseded by rmsf for the flexibility axis. Kept for benchmark comparison only.
        f"{prefix}_n_eff": effective_nconf(w),                   # nConf20 analog (raw; see note)
        f"{prefix}_kier_phi": kier_flexibility(mols[0]),         # Kier 2D flexibility (Begnini 2021)
    }
    for j in range(n_amide):
        out[f"{prefix}_cis_prob_{j}"] = round(float(cis_prob[j]), 3)
    # store internal for cross-solvent
    out[f"_{prefix}_cis_prob_vec"] = cis_prob.tolist()
    out[f"_{prefix}_dom_energy"] = float(ens["energies"][dom]) if np.isfinite(ens["energies"][dom]) else None
    return out


def cross_solvent(water: dict, apolar: dict, ap: str = "chcl3") -> dict:
    """Cross-solvent deltas: water minus the APOLAR leg (`ap` = its folder/column prefix,
    e.g. chcl3 or hexane). Delta column names stay prefix-free so downstream consumers
    (reports, notebooks, ml/) keep reading the same schema regardless of which apolar
    phase was used."""
    out = {}
    pairs = [
        ("delta_psa", "bw_psa"), ("delta_IMHB", "bw_IMHB"), ("delta_rg", "bw_rg"),
        ("delta_IMHB_bb", "bw_IMHB_bb"), ("delta_IMHB_res", "bw_IMHB_res"),
        ("delta_npr1", "bw_npr1"), ("delta_npr2", "bw_npr2"),
        ("delta_asphericity", "bw_asphericity"), ("delta_rmsf", "rmsf"),
        # v3 surface deltas (cross-solvent change in donor/acceptor exposure, hydrophobicity, amphipathicity)
        ("delta_SA_HD", "bw_SA_HD"),
        ("delta_SA_HA", "bw_SA_HA"),
        ("delta_hydrophobic_sasa", "bw_hydrophobic_sasa"),
        ("delta_amphi_moment", "bw_amphi_moment"),
    ]
    for out_key, feat in pairs:
        wv, mv = water.get(f"water_{feat}"), apolar.get(f"{ap}_{feat}")
        if wv is not None and mv is not None:
            out[out_key] = round(wv - mv, 3)

    # normalized delta PSA (Yu 2026): ΔPSA / mean total SASA
    if "delta_psa" in out:
        sasa = np.nanmean([water.get("water_sasa_total", np.nan),
                           apolar.get(f"{ap}_sasa_total", np.nan)])
        if np.isfinite(sasa) and sasa > 0:
            out["norm_delta_psa"] = round(out["delta_psa"] / sasa, 5)

    # ddG between dominant conformers across solvents (Hartree → kcal)
    ew, em = water.get("_water_dom_energy"), apolar.get(f"_{ap}_dom_energy")
    if ew is not None and em is not None:
        out["ddG_dom_kcal"] = round((ew - em) * 627.509, 2)

    # cis-amide switch: which bond changes most between solvents
    cw = water.get("_water_cis_prob_vec")
    cm = apolar.get(f"_{ap}_cis_prob_vec")
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


def _has_ensemble(d: Path | None) -> bool:
    """A usable leg = geometry (ensemble.sdf) present. Energies come from ensemble.xyz
    (current layout) or ensemble.json (legacy); load_ensemble handles either."""
    return bool(d) and (d / "ensemble.sdf").exists()


def process_compound(name: str, water_dir: Path | None, apolar_dir: Path | None,
                     ap: str = "chcl3", energies_csv: Path | None = None) -> dict:
    """One compound = the water leg + one apolar leg (`ap`: chcl3 or hexane).
    If `energies_csv` is given, both legs are Boltzmann-weighted by the solvated
    single-point populations from free_energy_calculator.py instead of raw GFN2."""
    row = {"compound": name, "apolar_solvent": ap}
    wd = {}
    md = {}
    if _has_ensemble(water_dir):
        wd = solvent_descriptors(water_dir, "water", energies_csv, "water")
        row.update({k: v for k, v in wd.items() if not k.startswith("_")})
    if _has_ensemble(apolar_dir):
        md = solvent_descriptors(apolar_dir, ap, energies_csv, ap)
        row.update({k: v for k, v in md.items() if not k.startswith("_")})
    if wd and md:
        row.update(cross_solvent(wd, md, ap))
        row["has_both_solvents"] = 1
    else:
        row["has_both_solvents"] = 0
    return row


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", action="append", default=[],
                   help="compound dir containing water/ and the apolar leg's subdir "
                        "(see --apolar). Works on results/conformers/<name>/ or "
                        "results/runs/run_*/ from crest_v3.2.py.")
    p.add_argument("--water-dir", action="append", default=[],
                   help="Directory holding a water-only ensemble (ensemble.sdf [+ .xyz/.json])")
    p.add_argument("--apolar", default="chcl3",
                   help="apolar leg's folder name = its column prefix (default: chcl3; "
                        "use hexane for the hexane transfer phase)")
    p.add_argument("--energies-csv", type=Path, required=True,
                   help="REQUIRED. free_energy_calculator.py per-conformer CSV. Descriptors are "
                        "Boltzmann-weighted by the SOLVATED CPCM-X populations from that run. "
                        "There is no ALPB/GFN2 fallback: if a leg's populations are missing the "
                        "run errors (stage 2 not run, or a broken energies path). Conformers "
                        "trimmed by --ewin correctly get weight 0.")
    p.add_argument("--name", action="append", default=[], help="Compound name (one per --run-dir/--water-dir)")
    p.add_argument("-o", "--out", default="results/ensemble_descriptors.csv", type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    ap = args.apolar
    jobs = []
    for i, rd in enumerate(args.run_dir):
        nm = args.name[i] if i < len(args.name) else Path(rd).name
        jobs.append((nm, Path(rd) / "water", Path(rd) / ap))
    offset = len(args.run_dir)
    for i, wdir in enumerate(args.water_dir):
        nm = args.name[offset + i] if offset + i < len(args.name) else Path(wdir).name
        jobs.append((nm, Path(wdir), None))

    print(f"weighting by CPCM-X solvated populations from: {args.energies_csv}")

    rows = []
    for nm, wdir, mdir in jobs:
        print(f"[{nm}] water={wdir} {ap}={mdir}")
        rows.append(process_compound(nm, wdir, mdir, ap, args.energies_csv))

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

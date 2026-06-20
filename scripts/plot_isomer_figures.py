# env: chameleon-calc
"""
plot_isomer_figures.py
----------------------
Publication figures for the non-diazirine DOPC R/S isomer pairs (plain xylene linker):
  pair A: 3-12-8-12  R vs S
  pair B: 3-12-10-12 R vs S

Per pair, from the intact CREST ensembles (results/conformers/), generates 3 figures
(SVG = Illustrator-editable + PNG) into results/figures/isomers/:

  reldiff_<pair>   : %|R-S| per descriptor. 2D/lipophilicity descriptors pinned at ~0
                     (blind to the stereocenter); 3D ensemble descriptors stick out.
  hbonds_<pair>    : native-unit box plots of the intramolecular H-bond family
                     (IMHB total, backbone vs residue, donors/acceptors) x R/S x water/mem.
  overlap3d_<pair> : robust-scaled (median/IQR) box plots of the continuous 3D descriptors,
                     water phase, R vs S — shows multi-axis divergence on one comparable axis.

Recomputes descriptors per conformer (phys_descriptors_v3) so the box plots show the real
ensemble distributions, not just Boltzmann means.

Usage:  python scripts/plot_isomer_figures.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Descriptors3D, rdMolDescriptors

sys.path.insert(0, str(Path(__file__).parent))
from phys_descriptors_v3 import (surface_descriptors_mol, imhb_descriptors_mol,
                                 backbone_hbond_atoms)

RDLogger.DisableLog("rdApp.*")

BASE = Path("results/conformers")
OUT = Path("results/figures/isomers")
OUT.mkdir(parents=True, exist_ok=True)
PAIRS = {
    "3-12-8-12":  {"R": BASE / "DOPC 3-12-8-12 R", "S": BASE / "DOPC 3-12-8-12 S"},
    "3-12-10-12": {"R": BASE / "3-12-10-12 R",     "S": BASE / "3-12-10-12 S"},
}
SOLVENTS = ["water", "mem"]
COL = {"R": "#d1495b", "S": "#30638e"}
HB_FAMILY = ["IMHB", "IMHB_bb", "IMHB_res", "IMHBD", "IMHBA"]
CONT_3D = ["rg", "asph", "sphe", "psa", "SA_HD", "SA_HA", "hydrophobic", "amphi"]
LBL = {"rg": "Rg", "asph": "asphericity", "sphe": "spherocity", "psa": "3D-PSA",
       "SA_HD": "SA_HD", "SA_HA": "SA_HA", "hydrophobic": "hydrophobic SASA",
       "amphi": "amphi. moment", "IMHB": "IMHB", "IMHB_bb": "IMHB backbone",
       "IMHB_res": "IMHB side-chain", "IMHBD": "IMHB donors", "IMHBA": "IMHB acceptors"}


# ── per-conformer descriptor table ────────────────────────────────────────────
def per_conf_df(run_dir: Path, isomer: str, pair: str) -> pd.DataFrame:
    frames = []
    for solv in SOLVENTS:
        jp, sp = run_dir / solv / "ensemble.json", run_dir / solv / "ensemble.sdf"
        if not jp.exists() or not sp.exists():
            continue
        confs = json.load(open(jp)).get("conformers", [])
        w = np.array([c.get("boltzmannweight", np.nan) for c in confs], float)
        psa = [c.get("psa", np.nan) for c in confs]
        mols = [m for m in Chem.SDMolSupplier(str(sp), removeHs=False, sanitize=True) if m]
        n = min(len(mols), len(confs))
        if n == 0:
            continue
        bb = backbone_hbond_atoms(mols[0])
        rows = []
        for i in range(n):
            m = mols[i]; cid = m.GetConformer().GetId()
            try:
                sd = surface_descriptors_mol(m, cid)
                ih = imhb_descriptors_mol(m, cid, bb)
                rows.append(dict(
                    pair=pair, isomer=isomer, solvent=solv, w=w[i], psa=sd["psa"],
                    rg=Descriptors3D.RadiusOfGyration(m, confId=cid),
                    asph=Descriptors3D.Asphericity(m, confId=cid),
                    sphe=Descriptors3D.SpherocityIndex(m, confId=cid),
                    SA_HD=sd["hbd_sasa"], SA_HA=sd["hba_sasa"],
                    hydrophobic=sd["hydrophobic_sasa"], amphi=sd["amphi_moment"],
                    IMHB=ih["imhb"], IMHB_bb=ih["imhb_bb"], IMHB_res=ih["imhb_res"],
                    IMHBD=ih["imhbd"], IMHBA=ih["imhba"]))
            except Exception:
                pass
        frames.append(pd.DataFrame(rows))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def bw_mean(df: pd.DataFrame, col: str) -> float:
    d = df.dropna(subset=[col, "w"])
    if d.empty or d["w"].sum() == 0:
        return float("nan")
    return float(np.average(d[col], weights=d["w"]))


def two_d(smiles: str) -> dict:
    m = Chem.MolFromSmiles(smiles)
    return {"MolWt": Descriptors.MolWt(m), "TPSA": Descriptors.TPSA(m),
            "cLogP": Crippen.MolLogP(m), "HBD": rdMolDescriptors.CalcNumHBD(m),
            "HBA": rdMolDescriptors.CalcNumHBA(m),
            "RotB": rdMolDescriptors.CalcNumRotatableBonds(m),
            "FrCSP3": rdMolDescriptors.CalcFractionCSP3(m)}


def _reldiff(a, b):
    den = (abs(a) + abs(b)) / 2
    return float("nan") if den == 0 else 100 * abs(a - b) / den


def _smiles(run_dir: Path) -> str:
    for solv in SOLVENTS:
        jp = run_dir / solv / "ensemble.json"
        if jp.exists():
            return json.load(open(jp)).get("smiles")
    return None


# ── figures ───────────────────────────────────────────────────────────────────
def fig_reldiff(pair, dfR, dfS, smiR, smiS):
    two_R, two_S = two_d(smiR), two_d(smiS)
    rows = [(f"2D: {k}", _reldiff(two_R[k], two_S[k]), "#9aa0a6") for k in two_R]
    fam_col = {"IMHB": "#2a9d8f", "IMHB_bb": "#2a9d8f", "IMHB_res": "#2a9d8f",
               "IMHBD": "#2a9d8f", "IMHBA": "#2a9d8f",
               "psa": "#e76f51", "SA_HD": "#e76f51", "SA_HA": "#e76f51",
               "hydrophobic": "#e76f51", "amphi": "#e76f51",
               "rg": "#e9c46a", "asph": "#e9c46a", "sphe": "#e9c46a"}
    wR, wS = dfR[dfR.solvent == "water"], dfS[dfS.solvent == "water"]
    for col in ["IMHB", "IMHB_bb", "IMHB_res", "IMHBD", "IMHBA",
                "psa", "SA_HD", "SA_HA", "hydrophobic", "amphi", "rg", "asph", "sphe"]:
        rows.append((f"3D: {LBL[col]}", _reldiff(bw_mean(wR, col), bw_mean(wS, col)),
                     fam_col[col]))
    rows = [r for r in rows if not np.isnan(r[1])]
    rows.sort(key=lambda r: r[1])
    labels, vals, colors = zip(*rows)
    fig, ax = plt.subplots(figsize=(7.5, 0.34 * len(rows) + 1))
    ax.barh(range(len(rows)), vals, color=colors)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("relative |R − S| difference (%, water phase)")
    ax.set_title(f"DOPC {pair}: 2D blind, 3D resolves the stereocenter")
    ax.axvline(0, color="k", lw=0.6)
    fig.tight_layout()
    _save(fig, f"reldiff_{pair}")


def _boxes(ax, groups, descr):
    """groups: list of (label, color, array). Draw a box per group."""
    data = [g[2][np.isfinite(g[2])] for g in groups]
    bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=False)
    for patch, g in zip(bp["boxes"], groups):
        patch.set_facecolor(g[1]); patch.set_alpha(0.65)
    for med in bp["medians"]:
        med.set_color("k")
    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels([g[0] for g in groups], fontsize=7)
    ax.set_title(descr, fontsize=9)


def fig_hbonds(pair, dfR, dfS):
    fig, axes = plt.subplots(1, len(HB_FAMILY), figsize=(2.5 * len(HB_FAMILY), 3.2))
    for ax, col in zip(axes, HB_FAMILY):
        groups = []
        for iso, df in (("R", dfR), ("S", dfS)):
            for solv in SOLVENTS:
                sub = df[df.solvent == solv]
                arr = sub[col].to_numpy(float) if not sub.empty else np.array([])
                groups.append((f"{iso}\n{('aq' if solv=='water' else 'mem')}",
                               COL[iso], arr))
        _boxes(ax, groups, LBL[col])
        ax.set_ylabel("count" if col == HB_FAMILY[0] else "")
    fig.suptitle(f"DOPC {pair}: intramolecular H-bonds (R/S × water/membrane)")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save(fig, f"hbonds_{pair}")


def fig_overlap3d(pair, dfR, dfS):
    """Robust-scaled (pooled median/IQR per descriptor) box plots, water phase, R vs S."""
    wR, wS = dfR[dfR.solvent == "water"], dfS[dfS.solvent == "water"]
    fig, ax = plt.subplots(figsize=(1.05 * len(CONT_3D) + 1.5, 4))
    pos, ticks, ticklabels = 1, [], []
    for col in CONT_3D:
        pooled = np.concatenate([wR[col].to_numpy(float), wS[col].to_numpy(float)])
        pooled = pooled[np.isfinite(pooled)]
        med = np.median(pooled)
        iqr = np.subtract(*np.percentile(pooled, [75, 25]))
        scale = iqr if iqr > 1e-9 else (np.std(pooled) or 1.0)
        rs = lambda a: (a[np.isfinite(a)] - med) / scale
        for j, (iso, df) in enumerate((("R", wR), ("S", wS))):
            bp = ax.boxplot(rs(df[col].to_numpy(float)), positions=[pos + j * 0.42],
                            widths=0.36, patch_artist=True, showfliers=False)
            bp["boxes"][0].set_facecolor(COL[iso]); bp["boxes"][0].set_alpha(0.65)
            bp["medians"][0].set_color("k")
        ticks.append(pos + 0.21); ticklabels.append(LBL[col]); pos += 1.4
    ax.axhline(0, color="grey", lw=0.6, ls="--")
    ax.set_xticks(ticks); ax.set_xticklabels(ticklabels, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("robust-scaled (median-centred, IQR units)")
    ax.set_title(f"DOPC {pair}: 3D descriptor distributions, water — R vs S")
    handles = [plt.Line2D([0], [0], color=COL[i], lw=6, alpha=0.65) for i in ("R", "S")]
    ax.legend(handles, ["R", "S"], loc="upper right", fontsize=8)
    fig.tight_layout()
    _save(fig, f"overlap3d_{pair}")


def _save(fig, name):
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {name}.svg / .png")


def main():
    for pair, dirs in PAIRS.items():
        print(f"[{pair}] loading ensembles…")
        dfR = per_conf_df(dirs["R"], "R", pair)
        dfS = per_conf_df(dirs["S"], "S", pair)
        if dfR.empty or dfS.empty:
            print(f"  !! missing data for {pair}, skipping")
            continue
        smiR, smiS = _smiles(dirs["R"]), _smiles(dirs["S"])
        print(f"  R: {len(dfR)} confs, S: {len(dfS)} confs")
        fig_reldiff(pair, dfR, dfS, smiR, smiS)
        fig_hbonds(pair, dfR, dfS)
        fig_overlap3d(pair, dfR, dfS)
    print(f"\nFigures → {OUT}")


if __name__ == "__main__":
    main()

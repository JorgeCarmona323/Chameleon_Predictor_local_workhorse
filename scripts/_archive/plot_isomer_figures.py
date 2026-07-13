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
    "3-12-8-12":  {"R": BASE / "DOPC 3-12-8-12/3-12-8-12 Xylene Linker/DOPC 3-12-8-12 R",
                   "S": BASE / "DOPC 3-12-8-12/3-12-8-12 Xylene Linker/DOPC 3-12-8-12 S"},
    "3-12-10-12": {"R": BASE / "WhC3/3-12-10-12 Xylene Linker/3-12-10-12 R",
                   "S": BASE / "WhC3/3-12-10-12 Xylene Linker/3-12-10-12 S"},
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
                    npr1=Descriptors3D.NPR1(m, confId=cid),   # I1/I3, disc/sphere axis
                    npr2=Descriptors3D.NPR2(m, confId=cid),   # I2/I3, rod axis
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
    _save(fig, pair, "fig1_reldiff")


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
                groups.append((f"{iso}\n{('aq' if solv=='water' else 'chcl3')}",
                               COL[iso], arr))
        _boxes(ax, groups, LBL[col])
        ax.set_ylabel("count" if col == HB_FAMILY[0] else "")
    fig.suptitle(f"DOPC {pair}: intramolecular H-bonds (R/S × water/chloroform)")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save(fig, pair, "fig2_hbonds")


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
    _save(fig, pair, "fig3_overlap3d")


# ── Figure 2 (NEW): validated-core descriptors, Begnini Fig-4 style ────────────
# Replaces the old hbonds + overlap3d figures. Rows = validated descriptors, columns =
# solvents (water | chloroform); within each panel R vs S. Per-conformer points sized by
# Boltzmann weight, minimum-energy conformer (MEC) ringed, weighted mean printed, and a
# Mann-Whitney R-vs-S p-value (Begnini 2021 used a Wilcoxon/Mann-Whitney test).
KEY_DESCRIPTORS = [
    ("psa", "SA 3D-PSA (Å²)"),
    ("rg", "radius of gyration (Å)"),
]   # backbone IMHB moved to its own population figure (discrete counts read poorly as box plots)


def _sig_stars(p):
    if p is None:
        return ""
    return ("****" if p <= 1e-4 else "***" if p <= 1e-3 else
            "**" if p <= 1e-2 else "*" if p <= 0.05 else "ns")


def _color_box(bp, iso):
    """Begnini-style: open box (transparent face, so dots show through) with a thin colored
    edge, colored median, colored whiskers/caps. bp is a single-box boxplot dict."""
    c = COL[iso]
    for b in bp["boxes"]:
        b.set_facecolor("none"); b.set_edgecolor(c); b.set_linewidth(1.6); b.set_zorder(3)
    for m in bp["medians"]:
        m.set_color(c); m.set_linewidth(1.8); m.set_zorder(3.1)
    for w in bp["whiskers"]:
        w.set_color(c); w.set_linewidth(1.2)
    for cp in bp["caps"]:
        cp.set_color(c); cp.set_linewidth(1.2)


def _representative(vals, wts, n=50):
    """Reduce a CREST ensemble (hundreds of conformers) to the n most Boltzmann-populated
    conformers, plotted as uniform circles (Begnini Fig 5: conformers as equal-size dots,
    not sized by weight). Returns (vals, wts) sorted by descending weight, truncated to n."""
    order = np.argsort(wts)[::-1][:n]
    return vals[order], wts[order]


def _key_panel(ax, dfR, dfS, col, solv):
    rng = np.random.default_rng(0)
    series, shown = {}, {}
    for iso, df in (("R", dfR), ("S", dfS)):
        sub = df[df.solvent == solv].dropna(subset=[col, "w"])
        vals, wts = sub[col].to_numpy(float), sub["w"].to_numpy(float)
        series[iso] = (vals, wts)
        shown[iso] = _representative(vals, wts) if len(vals) else (vals, wts)
    # box over the representative (populated) conformers, one call per isomer for its color
    for pos, iso in ((1, "R"), (2, "S")):
        if len(shown[iso][0]):
            _color_box(ax.boxplot([shown[iso][0]], positions=[pos], widths=0.42,
                                  patch_artist=True, showfliers=False), iso)
    means = {}
    for pos, iso in ((1, "R"), (2, "S")):
        vals, wts = shown[iso]
        if len(vals) == 0:
            continue
        x = rng.normal(pos, 0.08, size=len(vals))       # jitter cloud wider than the box
        ax.scatter(x, vals, s=26, color=COL[iso], alpha=0.6, edgecolor="white",
                   linewidth=0.3, zorder=2)             # uniform circles (Begnini Fig 5)
        mec = int(np.argmax(wts))                       # lowest-energy = highest weight
        ax.scatter([x[mec]], [vals[mec]], s=26, color=COL[iso], alpha=0.95, zorder=3)
        ax.scatter([x[mec]], [vals[mec]], s=95, facecolor="none",
                   edgecolor="#1f3b73", linewidth=1.5, zorder=4)   # MEC marked by a ring
        # weighted mean over the FULL ensemble (not just shown), the reported quantity
        fv, fw = series[iso]
        means[iso] = float(np.average(fv, weights=fw))
    p = None
    try:
        from scipy.stats import mannwhitneyu
        if len(series["R"][0]) and len(series["S"][0]):
            p = mannwhitneyu(series["R"][0], series["S"][0], alternative="two-sided").pvalue
    except Exception:
        p = None
    # significance bracket with end caps, above the data (Begnini style)
    y0, y1 = ax.get_ylim()
    span = y1 - y0
    yb = y1 + 0.05 * span
    ax.plot([1, 1, 2, 2], [yb - 0.02 * span, yb, yb, yb - 0.02 * span], color="k", lw=0.9)
    ax.text(1.5, yb + 0.005 * span, _sig_stars(p), ha="center", va="bottom", fontsize=10)
    ax.set_ylim(y0, yb + 0.12 * span)
    ax.set_xlim(0.5, 2.5)
    ax.set_xticks([1, 2]); ax.set_xticklabels(["R", "S"], fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8.5)
    return means


def _mean_row(ax, means):
    """Weighted-mean strip beneath a panel, Begnini's 'Mean' band."""
    from matplotlib.transforms import blended_transform_factory
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    ax.text(0.5, -0.17, "mean", transform=ax.transAxes, ha="left", va="top",
            fontsize=7.5, color="0.35", style="italic")
    for pos, iso in ((1, "R"), (2, "S")):
        if iso in means:
            ax.text(pos, -0.17, f"{means[iso]:.1f}", transform=tr, ha="center", va="top",
                    fontsize=8.5, color=COL[iso], fontweight="bold")


def fig_key_descriptors(pair, dfR, dfS):
    solv_cols = [("water", "water"), ("mem", "chloroform")]
    nrows, ncols = len(KEY_DESCRIPTORS), len(solv_cols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.8 * ncols, 4.2 * nrows), squeeze=False)
    for i, (col, label) in enumerate(KEY_DESCRIPTORS):
        for j, (solv_key, solv_lab) in enumerate(solv_cols):
            ax = axes[i][j]
            means = _key_panel(ax, dfR, dfS, col, solv_key)
            _mean_row(ax, means)
            if i == 0:
                ax.set_title(solv_lab, fontsize=12, fontweight="bold", pad=8)
            if j == 0:
                ax.set_ylabel(label, fontsize=10.5)
    fig.suptitle(f"DOPC {pair}: validated 3D descriptors, R vs S", fontsize=13, fontweight="bold")
    fig.text(0.5, 0.005, "open box = IQR (25-75th) with median; whiskers = 1.5×IQR; dots = 50 most "
             "Boltzmann-populated conformers (uniform); ringed = minimum-energy conformer (MEC); "
             "mean = Boltzmann-weighted over full ensemble; Wilcoxon rank-sum (Mann-Whitney U) R vs S "
             "(* ≤0.05, ** ≤0.01, *** ≤0.001, **** ≤0.0001)",
             ha="center", fontsize=7.5, color="0.30")
    fig.tight_layout(rect=[0, 0.035, 1, 0.95], h_pad=2.4, w_pad=2.0)
    _save(fig, pair, "fig2_key3d")


def fig_imhb(pair, dfR, dfS):
    """Backbone (transannular) IMHB as a Boltzmann-weighted population bar chart — a discrete
    count reads far better this way than as a box plot. One panel per solvent; for each isomer,
    bars give the fraction of the ensemble (by Boltzmann weight) holding 0,1,2,... transannular
    H-bonds, with the weighted-mean count marked."""
    solv_cols = [("water", "water"), ("mem", "chloroform")]
    allc = pd.concat([dfR, dfS]).dropna(subset=["IMHB_bb"])
    cmax = int(round(allc["IMHB_bb"].max())) if len(allc) else 4
    counts = list(range(0, cmax + 1))
    width = 0.38
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.9), squeeze=False, sharey=True)
    for j, (sk, sl) in enumerate(solv_cols):
        ax = axes[0][j]
        for off, iso, df in ((-width / 2, "R", dfR), (width / 2, "S", dfS)):
            sub = df[df.solvent == sk].dropna(subset=["IMHB_bb", "w"])
            if sub.empty:
                continue
            w = sub["w"].to_numpy(float); c = np.round(sub["IMHB_bb"].to_numpy(float))
            tot = w.sum()
            frac = [100 * w[c == k].sum() / tot if tot else 0 for k in counts]
            ax.bar(np.array(counts) + off, frac, width, color=COL[iso], alpha=0.85,
                   edgecolor="white", linewidth=0.6, label=iso, zorder=2)
            wm = float(np.average(c, weights=w)) if tot else np.nan
            ax.axvline(wm, color=COL[iso], ls="--", lw=1.4, alpha=0.8, zorder=3)
            ax.text(wm, 102, f"{wm:.1f}", color=COL[iso], ha="center", va="bottom",
                    fontsize=8.5, fontweight="bold")
        ax.set_title(sl, fontsize=12, fontweight="bold")
        ax.set_xticks(counts); ax.set_ylim(0, 108)
        ax.set_xlabel("backbone (transannular) IMHB count", fontsize=10)
        if j == 0:
            ax.set_ylabel("Boltzmann-weighted population (%)", fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=9)
    axes[0][0].legend(frameon=False, fontsize=10, loc="upper left")
    fig.suptitle(f"DOPC {pair}: backbone IMHB population, R vs S", fontsize=13, fontweight="bold")
    fig.text(0.5, 0.005, "bars = fraction of the ensemble (by Boltzmann weight) at each transannular "
             "backbone H-bond count; dashed line = weighted-mean count",
             ha="center", fontsize=7.5, color="0.30")
    fig.tight_layout(rect=[0, 0.04, 1, 0.93])
    _save(fig, pair, "fig3_imhb")


def _pmi_panel(ax, dfR, dfS, solv):
    """Normalized PMI triangle (Begnini SI Fig S4): NPR1 (x) vs NPR2 (y), vertices
    Rod (0,1) / Sphere (1,1) / Disc (0.5,0.5); conformers scattered, colored rug margins."""
    # triangle frame
    ax.plot([0, 1, 0.5, 0], [1, 1, 0.5, 1], color="0.25", lw=2.0, zorder=1)
    ax.text(0.02, 1.01, "Rod", ha="left", va="bottom", fontsize=15, color="0.25")
    ax.text(0.98, 1.01, "Sphere", ha="right", va="bottom", fontsize=15, color="0.25")
    ax.text(0.50, 0.515, "Disc", ha="center", va="bottom", fontsize=15, color="0.25")
    for pos, iso, df in ((1, "R", dfR), (2, "S", dfS)):
        sub = df[df.solvent == solv].dropna(subset=["npr1", "npr2", "w"])
        if sub.empty:
            continue
        v = sub[["npr1", "npr2"]].to_numpy(float)       # full set, for the rug margins
        wts = sub["w"].to_numpy(float)
        vr, wr = _representative(v, wts)                 # 50 most populated, for the dots
        ax.scatter(vr[:, 0], vr[:, 1], s=75, color=COL[iso],
                   alpha=0.6, edgecolor="white", linewidth=0.5, zorder=3, label=iso)
        # marginal rug: NPR1 along bottom, NPR2 along left, colored by isomer
        yb = 0.452 if iso == "R" else 0.438
        ax.plot(v[:, 0], np.full(len(v), yb), "|", color=COL[iso], ms=9, alpha=0.5, zorder=2)
        xb = 0.008 if iso == "R" else 0.022
        ax.plot(np.full(len(v), xb), v[:, 1], "_", color=COL[iso], ms=9, alpha=0.5, zorder=2)
        mec = int(np.argmax(wts))
        ax.scatter([v[mec, 0]], [v[mec, 1]], s=120, facecolor="none", edgecolor="#1f3b73",
                   linewidth=2.2, zorder=4)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(0.42, 1.06)
    ax.set_xlabel("NPR1 (I₁/I₃)", fontsize=17)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0]); ax.set_yticks([0.5, 0.75, 1.0])
    ax.tick_params(labelsize=14)
    ax.grid(True, color="0.92", lw=1.0, zorder=0)
    ax.set_axisbelow(True)


def fig_pmi(pair, dfR, dfS):
    solv_cols = [("water", "water"), ("mem", "chloroform")]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 6.6), squeeze=False)
    for j, (solv_key, solv_lab) in enumerate(solv_cols):
        ax = axes[0][j]
        _pmi_panel(ax, dfR, dfS, solv_key)
        ax.set_title(solv_lab, fontsize=20, fontweight="bold")
        if j == 0:
            ax.set_ylabel("NPR2 (I₂/I₃)", fontsize=17)
    h = [plt.Line2D([], [], marker="o", ls="", color=COL[i], label=i) for i in ("R", "S")]
    axes[0][0].legend(handles=h, loc="lower left", frameon=False, fontsize=15,
                      handletextpad=0.3, bbox_to_anchor=(0.0, 0.02))
    fig.suptitle(f"DOPC {pair}: conformational shape space (PMI), R vs S",
                 fontsize=20, fontweight="bold")
    fig.text(0.5, 0.005, "normalized principal moments of inertia; dots = 50 most Boltzmann-populated "
             "conformers (uniform); rug = full ensemble; ringed = minimum-energy conformer",
             ha="center", fontsize=12, color="0.30")
    fig.tight_layout(rect=[0, 0.04, 1, 0.92])
    _save(fig, pair, "fig4_pmi")


def _save(fig, pair, role):
    d = OUT / pair
    d.mkdir(parents=True, exist_ok=True)
    fig.savefig(d / f"{role}.svg", bbox_inches="tight")            # vector, resolution-independent
    fig.savefig(d / f"{role}.png", dpi=600, bbox_inches="tight")   # 600 dpi for publication raster
    plt.close(fig)
    print(f"  saved {pair}/{role}.svg / .png")


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
        fig_key_descriptors(pair, dfR, dfS)   # fig2: PSA + Rg box plots
        fig_imhb(pair, dfR, dfS)              # fig3: backbone IMHB population bars
        fig_pmi(pair, dfR, dfS)               # fig4: PMI shape-space triangle (Begnini SI S4)
    print(f"\nFigures → {OUT}")


if __name__ == "__main__":
    main()

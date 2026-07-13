# env: chameleon-calc
"""
plot_isomer_comparison.py
-------------------------
Two figures for the DOPC 3-12-8-12 R/S isomer story:

  1. box_hbonds.* : box plots of intramolecular H-bonds per conformer
     (R/S x water/mem) -- shows the ensemble-level structural difference directly.
  2. rel_diff_2d_vs_3d.* : relative %|R-S| difference per descriptor; 2D/lipophilicity
     descriptors sit at zero (blind to stereocenter), 3D ensemble descriptors stick out.

Inputs: ensemble.json (per-conformer) + results/2d_descriptors.csv + the 3D
ensemble descriptors computed inline. Output -> results/figures/isomers/
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("results/conformers")
OUT = Path("results/figures/isomers")
OUT.mkdir(parents=True, exist_ok=True)
CASES = {
    "R_water": BASE / "DOPC 3-12-8-12 R" / "water",
    "R_mem":   BASE / "DOPC 3-12-8-12 R" / "mem",
    "S_water": BASE / "DOPC 3-12-8-12 S" / "water",
    "S_mem":   BASE / "DOPC 3-12-8-12 S" / "mem",
}
C_R, C_S = "#d1495b", "#30638e"


def per_conf(case_dir):
    with open(case_dir / "ensemble.json") as f:
        d = json.load(f)
    hb = np.array([c["hbonds"] for c in d["conformers"]], float)
    psa = np.array([c["psa"] for c in d["conformers"]], float)
    w = np.array([c["boltzmannweight"] for c in d["conformers"]], float)
    return hb, psa, w / w.sum()


data = {k: per_conf(v) for k, v in CASES.items()}

# ---------- Figure 1: box plots of H-bonds per conformer ----------
fig, ax = plt.subplots(figsize=(8, 5))
order = ["R_water", "S_water", "R_mem", "S_mem"]
labels = ["R\nwater", "S\nwater", "R\nmembrane", "S\nmembrane"]
box_data = [data[k][0] for k in order]
bp = ax.boxplot(box_data, labels=labels, patch_artist=True, widths=0.6,
                medianprops=dict(color="black", lw=2), showfliers=True,
                flierprops=dict(marker="o", markersize=3, alpha=0.3))
for patch, k in zip(bp["boxes"], order):
    patch.set_facecolor(C_R if k.startswith("R") else C_S)
    patch.set_alpha(0.65)
# overlay Boltzmann-weighted mean as a diamond
for i, k in enumerate(order, 1):
    hb, _, w = data[k]
    ax.plot(i, (w * hb).sum(), "D", color="gold", markeredgecolor="black", markersize=9, zorder=5)
ax.set_ylabel("intramolecular H-bonds per conformer")
ax.set_title("DOPC 3-12-8-12 R vs S — H-bond distribution per ensemble\n"
             "(gold diamond = Boltzmann mean)", fontweight="bold")
ax.text(0.5, 0.02, "R opens in water (low, spread) - S stays closed (high) - both closed in membrane",
        transform=ax.transAxes, ha="center", fontsize=9, style="italic", color="grey")
fig.tight_layout()
fig.savefig(OUT / "box_hbonds.svg"); fig.savefig(OUT / "box_hbonds.png", dpi=160)
plt.close(fig)

# ---------- Figure 2: relative difference 2D vs 3D ----------
def reldiff(a, b):
    denom = (abs(a) + abs(b)) / 2
    return 0.0 if denom == 0 else abs(a - b) / denom * 100

# 2D descriptors (identical) from CSV
d2 = pd.read_csv("results/2d_descriptors.csv").set_index("name")
iso = [n for n in d2.index if "3-12-8-12" in n]
two_d = {}
for col in ["TPSA_2d", "MolLogP_Crippen", "MolWt", "NumHDonors", "NumHAcceptors",
            "FractionCSP3", "MolMR_Crippen", "LabuteASA"]:
    two_d[col] = reldiff(d2.loc[iso[0], col], d2.loc[iso[1], col])

# 3D ensemble descriptors -- read from the robust ensemble_descriptors.py output
ed = pd.read_csv("results/ensemble_descriptors_dopc_rs.csv").set_index("compound")
r3d, s3d = ed.loc["DOPC_R"], ed.loc["DOPC_S"]
three_d_cols = [
    "water_bw_psa", "water_bw_IMHB", "water_bw_rg", "water_bw_spherocity",
    "water_bw_asphericity", "water_p_dominant", "mem_bw_psa", "mem_bw_IMHB",
]
three_d = {c: reldiff(float(r3d[c]), float(s3d[c])) for c in three_d_cols}

fig, ax = plt.subplots(figsize=(8, 6))
names = list(two_d) + list(three_d)
vals = list(two_d.values()) + list(three_d.values())
colors = ["#9aa0a6"] * len(two_d) + [C_S] * len(three_d)
y = np.arange(len(names))[::-1]
ax.barh(y, vals, color=colors, edgecolor="white")
ax.set_yticks(y); ax.set_yticklabels(names)
ax.set_xlabel("relative difference  %|R - S|  /  mean")
ax.set_title("What distinguishes the isomers: 2D vs 3D descriptors", fontweight="bold")
ax.axvline(0, color="black", lw=0.8)
# group separators / legend
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#9aa0a6", label="2D / lipophilicity (identical -> 0%)"),
                   Patch(color=C_S, label="3D ensemble (Boltzmann)")], loc="upper right")
ax.set_xlim(0, 112)
for yi, v in zip(y, vals):
    ax.text(v + 1.5, yi, f"{v:.0f}%", va="center", fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "rel_diff_2d_vs_3d.svg"); fig.savefig(OUT / "rel_diff_2d_vs_3d.png", dpi=160)
plt.close(fig)

# ---------- Figure 3: 2D descriptor overlap (R and S coincide exactly) ----------
# canonical druglike / lipophilicity set; LogP is RDKit Crippen (WLOGP) -- the
# only logP RDKit ships. Any 2D logP (Crippen/aLogP/cLogP) is stereo-blind and
# would be identical for R/S, so the method choice does not affect the conclusion.
cols = ["MolWt", "TPSA_2d", "MolLogP_Crippen", "MolMR_Crippen",
        "NumHDonors", "NumHAcceptors", "NumRotatableBonds"]
labels2d = ["MolWt", "TPSA", "LogP\n(Crippen)", "MolMR",
            "HBD", "HBA", "RotBonds"]
rvals = [d2.loc[iso[0], c] for c in cols]
svals = [d2.loc[iso[1], c] for c in cols]
x = np.arange(len(cols))
fig, ax = plt.subplots(figsize=(9, 5))
ax.bar(x, rvals, width=0.6, color=C_R, alpha=0.6, label="DOPC_R", edgecolor="none")
ax.bar(x, svals, width=0.6, color=C_S, alpha=0.6, label="DOPC_S", edgecolor="none")
ax.set_yscale("symlog")          # handles the negative LogP and the large MolWt range
ax.axhline(0, color="black", lw=0.6)
ax.set_xticks(x); ax.set_xticklabels(labels2d)
ax.set_ylabel("descriptor value (symlog scale)")
ax.set_title("2D / lipophilicity descriptors: R and S overlap exactly\n"
             "(bars perfectly coincide -> identical for every descriptor)", fontweight="bold")
ax.legend(loc="upper right")
for xi, v in zip(x, rvals):
    ax.text(xi, v + (0.15 if v >= 0 else -0.35), f"{v:g}", ha="center",
            va="bottom" if v >= 0 else "top", fontsize=8, color="grey")
fig.tight_layout()
fig.savefig(OUT / "overlap_2d.svg"); fig.savefig(OUT / "overlap_2d.png", dpi=160)
plt.close(fig)

print("Saved box_hbonds, rel_diff_2d_vs_3d, overlap_2d ->", OUT)
print("\n2D relative diffs:", {k: round(v,1) for k,v in two_d.items()})
print("3D relative diffs:", {k: round(v,1) for k,v in three_d.items()})

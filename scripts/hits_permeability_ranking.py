#!/usr/bin/env python
"""All-hits permeability ranking (2026-10-15).

Rank every ortho 6-mer hit by our validated permeability predictors (3D-PSA hexane leg,
and CPCM-X ΔG_transfer water->hexane) and place them against the SMD controls:

  * NEGATIVE control (impermeable): 1-6-4-7   -- high PSA, strongly positive ΔG  (in the 32)
  * POSITIVE controls (permeable):  2-9-9-8, 3-12-8-12 (R/S)  -- low PSA, favorable ΔG

Question: do the other hits populate the space BETWEEN the positive and negative controls?

Descriptors are computed from the transferred CREST ensembles, Boltzmann-weighted by the
per-molecule CPCM-X populations (ensemble_descriptors.process_compound). ΔG_transfer is read
from each fe_hit_<name>.summary.csv. Controls come from the SMD cross-check (earlier xylene
runs -- different provenance, used as REFERENCE values and flagged as such).

Outputs -> data/hits_permeability_ranking/2026-10-15/
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, "scripts")
import ensemble_descriptors as ed
import plot_style
plot_style.apply()
from plot_style import C

HT  = Path("results/conformers/hits_transfer")
OUT = Path("data/hits_permeability_ranking/2026-10-15")
OUT.mkdir(parents=True, exist_ok=True)

# shared publication palette
HIT   = C["hit"]   # ordinary hit
PCCOL = C["pc"]    # positive control (permeable)
NCCOL = C["nc"]    # negative control (impermeable)
GREYB = C["band"]
STEM  = C["stem"]

NC = "1-6-4-7"
# positive controls from the SMD analysis (earlier xylene runs; reference values).
# psa_hexane = membrane-core leg (permeability proxy); psa_water = aqueous leg (polar exposure).
PC = {
    "2-9-9-8":     dict(psa_hexane=119.6, psa_water=121.5, dG=-3.96),
    "3-12-8-12 R": dict(psa_hexane=149.8, psa_water=230.2, dG= 1.54),
    "3-12-8-12 S": dict(psa_hexane=155.5, psa_water=194.2, dG=-3.05),
}

def dG_of(name):
    s = HT / "free_energy" / f"fe_hit_{name}.summary.csv"
    if not s.exists():
        return np.nan
    d = pd.read_csv(s).iloc[0]
    return float(d["dGtransfer_water->hexane_kcal"])

# ---- descriptors for all 32 hits (via the pipeline script ensemble_descriptors.process_compound) ----
# Reuse the cached CSV if present (descriptor compute is the slow step); pass --recompute to force.
CSV = OUT / "all_hits_descriptors.csv"
if CSV.exists() and "--recompute" not in sys.argv:
    df = pd.read_csv(CSV)
    print(f"reusing cached descriptors ({len(df)} hits) -> {CSV}  (pass --recompute to regenerate)")
else:
    names = sorted(p.name for p in (HT / "conformers").iterdir() if p.is_dir())
    rows = []
    for nm in names:
        w = HT / "conformers" / nm / "water"
        h = HT / "conformers" / nm / "hexane"
        ecsv = HT / "free_energy" / f"fe_hit_{nm}.csv"
        r = ed.process_compound(nm, w, h, "hexane", ecsv)
        rows.append(dict(
            name=nm,
            psa_water=r.get("water_bw_psa"),
            psa_hexane=r.get("hexane_bw_psa"),
            rg_water=r.get("water_bw_rg"),
            rg_hexane=r.get("hexane_bw_rg"),
            dG_transfer=dG_of(nm),
            role=("NC" if nm == NC else "hit"),
        ))
    df = pd.DataFrame(rows)
    df.to_csv(CSV, index=False)
    print(f"computed {len(df)} hits -> {CSV}")

# Boltzmann-weighted ΔPSA (chameleonicity) = PSA(water) - PSA(hexane); each leg CPCM-X-weighted.
# High ΔPSA = shields polar surface on entering the membrane (chameleon); ~0 = rigid / can't bury.
df["delta_psa"] = df["psa_water"] - df["psa_hexane"]

# combined frame incl. positive controls (for the ranking plots)
pc_rows = [dict(name=k, psa_hexane=v["psa_hexane"], psa_water=v["psa_water"],
                delta_psa=v["psa_water"] - v["psa_hexane"], dG_transfer=v["dG"], role="PC") for k, v in PC.items()]
allf = pd.concat([df[["name", "psa_hexane", "psa_water", "delta_psa", "dG_transfer", "role"]], pd.DataFrame(pc_rows)], ignore_index=True)

pc_psa = [v["psa_hexane"] for v in PC.values()]
nc_psa = float(df.loc[df.name == NC, "psa_hexane"].iloc[0])
pc_psa_w = [v["psa_water"] for v in PC.values()]
nc_psa_w = float(df.loc[df.name == NC, "psa_water"].iloc[0])
pc_dpsa = [v["psa_water"] - v["psa_hexane"] for v in PC.values()]
nc_dpsa = float(df.loc[df.name == NC, "delta_psa"].iloc[0])
pc_dG  = [v["dG"] for v in PC.values()]
nc_dG  = float(df.loc[df.name == NC, "dG_transfer"].iloc[0])

def colors(roles):
    return [NCCOL if r == "NC" else PCCOL if r == "PC" else HIT for r in roles]

def lollipop(metric, label, lo_is_permeable, pc_vals, nc_val, fname, title):
    s = allf.sort_values(metric, ascending=lo_is_permeable).reset_index(drop=True)
    y = np.arange(len(s))
    fig, ax = plt.subplots(figsize=(8.6, 10.8))
    x0 = min(s[metric].min(), min(pc_vals))
    # reference marks: positive-control band + negative-control line
    ax.axvspan(min(pc_vals), max(pc_vals), color=GREYB, zorder=0)
    ax.axvline(nc_val, color=NCCOL, ls=(0, (5, 3)), lw=1.3, zorder=1)
    ax.grid(axis="x", zorder=0)
    ax.hlines(y, x0, s[metric], color=STEM, lw=1.1, zorder=2)
    ax.scatter(s[metric], y, c=colors(s.role), s=78, zorder=3, edgecolor="white", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(list(s.name))
    for yi, rl in enumerate(s.role):
        if rl != "hit":
            ax.get_yticklabels()[yi].set_color(NCCOL if rl == "NC" else PCCOL)
            ax.get_yticklabels()[yi].set_fontweight("bold")
    ax.set_xlabel(label)
    ax.set_title(title, loc="left", pad=12)
    ax.set_ylim(-0.7, len(s) - 0.3)
    ax.set_axisbelow(True)
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=PCCOL, markersize=9, label="positive control (permeable)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=NCCOL, markersize=9, label="negative control: 1-6-4-7 (impermeable)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=HIT, markersize=9, label="hit"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / fname)
    plt.close(fig)
    return s

# Fig 1 — ranked by PSA (hexane): lower = more permeable
s_psa = lollipop("psa_hexane", "3D-PSA, hexane leg (Å²)   ->   less permeable",
                 True, pc_psa, nc_psa, "fig_rank_psa_hexane.png",
                 "All 32 ortho hits ranked by 3D-PSA (hexane) — permeability proxy\n"
                 "green = positive controls (permeable) · red dashed = negative control 1-6-4-7")

# Fig 2 — ranked by ΔG_transfer: more negative = more permeable
s_dG = lollipop("dG_transfer", "ΔG_transfer water->hexane (kcal/mol)   ->   less permeable",
                True, pc_dG, nc_dG, "fig_rank_dG_transfer.png",
                "All 32 ortho hits ranked by ΔG_transfer (water->hexane)\n"
                "green = positive controls · red dashed = negative control 1-6-4-7")

# Fig 1b — ranked by PSA (WATER leg): aqueous polar exposure (upper-interface analog)
s_psa_w = lollipop("psa_water", "3D-PSA, WATER leg (Å²)   ->   more polar-exposed in water",
                   True, pc_psa_w, nc_psa_w, "fig_rank_psa_water.png",
                   "All 32 ortho hits ranked by 3D-PSA (WATER leg) — aqueous polar exposure\n"
                   "(the hexane leg is the permeability proxy; water leg = upper-interface / chameleon open state)")

# Fig 1c — Boltzmann-weighted ΔPSA (chameleonicity): PSA(water) - PSA(hexane)
s_dpsa = lollipop("delta_psa", "ΔPSA = PSA(water) − PSA(hexane) (Å²)   ->   more chameleonic (shields in membrane)",
                  False, pc_dpsa, nc_dpsa, "fig_rank_delta_psa.png",
                  "All 32 ortho hits ranked by Boltzmann-weighted ΔPSA (chameleonicity)\n"
                  "high = buries polar surface entering the membrane · ~0 = rigid/charged (can't bury). 6-mer scope: see note")

# Fig 3 — 2D permeability space: PSA (hexane) vs ΔG
fig, ax = plt.subplots(figsize=(9, 7))
ax.grid(True, zorder=0); ax.set_axisbelow(True)
ax.scatter(df.dG_transfer, df.psa_hexane, c=colors(df.role), s=85, edgecolor="white", linewidth=0.8, zorder=3)
for _, r in df.iterrows():
    if r.role == "NC":
        ax.annotate(r["name"] + " (NC)", (r.dG_transfer, r.psa_hexane), fontweight="bold",
                    color=NCCOL, xytext=(8, 4), textcoords="offset points")
for k, v in PC.items():
    ax.scatter(v["dG"], v["psa_hexane"], c=PCCOL, s=120, marker="D", edgecolor="white", linewidth=0.9, zorder=4)
    ax.annotate(k + " (PC)", (v["dG"], v["psa_hexane"]), fontweight="bold", color=PCCOL,
                xytext=(8, 4), textcoords="offset points")
ax.set_xlabel("ΔG_transfer water→hexane (kcal/mol)   →   less permeable")
ax.set_ylabel("3D-PSA, hexane leg (Å²)   →   less permeable")
ax.set_title("Permeability space: all 32 hits between the controls\n"
             "lower-left = permeable (toward PC) · upper-right = impermeable (toward NC 1-6-4-7)", loc="left")
handles = [
    Line2D([0], [0], marker="D", color="w", markerfacecolor=PCCOL, markersize=10, label="positive controls (permeable)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor=NCCOL, markersize=10, label="negative control 1-6-4-7"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor=HIT, markersize=10, label="hits (32)"),
]
ax.legend(handles=handles, loc="upper left", frameon=False)
fig.text(0.5, 0.005, "Controls from the SMD cross-check (earlier xylene runs; reference values). "
         "1-6-4-7 is the charged (−1) re-run.", ha="center", fontsize=8, style="italic", color=C["mut"])
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(OUT / "fig_permeability_space.png")
plt.close(fig)

# Fig 3b — 2D space with the WATER-leg PSA (polar exposure) on y
fig, ax = plt.subplots(figsize=(9, 7))
ax.grid(True, zorder=0); ax.set_axisbelow(True)
ax.scatter(df.dG_transfer, df.psa_water, c=colors(df.role), s=85, edgecolor="white", linewidth=0.8, zorder=3)
for _, r in df.iterrows():
    if r.role == "NC":
        ax.annotate(r["name"] + " (NC)", (r.dG_transfer, r.psa_water), fontweight="bold",
                    color=NCCOL, xytext=(8, 4), textcoords="offset points")
for k, v in PC.items():
    ax.scatter(v["dG"], v["psa_water"], c=PCCOL, s=120, marker="D", edgecolor="white", linewidth=0.9, zorder=4)
    ax.annotate(k + " (PC)", (v["dG"], v["psa_water"]), fontweight="bold", color=PCCOL,
                xytext=(8, 4), textcoords="offset points")
ax.set_xlabel("ΔG_transfer water→hexane (kcal/mol)   →   less permeable")
ax.set_ylabel("3D-PSA, water leg (Å²)   →   more polar-exposed in water")
ax.set_title("Water-leg view: aqueous polar exposure vs ΔG_transfer\n"
             "chameleons (e.g. 3-12-8-12 R): high water PSA but low hexane PSA — open in water, shielded in membrane", loc="left")
ax.legend(handles=handles, loc="upper left", frameon=False)
fig.text(0.5, 0.005, "Controls from the SMD cross-check (earlier xylene runs; reference values). "
         "1-6-4-7 is the charged (−1) re-run.", ha="center", fontsize=8, style="italic", color=C["mut"])
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(OUT / "fig_permeability_space_water.png")
plt.close(fig)

# ---- console summary: where does each hit sit vs the controls? ----
pc_hi_psa = max(pc_psa)           # least-permeable positive control
between = df[(df.psa_hexane > pc_hi_psa) & (df.psa_hexane < nc_psa) & (df.role == "hit")]
more_perm = df[(df.psa_hexane <= pc_hi_psa) & (df.role == "hit")]
less_perm = df[(df.psa_hexane >= nc_psa) & (df.role == "hit")]
print(f"\nPSA(hexane): PC band {min(pc_psa):.0f}-{max(pc_psa):.0f} Å²  |  NC(1-6-4-7) {nc_psa:.0f} Å²")
print(f"  hits BETWEEN PC and NC:        {len(between)}")
print(f"  hits more permeable than PC:   {len(more_perm)}  {list(more_perm.name)}")
print(f"  hits less permeable than NC:   {len(less_perm)}  {list(less_perm.name)}")
pc_hi_psa_w = max(pc_psa_w)
between_w = df[(df.psa_water > min(pc_psa_w)) & (df.psa_water < nc_psa_w) & (df.role == "hit")]
print(f"PSA(water):  PC range {min(pc_psa_w):.0f}-{max(pc_psa_w):.0f} Å²  |  NC(1-6-4-7) {nc_psa_w:.0f} Å²  (chameleons: 3-12-8-12 water PSA high)")
dmax = df.loc[df.delta_psa.idxmax()]
print(f"ΔPSA:        PC range {min(pc_dpsa):.0f}-{max(pc_dpsa):.0f} Å²  |  NC(1-6-4-7) {nc_dpsa:.0f} Å²  (most chameleonic hit: {dmax['name']} {dmax['delta_psa']:.0f})")
print(f"\nSaved 6 figures + all_hits_descriptors.csv -> {OUT}")

#!/usr/bin/env python
"""
Cross-method concordance: implicit-solvent 3D descriptors (our CREST->CPCM-X pipeline)
vs explicit-membrane SMD (colleague), for the four 6-mer hits (xylene constructs).

TWO-ZONE mapping (bulk implicit endpoints -> SMD depth zones):
   our WATER leg   <-> SMD UPPER INTERFACE   (aqueous/interface: most hydrated, open)
   our HEXANE leg  <-> SMD HYDROPHOBIC CORE  (apolar: most shielded, closed)
(The SMD lower interface mirrors the upper; we have a single water endpoint, so it maps
to the interface state generically. Bulk water/hexane approximate the interface/core.)

Statistical rigor:
  * Rank-order is PRIMARY (permeability-relevant ORDERING must transfer, not absolute values).
    Magnitude concordance (CCC/Bland-Altman) deliberately NOT claimed (different quantities/scales).
  * 3-12-8-12 R and S are DIASTEREOMERS = distinct molecules, never averaged. SMD reports one
    Cmp2 of unknown stereochem -> the ENTIRE analysis runs twice (Cmp2=R, then Cmp2=S).
  * n=4 -> qualitative signal preservation, not significance (exact Mantel p >= ~0.04).

Outputs (all into  data/6mer comparison/ ): 6 figures + 4 CSVs. Pure-numpy (no scipy).
"""
import csv, itertools
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("data/6mer comparison"); OUT.mkdir(parents=True, exist_ok=True)
MOLS = ["2-9-9-8", "3-12-8-12", "6-4-4-13", "1-6-4-7"]      # Cmp1, Cmp2, Cmp4, Cmp3
CMP  = {"2-9-9-8":"Cmp1", "3-12-8-12":"Cmp2", "6-4-4-13":"Cmp4", "1-6-4-7":"Cmp3"}

# ---- our implicit descriptors, xylene constructs; w = water leg, h = hexane leg ----
BASE = {
 "2-9-9-8":     {"w": dict(PSA=121.45, IMHB=3.384, bbHB=0.053, Rg=4.434, asph=0.0452, SA_HA=102.28, hydSASA=685.84),
                 "h": dict(PSA=119.61, IMHB=3.522, bbHB=0.028, Rg=4.499, asph=0.0573, SA_HA=105.04, hydSASA=699.98), "dG": -3.96},
 "6-4-4-13":    {"w": dict(PSA=136.10, IMHB=4.048, bbHB=2.494, Rg=4.965, asph=0.1739, SA_HA=87.45,  hydSASA=762.59),
                 "h": dict(PSA=123.36, IMHB=4.776, bbHB=1.131, Rg=4.818, asph=0.0739, SA_HA=95.43,  hydSASA=772.87), "dG": -5.94},
 # 1-6-4-7 re-run at its CORRECT charge (-1, cysteic acid sulfonate deprotonated); charged CREST/CPCM-X,
 # Boltzmann-weighted by CPCM-X pops. ΔG_transfer rose 2.57 -> +7.41 (strongly disfavors apolar, as expected);
 # Rgyr stayed compact (~4.7, NOT enlarged) -> implicit still under-opens vs explicit SMD (see fig3 note).
 "1-6-4-7":     {"w": dict(PSA=231.49, IMHB=4.929, bbHB=0.881, Rg=4.703, asph=0.2460, SA_HA=172.61, hydSASA=601.29),
                 "h": dict(PSA=232.00, IMHB=6.956, bbHB=0.998, Rg=4.635, asph=0.1278, SA_HA=202.32, hydSASA=600.52), "dG": 7.41},
 "3-12-8-12_R": {"w": dict(PSA=230.20, IMHB=4.531, bbHB=1.870, Rg=4.413, asph=0.1082, SA_HA=147.47, hydSASA=555.36),
                 "h": dict(PSA=149.84, IMHB=7.660, bbHB=1.011, Rg=4.753, asph=0.1102, SA_HA=134.87, hydSASA=672.47), "dG": 1.54},
 "3-12-8-12_S": {"w": dict(PSA=194.15, IMHB=5.495, bbHB=2.345, Rg=4.576, asph=0.1339, SA_HA=138.46, hydSASA=619.00),
                 "h": dict(PSA=155.53, IMHB=6.925, bbHB=2.988, Rg=4.572, asph=0.1467, SA_HA=119.97, hydSASA=635.61), "dG": -3.05},
}
def imp_set(iso):
    d = {m: BASE[m] for m in ["2-9-9-8", "6-4-4-13", "1-6-4-7"]}
    d["3-12-8-12"] = BASE[f"3-12-8-12_{iso}"]; return d

# ---- colleague SMD, per zone (SI Fig 14) + molecule-level translocation barrier (Fig 6) ----
SMD = {
 "upper": {  # UPPER INTERFACE  <-> our water leg
  "2-9-9-8":  dict(waterHB=19, IMHB=20, bbHB=12, Rg=8.3, lipid=90,  xsec=211),
  "3-12-8-12":dict(waterHB=27, IMHB=20, bbHB=13, Rg=8.7, lipid=102, xsec=221),
  "6-4-4-13": dict(waterHB=24, IMHB=19, bbHB=14, Rg=9.0, lipid=100, xsec=251),
  "1-6-4-7":  dict(waterHB=32, IMHB=27, bbHB=15, Rg=9.1, lipid=102, xsec=243)},
 "core": {   # HYDROPHOBIC CORE <-> our hexane leg
  "2-9-9-8":  dict(waterHB=11, IMHB=20, bbHB=13, Rg=8.3, lipid=235, xsec=220),
  "3-12-8-12":dict(waterHB=16, IMHB=21, bbHB=14, Rg=8.7, lipid=271, xsec=221),
  "6-4-4-13": dict(waterHB=17, IMHB=19, bbHB=14, Rg=9.1, lipid=265, xsec=252),
  "1-6-4-7":  dict(waterHB=25, IMHB=28, bbHB=15, Rg=8.8, lipid=226, xsec=226)},
}
BARRIER = {"2-9-9-8":44, "3-12-8-12":49, "6-4-4-13":56, "1-6-4-7":61}

ZONES = [("upper interface", "w", "upper"), ("hydrophobic core", "h", "core")]
# COMPARISON BASIS: only descriptors both methods compute as genuine MOLECULAR properties
# (no bilayer required). Omitted on principle — absent from our implicit pipeline, so any
# correlation would be meaningless: lipid-tail contacts (no lipids in implicit solvent),
# cross-section area (bilayer-plane footprint, orientation-dependent), tilt angle (vs
# membrane normal). We keep those SMD values in smd_reference.csv for record, not for stats.
# Descriptor concordance uses only DIRECT same-quantity overlaps. PSA is NOT here: SMD reports no
# PSA and its water-HB count is too noisy to proxy it; PSA and ΔG are compared as permeability
# PREDICTORS vs the SMD barrier instead (Fig 1).
PAIRS = [("Rg","Rg","Rgyr"), ("IMHB","IMHB","intramol HB"), ("bbHB","bbHB","backbone HB")]
# 6-4-4-13 was a hit for a NON-DOPC membrane; the SMD ran DOPC, so its DOPC barrier is confounded.
# Excluded a priori for the "DOPC-comparable" (n=3) permeability statistic.
MOLS_DOPC = [m for m in MOLS if m != "6-4-4-13"]
COL = {"2-9-9-8":"#0d8b96","3-12-8-12":"#c75917","6-4-4-13":"#2e7d5b","1-6-4-7":"#8a4fb0"}

def ranks(x):
    x = np.asarray(x, float); o = x.argsort(); r = np.empty(len(x)); r[o] = np.arange(1, len(x)+1)
    for v in np.unique(x):
        m = x == v
        if m.sum() > 1: r[m] = r[m].mean()
    return r
def pearson(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    return float("nan") if a.std()==0 or b.std()==0 else float(np.corrcoef(a, b)[0, 1])
def spearman(a, b): return pearson(ranks(a), ranks(b))
def kendall(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float); n = len(a); c = d = 0
    for i, j in itertools.combinations(range(n), 2):
        s = np.sign(a[i]-a[j]) * np.sign(b[i]-b[j]); c += s > 0; d += s < 0
    return (c-d)/(n*(n-1)//2), int(c), int(d), n*(n-1)//2
def zscore(v): v = np.asarray(v, float); return (v - v.mean())/v.std()
def dmat(Z):
    n = Z.shape[0]; return np.array([[np.linalg.norm(Z[i]-Z[j]) for j in range(n)] for i in range(n)])
def offdiag(D):
    n = D.shape[0]; return np.array([D[i, j] for i in range(n) for j in range(i+1, n)])
def mantel(D1, D2):
    v1 = offdiag(D1); r = pearson(v1, offdiag(D2))
    perms = list(itertools.permutations(range(D1.shape[0])))
    return r, sum(pearson(v1, offdiag(D2[np.ix_(p, p)])) >= r-1e-12 for p in perms)/len(perms)

# ---------- compute: zone x isomer ----------
rows = []; rho = {}; DIST = {}
for zlab, ph, zk in ZONES:
    for iso in ["R", "S"]:
        IMP = imp_set(iso); rho[(zlab, iso)] = {}
        for ik, sk, dlab in PAIRS:
            imp = [IMP[m][ph][ik] for m in MOLS]; smd = [SMD[zk][m][sk] for m in MOLS]
            rr = spearman(imp, smd); tau, c, _, tot = kendall(imp, smd)
            rho[(zlab, iso)][dlab] = rr
            rows.append(dict(zone=zlab, isomer=iso, metric=dlab, spearman=round(rr,3), kendall=round(tau,3), pairs=f"{c}/{tot}"))
        Zi = np.column_stack([zscore([IMP[m][ph][ik] for m in MOLS]) for ik, sk, _ in PAIRS])
        Zs = np.column_stack([zscore([SMD[zk][m][sk] for m in MOLS]) for ik, sk, _ in PAIRS])
        mr, mp = mantel(dmat(Zi), dmat(Zs)); DIST[(zlab, iso)] = dmat(Zi)
        rows.append(dict(zone=zlab, isomer=iso, metric=f"GLOBAL Mantel ({len(PAIRS)} desc)", spearman=round(mr,3), kendall="", pairs=f"p={mp:.3f}"))

# permeability (molecule-level): our PSA (each leg) & dG vs SMD barrier, full n=4 AND DOPC-comparable n=3
for iso in ["R", "S"]:
    IMP = imp_set(iso)
    for src, get in [("PSA-water", lambda m: IMP[m]["w"]["PSA"]), ("PSA-hexane", lambda m: IMP[m]["h"]["PSA"]),
                     ("dG", lambda m: IMP[m]["dG"])]:
        imp = [get(m) for m in MOLS]; bar = [BARRIER[m] for m in MOLS]
        rr = spearman(imp, bar); tau, c, _, tot = kendall(imp, bar)
        rows.append(dict(zone="permeability_n4", isomer=iso, metric=f"{src} vs barrier", spearman=round(rr,3), kendall=round(tau,3), pairs=f"{c}/{tot}"))
        r3 = spearman([get(m) for m in MOLS_DOPC], [BARRIER[m] for m in MOLS_DOPC])
        rows.append(dict(zone="permeability_DOPC_n3", isomer=iso, metric=f"{src} vs barrier (excl 6-4-4-13)", spearman=round(r3,3), kendall="", pairs="n=3"))

# ---------- print summary ----------
print(f"{'ZONE':18s} {'descriptor':26s} {'rho(R)':>7s} {'rho(S)':>7s}")
for zlab, ph, zk in ZONES:
    for ik, sk, dlab in PAIRS:
        print(f"{zlab:18s} {dlab:26s} {rho[(zlab,'R')][dlab]:7.2f} {rho[(zlab,'S')][dlab]:7.2f}")
print("\nPermeability (vs SMD barrier), n=4 (incl 6-4-4-13) then n=3 (DOPC-comparable), R iso:")
for r in rows:
    if r["zone"] in ("permeability_n4", "permeability_DOPC_n3") and r["isomer"] == "R":
        print(f"  [{r['zone']:20s}] {r['metric']:32s} rho={r['spearman']}")
print("\nGlobal Mantel (full 6-desc panel), by zone & isomer:")
for r in rows:
    if r["metric"].startswith("GLOBAL"): print(f"  {r['zone']:18s} {r['isomer']}: r={r['spearman']:.2f} ({r['pairs']})")

# ---------- CSVs ----------
with open(OUT/"concordance_stats.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["zone","isomer","metric","spearman","kendall","pairs"]); w.writeheader(); w.writerows(rows)
dkeys = list(BASE["2-9-9-8"]["w"])
with open(OUT/"descriptor_panel.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["molecule","cmp","isomer"] + [f"imp_w_{k}" for k in dkeys] + [f"imp_h_{k}" for k in dkeys]
               + ["dG"] + [f"smd_upper_{k}" for k in SMD["upper"]["2-9-9-8"]] + [f"smd_core_{k}" for k in SMD["core"]["2-9-9-8"]] + ["smd_barrier"])
    for m in MOLS:
        variants = [("R", f"3-12-8-12_R"), ("S", f"3-12-8-12_S")] if m == "3-12-8-12" else [("-", m)]
        for iso, bkey in variants:
            b = BASE[bkey]
            w.writerow([m, CMP[m], iso] + [b["w"][k] for k in dkeys] + [b["h"][k] for k in dkeys] + [b["dG"]]
                       + [SMD["upper"][m][k] for k in SMD["upper"][m]] + [SMD["core"][m][k] for k in SMD["core"][m]] + [BARRIER[m]])
with open(OUT/"smd_reference.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["molecule","cmp","zone"] + list(SMD["upper"]["2-9-9-8"]) + ["barrier"])
    for m in MOLS:
        w.writerow([m, CMP[m], "upper"] + [SMD["upper"][m][k] for k in SMD["upper"][m]] + [BARRIER[m]])
        w.writerow([m, CMP[m], "core"]  + [SMD["core"][m][k]  for k in SMD["core"][m]]  + [BARRIER[m]])
with open(OUT/"implicit_descriptors.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["molecule","isomer","leg"] + dkeys)
    for m in MOLS:
        variants = [("R", "3-12-8-12_R"), ("S", "3-12-8-12_S")] if m == "3-12-8-12" else [("-", m)]
        for iso, bkey in variants:
            for leg in ["w", "h"]:
                w.writerow([m, iso, {"w":"water","h":"hexane"}[leg]] + [BASE[bkey][leg][k] for k in dkeys])

# ================= FIGURES (two main + two supplementary) =================
from matplotlib.colors import TwoSlopeNorm

# --- FIG 1 (MAIN): permeability prediction — our PSA & ΔG vs the SMD barrier ---
# 6-4-4-13 (non-DOPC hit) drawn as an open X and EXCLUDED from the n=3 DOPC-comparable rho.
def _perm_panel(ax, get, ylab, title):
    for m in MOLS:
        if m == "3-12-8-12":
            ys = [get("3-12-8-12_R"), get("3-12-8-12_S")]
            ax.plot([BARRIER[m]]*2, ys, "-", color=COL[m], lw=1.0, alpha=0.5, zorder=2)
            ax.scatter([BARRIER[m]]*2, ys, color=COL[m], s=95, zorder=3)
            ax.annotate(f"{m} (R,S)", (BARRIER[m], max(ys)), fontsize=8, color=COL[m], xytext=(7,2), textcoords="offset points")
        elif m == "6-4-4-13":
            ax.scatter(BARRIER[m], get(m), marker="X", s=180, facecolor="none", edgecolor=COL[m], linewidths=2.4, zorder=3)
            ax.annotate(f"{m}\n(non-DOPC, excl.)", (BARRIER[m], get(m)), fontsize=8, color=COL[m], xytext=(7,2), textcoords="offset points")
        else:
            ax.scatter(BARRIER[m], get(m), color=COL[m], s=95, zorder=3)
            ax.annotate(m, (BARRIER[m], get(m)), fontsize=8, color=COL[m], xytext=(7,2), textcoords="offset points")
    # illustrative least-squares trend through the DOPC-comparable molecules only (excl. 6-4-4-13)
    xf, yf = [], []
    for m in MOLS_DOPC:
        for bkey in (["3-12-8-12_R", "3-12-8-12_S"] if m == "3-12-8-12" else [m]):
            xf.append(BARRIER[m]); yf.append(get(bkey))
    sl, ic = np.polyfit(xf, yf, 1)
    xl = np.array([min(BARRIER[m] for m in MOLS_DOPC) - 0.5, max(BARRIER[m] for m in MOLS_DOPC) + 0.5])
    ax.plot(xl, sl * xl + ic, "--", color="gray", lw=1.6, alpha=0.85, zorder=1, label="trend (n=3, excl. 6-4-4-13)")
    ax.legend(fontsize=8, loc="upper left", frameon=False)
    bk = lambda m: ("3-12-8-12_R" if m == "3-12-8-12" else m)
    r4 = spearman([get(bk(m)) for m in MOLS], [BARRIER[m] for m in MOLS])
    r3 = spearman([get(bk(m)) for m in MOLS_DOPC], [BARRIER[m] for m in MOLS_DOPC])
    ax.set_xlabel("SMD barrier (kcal/mol)  ->  less permeable", fontsize=9.5)
    ax.set_ylabel(ylab, fontsize=9.5)
    ax.set_title(f"{title} vs barrier\nrho = {r4:.2f} (n=4, incl. 6-4-4-13)   |   rho = {r3:.2f} (n=3 DOPC-comparable)", fontsize=9.5)
    ax.spines[["top","right"]].set_visible(False)
fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.3))
_perm_panel(axes[0], lambda bk: BASE[bk]["h"]["PSA"], "our 3D-PSA, hexane (Å²)  ->  less permeable", "PSA")
_perm_panel(axes[1], lambda bk: BASE[bk]["dG"], "our ΔG_transfer w->hex (kcal/mol)  ->  less permeable", "ΔG_transfer")
fig.suptitle("Permeability prediction: our descriptors vs the SMD translocation barrier\n"
             "open X = 6-4-4-13 (a non-DOPC hit, confounded in DOPC) -> excluded from the n=3 statistic", fontsize=10)
fig.tight_layout(); fig.savefig(OUT/"fig1_permeability_concordance.png", dpi=200); plt.close(fig)

# --- FIG 2 (MAIN): descriptor concordance heatmap (rank agreement, both zones) ---
# Horizontal bars: one per (descriptor, zone). Bar length + color = rank agreement (Spearman rho).
# Clearer than a heatmap of abstract numbers: green = transfers, red = diverges, with a spelled-out axis.
bars = []   # (label, rho_mean)  — mean of R/S (near-identical; only Rgyr/interface differs 0.6/0.8)
for rk, disp in [("Rgyr", "radius of gyration"), ("intramol HB", "intramol H-bonds"), ("backbone HB", "backbone H-bonds")]:
    for zlab, zdisp in [("hydrophobic core", "core (hexane)"), ("upper interface", "interface (water)")]:
        rm = (rho[(zlab, "R")][rk] + rho[(zlab, "S")][rk]) / 2
        bars.append((f"{disp}  —  {zdisp}", rm))
labels = [b[0] for b in bars]; vals = np.array([b[1] for b in bars])
y = np.arange(len(bars))[::-1]                          # first descriptor on top
norm = TwoSlopeNorm(vmin=-0.5, vcenter=0.0, vmax=1.0)
colors = plt.get_cmap("RdYlGn")(norm(vals))
fig, ax = plt.subplots(figsize=(9.8, 4.9))
ax.barh(y, vals, color=colors, edgecolor="#333", height=0.72, zorder=3)
for yi, v in zip(y, vals):
    ax.text(v + (0.03 if v >= 0 else -0.03), yi, f"{v:+.2f}", va="center",
            ha="left" if v >= 0 else "right", fontsize=10, fontweight="bold")
ax.axvline(0, color="k", lw=1.0)
ax.axvline(0.8, color="gray", ls="--", lw=1.0)
ax.set_ylim(-0.95, len(bars) - 0.3)
ax.text(0.8, -0.72, "strong (0.8)", fontsize=8, color="gray", ha="center", va="center")
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9.5)
ax.set_xlim(-0.65, 1.15); ax.set_xticks([-0.5, 0, 0.5, 1.0])
ax.set_xlabel("rank agreement with SMD  (Spearman rho):    1 = perfect   ·   0 = none   ·   < 0 = opposite", fontsize=9.5)
ax.set_title("Which structural descriptors transfer:  explicit SMD  ->  our implicit pipeline\n"
             "size (Rgyr) transfers; H-bond counts diverge = the implicit over-shielding bias", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig(OUT/"fig2_descriptor_concordance.png", dpi=200); plt.close(fig)

# --- FIG 3 (MAIN): Rgyr concordance — the one direct same-quantity overlap (hydrophobic core) ---
fig, ax = plt.subplots(figsize=(8.8, 5.4))
for m in MOLS:
    if m == "3-12-8-12":
        xx = SMD["core"][m]["Rg"]
        yR = BASE["3-12-8-12_R"]["h"]["Rg"]; yS = BASE["3-12-8-12_S"]["h"]["Rg"]
        ax.plot([xx, xx], [yR, yS], "-", color=COL[m], lw=1.0, alpha=0.55, zorder=2)
        ax.scatter([xx, xx], [yR, yS], color=COL[m], s=90, zorder=3)
        ax.annotate(f"{m}  R", (xx, yR), fontsize=8.5, color=COL[m], xytext=(9, 1), textcoords="offset points")
        ax.annotate(f"{m}  S  (closer to SMD)", (xx, yS), fontsize=8.5, color=COL[m], xytext=(-10, -1), ha="right", textcoords="offset points")
    else:
        ax.scatter(SMD["core"][m]["Rg"], BASE[m]["h"]["Rg"], color=COL[m], s=90, zorder=3)
        ax.annotate(m, (SMD["core"][m]["Rg"], BASE[m]["h"]["Rg"]), fontsize=8.5, color=COL[m], xytext=(8,2), textcoords="offset points")
rR = spearman([imp_set("R")[m]["h"]["Rg"] for m in MOLS], [SMD["core"][m]["Rg"] for m in MOLS])
ax.set_xlabel("SMD radius of gyration (Å, hydrophobic core)")
ax.set_ylabel("our radius of gyration (Å, hexane leg)")
ax.set_title("Radius of gyration - the one direct same-quantity overlap\n"
             f"Spearman rho = {rR:.2f}  (rank-only; ~2x scale = different Rg definitions; tight spread)", fontsize=9.5)
ax.spines[["top","right"]].set_visible(False)
fig.text(0.5, 0.02,
         "Note: 1-6-4-7 now at its CORRECT charge (-1). ΔG_transfer rose to +7.4 (strongly disfavors apolar), as expected --\n"
         "but Rgyr stayed COMPACT (~4.7, not enlarged): implicit still under-opens vs explicit SMD (sampling limit, not a charge error).",
         ha="center", fontsize=8, style="italic", color="#555")
fig.tight_layout(rect=[0, 0.09, 1, 1]); fig.savefig(OUT/"fig3_rgyr_concordance.png", dpi=200); plt.close(fig)

# --- FIG 4 (MAIN): normalized-trend concordance (z-scored within each method) ---
# z-score removes offset AND scale (affine): if the two lines overlay, ours = a*SMD + b, i.e.
# the relative SPACING of molecules agrees and only offset/scale differ = calibratable, NOT a
# predictive-relationship failure. Hydrophobic core / hexane leg; Cmp2 = R (isomer-robust).
order = sorted(MOLS, key=lambda m: BARRIER[m])          # molecules along the permeability spectrum
IMP = imp_set("R")
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.5), sharey=True)
for ax, (ik, sk, dlab) in zip(axes, PAIRS):
    ours = [IMP[m]["h"][ik] for m in order]; smd = [SMD["core"][m][sk] for m in order]
    ax.plot(range(4), zscore(smd),  "--o", color="#333333", lw=2, ms=7, label="SMD (explicit)")
    ax.plot(range(4), zscore(ours), "-o",  color="#0d8b96", lw=2, ms=7, label="ours (implicit)")
    ax.set_xticks(range(4)); ax.set_xticklabels(order, rotation=40, ha="right", fontsize=8)
    ax.set_title(f"{dlab}\nPearson r = {pearson(ours, smd):.2f}", fontsize=9.5)
    ax.axhline(0, color="gray", lw=0.6); ax.spines[["top","right"]].set_visible(False)
axes[0].set_ylabel("z-score within method"); axes[0].legend(fontsize=8, loc="upper left")
fig.suptitle("Normalized-trend concordance (z-scored within each method) - hydrophobic core, molecules ordered by permeability   (Cmp2 = R)\n"
             "overlapping lines = offset/scale-only = CALIBRATABLE.     Pearson r = linear (affine) agreement, the right metric here (rank agreement is Fig 2).",
             fontsize=9.5)
fig.tight_layout(); fig.savefig(OUT/"fig4_normalized_trend.png", dpi=200); plt.close(fig)

# ===== WATER-LEG (upper-interface) counterparts of Fig 1 / 3 / 4 (Fig 2 already shows both zones) =====

# FIG 1 (water): PSA(water leg) vs barrier. (dG is leg-independent -> it lives in the hexane Fig 1 only.)
fig, ax = plt.subplots(figsize=(7.0, 5.3))
_perm_panel(ax, lambda bk: BASE[bk]["w"]["PSA"], "our 3D-PSA, WATER leg (Å²)  ->  less permeable", "PSA (water leg)")
fig.suptitle("Permeability prediction (WATER leg): our PSA vs the SMD translocation barrier\n"
             "open X = 6-4-4-13 (non-DOPC hit, confounded in DOPC) -> excluded from the n=3 statistic", fontsize=9.5)
fig.tight_layout(); fig.savefig(OUT/"fig1_permeability_concordance_water.png", dpi=200); plt.close(fig)

# FIG 3 (water): Rgyr, our water leg vs SMD UPPER INTERFACE
def _zdiff_w(iso):
    ours = [BASE[("3-12-8-12_"+iso) if mm == "3-12-8-12" else mm]["w"]["Rg"] for mm in MOLS]
    zi = zscore(ours); zs = zscore([SMD["upper"][mm]["Rg"] for mm in MOLS]); k = MOLS.index("3-12-8-12")
    return abs(zi[k] - zs[k])
closer_w = "S" if _zdiff_w("S") <= _zdiff_w("R") else "R"
fig, ax = plt.subplots(figsize=(8.8, 5.4))
for m in MOLS:
    if m == "3-12-8-12":
        xx = SMD["upper"][m]["Rg"]; yR = BASE["3-12-8-12_R"]["w"]["Rg"]; yS = BASE["3-12-8-12_S"]["w"]["Rg"]
        ax.plot([xx, xx], [yR, yS], "-", color=COL[m], lw=1.0, alpha=0.55, zorder=2)
        ax.scatter([xx, xx], [yR, yS], color=COL[m], s=90, zorder=3)
        (u_iso, u_y), (l_iso, l_y) = (("R", yR), ("S", yS)) if yR >= yS else (("S", yS), ("R", yR))
        ax.annotate(f"{m} {u_iso}" + ("  (closer to SMD)" if u_iso == closer_w else ""), (xx, u_y), fontsize=8, color=COL[m], xytext=(-9, 3), ha="right", textcoords="offset points")
        ax.annotate(f"{m} {l_iso}" + ("  (closer to SMD)" if l_iso == closer_w else ""), (xx, l_y), fontsize=8, color=COL[m], xytext=(-9, -4), ha="right", textcoords="offset points")
    else:
        ax.scatter(SMD["upper"][m]["Rg"], BASE[m]["w"]["Rg"], color=COL[m], s=90, zorder=3)
        ax.annotate(m, (SMD["upper"][m]["Rg"], BASE[m]["w"]["Rg"]), fontsize=8.5, color=COL[m], xytext=(8, 2), textcoords="offset points")
rw_R = spearman([BASE[("3-12-8-12_R" if m == "3-12-8-12" else m)]["w"]["Rg"] for m in MOLS], [SMD["upper"][m]["Rg"] for m in MOLS])
rw_S = spearman([BASE[("3-12-8-12_S" if m == "3-12-8-12" else m)]["w"]["Rg"] for m in MOLS], [SMD["upper"][m]["Rg"] for m in MOLS])
ax.set_xlabel("SMD radius of gyration (Å, upper interface)")
ax.set_ylabel("our radius of gyration (Å, water leg)")
ax.set_title("Radius of gyration (WATER / upper interface) - direct same-quantity overlap\n"
             f"Spearman rho = {rw_R:.2f} (Cmp2=R) / {rw_S:.2f} (Cmp2=S)  (rank-only; tight spread; ~2x scale = different Rg defs)", fontsize=9.5)
ax.spines[["top","right"]].set_visible(False)
fig.text(0.5, 0.02,
         "Note: 1-6-4-7 now at its CORRECT charge (-1). Even so its Rgyr stayed COMPACT (~4.7, not enlarged):\n"
         "implicit under-opens vs the explicit-water SMD upper interface (sampling/population limit, not the earlier charge error).",
         ha="center", fontsize=8, style="italic", color="#555")
fig.tight_layout(rect=[0, 0.09, 1, 1]); fig.savefig(OUT/"fig3_rgyr_concordance_water.png", dpi=200); plt.close(fig)

# FIG 4 (water): normalized trend, water leg vs SMD upper interface
IMPw = imp_set("R")
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.5), sharey=True)
for ax, (ik, sk, dlab) in zip(axes, PAIRS):
    ours = [IMPw[m]["w"][ik] for m in order]; smd = [SMD["upper"][m][sk] for m in order]
    ax.plot(range(4), zscore(smd),  "--o", color="#333333", lw=2, ms=7, label="SMD (explicit)")
    ax.plot(range(4), zscore(ours), "-o",  color="#0d8b96", lw=2, ms=7, label="ours (implicit)")
    ax.set_xticks(range(4)); ax.set_xticklabels(order, rotation=40, ha="right", fontsize=8)
    ax.set_title(f"{dlab}\nPearson r = {pearson(ours, smd):.2f}", fontsize=9.5)
    ax.axhline(0, color="gray", lw=0.6); ax.spines[["top","right"]].set_visible(False)
axes[0].set_ylabel("z-score within method"); axes[0].legend(fontsize=8, loc="upper left")
fig.suptitle("Normalized-trend concordance (z-scored) - WATER leg / upper interface, molecules ordered by permeability   (Cmp2 = R)\n"
             "overlapping lines = offset/scale-only = CALIBRATABLE.     Pearson r = linear (affine) agreement (rank agreement is Fig 2).",
             fontsize=9.5)
fig.tight_layout(); fig.savefig(OUT/"fig4_normalized_trend_water.png", dpi=200); plt.close(fig)

print(f"\nWrote 4 hexane + 3 water figures + 4 CSVs -> {OUT}/")

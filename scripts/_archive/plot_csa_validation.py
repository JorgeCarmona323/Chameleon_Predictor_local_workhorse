# env: chameleon-calc
#!/usr/bin/env python3
"""
plot_csa_validation.py
-----------------------
Publication-quality figures for the CsA water ensemble A1 validation report.

Produces 3 figures:
  1. Boltzmann weight distribution — bars colored by # A1 criteria satisfied
  2. A1 criteria pass rates — Boltzmann-weighted % per criterion
  3. Omega dihedral heatmap — all backbone amide bonds across 23 conformers

Output: results/figures/csa_validation_*.png

Usage:
  python scripts/plot_csa_validation.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
from rdkit import Chem

# ── Import analysis functions from validate_csa_water ─────────────────────────
sys.path.insert(0, str(Path(__file__).parent))
from validate_csa_water import (
    compute_omega_dihedrals, assign_ring_nhs, check_nh_hbonds,
    _ring_atom_set, _carbonyl_oxygens, _preceding_carbonyl_o,
    CIS_CUTOFF, DATA_DIR,
)

# ── Paths ──────────────────────────────────────────────────────────────────────
SDF_PATH  = DATA_DIR / "ensemble.sdf"
JSON_PATH = DATA_DIR / "ensemble.json"
OUT_DIR   = Path("results/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Style ──────────────────────────────────────────────────────────────────────
COLORS = {
    0: "#d1d5db",   # gray   — 0 criteria
    1: "#93c5fd",   # light blue — 1 criterion
    2: "#3b82f6",   # blue   — 2 criteria
    3: "#f97316",   # orange — 3 criteria
    4: "#16a34a",   # green  — all 4 (A1-like)
}
plt.rcParams.update({
    "font.family":    "sans-serif",
    "font.size":      10,
    "axes.spines.top":   False,
    "axes.spines.right": False,
})

# ── Load data ──────────────────────────────────────────────────────────────────
suppl = Chem.SDMolSupplier(str(SDF_PATH), removeHs=False, sanitize=True)
mols  = [m for m in suppl if m is not None]

with open(JSON_PATH) as fh:
    jdata = json.load(fh)
boltz = np.array([c["boltzmannweight"] for c in jdata["conformers"]], dtype=float)
boltz /= boltz.sum()

ref_mol    = mols[0]
nh_map     = assign_ring_nhs(ref_mol)
ring_atoms = _ring_atom_set(ref_mol)
acceptors  = _carbonyl_oxygens(ref_mol, ring_atoms)

ala7_expected_acc = _preceding_carbonyl_o(ref_mol, nh_map.get("Ala7", -1), ring_atoms) \
                    if "Ala7" in nh_map else None
abu2_adj_acc      = _preceding_carbonyl_o(ref_mol, nh_map.get("Abu2", -1), ring_atoms) \
                    if "Abu2" in nh_map else None

# ── Per-conformer analysis ─────────────────────────────────────────────────────
rows = []
all_omega_matrices = []   # for heatmap

for i, mol in enumerate(mols):
    omegas   = compute_omega_dihedrals(mol)
    n_cis    = sum(1 for o in omegas if abs(o["omega"]) < CIS_CUTOFF)
    hb       = check_nh_hbonds(mol, nh_map, ring_atoms, acceptors)

    abu2_hb = hb.get("Abu2") is not None
    val5_hb = hb.get("Val5") is not None
    ala7_hb = hb.get("Ala7") is not None

    abu2_info = hb.get("Abu2")
    ala7_info = hb.get("Ala7")
    abu2_lr   = (abu2_hb and abu2_adj_acc is not None
                 and abu2_info["a_idx"] != abu2_adj_acc)

    criteria = [
        n_cis == 1,      # cis amide
        abu2_hb,         # Abu2 H-bond
        ala7_hb,         # Ala7 H-bond
        not val5_hb,     # Val5 exposed
    ]
    n_pass  = sum(criteria)
    a1_like = (n_cis == 1 and abu2_hb and ala7_hb and not val5_hb)

    rows.append({
        "conf":   i + 1,
        "boltz":  float(boltz[i]),
        "n_cis":  n_cis,
        "n_pass": n_pass,
        "a1":     a1_like,
        "crit":   criteria,
    })
    all_omega_matrices.append([o["omega"] for o in omegas])

n_bonds = max(len(r) for r in all_omega_matrices)

# Pad rows with fewer bonds (shouldn't happen but safety)
for row in all_omega_matrices:
    while len(row) < n_bonds:
        row.append(float("nan"))
omega_matrix = np.array(all_omega_matrices)   # shape (23, n_bonds)

# ── Figure 1: Boltzmann weight distribution ───────────────────────────────────
fig1, ax = plt.subplots(figsize=(10, 4.5))

conf_ids = [r["conf"] for r in rows]
weights  = [r["boltz"] * 100 for r in rows]
bar_cols = [COLORS[r["n_pass"]] for r in rows]

bars = ax.bar(conf_ids, weights, color=bar_cols, edgecolor="white", linewidth=0.5)

# Annotate dominant conformer
ax.annotate(
    f"C1: {weights[0]:.1f}%\n(3/4 criteria)",
    xy=(1, weights[0]), xytext=(3, weights[0] + 2),
    arrowprops=dict(arrowstyle="->", color="#374151"),
    fontsize=8.5, color="#374151",
)

ax.set_xlabel("Conformer", fontsize=11)
ax.set_ylabel("Boltzmann weight (%)", fontsize=11)
ax.set_title("CsA Water Ensemble — Boltzmann Weight Distribution\n"
             "Color = number of A1 fingerprint criteria satisfied",
             fontsize=11, pad=10)
ax.set_xticks(conf_ids)
ax.set_xticklabels([str(c) for c in conf_ids], fontsize=7.5)
ax.set_xlim(0.2, len(rows) + 0.8)

legend_patches = [
    mpatches.Patch(color=COLORS[n], label=f"{n}/4 criteria")
    for n in sorted(COLORS)
    if any(r["n_pass"] == n for r in rows)
]
ax.legend(handles=legend_patches, loc="upper right", fontsize=9, framealpha=0.9)
ax.text(0.98, 0.60,
        "All conformers fully trans\n(no cis amide found)",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=9, color="#dc2626",
        bbox=dict(boxstyle="round,pad=0.3", fc="#fef2f2", ec="#dc2626", alpha=0.9))

plt.tight_layout()
p1 = OUT_DIR / "csa_validation_boltzmann.svg"
fig1.savefig(p1)
plt.close(fig1)
print(f"Saved: {p1}")

# ── Figure 2: A1 criteria pass rates ─────────────────────────────────────────
fig2, ax = plt.subplots(figsize=(7.5, 4))

criteria_labels = [
    "1. Cis amide\n   (MeVal11-MeBmt1)",
    "2. Abu2 NH\n   H-bonded",
    "3. Ala7 NH\n   H-bonded",
    "4. Val5 NH\n   solvent-exposed",
]
pass_rates = [
    sum(r["boltz"] for r in rows if r["crit"][i]) * 100
    for i in range(4)
]
bar_colors = ["#dc2626" if i == 0 else "#3b82f6" for i in range(4)]

y_pos = np.arange(4)
hbars = ax.barh(y_pos, pass_rates, color=bar_colors, height=0.55,
                edgecolor="white", linewidth=0.5)

# Value labels
for bar, val in zip(hbars, pass_rates):
    ax.text(val + 1, bar.get_y() + bar.get_height() / 2,
            f"{val:.1f}%", va="center", fontsize=10, fontweight="bold")

ax.axvline(100, color="#9ca3af", linestyle="--", linewidth=1, alpha=0.7)
ax.set_yticks(y_pos)
ax.set_yticklabels(criteria_labels, fontsize=10)
ax.set_xlabel("Boltzmann-weighted population satisfying criterion (%)", fontsize=10)
ax.set_title("CsA A1 Fingerprint — Criterion Satisfaction\n"
             "Reference: Limbach et al., JACS 2022",
             fontsize=11, pad=10)
ax.set_xlim(0, 115)

fail_patch = mpatches.Patch(color="#dc2626", label="FAIL (0.0%) — sampling gap")
pass_patch = mpatches.Patch(color="#3b82f6", label="PASS — H-bond network correct")
ax.legend(handles=[fail_patch, pass_patch], loc="lower right", fontsize=9)

ax.text(0.98, 0.18,
        "Root cause: CREST cannot cross\n~15–20 kcal/mol cis-trans barrier.\nFix: v3.3 multi-start enumeration.",
        transform=ax.transAxes, ha="right", va="bottom",
        fontsize=8.5, color="#374151",
        bbox=dict(boxstyle="round,pad=0.4", fc="#f9fafb", ec="#d1d5db", alpha=0.95))

plt.tight_layout()
p2 = OUT_DIR / "csa_validation_criteria.svg"
fig2.savefig(p2)
plt.close(fig2)
print(f"Saved: {p2}")

# ── Figure 3: Omega dihedral heatmap ──────────────────────────────────────────
fig3, ax = plt.subplots(figsize=(11, 6))

# Normalise: cis ~0°, trans ~±180° → color diverges from 0
norm  = TwoSlopeNorm(vmin=-180, vcenter=0, vmax=180)
cmap  = plt.cm.RdBu   # red = positive trans, blue = negative trans, white = cis

im = ax.imshow(omega_matrix, aspect="auto", cmap=cmap, norm=norm,
               interpolation="nearest")

# y-axis: conformer labels with Boltzmann weight, sized by weight
y_labels = [f"C{r['conf']}  ({r['boltz']*100:.1f}%)" for r in rows]
ax.set_yticks(np.arange(len(rows)))
ax.set_yticklabels(y_labels, fontsize=7.5)

# x-axis: bond index
bond_labels = [f"Bond\n{j+1}" for j in range(n_bonds)]
ax.set_xticks(np.arange(n_bonds))
ax.set_xticklabels(bond_labels, fontsize=8)

cbar = plt.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
cbar.set_label("Omega dihedral (°)", fontsize=10)
cbar.set_ticks([-180, -90, 0, 90, 180])

ax.set_title("CsA Water Ensemble — Backbone Omega Dihedrals\n"
             "All bonds trans (|ω| > 150°); cis (|ω| < 30°) would appear white",
             fontsize=11, pad=10)
ax.set_xlabel("Backbone amide bond index", fontsize=10)
ax.set_ylabel("Conformer (Boltzmann weight)", fontsize=10)

# Overlay: mark the MeVal11-MeBmt1 bond column if identifiable
# (We can't know the exact column without residue labeling, so add a note)
ax.text(0.01, 0.01,
        "White = cis (|ω| < 30°) — absent in all conformers\n"
        "Blue/Red = trans; MeVal11-MeBmt1 bond not yet labeled",
        transform=ax.transAxes, fontsize=8, va="bottom",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#d1d5db", alpha=0.9))

plt.tight_layout()
p3 = OUT_DIR / "csa_validation_omega_heatmap.svg"
fig3.savefig(p3)
plt.close(fig3)
print(f"Saved: {p3}")

print("\nAll figures saved to results/figures/")
print("  csa_validation_boltzmann.svg   — Boltzmann weight distribution")
print("  csa_validation_criteria.svg    — A1 criterion pass rates")
print("  csa_validation_omega_heatmap.svg — Omega dihedral heatmap")

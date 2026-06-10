"""
make_isomer_figures.py
----------------------
Figure pipeline for the DOPC 3-12-8-12 R/S isomer comparison.

Produces:
  1. PDB files of the dominant (highest Boltzmann weight) conformer for each
     isomer x solvent -> results/figures/isomers/{R,S}_{water,mem}_dominant.pdb
     (load these in PyMOL via scripts/isomer_figures.pml)
  2. Matplotlib ensemble-distribution figures:
     - intramolecular H-bond distribution (R vs S, water vs mem)
     - conformer population profile (R diffuse vs S concentrated)
     - PSA distribution per solvent
     -> results/figures/isomers/*.svg/.png

See docs/experiments/2026-06-05_dopc_rs_3d_vs_2d_descriptors.md
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

BASE = Path("results/conformers")
OUT = Path("results/figures/isomers")
OUT.mkdir(parents=True, exist_ok=True)

CASES = {
    "R_water": BASE / "DOPC 3-12-8-12 R" / "water",
    "R_mem":   BASE / "DOPC 3-12-8-12 R" / "mem",
    "S_water": BASE / "DOPC 3-12-8-12 S" / "water",
    "S_mem":   BASE / "DOPC 3-12-8-12 S" / "mem",
}
COLORS = {"R": "#d1495b", "S": "#30638e"}


def load(case_dir):
    with open(case_dir / "ensemble.json") as f:
        data = json.load(f)
    confs = data["conformers"]
    w = np.array([c["boltzmannweight"] for c in confs], dtype=float)
    w = w / w.sum()
    psa = np.array([c["psa"] for c in confs], dtype=float)
    hb = np.array([c["hbonds"] for c in confs], dtype=float)
    mols = [m for m in Chem.SDMolSupplier(str(case_dir / "ensemble.sdf"), removeHs=False) if m]
    return w, psa, hb, mols


data = {k: load(v) for k, v in CASES.items()}

# 1. Dominant conformer PDBs ---------------------------------------------------
for key, (w, psa, hb, mols) in data.items():
    dom = int(np.argmax(w))
    m = mols[dom]
    pdb = OUT / f"{key}_dominant.pdb"
    Chem.MolToPDBFile(m, str(pdb))
    print(f"{key}: dominant conf #{dom} (p={w[dom]:.3f}, HB={int(hb[dom])}, PSA={psa[dom]:.0f}) -> {pdb.name}")

# 2a. H-bond distribution (weighted) -------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for ax, solvent in zip(axes, ["water", "mem"]):
    for iso in ["R", "S"]:
        w, psa, hb, _ = data[f"{iso}_{solvent}"]
        bins = np.arange(-0.5, hb.max() + 1.5, 1)
        ax.hist(hb, bins=bins, weights=w, alpha=0.55, color=COLORS[iso],
                label=f"DOPC_{iso}", edgecolor="white")
        ax.axvline((w * hb).sum(), color=COLORS[iso], ls="--", lw=2)
    ax.set_title(f"{solvent.upper()}  (dashed = Boltzmann mean)")
    ax.set_xlabel("intramolecular H-bonds per conformer")
    ax.legend()
axes[0].set_ylabel("Boltzmann-weighted population")
fig.suptitle("DOPC 3-12-8-12 R vs S — intramolecular H-bond distribution", fontweight="bold")
fig.tight_layout()
fig.savefig(OUT / "hbond_distribution.svg")
fig.savefig(OUT / "hbond_distribution.png", dpi=160)
plt.close(fig)

# 2b. Conformer population profile (sorted weights) ----------------------------
fig, ax = plt.subplots(figsize=(7, 4.5))
for iso in ["R", "S"]:
    w, *_ = data[f"{iso}_water"]
    sw = np.sort(w)[::-1]
    ax.plot(np.arange(1, len(sw) + 1), np.cumsum(sw), color=COLORS[iso], lw=2.2,
            label=f"DOPC_{iso}  (p_dominant={w.max():.2f}, n={len(w)})")
ax.set_xscale("log")
ax.set_xlabel("conformer rank (by Boltzmann weight, log scale)")
ax.set_ylabel("cumulative Boltzmann population")
ax.set_title("Water ensemble concentration: R diffuse vs S concentrated", fontweight="bold")
ax.axhline(0.9, color="grey", ls=":", lw=1)
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "population_profile.svg")
fig.savefig(OUT / "population_profile.png", dpi=160)
plt.close(fig)

# 2c. PSA distribution ----------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for ax, solvent in zip(axes, ["water", "mem"]):
    for iso in ["R", "S"]:
        w, psa, hb, _ = data[f"{iso}_{solvent}"]
        ax.hist(psa, bins=20, weights=w, alpha=0.55, color=COLORS[iso],
                label=f"DOPC_{iso}", edgecolor="white")
        ax.axvline((w * psa).sum(), color=COLORS[iso], ls="--", lw=2)
    ax.set_title(f"{solvent.upper()}")
    ax.set_xlabel("3D polar SASA (Å²)")
    ax.legend()
axes[0].set_ylabel("Boltzmann-weighted population")
fig.suptitle("DOPC 3-12-8-12 R vs S — PSA distribution", fontweight="bold")
fig.tight_layout()
fig.savefig(OUT / "psa_distribution.svg")
fig.savefig(OUT / "psa_distribution.png", dpi=160)
plt.close(fig)

print(f"\nFigures + PDBs saved -> {OUT}")

# env: MACE
"""
benchmark_mace_vs_xtb.py
------------------------
Compare MACE-OFF23(M) vs GFN2-xTB+ALPB Boltzmann weights on the 23 CsA water
conformers from the existing cregen ensemble.

If the energy ranking agrees well (Pearson r > ~0.85 on relative energies), MACE-OFF
is trustworthy for this molecule class and the GPU pipeline is viable.

If they diverge, the most likely culprit is missing solvation in MACE-OFF (gas-phase
vs ALPB-water), which is expected and interpretable — see the bottom section for context.

Usage:
    pip install mace-torch ase scipy matplotlib
    python scripts/benchmark_mace_vs_xtb.py

Outputs:
    Console table comparing per-conformer energies and Boltzmann weights
    results/mace_vs_xtb_CsA_water.png  (scatter plots, if matplotlib available)
    results/mace_vs_xtb_CsA_water.csv  (full numerical comparison)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--model", default="MACE-OFF23(M)",
                    help="MACE model name or path to local .model file")
args = parser.parse_args()

# ── Paths ─────────────────────────────────────────────────────────────────────
def _find_csa_water_ensemble() -> tuple[Path, Path]:
    candidates = [
        Path("data/CREST_CsA_20260512"),
        *sorted(Path("results/runs").glob("*_CsA/water"), reverse=True),
    ]
    for d in candidates:
        xyz = d / "ensemble.xyz"
        jsn = d / "ensemble.json"
        if xyz.exists() and jsn.exists():
            return xyz, jsn
    sys.exit("ERROR: CsA water ensemble not found. Expected ensemble.xyz + ensemble.json in "
             "data/CREST_CsA_20260512/ or results/runs/*_CsA/water/")

XYZ_PATH, JSON_PATH = _find_csa_water_ensemble()
print(f"Using ensemble: {XYZ_PATH}")

# ── Parse multi-frame XYZ ─────────────────────────────────────────────────────
def parse_xyz_ensemble(path: Path) -> list:
    conformers = []
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        try:
            n_atoms = int(line)
        except ValueError:
            i += 1
            continue
        i += 1
        comment = lines[i].strip() if i < len(lines) else ""
        i += 1
        energy = None
        for tok in comment.replace("=", " ").replace(":", " ").split():
            try:
                energy = float(tok)
                break
            except ValueError:
                continue
        symbols, coords = [], []
        for _ in range(n_atoms):
            if i >= len(lines):
                break
            parts = lines[i].split()
            if len(parts) >= 4:
                symbols.append(parts[0])
                coords.append([float(p) for p in parts[1:4]])
            i += 1
        if len(symbols) == n_atoms:
            conformers.append((symbols, np.array(coords), energy))
    return conformers


# ── xTB reference ────────────────────────────────────────────────────────────
print("Loading xTB ensemble ...", flush=True)
conformers = parse_xyz_ensemble(XYZ_PATH)
print(f"  {len(conformers)} conformers loaded from {XYZ_PATH}")

xtb_energies = np.array([e for _, _, e in conformers], dtype=float)
KCAL = 627.509
RT   = 1.987e-3 * 298.15

xtb_kcal = xtb_energies * KCAL
xtb_rel  = xtb_kcal - xtb_kcal.min()
xtb_w    = np.exp(-xtb_rel / RT)
xtb_w   /= xtb_w.sum()

# Cross-check against stored JSON weights
with open(JSON_PATH) as f:
    jdata = json.load(f)
json_w = np.array([c["boltzmannweight"] for c in jdata["conformers"]])
json_mismatch = np.abs(xtb_w - json_w).max()
if json_mismatch > 0.01:
    print(f"  WARNING: max deviation from stored JSON weights = {json_mismatch:.3f}")
else:
    print(f"  xTB weights match stored JSON (max deviation {json_mismatch:.4f}) OK")
print(f"  xTB dominant conformer: C1 at {xtb_w[0]*100:.1f}%")
print(f"  xTB energy range: {xtb_rel.max():.2f} kcal/mol over {len(conformers)} conformers")

# ── MACE-OFF single-point energies ───────────────────────────────────────────
try:
    import torch
    from mace.calculators import mace_off
    from ase import Atoms
except ImportError:
    print(
        "\n[ERROR] mace-torch not installed.\n"
        "  Run:  pip install mace-torch ase\n"
        "  Then re-run this script."
    )
    sys.exit(1)

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"\nRunning MACE-OFF23(M) on {device} ...")
print("  (first run downloads ~200 MB model weights — cached after that)\n")

calc = mace_off(model=args.model, device=device, default_dtype="float64")

EV_TO_KCAL = 23.0609
mace_energies = []
for idx, (symbols, coords, _) in enumerate(conformers):
    atoms = Atoms(symbols=symbols, positions=coords)
    atoms.calc = calc
    e_ev = atoms.get_potential_energy()
    mace_energies.append(e_ev)
    e_kcal = e_ev * EV_TO_KCAL
    print(f"  Conf {idx+1:2d}/{len(conformers)}: {e_ev:.4f} eV  ({e_kcal:.2f} kcal/mol)",
          flush=True)

mace_energies = np.array(mace_energies)
mace_kcal = mace_energies * EV_TO_KCAL
mace_rel  = mace_kcal - mace_kcal.min()
mace_w    = np.exp(-mace_rel / RT)
mace_w   /= mace_w.sum()

# ── Comparison stats ─────────────────────────────────────────────────────────
from scipy.stats import pearsonr, spearmanr

r_energy,   _ = pearsonr(xtb_rel, mace_rel)
rho_weight, _ = spearmanr(xtb_w, mace_w)
r_weight,   _ = pearsonr(xtb_w, mace_w)

xtb_ranks  = np.argsort(np.argsort(-xtb_w))  + 1
mace_ranks = np.argsort(np.argsort(-mace_w)) + 1
rank_shift = np.abs(mace_ranks - xtb_ranks)

print("\n" + "=" * 65)
print("  MACE-OFF23(M) vs GFN2-xTB+ALPB  —  CsA water ensemble")
print("=" * 65)
print(f"  Pearson r  (relative energies):      {r_energy:+.3f}")
print(f"  Pearson r  (Boltzmann weights):      {r_weight:+.3f}")
print(f"  Spearman r (Boltzmann weights):      {rho_weight:+.3f}")
print(f"  Mean absolute rank shift:            {rank_shift.mean():.1f} positions")
print(f"  Max rank shift:                      {rank_shift.max():.0f} positions"
      f"  (C{np.argmax(rank_shift)+1})")
print()

print(f"  {'Conf':>5}  {'xTB w%':>7}  {'MACE w%':>7}  "
      f"{'xTB dE':>8}  {'MACE dE':>8}  {'dRank':>6}")
print(f"  {'':->5}  {'':->7}  {'':->7}  {'':->8}  {'':->8}  {'':->6}")
for i in np.argsort(-xtb_w):
    dr = int(mace_ranks[i]) - int(xtb_ranks[i])
    flag = "  <-- top MACE" if mace_ranks[i] == 1 else ""
    print(f"  C{i+1:>4d}  {xtb_w[i]*100:>6.1f}%  {mace_w[i]*100:>6.1f}%  "
          f"{xtb_rel[i]:>7.2f}  {mace_rel[i]:>7.2f}  {dr:>+6d}{flag}")

# ── Save CSV ──────────────────────────────────────────────────────────────────
import csv
out_dir = Path("results")
out_dir.mkdir(exist_ok=True)
csv_path = out_dir / "mace_vs_xtb_CsA_water.csv"
with open(csv_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["conformer", "xtb_energy_Ha", "xtb_rel_kcal", "xtb_weight",
                "mace_energy_eV", "mace_rel_kcal", "mace_weight",
                "xtb_rank", "mace_rank", "rank_shift"])
    for i in range(len(conformers)):
        w.writerow([
            f"C{i+1}",
            f"{xtb_energies[i]:.8f}", f"{xtb_rel[i]:.4f}", f"{xtb_w[i]:.6f}",
            f"{mace_energies[i]:.6f}", f"{mace_rel[i]:.4f}", f"{mace_w[i]:.6f}",
            int(xtb_ranks[i]), int(mace_ranks[i]), int(mace_ranks[i]-xtb_ranks[i]),
        ])
print(f"\n  CSV saved -> {csv_path}")

# ── Scatter plots ─────────────────────────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    fig.suptitle("MACE-OFF23(M) vs GFN2-xTB+ALPB  |  CsA 23 water conformers",
                 fontsize=11)

    # Panel 1: relative energies
    ax = axes[0]
    sc = ax.scatter(xtb_rel, mace_rel, c=xtb_w, cmap="viridis", s=55, zorder=3,
                    edgecolors="k", linewidths=0.4)
    lim = max(xtb_rel.max(), mace_rel.max()) * 1.05 + 0.2
    ax.plot([0, lim], [0, lim], "k--", alpha=0.35, lw=1, label="y=x")
    ax.set_xlabel("xTB relative energy (kcal/mol)")
    ax.set_ylabel("MACE-OFF relative energy (kcal/mol)")
    ax.set_title(f"Relative energies  |  r = {r_energy:.3f}")
    ax.set_xlim(-0.2, lim)
    ax.set_ylim(-0.2, lim)
    plt.colorbar(sc, ax=ax, label="xTB Boltzmann weight")
    for i, (x, y) in enumerate(zip(xtb_rel, mace_rel)):
        if xtb_w[i] > 0.06 or mace_w[i] > 0.06:
            ax.annotate(f"C{i+1}", (x, y), fontsize=7, ha="left",
                        xytext=(3, 3), textcoords="offset points")

    # Panel 2: Boltzmann weights
    ax = axes[1]
    ax.scatter(xtb_w, mace_w, c=xtb_w, cmap="viridis", s=55, zorder=3,
               edgecolors="k", linewidths=0.4)
    wlim = max(xtb_w.max(), mace_w.max()) * 1.1
    ax.plot([0, wlim], [0, wlim], "k--", alpha=0.35, lw=1)
    ax.set_xlabel("xTB Boltzmann weight")
    ax.set_ylabel("MACE-OFF Boltzmann weight")
    ax.set_title(f"Boltzmann weights  |  r = {r_weight:.3f}  rho = {rho_weight:.3f}")
    for i, (xw, mw) in enumerate(zip(xtb_w, mace_w)):
        if xw > 0.04 or mw > 0.04:
            ax.annotate(f"C{i+1}", (xw, mw), fontsize=7, ha="left",
                        xytext=(3, 3), textcoords="offset points")

    plt.tight_layout()
    png_path = out_dir / "mace_vs_xtb_CsA_water.png"
    plt.savefig(png_path, dpi=150)
    print(f"  Plot saved -> {png_path}")
except ImportError:
    print("  (matplotlib not available — install with: pip install matplotlib)")

# ── Interpretation ───────────────────────────────────────────────────────────
print("\n" + "-" * 65)
print("  Interpretation guide")
print("-" * 65)
print(f"  r (energies) = {r_energy:.3f}")
if r_energy > 0.90:
    verdict = "GOOD — MACE-OFF and xTB agree on conformer ordering."
    next_step = "Proceed with MACE-OFF as energy engine in the GPU pipeline."
elif r_energy > 0.75:
    verdict = "MODERATE — ranking mostly preserved; some conformers swap."
    next_step = ("Use MACE-OFF for screening but validate key conformers (A, C)\n"
                 "  against xTB or DFT before reporting final descriptors.")
else:
    verdict = "POOR — likely driven by missing solvation in MACE-OFF (gas phase)."
    next_step = ("MACE-OFF vacuum energies disagree with ALPB-solvated xTB.\n"
                 "  Options: (a) add OpenMM GBn2 implicit solvent on top of MACE,\n"
                 "           (b) use MACE for geometry quality only, xTB for ranking,\n"
                 "           (c) proceed anyway if rank shift < 3 for top conformers.")
print(f"\n  {verdict}")
print(f"\n  Next step: {next_step}")
print()
print("  Key caveat: MACE-OFF23 is gas-phase; xTB energies include ALPB(water).")
print("  Some systematic offset is expected. What matters most is whether the")
print("  RELATIVE ordering of A-like vs C-like conformers is preserved.")
print()

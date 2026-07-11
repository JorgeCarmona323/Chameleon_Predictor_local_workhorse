# env: rdkit_env   (RDKit prep + matplotlib compositing; PyMOL renders the .pml)
"""
make_isomer_3d_views.py
-----------------------
Publication 3D figures for the R/S isomer pairs. RDKit does the ensemble bookkeeping
(top-N weighted conformers, minimum-energy pick, common-frame Kabsch alignment) and writes
PDBs; PyMOL ray-traces them; matplotlib composites the two-panel (water | chloroform) figures.

Per pair, three two-panel figures (A = water, B = chloroform):
  Figure_4_<pair>.png   S isomer conformational ensemble (top-N overlaid)
  Figure_5_<pair>.png   R isomer conformational ensemble
  Figure_6_<pair>.png   R vs S minimum-energy conformers superimposed (exposes the stereocenter)

Everything is Kabsch-aligned to ONE reference (the R-water minimum-energy conformer, heavy
atoms), so all panels share a camera and the R/S overlay is a true superposition.

Run (rdkit_env):  python scripts/make_isomer_3d_views.py --pair 3-12-10-12
Outputs -> results/figures/isomers/Figure_{4,5,6}_<pair>.png  (panels + PDBs in .../isomers/3d/)
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Geometry import Point3D

RDLogger.DisableLog("rdApp.*")

BASE = Path("results/conformers")
FIG_ROOT = Path("results/figures/isomers")
# OUT_FIG (per-molecule final figures) and OUT (per-molecule 'work/' intermediates) are set
# per --pair in main(); see _set_dirs().
OUT_FIG = FIG_ROOT
OUT = FIG_ROOT / "work"


def _set_dirs(pair):
    global OUT_FIG, OUT
    OUT_FIG = FIG_ROOT / pair
    OUT = OUT_FIG / "work"
    OUT.mkdir(parents=True, exist_ok=True)

# nested Xylene-Linker (non-diazirine) ensemble folders, relative to BASE
PAIRS = {
    "3-12-8-12": {
        "R": "DOPC 3-12-8-12/3-12-8-12 Xylene Linker/DOPC 3-12-8-12 R",
        "S": "DOPC 3-12-8-12/3-12-8-12 Xylene Linker/DOPC 3-12-8-12 S",
    },
    "3-12-10-12": {
        "R": "WhC3/3-12-10-12 Xylene Linker/3-12-10-12 R",
        "S": "WhC3/3-12-10-12 Xylene Linker/3-12-10-12 S",
    },
}
SOLVENTS = ("water", "mem")               # data-folder keys on disk (kept as-is)
SOLV_LABEL = {"water": "water", "mem": "chloroform"}
SOLV_TOKEN = {"water": "water", "mem": "chcl3"}   # token used in OUTPUT filenames (no 'mem')
N_OVERLAY = 20          # thin, semi-transparent fan reads well across the full top-20 spread
RAY = 2400              # ray-trace panel size (px); high-res for crisp print/SVG embedding
R_COLOR, S_COLOR = "orange", "teal"
# Launch PyMOL as `python -m pymol` (the canonical headless entry) — this auto-loads the
# academic license at ~/.pymol/license.lic. The `Scripts\pymol.exe -cq` launcher does NOT,
# and falls back to watermarked evaluation mode.
PYMOL_PYTHON_CANDIDATES = [
    r"C:\ProgramData\pymol\python.exe",
]
LICENSE_CANDIDATES = [
    os.path.expanduser(r"~/.pymol/license.lic"),
    r"C:\Users\Admin\.pymol\license.lic",
    r"C:\Users\Admin\Downloads\pymol-edu-license.lic",
]


# ── ensemble IO ───────────────────────────────────────────────────────────────
def load(wd: Path):
    """Return (mols, weights, energies) aligned by index; energies in Hartree."""
    confs = json.load(open(wd / "ensemble.json"))["conformers"]
    e = np.array([c.get("totalenergy", np.nan) for c in confs], float)
    w = np.array([c.get("boltzmannweight", np.nan) for c in confs], float)
    mols = [m for m in Chem.SDMolSupplier(str(wd / "ensemble.sdf"), removeHs=False, sanitize=True) if m]
    n = min(len(mols), len(w), len(e))
    return mols[:n], w[:n], e[:n]


def heavy_idx(mol):
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]


def coords(mol):
    c = mol.GetConformer()
    return np.array([list(c.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())], float)


# ── Kabsch (align P onto Q on a heavy-atom subset; apply to all atoms) ─────────
def kabsch_apply(full, P_heavy, Q_heavy):
    Pc, Qc = P_heavy.mean(0), Q_heavy.mean(0)
    H = (P_heavy - Pc).T @ (Q_heavy - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return (full - Pc) @ R.T + Qc


def set_coords(mol, xyz):
    c = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        c.SetAtomPosition(i, Point3D(*map(float, xyz[i])))
    return mol


def write_pdb_multi(template, frames_xyz, path):
    base = Chem.Mol(template)
    base.RemoveAllConformers()
    for xyz in frames_xyz:
        conf = Chem.Conformer(base.GetNumAtoms())
        for i in range(base.GetNumAtoms()):
            conf.SetAtomPosition(i, Point3D(*map(float, xyz[i])))
        base.AddConformer(conf, assignId=True)
    w = Chem.PDBWriter(str(path))
    for cid in range(base.GetNumConformers()):
        w.write(base, confId=cid)
    w.close()


# ── prep: aligned PDBs for one pair ───────────────────────────────────────────
def prep_pair(pair):
    OUT.mkdir(parents=True, exist_ok=True)
    iso_d = PAIRS[pair]
    hidx = None
    ref = None  # global reference heavy-atom coords (R-water minimum)

    # establish reference from R-water minimum-energy conformer
    rmols, rw, re_ = load(BASE / iso_d["R"] / "water")
    hidx = heavy_idx(rmols[0])
    rmin = int(np.nanargmin(re_))
    ref = coords(rmols[rmin])[hidx]

    for iso, dname in iso_d.items():
        for solv in SOLVENTS:
            mols, wts, en = load(BASE / dname / solv)      # reads <dir>/mem/ on disk
            tag = f"{pair}_{iso}_{SOLV_TOKEN[solv]}"          # writes ..._chcl3_*.pdb

            imin = int(np.nanargmin(en))
            full = coords(mols[imin])
            # Fig 6 min: aligned to the GLOBAL reference (R-water min) so R and S superimpose
            write_pdb_multi(mols[imin], [kabsch_apply(full, full[hidx], ref)],
                            OUT / f"{tag}_min.pdb")

            # Fig 4/5 fan: per-solvent intra-fit — align to THIS ensemble's own lowest-energy
            # conformer, so each panel shows the solvent's true internal spread rather than a
            # cross-solvent (water-frame) alignment bias
            local_ref = full[hidx]
            top = [int(i) for i in np.argsort(wts)[::-1][:N_OVERLAY]]
            frames = [kabsch_apply(coords(mols[j]), coords(mols[j])[hidx], local_ref) for j in top]
            write_pdb_multi(mols[top[0]], frames, OUT / f"{tag}_overlay.pdb")

    print(f"  prepped aligned PDBs for {pair} -> {OUT}")


# ── PyMOL .pml ────────────────────────────────────────────────────────────────
HEAD = (
    "set multiplex, 0\n"               # load multi-MODEL PDB as states, not separate objects
    "bg_color white\n"
    "set ray_opaque_background, 0\n"
    "set orthoscopic, 1\n"
    "set ray_shadows, 0\n"
    "set antialias, 2\n"
    "set ambient_occlusion_mode, 1\n"
    "set stick_h_scale, 1\n"
)


def _fan_panel(pair, iso, solv, color):
    # load BOTH solvents and orient on the union so this figure's two panels share one camera,
    # then show only the target object. One panel per PyMOL process (avoids the session-
    # accumulation watermark trip in this build).
    # per-solvent intra-fit fans live in independent frames, so load + orient on this panel's
    # own object. Thin, semi-transparent sticks (no outline) so the conformational fan reads
    # as an overlaid bundle rather than an opaque blob.
    tok = SOLV_TOKEN[solv]
    return HEAD + (
        f"load {pair}_{iso}_{tok}_overlay.pdb, fan\n"
        "create mec, fan, 1, 1\n"          # state 1 = lowest-energy (MEC); split it off
        "remove hydrogens\nset ray_trace_mode, 0\n"
        # the rest of the ensemble: thin, semi-transparent fan
        "set all_states, on, fan\n"
        "set stick_radius, 0.07, fan\nset stick_transparency, 0.55, fan\n"
        # the minimum-energy conformer: thicker and opaque so it reads clearly on top
        "set stick_radius, 0.13, mec\nset stick_transparency, 0, mec\n"
        "orient fan\n"
        f"hide everything\nshow sticks, fan\nshow sticks, mec\n"
        f"color {color}, fan\ncolor {color}, mec\n"
        f"ray {RAY}, {RAY}\npng panel_{pair}_{iso}_{tok}.png, dpi=300\n"
    )


def _overlay_panel(pair, solv):
    a = "w" if solv == "water" else "m"
    R, S = f"R{a}", f"S{a}"
    tok = SOLV_TOKEN[solv]
    return HEAD + (
        f"load {pair}_R_water_min.pdb, Rw\nload {pair}_R_chcl3_min.pdb, Rm\n"
        f"load {pair}_S_water_min.pdb, Sw\nload {pair}_S_chcl3_min.pdb, Sm\n"
        "remove hydrogens\nset stick_transparency, 0\nset stick_radius, 0.16\n"
        "set ray_trace_mode, 1\nset ray_trace_color, grey30\norient\n"
        f"hide everything\nshow sticks, {R}\nshow sticks, {S}\n"
        f"color {R_COLOR}, {R}\ncolor {S_COLOR}, {S}\n"
        f"ray {RAY}, {RAY}\npng panel_{pair}_RvsS_{tok}.png, dpi=300\n"
    )


def write_pmls(pair):
    """One self-contained .pml per panel (rendered in its own process)."""
    jobs = []
    for solv in SOLVENTS:
        jobs.append((f"S_{SOLV_TOKEN[solv]}", _fan_panel(pair, "S", solv, S_COLOR)))      # Figure 4
    for solv in SOLVENTS:
        jobs.append((f"R_{SOLV_TOKEN[solv]}", _fan_panel(pair, "R", solv, R_COLOR)))      # Figure 5
    for solv in SOLVENTS:
        jobs.append((f"RvsS_{SOLV_TOKEN[solv]}", _overlay_panel(pair, solv)))             # Figure 6
    paths = []
    for name, body in jobs:
        p = OUT / f"_panel_{pair}_{name}.pml"
        p.write_text(body, encoding="utf-8")
        paths.append(p)
    return paths


def run_pymols(pmls):
    py = next((p for p in PYMOL_PYTHON_CANDIDATES if Path(p).exists()), None)
    if not py:
        sys.exit("PyMOL python not found; checked: " + ", ".join(PYMOL_PYTHON_CANDIDATES))
    # IMPORTANT: do NOT set PYMOL_LICENSE_FILE — PyMOL 3.x auto-detects the new-format
    # entitlement at ~/.pymol/license.lic; the legacy env var expects the OLD format and,
    # if set, overrides auto-detection and forces watermarked evaluation mode. Also render
    # ONE panel per process: many ray calls in a single session trip this build's watermark.
    lic = next((p for p in LICENSE_CANDIDATES if Path(p).exists()), None)
    print(f"  PyMOL license (auto-detected): {lic or 'NONE FOUND — will watermark'}")
    for pml in pmls:
        # cwd=OUT so the relative load/png paths in the .pml resolve there
        subprocess.run([py, "-m", "pymol", "-cq", pml.name], cwd=str(OUT), check=True)


# ── compositing ───────────────────────────────────────────────────────────────
def compose(pair):
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    def two_panel(panels, title, out):
        fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.7), dpi=300)
        for ax, (png, lab, sub) in zip(axes, panels):
            ax.imshow(mpimg.imread(str(OUT / png)))
            ax.axis("off")
            ax.set_title(sub, fontsize=10)
            ax.text(0.03, 0.97, lab, transform=ax.transAxes, fontsize=15,
                    fontweight="bold", va="top", ha="left")
        fig.suptitle(title, fontsize=11, fontweight="bold", y=0.99)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight")
        # SVG twin: layout/labels/titles are true vector; the ray-traced molecule panels
        # are embedded as full-res raster (3D ray-traces have no vector representation)
        fig.savefig(out.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {out}  (+ .svg)")

    # titles carry no "Figure N." prefix — the report caption supplies the number (avoids
    # double-numbering when embedded); the image keeps a descriptive standalone title
    figs = [
        ("fig5_S_ensemble.png", f"{pair} S isomer — conformational ensemble (top {N_OVERLAY})",
         [(f"panel_{pair}_S_water.png", "A", "water"), (f"panel_{pair}_S_chcl3.png", "B", "chloroform")]),
        ("fig6_R_ensemble.png", f"{pair} R isomer — conformational ensemble (top {N_OVERLAY})",
         [(f"panel_{pair}_R_water.png", "A", "water"), (f"panel_{pair}_R_chcl3.png", "B", "chloroform")]),
        ("fig7_RS_overlay.png", f"{pair} minimum-energy conformers — R (orange) vs S (teal)",
         [(f"panel_{pair}_RvsS_water.png", "A", "water"), (f"panel_{pair}_RvsS_chcl3.png", "B", "chloroform")]),
    ]
    for name, title, panels in figs:
        two_panel(panels, title, OUT_FIG / name)


def main():
    ap = argparse.ArgumentParser(description="Publication R/S 3D ensemble figures (PyMOL).")
    ap.add_argument("--pair", choices=list(PAIRS), required=True)
    ap.add_argument("--skip-render", action="store_true", help="reuse existing panel PNGs")
    args = ap.parse_args()

    _set_dirs(args.pair)
    prep_pair(args.pair)
    pmls = write_pmls(args.pair)
    if not args.skip_render:
        run_pymols(pmls)
    compose(args.pair)


if __name__ == "__main__":
    main()

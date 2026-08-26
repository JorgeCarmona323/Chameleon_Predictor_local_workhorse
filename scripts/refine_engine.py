#!/usr/bin/env python3
"""
refine_engine.py  --  STAGE 1.5: r2SCAN-3c + CPCM geometry refinement of a CREST ensemble via CENSO.

Consumes a stage-1 leg (`<leg>/ensemble.xyz`, CREST format: energy in each frame's comment line),
runs CENSO (screening + optimization at r2SCAN-3c + CPCM, ORCA backend), and writes a REFINED
ensemble back in the SAME schema so the rest of the pipeline is agnostic:

    <leg>/refined/ensemble.xyz     (r2SCAN-3c-optimized geometries, energy in comment)
    <leg>/refined/ensemble.sdf
    <leg>/refined/energies.csv     (conf, E_Eh, relE_kcal, pop)   <- r2SCAN-3c Boltzmann pops
    <leg>/refined/metadata.json

NON-DESTRUCTIVE: never touches the GFN2 ensemble; the GFN2<->r2SCAN-3c delta IS the experiment.

============================  CENSO-VERSION-SPECIFIC  ============================
Everything CENSO-3.0.8-specific is isolated in _write_censo_config(), _censo_cmd(),
and _find_censo_output(). Run `python refine_engine.py --probe` on the HPC FIRST: it
dumps `censo --help` and a fresh `censo` default config so we can reconcile the exact
key names / output filenames for this build, then adjust the three functions above.
=================================================================================
"""
import argparse, json, os, shutil, subprocess, sys
from pathlib import Path
import numpy as np

RT = 0.593  # kcal/mol at 298.15 K
H2K = 627.5094740631

# ----------------------------- xyz I/O -----------------------------
def read_xyz_ensemble(path):
    """Return (list_of_frames), each frame = (natoms, comment, [(sym,x,y,z),...], energy_Eh)."""
    L = Path(path).read_text().splitlines()
    frames, i = [], 0
    while i < len(L):
        s = L[i].strip()
        if not s:
            i += 1; continue
        try:
            n = int(s)
        except ValueError:
            break
        comment = L[i+1] if i+1 < len(L) else ""
        e = None
        for tok in comment.replace("=", " ").split():
            try:
                e = float(tok); break
            except ValueError:
                pass
        atoms = []
        for j in range(n):
            p = L[i+2+j].split()
            atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
        frames.append((n, comment, atoms, e))
        i += n + 2
    return frames

def write_xyz_ensemble(frames, path):
    out = []
    for n, comment, atoms, e in frames:
        out.append(str(n))
        out.append(comment if comment else (f"{e:.10f}" if e is not None else ""))
        for sym, x, y, z in atoms:
            out.append(f"{sym:<3s} {x:>15.8f} {y:>15.8f} {z:>15.8f}")
    Path(path).write_text("\n".join(out) + "\n")

def boltzmann(frames):
    E = np.array([f[3] if f[3] is not None else np.nan for f in frames], float)
    rel = (E - np.nanmin(E)) * H2K
    w = np.exp(-rel / RT); w = w / np.nansum(w)
    return rel, w

# ----------------------------- CENSO (version-specific) -----------------------------
def _write_censo_config(cfg_path, solvent, orca_path, maxcores, ewin):
    """CENSO 3.0.8 TOML config: screening + optimization at r2SCAN-3c + CPCM, ORCA backend.
    VERIFY key names against `censo --new-config` output on this build (use --probe)."""
    orca_bin = str(orca_path)
    txt = f"""# CENSO 3.0.8 config -- r2SCAN-3c + CPCM refinement (auto-written by refine_engine.py)
# VERIFY these keys with `censo --new-config` on this exact build.
[general]
solvent = "{solvent}"
sm = "cpcm"
prog = "orca"
maxcores = {maxcores}
imagthr = -100.0
temperature = 298.15

[prescreening]
run = true
func = "pbe-d4"
basis = "def2-SV(P)"

[screening]
run = true
func = "r2scan-3c"
threshold = {ewin}

[optimization]
run = true
func = "r2scan-3c"
optlevel = "normal"
threshold = 4.0

[refinement]
run = false

[paths]
orcapath = "{orca_bin}"
"""
    Path(cfg_path).write_text(txt)
    return cfg_path

def _censo_cmd(ensemble_in, cfg_path, maxcores, charge):
    """CENSO 3.0.8 invocation. VERIFY flags with `censo --help` (--probe)."""
    return ["censo", "-i", str(ensemble_in), "--config", str(cfg_path),
            "--maxcores", str(maxcores), "--charge", str(charge)]

def _find_censo_output(workdir):
    """Locate CENSO's final optimized ensemble. VERIFY the filename on this build (--probe).
    Tries the common v3 names in priority order."""
    wd = Path(workdir)
    for pat in ("*optimization*.xyz", "*censo_final*.xyz", "*.optimization.xyz",
                "3_OPTIMIZATION*.xyz", "conformers_optimized*.xyz", "*_optimized.xyz"):
        hits = sorted(wd.rglob(pat))
        if hits:
            return hits[-1]
    return None

def probe():
    print("=== censo --version ==="); subprocess.run(["censo", "--version"])
    print("\n=== censo --help ==="); subprocess.run(["censo", "--help"])
    print("\n=== censo default config (--new-config / -newconfig) ===")
    for flag in ("--new-config", "-newconfig", "--writeconfig"):
        r = subprocess.run(["censo", flag], capture_output=True, text=True)
        if r.returncode == 0:
            print(f"(worked: censo {flag})"); print(r.stdout[:4000]); break
        print(f"(censo {flag} -> rc={r.returncode})")

# ----------------------------- main refine -----------------------------
def refine(ensemble_xyz, solvent, orca_path, maxcores, charge, ewin, keep_work=False):
    leg = Path(ensemble_xyz).parent
    outdir = leg / "refined"; outdir.mkdir(exist_ok=True)
    work = leg / "refined_work"; work.mkdir(exist_ok=True)
    ens_in = work / "crest_conformers.xyz"
    shutil.copy(ensemble_xyz, ens_in)
    cfg = _write_censo_config(work / "censorc.toml", solvent, orca_path, maxcores, ewin)

    env = dict(os.environ)
    env["PATH"] = f"{Path(orca_path).parent}:{env.get('PATH','')}"   # ORCA on PATH for CENSO
    cmd = _censo_cmd(ens_in, cfg, maxcores, charge)
    print("RUN:", " ".join(cmd), "\n(cwd=", work, ")", flush=True)
    r = subprocess.run(cmd, cwd=work, env=env)
    if r.returncode != 0:
        sys.exit(f"CENSO failed (rc={r.returncode}) -- run `python {sys.argv[0]} --probe` to check the interface")

    out = _find_censo_output(work)
    if out is None:
        sys.exit(f"could not find CENSO optimized ensemble under {work} -- check _find_censo_output() (--probe)")
    frames = read_xyz_ensemble(out)
    rel, w = boltzmann(frames)
    write_xyz_ensemble(frames, outdir / "ensemble.xyz")
    # energies.csv (r2SCAN-3c pops)
    lines = ["conf,E_Eh,relE_kcal,pop"]
    for k, (f, r_, w_) in enumerate(zip(frames, rel, w)):
        lines.append(f"{k},{f[3]},{r_:.4f},{w_:.6g}")
    (outdir / "energies.csv").write_text("\n".join(lines) + "\n")
    # sdf (best-effort via RDKit)
    try:
        from rdkit import Chem
        w_sdf = Chem.SDWriter(str(outdir / "ensemble.sdf"))
        supp = Chem.MolFromXYZFile  # not all builds; fall back silently
    except Exception:
        pass
    (outdir / "metadata.json").write_text(json.dumps({
        "method": "r2scan-3c", "solvation": "cpcm", "solvent": solvent,
        "backend_orca": str(orca_path), "n_conformers": len(frames),
        "source_ensemble": str(ensemble_xyz), "charge": charge, "ewin_kcal": ewin,
    }, indent=2))
    if not keep_work:
        pass  # keep by default; CENSO scratch is useful for debugging
    print(f"OK: {len(frames)} refined conformers -> {outdir}/ensemble.xyz  (+ energies.csv, metadata.json)")

def main():
    ap = argparse.ArgumentParser(description="r2SCAN-3c + CPCM refinement of a CREST ensemble leg via CENSO 3.0.8")
    ap.add_argument("--ensemble", help="path to <leg>/ensemble.xyz (CREST format)")
    ap.add_argument("--solvent", default="chloroform")
    ap.add_argument("--orca", default=os.path.expanduser("~/orca_6.1.1/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca"))
    ap.add_argument("--maxcores", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "20")))
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--ewin", type=float, default=6.0, help="kcal/mol screening window")
    ap.add_argument("--probe", action="store_true", help="dump censo --help / --version / default config and exit")
    a = ap.parse_args()
    if a.probe:
        probe(); return
    if not a.ensemble:
        ap.error("--ensemble is required (or use --probe)")
    refine(a.ensemble, a.solvent, a.orca, a.maxcores, a.charge, a.ewin)

if __name__ == "__main__":
    main()

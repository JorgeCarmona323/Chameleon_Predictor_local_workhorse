#!/usr/bin/env python3
"""
refine_engine.py  --  STAGE 1.5: r2SCAN-3c + CPCM refinement of a CREST/GOAT ensemble via CENSO 3.0.8.

Consumes a stage-1 leg (`<leg>/ensemble.xyz`, CREST/GOAT format: energy in each frame's comment),
runs CENSO screening + optimization at r2SCAN-3c + CPCM (ORCA backend), and writes a REFINED
ensemble back in the SAME schema so the rest of the pipeline is agnostic:

    <leg>/refined/ensemble.xyz     (r2SCAN-3c-optimized geometries)
    <leg>/refined/energies.csv     (conf, E_Eh, relE_kcal, pop)
    <leg>/refined/2_OPTIMIZATION.{xyz,json,out}   (raw CENSO outputs, for exact re-weighting)
    <leg>/refined/metadata.json

NON-DESTRUCTIVE: never touches the input ensemble.

Verified against the CENSO 3.0.8 source (reference repos/CENSO):
  - config is INI (~/.censo2rc style); parts enabled by CLI flags (-P/-S/-O/-R)
  - prog=orca, sm=cpcm (OrcaSolvMod), func=r2scan-3c, basis=def2-mtzvpp
  - [paths] orca=<real orca> fixes the /usr/bin/orca screen-reader auto-detect
  - solvent identifier is CENSO's key (chcl3, h2o, dmso, hexane); CENSO maps it per backend
  - optimized ensemble is written as 2_OPTIMIZATION.xyz
"""
import argparse, json, os, shutil, subprocess, sys
from pathlib import Path
# numpy imported lazily in boltzmann() so --probe works in a minimal env

RT = 0.593           # kcal/mol at 298.15 K
H2K = 627.5094740631

# our solvent name -> CENSO solvent identifier (key in censo assets/solvents.json)
_CENSO_SOLV = {"chloroform": "chcl3", "chcl3": "chcl3", "water": "h2o", "h2o": "h2o",
               "dmso": "dmso", "hexane": "hexane"}

# ----------------------------- xyz I/O -----------------------------
def read_xyz_ensemble(path):
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
        atoms = [(p[0], float(p[1]), float(p[2]), float(p[3]))
                 for p in (L[i+2+j].split() for j in range(n))]
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
    import numpy as np
    E = np.array([f[3] if f[3] is not None else np.nan for f in frames], float)
    if np.all(np.isnan(E)):
        w = np.full(len(frames), 1.0 / max(len(frames), 1)); return np.zeros(len(frames)), w
    rel = (E - np.nanmin(E)) * H2K
    w = np.exp(-rel / RT); w = w / np.nansum(w)
    return rel, w

# ----------------------------- CENSO (3.0.8) -----------------------------
def _write_censo_config(cfg_path, solvent_id, orca_path, xtb_path, ewin, evaluate_rrho):
    """INI config for CENSO 3.0.8: screening + optimization at r2SCAN-3c + CPCM, ORCA backend."""
    txt = f"""[general]
temperature = 298.15
evaluate_rrho = {evaluate_rrho}
sm_rrho = alpb
imagthr = -100.0
sthr = 50.0
solvent = {solvent_id}
gas_phase = False
copy_mo = True
balance = True
ignore_failed = True

[prescreening]
prog = orca
func = pbe-d3
basis = def2-sv(p)
gfnv = gfn2
threshold = 4.0
template = False

[screening]
prog = orca
func = r2scan-3c
basis = def2-mtzvpp
sm = cpcm
gfnv = gfn2
threshold = 3.5
gsolv_included = False
template = False

[optimization]
prog = orca
func = r2scan-3c
basis = def2-mtzvpp
sm = cpcm
gfnv = gfn2
optcycles = 8
maxcyc = 200
optlevel = normal
threshold = {ewin}
gradthr = 0.01
hlow = 0.01
macrocycles = True
constrain = False
xtb_opt = True
template = False

[paths]
orca = {orca_path}
xtb = {xtb_path}
tm =
cosmotherm =
cosmorssetup =
"""
    Path(cfg_path).write_text(txt)
    return cfg_path

def probe():
    print("=== censo --version ==="); subprocess.run(["censo", "--version"])
    print("\n=== censo --help ==="); subprocess.run(["censo", "--help"])
    print("\n=== censo --new-config (needs real ORCA on PATH) ===")
    r = subprocess.run(["censo", "--new-config"], capture_output=True, text=True)
    print(f"(rc={r.returncode})"); print((r.stdout or "")[:2000]); print((r.stderr or "")[:2000])

# ----------------------------- main refine -----------------------------
def refine(ensemble_xyz, solvent, orca_path, xtb_path, maxcores, charge, ewin, nconf, evaluate_rrho):
    solvent_id = _CENSO_SOLV.get(solvent.lower(), solvent.lower())
    leg = Path(ensemble_xyz).parent
    outdir = leg / "refined"; outdir.mkdir(exist_ok=True)
    work = leg / "refined_work"; work.mkdir(exist_ok=True)
    ens_in = work / "crest_conformers.xyz"
    shutil.copy(ensemble_xyz, ens_in)
    cfg = _write_censo_config(work / "censo2rc", solvent_id, orca_path, xtb_path, ewin, evaluate_rrho)

    cmd = ["censo", "-i", ens_in.name, "--inprc", "censo2rc", "-S", "-O",
           "-c", str(charge), "--maxcores", str(maxcores)]
    if nconf and nconf > 0:
        cmd += ["-n", str(nconf)]
    env = dict(os.environ)
    env["PATH"] = f"{Path(orca_path).parent}:{env.get('PATH','')}"   # real ORCA ahead of /usr/bin/orca
    print("RUN:", " ".join(cmd), "  (cwd=", work, ")", flush=True)
    r = subprocess.run(cmd, cwd=work, env=env)
    if r.returncode != 0:
        sys.exit(f"CENSO failed (rc={r.returncode}). Check {work}/censo.log")

    opt_xyz = work / "2_OPTIMIZATION.xyz"
    if not opt_xyz.is_file():
        sys.exit(f"no 2_OPTIMIZATION.xyz under {work} -- CENSO optimization did not complete")
    frames = read_xyz_ensemble(opt_xyz)
    rel, w = boltzmann(frames)
    write_xyz_ensemble(frames, outdir / "ensemble.xyz")
    lines = ["conf,E_Eh,relE_kcal,pop"]
    for k, (f, r_, w_) in enumerate(zip(frames, rel, w)):
        lines.append(f"{k},{f[3]},{r_:.4f},{w_:.6g}")
    (outdir / "energies.csv").write_text("\n".join(lines) + "\n")
    for extra in ("2_OPTIMIZATION.xyz", "2_OPTIMIZATION.json", "2_OPTIMIZATION.out"):
        p = work / extra
        if p.is_file(): shutil.copy(p, outdir / extra)
    (outdir / "metadata.json").write_text(json.dumps({
        "method": "r2scan-3c", "solvation": "cpcm", "solvent": solvent_id,
        "backend_orca": str(orca_path), "n_conformers": len(frames),
        "source_ensemble": str(ensemble_xyz), "charge": charge, "ewin_kcal": ewin,
        "evaluate_rrho": evaluate_rrho, "nconf_cap": nconf,
    }, indent=2))
    print(f"OK: {len(frames)} refined conformers -> {outdir}/ensemble.xyz  (+ energies.csv, raw CENSO outputs)")

def main():
    ap = argparse.ArgumentParser(description="r2SCAN-3c + CPCM refinement of a CREST/GOAT ensemble via CENSO 3.0.8")
    ap.add_argument("--ensemble", help="path to <leg>/ensemble.xyz")
    ap.add_argument("--solvent", default="chloroform", help="chloroform/water/dmso/hexane (mapped to CENSO id)")
    ap.add_argument("--orca", default=os.path.expanduser("~/orca_6.1.1/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca"))
    ap.add_argument("--xtb", default=os.environ.get("XTBEXE", ""), help="xtb path (default $XTBEXE)")
    ap.add_argument("--maxcores", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "20")))
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--ewin", type=float, default=4.0, help="optimization ΔG window (kcal/mol)")
    ap.add_argument("--nconf", type=int, default=150, help="cap input to first N (lowest-GFN2) conformers; 0=all")
    ap.add_argument("--rrho", action="store_true", help="evaluate mRRHO thermal corrections (slower; off by default)")
    ap.add_argument("--probe", action="store_true", help="dump censo --version/--help/--new-config and exit")
    a = ap.parse_args()
    if a.probe:
        probe(); return
    if not a.ensemble:
        ap.error("--ensemble is required (or use --probe)")
    xtb = a.xtb or shutil.which("xtb") or ""
    if not xtb:
        sys.exit("no xtb path (set --xtb or $XTBEXE)")
    refine(a.ensemble, a.solvent, a.orca, xtb, a.maxcores, a.charge, a.ewin, a.nconf, str(bool(a.rrho)))

if __name__ == "__main__":
    main()

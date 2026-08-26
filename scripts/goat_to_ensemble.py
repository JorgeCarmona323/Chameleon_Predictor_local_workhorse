#!/usr/bin/env python3
"""Convert an ORCA GOAT final ensemble into our ensemble schema (ensemble.xyz + energies.csv +
metadata.json), so GOAT output drops straight into refine_engine.py / the descriptor scripts.

Reuses refine_engine's xyz I/O + Boltzmann weighting (DRY)."""
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from refine_engine import read_xyz_ensemble, write_xyz_ensemble, boltzmann

def find_final_ensemble(d):
    """ORCA GOAT writes <base>.finalensemble.xyz. VERIFY the name on this ORCA build."""
    d = Path(d)
    for pat in ("*.finalensemble.xyz", "*finalensemble*.xyz", "*.globalminimum.xyz", "*.ensemble.xyz"):
        hits = sorted(d.glob(pat))
        if hits:
            return hits[-1]
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--goat-dir", required=True, help="dir holding ORCA GOAT output (*.finalensemble.xyz)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--solvent", default="")
    ap.add_argument("--charge", type=int, default=0)
    a = ap.parse_args()
    fe = find_final_ensemble(a.goat_dir)
    if fe is None:
        sys.exit(f"no GOAT final ensemble under {a.goat_dir} -- check find_final_ensemble() vs this ORCA build")
    frames = read_xyz_ensemble(fe)
    rel, w = boltzmann(frames)
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    write_xyz_ensemble(frames, out / "ensemble.xyz")
    lines = ["conf,E_Eh,relE_kcal,pop"]
    for k, (f, r_, w_) in enumerate(zip(frames, rel, w)):
        lines.append(f"{k},{f[3]},{r_:.4f},{w_:.6g}")
    (out / "energies.csv").write_text("\n".join(lines) + "\n")
    (out / "metadata.json").write_text(json.dumps({
        "method": "goat-gfn2-xtb", "solvation": "alpb", "solvent": a.solvent,
        "n_conformers": len(frames), "source": str(fe), "charge": a.charge,
    }, indent=2))
    print(f"OK: {len(frames)} GOAT conformers -> {out}/ensemble.xyz  (+ energies.csv, metadata.json)")

if __name__ == "__main__":
    main()

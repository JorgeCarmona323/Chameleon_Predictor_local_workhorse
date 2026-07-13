# env: chameleon-calc
"""
build_ensemble_from_crest.py
----------------------------
Turn a raw CREST `crest_conformers.xyz` (coords only, no bonds) into the
`ensemble.sdf` + `ensemble.json` pair the analysis scripts expect
(validate_csa_water.py, ensemble_descriptors.py). Useful for a CREST run that
was stopped before the pipeline's post-processing wrote those files.

Bonds/topology come from a template SDF (a completed ensemble of the SAME
molecule, same atom order). Boltzmann weights are computed from the per-conformer
energies in the xyz comment lines (GFN2 Hartree) at 298.15 K.

Usage:
  python scripts/build_ensemble_from_crest.py \
      --xyz results/conformers/CSA_v2_water/crest_conformers.xyz \
      --template data/CREST_CsA_20260512/ensemble.sdf \
      --outdir results/conformers/CSA_v2_water
"""
import argparse
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Conformer

RT_KCAL = 0.5924096   # kcal/mol at 298.15 K
HARTREE_KCAL = 627.5094740631


def parse_xyz(path):
    lines = Path(path).read_text().splitlines()
    confs, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        comment = lines[i + 1] if i + 1 < len(lines) else ""
        try:
            e = float(comment.split()[0])
        except Exception:
            e = np.nan
        block = lines[i + 2: i + 2 + n]
        els = [b.split()[0] for b in block]
        xyz = np.array([[float(v) for v in b.split()[1:4]] for b in block])
        confs.append((els, xyz, e))
        i += 2 + n
    return confs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xyz", required=True)
    ap.add_argument("--template", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    confs = parse_xyz(args.xyz)
    print(f"parsed {len(confs)} conformers from {args.xyz}")

    tmpl = next(m for m in Chem.SDMolSupplier(args.template, removeHs=False, sanitize=True) if m)
    tel = [a.GetSymbol() for a in tmpl.GetAtoms()]
    cel = confs[0][0]
    if tel != cel:
        print(f"ATOM ORDER MISMATCH: template {len(tel)} vs xyz {len(cel)} atoms")
        for k, (a, b) in enumerate(zip(tel, cel)):
            if a != b:
                print(f"  first diff at atom {k}: template={a} xyz={b}")
                break
        raise SystemExit(1)
    print(f"atom order matches ({len(tel)} atoms)")

    E = np.array([c[2] for c in confs], float) * HARTREE_KCAL
    Erel = E - np.nanmin(E)
    w = np.exp(-Erel / RT_KCAL)
    w = w / np.nansum(w)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(outdir / "ensemble.sdf"))
    for els, xyz, e in confs:
        m = Chem.Mol(tmpl)
        m.RemoveAllConformers()
        conf = Conformer(m.GetNumAtoms())
        for k, (x, y, z) in enumerate(xyz):
            conf.SetAtomPosition(k, (float(x), float(y), float(z)))
        m.AddConformer(conf, assignId=True)
        writer.write(m)
    writer.close()

    entries = [{"totalenergy": float(c[2]), "boltzmannweight": float(wi)}
               for c, wi in zip(confs, w)]
    json.dump({"conformers": entries}, open(outdir / "ensemble.json", "w"))
    print(f"wrote {outdir/'ensemble.sdf'} + ensemble.json  "
          f"(min-E conf weight {w.max():.3f}, {int((w > 0.001).sum())} conf > 0.1%)")


if __name__ == "__main__":
    main()

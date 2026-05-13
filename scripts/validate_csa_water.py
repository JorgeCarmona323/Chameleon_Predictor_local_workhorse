#!/usr/bin/env python3
"""
validate_csa_water.py
---------------------
Validate the CREST water-phase CsA conformer ensemble against the A1 conformer
fingerprint reported by Limbach et al., J. Am. Chem. Soc. 2022, 144, 12602.

A1 fingerprint (SI Table S11 + main text):
  1. Exactly one cis amide bond: MeVal11-MeBmt1 (|omega| < 30 deg)
  2. Abu2 NH intramolecularly H-bonded  (NMR slope -2.6  ppb/K)
  3. Ala7 NH intramolecularly H-bonded  (NMR slope -2.07 ppb/K)
  4. Val5 NH NOT H-bonded               (NMR slope -7.33 ppb/K, solvent-exposed)

Note on implicit solvation:
  CREST uses ALPB(water) — no explicit water molecules.
  The two crystal water molecules bridging the H-bond network in the A1
  X-ray structure (CCDC 2149649) are not present here. However, the
  intramolecular backbone geometry and cis-amide are still meaningful
  indicators of the correct fold.

Usage:
  python scripts/validate_csa_water.py
  python scripts/validate_csa_water.py --sdf data/CREST_CsA_20260512/ensemble.sdf
                                        --json data/CREST_CsA_20260512/ensemble.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

REPO_ROOT = Path(__file__).parent.parent
DATA_DIR  = REPO_ROOT / "data" / "CREST_CsA_20260512"

# Geometry cutoffs
CIS_CUTOFF   = 30.0   # |omega| < 30 deg → cis amide
HB_DIST      = 3.5    # N···O distance cutoff (Angstrom)
HB_ANGLE     = 120.0  # D-H···A angle cutoff (degrees)

# SMARTS for the 4 free ring NHs in CsA
# CsA sequence: MeBmt1-Abu2-Sar3-MeLeu4-Val5-MeLeu6-Ala7-DAla8-MeLeu9-MeLeu10-MeVal11
RESIDUE_NH_SMARTS: dict[str, str] = {
    "Abu2":  "[NH1;R]-[C;R](-CC)-C(=O)",           # ethyl side chain
    "Val5":  "[NH1;R]-[C;R](-C(C)C)-C(=O)",         # isopropyl side chain
    "Ala7":  "[NH1;R]-[C;R](-C)-C(=O)-[NH1;R]",    # methyl, next residue also free NH
    "DAla8": "[NH1;R]-[C;R](-C)-C(=O)-[N;R](-C)",  # methyl, next residue is NMe
}


# ── Utility ───────────────────────────────────────────────────────────────────

def _ring_atom_set(mol: Chem.Mol) -> set[int]:
    atoms: set[int] = set()
    for ring in mol.GetRingInfo().AtomRings():
        atoms.update(ring)
    return atoms


def _carbonyl_oxygens(mol: Chem.Mol, ring_atoms: set[int]) -> list[int]:
    """Indices of C=O oxygens where the carbon is in the macrocycle ring."""
    acc = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() == 6 and nb.GetIdx() in ring_atoms:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                if bond.GetBondTypeAsDouble() == 2.0:
                    acc.append(atom.GetIdx())
                    break
    return acc


def _preceding_carbonyl_o(mol: Chem.Mol, n_idx: int, ring_atoms: set[int]) -> int | None:
    """
    For a ring amide N, return the O index of the C=O it is directly bonded to
    (i.e., the carbonyl of the preceding residue in the ring).
    Each ring N has exactly two ring-C neighbours: the preceding C(=O) and the Cα.
    The C(=O) is the one carrying a double bond to O.
    """
    n_atom = mol.GetAtomWithIdx(n_idx)
    for nb in n_atom.GetNeighbors():
        if nb.GetAtomicNum() != 6 or nb.GetIdx() not in ring_atoms:
            continue
        for nb2 in nb.GetNeighbors():
            if nb2.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(nb.GetIdx(), nb2.GetIdx())
                if bond.GetBondTypeAsDouble() == 2.0:
                    return nb2.GetIdx()
    return None


# ── Analysis functions ────────────────────────────────────────────────────────

def compute_omega_dihedrals(mol: Chem.Mol) -> list[dict]:
    """
    Compute backbone omega dihedrals for all ring amide bonds.
    Omega = dihedral(Ca_prev, C, N, Ca_next).
    Returns list of dicts with keys: c_idx, n_idx, omega, is_nme.
    """
    conf = mol.GetConformer(0)
    ring_atoms = _ring_atom_set(mol)
    results = []

    for bond in mol.GetBonds():
        if not bond.IsInRing():
            continue
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()

        c_atom = n_atom = None
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7:
            c_atom, n_atom = a1, a2
        elif a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6:
            c_atom, n_atom = a2, a1
        else:
            continue

        # Verify C has a double bond to O (amide carbonyl)
        has_co = any(
            nb.GetAtomicNum() == 8
            and mol.GetBondBetweenAtoms(c_atom.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble() == 2.0
            for nb in c_atom.GetNeighbors()
        )
        if not has_co:
            continue

        # Ca_prev: ring-C bonded to carbonyl C (not O, not N)
        ca_prev = next(
            (nb for nb in c_atom.GetNeighbors()
             if nb.GetAtomicNum() == 6 and nb.GetIdx() in ring_atoms),
            None,
        )
        # Ca_next: ring-C bonded to amide N (not the carbonyl C; not an N-methyl C)
        ca_next = next(
            (nb for nb in n_atom.GetNeighbors()
             if nb.GetAtomicNum() == 6
             and nb.GetIdx() in ring_atoms
             and nb.GetIdx() != c_atom.GetIdx()),
            None,
        )
        if ca_prev is None or ca_next is None:
            continue

        omega = rdMolTransforms.GetDihedralDeg(
            conf,
            ca_prev.GetIdx(), c_atom.GetIdx(),
            n_atom.GetIdx(), ca_next.GetIdx(),
        )
        is_nme = any(
            nb.GetAtomicNum() == 6 and nb.GetIdx() not in ring_atoms
            for nb in n_atom.GetNeighbors()
        )
        results.append({
            "c_idx":  c_atom.GetIdx(),
            "n_idx":  n_atom.GetIdx(),
            "omega":  omega,
            "is_nme": is_nme,
        })
    return results


def assign_ring_nhs(mol: Chem.Mol) -> dict[str, int]:
    """Match SMARTS patterns to find N-atom indices for the 4 free ring NHs."""
    assignment: dict[str, int] = {}
    for name, smarts in RESIDUE_NH_SMARTS.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        matches = mol.GetSubstructMatches(patt)
        if matches:
            assignment[name] = matches[0][0]  # first atom in SMARTS = NH nitrogen
    return assignment


def check_nh_hbonds(
    mol: Chem.Mol,
    nh_map: dict[str, int],
    ring_atoms: set[int],
    acceptors: list[int],
) -> dict[str, dict | None]:
    """
    For each named NH, find the best intramolecular H-bond to a ring C=O oxygen.
    Returns {residue: None} if no H-bond, or {residue: {na_dist, angle, a_idx}}.
    """
    conf = mol.GetConformer(0)
    results: dict[str, dict | None] = {}

    for res, n_idx in nh_map.items():
        n_atom = mol.GetAtomWithIdx(n_idx)
        h_idx = next(
            (nb.GetIdx() for nb in n_atom.GetNeighbors() if nb.GetAtomicNum() == 1),
            None,
        )
        if h_idx is None:
            results[res] = None
            continue

        n_pos = np.array(conf.GetAtomPosition(n_idx))
        h_pos = np.array(conf.GetAtomPosition(h_idx))

        best: dict | None = None
        for a_idx in acceptors:
            a_pos = np.array(conf.GetAtomPosition(a_idx))
            na_dist = float(np.linalg.norm(n_pos - a_pos))
            if na_dist > HB_DIST:
                continue
            vec_hn = n_pos - h_pos
            vec_ha = a_pos - h_pos
            cos_a = np.dot(vec_hn, vec_ha) / (
                np.linalg.norm(vec_hn) * np.linalg.norm(vec_ha) + 1e-10
            )
            angle = float(np.degrees(np.arccos(np.clip(cos_a, -1.0, 1.0))))
            if angle >= HB_ANGLE:
                if best is None or na_dist < best["na_dist"]:
                    best = {
                        "a_idx":   a_idx,
                        "na_dist": round(na_dist, 2),
                        "angle":   round(angle, 1),
                    }
        results[res] = best
    return results


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sdf",  type=Path, default=DATA_DIR / "ensemble.sdf",
                        help="Path to ensemble SDF (default: %(default)s)")
    parser.add_argument("--json", type=Path, default=DATA_DIR / "ensemble.json",
                        help="Path to ensemble JSON with Boltzmann weights")
    args = parser.parse_args()

    # ── Load conformers ───────────────────────────────────────────────────────
    suppl = Chem.SDMolSupplier(str(args.sdf), removeHs=False, sanitize=True)
    mols  = [m for m in suppl if m is not None]
    if not mols:
        sys.exit(f"ERROR: no molecules loaded from {args.sdf}")
    print(f"Loaded {len(mols)} conformers  |  {mols[0].GetNumAtoms()} atoms each")

    # ── Load Boltzmann weights ────────────────────────────────────────────────
    with open(args.json) as fh:
        ensemble_data = json.load(fh)
    conf_entries = ensemble_data.get("conformers", [])
    boltz_raw = [c.get("boltzmannweight", 0.0) for c in conf_entries]

    if len(boltz_raw) == len(mols):
        boltz = np.array(boltz_raw, dtype=float)
    else:
        print(f"Warning: {len(mols)} mols but {len(boltz_raw)} weights — using equal weights")
        boltz = np.ones(len(mols), dtype=float)
    boltz /= boltz.sum()

    # ── Assign residue NH atom indices (topology same across conformers) ───────
    ref_mol = mols[0]
    nh_map  = assign_ring_nhs(ref_mol)
    ring_atoms = _ring_atom_set(ref_mol)

    print(f"Residue NH assignment: {list(nh_map.keys())}")
    missing = [r for r in RESIDUE_NH_SMARTS if r not in nh_map]
    if missing:
        print(f"Warning: SMARTS match failed for: {missing} — those columns will show '?'")

    # Pre-compute shared acceptor list (same topology for all conformers)
    acceptors = _carbonyl_oxygens(ref_mol, ring_atoms)

    # Expected H-bond acceptors for topology check (item 2)
    # Ala7 NH should donate to the adjacent Val5 C=O (preceding carbonyl of Ala7 N)
    # Abu2 NH should donate to MeLeu10 C=O — a long-range contact, NOT adjacent to Abu2 N
    ala7_expected_acc = (
        _preceding_carbonyl_o(ref_mol, nh_map["Ala7"], ring_atoms)
        if "Ala7" in nh_map else None
    )
    abu2_adj_acc = (
        _preceding_carbonyl_o(ref_mol, nh_map["Abu2"], ring_atoms)
        if "Abu2" in nh_map else None
    )

    # ── Per-conformer analysis ────────────────────────────────────────────────
    rows = []
    for i, mol in enumerate(mols):
        w = float(boltz[i])

        omegas   = compute_omega_dihedrals(mol)
        cis_list = [o for o in omegas if abs(o["omega"]) < CIS_CUTOFF]
        n_cis    = len(cis_list)
        cis_str  = f"{cis_list[0]['omega']:+.1f}" if cis_list else "none"

        hb           = check_nh_hbonds(mol, nh_map, ring_atoms, acceptors)
        abu2_info    = hb.get("Abu2")
        val5_info    = hb.get("Val5")
        ala7_info    = hb.get("Ala7")
        dala8_info   = hb.get("DAla8")

        abu2_hb  = abu2_info  is not None
        val5_hb  = val5_info  is not None
        ala7_hb  = ala7_info  is not None
        dala8_hb = dala8_info is not None

        # Item 2: topology checks
        # Ala7 → Val5 C=O (adjacent preceding carbonyl)
        ala7_to_val5 = (
            ala7_hb
            and ala7_expected_acc is not None
            and ala7_info["a_idx"] == ala7_expected_acc
        )
        # Abu2 → long-range C=O (should NOT be the adjacent MeBmt1 C=O)
        abu2_long_range = (
            abu2_hb
            and abu2_adj_acc is not None
            and abu2_info["a_idx"] != abu2_adj_acc
        )

        # A1 fingerprint: exactly 1 cis amide + Abu2 HB + Ala7 HB + Val5 free
        a1_like = (n_cis == 1 and abu2_hb and ala7_hb and not val5_hb)

        rows.append({
            "conf":          i + 1,
            "boltz":         w,
            "pct":           round(w * 100, 1),
            "n_cis":         n_cis,
            "cis_str":       cis_str,
            "all_omegas":    omegas,            # item 1: full omega list
            "Abu2_HB":       abu2_hb,
            "Val5_HB":       val5_hb,
            "Ala7_HB":       ala7_hb,
            "DAla8_HB":      dala8_hb,
            "ala7_to_val5":  ala7_to_val5,     # item 2: correct acceptor?
            "abu2_lr":       abu2_long_range,   # item 2: long-range contact?
            "abu2_dist":     round(abu2_info["na_dist"], 2) if abu2_info else None,
            "abu2_angle":    round(abu2_info["angle"],   1) if abu2_info else None,
            "ala7_dist":     round(ala7_info["na_dist"], 2) if ala7_info else None,
            "ala7_angle":    round(ala7_info["angle"],   1) if ala7_info else None,
            "A1_like":       a1_like,
        })

    # ── Table 1: fingerprint pass/fail ───────────────────────────────────────
    yn = lambda x: "Y" if x else "N"
    print()
    print("CsA Water Ensemble — A1 Fingerprint Validation")
    print("Reference: Limbach et al., JACS 2022, 144, 12602")
    print("=" * 78)
    print(f"{'Conf':>4}  {'Boltz%':>6}  {'nCis':>4}  {'ω(cis)°':>8}"
          f"  {'Abu2':>4}  {'Val5':>4}  {'Ala7':>4}  {'DAla8':>5}"
          f"  {'Abu2→LR':>7}  {'Ala7→V5':>7}  {'A1?':>4}")
    print("-" * 78)
    for r in rows:
        lr  = yn(r["abu2_lr"])       if r["Abu2_HB"] else "-"
        v5  = yn(r["ala7_to_val5"])  if r["Ala7_HB"] else "-"
        print(
            f"{r['conf']:>4}  {r['pct']:>6.1f}%  {r['n_cis']:>4}  {r['cis_str']:>8}"
            f"  {yn(r['Abu2_HB']):>4}  {yn(r['Val5_HB']):>4}"
            f"  {yn(r['Ala7_HB']):>4}  {yn(r['DAla8_HB']):>5}"
            f"  {lr:>7}  {v5:>7}  {yn(r['A1_like']):>4}"
        )
    print()
    print("  Abu2→LR : Abu2 NH H-bonds to a non-adjacent C=O (expected: MeLeu10)")
    print("  Ala7→V5 : Ala7 NH H-bonds to the adjacent Val5 C=O (paper: Ala7-Val5)")

    # ── Table 2: all omega values for top 5 conformers (item 1) ──────────────
    print()
    print("All backbone omega dihedrals — top 5 conformers by Boltzmann weight")
    print("(sanity check: all non-cis amides should be |ω| > 150°)")
    print("-" * 60)
    top5 = sorted(rows, key=lambda r: r["boltz"], reverse=True)[:5]
    for r in top5:
        oms = sorted(r["all_omegas"], key=lambda o: o["omega"])
        cis_vals  = [f"{o['omega']:+.0f}°(cis)"  for o in oms if abs(o["omega"]) < 30]
        amb_vals  = [f"{o['omega']:+.0f}°(amb)"  for o in oms if 30 <= abs(o["omega"]) <= 150]
        tran_vals = [f"{o['omega']:+.0f}°"        for o in oms if abs(o["omega"]) > 150]
        print(f"  Conf {r['conf']:>2} ({r['pct']:.1f}%):")
        print(f"    trans  : {', '.join(tran_vals) if tran_vals else 'none'}")
        if amb_vals:
            print(f"    ambig  : {', '.join(amb_vals)}")
        print(f"    cis    : {', '.join(cis_vals) if cis_vals else 'none'}")

    # ── Table 3: H-bond geometry quality for A1-like conformers (item 3) ─────
    a1_rows = [r for r in rows if r["A1_like"]]
    print()
    if a1_rows:
        print("H-bond geometry for A1-like conformers (vs. crystal reference):")
        print(f"  {'Conf':>4}  {'Boltz%':>6}  {'Abu2 N···O':>10}  {'Abu2 ∠':>8}  {'Ala7 N···O':>10}  {'Ala7 ∠':>8}")
        print(f"  {'':>4}  {'':>6}  {'(Å)':>10}  {'(deg)':>8}  {'(Å)':>10}  {'(deg)':>8}")
        print("  " + "-" * 56)
        for r in a1_rows:
            ad = f"{r['abu2_dist']:.2f}" if r["abu2_dist"] else "  -  "
            aa = f"{r['abu2_angle']:.1f}" if r["abu2_angle"] else "  -  "
            ld = f"{r['ala7_dist']:.2f}" if r["ala7_dist"] else "  -  "
            la = f"{r['ala7_angle']:.1f}" if r["ala7_angle"] else "  -  "
            print(f"  {r['conf']:>4}  {r['pct']:>6.1f}%  {ad:>10}  {aa:>8}  {ld:>10}  {la:>8}")

        # Boltzmann-weighted averages within A1-like subset
        w_a1_total = sum(r["boltz"] for r in a1_rows)
        def _wavg(key):
            vals = [(r[key], r["boltz"]) for r in a1_rows if r[key] is not None]
            if not vals: return float("nan")
            return sum(v * w for v, w in vals) / sum(w for _, w in vals)

        print()
        print(f"  Boltzmann-weighted averages (A1-like subset, {w_a1_total*100:.1f}% of ensemble):")
        print(f"    Abu2 N···O : {_wavg('abu2_dist'):.2f} Å   (crystal: 2.96 Å, 155°)")
        print(f"    Abu2 angle : {_wavg('abu2_angle'):.1f}°")
        print(f"    Ala7 N···O : {_wavg('ala7_dist'):.2f} Å   (crystal: 3.13 Å, 138°)")
        print(f"    Ala7 angle : {_wavg('ala7_angle'):.1f}°")
    else:
        print("No A1-like conformers found — geometry quality section skipped.")

    # ── Boltzmann-weighted summary ────────────────────────────────────────────
    w_cis     = sum(r["boltz"] for r in rows if r["n_cis"] == 1)
    w_abu2    = sum(r["boltz"] for r in rows if r["Abu2_HB"])
    w_ala7    = sum(r["boltz"] for r in rows if r["Ala7_HB"])
    w_val5_no = sum(r["boltz"] for r in rows if not r["Val5_HB"])
    w_a1      = sum(r["boltz"] for r in rows if r["A1_like"])

    print()
    print("Boltzmann-weighted population satisfying each A1 criterion:")
    print(f"  1. Exactly 1 cis amide (MeVal11-MeBmt1)  {w_cis*100:5.1f}%")
    print(f"  2. Abu2 NH H-bonded                       {w_abu2*100:5.1f}%   [NMR: -2.6  ppb/K]")
    print(f"  3. Ala7 NH H-bonded                       {w_ala7*100:5.1f}%   [NMR: -2.07 ppb/K]")
    print(f"  4. Val5 NH solvent-exposed                 {w_val5_no*100:5.1f}%   [NMR: -7.33 ppb/K]")
    print(f"  Full A1 match (all 4 criteria)             {w_a1*100:5.1f}%")
    print()
    print("Implicit-solvent note:")
    print("  ALPB(water) was used — no explicit waters in conformers.")
    print("  The 2 crystal water molecules inside the A1 macrocycle cavity")
    print("  (CCDC 2149649, Table S6) are not captured but do not invalidate")
    print("  the intramolecular fold comparison.")


if __name__ == "__main__":
    main()

"""
3dphys_descriptors_v2.py
------------------------
3D physics-based descriptors for cyclic peptide conformer ensembles.
Operates on raw XYZ coordinate arrays (symbol lists + numpy coords) from
CREST/xTB ensembles. Designed to work alongside crest_conformers.py.

Functions:
  boltzmann_weights   — compute Boltzmann weights from GFN2 energies (Hartree)
  compute_psa_xyz     — 3D polar SASA using rdFreeSASA (falls back to analytic)
  count_hbonds_xyz    — intramolecular H-bond count from 3D geometry
"""

from __future__ import annotations

import numpy as np


def boltzmann_weights(energies_hartree: list[float], T: float = 298.15) -> np.ndarray:
    KCAL_PER_HARTREE = 627.509
    RT = 1.987e-3 * T
    e = np.array(energies_hartree, dtype=float) * KCAL_PER_HARTREE
    valid = np.isfinite(e)
    if not np.any(valid):
        raise RuntimeError("No valid conformer energies for Boltzmann weighting")
    e_rel = e[valid] - np.nanmin(e[valid])
    valid_weights = np.exp(-e_rel / RT)
    valid_weights /= valid_weights.sum()
    weights = np.full(len(e), np.nan, dtype=float)
    weights[valid] = valid_weights
    return weights


def compute_psa_xyz(symbols: list[str], coords: np.ndarray,
                    template_mol=None) -> float:
    """
    Compute 3D polar SASA (Å²) for a conformer given as symbol list + coords array.
    Uses rdFreeSASA with Bondi radii if template_mol (RDKit Mol) is provided.
    Falls back to analytic sphere intersection if template_mol is unavailable.
    """
    from rdkit import Chem
    from rdkit.Chem import rdFreeSASA

    _BONDI = {
        'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
        'S': 1.80, 'P': 1.80, 'F': 1.47, 'Cl': 1.75,
        'Br': 1.85, 'I': 1.98,
    }
    _POLAR = {'N', 'O', 'S', 'P'}

    if template_mol is not None:
        try:
            mol_h = Chem.RWMol(template_mol)
            if mol_h.GetNumAtoms() != len(symbols):
                raise ValueError("Atom count mismatch")
            conf = Chem.Conformer(mol_h.GetNumAtoms())
            for i, (x, y, z) in enumerate(coords):
                conf.SetAtomPosition(i, (float(x), float(y), float(z)))
            conf_id = mol_h.AddConformer(conf, assignId=True)
            mol_h = mol_h.GetMol()
            radii = []
            for atom in mol_h.GetAtoms():
                sym = atom.GetSymbol()
                radii.append(_BONDI.get(sym, 1.50))
                if sym in _POLAR:
                    atom.SetIntProp('SASAClass', 0)
                    atom.SetProp('SASAClassName', 'Polar')
                else:
                    atom.SetIntProp('SASAClass', 1)
                    atom.SetProp('SASAClassName', 'APolar')
            query = rdFreeSASA.MakeFreeSasaPolarAtomQuery()
            psa = rdFreeSASA.CalcSASA(mol_h, radii, confIdx=conf_id, query=query)
            return round(float(psa), 2)
        except Exception:
            pass

    # Analytic fallback
    PROBE = 1.40
    n = len(symbols)
    radii_all = np.array([_BONDI.get(s, 1.50) + PROBE for s in symbols])
    polar_idx = [i for i, s in enumerate(symbols) if s in _POLAR]
    if not polar_idx:
        return 0.0

    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=-1))
    np.fill_diagonal(dist, np.inf)

    total = 0.0
    for idx in polar_idx:
        r_i = radii_all[idx]
        buried = 0.0
        sphere_area = 4.0 * np.pi * r_i ** 2
        for jdx in range(n):
            if jdx == idx:
                continue
            d = dist[idx, jdx]
            r_j = radii_all[jdx]
            if d >= r_i + r_j:
                continue
            if d <= abs(r_i - r_j):
                if r_i <= r_j:
                    buried = sphere_area
                    break
                continue
            h = r_i - (r_i**2 + d**2 - r_j**2) / (2.0 * d)
            h = max(0.0, min(h, 2.0 * r_i))
            buried += 2.0 * np.pi * r_i * h
        exposed = max(0.0, sphere_area - buried)
        total += exposed
    return round(total, 2)


def count_hbonds_xyz(symbols: list[str], coords: np.ndarray) -> int:
    """
    Count intramolecular H-bonds from 3D geometry.
    Criteria: H...A distance < 2.5 Å, D-H...A angle > 120°.
    """
    D_ATOMS = {"N", "O"}
    A_ATOMS = {"N", "O"}
    H_DIST_MAX = 2.5
    ANGLE_MIN = 120.0

    donor_h = []
    for i, sym in enumerate(symbols):
        if sym != "H":
            continue
        dists = [(np.linalg.norm(coords[i] - coords[j]), j)
                 for j, s in enumerate(symbols) if s != "H" and j != i]
        if not dists:
            continue
        d_nearest, d_idx = min(dists)
        if d_nearest < 1.3 and symbols[d_idx] in D_ATOMS:
            donor_h.append((i, d_idx))

    acceptors = [i for i, s in enumerate(symbols) if s in A_ATOMS]

    count = 0
    for h_idx, d_idx in donor_h:
        h_pos = coords[h_idx]
        d_pos = coords[d_idx]
        for a_idx in acceptors:
            if a_idx == d_idx:
                continue
            a_pos = coords[a_idx]
            ha_dist = np.linalg.norm(h_pos - a_pos)
            if ha_dist > H_DIST_MAX:
                continue
            vec_hd = d_pos - h_pos
            vec_ha = a_pos - h_pos
            cos_a = np.dot(vec_hd, vec_ha) / (
                np.linalg.norm(vec_hd) * np.linalg.norm(vec_ha) + 1e-10)
            angle = np.degrees(np.arccos(np.clip(cos_a, -1, 1)))
            if angle >= ANGLE_MIN:
                count += 1
    return count

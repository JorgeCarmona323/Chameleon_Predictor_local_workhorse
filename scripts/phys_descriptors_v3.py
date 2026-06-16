# env: chameleon-calc
"""
phys_descriptors_v3.py
----------------------
3D physics-based descriptors for cyclic peptide conformer ensembles — v3.

ADDITIVE SUPERSET of phys_descriptors_v2. Nothing from v2 is changed or removed:
the v2 functions are imported and re-exported here so existing callers can switch
`import phys_descriptors_v2` -> `import phys_descriptors_v3` with no behavior change.

v3 adds the literature-validated surface/shape descriptors we were missing
(see docs/experiments/2026-06-13_descriptor_literature_review.md). All are
computable from ensembles we already have — no new CREST runs:

  surface_descriptors_mol / surface_descriptors_xyz
      psa                 — 3D polar SASA  (N/O/S/P)            [== v2 compute_psa_xyz target]
      hbd_sasa            — SA_HD: SASA of H-bond DONOR H's     [Rzepiela 2022, top descriptor]
      hba_sasa            — SA_HA: SASA of H-bond ACCEPTOR atoms (N/O)
      hydrophobic_sasa    — ASA_H: SASA of apolar atoms (C, H-on-C) [Rzepiela 2022]
      total_sasa          — full molecular SASA
      amphi_moment        — Å separation of polar vs apolar SASA centroids
                            [García Jiménez 2024, integy/amphipathic moment]
  imhb_descriptors_mol    — intramolecular H-bond breakdown (geometric, same criterion
                            as count_hbonds_xyz): imhb (total), imhbd/imhba (engaged
                            donors/acceptors), imhb_bb/imhb_res (backbone vs side-chain)
  effective_nconf         — exp(Shannon entropy of Boltzmann weights):
                            in-ensemble flexibility analog of nConf20 [Wicker & Cooper]

The rationale (keep everything, add the missing, benchmark to decide) is recorded
in the literature-review doc. v3 is purely additive so the in-house ML benchmark
can train on the full column set and let feature importance pick the final set.

Re-exported from v2: boltzmann_weights, compute_psa_xyz, count_hbonds_xyz.
"""

from __future__ import annotations

import numpy as np

# Re-export v2 unchanged so v3 is a drop-in superset.
from phys_descriptors_v2 import (  # noqa: F401
    boltzmann_weights,
    compute_psa_xyz,
    count_hbonds_xyz,
)

# Bondi van der Waals radii (Å) — same table as v2.
_BONDI = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52,
    "S": 1.80, "P": 1.80, "F": 1.47, "Cl": 1.75,
    "Br": 1.85, "I": 1.98,
}
_POLAR = {"N", "O", "S", "P"}   # polar heavy atoms (PSA atoms)
_DONOR_HEAVY = {"N", "O"}        # H-bond donor heavy atoms
_ACCEPTOR = {"N", "O"}           # H-bond acceptor heavy atoms

_HB_DIST_MAX = 2.5               # H...acceptor distance cutoff (Å)
_HB_ANGLE_MIN = 120.0            # D-H...A angle cutoff (deg)


# ── per-atom SASA ─────────────────────────────────────────────────────────────
def _per_atom_sasa(mol, conf_id: int) -> np.ndarray | None:
    """Return per-atom SASA (Å²) array indexed by atom, or None on failure.

    Computes the full molecular SASA once with rdFreeSASA (Bondi radii); rdFreeSASA
    stores each atom's contribution as the double property 'SASA', which we read back
    and subset by atom class. One CalcSASA call serves PSA, SA_HBD, ASA_H and the
    amphipathic moment — they are all partitions of the same per-atom surface.
    """
    from rdkit.Chem import rdFreeSASA
    try:
        radii = [_BONDI.get(a.GetSymbol(), 1.50) for a in mol.GetAtoms()]
        rdFreeSASA.CalcSASA(mol, radii, confIdx=conf_id)
        return np.array([float(a.GetPropsAsDict().get("SASA", 0.0))
                         for a in mol.GetAtoms()], dtype=float)
    except Exception:
        return None


def _atom_classes(mol):
    """Per-atom boolean masks: (is_polar_heavy, is_donor_H, is_apolar, is_acceptor).

    donor H  = H bonded to N or O.
    acceptor = N/O heavy atoms.
    apolar   = C, or H bonded to C.
    polar    = N/O/S/P heavy atoms (PSA set).
    """
    n = mol.GetNumAtoms()
    is_polar = np.zeros(n, dtype=bool)
    is_donor_h = np.zeros(n, dtype=bool)
    is_apolar = np.zeros(n, dtype=bool)
    is_acceptor = np.zeros(n, dtype=bool)
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        sym = atom.GetSymbol()
        if sym in _POLAR:
            is_polar[i] = True
            if sym in _ACCEPTOR:
                is_acceptor[i] = True
        elif sym == "C":
            is_apolar[i] = True
        elif sym == "H":
            nbrs = atom.GetNeighbors()
            if nbrs:
                hsym = nbrs[0].GetSymbol()
                if hsym in _DONOR_HEAVY:
                    is_donor_h[i] = True
                elif hsym == "C":
                    is_apolar[i] = True
    return is_polar, is_donor_h, is_apolar, is_acceptor


def surface_descriptors_mol(mol, conf_id: int = -1) -> dict:
    """v3 surface/shape descriptors for one RDKit conformer (mol must have Hs + bonds).

    Returns dict with keys: psa, hbd_sasa (SA_HD), hba_sasa (SA_HA),
    hydrophobic_sasa, total_sasa, amphi_moment. NaN if the SASA calculation fails.
    """
    if conf_id == -1:
        conf_id = mol.GetConformer().GetId()
    nan = float("nan")
    out = {"psa": nan, "hbd_sasa": nan, "hba_sasa": nan, "hydrophobic_sasa": nan,
           "total_sasa": nan, "amphi_moment": nan}

    sasa = _per_atom_sasa(mol, conf_id)
    if sasa is None:
        return out

    is_polar, is_donor_h, is_apolar, is_acceptor = _atom_classes(mol)
    out["total_sasa"] = round(float(sasa.sum()), 2)
    out["psa"] = round(float(sasa[is_polar].sum()), 2)
    out["hbd_sasa"] = round(float(sasa[is_donor_h].sum()), 2)   # SA_HD: donor-H surface
    out["hba_sasa"] = round(float(sasa[is_acceptor].sum()), 2)  # SA_HA: acceptor-atom surface
    out["hydrophobic_sasa"] = round(float(sasa[is_apolar].sum()), 2)

    # Amphipathic moment: distance between SASA-weighted centroids of polar vs
    # apolar surface. Large = polar/nonpolar surface is spatially segregated.
    coords = mol.GetConformer(conf_id).GetPositions()
    out["amphi_moment"] = round(_amphi_moment(coords, sasa, is_polar, is_apolar), 3)
    return out


def _amphi_moment(coords, sasa, is_polar, is_apolar) -> float:
    wp = sasa * is_polar
    wa = sasa * is_apolar
    sp, sa = wp.sum(), wa.sum()
    if sp <= 0 or sa <= 0:
        return float("nan")
    c_polar = (coords * wp[:, None]).sum(axis=0) / sp
    c_apolar = (coords * wa[:, None]).sum(axis=0) / sa
    return float(np.linalg.norm(c_polar - c_apolar))


def surface_descriptors_xyz(symbols: list[str], coords: np.ndarray,
                            template_mol) -> dict:
    """XYZ-array parity wrapper for surface_descriptors_mol.

    Needs a template RDKit Mol (with bonds) to classify donors; embeds `coords`
    as a conformer on a copy of it. Returns NaNs if no usable template is given.
    """
    from rdkit import Chem
    nan = float("nan")
    if template_mol is None:
        return {"psa": nan, "hbd_sasa": nan, "hba_sasa": nan, "hydrophobic_sasa": nan,
                "total_sasa": nan, "amphi_moment": nan}
    mol = Chem.RWMol(template_mol)
    if mol.GetNumAtoms() != len(symbols):
        return {"psa": nan, "hbd_sasa": nan, "hba_sasa": nan, "hydrophobic_sasa": nan,
                "total_sasa": nan, "amphi_moment": nan}
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, (float(x), float(y), float(z)))
    mol = mol.GetMol()
    cid = mol.AddConformer(conf, assignId=True)
    return surface_descriptors_mol(mol, cid)


# ── intramolecular H-bond breakdown ───────────────────────────────────────────
def macrocycle_atoms(mol) -> set:
    """Atom indices of the largest ring = the macrocycle backbone.
    Side-chain rings (thiophene, benzene, proline) are small and excluded."""
    rings = mol.GetRingInfo().AtomRings()
    return set(max(rings, key=len)) if rings else set()


def backbone_hbond_atoms(mol) -> set:
    """Atoms counted as 'backbone' for H-bond classification: the macrocycle ring PLUS
    carbonyl/hydroxyl O's directly bonded to a ring atom. The amide C=O oxygen points
    OUT of the ring (it is not itself a ring atom), so without this a transannular
    backbone N-H...O=C would be misclassified as side-chain."""
    ring = macrocycle_atoms(mol)
    bb = set(ring)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "O" and any(nb.GetIdx() in ring for nb in atom.GetNeighbors()):
            bb.add(atom.GetIdx())
    return bb


def imhb_descriptors_mol(mol, conf_id: int = -1, backbone_atoms: set | None = None) -> dict:
    """Intramolecular H-bond breakdown for one conformer (mol must have Hs + bonds).

    Same geometric criterion as count_hbonds_xyz (H...A < 2.5 Å, D-H...A > 120°),
    but delineated by role and location:
      imhb     — total intramolecular H-bonds (== bw_hb predecessor)
      imhbd    — distinct DONOR H's engaged in >=1 IMHB
      imhba    — distinct ACCEPTOR atoms engaged in >=1 IMHB
      imhb_bb  — IMHBs where BOTH partners are macrocycle-ring (backbone/transannular)
      imhb_res — IMHBs involving >=1 side-chain (residue) atom
    imhb_bb + imhb_res == imhb.
    """
    if conf_id == -1:
        conf_id = mol.GetConformer().GetId()
    if backbone_atoms is None:
        backbone_atoms = backbone_hbond_atoms(mol)
    coords = mol.GetConformer(conf_id).GetPositions()

    donors, acceptors = [], []
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym in _ACCEPTOR:
            acceptors.append(atom.GetIdx())
        elif sym == "H":
            nbrs = atom.GetNeighbors()
            if nbrs and nbrs[0].GetSymbol() in _DONOR_HEAVY:
                donors.append((atom.GetIdx(), nbrs[0].GetIdx()))

    imhb = bb = res = 0
    donor_set, acc_set = set(), set()
    for h_idx, d_idx in donors:
        h, d = coords[h_idx], coords[d_idx]
        for a_idx in acceptors:
            if a_idx == d_idx:
                continue
            a = coords[a_idx]
            if np.linalg.norm(h - a) > _HB_DIST_MAX:
                continue
            v1, v2 = d - h, a - h
            cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
            if np.degrees(np.arccos(np.clip(cos, -1, 1))) < _HB_ANGLE_MIN:
                continue
            imhb += 1
            donor_set.add(h_idx)
            acc_set.add(a_idx)
            if d_idx in backbone_atoms and a_idx in backbone_atoms:
                bb += 1
            else:
                res += 1
    return {"imhb": imhb, "imhbd": len(donor_set), "imhba": len(acc_set),
            "imhb_bb": bb, "imhb_res": res}


# ── ensemble flexibility ──────────────────────────────────────────────────────
def effective_nconf(weights) -> float:
    """Effective number of populated conformers = exp(Shannon entropy of weights).

    In-ensemble analog of nConf20 (Wicker & Cooper): n_eff = 1 for a single
    dominant conformer, up to N for a perfectly flat ensemble. Complements
    p_dominant — p_dominant gives the top peak; n_eff gives how many states matter.
    """
    w = np.asarray(weights, dtype=float)
    w = w[np.isfinite(w) & (w > 0)]
    if w.size == 0:
        return float("nan")
    w = w / w.sum()
    h = -(w * np.log(w)).sum()
    return round(float(np.exp(h)), 3)

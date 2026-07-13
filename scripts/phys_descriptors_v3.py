# env: chameleon-calc
"""
phys_descriptors_v3.py
----------------------
3D physics-based descriptors for cyclic peptide conformer ensembles — v3.

Canonical 3D-descriptor library. The former phys_descriptors_v2 base functions
(boltzmann_weights, compute_psa_xyz, count_hbonds_xyz) are folded in below unchanged;
v2 itself has been retired to scripts/_archive/. Import everything from here.

v3 adds the literature-validated surface/shape descriptors we were missing
(see docs/experiments/2026-06-13_descriptor_literature_review.md). All are
computable from ensembles we already have — no new CREST runs:

  surface_descriptors_mol / surface_descriptors_xyz
      psa                 — 3D PSA: SASA over N/O + polar H + oxidized S [Ono 2019/Begnini 2021 def;
                            reduced S (thiophene/thioether/thiol) excluded per Ertl TPSA convention]
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
  weighted_rmsf           — Boltzmann-weighted RMSF (threshold-free ensemble flexibility, Å)
  kier_flexibility        — Kier Φ = κ₁·κ₂/N_heavy (2D molecular flexibility) [Begnini 2021]

The rationale (keep everything, add the missing, benchmark to decide) is recorded
in the literature-review doc. v3 is purely additive so the in-house ML benchmark
can train on the full column set and let feature importance pick the final set.

Re-exported from v2: boltzmann_weights, compute_psa_xyz, count_hbonds_xyz.
"""

from __future__ import annotations

import numpy as np

# ── Base descriptors (folded in from the former phys_descriptors_v2) ──────────
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
    # 3D PSA — Ono 2019 / Begnini 2021 definition: SASA over N/O acceptors + polar H attached
    # to N/O, + OXIDIZED sulfur only (S=O: sulfoxide/sulfone/sulfonamide). Reduced S (thiophene,
    # thioether, thiol, disulfide) contributes 0 — the Ertl TPSA sulfur convention.
    oxidized_s = np.array([
        a.GetSymbol() == "S" and any(
            b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(a).GetSymbol() == "O"
            for b in a.GetBonds())
        for a in mol.GetAtoms()])
    out["total_sasa"] = round(float(sasa.sum()), 2)
    out["psa"] = round(float(sasa[is_acceptor | is_donor_h | oxidized_s].sum()), 2)
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


def weighted_rmsf(mols, weights, heavy_only: bool = True) -> float:
    """Boltzmann-weighted root-mean-square fluctuation (Å) of the ensemble.

    The clean, THRESHOLD-FREE, count-free flexibility descriptor: aligns all conformers
    to the highest-weight one, then measures the weighted geometric spread about the
    weighted-mean structure. Unlike p_dominant / n_eff it does not depend on how finely
    CREST discretizes a basin (near-duplicate conformers sit on the mean and contribute
    ~0), nor on a clustering RMSD threshold. High = floppy, low = rigid.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    w = np.asarray(weights, dtype=float)
    keep = np.where(np.isfinite(w) & (w > 0))[0]
    if keep.size < 2:
        return float("nan")
    ref = int(keep[np.argmax(w[keep])])           # highest-weight conformer = alignment frame
    order = [ref] + [int(i) for i in keep if i != ref]

    base = Chem.Mol(mols[ref]); base.RemoveAllConformers()
    for i in order:
        base.AddConformer(Chem.Conformer(mols[i].GetConformer()), assignId=True)
    atom_ids = [a.GetIdx() for a in base.GetAtoms()
                if (a.GetAtomicNum() > 1 or not heavy_only)]
    try:
        rdMolAlign.AlignMolConformers(base, atomIds=atom_ids)
    except Exception:
        return float("nan")

    P = np.array([c.GetPositions()[atom_ids] for c in base.GetConformers()])  # (N, A, 3)
    ww = w[order]; ww = ww / ww.sum()
    mean = (ww[:, None, None] * P).sum(axis=0)                                 # (A, 3)
    msf = (ww[:, None] * ((P - mean) ** 2).sum(axis=2)).sum(axis=0)            # (A,)
    return round(float(np.sqrt(msf.mean())), 3)


def kier_flexibility(mol) -> float:
    """Kier molecular flexibility index Φ = κ₁·κ₂ / N, with κ₁, κ₂ on the hydrogen-suppressed
    graph and N the TOTAL atom count (including H) — MOE's `KierFlex` normalization. This is the
    scale on which Begnini 2021 / Poongavanam 2021 report Φ ≈ 7–9 for macrocyclic peptides and
    the Φ ≲ 10 ceiling below which conformational sampling reliably ranks macrocycle permeability.
    (The original Kier 1989 heavy-atom normalization runs ~1.8× higher and is NOT comparable to
    that threshold.) Topological (2D) — identical for stereoisomers; varies with the backbone
    graph (azetidine vs sarcosine)."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    mh = Chem.RemoveHs(mol)
    n = Chem.AddHs(mh).GetNumAtoms()
    if n == 0:
        return float("nan")
    return round(float(rdMolDescriptors.CalcKappa1(mh) * rdMolDescriptors.CalcKappa2(mh) / n), 3)

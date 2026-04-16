"""
cremp_deltapsa.py
-----------------
Compute ΔPSA from CREMP pre-computed conformational ensembles.

For each compound in pickle.tar.gz:
  1. Load the RDKit mol (rd_mol) — already has N unique conformers embedded
  2. Compute polar SASA for every conformer using Bondi radii (same as conformer_engine.py)
  3. Identify aqueous conformer (max PSA) and membrane conformer (min PSA)
  4. Compute Δ descriptors and Boltzmann-weighted PSA

Output:
  results/cremp_deltapsa.csv — one row per CREMP compound, SMILES-joinable to feature_matrix.csv

Usage:
  python scripts/cremp_deltapsa.py
  python scripts/cremp_deltapsa.py --pickle dependencies/pickle.tar.gz --outdir results
"""

import argparse
import tarfile
import pickle
import warnings
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFreeSASA, rdMolDescriptors
from tqdm import tqdm

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

# ── Bondi radii + polar atom definition (identical to conformer_engine.py) ────
_BONDI = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
    'S': 1.80, 'P': 1.80, 'F': 1.47, 'Cl': 1.75,
    'Br': 1.85, 'I': 1.98,
}
_POLAR_ELEMENTS = {'N', 'O', 'S', 'P'}


def _polar_sasa(mol: Chem.Mol, conf_id: int) -> float:
    """
    Polar SASA for one conformer using Bondi radii.
    Exact copy of the function in conformer_engine.py — do not diverge.
    Returns np.nan on failure.
    """
    try:
        radii = []
        for atom in mol.GetAtoms():
            sym = atom.GetSymbol()
            radii.append(_BONDI.get(sym, 1.50))
            if sym in _POLAR_ELEMENTS:
                atom.SetIntProp('SASAClass', 0)
                atom.SetProp('SASAClassName', 'Polar')
            else:
                atom.SetIntProp('SASAClass', 1)
                atom.SetProp('SASAClassName', 'APolar')
        query = rdFreeSASA.MakeFreeSasaPolarAtomQuery()
        return round(rdFreeSASA.CalcSASA(mol, radii, confIdx=conf_id, query=query), 4)
    except Exception:
        return np.nan


def _total_sasa(mol: Chem.Mol, conf_id: int) -> float:
    """Total SASA (all atoms) for normalizing ΔPSA per Yu et al. 2026."""
    try:
        radii = [_BONDI.get(atom.GetSymbol(), 1.50) for atom in mol.GetAtoms()]
        return round(rdFreeSASA.CalcSASA(mol, radii, confIdx=conf_id), 4)
    except Exception:
        return np.nan


def process_pickle(name: str, obj: dict, max_confs: int = 100) -> dict:
    """
    Compute ΔPSA descriptors for one CREMP compound.
    name      — filename stem (e.g. 'meA.L.bHph.P.L.F')
    obj       — deserialized pickle dict
    max_confs — max conformers to evaluate (sampled uniformly; 100 ≈ 6% of 1774
                which is sufficient to find near-minimum PSA)
    """
    compound_id = Path(name).stem

    mol = obj.get('rd_mol')
    smiles = obj.get('smiles', '')
    conformers_meta = obj.get('conformers', [])  # list of dicts with boltzmannweight

    if mol is None or mol.GetNumConformers() == 0:
        return {'compound_id': compound_id, 'smiles': smiles, 'error': 'no_rdmol'}

    n_confs = mol.GetNumConformers()

    # ── Sample conformers ─────────────────────────────────────────────────────
    # Strategy: take the lowest-energy conformers (first in list, relativeenergy=0
    # is set 1) PLUS a uniform sample across the rest. This ensures we catch the
    # most collapsed CHCl₃ conformers (low energy = most populated in membrane env)
    # while also sampling conformational diversity for psa3d_std.
    if n_confs <= max_confs:
        sample_ids = list(range(n_confs))
    else:
        # Always include first 20 (lowest energy) + uniform sample of the rest
        n_top = min(20, max_confs // 2)
        n_rand = max_confs - n_top
        top_ids = list(range(n_top))
        rest = list(range(n_top, n_confs))
        step = max(1, len(rest) // n_rand)
        rand_ids = rest[::step][:n_rand]
        sample_ids = sorted(set(top_ids + rand_ids))

    # ── Per-conformer polar SASA ──────────────────────────────────────────────
    psas = []
    for i in sample_ids:
        psas.append(_polar_sasa(mol, i))

    psas = np.array(psas, dtype=float)
    valid_mask = ~np.isnan(psas)

    if valid_mask.sum() == 0:
        return {'compound_id': compound_id, 'smiles': smiles, 'error': 'psa_failed'}

    psas_valid = psas[valid_mask]
    valid_indices = np.where(valid_mask)[0]

    aq_psas_idx  = valid_indices[np.argmax(psas_valid)]   # index into psas[] / sample_ids[]
    mem_psas_idx = valid_indices[np.argmin(psas_valid)]

    aq_conf_id  = sample_ids[int(aq_psas_idx)]   # actual RDKit conformer ID
    mem_conf_id = sample_ids[int(mem_psas_idx)]

    aq_psa  = float(psas[aq_psas_idx])
    mem_psa = float(psas[mem_psas_idx])

    # ── Boltzmann-weighted PSA (sampled subset only) ──────────────────────────
    bw_psa = np.nan
    if len(conformers_meta) == n_confs:
        weights_all = np.array([conformers_meta[i].get('boltzmannweight', 0.0)
                                for i in sample_ids], dtype=float)
        w_subset = weights_all[valid_mask]
        w_sum = w_subset.sum()
        if w_sum > 0:
            bw_psa = float(np.sum(w_subset * psas_valid) / w_sum)

    # ── Normalized ΔPSA (Yu et al. 2026) — ΔPSA / SASA_aq_total ─────────────
    # Denominator is total SASA of the aq conformer (Yu 2026 definition).
    # Note: bw_psa3d is computed over ~100 sampled conformers (59% of total BW
    # weight for typical compounds). Sampling systematically overestimates
    # bw_psa3d by ~3-4 Å² relative to the full ensemble, because the tail of
    # medium-energy conformers is underrepresented. Accepted approximation.
    total_aq  = _total_sasa(mol, aq_conf_id)
    norm_delta_psa = np.nan
    if not np.isnan(total_aq) and total_aq > 0:
        norm_delta_psa = round((aq_psa - mem_psa) / total_aq, 6)

    # ── 2D TPSA (reference) ───────────────────────────────────────────────────
    tpsa_2d = float(rdMolDescriptors.CalcTPSA(mol))

    return {
        'compound_id':       compound_id,
        'smiles':            smiles,
        'n_confs':           n_confs,
        'temperature_K':     obj.get('temperature', np.nan),
        'aq_psa3d':          aq_psa,
        'mem_psa3d':         mem_psa,
        'delta_psa3d':       round(aq_psa - mem_psa, 4),
        'psa3d_std':         round(float(np.std(psas_valid)), 4),
        'psa3d_spread':      round(float(psas_valid.max() - psas_valid.min()), 4),
        'bw_psa3d':          round(bw_psa, 4) if not np.isnan(bw_psa) else np.nan,
        'norm_delta_psa':    norm_delta_psa,   # ΔPSA / SASA_aq (Yu 2026)
        'tpsa_2d':           tpsa_2d,
        'ensemble_energy':   obj.get('ensembleenergy', np.nan),
        'pop_lowest_pct':    obj.get('poplowestpct', np.nan),
        'unique_confs':      obj.get('uniqueconfs', np.nan),
        'error':             None,
    }


def run(pickle_path: str, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / 'cremp_deltapsa.csv'

    tf = tarfile.open(pickle_path)
    members = [m for m in tf.getmembers() if m.name.endswith('.pickle')]
    print(f"Found {len(members)} pickle files in {pickle_path}")

    results = []
    failed = 0

    for member in tqdm(members, desc="CREMP ΔPSA"):
        try:
            f = tf.extractfile(member)
            obj = pickle.load(f)
            row = process_pickle(member.name, obj)
        except Exception as e:
            row = {'compound_id': Path(member.name).stem, 'error': str(e)}
            failed += 1
        results.append(row)

    tf.close()

    df = pd.DataFrame(results)
    n_ok = df['error'].isna().sum() if 'error' in df.columns else len(df)
    print(f"\nSucceeded: {n_ok} / {len(df)}")
    print(f"Failed:    {failed}")
    if failed > 0:
        print(df[df['error'].notna()][['compound_id', 'error']].to_string())

    df.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")

    # Quick summary
    ok = df[df['error'].isna()].copy() if 'error' in df.columns else df
    print(f"\nΔPSA summary (n={len(ok)}):")
    print(ok['delta_psa3d'].describe().round(2).to_string())


def parse_args():
    parser = argparse.ArgumentParser(description="Compute ΔPSA from CREMP pickle ensembles")
    parser.add_argument('--pickle', default='dependencies/pickle.tar.gz')
    parser.add_argument('--outdir', default='results')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    run(args.pickle, Path(args.outdir))

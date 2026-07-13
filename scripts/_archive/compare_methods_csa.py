# env: chameleon-calc
"""
compare_methods_csa.py
----------------------
Three-way CsA 3D-descriptor comparison: RDKit vacuum (feature_matrix.csv) vs
CREST V1 implicit-solvent ensemble vs experimental crystal structures.

Each method's "open" and "closed" states:
  - RDKit vacuum / CREST:  max-PSA conformer (open) and min-PSA conformer (closed)
  - Experimental:          A1 crystal (open, aqueous) and C1 crystal (closed)

Output: results/csa_threeway_descriptors.csv
See docs/experiments/2026-06-09_csa_threeway_method_comparison.md
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

def rg(c):
    c = c - c.mean(0)
    return float(np.sqrt((c**2).sum(1).mean()))

# --- RDKit vacuum (feature_matrix.csv, CsA = first row) ---
fm = pd.read_csv("results/archive/feature_matrix.csv", low_memory=False)
csa = fm.iloc[0]
rdkit_vac = {
    "open_psa": round(float(csa["aq_psa3d"]), 1),
    "closed_psa": round(float(csa["mem_psa3d"]), 1),
    "delta_psa": round(float(csa["delta_psa3d"]), 1),
    "open_hb": float(csa["aq_hb_count"]),
    "closed_hb": float(csa["mem_hb_count"]),
    "open_Rg": round(float(csa["aq_Rg"]), 2),
    "closed_Rg": round(float(csa["mem_Rg"]), 2),
    "boltzmann_psa": np.nan,
}

# --- CREST V1 implicit-solvent ensemble ---
with open("data/CREST_CsA_20260512/ensemble.json") as f:
    data = json.load(f)
w = np.array([c["boltzmannweight"] for c in data["conformers"]]); w /= w.sum()
psa = np.array([c["psa"] for c in data["conformers"]])
hb = np.array([c["hbonds"] for c in data["conformers"]])
sup = [m for m in Chem.SDMolSupplier("data/CREST_CsA_20260512/ensemble.sdf", removeHs=False) if m]
rgs = np.array([rg(m.GetConformer().GetPositions()) for m in sup[:len(w)]])
imax, imin = int(psa.argmax()), int(psa.argmin())
crest = {
    "open_psa": round(float(psa[imax]), 1),
    "closed_psa": round(float(psa[imin]), 1),
    "delta_psa": round(float(psa[imax] - psa[imin]), 1),
    "open_hb": int(hb[imax]),
    "closed_hb": int(hb[imin]),
    "open_Rg": round(float(rgs[imax]), 2),
    "closed_Rg": round(float(rgs[imin]), 2),
    "boltzmann_psa": round(float(np.dot(w, psa)), 1),
}

# --- Experimental crystal structures (computed 2026-06-05, rdFreeSASA) ---
exp = {
    "open_psa": 137.5,    # A1 aqueous (X-ray, CCDC 2149649)
    "closed_psa": 95.9,   # C1 closed (DEKSAN, CCDC 1138505)
    "delta_psa": round(137.5 - 95.9, 1),
    "open_hb": 2, "closed_hb": 4,
    "open_Rg": 6.15, "closed_Rg": 6.42,
    "boltzmann_psa": np.nan,
}

df = pd.DataFrame({"RDKit_vacuum": rdkit_vac, "CREST_V1_implicit": crest, "Experimental_crystal": exp})
df = df.loc[["open_psa", "closed_psa", "delta_psa", "boltzmann_psa",
             "open_hb", "closed_hb", "open_Rg", "closed_Rg"]]
out = Path("results/csa_threeway_descriptors.csv")
df.to_csv(out)
print(df.to_string())
print(f"\nSaved -> {out}")
print("open = max-PSA conf (RDKit/CREST) or A1 crystal (exp); closed = min-PSA conf or C1 crystal")

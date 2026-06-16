# env: MACE
"""
mace_tip3p_feasibility.py
--------------------------
Feasibility test: MACE-OFF23(M) + explicit TIP3P water MD on CsA C1 conformer.

Hybrid ML/MM setup:
  - MACE-OFF23(M)  : CsA intramolecular forces (bonds, angles, torsions, vdW)
  - GAFF2          : CsA partial charges + LJ for solute-water cross interactions
  - TIP3P          : explicit water (water-water + water-solute non-bonded)
  - OpenMM         : MD engine, PME electrostatics, Langevin thermostat, MC barostat

Steps:
  1. Load C1 conformer from ensemble.sdf (A-like, 46.3% Boltzmann weight in xTB water)
  2. Assign GAFF2 parameters + partial charges to CsA
  3. Solvate in TIP3P water box
  4. Create mixed ML/MM system (MACE solute + TIP3P water)
  5. Energy minimize
  6. NVT equilibration at 300 K
  7. NPT production at 300 K, 1 atm
  8. Report step rate -> project cost to 50 ns x N solvents x library scale

Key outputs:
  - Console: step rate, density, projected cost breakdown
  - results/mace_tip3p_feasibility.csv  (step, energy, density, Rg)

Required setup (one-time, in MACE env):
  pip install openmm openmmml openmmforcefields "openff-toolkit>=0.16"

Usage:
  python scripts/mace_tip3p_feasibility.py
  python scripts/mace_tip3p_feasibility.py --model models/MACE-OFF23_medium.model
  python scripts/mace_tip3p_feasibility.py --steps 100000 --padding 10
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np

# ── Args ───────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--model", default="mace-off23-medium",
                    help="openmmml model name or path to local .model file "
                         "(default: mace-off23-medium)")
parser.add_argument("--steps", type=int, default=50_000,
                    help="NPT production steps at 2 fs (default 50000 = 100 ps)")
parser.add_argument("--equil-steps", type=int, default=10_000,
                    help="NVT equilibration steps (default 10000 = 20 ps)")
parser.add_argument("--padding", type=float, default=12.0,
                    help="Water box padding in Angstrom (default 12)")
args = parser.parse_args()

# ── Package checks ─────────────────────────────────────────────────────────────
missing = []
for pkg, install in [("openmm", "openmm"), ("openmmml", "openmmml"),
                     ("openmmforcefields", "openmmforcefields"),
                     ("openff.toolkit", "openff-toolkit>=0.16"),
                     ("rdkit", "rdkit")]:
    try:
        __import__(pkg.replace(".", "/").replace("/", "."))
    except ImportError:
        missing.append(f"  pip install {install}")
if missing:
    sys.exit("Missing packages — run:\n" + "\n".join(missing))

import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmmml import MLPotential
from openff.toolkit import Molecule as OFFMolecule, Topology as OFFTopology
from openmmforcefields.generators import GAFFTemplateGenerator
from rdkit import Chem

# ── Find C1 conformer SDF ──────────────────────────────────────────────────────
def _find_sdf() -> Path:
    candidates = [
        *sorted(Path("results/runs").glob("*_CsA/water/ensemble.sdf"), reverse=True),
        Path("data/CREST_CsA_20260512/ensemble.sdf"),
    ]
    for p in candidates:
        if p.exists():
            return p
    sys.exit("ERROR: ensemble.sdf not found. Expected in results/runs/*_CsA/water/")

sdf_path = _find_sdf()
print(f"\nLoading CsA C1 conformer: {sdf_path}")

# ── Load first conformer (C1 = highest Boltzmann weight) ──────────────────────
suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
rdmol = next((m for m in suppl if m is not None), None)
if rdmol is None:
    sys.exit("ERROR: Could not parse any conformer from SDF.")

n_solute = rdmol.GetNumAtoms()
print(f"  Atoms in CsA: {n_solute}")

# Positions: RDKit stores in Angstrom
conf = rdmol.GetConformer(0)
solute_pos_ang = np.array(
    [[conf.GetAtomPosition(i).x,
      conf.GetAtomPosition(i).y,
      conf.GetAtomPosition(i).z] for i in range(n_solute)]
)

# ── OpenFF molecule + partial charges ─────────────────────────────────────────
print("Assigning partial charges ...")
off_mol = OFFMolecule.from_rdkit(rdmol, allow_undefined_stereo=True)

# Charge method priority:
#   1. NAGL (neural net, reproduces AM1-BCC accuracy, no atom size limit, fast)
#   2. AM1-BCC (semi-empirical QM, standard for GAFF; slow/unreliable >150 atoms)
#   3. Gasteiger (empirical fallback, fast, less accurate for solute-water interactions)
# For production conformer sampling in water, NAGL or AM1-BCC is required.
# Gasteiger charges will underestimate solute-water electrostatics.
def _get_nagl_methods() -> list[str]:
    try:
        import openff.nagl_models as nm
        models = [str(m) for m in nm.list_available_nagl_models()
                  if "am1bcc" in str(m).lower()]
        def _rank(p):
            return 2 if "alpha" in p else 1 if "rc" in p else 0
        return sorted(models, key=_rank)
    except Exception:
        return []

charge_assigned = False
nagl_methods = _get_nagl_methods()
methods = nagl_methods + ["am1bcc", "gasteiger"]

for method in methods:
    try:
        off_mol.assign_partial_charges(method)
        print(f"  Charge method: {method}")
        if method == "gasteiger":
            print("  WARNING: Gasteiger charges will underestimate solute-water "
                  "electrostatics. Install openff-nagl or use AM1-BCC for production.")
        charge_assigned = True
        break
    except Exception as e:
        print(f"  {method} failed: {type(e).__name__}")

if not charge_assigned:
    sys.exit("ERROR: Could not assign partial charges. "
             "Install AmberTools (conda install -c conda-forge ambertools) "
             "or OpenEye.")

# ── GAFF2 template generator ──────────────────────────────────────────────────
gaff = GAFFTemplateGenerator(molecules=[off_mol], forcefield="gaff-2.11")

# ── OpenMM force field ────────────────────────────────────────────────────────
ff = app.ForceField("tip3p.xml")
ff.registerTemplateGenerator(gaff.generator)

# ── Build OpenMM topology from OpenFF ─────────────────────────────────────────
off_top = OFFTopology.from_molecules([off_mol])
omm_top = off_top.to_openmm()
positions_nm = solute_pos_ang * 0.1  # Å -> nm

# ── Solvate ───────────────────────────────────────────────────────────────────
pad = args.padding * unit.angstrom
print(f"\nSolvating with TIP3P ({args.padding:.0f} Å padding) ...")
modeller = app.Modeller(omm_top, positions_nm * unit.nanometer)
modeller.addSolvent(ff, model="tip3p", padding=pad)
n_water = (modeller.topology.getNumAtoms() - n_solute) // 3
n_total = modeller.topology.getNumAtoms()
print(f"  Water molecules: {n_water}")
print(f"  Total atoms:     {n_total}")

# ── MM system ─────────────────────────────────────────────────────────────────
mm_system = ff.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=10.0 * unit.angstrom,
    constraints=app.HBonds,
)

# ── MACE ML potential ─────────────────────────────────────────────────────────
print(f"\nLoading MACE potential: {args.model}")
model_path = Path(args.model)
if model_path.exists():
    # Local model file — try mace backend with path
    try:
        potential = MLPotential("mace", modelPath=str(model_path))
    except TypeError:
        potential = MLPotential("mace", model_path=str(model_path))
else:
    potential = MLPotential(args.model)

# Create hybrid ML(solute) / MM(water) system
solute_atoms = list(range(n_solute))
print("Building mixed ML/MM system ...")
print("  (MACE handles CsA intramolecular; GAFF2+TIP3P handles solute-water)")
ml_system = potential.createMixedSystem(
    modeller.topology, mm_system, solute_atoms, interpolate=False
)
print("  Mixed system OK")

# ── Integrator + CUDA platform ────────────────────────────────────────────────
integrator = mm.LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picosecond,
    0.002 * unit.picoseconds,
)

try:
    platform = mm.Platform.getPlatformByName("CUDA")
    props = {"CudaPrecision": "mixed"}
    print("  Platform: CUDA")
except Exception:
    platform = mm.Platform.getPlatformByName("CPU")
    props = {}
    print("  WARNING: CUDA not available, falling back to CPU")

simulation = app.Simulation(
    modeller.topology, ml_system, integrator, platform, props
)
simulation.context.setPositions(modeller.positions)

# ── Energy minimize ────────────────────────────────────────────────────────────
print("\n[1/3] Energy minimization ...")
t0 = time.perf_counter()
simulation.minimizeEnergy(maxIterations=2000)
state = simulation.context.getState(getEnergy=True)
e_min = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
print(f"  Minimized energy : {e_min:,.1f} kcal/mol  ({time.perf_counter()-t0:.1f} s)")

# ── NVT equilibration ─────────────────────────────────────────────────────────
print(f"\n[2/3] NVT equilibration ({args.equil_steps} steps = "
      f"{args.equil_steps * 0.002:.1f} ps) ...")
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
t0 = time.perf_counter()
simulation.step(args.equil_steps)
equil_elapsed = time.perf_counter() - t0
equil_rate = args.equil_steps / equil_elapsed
state = simulation.context.getState(getEnergy=True)
e_nvt = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
print(f"  Step rate        : {equil_rate:.1f} steps/s")
print(f"  Energy after NVT : {e_nvt:,.1f} kcal/mol")

# ── NPT production ─────────────────────────────────────────────────────────────
ml_system.addForce(
    mm.MonteCarloBarostat(1.0 * unit.atmospheres, 300 * unit.kelvin, 25)
)
simulation.context.reinitialize(preserveState=True)

print(f"\n[3/3] NPT production ({args.steps} steps = {args.steps * 0.002:.1f} ps) ...")
report_every = max(1, args.steps // 20)
records = []
t0 = time.perf_counter()

for i in range(0, args.steps, report_every):
    simulation.step(report_every)
    state = simulation.context.getState(
        getEnergy=True, getPositions=True, enforcePeriodicBox=True
    )
    E = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
    vol_nm3 = float(np.linalg.det(box))

    # Approximate density: use water mass only (dominant contribution)
    mass_g = n_water * 18.015 / 6.02214e23
    density = mass_g / (vol_nm3 * 1e-21)

    # Radius of gyration of solute (conformational indicator)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    sol_pos = pos[:n_solute]
    com = sol_pos.mean(axis=0)
    rg = float(np.sqrt(((sol_pos - com) ** 2).sum(axis=1).mean()))

    step = i + report_every
    records.append((step, E, density, rg))
    print(f"  Step {step:7d}: E = {E:12,.1f} kcal/mol  "
          f"rho = {density:.3f} g/cm3  Rg = {rg:.3f} nm", flush=True)

prod_elapsed = time.perf_counter() - t0
prod_rate = args.steps / prod_elapsed

# ── Save CSV ───────────────────────────────────────────────────────────────────
out_path = Path("results/mace_tip3p_feasibility.csv")
out_path.parent.mkdir(exist_ok=True)
with open(out_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["step", "energy_kcal_mol", "density_g_cm3", "rg_nm"])
    w.writerows(records)
print(f"\n  CSV -> {out_path}")

# ── Cost projection ────────────────────────────────────────────────────────────
STEPS_PER_NS   = 500_000   # at 2 fs timestep
TARGET_NS      = 50.0
N_SOLVENTS     = 2         # water + CHCl3
N_STARTS       = 4         # multi-start cis/trans seeds (parallel per compound)
N_COMPOUNDS    = 50

hours_per_run  = (TARGET_NS * STEPS_PER_NS) / prod_rate / 3600
hours_per_cpd  = hours_per_run * N_SOLVENTS * N_STARTS   # if runs sequential
hours_parallel = hours_per_run                             # if N_SOLVENTS*N_STARTS GPUs per cpd
hours_library_seq  = hours_per_cpd  * N_COMPOUNDS
hours_library_par  = hours_parallel * N_COMPOUNDS          # best case: all parallel

avg_density = float(np.mean([r[2] for r in records]))
density_ok  = 0.95 < avg_density < 1.05

print("\n" + "=" * 68)
print("  Results summary")
print("=" * 68)
print(f"  System          : {n_solute} solute atoms + {n_water} water molecules "
      f"({n_total} total)")
print(f"  Equil rate      : {equil_rate:8.1f} steps/s")
print(f"  Production rate : {prod_rate:8.1f} steps/s")
print(f"  Avg NPT density : {avg_density:.3f} g/cm3  "
      f"({'OK' if density_ok else 'WARNING: off from 1.00'})")
print()
print(f"  -- Projected cost ({TARGET_NS:.0f} ns production at 2 fs timestep) --")
print(f"  Per run (1 start, 1 solvent) : {hours_per_run:.1f} h")
print(f"  Per compound ({N_SOLVENTS} solvents, {N_STARTS} starts, serial) : "
      f"{hours_per_cpd:.1f} h  ({hours_per_cpd/24:.1f} days)")
print(f"  Per compound ({N_SOLVENTS*N_STARTS} GPUs in parallel)           : "
      f"{hours_parallel:.1f} h")
print(f"  Full library ({N_COMPOUNDS} compounds, serial per cpd)   : "
      f"{hours_library_seq:.0f} h  ({hours_library_seq/24:.0f} days, {N_SOLVENTS*N_STARTS} GPUs each)")
print()

if prod_rate >= 100:
    verdict = "FEASIBLE — explicit TIP3P + MACE viable at library scale with GPU cluster."
    detail  = f"  {hours_per_cpd:.1f} h/compound with {N_SOLVENTS*N_STARTS} GPUs in parallel."
elif prod_rate >= 20:
    verdict = "MARGINAL — feasible for small pilot libraries (<15 compounds)."
    detail  = ("  Consider shorter sampling (10-20 ns) or GBn2 implicit solvent "
               "as a faster approximation.")
else:
    verdict = "SLOW — explicit TIP3P + MACE not practical at library scale."
    detail  = ("  Options: (a) GBn2 implicit solvent benchmark, "
               "(b) GFN-FF for sampling + MACE for re-ranking, "
               "(c) reduce to ~5 ns and accept reduced sampling quality.")

print(f"  Verdict: {verdict}")
print(f"  {detail}")
print("=" * 68)

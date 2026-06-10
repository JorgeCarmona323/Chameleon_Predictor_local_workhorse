# isomer_figures.pml
# PyMOL script for the DOPC 3-12-8-12 R/S isomer story.
# Run from repo root:  pymol scripts/isomer_figures.pml
# Requires the PDBs from scripts/make_isomer_figures.py (results/figures/isomers/).
#
# Story: R opens in water (few intramolecular H-bonds) but closes in membrane;
#        S is pre-organized -- closed (H-bonded) in both environments.

reinitialize
bg_color white
set ray_shadows, 0
set ray_opaque_background, 0
set cartoon_transparency, 0.0
set dash_color, gold
set dash_gap, 0.35
set dash_width, 3.0
set label_size, 18

# --- load the four dominant conformers ---
load results/figures/isomers/R_water_dominant.pdb, R_water
load results/figures/isomers/R_mem_dominant.pdb,   R_mem
load results/figures/isomers/S_water_dominant.pdb, S_water
load results/figures/isomers/S_mem_dominant.pdb,   S_mem

# --- common style: licorice backbone, thin, polar atoms colored ---
hide everything
show sticks
set stick_radius, 0.13
util.cbaw            # color by atom, carbons white
color salmon, R_water and elem C
color salmon, R_mem and elem C
color skyblue, S_water and elem C
color skyblue, S_mem and elem C

# --- intramolecular H-bonds (polar contacts, needs the explicit H in the PDB) ---
python
for obj in ["R_water", "R_mem", "S_water", "S_mem"]:
    cmd.distance(f"hb_{obj}", f"{obj} and (elem N+O)", f"{obj} and (elem N+O)", 3.5, mode=2)
    cmd.hide("labels", f"hb_{obj}")
python end

# --- grid view: all four side by side ---
set grid_mode, 1
set grid_slot, 1, R_water
set grid_slot, 2, S_water
set grid_slot, 3, R_mem
set grid_slot, 4, S_mem

orient
zoom all, 3

# --- render the 2x2 comparison (top: water, bottom: membrane) ---
set ray_trace_mode, 1
ray 1600, 1600
png results/figures/isomers/pymol_2x2_RS_water_mem.png, dpi=300

# --- headline single panel: R_water vs S_water (the H-bond difference) ---
set grid_mode, 0
disable all
enable R_water
enable S_water
# overlay them aligned for direct shape comparison
align S_water, R_water
orient R_water or S_water
zoom (R_water or S_water), 2
ray 1600, 1200
png results/figures/isomers/pymol_R_vs_S_water_overlay.png, dpi=300

# Notes printed to PyMOL log
print "Rendered: pymol_2x2_RS_water_mem.png (R/S x water/mem, gold dashes = intramolecular H-bonds)"
print "Rendered: pymol_R_vs_S_water_overlay.png (R salmon vs S skyblue, aligned)"
print "Story: R_water has fewer H-bonds (open); S_water and all membrane forms are H-bonded (closed)."

# visualize_csa_vmd.tcl
# ----------------------
# VMD visualization for the CsA water conformer ensemble.
# Highlights the A1 fingerprint: cis-amide region and NH H-bond donors.
#
# Reference: Limbach et al., JACS 2022, 144, 12602
#
# HOW TO USE:
#   Step 1) In VMD: File > New Molecule
#           Load:  data/CREST_CsA_20260512/ensemble.xyz   (type: xyz)
#           This creates mol 0 with 23 frames (one per conformer).
#
#   Step 2) Open VMD Tk Console (Extensions > Tk Console) and run:
#           source scripts/visualize_csa_vmd.tcl
#
#   Step 3) Optional — add A1 crystal reference:
#           File > New Molecule
#           Load:  data/experimental_structure_references_CsA/CsA_A1_xray_CCDC2149649.cif
#           Then align with the ensemble backbone (see bottom of this file).
#
# KEY: Run validate_csa_water.py first to know which frames are A1-like.
#      Frame 0 = conformer 1 (highest Boltzmann weight, likely A1-like).

# ── Display settings ──────────────────────────────────────────────────────────
display projection Orthographic
display depthcue off
color Display Background white
axes location off

# ── Ensemble representations (mol 0) ─────────────────────────────────────────
set ens 0

# Clear default representation
mol delrep 0 $ens

# Rep 0: Backbone tube — all heavy atoms, thin, gray
mol addrep $ens
mol modstyle   0 $ens Tube 0.15 12
mol modcolor   0 $ens ColorID 2
mol modselect  0 $ens "noh"
mol modmaterial 0 $ens Transparent

# Rep 1: Ring nitrogens — large VDW, colored by type
#   Free NH (Abu2, Val5, Ala7, D-Ala8): blue
#   N-methyl: cyan
mol addrep $ens
mol modstyle   1 $ens VDW 0.5 12
mol modcolor   1 $ens ColorID 0
mol modselect  1 $ens "name N and within 1.15 of name H"
mol modmaterial 1 $ens Glossy

# Rep 2: Carbonyl oxygens — red VDW (H-bond acceptors)
mol addrep $ens
mol modstyle   2 $ens VDW 0.4 12
mol modcolor   2 $ens ColorID 1
mol modselect  2 $ens "name O"
mol modmaterial 2 $ens Glossy

# Rep 3: All bonds as thin licorice for context
mol addrep $ens
mol modstyle   3 $ens Licorice 0.08 10 10
mol modcolor   3 $ens ColorID 8
mol modselect  3 $ens "noh"
mol modmaterial 3 $ens Transparent

# ── Animation ─────────────────────────────────────────────────────────────────
animate speed 0.15
animate style loop

# ── Console output ────────────────────────────────────────────────────────────
puts ""
puts "=== CsA Water Ensemble | Limbach JACS 2022 A1 Validation ==="
puts ""
puts "Loaded ensemble: mol $ens  (23 frames)"
puts "  Frame 0  = Conformer 1  (highest Boltzmann weight, ~46%)"
puts "  Frame 22 = Conformer 23 (lowest weight)"
puts ""
puts "Color scheme:"
puts "  Blue  : Free NH donors  (Abu2, Val5, Ala7, D-Ala8)"
puts "  Red   : Carbonyl oxygens (H-bond acceptors)"
puts "  Gray  : Backbone"
puts ""
puts "A1 fingerprint to look for (run validate_csa_water.py for quantitative results):"
puts "  - Cis amide at MeVal11-MeBmt1: the N-methyl amide where the"
puts "    two adjacent Ca-H protons point toward each other"
puts "  - Abu2 (blue, residue 2) NH close to a red C=O on the far side of the ring"
puts "  - Ala7 (blue, residue 7) NH close to a red C=O"
puts "  - Val5 (blue, residue 5) NH should be open / far from acceptors"
puts ""
puts "To show a single conformer (e.g. frame 0):"
puts "  animate goto 0"
puts "  animate pause"
puts ""
puts "To measure a dihedral in VMD:"
puts "  Mouse > Label > Dihedral, then click 4 atoms"
puts ""

# ── Crystal structure alignment (uncomment after loading CIF as mol 1) ────────
# After loading CsA_A1_xray_CCDC2149649.cif as mol 1, run these lines:
#
# set ref  [atomselect 1 "noh"]
# set mob  [atomselect 0 "noh" frame 0]
# set M    [measure fit $mob $ref]
# set all  [atomselect 0 "all"]
# $all move $M
# puts "Ensemble aligned to A1 crystal structure (CCDC 2149649)"

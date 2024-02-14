# define input structure
structure="1v4s_Nmut_capped_A173F"

# generate topology from pdb file without HETATMs
$gmx pdb2gmx -f ../_input_structures/$structure.pdb -o processed_original.gro -ff 'a99SBdisp'<<EOF
1 
EOF

# rename topology file 
mv topol.top topol_original.top

# add Na coordinates from WT structure file using add_crystal_ion.py
# the script adds ion coordinates to processed.gro, changes the total number of
# atoms in the system at the top of the .gro file, and renumbers the atom number of
# the ion. Also, NA is added to topol.top under [ molecules ].
# remember to change Na coordinates for 1v4t structures
python add_crystal_ion.py

# define the simulation box (with 1.1 nm to box boundaries), center protein in box
$gmx editconf -f processed.gro -o processed_box.gro -c -d 1.1 -bt dodecahedron

# solvate with disp water
$gmx solvate -cp processed_box.gro -cs a99SBdisp.ff/a99SBdisp_water -o processed_box_solv.gro -p topol.top

# generate .tpr file to be used as input for genion
$gmx grompp -f ../../mdp_files/min.mdp -c processed_box_solv.gro -p topol.top -o ions.tpr -maxwarn 1

# add ions
$gmx genion -s ions.tpr -o processed_box_solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15<<EOF
15
EOF

# generate ca position restraint files
$gmx make_ndx -f processed.gro -o ca_index.ndx<<EOF
keep 3
q
EOF
$gmx genrestr -f processed.gro -n ca_index.ndx -o ca_posre.itp
#!/bin/bash

# set path to mdp files
mdp_path="../../mdp_files"

# generate tpr file for run
gmx_mpi grompp -f $mdp_path/em.mdp -c ../../md_input_files/1v4s_A173F_setup/processed_box_solv_ions.gro -p ../../md_input_files/1v4s_A173F_setup/topol.top -o em.tpr 

# run minimization
gmx_mpi mdrun -s em.tpr -v -deffnm em

# plot potential energy as function of minimization step
gmx_mpi energy -f em.edr -o em.xvg <<EOF
10 0
EOF
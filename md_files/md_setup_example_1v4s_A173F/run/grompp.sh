#!/bin/bash

# define paths to inputs
input_structure_path="../../../md_input_files/1v4s_A173F_setup"
mdp_path="../../../mdp_files"
eq_out_path="../../eq"

# prepare input for mdrun
for run in run1 run2 run3
do
cd $run
gmx_mpi grompp -f $mdp_path/run.mdp -c $eq_out_path/$run/npt2.gro -t $eq_out_path/$run/npt2.cpt -p $input_structure_path/topol.top -o run.tpr
cd ..
done
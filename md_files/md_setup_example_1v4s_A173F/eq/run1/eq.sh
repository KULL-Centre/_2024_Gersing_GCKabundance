#!/bin/bash

# define paths to inputs
input_structure_path="../../../md_input_files/1v4s_A173F_setup"
mdp_path="../../../mdp_files"

# generate new input from energy minimized structure
# apply position restraints to all protein heavy atoms with position restraint origins defined by em.gro coordinates
gmx_mpi grompp -f $mdp_path/eq_NVT1.mdp -c ../../em/em.gro -r ../../em/em.gro -p $input_structure_path/topol.top -o nvt1.tpr
gmx_mpi mdrun -s nvt1.tpr -deffnm nvt1 -nb gpu -pme gpu -npme 0 -notunepme -pin on -ntomp 18 -maxh 23.9 -v

# generate new input from first NVT run, now with position restraints only on the CA atoms
gmx_mpi grompp -f $mdp_path/eq_NVT2.mdp -c nvt1.gro -r nvt1.gro -t nvt1.cpt -p $input_structure_path/topol.top -o nvt2.tpr
gmx_mpi mdrun -s nvt2.tpr -deffnm nvt2 -nb gpu -pme gpu -npme 0 -notunepme -pin on -ntomp 18 -maxh 23.9 -v

# plot T and P for NVT runs
for traj in nvt1 nvt2
do
# plot T
gmx_mpi energy -f $traj.edr -o plot/$traj-temp.xvg<<EOF
16 0
EOF
# plot P
gmx_mpi energy -f $traj.edr -o plot/$traj-pressure.xvg<<EOF
18 0
EOF
done

# add pressure equilibration, now without any position restraints
gmx_mpi grompp -f $mdp_path/eq_NPT1.mdp -c nvt2.gro -t nvt2.cpt -p $input_structure_path/topol.top -o npt1.tpr
gmx_mpi mdrun -s npt1.tpr -deffnm npt1 -nb gpu -pme gpu -npme 0 -notunepme -pin on -ntomp 18 -maxh 23.9 -v

# add pressure equilibration, now without any position restraints
gmx_mpi grompp -f $mdp_path/eq_NPT2.mdp -c npt1.gro -t npt1.cpt -p $input_structure_path/topol.top -o npt2.tpr
gmx_mpi mdrun -s npt2.tpr -deffnm npt2 -nb gpu -pme gpu -npme 0 -notunepme -pin on -ntomp 18 -maxh 23.9 -v

# plot T and P for NPT runs
for traj in npt1 npt2
do
# plot T
gmx_mpi energy -f $traj.edr -o plot/$traj-temp.xvg<<EOF
15 0
EOF
# plot P
gmx_mpi energy -f $traj.edr -o plot/$traj-pressure.xvg<<EOF
17 0
EOF
done

# convert trajectories for visualisation 
for traj in nvt1 nvt2 npt1 npt2
do 
gmx_mpi trjconv -f $traj.xtc -s $traj.tpr -o $traj\_trjconv.xtc -pbc mol -ur compact<<EOF
0
EOF
done
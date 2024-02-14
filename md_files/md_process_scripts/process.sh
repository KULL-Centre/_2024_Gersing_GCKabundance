#!/bin/bash

for run in run1 run2 run3
do 

mkdir $run/plot

# plot thermodynamic properties from .edr files

gmx_mpi energy -f $run/run.edr -o $run/plot/potential.xvg<<EOF
11 0 
EOF

gmx_mpi energy -f $run/run.edr -o $run/plot/kinetic.xvg<<EOF
12 0
EOF

gmx_mpi energy -f $run/run.edr -o $run/plot/totalene.xvg<<EOF
13 0
EOF

gmx_mpi energy -f $run/run.edr -o $run/plot/temp.xvg<<EOF
15 0
EOF

gmx_mpi energy -f $run/run.edr -o $run/plot/pressure.xvg<<EOF
17 0
EOF

gmx_mpi energy -f $run/run.edr -o $run/plot/volume.xvg<<EOF
22 0
EOF

gmx_mpi energy -f $run/run.edr -o $run/plot/density.xvg<<EOF
23 0 
EOF

mkdir $run/processed_traj

# center protein and make molecules whole
# write coordinates for entire system to be able to keep track of ions
gmx_mpi trjconv -f $run/run.xtc -s $run/run.tpr -o $run/processed_traj/run_center.xtc -pbc mol -center -ur compact<<EOF
1 0
EOF

gmx_mpi trjconv -f $run/run.xtc -s $run/run.tpr -dump 0 -o $run/processed_traj/run_center.gro -pbc mol -center -ur compact<<EOF
1 0
EOF

# center protein and make molecules whole
# write coordinates for protein only
gmx_mpi trjconv -f $run/run.xtc -s $run/run.tpr -o $run/processed_traj/run_center_prot.xtc -pbc mol -center -ur compact<<EOF
1 1
EOF

gmx_mpi trjconv -f $run/run.xtc -s $run/run.tpr -dump 0 -o $run/processed_traj/run_center_prot.gro -pbc mol -center -ur compact<<EOF
1 1
EOF

# calc min dist to periodic image 
gmx_mpi mindist -f $run/run.xtc -s $run/run.tpr -pbc yes -od $run/plot/min_periodic_dist.xvg -pi<<EOF
1
EOF

done
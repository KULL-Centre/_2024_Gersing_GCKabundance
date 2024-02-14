#!/bin/bash

# run simulation
gmx_mpi mdrun -s run.tpr -deffnm run -cpi run.cpt -nb gpu -pme gpu -npme 0 -notunepme -pin on -ntomp 18 -maxh 23.9 -v
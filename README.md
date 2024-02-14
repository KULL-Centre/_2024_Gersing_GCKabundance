# _2024_Gersing_GCKabundance

This repository contains the Jupyter Notebooks and input data to generate the plots found in "Characterizing glucokinase variant mechanisms using a multiplexed abundance assay" by Sarah Gersing, Thea K Schulze, Matteo Cagiada, Amelie Stein, Frederick P Roth, Kresten Lindorff-Larsen and Rasmus Hartmann-Petersen.

The preprinted paper can be found at https://doi.org/10.1101/2023.05.24.542036

Simulation data can be found at https://doi.org/10.17894/ucph.4e38c597-fae2-4654-a48a-8e714263c5a1

Content of directories in repository:

`md_files/mdp_files/` contains mdp files used to run the MD simulations. 

`md_files/input_files/` contains MD simulation starting structures and scripts for generating system starting configurations for simulations based on these structures.

`md_files/process_scripts/` contains the script used to processed raw trajectories prior to other analyses. 

`md_files/md_setup_example_1v4s_A173F/` contains example scripts for launching energy minimisation, equilibration and production runs.

`md_files/analysis_scripts/` contains scripts used to calculate RMSDs, RMSFs and other trajectory observables reported in the paper. 

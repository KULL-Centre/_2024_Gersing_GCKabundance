import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mdtraj as md
import itertools

def calc_angles_dists(variant, conf, run, data_dir):
    """
    
    Calculate cleft/opening angles between small and large domain for specified 
    GCK variant, and calculate distance between residues 159 and 452.
    Observables are calculated for each frame in input trajectories. 

    Angle definition 1 and the 159-452 distance were
    previously introduced for analysis of GCK simulations 
    in the following work:
    "Zhang J, Li C, Chen K, Zhu W, Shen X, Jiang H. 
    Conformational transition pathway in the allosteric 
    process of human glucokinase. Proceedings of the 
    National Academy of Sciences. 2006 Sep; 103(36):13368–13373."
    https://doi.org/10.1073/pnas.0605738103
    
    """
    # load trajectory
    traj = md.load(f"{data_dir}/{variant}/{conf}/run/{run}/run_center_prot.xtc",
                   top=f"{data_dir}/{variant}/{conf}/run/{run}/run_center_prot.gro")

    # define time for each frame in ps and ns
    # (we save a frame every 100th ps)
    t_ps = np.arange(0,traj.n_frames*100,100)
    t_ns = t_ps / 1000

    # get atom numbers of CA atoms to be used for angle calculation
    angle_ca_1 = traj.top.select('(residue == 229) and (name == CA)')[0]
    angle_ca_2 = traj.top.select('(residue == 233) and (name == CA)')[0]
    angle_ca_3_1 = traj.top.select('(residue == 169) and (name == CA)')[0] 
    angle_ca_3_2 = traj.top.select('(residue == 109) and (name == CA)')[0] 

    # define angles
    angle_definition_3_1 = np.array([[angle_ca_1, angle_ca_2, angle_ca_3_1]])
    angle_definition_3_2 = np.array([[angle_ca_1, angle_ca_2, angle_ca_3_2]])

    # calculate angle in each frame using angle definition 1
    angles_traj_1 = md.compute_angles(traj, angle_indices=angle_definition_3_1, periodic=False) 
    angles_traj_1 = angles_traj_1 * (180/np.pi) # from radians to degrees

    # calculate angle in each frame using angle definition 2
    angles_traj_2 = md.compute_angles(traj, angle_indices=angle_definition_3_2, periodic=False)
    angles_traj_2 = angles_traj_2 * (180/np.pi) # from radians to degrees

    # get atom numbers for CA atoms to be used for distance calculation
    dist_ca_1 = traj.top.select('(residue == 159) and (name == CA)')[0]
    dist_ca_2 = traj.top.select('(residue == 452) and (name == CA)')[0]

    # calculate distance in each frame
    dist = md.compute_distances(traj, atom_pairs=np.array([[dist_ca_1,dist_ca_2]]), periodic=False)

    # save to txt files
    np.savetxt(f'../analysis_output/angles+dists/{variant}_{conf}_{run}_angle1.txt', np.column_stack([t_ns,angles_traj_1.flatten()]))
    np.savetxt(f'../analysis_output/angles+dists/{variant}_{conf}_{run}_angle2.txt', np.column_stack([t_ns,angles_traj_2.flatten()]))
    np.savetxt(f'../analysis_output/angles+dists/{variant}_{conf}_{run}_dist.txt', np.column_stack([t_ns,dist.flatten()]))
    
if __name__ == "__main__":
    
    # give name of directory with traj data
    #data_dir = 

    variants = ["gck_wt","gck_V455M","gck_G175E","gck_D158A","gck_A173F","gck_G162Q"]
    conformations = ["closed_conf","superopen_conf"]
    runs = ["run1","run2","run3"]
    
    # loop over all variants, conformations and runs to calculate
    # angles and distance for every individual trajectory
    for variant, conf, run in itertools.product(variants,conformations, runs):
        print(variant, conf)
        calc_angles_dists(variant, conf, run, data_dir)
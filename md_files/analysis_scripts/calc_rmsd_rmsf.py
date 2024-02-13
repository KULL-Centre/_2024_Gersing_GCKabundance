import numpy as np
import pandas as pd
import mdtraj as md
import itertools
from Bio.PDB.DSSP import make_dssp_dict

def define_domains():
    """
    Define residues in large domain (ld), small domain (sd) and loop region (loop)
    All numbers are residue numbers, not indices.
    """
    # for super-open structure (excluding hinge regions)
    loop_res_arr_open = np.arange(151,179+1)
    ld_res_arr_open = np.concatenate((np.arange(15,45+1), np.arange(205,440+1)))
    sd_loop_res_arr_open = np.concatenate((np.arange(74,201+1),np.arange(448,459+1))) 
    sd_res_arr_open = np.setdiff1d(sd_loop_res_arr_open,loop_res_arr_open)
    domain_open_resi_dict = {"large_domain":ld_res_arr_open, "small_domain":sd_res_arr_open, "loop":loop_res_arr_open}
    
    # for the closed structure (excluding hinge regions)
    loop_res_arr_closed = np.arange(151,179+1)
    ld_res_arr_closed = np.concatenate((np.arange(14,45+1), np.arange(205,440+1)))
    sd_loop_res_arr_closed = np.concatenate((np.arange(72,203+1),np.arange(444,459+1))) 
    sd_res_arr_closed = np.setdiff1d(sd_loop_res_arr_closed,loop_res_arr_closed)
    domain_closed_resi_dict = {"large_domain":ld_res_arr_closed, "small_domain":sd_res_arr_closed, "loop":loop_res_arr_closed}

    define_domains_dict = {"closed_conf":domain_closed_resi_dict,"superopen_conf":domain_open_resi_dict}
    
    return define_domains_dict

def calc_rmsd(target_traj, ref_traj, atom_idx):
    """
    Calculate rmsd of frames in target_traj with respect to first 
    frame of ref_traj using the mdtraj function rmsd. 
    Rmsd is calculated only for the part of the structure 
    that is defined by atom indices in atom_idx. 
    atom_idx is the atom_idx used for both target and ref trajs.
    """
    rmsd = md.rmsd(target_traj, ref_traj, frame=0, atom_indices=atom_idx, 
                   ref_atom_indices=atom_idx, precentered=False)
    return rmsd

def get_atom_idx_domain(traj, domain_resi_dict, atom_group="all"):
    """
    Create dictionary with atom indices for 
    residues in the different GCK domains defined 
    in domain_res_dict. 
    """    
    # get topology information in dataframe format
    top_df, _ = traj.top.to_dataframe()
    
    # slice topology to only contain carbon alpha atoms
    if atom_group == "carbon-alpha":
        top_df = top_df.iloc[traj.top.select("protein and backbone and name CA")]
    
    # slice topology to only contain backbone heavy atoms
    if atom_group == "backbone-heavy":
        top_df = top_df.iloc[traj.top.select("protein and backbone and not name H")]
        
    # extract residue information from topology dataframe
    top_res_arr = top_df.resSeq.values
    
    # get indices of atoms of residues in each domain
    sd_atom_idx = np.array(top_df.iloc[np.isin(top_res_arr, domain_resi_dict["small_domain"])].index)
    ld_atom_idx = np.array(top_df.iloc[np.isin(top_res_arr, domain_resi_dict["large_domain"])].index)
    loop_atom_idx = np.array(top_df.iloc[np.isin(top_res_arr, domain_resi_dict["loop"])].index)
    
    # create dict with all indices
    domain_atom_idx_dict = {"large_domain":ld_atom_idx, "small_domain":sd_atom_idx, "loop":loop_atom_idx}

    return domain_atom_idx_dict
    
def calc_rmsd_all_domains(variant, conf, run, data_dir):
    """
    Calculate RMSD for each frame in trajectory to (1) the first
    production run structure and (2) the initial structure/crystal
    structure. RMSD is calculated per domain, with domain definitions
    given in domain_resi_dict. 
    """
    # load trajectory 
    traj = md.load(f"{data_dir}/{variant}/{conf}/run/{run}/run_center_prot.xtc",
                   top=f"{data_dir}/{variant}/{conf}/run/{run}/run_center_prot.gro")

    # load matching crystal structure (without Na)
    pdbid_dict = {"closed_conf":"1v4s","superopen_conf":"1v4t"}
    pdbid = pdbid_dict[conf]
    if variant[4:] != "wt":
        initial_struc = md.load(f"../md_intput_files/{pdbid}_{variant[4:]}_setup/processed_original.gro") 
    else:
        initial_struc = md.load(f"../md_intput_files/{pdbid}_setup/processed_original.gro") 

    # get dictionary with domain definitions
    define_domains_dict = define_domains()
    domain_resi_dict = define_domains_dict[conf]

    # get atom indices for each domain
    domain_atom_idx_dict = get_atom_idx_domain(traj, domain_resi_dict, atom_group="backbone-heavy")

    # calculate rmsd to initial structure (crystal structure)
    # and to first structure in production run for all three domains    
    for domain in ['small_domain','large_domain','loop']:

        # get atom indices for given domain
        atom_idx_arr = domain_atom_idx_dict[domain]
        
        # calc rmsd to production run starting structure
        rmsd_to_start = calc_rmsd(traj, traj, atom_idx_arr)
      
        # calc rmsd to input crystal structure
        rmsd_to_crystal = calc_rmsd(traj, initial_struc, atom_idx_arr)
        
        # save
        outdir = "rmsd_all_res"
        np.save(f'../analysis_output/{outdir}/{variant}_{conf}_{run}_{domain}_start.npy', rmsd_to_start)
        np.save(f'../analysis_output/{outdir}/{variant}_{conf}_{run}_{domain}_crystal.npy', rmsd_to_crystal)

def calc_rmsf(variant, conf, run, data_dir):
    """
    Calculate RMSFs for specified GCK variant
    by doing structural alignment and RMSF 
    calculations separately for each domain in GCK.
    """
    # load trajectory
    traj = md.load(f"{data_dir}/{variant}/{conf}/run/{run}/run_center_prot.xtc",
                   top=f"{data_dir}/{variant}/{conf}/run/{run}/run_center_prot.gro")
    
    # define empty rmsd array
    rmsf_arr = np.zeros(465)
    rmsf_arr[:] = np.nan
    
    # get topology information in dataframe format
    top_df, _ = traj.top.to_dataframe()
    
    # keep only ca atoms
    ca_top_df = top_df[top_df.name == 'CA']
    
    # extract residue information from topology dataframe
    ca_top_res_arr = ca_top_df.resSeq.values
        
    # get dictionary with domain definitions
    define_domains_dict = define_domains()    
    domain_resi_dict = define_domains_dict[conf]
    
    # calc rmsf for each domain independently (doing alignment on domain-basis)
    for domain in ["large_domain","small_domain","loop"]:
        
        # get CA atom indices in top
        domain_res_arr = domain_resi_dict[domain] 
        ca_atom_idx = np.array(ca_top_df[np.isin(ca_top_res_arr, domain_res_arr)].index) 
        
        # slice trajectory to contain only CA atoms for target domain
        ca_traj = traj.atom_slice(ca_atom_idx)
        
        # align frames on CA atoms 
        ca_sp_traj = ca_traj.superpose(reference=ca_traj, frame=0)
        
        # calc average structure based on aligned frames
        avg_xyz = np.mean(ca_sp_traj.xyz, axis=0) 
        
        # store average structure coordinates as ca_avg_traj trajectory
        ca_avg_traj = traj.atom_slice(ca_atom_idx)
        ca_avg_traj.xyz = avg_xyz
        
        # calc rmsf to average structure using only CA atoms
        ca_rmsf = md.rmsf(target = ca_sp_traj, reference = ca_avg_traj, frame = 0)
            
        # add result to overall rmsf array
        ca_domain_df, _ = ca_traj.top.to_dataframe()
        rmsf_arr[ca_domain_df.resSeq.values - 1] = ca_rmsf

    outdir = "rmsf_all_res"
    np.save(f'../analysis_output/{outdir}/{variant}_{conf}_{run}.npy', rmsf_arr)

if __name__ == "__main__":
    
    # give name of directory with traj data
    #data_dir =

    variants = ["gck_wt","gck_V455M","gck_G175E","gck_D158A","gck_A173F","gck_G162Q"]
    conformations = ["closed_conf","superopen_conf"]
    runs = ["run1","run2","run3"]
    
    # loop over all variants, conformations and runs to calculate
    # RMSFs for every individual trajectory
    for variant, conf, run in itertools.product(variants,conformations, runs):
        print(variant, conf, run)
        calc_rmsd_all_domains(variant, conf, run, data_dir)
        calc_rmsf(variant, conf, run, data_dir)


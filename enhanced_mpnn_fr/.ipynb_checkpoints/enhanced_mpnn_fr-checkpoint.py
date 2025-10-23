#!/usr/bin/env python

import os
import sys
import glob
import shutil
import argparse
import subprocess
from collections import OrderedDict
import time
import json

# import torch
from pyrosetta import *
from pyrosetta.rosetta import *

# Initialize PyRosetta with required flags
init("-beta_nov16 -in:file:silent_struct_type binary")

# Import silent_tools if needed
sys.path.insert(0, '/home/achen918/software/silent_tools/')
import silent_tools  

def redesign_res_make_cst(pdb_path, cst_dir, omit_AA_json_dir):
    """
    Generates constraints and identifies redesigned residues.
    """
        
    first_peptide_size = 0
    second_peptide_size = 0
    first_flank_size = 0
    second_flank_size = 0
    
    pose = pose_from_pdb(pdb_path)
    pdb_info = pose.pdb_info()
    chain_b_start = pdb_info.number(pose.chain_begin(2))
    chain_a_end = pdb_info.number(pose.chain_end(1))

    # Residues in the diffused region in chain B
    diffused_residue_numbers_2 = [i for i in range(pose.chain_begin(2) + first_peptide_size, pose.chain_end(2) - second_peptide_size + 1)]
    redesigned_residues = ["B" + str(i + chain_b_start - chain_a_end - 1) for i in diffused_residue_numbers_2]

    # Residues in chain A
    chain_a_start = pdb_info.number(pose.chain_begin(1))
    diffused_residue_numbers_1 = [i for i in range(pose.chain_begin(1) + first_flank_size, pose.chain_end(1) - second_flank_size + 1)]
    zinc_residue = pose.chain_begin(3)
    
    his_residues_temp = []
    his_residues = []
    glu_residues_temp = []
    glu_residues = []
    tyr_residue = None
    for i in diffused_residue_numbers_1:
        res_name = pose.residue(i).name3()
        if res_name not in ["HIS", "GLU", "TYR"]:
            redesigned_residues.append("A" + str(i + chain_a_start - 1))
        elif res_name == "HIS":
            his_residues_temp.append(i)
        elif res_name == "GLU":
            glu_residues_temp.append(i)
        elif res_name == "TYR":
            tyr_residue = i
    for glu_residue in glu_residues_temp:
        if (pose.residue(glu_residue-1).name3() == "HIS") and (pose.residue(glu_residue+3).name3() == "HIS"):
            his_residues.append(glu_residue-1)
            his_residues.append(glu_residue+3)
            glu_residues.append(glu_residue)
            break
            
    glu_residues.append([glu for glu in glu_residues_temp if glu != glu_residue][0])
    
    print(f"histidines: {his_residues}")
    print(f"tyrosine: {tyr_residue}")
    print(f"glutamate: {glu_residues}")
        

    # Constraint file creation
    ### Copy the input of invrotzyme over to here!! At least make sure the Zinc-His interactions are all included
    ### (6x3 = 18 interactions at least.)
    ### add constraints for tyrosine: 1. OH coplanar with phenyl ring (dihedral CE1, CZ, OH, target O in CO, absolute value close to 0. ) 2. angle CA,O, O in target CO close to 103 (native) 
    if len(his_residues) == 2 and len(glu_residues) == 2 and tyr_residue:
        chain_B_end = pose.chain_end(2)
        cst_content = f"""\
AtomPair NE2 {his_residues[0]} ZN {chain_B_end+1} HARMONIC 1.92 0.1
Angle O {chain_B_end - 4} ZN {chain_B_end+1} NE2 {his_residues[0]} HARMONIC 2.07 0.1
Angle ZN {chain_B_end+1} NE2 {his_residues[0]} CE1 {his_residues[0]} HARMONIC 2.50 0.1
Dihedral C {chain_B_end - 4} O {chain_B_end - 4} ZN {chain_B_end+1} NE2 {his_residues[0]} HARMONIC 1.68 0.1
Dihedral O {chain_B_end - 4} ZN {chain_B_end+1} NE2 {his_residues[0]} CE1 {his_residues[0]} HARMONIC 1.44 0.1
Dihedral ZN {chain_B_end+1} NE2 {his_residues[0]} CE1 {his_residues[0]} ND1 {his_residues[0]} HARMONIC 2.99 0.1
AtomPair NE2 {his_residues[1]} ZN {chain_B_end+1} HARMONIC 2.02 0.1
Angle O {chain_B_end - 4} ZN {chain_B_end+1} NE2 {his_residues[1]} HARMONIC 2.31 0.1
Angle ZN {chain_B_end+1} NE2 {his_residues[1]} CE1 {his_residues[1]} HARMONIC 2.34 0.1
Dihedral C {chain_B_end - 4} O {chain_B_end - 4} ZN {chain_B_end+1} NE2 {his_residues[1]} HARMONIC -1.56 0.1
Dihedral O {chain_B_end - 4} ZN {chain_B_end+1} NE2 {his_residues[1]} CE1 {his_residues[1]} HARMONIC -1.57 0.1
Dihedral ZN {chain_B_end+1} NE2 {his_residues[1]} CE1 {his_residues[1]} ND1 {his_residues[1]} HARMONIC -2.81 0.1
AtomPair OE1 {glu_residues[1]} ZN {chain_B_end+1} HARMONIC 2.23 0.1
Angle O {chain_B_end - 4} ZN {chain_B_end+1} OE1 {glu_residues[1]} HARMONIC 1.38 0.1
Angle ZN {chain_B_end+1} OE1 {glu_residues[1]} CD {glu_residues[1]} HARMONIC 1.97 0.1
Dihedral C {chain_B_end - 4} O {chain_B_end - 4} ZN {chain_B_end+1} OE1 {glu_residues[1]} HARMONIC -3.01 0.1
Dihedral O {chain_B_end - 4} ZN {chain_B_end+1} OE1 {glu_residues[1]} CD {glu_residues[1]} HARMONIC 1.05 0.1
Dihedral ZN {chain_B_end+1} OE1 {glu_residues[1]} CD {glu_residues[1]} CG {glu_residues[1]} HARMONIC -2.99 0.1
AtomPair O {chain_B_end - 4} OH {tyr_residue} HARMONIC 1.9 0.1
Angle C {chain_B_end - 4} O {chain_B_end - 4} OH {tyr_residue} HARMONIC 2.86 0.1
Angle O {chain_B_end - 4} OH {tyr_residue}  CZ {tyr_residue} HARMONIC 2.29 0.1
Dihedral CA {chain_B_end - 4} C {chain_B_end - 4} O {chain_B_end - 4} OH {tyr_residue}  HARMONIC -1.39 0.1
Dihedral C {chain_B_end - 4} O {chain_B_end - 4} OH {tyr_residue}  CZ {tyr_residue} HARMONIC -2.04 0.1
Dihedral O {chain_B_end - 4} OH {tyr_residue}  CZ {tyr_residue} CE1 {tyr_residue} HARMONIC 2.92 0.1
AtomPair C {chain_B_end - 4} OE1 {glu_residues[0]} HARMONIC 5.1 0.1
AtomPair C {chain_B_end - 4} OE2 {glu_residues[0]} HARMONIC 4.8 0.1
AtomPair ZN {chain_B_end+1} OE1 {glu_residues[0]} HARMONIC 4.3 0.1
AtomPair ZN {chain_B_end+1} OE2 {glu_residues[0]} HARMONIC 4.4 0.1
AtomPair C {chain_B_end - 4} ZN {chain_B_end+1} HARMONIC 3.7 0.1
AtomPair O {chain_B_end - 4} ZN {chain_B_end+1} HARMONIC 3.5 0.1
"""
        cst_file_path = os.path.join(cst_dir, os.path.basename(pdb_path).replace('.pdb', '.cst'))
        with open(cst_file_path, 'w') as cst_file:
            cst_file.write(cst_content)
    else:
        print(f"Skipping {pdb_path} due to missing key residues.")
    
    ### json file creation for banning DHE
    omit_AA_json_path = os.path.join(omit_AA_json_dir, os.path.basename(pdb_path).replace('.pdb', '.json'))
    # Define the dictionary
    data = {
        f"B{pdb_info.number(pose.chain_end(2)) - 1}": "DE",
        f"B{pdb_info.number(pose.chain_end(2)) - 2}": "DE",
        f"B{pdb_info.number(pose.chain_end(2)) - 3}": "DE",
        f"B{pdb_info.number(pose.chain_end(2))}": "DE",
        
    }
    # Write the dictionary to a JSON file
    with open(omit_AA_json_path, "w") as f:
        json.dump(data, f, indent=4)




    return " ".join(redesigned_residues)

def write_mpnn_cmd(pdb_path, 
                   output_folder, 
                   redesigned_residues,
                   omit_AA_json_path):
    """
    Generate the MPNN command string.
    """
    
    command = (
        f"/usr/bin/singularity exec /software/containers/mlfold.sif python /databases/mpnn/fused_mpnn/run.py "
        f"--model_type ligand_mpnn "
        f"--pdb_path {pdb_path} "
        f"--out_folder \"{output_folder}\" "
        f"--redesigned_residues \"{redesigned_residues}\" "
        f"--batch_size 1 "
        f"--number_of_batches 1 "
        f"--temperature 0.1 "
        f"--omit_AA CH "
        f"--omit_AA_per_residue {omit_AA_json_path} "
        f"--chains_to_design A,B "
        f"--pack_side_chains 1 "
        f"--enhance plddt_16_20240910-b65a33eb"
    )
    return command

def relax_pose(pose):
    """
    Apply FastRelax to the given pose.
    """
    FastRelax().apply(pose)
    return pose


def enhanced_mpnnfr(pdb_path, 
                    working_dir, 
                    num_relax_cycles, 
                    save_relaxed_cycles, 
                    cst_dir,
                    rosetta_xml_path,
                    omit_AA_json_dir):
    """
    iterated enhanced mpnn and fr for sequence design
    """
    ## preparation before starting the first cycle
    pdb_file = os.path.basename(pdb_path)
    input_pdb_name = pdb_file.replace('.pdb', '')
    # Create the output folder for this PDB file
    output_folder = os.path.join(working_dir, 'mpnn_output', input_pdb_name)
    os.makedirs(output_folder, exist_ok=True)

    ## prepare writing cmd, make cst files
    redesigned_residues = redesign_res_make_cst(pdb_path, cst_dir, omit_AA_json_dir)   

    ## Set up fast relax
    # Construct the full path to the constraint file using the provided directory and the PDB identifier
    cst_file = os.path.join(cst_dir , f"{input_pdb_name}.cst")
    omit_AA_json_path = os.path.join(omit_AA_json_dir, os.path.basename(pdb_path).replace('.pdb', '.json'))
    
    # Modify the XML file with the correct constraint file

    with open(rosetta_xml_path, 'r') as file:
        xml_content = file.read()
    
    # Replace the placeholder with the actual cst file path
    xml_content = xml_content.replace('{1}', cst_file)
    
    # Write the modified XML to a temporary file
    temp_xml_path = working_dir+ "/TempRosettaFastRelaxUtil.xml"
    with open(temp_xml_path, 'w') as file:
        file.write(xml_content)   
    objs = protocols.rosetta_scripts.XmlObjects.create_from_file(temp_xml_path)
    
    # Load the movers we will need
    FastRelax = objs.get_mover('FastRelax')

    print(f"designing {input_pdb_name}")
    tmp_pdb_file = os.path.join(output_folder, "tmp.pdb") ## this is where the relaxed structures are temperaorily put. 
    shutil.copy(pdb_path, tmp_pdb_file)
    ## enhanced mpnn fr cycles 
    for i,cycle in enumerate(range(num_relax_cycles)):
        cmd = write_mpnn_cmd(tmp_pdb_file, output_folder, redesigned_residues, omit_AA_json_path) ## call redesigned residues and cst file writing inside write mpnn cmd
        print(cmd)
        for file in glob.glob(output_folder + "/*/*.pdb"): # clean up intermediate enhanced mpnn files
            os.remove(file)
        for file in glob.glob(output_folder + "/*/*.fa"): # clean up intermediate enhanced mpnn folder
            os.remove(file)
        result = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True) ## do not use shell = True
        print(result.stdout)  # Print command output
        
        time.sleep(5)  # Wait 5 seconds before checking for files to ensure enhancedmpnn file writing is done.
        packed_pdbs = glob.glob(output_folder + "/packed/*.pdb")
        if not packed_pdbs:
            raise FileNotFoundError(f"No packed PDB files found in {output_folder}/packed/")
        pdb_enhanced_mpnn_packed = packed_pdbs[0] ## only make 1 seq with enhanced mpnn
        pose = pose_from_pdb(pdb_enhanced_mpnn_packed)
        FastRelax.apply(pose)
        if pose:
            pose.dump_pdb(tmp_pdb_file)
        else:
            print("Error: Relaxed pose is None, skipping final PDB save.")

        if save_relaxed_cycles:
            shutil.copy(tmp_pdb_file, output_folder + f"/{input_pdb_name}_mpnnfr_cycle_{i}.pdb")

    pose.dump_pdb(output_folder + f"/{input_pdb_name}_mpnnfr_final.pdb")
    os.remove(tmp_pdb_file)

def main():
    parser = argparse.ArgumentParser(description="Run enhanced MPNN with FastRelax.")
    parser.add_argument("save_relaxed_cycles", type=bool, help="Save intermediate relaxed structures.")
    parser.add_argument("num_relax_cycles", type=int, help="Number of relax cycles.")
    parser.add_argument("rosetta_xml_path", type=str, help="Rosetta temporary xml path")
    parser.add_argument("input_pdb_dir", type=str, help="Input PDB directory.")
    parser.add_argument("cst_dir", type=str, help="Constraint file directory.")
    parser.add_argument("working_dir", type=str, help="Working directory.")
    parser.add_argument("start_index", type=int, help="Start index for PDB files.")
    parser.add_argument("end_index", type=int, help="End index for PDB files.")
    parser.add_argument("omit_AA_json_dir",type=str, help="json path for omitting AA")

    args = parser.parse_args()

    os.makedirs(args.cst_dir, exist_ok=True)
    os.makedirs(args.omit_AA_json_dir, exist_ok=True)
    omit_AA_json_dir = args.omit_AA_json_dir
    cst_dir = args.cst_dir
    input_pdbs = sorted(glob.glob(os.path.join(args.input_pdb_dir, "*.pdb")))[args.start_index:args.end_index + 1]

    for pdb_path in input_pdbs:
        print(f"found {pdb_path}, processing...")
        enhanced_mpnnfr(pdb_path, 
                        args.working_dir, 
                        args.num_relax_cycles, 
                        args.save_relaxed_cycles, 
                        args.cst_dir,
                        args.rosetta_xml_path,
                        args.omit_AA_json_dir)

if __name__ == "__main__":
    main()

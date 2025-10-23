import os
import glob
import json
import re
import pandas as pd
import numpy as np
import argparse
import fnmatch

def process_pdb_files(output_act_dir, start_idx, end_idx, combined_output_file):
    # List to accumulate row dictionaries
    data_list = []

    # Get all descriptions (folders) in the output_act_dir
    descriptions = sorted(os.listdir(output_act_dir))

    # Process only the range specified by start_idx and end_idx
    for description in descriptions[start_idx-1:end_idx]:
        desc_dir = os.path.join(output_act_dir, description)
        input_json_folder = output_act_dir + "/../json_files"
        
        # List all files in the directory
        file_found = False
        for filename in os.listdir(input_json_folder):
            # Check if the filename (lowercased) matches the directory name with ".json"
            if fnmatch.fnmatch(filename.lower(), description.lower() + ".json"):
                input_json = os.path.join(input_json_folder, filename)
                file_found = True
                break  # Stop after finding the first match     
        if file_found is False:
            print(f"No matching JSON file found for {description} at {os.path.join(input_json_folder, filename)}.")
        
        with open(input_json, 'r') as f:
            seq_data = json.load(f)
        chainA_len, chainB_len = [len(entry["protein"]["sequence"]) for entry in seq_data[0]["sequences"][:2]]

        if not os.path.isdir(desc_dir):
            continue  # Skip if not a directory

        plddt_pattern = os.path.join(desc_dir, f"seed-1_sample-*/confidences.json")
        plddt_files = glob.glob(plddt_pattern)

        for plddt_file in plddt_files:

            match = re.search(r"sample-(\d+)/confidences\.json", plddt_file)
            if not match:
                continue
            sample_num = match.group(1)

            try:
                with open(plddt_file, 'r') as f:
                    plddt_data = json.load(f)
                atom_plddts = plddt_data.get("atom_plddts", [])
                if not atom_plddts:
                    print(f"No 'atom_plddts' data in {plddt_file}")
                    continue
                avg_plddt = sum(atom_plddts) / len(atom_plddts)
                
                pdb_path = os.path.join(desc_dir, f"seed-1_sample-{sample_num}/model.pdb")

                pae_matrix = np.array(plddt_data.get("pae", []))

                pae_00 = pae_matrix[:chainA_len, :chainA_len].mean()
                pae_01 = pae_matrix[:chainA_len, chainA_len:chainA_len+chainB_len].mean()
                pae_02 = pae_matrix[:chainA_len, -1].mean()
                pae_10 = pae_matrix[chainA_len:chainA_len+chainB_len, :chainA_len].mean()
                pae_11 = pae_matrix[chainA_len:chainA_len+chainB_len, chainA_len:chainA_len+chainB_len].mean()
                pae_12 = pae_matrix[chainA_len:chainA_len+chainB_len, -1].mean()
                pae_20 = pae_matrix[-1, :chainA_len].mean()
                pae_21 = pae_matrix[-1, chainA_len:chainA_len+chainB_len].mean()
                pae_22 = pae_matrix[-1, -1].mean()
                
                pae_chain_matrix = [[pae_00, pae_01, pae_02],
                               [pae_10, pae_11, pae_12],
                               [pae_20, pae_21, pae_22]]
    
                pae_interface = pae_matrix[chainA_len:chainA_len+chainB_len, :chainA_len].mean()/2+\
                                 pae_matrix[:chainA_len, chainA_len:chainA_len+chainB_len].mean()/2

                data_list.append({
                "description": description,
                "pdb_path": pdb_path,
                "complex_plddt": avg_plddt,
                "per_chain_plddt_A": np.array(atom_plddts)[:chainA_len].mean(),
                "per_chain_plddt_B": np.array(atom_plddts)[chainA_len:chainA_len+chainB_len].mean(),
                "per_chain_plddt_C": np.array(atom_plddts)[-1].mean(),
                "pae_chain_matrix":pae_chain_matrix,
                "pae_interface": pae_interface,
                "chainA_len": chainA_len,
                "chainB_len": chainB_len
            })
            except Exception as e:
                print(f"Error processing PLDDT file {plddt_file}: {e}")
                continue

    # Save results to CSV
    df = pd.DataFrame(data_list)
    df.to_csv(combined_output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process PDB files in parallel.")
    parser.add_argument("output_act_dir", help="Directory containing output folders")
    parser.add_argument("start_idx", type=int, help="Start index for processing")
    parser.add_argument("end_idx", type=int, help="End index for processing")
    parser.add_argument("combined_output_file", help="Path to save the combined CSV file")
    args = parser.parse_args()

    process_pdb_files(args.output_act_dir, args.start_idx, args.end_idx, args.combined_output_file)

    # Print a summary of the job
    num_pdbs_processed = args.end_idx - args.start_idx + 1
    print(f"Processed {num_pdbs_processed} PDBs from index {args.start_idx} to {args.end_idx} in this job.")

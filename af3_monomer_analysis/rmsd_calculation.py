import pandas as pd
from Bio.PDB import PDBParser, Superimposer ### use the PPI_design kernel for this 
import numpy as np
import os
from Bio import PDB  ### use the PPI_design kernel for this 
import sys
from multiprocessing import Pool, cpu_count
import fnmatch
import glob


parser = PDB.PDBParser(QUIET=True) ## do not remove

# Function to extract residue numbers for HEHHY, IR from the template structure
def extract_residue_numbers(template_pdb_path):
    """
    input: path of the template pdb from filtered diffusion outcome.
    output: tuple which contains: list of three histidine residue id and four integers.
    """
    structure = parser.get_structure('template', template_pdb_path)
    his_residues = []
    glu_residues = []
    
    # hardcode fixed flank sizes for kejia's binders, outside of for loop
    first_flank_size = 0
    second_flank_size = 0
    first_peptide_size = 0
    second_peptide_size = 0
    
    chain_a = structure[0]['A']
    chain_b = structure[0]['B']
    diffused_residues_A = [res for i,res in enumerate(chain_a) if (i>=first_flank_size) and i<len(chain_a)-second_flank_size]

    diffused_residues = diffused_residues_A
    
    for residue in diffused_residues:
        if not residue.id[0] == ' ':  # Skip hetero residues
            continue
        res_name = residue.get_resname()
        res_id = residue.get_id()[1]
        if res_name == 'HIS':
            his_residues.append(res_id)
        elif res_name == 'GLU':
            glu_residues.append(res_id)
        elif res_name == 'TYR':
            tyr_residue = res_id


            
    zinc = next(structure[0]['C'].get_residues()).get_id()[1]  
    # zinc = None

    if (len(his_residues)==2) and (len(glu_residues)==2) and zinc:
        # print(f"found all residues for {template_pdb_path}")
        return his_residues, glu_residues, tyr_residue, zinc



# Function to extract atom coordinates for a given residue number and atom name
def get_atom_coordinates(structure, chain_A_length, residue_number, atom_name):
    try:
        if atom_name == 'ZN1': 
            residue = next(structure[0]['C'].get_residues())
        elif residue_number <= chain_A_length:
            chain = structure[0]['A']
            residue = chain[residue_number]
        elif residue_number > chain_A_length:
            chain = structure[0]['B']
            residue = chain[residue_number-chain_A_length]
        # print(f"{atom_name } found in {residue}")
        atom = residue[atom_name]
       
        return atom.coord
    except KeyError:
        print(f"Atom {atom_name} or residue {residue_number} not found.")
        return None


## structure-oracle agreement

def calculate_rmsd(structure1, structure2, atom_selector):
    # Superimpose and calculate RMSD between selected atoms in structures
    super_imposer = Superimposer()
    atoms1 = [atom for atom in structure1.get_atoms() if atom_selector(atom)]
    atoms2 = [atom for atom in structure2.get_atoms() if atom_selector(atom)]

    if len(atoms1) != len(atoms2):
        raise ValueError(f"Atom count mismatch: {len(atoms1)} in structure1 vs {len(atoms2)} in structure2")
    
    super_imposer.set_atoms(atoms1, atoms2)
    super_imposer.apply(structure2.get_atoms())
    return super_imposer.rms

def calculate_side_chain_rmsd(structure_pred, structure_temp, residue_numbers):

    def side_chain_selector(atom):
        residue_id = atom.get_parent().get_id()[1]
        chain_id = atom.get_parent().get_parent().id
        return (atom.get_name() not in ['N', 'CA', 'C', 'O'] and
                residue_id in residue_numbers and
                chain_id == 'A' and
                not atom.get_name().startswith('H'))  # Exclude hydrogens
    
    return calculate_rmsd(structure_pred, structure_temp, side_chain_selector)

def calculate_cata_res_bb_rmsd(structure_pred, structure_temp, residue_numbers):

    def side_chain_selector(atom):
        residue_id = atom.get_parent().get_id()[1]
        chain_id = atom.get_parent().get_parent().id
        return (atom.get_name()in ['N', 'CA', 'C', 'O'] and
                residue_id in residue_numbers and
                chain_id == 'A' and
                not atom.get_name().startswith('H'))  # Exclude hydrogens
    
    return calculate_rmsd(structure_pred, structure_temp, side_chain_selector)



def calculate_ca_rmsd(structure_pred, structure_temp):
    
    def ca_selector(atom):
        return atom.get_name() == 'CA' and atom.get_parent().get_parent().id == 'A'

    return calculate_rmsd(structure_pred, structure_temp, ca_selector)




def main(input_csv, template_pdb_dir, output_csv, start_index, end_index, parser):
    data = pd.read_csv(input_csv)
    
    total_files = len(data)
    print(f'found {total_files} files')
    chunk = data.iloc[start_index:end_index]
    print(f'processing rows from {start_index} to {end_index}')
 
    results = []
    for i, (_, row) in enumerate(chunk.iterrows(), start=start_index):
        result = process_row((row, template_pdb_dir, parser))
        if result:
            results.append(result)
        if (i + 1) % 100 == 0 or (i + 1) == end_index:
            print(f'Processed {i + 1} files out of {total_files} files')

    distance_df = pd.DataFrame(results)
    distance_df.to_csv(output_csv, index=False)
    print(f"Updated file with key residue distances saved to {output_csv}")


    

def process_row(args):
    row, template_pdb_dir, parser = args

    pdb_dir = row['pdb_path']
    description = os.path.basename(os.path.dirname(os.path.dirname(pdb_dir)))
    template_pdb_name_lower = description.split("_mpnnfr")[0]
  
        # List all files in the directory
    file_found = False
    # Expand the wildcard to get all matching directories
    expanded_dirs = glob.glob(template_pdb_dir)
    
    file_found=False
    for dir_path in expanded_dirs:
        # print(dir_path)
        for filename in os.listdir(dir_path):
            # print(filename.lower())
            # print(template_pdb_name_lower + ".pdb")
             # Check if the filename (lowercased) matches the directory name with ".json"
            if fnmatch.fnmatch(filename.lower(), template_pdb_name_lower + ".pdb"):
                template_pdb_path = os.path.join(template_pdb_dir, dir_path, filename)
                file_found=True
                break  # Stop after finding the first match     
        if file_found:
            break
                


    if not os.path.exists(template_pdb_path):
        print(f"Template PDB file missing for {description}")
        return None
    if not os.path.exists(pdb_dir):
        print(f"PDB file missing for {description}")
        return None
    
    pred_structure = parser.get_structure('chai', pdb_dir)
    temp_structure = parser.get_structure('temp', template_pdb_path)
    


    try:
      

        ### calculate structure_prediction_similarity
        rmsd_row={}
        rmsd_row.update(row) ## add row information pdb path, description, etc
        rmsd_row['bb_path'] = template_pdb_path
        all_res_ca_rmsd = calculate_ca_rmsd(temp_structure, pred_structure)
        # Calculate side-chain RMSD for key residues
        
        his_residues, glu_residues, tyr_residue, zinc = extract_residue_numbers(template_pdb_path)
        cata_res_sc_rmsd= calculate_side_chain_rmsd(temp_structure, pred_structure, his_residues+glu_residues+[tyr_residue])
        cata_res_bb_rmsd = calculate_cata_res_bb_rmsd(temp_structure, pred_structure, his_residues+glu_residues+[tyr_residue])
        # Add results to the row
        rmsd_row['all_res_ca_rmsd'] = all_res_ca_rmsd 
        rmsd_row['cata_res_sc_rmsd'] = cata_res_sc_rmsd
        rmsd_row['cata_res_bb_rmsd'] = cata_res_bb_rmsd


        return rmsd_row

    except Exception as e:
        print(f"Error processing {description}: {e}")
        return None




if __name__ == "__main__":
    input_csv = sys.argv[1]
    template_pdb_dir = sys.argv[2]
    output_csv = sys.argv[3]
    start_index = int(sys.argv[4])
    end_index = int(sys.argv[5])
    parser = PDB.PDBParser(QUIET=True) ## do not remove
    main(input_csv, template_pdb_dir, output_csv, start_index, end_index, parser)

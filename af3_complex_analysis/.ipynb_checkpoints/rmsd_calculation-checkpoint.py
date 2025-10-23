##### 250617 all catalytically relevant distances, angles, dihedrals measured based on 4qhp

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
    
    for idx,residue in enumerate(diffused_residues):
        if not residue.id[0] == ' ':  # Skip hetero residues
            continue
        res_name = residue.get_resname()
        res_id = residue.get_id()[1]
        if res_name == 'HIS':
            if diffused_residues[idx+1].get_resname() == "GLU":
                his_residues.append(res_id)
                his_residues.append(res_id + 4)
                glu_residues.append(res_id + 1) ## first glu glu_residues[0] is base, second glu glu_residues[1] is zinc chelator
                
                
        elif res_name == 'GLU':
            if diffused_residues[idx-1].get_resname() != "HIS":
                zinc_chelating_glu_res_id = res_id 
        elif res_name == 'TYR':
            tyr_residue = res_id  
            
    glu_residues.append(zinc_chelating_glu_res_id)


        
            
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

def calculate_ca_rmsd(structure_pred, structure_temp):
    
    def ca_selector(atom):
        chain_id = atom.get_parent().get_parent().id
        return (atom.get_name() == 'CA' and
               chain_id == 'A')

    return calculate_rmsd(structure_pred, structure_temp, ca_selector)



        
def calculate_three_angle_rmsd(angle1, angle2, angle3):
    """
    calculate the three angle root mean square difference from the ideal values of ~ 120 degrees
    input all in degrees
    """

    return np.sqrt((angle1-125.8) **2 + (angle2-134.8)**2 + (angle3-107.7)**2) / np.sqrt(3)

def calculate_three_dihedral_rmsd(dihedral1, dihedral2, dihedral3):
    """
    calculate the three dihedral root mean square difference from the ideal values of ~ 165 degrees, measured from native motif
    input all in degrees
    """

    return np.sqrt((np.abs(dihedral1)-169.8) **2 + (np.abs(dihedral2)-174.4)**2 + (np.abs(dihedral3)-12.2)**2) / np.sqrt(3)
    

def calculate_angle(atom1, atom2, atom3):
    """
    Calculate the angle between three atoms.
    
    Parameters:
    atom1, atom2, atom3: numpy arrays of shape (3,) containing x, y, z coordinates
    
    Returns:
    angle in degrees
    """
    # Calculate vectors
    v1 = atom1 - atom2
    v2 = atom3 - atom2
    
    # Calculate dot product
    dot_product = np.dot(v1, v2)
    
    # Calculate magnitudes
    v1_mag = np.linalg.norm(v1)
    v2_mag = np.linalg.norm(v2)
    
    # Calculate cosine of the angle
    cos_angle = dot_product / (v1_mag * v2_mag)
    
    # Ensure the value is within [-1, 1] to avoid domain errors
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    
    # Calculate angle in radians
    angle_rad = np.arccos(cos_angle)
    
    # Convert to degrees
    angle_deg = np.degrees(angle_rad)
    
    return angle_deg

def calculate_dihedral(atom1, atom2, atom3, atom4):
    """
    Calculate the dihedral angle between four atoms.
    
    Parameters:
    atom1, atom2, atom3, atom4: numpy arrays of shape (3,) containing x, y, z coordinates
    
    Returns:
    dihedral angle in degrees
    """
    # Calculate vectors
    b1 = atom2 - atom1
    b2 = atom3 - atom2
    b3 = atom4 - atom3

    # Calculate normal vectors
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    # Normalize normal vectors
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)

    # Calculate the angle between normal vectors
    x = np.dot(n1, n2)
    y = np.dot(np.cross(n1, b2/np.linalg.norm(b2)), n2)
    
    # Calculate dihedral angle
    angle = np.arctan2(y, x)
    
    # Convert to degrees
    angle_deg = np.degrees(angle)
    
    return angle_deg


def rmsd_important_distances(distances_pred,
                             distances_ideal = np.array([3.2,3.2,3.2,4,4,4,4,4,4])):
    """
    calculate difference from ideal distances
    ========================================
    input:
    output:
    """
    distances = np.array(distances_pred)
    return np.sqrt(((distances - distances_ideal)**2).mean())



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
    description = pdb_dir.split("/output/")[1].split("/seed-1_")[0]
    template_pdb_name_lower = description.split("_mpnnfr_")[0].removeprefix("af3_output_")
    
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



    try:
        his_residues, glu_residues, tyr_residue, zinc = extract_residue_numbers(template_pdb_path)
#### fix this part to make it more general.....        
        temp_structure = parser.get_structure('temp', template_pdb_path)
        chain_B_first = next(temp_structure[0]['B'].get_residues()).get_id()[1]
        chain_B_length = len([res for res in temp_structure[0]['B'].get_residues()])
        chain_A_length = len([res for res in temp_structure[0]['A'] if res.id[0] == ' '])
        chai_structure = parser.get_structure('chai', pdb_dir)

        cleavage_site = chain_B_length - 3

        distances = {}
        his0 = get_atom_coordinates(chai_structure, chain_B_first - 1, his_residues[0], 'NE2')
        his1 = get_atom_coordinates(chai_structure, chain_B_first - 1, his_residues[1], 'NE2')

        co_C = get_atom_coordinates(chai_structure, chain_B_first - 1, cleavage_site-1 +chain_B_first, 'C')
        co_O = get_atom_coordinates(chai_structure, chain_B_first - 1, cleavage_site-1 +chain_B_first, 'O')
        leaving_N = get_atom_coordinates(chai_structure, chain_B_first - 1, cleavage_site +chain_B_first, 'N')


        glu_base_O1 = get_atom_coordinates(chai_structure, chain_B_first - 1, glu_residues[0], 'OE1')# base
        glu_base_O2 = get_atom_coordinates(chai_structure, chain_B_first - 1, glu_residues[0], 'OE2')# base
        
        glu_chel_O1 = get_atom_coordinates(chai_structure, chain_B_first - 1, glu_residues[1], 'OE1')# zinc chelator
        glu_chel_O2 = get_atom_coordinates(chai_structure, chain_B_first - 1, glu_residues[1], 'OE2')# zinc chelator
        
        tyr_OH = get_atom_coordinates(chai_structure, chain_B_first - 1, tyr_residue, 'OH')
        
        zinc_ZN = get_atom_coordinates(chai_structure, chain_B_first - 1, zinc, 'ZN1')

        distances['zinc_his0'] = (zinc_ZN,his0)
        distances['zinc_his1'] = (zinc_ZN,his1)
        distances['zinc_glu_chel_O1'] = (zinc_ZN, glu_chel_O1)
        distances['zinc_glu_chel_O2'] = (zinc_ZN, glu_chel_O2)


    
        distances['co_glu_base_1'] = (glu_base_O1, co_C) #base
        distances['co_glu_base_2'] = (glu_base_O2, co_C) #base
        distances['zinc_sub'] = (co_C, zinc_ZN)
        
        distances['co_tyr'] = (co_O, tyr_OH)
        
        distances['leaving_n_glu_1']= (glu_base_O1, leaving_N) # base
        distances['leaving_n_glu_2']= (glu_base_O2, leaving_N) # base
        
        


        distance_row = {}
        distance_row.update(row)
        for key, (atom1_coord, atom2_coord) in distances.items():
            if atom1_coord is not None and atom2_coord is not None:
                distance = np.linalg.norm(atom1_coord - atom2_coord)
                distance_row[key] = distance
            else:
                distance_row[key] = None
        distance_row['zinc_glu_chel_O'] = min(distance_row['zinc_glu_chel_O1'], distance_row['zinc_glu_chel_O2']) 
        if distance_row['zinc_glu_chel_O'] == distance_row['zinc_glu_chel_O1']:
            glu_chel_Oclose = glu_chel_O1
            glu_chel_Ofar = glu_chel_O2
            
        else:
            glu_chel_Oclose = glu_chel_O2
            glu_chel_Ofar = glu_chel_O1
        
        
        distance_row['rmsd_zinc_chelation'] = rmsd_important_distances(
            [distance_row['zinc_his0'], distance_row['zinc_his1'],distance_row['zinc_glu_chel_O']], 
            distances_ideal = np.array([2.1,2.1,1.9])
        )

        distance_row['leaving_n_glu'] = min(distance_row['leaving_n_glu_1'], distance_row['leaving_n_glu_2'])

        his0_CE1 = get_atom_coordinates(chai_structure, chain_B_first - 1, his_residues[0], 'CE1')
        his0_ND1 = get_atom_coordinates(chai_structure, chain_B_first - 1, his_residues[0], 'ND1')

        his1_CE1 = get_atom_coordinates(chai_structure, chain_B_first - 1, his_residues[1], 'CE1')
        his1_ND1 = get_atom_coordinates(chai_structure, chain_B_first - 1, his_residues[1], 'ND1')

        glu_chel_CD = get_atom_coordinates(chai_structure, chain_B_first - 1, glu_residues[1], 'CD') # zinc chelator

        angle_row={}
        angle_row["angle_zinc_his_0"]=calculate_angle(zinc_ZN, his0, his0_CE1) ##his0 is actually NE2
        angle_row["angle_zinc_his_1"]=calculate_angle(zinc_ZN, his1, his1_CE1) ##his0 is actually NE2
        angle_row["angle_zinc_glu"]=calculate_angle(zinc_ZN, glu_chel_Oclose, glu_chel_CD) ##his0 is actually NE2
        angle_row["three_angle_rmsd"]= calculate_three_angle_rmsd(angle_row["angle_zinc_his_0"],
                                                                  angle_row["angle_zinc_his_1"],
                                                                  angle_row["angle_zinc_glu"]) ### perfect angle is hard coded in the calculate_three_angle_rmsd function
        

        dihedral_row={}
        dihedral_row["dihedral_zinc_his_0"]=calculate_dihedral(zinc_ZN, his0, his0_CE1, his0_ND1) ##his0 is actually NE2
        dihedral_row["dihedral_zinc_his_1"]=calculate_dihedral(zinc_ZN, his1, his1_CE1, his1_ND1) ##his0 is actually NE2
        dihedral_row["dihedral_zinc_glu"]=calculate_dihedral(zinc_ZN, glu_chel_Oclose, glu_chel_CD, glu_chel_Ofar) ##his0 is actually NE2
        dihedral_row["three_dihedral_rmsd"]= calculate_three_dihedral_rmsd(dihedral_row["dihedral_zinc_his_0"],
                                                                  dihedral_row["dihedral_zinc_his_1"],
                                                                  dihedral_row["dihedral_zinc_glu"]) ### perfect angle is hard coded in the 
        
        

        ### calculate structure_prediction_similarity
        rmsd_row={}
        ca_rmsd = calculate_ca_rmsd(temp_structure, chai_structure)
        # Calculate side-chain RMSD for key residues
        sc_rmsd = calculate_side_chain_rmsd(temp_structure, chai_structure, his_residues+glu_residues+[tyr_residue])
        # Add results to the row
        rmsd_row['ca_rmsd'] = ca_rmsd
        rmsd_row['sc_rmsd'] = sc_rmsd

        ### opening angle
        opening_angle_row={}
        opening_angle_row["largest_opening_angle"] = opening_angle(chai_structure)

        ### enzyme size
        enzyme_length_row = {}
        enzyme_length_row["enzyme_length"] = chain_A_length
        
        return distance_row | angle_row | dihedral_row | rmsd_row | opening_angle_row #| enzyme_length_row

    except Exception as e:
        print(f"Error processing {description}: {e}")
        return None


def get_backbone_coords(structure, chain_id, residue_begin_end=(None,None)):
    """Extract Cα, C, and N coordinates from the last 125 residues of a given chain."""
    coords = []
    chain = structure[0][chain_id]
    residues = list(chain)[residue_begin_end[0]:residue_begin_end[1]]  # Get the last 125 residues
    
    for res in residues:
        for atom_name in ['CA', 'C', 'N']:  # Include backbone atoms
            if atom_name in res:
                coords.append(res[atom_name].get_coord())
    
    return np.array(coords)

def best_fit_line(coords):
    """Compute the best-fit 1D vector minimizing distance to all points."""
    centroid = np.mean(coords, axis=0)
    _, _, Vt = np.linalg.svd(coords - centroid)
    direction = Vt[0]  # First principal component
    return centroid, direction

def project_to_plane(coords, plane_point, plane_normal):
    """Project points onto a plane given by a point and a normal vector."""
    plane_normal = plane_normal / np.linalg.norm(plane_normal)
    projections = []
    for point in coords:
        vector_to_plane = point - plane_point
        distance = np.dot(vector_to_plane, plane_normal)
        projection = point - distance * plane_normal
        projections.append(projection)
    return np.array(projections)

def compute_largest_angle(center, points):
    """Find the largest angle between neighboring vectors from center to points."""
    vectors = points - center  # Compute vectors from center to points
    angles = np.arctan2(vectors[:, 1], vectors[:, 0])  # Compute angles w.r.t x-axis
    sorted_indices = np.argsort(angles)  # Sort points by angle
    sorted_angles = angles[sorted_indices]  # Get sorted angles

    # Compute angle differences between consecutive points
    angle_diffs = np.diff(sorted_angles)
    angle_diffs = np.append(angle_diffs, 2 * np.pi + sorted_angles[0] - sorted_angles[-1])  # Wrap around

    largest_angle = np.max(angle_diffs)  # Find the largest angle
    return np.degrees(largest_angle)  # Convert to degrees

def opening_angle(structure):
    """
    calculate the larges open angle after projecting all atoms to a plane perpendicular to the peptide, then connecting the center of peptide to each chain A bb atom. 
    """


    # Step 1: Find best-fit 1D vector for chain B
    ca_chain_B = get_backbone_coords(structure, 'B', residue_begin_end=(0,None))  # Get all backbone atoms
    center_B, best_vector = best_fit_line(ca_chain_B)

    # Step 2: Project all backbone atoms onto the perpendicular plane
    backbone_chain_A = get_backbone_coords(structure, 'A', residue_begin_end=(0,None))  # Get backbone atoms
    projected_B = project_to_plane(ca_chain_B, center_B, best_vector)
    projected_A = project_to_plane(backbone_chain_A, center_B, best_vector)

    # Step 3: Find a 2D coordinate system on the plane
    v1 = np.cross(best_vector, [1, 0, 0])
    if np.linalg.norm(v1) < 1e-6:
        v1 = np.cross(best_vector, [0, 1, 0])
    v1 /= np.linalg.norm(v1)  # Normalize
    v2 = np.cross(best_vector, v1)  # Second basis vector

    # Convert 3D projected points to 2D coordinates
    proj_B_2D = np.array([[np.dot(p - center_B, v1), np.dot(p - center_B, v2)] for p in projected_B])
    proj_A_2D = np.array([[np.dot(p - center_B, v1), np.dot(p - center_B, v2)] for p in projected_A])

    # Step 4: Compute the largest angle
    largest_angle = compute_largest_angle(np.array([0, 0]), proj_A_2D)  # Center is at (0,0)

    # # Step 5: Visualize the projection in 2D
    # plt.figure(figsize=(7, 7))
    # plt.scatter(proj_B_2D[:, 0], proj_B_2D[:, 1], color='blue', label="Projected Chain B (Peptide)")
    # plt.scatter(proj_A_2D[:, 0], proj_A_2D[:, 1], color='red', label="Projected Chain A (Protease)")
    # plt.scatter(0, 0, color='black', marker='x', s=100, label="Center of Chain B")

    # # Draw lines from center to projected A atoms
    # for p in proj_A_2D:
    #     plt.plot([0, p[0]], [0, p[1]], 'gray', alpha=0.5)

    # plt.xlabel("Projected X")
    # plt.ylabel("Projected Y")
    # plt.title(f"2D Projection of Backbone Atoms onto Plane\nLargest Angle = {largest_angle:.2f}°")
    # plt.legend()
    # plt.show()

    # print(f"Largest angle between neighboring vectors: {largest_angle:.2f}°")

    return largest_angle






if __name__ == "__main__":
    input_csv = sys.argv[1]
    template_pdb_dir = sys.argv[2]
    output_csv = sys.argv[3]
    start_index = int(sys.argv[4])
    end_index = int(sys.argv[5])
    parser = PDB.PDBParser(QUIET=True) ## do not remove
    main(input_csv, template_pdb_dir, output_csv, start_index, end_index, parser)

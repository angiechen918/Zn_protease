#!/usr/bin/env python3
import sys
sys.path.append('/home/achen918/scripts_from_others/from_hojae/util') #Copy this folder to your home directory
import os
import argparse
import glob
import csv
from pdb_util import parse_pdb
from calc_dssp import annotate_sse
import torch
from pyrosetta import *
from pyrosetta.rosetta.core.pose import *

from pyrosetta.rosetta import *

import pandas as pd
from Bio.PDB import PDBParser, Superimposer ### use the PPI_design kernel for this 
import numpy as np
from Bio import PDB  ### use the PPI_design kernel for this 

# Initialize PyRosetta with muted output
options = "--mute all"
init(extra_options=options)
# Hardcoded for Kejia's IRP, PLP flanks
first_flank_size = 0
second_flank_size = 0

first_peptide_size = 0
second_peptide_size = 0

def qc_bb(pose, chain_break_thresh = 1.7, floating_res_thresh = 1.7):
    """
    check the quality of generated bb.
    quality indicator 1: cat_residue_in_peptide
    quality indicator 2: chain_break
    quality indicator 3: floating_sidechain
    
    Input: pyrosetta.pose object
    Output: tuple of three booleans: floating_sidechain, cat_residue_in_peptide, chain_break
    """
    
    diffused_residue_numbers_1 = [i for i in range(pose.chain_begin(1) + first_flank_size, pose.chain_end(1) - second_flank_size)]
    diffused_residue_numbers_2 = [i for i in range(pose.chain_begin(2) + first_peptide_size, pose.chain_end(2) - second_peptide_size)]

    ### check cat residue picked up in wrong chain  
    cat_residue_in_peptide = False
    for i in diffused_residue_numbers_2:
        if pose.residue(i).name3() in ["HIS","GLU","TYR"]:
            cat_residue_in_peptide = True
            break

    ### check chain break
    chain_break = False
    for chain_num in [1,2]:
        for i in range(pose.chain_begin(chain_num), pose.chain_end(chain_num)): ## in total num_res - 1 peptide bonds
            carbonyl_C_i_xyz = pose.residue(i).xyz('C')
            nitrogen_iplus1_xyz = pose.residue(i+1).xyz('N')
            # Assuming vector1 and vector2 have been populated with coordinates
            if carbonyl_C_i_xyz.distance(nitrogen_iplus1_xyz)> chain_break_thresh:
                chain_break = True
                break
        if chain_break:
                break

    ### check floating side chains
    floating_sidechain = False
    for i in range(pose.chain_begin(1), pose.chain_end(1)+1): ## check all residues
        if pose.residue(i).name3() in ["HIS","GLU","TYR"]: ### guideposted residues, no check Tyr 
            CA_xyz = pose.residue(i).xyz('CA')
            CB_xyz = pose.residue(i).xyz('CB')
            if CA_xyz.distance(CB_xyz)> floating_res_thresh:
                floating_sidechain = True
                break
    return floating_sidechain, cat_residue_in_peptide, chain_break



## This section checks if the newly grown region itself wraps around peptide and remvoe the ones that do wrap. 
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
    ca_chain_B = get_backbone_coords(structure, 'B', residue_begin_end=(None,None))  # Get all backbone atoms
    center_B, best_vector = best_fit_line(ca_chain_B)

    # Step 2: Project all backbone atoms onto the perpendicular plane
    backbone_chain_A = get_backbone_coords(structure, 'A', residue_begin_end=(None,None))  # Get backbone atoms
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

    
def process_pdb(pdb):
    
    parse = parse_pdb(pdb)
    chain_a_indices = [i for i, chain in enumerate(parse['pdb_idx']) if chain[0] == 'A']
    if not chain_a_indices:
        return None  # Skip if no chain A residues

    ca_xyz_chain_a = parse["xyz"][chain_a_indices, 1, :]  # Grab Ca xyz coordinate of chain A
    dssp = annotate_sse(ca_xyz_chain_a)
    sse_portion = dssp.sum(dim=0) / dssp.shape[0]
    
    dssp_sequence = ''.join(['H' if x == 0 else 'S' if x == 1 else 'L' if x == 2 else 'X' for x in dssp.argmax(dim=1).tolist()])
    
    pose = pose_from_pdb(pdb)
    floating_sidechain, cat_residue_in_peptide, chain_break = qc_bb(pose)
    bioparser = PDB.PDBParser(QUIET=True) ## do not remove
    structure = bioparser.get_structure('bb', pdb)
    largest_angle= opening_angle(structure)

    
    file_name = os.path.basename(pdb)
    return [file_name, 
            sse_portion[0].item(), 
            sse_portion[1].item(), 
            sse_portion[2].item(), 
            sse_portion[3].item(),
            dssp_sequence,
            floating_sidechain, 
            cat_residue_in_peptide, 
            chain_break,
            largest_angle,
            pdb]

def main(csv_file, bb_dir, start=None, end=None):

    all_pdbs = sorted(glob.glob(bb_dir+ "/*.pdb"))[start-1:end]
    
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['file_name', 'H', 'S', 'L', 'X', 'dssp_sequence', 'floating_sidechain', 'cat_residue_in_peptide', 'chain_break', 'largest_opening_angle', 'file_path'])
        
        for pdb in all_pdbs:
            result = process_pdb(pdb)
            if result:
                writer.writerow(result)
        
    print(f'Secondary structure portions and DSSP sequences written to {csv_file}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process PDB files and calculate secondary structure information.")
    parser.add_argument("csv_file", help="Path to the output CSV file")
    parser.add_argument("bb_dir", help="Directory containing PDB files to process")
    parser.add_argument("--start", type=int, help="Start index of PDBs to process")
    parser.add_argument("--end", type=int, help="End index of PDBs to process")
    args = parser.parse_args()

    main(args.csv_file, args.bb_dir, args.start, args.end)

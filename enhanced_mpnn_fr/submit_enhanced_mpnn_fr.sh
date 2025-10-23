#!/bin/bash
#SBATCH -p cpu
#SBATCH --mem=4g
#SBATCH -c 1
#SBATCH -J enhanced_mpnn_fr
#SBATCH -o /net/scratch/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/enhanced_mpnn_fr/log/enhanced_mpnn_fr_%A_%a.out
#SBATCH -t 02:00:00
#SBATCH -a 0-5225

# Set variables
INPUT_PDB_DIR="/net/scratch/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/filter_bb/filtered_bbs"
CST_DIR="/net/scratch/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/enhanced_mpnn_fr/cst_files"
WORKING_DIR="/net/scratch/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/enhanced_mpnn_fr"
ROSETTA_XML="/home/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/enhanced_mpnn_fr/RosettaFastRelaxUtil_cst.xml"
OMIT_AA_JSON_DIR="/home/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/enhanced_mpnn_fr/ban_DEH_from_zinc"

# Calculate indices
START_INDEX=$((SLURM_ARRAY_TASK_ID * 1))
END_INDEX=$((START_INDEX + 1 - 1))

# Ensure END_INDEX doesn't exceed total PDB count
if [ $END_INDEX -ge 5226 ]; then
    END_INDEX=$((5226 - 1))
fi

source activate my-pyrosetta 
# Run the Python script
python /home/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/enhanced_mpnn_fr/enhanced_mpnn_fr.py True 15 "$ROSETTA_XML" "$INPUT_PDB_DIR" "$CST_DIR" "$WORKING_DIR" $START_INDEX $END_INDEX "$OMIT_AA_JSON_DIR"

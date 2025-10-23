#!/bin/bash
#SBATCH -p cpu
#SBATCH --mem=1g
#SBATCH -c 1
#SBATCH -J calculate_geometry
#SBATCH -o /net/scratch/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/af3_complex_analysis//log_catalytic_site_geometry/process_data_%A_%a.out
#SBATCH -t 02:00:00
#SBATCH -a 0-1

# Set variables
INPUT_CSV="/net/scratch/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/af3_complex_analysis//folding_out/af3_data.csv"
TEMPLATE_PDB_DIR="/net/scratch/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/filter_bb/filtered_bbs"
OUTPUT_DIR="/net/scratch/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/af3_complex_analysis//rmsd_calculation_output"

# Calculate chunk size and indices
TOTAL_ROWS=$(wc -l < "$INPUT_CSV")
CHUNK_SIZE=$(( (TOTAL_ROWS + SLURM_ARRAY_TASK_COUNT - 1) / SLURM_ARRAY_TASK_COUNT ))
START_INDEX=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))
END_INDEX=$(( START_INDEX + CHUNK_SIZE ))

# Ensure END_INDEX doesn't exceed TOTAL_ROWS
if [ $END_INDEX -gt $TOTAL_ROWS ]; then
    END_INDEX=$TOTAL_ROWS
fi

# Run the Python script
/software/containers/PPI_design.sif /home/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/af3_complex_analysis/convert_modify_pdb.py "$INPUT_CSV" $START_INDEX $END_INDEX
/software/containers/PPI_design.sif /home/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/af3_complex_analysis/rmsd_calculation.py "$INPUT_CSV" "$TEMPLATE_PDB_DIR" "${OUTPUT_DIR}/output_${SLURM_ARRAY_TASK_ID}.csv" $START_INDEX $END_INDEX

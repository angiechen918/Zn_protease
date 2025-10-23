#!/bin/bash
#SBATCH -p cpu
#SBATCH --mem=4g
#SBATCH -c 1
#SBATCH -J collect_data_job
#SBATCH -o /net/scratch/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/af3_monomer_analysis/folding_out/logs/out.log%a
#SBATCH -t 01:00:00
#SBATCH -a 1-674

# Create a directory for individual CSV outputs
output_dir="/net/scratch/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/af3_monomer_analysis/folding_out/logs/../parallel_csvs"
mkdir -p $output_dir


# Calculate start and end indices for this job
start_idx=$(( ($SLURM_ARRAY_TASK_ID - 1) * 100 + 1 ))
end_idx=$(( $SLURM_ARRAY_TASK_ID * 100 ))

# Run the Python script with the calculated indices
/software/containers/PPI_design.sif /home/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/af3_monomer_analysis/collect_data.py "/net/scratch/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/af3_monomer/output" $start_idx $end_idx "$output_dir/output_$SLURM_ARRAY_TASK_ID.csv"

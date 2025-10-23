#!/bin/bash
#SBATCH -p cpu
#SBATCH --mem=16g
#SBATCH -c 1
#SBATCH -J benchmarl_bb
#SBATCH -o /net/scratch/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion//filter_bb_strand_pair_substrate/log/out.log%a
#SBATCH -t 00:30:00
#SBATCH -a 1-1754

# Calculate start and end indices for this job
start_idx=$(( ($SLURM_ARRAY_TASK_ID - 1) * 50 + 1 ))
end_idx=$(( $SLURM_ARRAY_TASK_ID * 50 ))

# Create a temporary directory for this job
temp_dir="/net/scratch/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion//filter_bb_strand_pair_substrate/parallel_benchmark/tmp/$SLURM_JOB_ID"
mkdir -p $temp_dir

# Process the subset of PDBs for this job
/software/containers/PPI_design.sif /home/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/filter_bb/benchmark.py $temp_dir/output_$SLURM_ARRAY_TASK_ID.csv "/net/scratch/achen918/catalytic_binders/250615_motif_aminopeptidase_N_no_guideposting_free_diffusion/diffusion" --start $start_idx --end $end_idx


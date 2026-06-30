#!/bin/bash
#SBATCH --job-name=parallel_r_jobs    # Job name
# #SBATCH --output=output_%A_%a.txt     # Standard output file (unique for each task)
# #SBATCH --error=error_%A_%a.txt       # Standard error file (unique for each task)
#SBATCH --partition=shared            # Partition name
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --array=1-5                   # Array job with 6 tasks
#SBATCH --mem=20G
#SBATCH --time=48:00:00                # Maximum runtime (D-HH:MM:SS)
module load R 

# Array of R scripts
scripts=("exposure_BMI.R" "exposure_DBP.R" "exposure_HDL.R" "exposure_LDL.R" "exposure_SBP.R" "exposure_logTG.R")

# Run the R script corresponding to the array task ID
Rscript ${scripts[$SLURM_ARRAY_TASK_ID-1]}

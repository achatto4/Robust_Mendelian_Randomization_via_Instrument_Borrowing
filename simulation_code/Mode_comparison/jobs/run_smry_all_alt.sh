#!/bin/bash
#SBATCH --job-name=run_R_scripts
# #SBATCH --output=logs/output_%A_%a.out
# #SBATCH --error=logs/error_%A_%a.err
#SBATCH --array=1-315
#SBATCH --time=48:00:00   # Adjust the time based on expected runtime
#SBATCH --mem=10G          # Adjust memory requirements as needed
#SBATCH --cpus-per-task=4 # Adjust number of CPUs as needed

module load R

# Read the parameter combination for this task
params=$(sed -n "${SLURM_ARRAY_TASK_ID}p" parameter_combinations_alt.txt)
set -- $params
est_theta=$1
thetaU=$2
N=$3
prop_invalid=$4
phi=$5
overlap=$6
theta_alt=$7
NxNy_alt_ratio=$8

# Run the R scripts in parallel with the parameters
Rscript ThetaU_Nalt.R $est_theta $thetaU $N $prop_invalid $phi $overlap $theta_alt $NxNy_alt_ratio

# Wait for all background processes to complete
wait

echo "Task $SLURM_ARRAY_TASK_ID with parameters $params completed."

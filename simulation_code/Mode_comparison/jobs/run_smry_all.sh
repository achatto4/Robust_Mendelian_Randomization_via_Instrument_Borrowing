#!/bin/bash
#SBATCH --job-name=run_R_scripts
# #SBATCH --output=logs/output_%A_%a.out
# #SBATCH --error=logs/error_%A_%a.err
#SBATCH --array=1-8640
#SBATCH --time=100:00:00   # Adjust the time based on expected runtime
#SBATCH --mem=10G          # Adjust memory requirements as needed
#SBATCH --cpus-per-task=4 # Adjust number of CPUs as needed

module load R

# Read the parameter combination for this task
params=$(sed -n "${SLURM_ARRAY_TASK_ID}p" parameter_combinations.txt)
set -- $params
est_theta=$1
thetaU=$2
N=$3
prop_invalid=$4
phi=$5
overlap=$6


# Run the R scripts in parallel with the parameters
Rscript bal_Inside_mode.R $est_theta $thetaU $N $prop_invalid $phi $overlap &
Rscript dir_noInside_mode.R $est_theta $thetaU $N $prop_invalid $phi $overlap &
Rscript bal_noInside_mode.R $est_theta $thetaU $N $prop_invalid $phi $overlap &
Rscript dir_Inside_mode.R $est_theta $thetaU $N $prop_invalid $phi $overlap &

# Wait for all background processes to complete
wait

echo "Task $SLURM_ARRAY_TASK_ID with parameters $params completed."

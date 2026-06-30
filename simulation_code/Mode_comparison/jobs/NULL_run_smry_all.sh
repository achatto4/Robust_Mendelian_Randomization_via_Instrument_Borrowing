#!/bin/bash
#SBATCH --job-name=NULL_run_R_scripts
# #SBATCH --output=logs/null_output_%A_%a.out
# #SBATCH --error=logs/null_error_%A_%a.err
#SBATCH --array=1-400
#SBATCH --time=100:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4

# Avoid module-script failure on clusters where ADDR2LINE is referenced before assignment.
export ADDR2LINE=${ADDR2LINE:-addr2line}

module load R

param_file="parameter_combinations_NULL.txt"

if [[ ! -f "$param_file" ]]; then
  echo "Missing $param_file. Run: bash NULL_param_comb.sh"
  exit 1
fi

params=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$param_file")
if [[ -z "$params" ]]; then
  echo "No parameter row for task ${SLURM_ARRAY_TASK_ID} in $param_file"
  exit 0
fi

set -- $params
theta_alt=$1
thetaU=$2
N=$3
prop_invalid=$4
phi=$5
overlap=$6

Rscript NULL_bal_Inside_mode.R "$theta_alt" "$thetaU" "$N" "$prop_invalid" "$phi" "$overlap" &
Rscript NULL_dir_noInside_mode.R "$theta_alt" "$thetaU" "$N" "$prop_invalid" "$phi" "$overlap" &
Rscript NULL_bal_noInside_mode.R "$theta_alt" "$thetaU" "$N" "$prop_invalid" "$phi" "$overlap" &
Rscript NULL_dir_Inside_mode.R "$theta_alt" "$thetaU" "$N" "$prop_invalid" "$phi" "$overlap" &

wait

echo "Task ${SLURM_ARRAY_TASK_ID} completed with params: $params"

#!/bin/bash
#SBATCH --job-name=POWER_run_R_scripts
# #SBATCH --output=logs/power_output_%A_%a.out
# #SBATCH --error=logs/power_error_%A_%a.err
#SBATCH --array=1-4000
#SBATCH --time=100:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4

# Avoid module-script failure on clusters where ADDR2LINE is referenced before assignment.
export ADDR2LINE=${ADDR2LINE:-addr2line}

module load R

param_file="parameter_combinations_POWER.txt"

if [[ ! -f "$param_file" ]]; then
  echo "Missing $param_file. Run: bash POWER_param_comb.sh"
  exit 1
fi

params=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$param_file")
if [[ -z "$params" ]]; then
  echo "No parameter row for task ${SLURM_ARRAY_TASK_ID} in $param_file"
  exit 0
fi

set -- $params
theta=$1
theta_alt=$2
thetaU=$3
N=$4
prop_invalid=$5
phi=$6
overlap=$7
alt_ratio=$8

Rscript ../scenarios/POWER_dir_noInside_mode.R "$theta" "$theta_alt" "$thetaU" "$N" "$prop_invalid" "$phi" "$overlap" "$alt_ratio"

echo "Task ${SLURM_ARRAY_TASK_ID} completed with params: $params"

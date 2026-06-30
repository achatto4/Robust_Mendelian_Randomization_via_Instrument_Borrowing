#!/bin/bash

# Power simulation parameter grid for POWER_*_mode.R
# Argument order expected by POWER scripts:
# theta_idx theta_alt_idx thetaU_idx N_idx prop_invalid_idx phi_idx overlap_idx alt_ratio_idx

theta_values=(1 2)            # maps to c(-0.1, 0.1)
theta_alt_values=(1 2 3 4 5)  # maps to c(-0.2, -0.1, 0, 0.1, 0.2)
thetaU_values=(1)             # fixed at 0.3 in code
N_values=(3)                  # fixed at 1e5 in code
prop_invalid_values=(3)       # fixed at 0.5 in code
phi_values=(1)                # fixed at phi=1 in code
overlap_values=(4)            # fixed at 0.75 in code
alt_ratio_values=(1 2 3 4 5)  # maps to c(0.2, 0.5, 1, 2, 5)

param_file="parameter_combinations_POWER.txt"
> "$param_file"

for theta in "${theta_values[@]}"; do
  for theta_alt in "${theta_alt_values[@]}"; do
    for thetaU in "${thetaU_values[@]}"; do
      for N in "${N_values[@]}"; do
        for prop_invalid in "${prop_invalid_values[@]}"; do
          for phi in "${phi_values[@]}"; do
            for overlap in "${overlap_values[@]}"; do
              for alt_ratio in "${alt_ratio_values[@]}"; do
                echo "$theta $theta_alt $thetaU $N $prop_invalid $phi $overlap $alt_ratio" >> "$param_file"
              done
            done
          done
        done
      done
    done
  done
done

echo "Parameter combinations file created: $param_file"
echo "Total combinations: $(wc -l < "$param_file")"

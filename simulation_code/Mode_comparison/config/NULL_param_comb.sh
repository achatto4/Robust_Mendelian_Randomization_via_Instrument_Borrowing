#!/bin/bash

# Null simulation parameter grid for NULL_*_mode.R
# Argument order expected by NULL scripts:
# theta_alt_idx thetaU_idx N_idx prop_invalid_idx phi_idx overlap_idx

theta_alt_values=(1 2 3 4 5)      # maps to c(-0.2, -0.1, 0, 0.1, 0.2)
thetaU_values=(1)             # fixed at 0.3 in code
N_values=(1 2 3 4)           # maps to c(5e4, 8e4, 1e5, 2e5)
prop_invalid_values=(1 2 3 4)
phi_values=(1)                # fixed at phi=1 in code
overlap_values=(1 2 3 4 5)

param_file="parameter_combinations_NULL.txt"
> "$param_file"

for theta_alt in "${theta_alt_values[@]}"; do
  for thetaU in "${thetaU_values[@]}"; do
    for N in "${N_values[@]}"; do
      for prop_invalid in "${prop_invalid_values[@]}"; do
        for phi in "${phi_values[@]}"; do
          for overlap in "${overlap_values[@]}"; do
            echo "$theta_alt $thetaU $N $prop_invalid $phi $overlap" >> "$param_file"
          done
        done
      done
    done
  done
done

echo "Parameter combinations file created: $param_file"
echo "Total combinations: $(wc -l < "$param_file")"

#!/bin/bash

# Define parameter values
est_theta_values=(1 2 3)
thetaU_values=(1)
N_values=(2)
prop_invalid_values=(3)
phi_values=(1)
overlap_values=(3)
theta_alt_vec=({1..21})
NxNy_alt_ratio_vec=({1..5})


# Output file for storing parameter combinations
param_file="parameter_combinations_alt.txt"
> $param_file

# Generate all combinations and write them to the file
for est_theta in "${est_theta_values[@]}"; do
    for thetaU in "${thetaU_values[@]}"; do
        for N in "${N_values[@]}"; do
            for prop_invalid in "${prop_invalid_values[@]}"; do
                for phi in "${phi_values[@]}"; do
                    for overlap in "${overlap_values[@]}"; do
                        for theta_alt in "${theta_alt_vec[@]}"; do
                        for NxNy_alt_ratio in "${NxNy_alt_ratio_vec[@]}"; do
                        echo "$est_theta $thetaU $N $prop_invalid $phi $overlap $theta_alt $NxNy_alt_ratio" >> $param_file
                        done
                        done
                    done
                done
            done
        done
    done
done

echo "Parameter combinations file created: $param_file"

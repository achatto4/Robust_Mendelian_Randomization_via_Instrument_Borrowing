#!/bin/bash

# Define parameter values
est_theta_values=(1 2 3)
thetaU_values=(1 2)
N_values=(1 2 3 5)
prop_invalid_values=(1 2 3 4)
phi_values=({1..30})
overlap_values=(1 2 3)

# Output file for storing parameter combinations
param_file="parameter_combinations.txt"
> $param_file

# Generate all combinations and write them to the file
for est_theta in "${est_theta_values[@]}"; do
    for thetaU in "${thetaU_values[@]}"; do
        for N in "${N_values[@]}"; do
            for prop_invalid in "${prop_invalid_values[@]}"; do
                for phi in "${phi_values[@]}"; do
                    for overlap in "${overlap_values[@]}"; do
                        echo "$est_theta $thetaU $N $prop_invalid $phi $overlap" >> $param_file
                    done
                done
            done
        done
    done
done

echo "Parameter combinations file created: $param_file"

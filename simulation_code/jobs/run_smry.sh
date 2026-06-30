#!/bin/bash

module load R

# Define input arguments for all four scenarios
est_theta=2
thetaU=1
N=2
prop_invalid=3 
overlap=2

# Run all four R scripts in parallel, passing the same arguments to each script
Rscript bal_Inside_paper.R $est_theta $thetaU $N $prop_invalid $overlap &
Rscript dir_noInside_paper.R $est_theta $thetaU $N $prop_invalid $overlap &
Rscript bal_noInside_paper.R $est_theta $thetaU $N $prop_invalid $overlap &
Rscript dir_Inside_paper.R $est_theta $thetaU $N $prop_invalid $overlap &

# Wait for all background processes to complete
wait

echo "All scripts have completed."

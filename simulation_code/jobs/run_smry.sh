#!/bin/bash
# Run one simulation scenario (directional pleiotropy, InSIDE violated).
# Arguments are indices into the parameter vectors defined at the top of the
# scenario script: est_theta thetaU N prop_invalid overlap.

module load R

est_theta=2
thetaU=1
N=2
prop_invalid=3
overlap=2

Rscript ../scenarios/dir_noInside_paper.R $est_theta $thetaU $N $prop_invalid $overlap

echo "Done."

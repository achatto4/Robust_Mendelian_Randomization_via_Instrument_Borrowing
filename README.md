# Improving Mendelian Randomization Analysis by Instrument Borrowing from Auxiliary Outcome Traits

Code accompanying the paper **"Improving Mendelian Randomization Analysis by
Instrument Borrowing from Auxiliary Outcome Traits"** by Anagh Chattopadhyay and
Nilanjan Chatterjee.

This repository holds the **simulation** and **real-data analysis** code. The core
instrument-borrowing methods (coheterogeneity screening, IB-Mode, IB-PRESSO) live
in a separate R package, **[IBMR](https://github.com/achatto4/IBMR)**, which this
code depends on.

## Dependencies

```r
# install.packages("devtools")
devtools::install_github("achatto4/IBMR")
```

Other R packages used: `TwoSampleMR`, `MendelianRandomization`, `MRMix`,
`MRPRESSO`, `MR2`, `data.table`, `dplyr`, `tidyr`, `ggplot2`, `MASS`,
`penalized`, `ks`, `ggridges`, `patchwork`, `scales`.

## Contents

Simulation scripts are in [`simulation_code/`](simulation_code); real-data
analysis scripts are in [`RDA_code/`](RDA_code). GWAS summary statistics are not
included; sources are listed in the paper's Data Availability.

The paper reports four simulation scenarios, crossing directional against
balanced pleiotropy with the InSIDE assumption violated or satisfied. This
repository provides the directional, InSIDE-violated scenario
(`dir_noInside_paper.R`), which is the primary setting in the paper, together
with the MR2 comparison (`dir_noInside_MR2.R`). Balanced pleiotropy is obtained
by setting the mean of the direct-effect draws (`alpha`, `alpha_alt`) to zero;
the InSIDE-satisfied variants additionally remove the confounder-mediated
component `phi`, so they are not reachable by a parameter change alone.

Real-data analysis is illustrated with one representative exposure
(`exposure_BMI.R`); the remaining exposures follow the same structure, differing
only in the trait read in and its declared primary and auxiliary outcomes.

## License

MIT — see [LICENSE](LICENSE).

## Citation

Chattopadhyay A, Chatterjee N. *Improving Mendelian Randomization Analysis by Instrument Borrowing from Auxiliary Outcome Traits.*

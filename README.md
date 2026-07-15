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
analysis scripts are in [`RDA_code/`](RDA_code). Scripts are provided for the
directional-pleiotropy, InSIDE-violated scenario reported in the paper (other
scenarios are obtained by changing the scenario parameters), with a
representative real-data exposure analysis (`exposure_BMI.R`). GWAS summary
statistics are not included; sources are listed in the paper's Data Availability.

## Citation

Chattopadhyay A, Chatterjee N. *Improving Mendelian Randomization Analysis by Instrument Borrowing from Auxiliary Outcome Traits.*

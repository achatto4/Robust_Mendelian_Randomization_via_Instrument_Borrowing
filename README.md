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

## Structure

- `RDA_code/` — real-data analysis
  - `exposure_*.R` — one script per exposure (BMI, LDL, SBP, …, vitamin D)
  - `ib_mr_methods_pipeline.R` — shared pipeline sourced by the exposure scripts;
    calls `IBMR::IBMODE`, `IBMR::IBPRESSO`, and the standard MR methods
  - `causalcompMorrison.R` — RDA figures; `table_results.R` — RDA tables
  - `mvmode.sh` — cluster job
- `simulation_code/`
  - `scenarios/` — main paper simulations (balanced/directional × inside/no-inside)
  - `plots/` — figure-generation / publication plotting scripts
  - `jobs/` — batch (SLURM) helpers
  - `Mode_comparison/` — IB-Mode type-I-error calibration and power scenarios
    (`NULL_*` calibration, `POWER_*` power, with the corresponding `*_view.R` plotters)
  - `simul_trait_n.R` — finite-sample validation of the coheterogeneity statistic
    and its standard error (generates the Figure 2 coverage plot); calls
    `IBMR::coheterogeneity_Q`

## Usage

```bash
Rscript simulation_code/scenarios/bal_Inside_paper.R   # a simulation scenario
```

Real-data scripts (`RDA_code/`) require GWAS summary statistics not included here; set the input paths for your environment.

## Citation

Chattopadhyay A, Chatterjee N. *Improving Mendelian Randomization Analysis by Instrument Borrowing from Auxiliary Outcome Traits.*

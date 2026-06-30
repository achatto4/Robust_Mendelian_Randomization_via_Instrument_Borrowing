# Robust Mendelian Randomization via Instrument Borrowing

Code accompanying the paper **"Robust Mendelian Randomization via Instrument
Borrowing"** by Anagh Chattopadhyay and Nilanjan Chatterjee.

This repository holds the **simulation** and **real-data analysis** code. The core
instrument-borrowing methods (coheterogeneity screening, IB-Mode, IB-PRESSO) live
in a separate R package, **[IBMR](https://github.com/achatto4/IBMR)**, which this
code depends on.

## Dependencies

```r
# install.packages("devtools")
devtools::install_github("achatto4/IBMR")
```

Other R packages used: `data.table`, `dplyr`, `tidyr`, `ggplot2`, `MASS`,
`MendelianRandomization`, `MRMix`, `penalized`, `ks`, `MRPRESSO`, `ggridges`,
`patchwork`, `scales`.

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
  - `Mode_comparison/` — IB-Mode vs IB-PRESSO type-I-error & power scenarios
  - `Cohet_theory_validation.R` — finite-sample validation of the coheterogeneity SE

## Usage

```bash
Rscript RDA_code/exposure_BMI.R                       # a real-data exposure analysis
Rscript simulation_code/scenarios/bal_Inside_paper.R # a main simulation scenario
```

> Real-data scripts read GWAS summary statistics from paths set for the original
> compute environment — update these for your setup. For the paper's overlap-aware
> IB-Mode standard errors, pass the cross-trait LDSC intercept for each
> (primary, auxiliary) outcome pair via `ldsc_intercept=` in `run_ibmode()` /
> `harmonize_and_evaluate()`; the default (`NULL`) uses the independent bootstrap.

## Citation

Chattopadhyay A, Chatterjee N. *Robust Mendelian Randomization via Instrument Borrowing.*

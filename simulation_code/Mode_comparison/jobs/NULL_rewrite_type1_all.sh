#!/bin/bash
set -euo pipefail

# Usage:
#   bash NULL_rewrite_type1_all.sh [results_dir] [alpha_csv] [update_rda]
# Example:
#   bash NULL_rewrite_type1_all.sh ../../../results/simulation_results/mode_comp 0.001,0.005,0.01,0.05,0.1 false

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RESULTS_DIR="${1:-../../../results/simulation_results/mode_comp}"
ALPHAS="${2:-0.001,0.005,0.01,0.05,0.1}"
UPDATE_RDA="${3:-false}"

Rscript "$SCRIPT_DIR/NULL_rewrite_type1_from_rda.R" "$RESULTS_DIR" "$ALPHAS" "$UPDATE_RDA"

#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
mkdir -p figures .work
Rscript scripts/01_bodyweight.R
Rscript scripts/02_HE_quantification.R
if [[ -s data/estrous_cycle.csv ]]; then
  Rscript scripts/03_estrous_cycle.R
else
  echo 'SKIP Fig1C/D: no verified estrous observation data supplied.'
fi
if ! bash check_inputs.sh; then
  echo 'STOP: missing omics inputs; Figure 1 scripts have completed, but omics cannot be reproduced.' >&2
  exit 2
fi
Rscript scripts/04_Fig2_single_treatments.R
Rscript scripts/05_Fig3_combination.R
Rscript scripts/06_Fig4_ovarian_RNA.R
Rscript scripts/07_Fig4G_cross_compartment.R
echo 'Finished: output files in figures/; analysis/QC in .work/'

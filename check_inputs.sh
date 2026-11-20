#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
required=(
  'data/bodyweight_long.csv'
  'data/HE_raw_counts.csv'
  'data/raw_snapshot/microbiome_otu_table.tsv'
  'data/raw_snapshot/fecal_all_metab_data.tsv'
  'data/raw_snapshot/serum_all_metab_data.tsv'
  'data/raw_snapshot/rna_gene_count.tsv'
  'data/raw_snapshot/rna_gene_tpm_annotation.tsv'
  'data/taxonomy/seqs_tax_assignments.txt'
)
missing=0
for f in "${required[@]}"; do
  if [[ -s "$ROOT/$f" ]]; then printf '[OK] %s\n' "$f"; else printf '[MISSING] %s\n' "$f"; missing=1; fi
done
if [[ -s "$ROOT/data/estrous_cycle.csv" ]]; then
  echo '[OK] observed estrous-cycle data present; confirm provenance before run'
else
  echo '[OPTIONAL] data/estrous_cycle.csv absent: Fig1C/D not generated'
fi
exit "$missing"

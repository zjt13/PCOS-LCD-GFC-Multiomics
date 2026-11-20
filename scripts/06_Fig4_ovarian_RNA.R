#!/usr/bin/env Rscript
# Reproducible PCOS multi-omics figure pipeline
# Input: data/raw_snapshot/*.tsv and data/taxonomy/seqs_tax_assignments.txt
# Output: figures/*.pdf and *.png; intermediate analysis and QC: .work/
# Run from any directory using Rscript scripts/<this-script>.R
# Individual plot themes are preserved. Figure sizes follow original x0.70 export.
# IMPORTANT: inspect README.md before running, including sample identity and QC notes.

# ---- Portable paths resolved from the location of this script ----
args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_full, value = TRUE)
if (length(script_arg) != 1L) stop("Run using Rscript scripts/<script>.R")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
ROOT <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
REPO_ROOT <- ROOT
DATA_DIR <- file.path(ROOT, "data", "raw_snapshot")
TAXONOMY_INPUT <- file.path(ROOT, "data", "taxonomy", "seqs_tax_assignments.txt")
FIGURE_DIR <- file.path(ROOT, "figures")
WORK_DIR <- file.path(ROOT, ".work")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)

FINAL_FIGURE_MAP <- c("FigS4A_OvaryTranscriptome_PCA"="Fig4A", "Fig4A_DEG_Counts"="Fig4B", "Fig4B_Venn_TreatmentVsPCOS_AllDEG"="Fig4C", "Fig4C_GO_GSEA_TreatmentVsPCOS"="Fig4D", "Fig4D_CombineBeyondBoth_Pathways"="Fig4E", "Fig4F_CircadianInteraction_Genes"="Fig4F", "Fig4H_OriginalFourGenera_CLR"="Fig4H")
# Individual plot theme() settings are preserved during export.
# No global font-size overrides are applied in save_plot().

options(stringsAsFactors = FALSE)
set.seed(20260930)

# ======================================================================
# 0. Paths
# ======================================================================

RAW_DIR <- DATA_DIR

COUNT_FILE <- file.path(
  RAW_DIR,
  "rna_gene_count.tsv"
)

TPM_FILE <- file.path(
  RAW_DIR,
  "rna_gene_tpm_annotation.tsv"
)

MICROBIOME_FILE <- file.path(
  RAW_DIR,
  "microbiome_otu_table.tsv"
)

FECAL_FILE <- file.path(
  RAW_DIR,
  "fecal_all_metab_data.tsv"
)

TAXONOMY_FILE <- TAXONOMY_INPUT

FIG3_TABLE_DIR <- file.path(WORK_DIR,"Fig3","Tables")

FIG3_KEY_GENUS <- file.path(
  FIG3_TABLE_DIR,
  "MAIN_Key_Genera.csv"
)

FIG3_KEY_METAB <- file.path(
  FIG3_TABLE_DIR,
  "MAIN_Key_Annotated_Metabolites.csv"
)

# OPTIONAL reference file for disease-normalization analysis.
#
# Required columns:
#   Gene_name
#   Disease_logFC
#
# Disease_logFC must be PCOS - Healthy/Control:
#   positive = higher in PCOS
#   negative = lower in PCOS
#
# The script NEVER fabricates this direction from the treatment-only data.
PCOS_REFERENCE_SIGNATURE <- file.path(
  RAW_DIR,
  "pcos_reference_signature.tsv"
)

OUT_DIR <- file.path(WORK_DIR,"Fig4")
MAIN_DIR <- file.path(OUT_DIR, "Main")
SUPP_DIR <- file.path(OUT_DIR, "Supplementary")
TABLE_DIR <- file.path(OUT_DIR, "Tables")
GO_DIR <- file.path(TABLE_DIR, "GO")
QC_DIR <- file.path(OUT_DIR, "QC")

dir.create(MAIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GO_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

# ======================================================================
# 1. Packages
# ======================================================================

cran_pkgs <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "ggrepel",
  "pheatmap",
  "tibble",
  "statmod",
  "ggVennDiagram"
)

missing_cran <- cran_pkgs[
  !vapply(
    cran_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_cran) > 0) {
  install.packages(
    missing_cran,
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

bioc_pkgs <- c(
  "DESeq2",
  "clusterProfiler",
  "org.Mm.eg.db"
)

missing_bioc <- bioc_pkgs[
  !vapply(
    bioc_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_bioc) > 0) {
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggrepel)
library(pheatmap)
library(tibble)

# ======================================================================
# 2. Global settings
# ======================================================================

GROUP_LEVELS <- c(
  "PCOS",
  "LCD",
  "GFC",
  "Combine"
)

GROUP_COLS <- c(
  "PCOS" = "#B18459",
  "LCD" = "#8E44AD",
  "GFC" = "#F28E2B",
  "Combine" = "#E76F51"
)

GROUP_SHAPES <- c(
  "PCOS" = 16,
  "LCD" = 15,
  "GFC" = 17,
  "Combine" = 18
)

FDR_CUTOFF <- 0.05
LOGFC_CUTOFF <- 1.00
NOMINAL_P_CUTOFF <- 0.05
GO_MIN_PRIMARY_GENES <- 8

PAIRWISE_CONTRASTS <- c(
  "LCD_vs_PCOS",
  "GFC_vs_PCOS",
  "Combine_vs_PCOS",
  "LCD_vs_GFC",
  "Combine_vs_LCD",
  "Combine_vs_GFC"
)

ALL_ANALYSIS_CONTRASTS <- c(
  PAIRWISE_CONTRASTS,
  "Interaction_LCDxGFC"
)

# ======================================================================
# 3. Helper functions
# ======================================================================

save_plot <- function(p, filename, dir=NULL, width=6, height=4.5, dpi=400, ...) {
  nm <- FINAL_FIGURE_MAP[filename]
  if (length(nm)==0L || is.na(nm) || !nzchar(nm)) return(invisible(NULL))
  w <- as.numeric(width)*0.70; h <- as.numeric(height)*0.70
  ggplot2::ggsave(file.path(FIGURE_DIR,paste0(nm,".pdf")),p,
                  width=w,height=h,units="in",bg="white",limitsize=FALSE)
  ggplot2::ggsave(file.path(FIGURE_DIR,paste0(nm,".png")),p,
                  width=w,height=h,units="in",dpi=400,bg="white",limitsize=FALSE)
  message("Saved manuscript panel: ",nm)
  invisible(NULL)
}


group_from_sample <- function(x) {

  prefix <- sub(
    "_.*$",
    "",
    x
  )

  dplyr::recode(
    prefix,
    "MX" = "PCOS",
    "DT" = "LCD",
    "GZFL" = "GFC",
    "LH" = "Combine",
    "QC" = "QC",
    .default = NA_character_
  )
}

is_rna_sample <- function(x) {
  grepl(
    "^(MX|DT|GZFL|LH)_[0-9]+$",
    x
  )
}

sample_number <- function(x) {
  suppressWarnings(
    as.numeric(
      sub("^.*_", "", x)
    )
  )
}

sample_display <- function(x) {

  g <- group_from_sample(x)
  n <- sample_number(x)

  ifelse(
    !is.na(g) & !is.na(n),
    paste0(g, "_", n),
    x
  )
}

safe_numeric_matrix <- function(
  df,
  sample_cols,
  row_ids
) {

  m <- as.matrix(
    df[
      ,
      sample_cols,
      drop = FALSE
    ]
  )

  suppressWarnings(
    storage.mode(m) <- "numeric"
  )

  rownames(m) <- make.unique(
    as.character(row_ids)
  )

  m
}

format_contrast <- function(x) {

  dplyr::recode(
    x,
    "LCD_vs_PCOS" = "LCD vs PCOS",
    "GFC_vs_PCOS" = "GFC vs PCOS",
    "Combine_vs_PCOS" = "Combine vs PCOS",
    "LCD_vs_GFC" = "LCD vs GFC",
    "Combine_vs_LCD" = "Combine vs LCD",
    "Combine_vs_GFC" = "Combine vs GFC",
    "Interaction_LCDxGFC" = "LCD × GFC interaction",
    .default = x
  )
}

ratio_to_numeric <- function(x) {

  sapply(
    strsplit(
      as.character(x),
      "/",
      fixed = TRUE
    ),
    function(z) {

      if (length(z) != 2) {
        return(NA_real_)
      }

      as.numeric(z[1]) /
        as.numeric(z[2])
    }
  )
}

# ======================================================================
# 4. Read counts + TPM annotation
# ======================================================================

if (!file.exists(COUNT_FILE)) {
  stop(
    "Cannot find count file: ",
    COUNT_FILE
  )
}

if (!file.exists(TPM_FILE)) {
  stop(
    "Cannot find TPM annotation file: ",
    TPM_FILE
  )
}

count_df <- readr::read_tsv(
  COUNT_FILE,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

names(count_df)[1] <- "Gene_id"

count_sample_cols <- names(count_df)[
  vapply(
    names(count_df),
    is_rna_sample,
    logical(1)
  )
]

if (length(count_sample_cols) == 0) {
  stop(
    "No RNA sample columns detected."
  )
}

group_vec <- factor(
  group_from_sample(
    count_sample_cols
  ),
  levels = GROUP_LEVELS
)

if (any(is.na(group_vec))) {
  stop(
    "Some RNA sample names could not be mapped to groups."
  )
}

count_mat <- safe_numeric_matrix(
  count_df,
  count_sample_cols,
  count_df$Gene_id
)

count_mat[
  is.na(count_mat)
] <- 0

count_mat[
  count_mat < 0
] <- 0

# TPM / annotation table
tpm_df <- readr::read_tsv(
  TPM_FILE,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

if (ncol(tpm_df) < 4) {
  stop(
    "rna_gene_tpm_annotation.tsv must contain at least four annotation columns."
  )
}

names(tpm_df)[1:4] <- c(
  "Gene_id",
  "Gene_name",
  "Gene_description",
  "Bio_type"
)

annotation_df <- tpm_df %>%
  dplyr::select(
    Gene_id,
    Gene_name,
    Gene_description,
    Bio_type
  ) %>%
  distinct(
    Gene_id,
    .keep_all = TRUE
  )

annotation_df$Gene_name <- as.character(
  annotation_df$Gene_name
)

annotation_df$Gene_name[
  is.na(annotation_df$Gene_name) |
    annotation_df$Gene_name == ""
] <- annotation_df$Gene_id[
  is.na(annotation_df$Gene_name) |
    annotation_df$Gene_name == ""
]

tpm_sample_cols <- names(tpm_df)[
  vapply(
    names(tpm_df),
    is_rna_sample,
    logical(1)
  )
]

# ======================================================================
# 5. DESeq2 primary differential-expression analysis
#
# The company report states:
#   DESeq2, |log2FC| >= 1, padj <= 0.05
#
# The supplied count matrix contains some fractional expected counts.
# DESeq2 requires integer counts, so this script rounds to the nearest
# non-negative integer. The fraction of non-integer entries is written
# to QC/count_rounding_QC.csv.
# ======================================================================

fractional_fraction <- mean(
  abs(
    count_mat -
      round(count_mat)
  ) > 1e-8
)

count_mat_deseq <- round(
  count_mat
)

count_mat_deseq[
  count_mat_deseq < 0
] <- 0

write_csv(
  data.frame(
    Metric = c(
      "Fraction_non_integer_input_counts",
      "Minimum_input_count",
      "Maximum_input_count",
      "Genes_input",
      "Samples_input"
    ),
    Value = c(
      fractional_fraction,
      min(count_mat, na.rm = TRUE),
      max(count_mat, na.rm = TRUE),
      nrow(count_mat),
      ncol(count_mat)
    )
  ),
  file.path(
    QC_DIR,
    "count_rounding_QC.csv"
  )
)

sample_meta <- data.frame(
  Sample = count_sample_cols,
  Group = factor(
    group_vec,
    levels = GROUP_LEVELS
  ),
  LCD = ifelse(
    group_vec %in% c(
      "LCD",
      "Combine"
    ),
    1,
    0
  ),
  GFC = ifelse(
    group_vec %in% c(
      "GFC",
      "Combine"
    ),
    1,
    0
  ),
  stringsAsFactors = FALSE
)

rownames(sample_meta) <- sample_meta$Sample

# 2x2 factorial design:
# PCOS    LCD=0 GFC=0
# LCD     LCD=1 GFC=0
# GFC     LCD=0 GFC=1
# Combine LCD=1 GFC=1
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_mat_deseq,
  colData = sample_meta,
  design = ~ LCD + GFC + LCD:GFC
)

# Conservative minimal prefilter only; DESeq2 still performs independent
# filtering at the result stage.
dds <- dds[
  rowSums(
    DESeq2::counts(dds) >= 10
  ) >= 2,
]

dds <- DESeq2::DESeq(
  dds,
  quiet = TRUE
)

rn <- DESeq2::resultsNames(dds)

main_lcd <- rn[
  rn == "LCD"
]

main_gfc <- rn[
  rn == "GFC"
]

interaction_name <- rn[
  grepl(
    "LCD.*GFC|GFC.*LCD",
    rn
  )
]

if (
  length(main_lcd) != 1 ||
    length(main_gfc) != 1 ||
    length(interaction_name) != 1
) {

  stop(
    paste0(
      "Could not identify DESeq2 coefficients. resultsNames(dds) = ",
      paste(
        rn,
        collapse = ", "
      )
    )
  )
}

make_numeric_contrast <- function(
  lcd = 0,
  gfc = 0,
  interaction = 0
) {

  z <- rep(
    0,
    length(rn)
  )

  names(z) <- rn

  z[main_lcd] <- lcd
  z[main_gfc] <- gfc
  z[interaction_name] <- interaction

  unname(z)
}

deseq_contrasts <- list(

  LCD_vs_PCOS =
    make_numeric_contrast(
      lcd = 1
    ),

  GFC_vs_PCOS =
    make_numeric_contrast(
      gfc = 1
    ),

  Combine_vs_PCOS =
    make_numeric_contrast(
      lcd = 1,
      gfc = 1,
      interaction = 1
    ),

  LCD_vs_GFC =
    make_numeric_contrast(
      lcd = 1,
      gfc = -1
    ),

  Combine_vs_LCD =
    make_numeric_contrast(
      gfc = 1,
      interaction = 1
    ),

  Combine_vs_GFC =
    make_numeric_contrast(
      lcd = 1,
      interaction = 1
    ),

  Interaction_LCDxGFC =
    make_numeric_contrast(
      interaction = 1
    )
)

get_deseq_table <- function(
  contrast_vector
) {

  res <- DESeq2::results(
    dds,
    contrast = contrast_vector,
    alpha = FDR_CUTOFF,
    independentFiltering = TRUE
  )

  tb <- as.data.frame(
    res
  ) %>%
    tibble::rownames_to_column(
      "Gene_id"
    ) %>%
    dplyr::transmute(
      Gene_id = Gene_id,
      baseMean = baseMean,
      logFC = log2FoldChange,
      lfcSE = lfcSE,
      stat = stat,
      P.Value = pvalue,
      adj.P.Val = padj
    ) %>%
    dplyr::left_join(
      annotation_df,
      by = "Gene_id"
    ) %>%
    dplyr::relocate(
      Gene_id,
      Gene_name,
      Gene_description,
      Bio_type
    )

  tb
}

de_list <- lapply(
  deseq_contrasts,
  get_deseq_table
)

# Convenience object used in later sections.
interaction_tt <- de_list$Interaction_LCDxGFC

for (nm in names(de_list)) {

  write_csv(
    de_list[[nm]],
    file.path(
      TABLE_DIR,
      paste0(
        "DE_",
        nm,
        "_DESeq2.csv"
      )
    )
  )
}

# Variance-stabilized expression for PCA, heatmaps and module scoring.
vsd <- DESeq2::vst(
  dds,
  blind = FALSE
)

norm_expr <- SummarizedExperiment::assay(
  vsd
)

# Save DESeq2 metadata.
writeLines(
  c(
    paste0(
      "DESeq2 resultsNames: ",
      paste(
        rn,
        collapse = "; "
      )
    ),
    paste0(
      "Primary threshold: padj <= ",
      FDR_CUTOFF,
      " and |log2FC| >= ",
      LOGFC_CUTOFF
    ),
    paste0(
      "Fraction of non-integer input counts rounded: ",
      signif(
        fractional_fraction,
        5
      )
    )
  ),
  file.path(
    QC_DIR,
    "DESeq2_analysis_notes.txt"
  )
)

# ======================================================================
# 7. DEG sensitivity analysis
# ======================================================================

threshold_definitions <- data.frame(
  Threshold = c(
    "Company_DESeq2_FDR0.05_FC1",
    "FDR0.10_FC0.58",
    "NominalP0.05_FC1",
    "FDR0.05_NoFC"
  ),
  stringsAsFactors = FALSE
)

count_deg_by_rule <- function(
  tt,
  rule
) {

  if (
    rule ==
      "Company_DESeq2_FDR0.05_FC1"
  ) {

    sig <- tt$adj.P.Val < 0.05 &
      abs(tt$logFC) >= 1

  } else if (
    rule ==
      "FDR0.10_FC0.58"
  ) {

    sig <- tt$adj.P.Val < 0.10 &
      abs(tt$logFC) >= 0.58

  } else if (
    rule ==
      "NominalP0.05_FC1"
  ) {

    sig <- tt$P.Value < 0.05 &
      abs(tt$logFC) >= 1

  } else {

    sig <- tt$adj.P.Val < 0.05
  }

  data.frame(
    Up = sum(
      sig &
        tt$logFC > 0,
      na.rm = TRUE
    ),
    Down = sum(
      sig &
        tt$logFC < 0,
      na.rm = TRUE
    )
  )
}

sensitivity_rows <- list()
ks <- 1

for (nm in ALL_ANALYSIS_CONTRASTS) {

  for (rule in threshold_definitions$Threshold) {

    tmp <- count_deg_by_rule(
      de_list[[nm]],
      rule
    )

    sensitivity_rows[[ks]] <- data.frame(
      Contrast = nm,
      Threshold = rule,
      Up = tmp$Up,
      Down = tmp$Down,
      Total = tmp$Up + tmp$Down,
      stringsAsFactors = FALSE
    )

    ks <- ks + 1
  }
}

deg_sensitivity <- bind_rows(
  sensitivity_rows
)

write_csv(
  deg_sensitivity,
  file.path(
    TABLE_DIR,
    "DEG_threshold_sensitivity.csv"
  )
)

primary_counts <- deg_sensitivity %>%
  filter(
    Threshold ==
      "Company_DESeq2_FDR0.05_FC1"
  )


company_report_counts <- data.frame(
  Contrast = c(
    "LCD_vs_GFC",
    "LCD_vs_PCOS",
    "GFC_vs_PCOS",
    "Combine_vs_LCD",
    "Combine_vs_GFC",
    "Combine_vs_PCOS"
  ),
  Company_Total = c(
    489,
    541,
    239,
    129,
    1443,
    1130
  ),
  Company_Up = c(
    452,
    450,
    171,
    48,
    1265,
    1022
  ),
  Company_Down = c(
    37,
    91,
    68,
    81,
    178,
    108
  ),
  stringsAsFactors = FALSE
)

deg_company_compare <- primary_counts %>%
  filter(
    Contrast %in%
      company_report_counts$Contrast
  ) %>%
  dplyr::select(
    Contrast,
    ThisRun_Up = Up,
    ThisRun_Down = Down,
    ThisRun_Total = Total
  ) %>%
  dplyr::left_join(
    company_report_counts,
    by = "Contrast"
  ) %>%
  mutate(
    Delta_Total =
      ThisRun_Total -
        Company_Total,
    Ratio_Total =
      ifelse(
        Company_Total > 0,
        ThisRun_Total /
          Company_Total,
        NA_real_
      )
  )

write_csv(
  deg_company_compare,
  file.path(
    TABLE_DIR,
    "DEG_counts_vs_company_report.csv"
  )
)

# ======================================================================
# 8. Fig4A | DEG numbers
# ======================================================================

contrast_order <- c(
  "LCD_vs_PCOS",
  "GFC_vs_PCOS",
  "Combine_vs_PCOS",
  "LCD_vs_GFC",
  "Combine_vs_LCD",
  "Combine_vs_GFC",
  "Interaction_LCDxGFC"
)

deg_plot_df <- primary_counts %>%
  dplyr::select(
    Contrast,
    Up,
    Down
  ) %>%
  pivot_longer(
    cols = c(
      Up,
      Down
    ),
    names_to = "Direction",
    values_to = "N"
  ) %>%
  mutate(
    Contrast = factor(
      Contrast,
      levels = contrast_order
    ),
    Nplot = ifelse(
      Direction == "Down",
      -N,
      N
    ),
    Label = ifelse(
      N == 0,
      "",
      as.character(N)
    )
  )

p_A <- ggplot(
  deg_plot_df,
  aes(
    x = Contrast,
    y = Nplot,
    fill = Direction
  )
) +

  geom_col(
    width = 0.68
  ) +

  geom_hline(
    yintercept = 0,
    linewidth = 0.45
  ) +

  geom_text(
    aes(
      label = Label
    ),
    position = position_stack(
      vjust = 0.5
    ),
    size = 3.6
  ) +

  scale_fill_manual(
    values = c(
      "Up" = "#C94C4C",
      "Down" = "#3B73B9"
    ),
    breaks = c(
      "Down",
      "Up"
    )
  ) +

  scale_y_continuous(
    labels = abs,
    expand = expansion(
      mult = c(
        0.08,
        0.10
      )
    )
  ) +

  scale_x_discrete(
    labels = c(
      "LCD_vs_PCOS" = "LCD\nvs PCOS",
      "GFC_vs_PCOS" = "GFC\nvs PCOS",
      "Combine_vs_PCOS" = "Combine\nvs PCOS",
      "LCD_vs_GFC" = "LCD\nvs GFC",
      "Combine_vs_LCD" = "Combine\nvs LCD",
      "Combine_vs_GFC" = "Combine\nvs GFC",
      "Interaction_LCDxGFC" = "LCD×GFC\ninteraction"
    )
  ) +

  labs(
    x = NULL,
    y = "Number of DE genes",
    fill = NULL
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    legend.position = "top",
    axis.text.x = element_text(
      size = 9.5,
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    ),
    axis.title.y = element_text(
      face = "bold"
    )
  )

save_plot(
  p_A,
  "Fig4A_DEG_Counts",
  MAIN_DIR,
  width = 8.2,
  height = 5.1
)

# ======================================================================
# 9. Fig4B | Venn diagram of treatment-vs-PCOS DEGs
#
# Primary Venn uses stable Gene_id (not unique gene symbols) under the same threshold:
#   padj < 0.05 AND |log2FC| >= 1
# Separate Up-only and Down-only Venn diagrams are also written to
# Supplementary.
# ======================================================================

# Use unique stable Gene_id, matching the counting unit of Fig. 4B.
# Using unique(Gene_name) collapsed two distinct LCD gene IDs into
# one symbol and yielded 768 rather than 770 in the original Fig. 4C.
get_primary_deg_ids <- function(
  contrast_name,
  direction = c("All", "Up", "Down")
) {
  direction <- match.arg(direction)
  x <- de_list[[contrast_name]] %>%
    dplyr::filter(
      !is.na(adj.P.Val),
      adj.P.Val < FDR_CUTOFF,
      abs(logFC) >= LOGFC_CUTOFF,
      !is.na(Gene_id),
      Gene_id != ""
    )
  if (direction == "Up") x <- dplyr::filter(x, logFC > 0)
  if (direction == "Down") x <- dplyr::filter(x, logFC < 0)
  unique(as.character(x$Gene_id))
}

venn_all <- list(
  LCD = get_primary_deg_ids(
    "LCD_vs_PCOS",
    "All"
  ),
  GFC = get_primary_deg_ids(
    "GFC_vs_PCOS",
    "All"
  ),
  Combine = get_primary_deg_ids(
    "Combine_vs_PCOS",
    "All"
  )
)

venn_up <- list(
  LCD = get_primary_deg_ids(
    "LCD_vs_PCOS",
    "Up"
  ),
  GFC = get_primary_deg_ids(
    "GFC_vs_PCOS",
    "Up"
  ),
  Combine = get_primary_deg_ids(
    "Combine_vs_PCOS",
    "Up"
  )
)

venn_down <- list(
  LCD = get_primary_deg_ids(
    "LCD_vs_PCOS",
    "Down"
  ),
  GFC = get_primary_deg_ids(
    "GFC_vs_PCOS",
    "Down"
  ),
  Combine = get_primary_deg_ids(
    "Combine_vs_PCOS",
    "Down"
  )
)

# Export exact Venn memberships for transparent interpretation.
venn_membership_table <- function(sets) {
  gene_ids <- sort(unique(unlist(sets, use.names = FALSE)))
  data.frame(
    Gene_id = gene_ids,
    LCD = gene_ids %in% sets$LCD,
    GFC = gene_ids %in% sets$GFC,
    Combine = gene_ids %in% sets$Combine,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(annotation_df %>%
      dplyr::select(Gene_id, Gene_name, Gene_description) %>%
      dplyr::distinct(Gene_id, .keep_all = TRUE),
      by = "Gene_id") %>%
    dplyr::mutate(
      Pattern = paste0(
        ifelse(LCD, "LCD;", ""),
        ifelse(GFC, "GFC;", ""),
        ifelse(Combine, "Combine", "")
      ),
      Pattern = sub(";$", "", Pattern)
    ) %>%
    dplyr::relocate(Gene_id, Gene_name, Gene_description)
}

# Check gene-ID Venn totals against the very same primary DEG-count table
# used for the bar plot. If this fails, do not export an inconsistent figure.
venn_set_mapping <- c(
  LCD = "LCD_vs_PCOS",
  GFC = "GFC_vs_PCOS",
  Combine = "Combine_vs_PCOS"
)
venn_count_audit <- data.frame(
  Treatment = names(venn_set_mapping),
  Contrast = unname(venn_set_mapping),
  VennGeneIDs = vapply(venn_all, length, integer(1)),
  BarDEGs = primary_counts$Total[match(unname(venn_set_mapping),
                                      primary_counts$Contrast)],
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(Delta = VennGeneIDs - BarDEGs)
readr::write_csv(
  venn_count_audit,
  file.path(QC_DIR, "Fig4BC_VennGeneID_vs_DEGCount_Audit.csv")
)
if (anyNA(venn_count_audit$Delta) || any(venn_count_audit$Delta != 0)) {
  stop("Fig4B vs Fig4C mismatch remains; inspect QC/Fig4BC_VennGeneID_vs_DEGCount_Audit.csv")
}

# Exclusive LCD∩Combine excludes GFC; inclusive includes the triple overlap.
venn_tbl_all <- venn_membership_table(venn_all)
venn_overlap_audit <- data.frame(
  Definition = c("LCD_Only", "LCD_and_Combine_Exclusive", 
                 "LCD_and_GFC_Exclusive", "All_Three", "LCD_and_Combine_Inclusive"),
  N_GeneIDs = c(
    sum(venn_tbl_all$LCD & !venn_tbl_all$GFC & !venn_tbl_all$Combine),
    sum(venn_tbl_all$LCD & !venn_tbl_all$GFC & venn_tbl_all$Combine),
    sum(venn_tbl_all$LCD & venn_tbl_all$GFC & !venn_tbl_all$Combine),
    sum(venn_tbl_all$LCD & venn_tbl_all$GFC & venn_tbl_all$Combine),
    sum(venn_tbl_all$LCD & venn_tbl_all$Combine)
  )
)
readr::write_csv(
  venn_overlap_audit,
  file.path(QC_DIR, "Fig4C_ExclusiveInclusive_Overlap_Audit.csv")
)

write_csv(
  venn_membership_table(
    venn_all
  ),
  file.path(
    TABLE_DIR,
    "Venn_TreatmentVsPCOS_AllDEG_membership.csv"
  )
)

write_csv(
  venn_membership_table(
    venn_up
  ),
  file.path(
    TABLE_DIR,
    "Venn_TreatmentVsPCOS_UpDEG_membership.csv"
  )
)

write_csv(
  venn_membership_table(
    venn_down
  ),
  file.path(
    TABLE_DIR,
    "Venn_TreatmentVsPCOS_DownDEG_membership.csv"
  )
)

make_venn_plot <- function(
  sets,
  title_text
) {

  ggVennDiagram::ggVennDiagram(
    sets,
    label_alpha = 0,
    edge_size = 0.8,
    set_size = 4.4,
    label = "count"
  ) +
    scale_fill_gradient(
      low = "white",
      high = "grey65"
    ) +
    labs(
      title = title_text
    ) +
    theme_void(
      base_size = 14
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 13
      ),
      legend.position = "none"
    )
}

p_B <- make_venn_plot(
  venn_all,
  "DEGs vs PCOS"
)

save_plot(
  p_B,
  "Fig4B_Venn_TreatmentVsPCOS_AllDEG",
  MAIN_DIR,
  width = 5.8,
  height = 5.2
)

p_S4D <- make_venn_plot(
  venn_up,
  "Upregulated DEGs vs PCOS"
)

save_plot(
  p_S4D,
  "FigS4D_Venn_TreatmentVsPCOS_UpDEG",
  SUPP_DIR,
  width = 5.8,
  height = 5.2
)

p_S4E <- make_venn_plot(
  venn_down,
  "Downregulated DEGs vs PCOS"
)

save_plot(
  p_S4E,
  "FigS4E_Venn_TreatmentVsPCOS_DownDEG",
  SUPP_DIR,
  width = 5.8,
  height = 5.2
)

# ======================================================================
# 10. GO-BP pathway analysis
#
# MAIN figures use GSEA on the full DESeq2 Wald-statistic ranking.
# This avoids the major weakness of DEG-ORA when one treatment (e.g. GFC)
# has very few threshold-passing DEGs.
#
# ORA is still produced in Tables/GO_ORA as a supplementary sensitivity
# analysis, but it is not used to define the main mechanistic story.
# ======================================================================

GO_ORA_DIR <- file.path(
  GO_DIR,
  "ORA"
)

GO_GSEA_DIR <- file.path(
  GO_DIR,
  "GSEA"
)

dir.create(
  GO_ORA_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  GO_GSEA_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ----------------------------------------------------------------------
# 10.1 GO theme classifier
#
# This does NOT change the enrichment statistics. It only provides a
# transparent way to choose a readable, biologically interpretable subset
# for the main figure. All unfiltered terms are exported.
# ----------------------------------------------------------------------

assign_go_theme <- function(
  description
) {

  x <- tolower(
    description
  )

  dplyr::case_when(

    grepl(
      "circadian|rhythmic|clock",
      x
    ) ~
      "Circadian rhythm",

    grepl(
      "oocyte|ovari|follic|ovulation|reproduct|female gonad|lutein|steroid hormone|gonadotrop",
      x
    ) ~
      "Reproduction / follicle",

    grepl(
      "autoph|apopt|programmed cell death|cell death",
      x
    ) ~
      "Autophagy / apoptosis",

    grepl(
      "oxidative stress|redox|reactive oxygen|inflamm|cytokine|immune|nf.?k|tumor necrosis|interleukin",
      x
    ) ~
      "Redox / inflammation",

    grepl(
      "mitochond|respir|oxidative phosphorylation|atp|electron transport|energy derivation",
      x
    ) ~
      "Mitochondria / energy",

    grepl(
      "lipid|fatty acid|cholesterol|sterol|steroid biosynth|glucose|glycol|insulin|carbohydrate",
      x
    ) ~
      "Metabolism / steroid-lipid",

    grepl(
      "cilium|cilia|axonem|flagell|microtubule|actin|cytoskeleton",
      x
    ) ~
      "Cilium / cytoskeleton",

    grepl(
      "extracellular matrix|collagen|matrix organization|cell adhesion",
      x
    ) ~
      "ECM / adhesion",

    TRUE ~
      "Other"
  )
}

theme_order <- c(
  "Reproduction / follicle",
  "Circadian rhythm",
  "Autophagy / apoptosis",
  "Redox / inflammation",
  "Mitochondria / energy",
  "Metabolism / steroid-lipid",
  "ECM / adhesion",
  "Cilium / cytoskeleton",
  "Other"
)

# ----------------------------------------------------------------------
# 10.2 GO-GSEA on DESeq2 Wald-statistic ranking
# ----------------------------------------------------------------------

make_gene_ranking <- function(
  tt
) {

  r <- tt %>%
    dplyr::filter(
      !is.na(Gene_name),
      Gene_name != "",
      is.finite(stat)
    ) %>%
    dplyr::group_by(
      Gene_name
    ) %>%
    dplyr::slice_max(
      order_by = abs(stat),
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup()

  gene_list <- r$stat
  names(gene_list) <- r$Gene_name

  sort(
    gene_list,
    decreasing = TRUE
  )
}

run_gsego_bp <- function(
  tt,
  contrast_name
) {

  gene_list <- make_gene_ranking(
    tt
  )

  if (length(gene_list) < 100) {
    return(NULL)
  }

  eg <- tryCatch(
    clusterProfiler::gseGO(
      geneList = gene_list,
      OrgDb = org.Mm.eg.db::org.Mm.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      eps = 0,
      verbose = FALSE,
      seed = TRUE
    ),
    error = function(e) {
      message(
        "gseGO failed for ",
        contrast_name,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (
    is.null(eg) ||
      nrow(
        as.data.frame(eg)
      ) == 0
  ) {
    return(NULL)
  }

  tb <- as.data.frame(
    eg
  ) %>%
    dplyr::mutate(
      Contrast = contrast_name,
      Theme = assign_go_theme(
        Description
      ),
      Direction = ifelse(
        NES >= 0,
        "Positive enrichment",
        "Negative enrichment"
      )
    )

  tb
}

gsea_list <- list()

for (nm in ALL_ANALYSIS_CONTRASTS) {

  x <- run_gsego_bp(
    de_list[[nm]],
    nm
  )

  if (!is.null(x)) {

    gsea_list[[nm]] <- x

    write_csv(
      x,
      file.path(
        GO_GSEA_DIR,
        paste0(
          "GSEA_GO_BP_",
          nm,
          ".csv"
        )
      )
    )
  }
}

if (length(gsea_list) > 0) {

  gsea_all <- bind_rows(
    gsea_list
  )

  write_csv(
    gsea_all,
    file.path(
      TABLE_DIR,
      "GO_BP_GSEA_all_contrasts.csv"
    )
  )

} else {

  gsea_all <- data.frame()
}

# ----------------------------------------------------------------------
# 10.3 Representative GO-term selector for MAIN figures
#
# Selection rules:
#   1. FDR < 0.10.
#   2. Prefer biologically interpretable PCOS/ovary themes.
#   3. Keep at most two terms per theme to prevent a redundant cilium or
#      transport family from visually dominating the figure.
#   4. Always add the top two terms by FDR even when they are "Other";
#      therefore the main figure is not silently cherry-picked.
# ----------------------------------------------------------------------

select_representative_gsea_terms <- function(
  df,
  contrasts,
  max_terms = 14,
  n_per_theme = 2
) {

  z <- df %>%
    dplyr::filter(
      Contrast %in%
        contrasts,
      !is.na(p.adjust),
      p.adjust < 0.10
    )

  if (nrow(z) == 0) {
    return(
      data.frame()
    )
  }

  term_summary <- z %>%
    dplyr::group_by(
      ID,
      Description,
      Theme
    ) %>%
    dplyr::summarise(
      MinFDR = min(
        p.adjust,
        na.rm = TRUE
      ),
      MaxAbsNES = max(
        abs(NES),
        na.rm = TRUE
      ),
      Ncontrasts = dplyr::n_distinct(
        Contrast
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Theme = factor(
        Theme,
        levels = theme_order
      )
    )

  preferred <- term_summary %>%
    dplyr::filter(
      Theme != "Other"
    ) %>%
    dplyr::group_by(
      Theme
    ) %>%
    dplyr::arrange(
      MinFDR,
      dplyr::desc(Ncontrasts),
      dplyr::desc(MaxAbsNES),
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(
      n = n_per_theme
    ) %>%
    dplyr::ungroup()

  top_unbiased <- term_summary %>%
    dplyr::arrange(
      MinFDR,
      dplyr::desc(MaxAbsNES)
    ) %>%
    dplyr::slice_head(
      n = 2
    )

  selected <- dplyr::bind_rows(
    preferred,
    top_unbiased
  ) %>%
    dplyr::distinct(
      ID,
      .keep_all = TRUE
    ) %>%
    dplyr::arrange(
      Theme,
      MinFDR,
      dplyr::desc(MaxAbsNES)
    ) %>%
    dplyr::slice_head(
      n = max_terms
    )

  selected
}

make_gsea_heatmap <- function(
  all_df,
  selected_df,
  contrasts,
  contrast_labels,
  filename,
  width,
  dir = MAIN_DIR
) {

  if (
    nrow(selected_df) == 0
  ) {
    return(NULL)
  }

  plot_df <- tidyr::expand_grid(
    ID = selected_df$ID,
    Contrast = contrasts
  ) %>%
    dplyr::left_join(
      selected_df %>%
        dplyr::select(
          ID,
          Description,
          Theme
        ),
      by = "ID"
    ) %>%
    dplyr::left_join(
      all_df %>%
        dplyr::filter(
          Contrast %in%
            contrasts
        ) %>%
        dplyr::select(
          ID,
          Contrast,
          NES,
          p.adjust
        ),
      by = c(
        "ID",
        "Contrast"
      )
    ) %>%
    dplyr::mutate(
      NES = ifelse(
        is.na(NES),
        0,
        NES
      ),
      Sig = dplyr::case_when(
        !is.na(p.adjust) &
          p.adjust < 0.05 ~
          "**",
        !is.na(p.adjust) &
          p.adjust < 0.10 ~
          "*",
        TRUE ~
          ""
      ),
      TermLabel = stringr::str_wrap(
        Description,
        width = 38
      ),
      Theme = factor(
        Theme,
        levels = theme_order
      ),
      Contrast = factor(
        Contrast,
        levels = contrasts
      )
    )

  # Preserve the selected order inside each theme.
  y_order <- selected_df %>%
    dplyr::mutate(
      Theme = factor(
        Theme,
        levels = theme_order
      )
    ) %>%
    dplyr::arrange(
      Theme,
      MinFDR
    ) %>%
    dplyr::pull(
      Description
    )

  y_order <- unique(
    stringr::str_wrap(
      y_order,
      width = 38
    )
  )

  plot_df$TermLabel <- factor(
    plot_df$TermLabel,
    levels = rev(y_order)
  )

  nes_lim <- max(
    abs(
      plot_df$NES
    ),
    na.rm = TRUE
  )

  if (
    !is.finite(nes_lim) ||
      nes_lim < 1
  ) {
    nes_lim <- 1
  }

  p <- ggplot(
    plot_df,
    aes(
      x = Contrast,
      y = TermLabel,
      fill = NES
    )
  ) +

    geom_tile(
      colour = "white",
      linewidth = 0.55
    ) +

    geom_text(
      aes(label = Sig),
      size = 3.8
    ) +

    facet_grid(
      Theme ~ .,
      scales = "free_y",
      space = "free_y",
      switch = "y",
      drop = TRUE
    ) +

    scale_fill_gradient2(
      low = "#3B73B9",
      mid = "white",
      high = "#C94C4C",
      midpoint = 0,
      limits = c(
        -nes_lim,
        nes_lim
      ),
      name = "NES"
    ) +

    scale_x_discrete(
      labels = contrast_labels
    ) +

    labs(
      x = NULL,
      y = NULL,
      caption = "* FDR < 0.10; ** FDR < 0.05"
    ) +

    theme_classic(
      base_size = 14
    ) +

    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(
        angle = 90,
        face = "bold",
        size = 9.3
      ),
      axis.text.x = element_text(
        size = 10,
        face = "bold",
        colour = "black"
      ),
      axis.text.y = element_text(
        size = 9.2,
        colour = "black",
        lineheight = 0.95
      ),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "right",
      plot.caption = element_text(
        size = 9,
        hjust = 0
      ),
      plot.margin = margin(
        8,
        8,
        8,
        14
      )
    )

  n_terms <- length(
    unique(
      plot_df$TermLabel
    )
  )

  save_plot(
    p,
    filename,
    dir,
    width = width,
    height = max(
      5.4,
      2.8 +
        0.34 *
        n_terms
    )
  )

  p
}

# ======================================================================
# 11. Fig4C | GO-GSEA: treatment responses vs PCOS
# ======================================================================

treatment_vs_pcos <- c(
  "LCD_vs_PCOS",
  "GFC_vs_PCOS",
  "Combine_vs_PCOS"
)

if (nrow(gsea_all) > 0) {

  gsea_treat_selected <- select_representative_gsea_terms(
    gsea_all,
    treatment_vs_pcos,
    max_terms = 14,
    n_per_theme = 2
  )

  write_csv(
    gsea_treat_selected,
    file.path(
      TABLE_DIR,
      "Fig4C_GO_GSEA_SelectedTerms.csv"
    )
  )

  p_C <- make_gsea_heatmap(
    all_df = gsea_all,
    selected_df = gsea_treat_selected,
    contrasts = treatment_vs_pcos,
    contrast_labels = c(
      "LCD_vs_PCOS" = "LCD",
      "GFC_vs_PCOS" = "GFC",
      "Combine_vs_PCOS" = "Combine"
    ),
    filename =
      "Fig4C_GO_GSEA_TreatmentVsPCOS",
    width = 8.8
  )
}

# ======================================================================
# 12. Fig4D | Pathways in which Combine extends beyond BOTH monotherapies
#
# This panel does NOT claim that the direction is beneficial.
# It asks a narrower, fully identifiable question from the current data:
#
#   Is the Combine-vs-PCOS pathway effect extended in the SAME direction
#   relative to both LCD and GFC?
#
# Strict evidence:
#   Combine vs PCOS FDR < 0.10
#   Combine vs LCD  FDR < 0.10
#   Combine vs GFC  FDR < 0.10
#   signs of all three NES values are the same
#
# If no pathway meets the strict rule, an exploratory fallback uses:
#   Combine vs PCOS FDR < 0.10
#   Combine vs LCD/GFC nominal P < 0.05
# and is labeled explicitly.
# ======================================================================

combine_pathway_table <- data.frame()

if (nrow(gsea_all) > 0) {

  g_c <- gsea_all %>%
    dplyr::filter(
      Contrast == "Combine_vs_PCOS"
    ) %>%
    dplyr::select(
      ID,
      Description,
      Theme,
      NES_C_vs_PCOS = NES,
      P_C_vs_PCOS = pvalue,
      FDR_C_vs_PCOS = p.adjust,
      Core_C_vs_PCOS = core_enrichment
    )

  g_l <- gsea_all %>%
    dplyr::filter(
      Contrast == "Combine_vs_LCD"
    ) %>%
    dplyr::select(
      ID,
      NES_C_vs_LCD = NES,
      P_C_vs_LCD = pvalue,
      FDR_C_vs_LCD = p.adjust
    )

  g_g <- gsea_all %>%
    dplyr::filter(
      Contrast == "Combine_vs_GFC"
    ) %>%
    dplyr::select(
      ID,
      NES_C_vs_GFC = NES,
      P_C_vs_GFC = pvalue,
      FDR_C_vs_GFC = p.adjust
    )

  combine_pathway_table <- g_c %>%
    dplyr::inner_join(
      g_l,
      by = "ID"
    ) %>%
    dplyr::inner_join(
      g_g,
      by = "ID"
    ) %>%
    dplyr::mutate(

      SameDirection =
        sign(NES_C_vs_PCOS) ==
          sign(NES_C_vs_LCD) &
        sign(NES_C_vs_PCOS) ==
          sign(NES_C_vs_GFC),

      StrictBeyondBoth =
        !is.na(FDR_C_vs_PCOS) &
        !is.na(FDR_C_vs_LCD) &
        !is.na(FDR_C_vs_GFC) &
        FDR_C_vs_PCOS < 0.10 &
        FDR_C_vs_LCD < 0.10 &
        FDR_C_vs_GFC < 0.10 &
        SameDirection,

      ExploratoryBeyondBoth =
        !is.na(FDR_C_vs_PCOS) &
        FDR_C_vs_PCOS < 0.10 &
        !is.na(P_C_vs_LCD) &
        P_C_vs_LCD < 0.05 &
        !is.na(P_C_vs_GFC) &
        P_C_vs_GFC < 0.05 &
        SameDirection,

      EvidenceMode = dplyr::case_when(
        StrictBeyondBoth ~
          "Strict FDR in all 3 contrasts",
        ExploratoryBeyondBoth ~
          "Exploratory: Combine FDR + direct nominal P",
        TRUE ~
          "Not beyond both"
      ),

      Direction = ifelse(
        NES_C_vs_PCOS > 0,
        "Higher in Combine",
        "Lower in Combine"
      ),

      PriorityScore =
        -log10(
          pmax(
            FDR_C_vs_PCOS,
            1e-300
          )
        ) +
        -log10(
          pmax(
            ifelse(
              StrictBeyondBoth,
              pmax(
                FDR_C_vs_LCD,
                FDR_C_vs_GFC
              ),
              pmax(
                P_C_vs_LCD,
                P_C_vs_GFC
              )
            ),
            1e-300
          )
        )
    )

  write_csv(
    combine_pathway_table,
    file.path(
      TABLE_DIR,
      "Combine_BeyondBoth_Pathways_All.csv"
    )
  )

  strict_path <- combine_pathway_table %>%
    dplyr::filter(
      StrictBeyondBoth
    )

  if (nrow(strict_path) > 0) {

    combine_path_plot <- strict_path

  } else {

    combine_path_plot <- combine_pathway_table %>%
      dplyr::filter(
        ExploratoryBeyondBoth
      )
  }

  if (nrow(combine_path_plot) > 0) {

    # Main display: biologically interpretable terms first, but retain the
    # top two unbiased terms to avoid hidden cherry-picking.
    interpretable_path <- combine_path_plot %>%
      dplyr::filter(
        Theme != "Other"
      ) %>%
      dplyr::group_by(
        Theme
      ) %>%
      dplyr::arrange(
        dplyr::desc(PriorityScore),
        .by_group = TRUE
      ) %>%
      dplyr::slice_head(
        n = 2
      ) %>%
      dplyr::ungroup()

    unbiased_path <- combine_path_plot %>%
      dplyr::arrange(
        dplyr::desc(PriorityScore)
      ) %>%
      dplyr::slice_head(
        n = 2
      )

    combine_path_plot <- dplyr::bind_rows(
      interpretable_path,
      unbiased_path
    ) %>%
      dplyr::distinct(
        ID,
        .keep_all = TRUE
      ) %>%
      dplyr::arrange(
        Direction,
        dplyr::desc(PriorityScore)
      ) %>%
      dplyr::slice_head(
        n = 12
      ) %>%
      dplyr::mutate(
        TermLabel = stringr::str_wrap(
          Description,
          width = 38
        ),
        Direction = factor(
          Direction,
          levels = c(
            "Higher in Combine",
            "Lower in Combine"
          )
        )
      )

    write_csv(
      combine_path_plot,
      file.path(
        TABLE_DIR,
        "Fig4D_CombineBeyondBoth_SelectedPathways.csv"
      )
    )

    # FINAL Fig4E: one size legend only. Red/blue fills show the NES
    # direction, while EvidenceMode remains available in the exported table.
    # The complete screening logic and prioritization scores are unchanged.
    p_D_main <- ggplot(
      combine_path_plot,
      aes(
        x = NES_C_vs_PCOS,
        y = reorder(TermLabel, NES_C_vs_PCOS),
        size = PriorityScore,
        fill = Direction
      )
    ) +
      geom_vline(
        xintercept = 0,
        colour = "grey75",
        linewidth = 0.45
      ) +
      geom_point(
        shape = 21,              # same filled circle for all terms
        colour = "grey35",       # outline stays constant
        stroke = 0.35,
        alpha = 0.95
      ) +
      scale_fill_manual(
        values = c(
          "Higher in Combine" = "#C94C4C",
          "Lower in Combine"  = "#3B73B9"
        ),
        guide = "none"           # retain colors; remove direction legend
      ) +
      scale_size_continuous(
        range = c(4.5, 9.0),
        name = "Prioritization score",
        breaks = scales::pretty_breaks(n = 3)
      ) +
      guides(
        size = guide_legend(
          direction = "horizontal",
          title.position = "left",
          override.aes = list(
            fill = "grey65",
            colour = "grey35",
            shape = 21,
            alpha = 1
          )
        ),
        fill = "none",
        shape = "none"
      ) +
      labs(
        x = "NES: Combine vs PCOS",
        y = NULL,
        caption = NULL
      ) +
      theme_classic(base_size = 14) +
      theme(
        axis.text.x = element_text(
          size = 11,
          colour = "black"
        ),
        axis.text.y = element_text(
          size = 9.5,
          colour = "black",
          lineheight = 0.94
        ),
        axis.title.x = element_text(
          size = 12,
          face = "bold"
        ),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 10.5),
        legend.text = element_text(size = 10),
        legend.key.width = grid::unit(0.55, "cm"),
        legend.key.height = grid::unit(0.55, "cm"),
        legend.spacing.x = grid::unit(0.12, "cm"),
        legend.margin = margin(t = 1, b = 0),
        plot.caption = element_blank(),
        plot.margin = margin(6, 8, 5, 6)
      )

    save_plot(
      p_D_main,
      "Fig4D_CombineBeyondBoth_Pathways",
      MAIN_DIR,
      width = 9.0,
      height = max(
        5.0,
        2.8 +
          0.38 *
          nrow(
            combine_path_plot
          )
      )
    )

  } else {

    writeLines(
      "No pathway met either the strict or exploratory Combine-beyond-both definition.",
      file.path(
        QC_DIR,
        "Fig4D_No_CombineBeyondBoth_Pathway.txt"
      )
    )
  }
}

# ======================================================================
# 13. Fig4E | Genes in which Combine extends beyond BOTH monotherapies
#
# Strict gene-level definition:
#   Combine vs PCOS: FDR < 0.05 and |log2FC| >= 1
#   Combine vs LCD : FDR < 0.10 and |log2FC| >= 0.58
#   Combine vs GFC : FDR < 0.10 and |log2FC| >= 0.58
#   all three effects have the same sign
#
# Again: "beyond both" means a stronger treatment response away from PCOS,
# NOT proven normalization toward a healthy state.
# ======================================================================

cpc <- de_list$Combine_vs_PCOS %>%
  dplyr::select(
    Gene_id,
    Gene_name,
    Gene_description,
    FC_C_vs_PCOS = logFC,
    P_C_vs_PCOS = P.Value,
    FDR_C_vs_PCOS = adj.P.Val
  )

cld <- de_list$Combine_vs_LCD %>%
  dplyr::select(
    Gene_id,
    FC_C_vs_LCD = logFC,
    P_C_vs_LCD = P.Value,
    FDR_C_vs_LCD = adj.P.Val
  )

cgf <- de_list$Combine_vs_GFC %>%
  dplyr::select(
    Gene_id,
    FC_C_vs_GFC = logFC,
    P_C_vs_GFC = P.Value,
    FDR_C_vs_GFC = adj.P.Val
  )

combine_beyond_genes <- cpc %>%
  dplyr::inner_join(
    cld,
    by = "Gene_id"
  ) %>%
  dplyr::inner_join(
    cgf,
    by = "Gene_id"
  ) %>%
  dplyr::mutate(

    SameDirection =
      sign(FC_C_vs_PCOS) ==
        sign(FC_C_vs_LCD) &
      sign(FC_C_vs_PCOS) ==
        sign(FC_C_vs_GFC),

    StrictBeyondBoth =
      !is.na(FDR_C_vs_PCOS) &
      !is.na(FDR_C_vs_LCD) &
      !is.na(FDR_C_vs_GFC) &
      FDR_C_vs_PCOS < 0.05 &
      abs(FC_C_vs_PCOS) >= 1 &
      FDR_C_vs_LCD < 0.10 &
      abs(FC_C_vs_LCD) >= 0.58 &
      FDR_C_vs_GFC < 0.10 &
      abs(FC_C_vs_GFC) >= 0.58 &
      SameDirection,

    Direction = ifelse(
      FC_C_vs_PCOS > 0,
      "Higher in Combine",
      "Lower in Combine"
    ),

    WorstDirectFDR = pmax(
      FDR_C_vs_LCD,
      FDR_C_vs_GFC,
      na.rm = TRUE
    ),

    MinDirectMagnitude = pmin(
      abs(FC_C_vs_LCD),
      abs(FC_C_vs_GFC)
    )
  ) %>%
  dplyr::filter(
    StrictBeyondBoth
  ) %>%
  dplyr::arrange(
    WorstDirectFDR,
    dplyr::desc(
      MinDirectMagnitude
    ),
    FDR_C_vs_PCOS
  )

write_csv(
  combine_beyond_genes,
  file.path(
    TABLE_DIR,
    "Combine_BeyondBoth_Genes_Strict.csv"
  )
)

# Map genes to selected Combine-beyond-both pathways when possible.
combine_core_genes <- character(0)

if (
  exists("combine_path_plot") &&
    nrow(combine_path_plot) > 0
) {

  combine_core_genes <- unique(
    unlist(
      strsplit(
        paste(
          combine_path_plot$Core_C_vs_PCOS,
          collapse = "/"
        ),
        "/",
        fixed = TRUE
      )
    )
  )
}

known_ovary_story_genes <- unique(
  c(
    "Star",
    "Cyp11a1",
    "Cyp17a1",
    "Cyp19a1",
    "Hsd3b1",
    "Fshr",
    "Lhcgr",
    "Pgr",
    "Ptgs2",
    "Areg",
    "Ereg",
    "Adamts1",
    "Mtor",
    "Rptor",
    "Rictor",
    "Becn1",
    "Atg5",
    "Atg7",
    "Sqstm1",
    "Nfkb1",
    "Rela",
    "Nfkbia",
    "Nfe2l2",
    "Hmox1",
    "Nqo1",
    "Sod2",
    "Pparg",
    "Ppara",
    "Ppargc1a",
    "Cpt1a",
    "Tfam"
  )
)

combine_beyond_genes <- combine_beyond_genes %>%
  dplyr::mutate(
    InSelectedPathwayCore =
      Gene_name %in%
        combine_core_genes,

    KnownOvaryMechanism =
      Gene_name %in%
        known_ovary_story_genes,

    SelectionTier = dplyr::case_when(
      InSelectedPathwayCore &
        KnownOvaryMechanism ~ 1,
      InSelectedPathwayCore ~ 2,
      KnownOvaryMechanism ~ 3,
      TRUE ~ 4
    )
  ) %>%
  dplyr::arrange(
    SelectionTier,
    WorstDirectFDR,
    dplyr::desc(
      MinDirectMagnitude
    ),
    FDR_C_vs_PCOS
  )

write_csv(
  combine_beyond_genes,
  file.path(
    TABLE_DIR,
    "Combine_BeyondBoth_Genes_Ranked.csv"
  )
)

if (nrow(combine_beyond_genes) > 0) {

  heat_gene_tbl <- combine_beyond_genes %>%
    dplyr::filter(
      !is.na(Gene_name),
      Gene_name != ""
    ) %>%
    dplyr::distinct(
      Gene_name,
      .keep_all = TRUE
    ) %>%
    dplyr::slice_head(
      n = 24
    )

  heat_ids <- intersect(
    heat_gene_tbl$Gene_id,
    rownames(
      norm_expr
    )
  )

  if (length(heat_ids) >= 2) {

    heat_beyond <- norm_expr[
      heat_ids,
      ,
      drop = FALSE
    ]

    heat_beyond <- t(
      scale(
        t(
          heat_beyond
        )
      )
    )

    heat_beyond[
      !is.finite(
        heat_beyond
      )
    ] <- 0

    labels_beyond <- heat_gene_tbl$Gene_name[
      match(
        rownames(
          heat_beyond
        ),
        heat_gene_tbl$Gene_id
      )
    ]

    labels_beyond[
      is.na(labels_beyond) |
        labels_beyond == ""
    ] <- rownames(
      heat_beyond
    )[
      is.na(labels_beyond) |
        labels_beyond == ""
    ]

    rownames(
      heat_beyond
    ) <- make.unique(
      labels_beyond
    )

    sample_order_beyond <- sample_meta %>%
      dplyr::mutate(
        SampleNumber =
          sample_number(
            Sample
          ),
        DisplaySample =
          sample_display(
            Sample
          )
      ) %>%
      dplyr::arrange(
        Group,
        SampleNumber
      )

    heat_beyond <- heat_beyond[
      ,
      sample_order_beyond$Sample,
      drop = FALSE
    ]

    colnames(
      heat_beyond
    ) <- sample_order_beyond$DisplaySample

    anno_beyond <- data.frame(
      Group = factor(
        as.character(
          sample_order_beyond$Group
        ),
        levels = GROUP_LEVELS
      ),
      row.names =
        sample_order_beyond$DisplaySample
    )

    anno_colors_beyond <- list(
      Group = GROUP_COLS
    )

    gaps_beyond <- cumsum(
      as.integer(
        table(
          factor(
            sample_order_beyond$Group,
            levels = GROUP_LEVELS
          )
        )
      )
    )

    gaps_beyond <- gaps_beyond[
      -length(
        gaps_beyond
      )
    ]

    pdf(
      file.path(tempdir(),"Fig4E_CombineBeyondBoth_GeneHeatmap.pdf"),
      width = 9.6,
      height = max(
        5.8,
        2.8 +
          0.22 *
          nrow(
            heat_beyond
          )
      )
    )

    pheatmap::pheatmap(
      heat_beyond,
      annotation_col =
        anno_beyond,
      annotation_colors =
        anno_colors_beyond,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      gaps_col = gaps_beyond,
      show_colnames = TRUE,
      show_rownames = TRUE,
      fontsize = 11,
      fontsize_row = 9.3,
      fontsize_col = 9.5,
      angle_col = 45,
      border_color = NA,
      main = ""
    )

    dev.off()

    png(
      file.path(tempdir(),"Fig4E_CombineBeyondBoth_GeneHeatmap.png"),
      width = 9.6,
      height = max(
        5.8,
        2.8 +
          0.22 *
          nrow(
            heat_beyond
          )
      ),
      units = "in",
      res = 400
    )

    pheatmap::pheatmap(
      heat_beyond,
      annotation_col =
        anno_beyond,
      annotation_colors =
        anno_colors_beyond,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      gaps_col = gaps_beyond,
      show_colnames = TRUE,
      show_rownames = TRUE,
      fontsize = 11,
      fontsize_row = 9.3,
      fontsize_col = 9.5,
      angle_col = 45,
      border_color = NA,
      main = ""
    )

    dev.off()
  }

} else {

  writeLines(
    "No gene met the strict Combine-beyond-both definition.",
    file.path(
      QC_DIR,
      "Fig4E_No_Strict_CombineBeyondBoth_Gene.txt"
    )
  )
}

# ======================================================================
# 14. Fig4F | OPTIONAL PCOS disease-signature reversal
#
# This is the ONLY panel designed to answer whether a treatment response
# moves opposite to a PCOS-vs-healthy disease signature.
#
# It is generated only when pcos_reference_signature.tsv exists.
#
# Required:
#   Gene_name
#   Disease_logFC  (PCOS - Healthy)
#
# ReversalScore = -Spearman(Disease_logFC, Treatment_logFC)
#   positive -> treatment tends to oppose the disease signature
#   zero     -> no systematic reversal
#   negative -> treatment tends to reinforce the disease signature
# ======================================================================

if (file.exists(
  PCOS_REFERENCE_SIGNATURE
)) {

  ref_sig <- readr::read_tsv(
    PCOS_REFERENCE_SIGNATURE,
    show_col_types = FALSE,
    progress = FALSE
  )

  required_ref_cols <- c(
    "Gene_name",
    "Disease_logFC"
  )

  if (!all(
    required_ref_cols %in%
      names(ref_sig)
  )) {

    stop(
      paste0(
        "pcos_reference_signature.tsv must contain: ",
        paste(
          required_ref_cols,
          collapse = ", "
        )
      )
    )
  }

  ref_sig <- ref_sig %>%
    dplyr::filter(
      !is.na(Gene_name),
      Gene_name != "",
      is.finite(
        Disease_logFC
      )
    ) %>%
    dplyr::group_by(
      Gene_name
    ) %>%
    dplyr::slice_max(
      order_by = abs(
        Disease_logFC
      ),
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup()

  reversal_contrasts <- c(
    "LCD_vs_PCOS",
    "GFC_vs_PCOS",
    "Combine_vs_PCOS"
  )

  reversal_rows <- list()
  kr <- 1

  set.seed(
    20260930
  )

  for (nm in reversal_contrasts) {

    ttb <- de_list[[nm]] %>%
      dplyr::filter(
        !is.na(Gene_name),
        Gene_name != "",
        is.finite(logFC)
      ) %>%
      dplyr::group_by(
        Gene_name
      ) %>%
      dplyr::slice_max(
        order_by = abs(stat),
        n = 1,
        with_ties = FALSE
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(
        Gene_name,
        Treatment_logFC = logFC
      )

    zz <- ref_sig %>%
      dplyr::inner_join(
        ttb,
        by = "Gene_name"
      )

    if (nrow(zz) < 30) {
      next
    }

    rho <- suppressWarnings(
      cor(
        zz$Disease_logFC,
        zz$Treatment_logFC,
        method = "spearman",
        use = "complete.obs"
      )
    )

    reversal_score <- -rho

    fraction_opposite <- mean(
      sign(
        zz$Disease_logFC
      ) ==
        -sign(
          zz$Treatment_logFC
        ),
      na.rm = TRUE
    )

    boot_scores <- replicate(
      1000,
      {
        idx <- sample(
          seq_len(
            nrow(zz)
          ),
          replace = TRUE
        )

        -suppressWarnings(
          cor(
            zz$Disease_logFC[idx],
            zz$Treatment_logFC[idx],
            method = "spearman",
            use = "complete.obs"
          )
        )
      }
    )

    ci <- stats::quantile(
      boot_scores,
      probs = c(
        0.025,
        0.975
      ),
      na.rm = TRUE
    )

    reversal_rows[[kr]] <- data.frame(
      Contrast = nm,
      N_genes = nrow(zz),
      ReversalScore = reversal_score,
      CI_low = unname(
        ci[1]
      ),
      CI_high = unname(
        ci[2]
      ),
      FractionOpposite =
        fraction_opposite,
      stringsAsFactors = FALSE
    )

    kr <- kr + 1
  }

  reversal_summary <- dplyr::bind_rows(
    reversal_rows
  )

  write_csv(
    reversal_summary,
    file.path(
      TABLE_DIR,
      "Fig4F_PCOS_DiseaseSignature_Reversal.csv"
    )
  )

  if (nrow(reversal_summary) > 0) {

    reversal_summary <- reversal_summary %>%
      dplyr::mutate(
        Treatment = dplyr::recode(
          Contrast,
          "LCD_vs_PCOS" = "LCD",
          "GFC_vs_PCOS" = "GFC",
          "Combine_vs_PCOS" = "Combine"
        ),
        Treatment = factor(
          Treatment,
          levels = c(
            "LCD",
            "GFC",
            "Combine"
          )
        )
      )

    p_F_main <- ggplot(
      reversal_summary,
      aes(
        x = Treatment,
        y = ReversalScore,
        colour = Treatment
      )
    ) +

      geom_hline(
        yintercept = 0,
        colour = "grey70",
        linewidth = 0.45
      ) +

      geom_errorbar(
        aes(
          ymin = CI_low,
          ymax = CI_high
        ),
        width = 0.08,
        linewidth = 0.7
      ) +

      geom_point(
        size = 3.5
      ) +

      scale_colour_manual(
        values = GROUP_COLS[
          c(
            "LCD",
            "GFC",
            "Combine"
          )
        ],
        guide = "none"
      ) +

      labs(
        x = NULL,
        y = "PCOS disease-signature reversal score",
        caption = "Positive values indicate treatment effects opposing an external PCOS-vs-healthy transcriptional signature."
      ) +

      theme_classic(
        base_size = 14
      ) +

      theme(
        axis.text = element_text(
          colour = "black"
        ),
        axis.text.x = element_text(
          face = "bold"
        ),
        axis.title.y = element_text(
          face = "bold"
        ),
        plot.caption = element_text(
          size = 8.8,
          hjust = 0
        )
      )

    save_plot(
      p_F_main,
      "Fig4F_PCOS_DiseaseSignature_Reversal",
      MAIN_DIR,
      width = 5.6,
      height = 4.8
    )
  }

} else {

  writeLines(
    c(
      "Fig4F disease-signature reversal was NOT generated.",
      "",
      "Reason: no healthy-control/reference PCOS signature is available in the treatment-only RNA-seq dataset.",
      "",
      "To enable the analysis, create:",
      PCOS_REFERENCE_SIGNATURE,
      "",
      "Required columns:",
      "Gene_name",
      "Disease_logFC",
      "",
      "Disease_logFC must be PCOS - Healthy/Control.",
      "Positive = up in PCOS; negative = down in PCOS.",
      "",
      "Until such a reference is provided, 'Combine is farther from PCOS' must NOT be interpreted as 'Combine is closer to healthy'."
    ),
    file.path(
      QC_DIR,
      "Fig4F_REFERENCE_SIGNATURE_REQUIRED.txt"
    )
  )
}

# ======================================================================
# 15. Supplementary | direct treatment GO-GSEA

# ======================================================================

direct_contrasts <- c(
  "LCD_vs_GFC",
  "Combine_vs_LCD",
  "Combine_vs_GFC"
)

if (nrow(gsea_all) > 0) {

  gsea_direct_selected <- select_representative_gsea_terms(
    gsea_all,
    direct_contrasts,
    max_terms = 12,
    n_per_theme = 2
  )

  write_csv(
    gsea_direct_selected,
    file.path(
      TABLE_DIR,
      "Fig4D_GO_GSEA_SelectedTerms.csv"
    )
  )

  p_D <- make_gsea_heatmap(
    all_df = gsea_all,
    selected_df = gsea_direct_selected,
    contrasts = direct_contrasts,
    contrast_labels = c(
      "LCD_vs_GFC" = "LCD\nvs GFC",
      "Combine_vs_LCD" = "Combine\nvs LCD",
      "Combine_vs_GFC" = "Combine\nvs GFC"
    ),
    filename =
      "FigS4D_GO_GSEA_DirectTreatmentDifferences",
    width = 8.9,
    dir = SUPP_DIR
  )
}

# ======================================================================
# 16. Supplementary | GO-GSEA: LCD x GFC interaction
#
# The old ORA identified "rhythmic process". GO-GSEA now asks whether
# this signal is supported by the full ranked transcriptome rather than
# only threshold-passing genes.
# ======================================================================

gsea_int_selected <- data.frame()

if (nrow(gsea_all) > 0) {

  gsea_int <- gsea_all %>%
    dplyr::filter(
      Contrast ==
        "Interaction_LCDxGFC",
      !is.na(p.adjust),
      p.adjust < 0.10
    ) %>%
    dplyr::arrange(
      p.adjust,
      dplyr::desc(abs(NES))
    )

  if (nrow(gsea_int) > 0) {

    # Prioritize interpretable themes but preserve top unbiased results.
    gsea_int_selected <- select_representative_gsea_terms(
      gsea_all,
      "Interaction_LCDxGFC",
      max_terms = 12,
      n_per_theme = 3
    )

    int_plot_df <- gsea_int %>%
      dplyr::filter(
        ID %in%
          gsea_int_selected$ID
      ) %>%
      dplyr::mutate(
        TermLabel = stringr::str_wrap(
          Description,
          width = 40
        ),
        Direction = ifelse(
          NES >= 0,
          "Positive interaction",
          "Negative interaction"
        ),
        Strength = pmin(
          -log10(
            pmax(
              p.adjust,
              1e-300
            )
          ),
          12
        )
      ) %>%
      dplyr::arrange(
        NES
      )

    int_plot_df$TermLabel <- factor(
      int_plot_df$TermLabel,
      levels = int_plot_df$TermLabel
    )

    write_csv(
      int_plot_df,
      file.path(
        TABLE_DIR,
        "Fig4E_GO_GSEA_Interaction_SelectedTerms.csv"
      )
    )

    p_E <- ggplot(
      int_plot_df,
      aes(
        x = NES,
        y = TermLabel,
        size = Strength,
        fill = Direction
      )
    ) +

      geom_vline(
        xintercept = 0,
        colour = "grey70",
        linewidth = 0.45
      ) +

      geom_point(
        shape = 21,
        colour = "grey25",
        stroke = 0.35,
        alpha = 0.92
      ) +

      scale_fill_manual(
        values = c(
          "Positive interaction" = "#C94C4C",
          "Negative interaction" = "#3B73B9"
        )
      ) +

      scale_size_continuous(
        range = c(
          3.0,
          7.0
        ),
        name = expression(
          -log[10](FDR)
        )
      ) +

      labs(
        x = "Normalized enrichment score (NES)",
        y = NULL,
        fill = NULL,
        caption =
          "Interaction = Combine − LCD − GFC + PCOS"
      ) +

      theme_classic(
        base_size = 14
      ) +

      theme(
        legend.position = "right",
        axis.text.x = element_text(
          colour = "black"
        ),
        axis.text.y = element_text(
          size = 9.7,
          colour = "black",
          lineheight = 0.95
        ),
        axis.title.x = element_text(
          face = "bold"
        ),
        plot.caption = element_text(
          size = 9,
          hjust = 0
        ),
        plot.margin = margin(
          8,
          10,
          8,
          8
        )
      )

    save_plot(
      p_E,
      "FigS4E_GO_GSEA_LCDxGFC_Interaction",
      SUPP_DIR,
      width = 8.2,
      height = max(
        4.8,
        2.8 +
          0.38 *
          nrow(
            int_plot_df
          )
      )
    )
  }
}

# ----------------------------------------------------------------------
# 16.1 Supplementary ORA sensitivity analysis
#
# Require at least 3 genes supporting a GO term for plotting. This prevents
# one-gene terms (e.g. isolated neurotransmitter-transport annotations) from
# being presented as strong pathway-level findings.
# ----------------------------------------------------------------------

background_symbols <- annotation_df$Gene_name[
  match(
    rownames(norm_expr),
    annotation_df$Gene_id
  )
]

background_symbols <- unique(
  background_symbols[
    !is.na(background_symbols) &
      background_symbols != ""
  ]
)

select_genes_for_ora <- function(
  tt,
  direction
) {

  if (direction == "Up") {

    x <- tt %>%
      dplyr::filter(
        !is.na(adj.P.Val),
        adj.P.Val < FDR_CUTOFF,
        logFC >= LOGFC_CUTOFF,
        !is.na(Gene_name),
        Gene_name != ""
      )

  } else {

    x <- tt %>%
      dplyr::filter(
        !is.na(adj.P.Val),
        adj.P.Val < FDR_CUTOFF,
        logFC <= -LOGFC_CUTOFF,
        !is.na(Gene_name),
        Gene_name != ""
      )
  }

  unique(
    x$Gene_name
  )
}

ora_list <- list()

for (nm in ALL_ANALYSIS_CONTRASTS) {

  for (direction in c(
    "Up",
    "Down"
  )) {

    genes <- select_genes_for_ora(
      de_list[[nm]],
      direction
    )

    if (length(genes) >= 5) {

      ego <- tryCatch(
        clusterProfiler::enrichGO(
          gene = genes,
          universe = background_symbols,
          OrgDb =
            org.Mm.eg.db::org.Mm.eg.db,
          keyType = "SYMBOL",
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = 1,
          qvalueCutoff = 1,
          readable = TRUE
        ),
        error = function(e) NULL
      )

      if (!is.null(ego)) {

        tb <- as.data.frame(
          ego
        )

        if (nrow(tb) > 0) {

          tb <- tb %>%
            dplyr::mutate(
              Contrast = nm,
              Direction = direction,
              Theme = assign_go_theme(
                Description
              )
            )

          ora_list[[
            paste(
              nm,
              direction,
              sep = "__"
            )
          ]] <- tb

          write_csv(
            tb,
            file.path(
              GO_ORA_DIR,
              paste0(
                "ORA_GO_BP_",
                nm,
                "_",
                direction,
                ".csv"
              )
            )
          )
        }
      }
    }
  }
}

if (length(ora_list) > 0) {

  ora_all <- bind_rows(
    ora_list
  )

  write_csv(
    ora_all,
    file.path(
      TABLE_DIR,
      "GO_BP_ORA_all_contrasts_supplementary.csv"
    )
  )
}

# 13. Strict Combine-specific up/down gene sets
#
# Strict down:
#   Combine vs LCD < -0.58, FDR < 0.10
#   AND
#   Combine vs GFC < -0.58, FDR < 0.10
#
# Strict up is analogous.
# ======================================================================

cl <- de_list$Combine_vs_LCD %>%
  dplyr::select(
    Gene_id,
    Gene_name,
    Gene_description,
    logFC_Combine_LCD = logFC,
    FDR_Combine_LCD = adj.P.Val
  )

cg <- de_list$Combine_vs_GFC %>%
  dplyr::select(
    Gene_id,
    logFC_Combine_GFC = logFC,
    FDR_Combine_GFC = adj.P.Val
  )

int_small <- interaction_tt %>%
  dplyr::select(
    Gene_id,
    Interaction_logFC = logFC,
    Interaction_FDR = adj.P.Val
  )

combine_specific <- cl %>%
  inner_join(
    cg,
    by = "Gene_id"
  ) %>%
  left_join(
    int_small,
    by = "Gene_id"
  )

combine_specific_down <- combine_specific %>%
  filter(
    FDR_Combine_LCD < FDR_CUTOFF,
    FDR_Combine_GFC < FDR_CUTOFF,
    logFC_Combine_LCD <=
      -LOGFC_CUTOFF,
    logFC_Combine_GFC <=
      -LOGFC_CUTOFF
  ) %>%
  arrange(
    pmax(
      FDR_Combine_LCD,
      FDR_Combine_GFC
    ),
    Interaction_FDR
  )

combine_specific_up <- combine_specific %>%
  filter(
    FDR_Combine_LCD < FDR_CUTOFF,
    FDR_Combine_GFC < FDR_CUTOFF,
    logFC_Combine_LCD >=
      LOGFC_CUTOFF,
    logFC_Combine_GFC >=
      LOGFC_CUTOFF
  ) %>%
  arrange(
    pmax(
      FDR_Combine_LCD,
      FDR_Combine_GFC
    ),
    Interaction_FDR
  )

write_csv(
  combine_specific_down,
  file.path(
    TABLE_DIR,
    "Combine_specific_downregulated_genes.csv"
  )
)

write_csv(
  combine_specific_up,
  file.path(
    TABLE_DIR,
    "Combine_specific_upregulated_genes.csv"
  )
)

# Formal GO for strict Combine-specific sets (no nominal fallback)
run_strict_go <- function(
  genes,
  label
) {

  genes <- unique(
    genes[
      !is.na(genes) &
        genes != ""
    ]
  )

  if (length(genes) < 5) {
    return(NULL)
  }

  ego <- tryCatch(
    clusterProfiler::enrichGO(
      gene = genes,
      universe = background_symbols,
      OrgDb = org.Mm.eg.db::org.Mm.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(ego)) {
    return(NULL)
  }

  tb <- as.data.frame(
    ego
  )

  if (nrow(tb) == 0) {
    return(NULL)
  }

  tb$GeneSet <- label
  tb
}

go_combine_down <- run_strict_go(
  combine_specific_down$Gene_name,
  "Combine-specific down"
)

go_combine_up <- run_strict_go(
  combine_specific_up$Gene_name,
  "Combine-specific up"
)

combine_go_list <- Filter(
  Negate(is.null),
  list(
    go_combine_down,
    go_combine_up
  )
)

if (length(combine_go_list) > 0) {

  combine_go <- bind_rows(
    combine_go_list
  )

  write_csv(
    combine_go,
    file.path(
      TABLE_DIR,
      "GO_CombineSpecific_Strict.csv"
    )
  )

  combine_go_plot <- combine_go %>%
    filter(
      p.adjust < 0.10,
      Count >= 3
    ) %>%
    group_by(
      GeneSet
    ) %>%
    arrange(
      p.adjust,
      desc(Count),
      .by_group = TRUE
    ) %>%
    slice_head(
      n = 10
    ) %>%
    ungroup() %>%
    mutate(
      TermLabel = stringr::str_wrap(
        Description,
        width = 40
      ),
      Strength = -log10(
        pmax(
          p.adjust,
          1e-300
        )
      )
    )

  if (nrow(combine_go_plot) == 0) {

    writeLines(
      "No robust strict Combine-specific GO-BP term remained after requiring FDR < 0.10 and Count >= 3. The previous Count=1 GO terms should not be used as pathway-level evidence.",
      file.path(
        QC_DIR,
        "StrictCombineSpecific_GO_NOT_ROBUST.txt"
      )
    )
  }

  if (nrow(combine_go_plot) > 0) {

    p_S4D <- ggplot(
      combine_go_plot,
      aes(
        x = Strength,
        y = TermLabel,
        size = Count,
        fill = GeneSet
      )
    ) +

      geom_point(
        shape = 21,
        colour = "grey25",
        alpha = 0.90
      ) +

      facet_wrap(
        ~ GeneSet,
        ncol = 2,
        scales = "free_y"
      ) +

      scale_fill_manual(
        values = c(
          "Combine-specific up" = "#C94C4C",
          "Combine-specific down" = "#3B73B9"
        )
      ) +

      labs(
        x = expression(
          -log[10](FDR)
        ),
        y = NULL,
        size = "Genes",
        fill = NULL
      ) +

      theme_classic(
        base_size = 14
      ) +

      theme(
        strip.background = element_blank(),
        strip.text = element_text(
          face = "bold"
        ),
        axis.text.y = element_text(
          size = 9.5,
          colour = "black"
        )
      )

    save_plot(
      p_S4D,
      "FigS4H_GO_StrictCombineSpecific",
      SUPP_DIR,
      width = 9.8,
      height = 6.5
    )
  }
}

# ======================================================================
# 18. Circadian interaction representative genes
#
# The current data show a coherent circadian-clock interaction program.
# Therefore the main representative panel prioritizes canonical clock genes
# that are also supported by the actual interaction statistics.
#
# No gene is forced into the figure if it lacks interaction evidence.
# ======================================================================

interaction_gsea_core_genes <- character(0)

if (
  exists("gsea_int_selected") &&
    nrow(gsea_int_selected) > 0 &&
    nrow(gsea_all) > 0
) {

  core_rows <- gsea_all %>%
    dplyr::filter(
      Contrast ==
        "Interaction_LCDxGFC",
      ID %in%
        gsea_int_selected$ID
    )

  if (
    "core_enrichment" %in%
      names(core_rows)
  ) {

    interaction_gsea_core_genes <- unique(
      unlist(
        strsplit(
          paste(
            core_rows$core_enrichment,
            collapse = "/"
          ),
          "/",
          fixed = TRUE
        )
      )
    )
  }
}

clock_priority <- c(
  "Bmal1",
  "Arntl",
  "Npas2",
  "Per3",
  "Per2",
  "Dbp",
  "Nr1d2",
  "Tef",
  "Hlf",
  "Clock",
  "Cry1",
  "Cry2",
  "Rbm3"
)

secondary_story_genes <- c(
  "Mtor",
  "Rptor",
  "Becn1",
  "Atg5",
  "Atg7",
  "Sqstm1",
  "Star",
  "Cyp11a1",
  "Cyp17a1",
  "Cyp19a1",
  "Lhcgr",
  "Fshr",
  "Nfkb1",
  "Rela",
  "Nfkbia",
  "Nfe2l2",
  "Hmox1",
  "Nqo1",
  "Dcn",
  "Enpp2"
)

interaction_ranked <- interaction_tt %>%
  dplyr::mutate(
    PrimarySig =
      !is.na(adj.P.Val) &
      adj.P.Val < FDR_CUTOFF &
      abs(logFC) >= LOGFC_CUTOFF,

    ClockRank = match(
      Gene_name,
      clock_priority
    ),

    ClockRank = ifelse(
      is.na(ClockRank),
      999,
      ClockRank
    ),

    IsClock =
      Gene_name %in%
      clock_priority,

    InInteractionGO =
      Gene_name %in%
      interaction_gsea_core_genes,

    IsSecondaryStory =
      Gene_name %in%
      secondary_story_genes,

    SelectionTier = dplyr::case_when(

      PrimarySig &
        IsClock &
        InInteractionGO ~ 1,

      PrimarySig &
        IsClock ~ 2,

      PrimarySig &
        InInteractionGO ~ 3,

      PrimarySig &
        IsSecondaryStory ~ 4,

      PrimarySig ~ 5,

      TRUE ~ 9
    ),

    SelectionReason = dplyr::case_when(

      SelectionTier == 1 ~
        "Significant interaction + canonical clock gene + interaction GO core",

      SelectionTier == 2 ~
        "Significant interaction + canonical circadian-clock gene",

      SelectionTier == 3 ~
        "Significant interaction + interaction GO core gene",

      SelectionTier == 4 ~
        "Significant interaction + PCOS/ovary mechanism gene",

      SelectionTier == 5 ~
        "Strong remaining interaction gene",

      TRUE ~
        "Lower priority"
    )
  ) %>%
  dplyr::arrange(
    SelectionTier,
    ClockRank,
    adj.P.Val,
    dplyr::desc(abs(logFC))
  )

representative_genes <- interaction_ranked %>%
  dplyr::filter(
    SelectionTier < 9,
    !is.na(Gene_name),
    Gene_name != ""
  ) %>%
  dplyr::distinct(
    Gene_name,
    .keep_all = TRUE
  ) %>%
  dplyr::slice_head(
    n = 6
  )

if (nrow(representative_genes) < 6) {

  fallback <- interaction_ranked %>%
    dplyr::filter(
      !is.na(Gene_name),
      Gene_name != "",
      !Gene_name %in%
        representative_genes$Gene_name
    ) %>%
    dplyr::distinct(
      Gene_name,
      .keep_all = TRUE
    ) %>%
    dplyr::slice_head(
      n =
        6 -
        nrow(representative_genes)
    )

  representative_genes <- dplyr::bind_rows(
    representative_genes,
    fallback
  )
}

representative_genes <- representative_genes %>%
  dplyr::slice_head(
    n = 6
  )

write_csv(
  representative_genes,
  file.path(
    TABLE_DIR,
    "Fig4F_CircadianInteractionGenes_Selected.csv"
  )
)

# 19. Supplementary | interaction volcano
# ======================================================================

interaction_plot_df <- interaction_ranked %>%
  mutate(
    Pplot = pmax(
      P.Value,
      1e-300
    ),
    Direction = case_when(
      PrimarySig &
        logFC > 0 ~
        "Positive interaction",
      PrimarySig &
        logFC < 0 ~
        "Negative interaction",
      TRUE ~
        "Not significant"
    )
  )

label_interaction <- interaction_plot_df %>%
  dplyr::filter(
    Gene_name %in%
      representative_genes$Gene_name |
      (
        PrimarySig &
          (
            IsClock |
              InInteractionGO
          )
      )
  ) %>%
  arrange(
    SelectionTier,
    adj.P.Val,
    desc(abs(logFC))
  ) %>%
  distinct(
    Gene_name,
    .keep_all = TRUE
  ) %>%
  slice_head(
    n = 12
  )

volcano_zoom <- as.numeric(
  stats::quantile(
    abs(
      interaction_plot_df$logFC[
        is.finite(
          interaction_plot_df$logFC
        )
      ]
    ),
    probs = 0.99,
    na.rm = TRUE
  )
)

if (!is.finite(volcano_zoom)) {
  volcano_zoom <- 6
}

volcano_zoom <- min(
  max(
    volcano_zoom,
    5
  ),
  10
)

p_E <- ggplot(
  interaction_plot_df,
  aes(
    x = logFC,
    y = -log10(Pplot)
  )
) +

  geom_vline(
    xintercept = 0,
    colour = "grey70",
    linewidth = 0.45
  ) +

  geom_point(
    data = interaction_plot_df %>%
      filter(
        Direction ==
          "Not significant"
      ),
    colour = "grey72",
    size = 1.0,
    alpha = 0.22
  ) +

  geom_point(
    data = interaction_plot_df %>%
      filter(
        Direction !=
          "Not significant"
      ),
    aes(
      colour = Direction
    ),
    size = 1.8,
    alpha = 0.75
  ) +

  ggrepel::geom_text_repel(
    data = label_interaction,
    aes(
      label = Gene_name,
      colour = Direction
    ),
    size = 3.25,
    max.overlaps = Inf,
    box.padding = 0.55,
    point.padding = 0.30,
    min.segment.length = 0,
    segment.colour = "grey55",
    segment.size = 0.35,
    seed = 20260930,
    show.legend = FALSE
  ) +

  scale_colour_manual(
    values = c(
      "Positive interaction" = "#C94C4C",
      "Negative interaction" = "#3B73B9"
    )
  ) +

  coord_cartesian(
    xlim = c(
      -volcano_zoom,
      volcano_zoom
    )
  ) +

  labs(
    x = "LCD × GFC interaction effect\n(Combine − LCD − GFC + PCOS)",
    y = expression(
      -log[10](italic(P))
    ),
    colour = NULL
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    legend.position = "top",
    axis.text = element_text(
      colour = "black"
    ),
    axis.title = element_text(
      face = "bold"
    )
  )

save_plot(
  p_E,
  "FigS4F_LCDxGFC_Interaction_Volcano",
  SUPP_DIR,
  width = 6.5,
  height = 5.3
)


p_S4G <- ggplot(
  interaction_plot_df,
  aes(
    x = logFC,
    y = -log10(Pplot)
  )
) +

  geom_vline(
    xintercept = 0,
    colour = "grey70",
    linewidth = 0.45
  ) +

  geom_point(
    data = interaction_plot_df %>%
      dplyr::filter(
        Direction ==
          "Not significant"
      ),
    colour = "grey72",
    size = 0.9,
    alpha = 0.20
  ) +

  geom_point(
    data = interaction_plot_df %>%
      dplyr::filter(
        Direction !=
          "Not significant"
      ),
    aes(
      colour = Direction
    ),
    size = 1.5,
    alpha = 0.70
  ) +

  scale_colour_manual(
    values = c(
      "Positive interaction" = "#C94C4C",
      "Negative interaction" = "#3B73B9"
    )
  ) +

  labs(
    x = "LCD × GFC interaction effect",
    y = expression(
      -log[10](italic(P))
    ),
    colour = NULL
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    legend.position = "top",
    axis.text = element_text(
      colour = "black"
    ),
    axis.title = element_text(
      face = "bold"
    )
  )

save_plot(
  p_S4G,
  "FigS4G_LCDxGFC_Interaction_Volcano_FullRange",
  SUPP_DIR,
  width = 6.3,
  height = 5.0
)

# ======================================================================
# 20. Fig4G | Circadian interaction genes (mechanism, NOT efficacy)
# ======================================================================

if (length(tpm_sample_cols) > 0) {

  tpm_mat <- safe_numeric_matrix(
    tpm_df,
    tpm_sample_cols,
    tpm_df$Gene_id
  )

  tpm_mat[
    is.na(tpm_mat)
  ] <- 0

  expr_rows <- list()
  ke <- 1

  for (i in seq_len(
    nrow(representative_genes)
  )) {

    gid <- representative_genes$Gene_id[i]
    gname <- representative_genes$Gene_name[i]

    if (gid %in%
        rownames(tpm_mat)) {

      vals <- log2(
        tpm_mat[
          gid,
          ,
          drop = TRUE
        ] + 1
      )

      expr_rows[[ke]] <- data.frame(
        Sample = names(vals),
        Gene = gname,
        Expression = as.numeric(vals),
        stringsAsFactors = FALSE
      )

      ke <- ke + 1
    }
  }

  representative_expr <- bind_rows(
    expr_rows
  )

  y_lab <- "log2(TPM + 1)"

} else {

  expr_rows <- list()
  ke <- 1

  for (i in seq_len(
    nrow(representative_genes)
  )) {

    gid <- representative_genes$Gene_id[i]
    gname <- representative_genes$Gene_name[i]

    if (gid %in%
        rownames(norm_expr)) {

      vals <- norm_expr[
        gid,
        ,
        drop = TRUE
      ]

      expr_rows[[ke]] <- data.frame(
        Sample = names(vals),
        Gene = gname,
        Expression = as.numeric(vals),
        stringsAsFactors = FALSE
      )

      ke <- ke + 1
    }
  }

  representative_expr <- bind_rows(
    expr_rows
  )

  y_lab <- "Voom log-expression"
}

representative_expr <- representative_expr %>%
  mutate(
    Group = factor(
      group_from_sample(Sample),
      levels = GROUP_LEVELS
    ),
    GeneDisplay = stringr::str_wrap(
      Gene,
      width = 18
    ),
    GeneDisplay = factor(
      GeneDisplay,
      levels = unique(
        stringr::str_wrap(
          representative_genes$Gene_name,
          width = 18
        )
      )
    )
  )

p_F <- ggplot(
  representative_expr,
  aes(
    x = Group,
    y = Expression,
    fill = Group
  )
) +

  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    alpha = 0.28,
    colour = "grey35"
  ) +

  geom_jitter(
    aes(
      colour = Group
    ),
    width = 0.08,
    size = 2.2,
    alpha = 0.90
  ) +

  facet_wrap(
    ~ GeneDisplay,
    scales = "free_y",
    ncol = 3,
    drop = FALSE
  ) +

  scale_fill_manual(
    values = GROUP_COLS,
    guide = "none"
  ) +

  scale_colour_manual(
    values = GROUP_COLS,
    guide = "none"
  ) +

  labs(
    x = NULL,
    y = y_lab,
    caption = "Non-additive circadian interaction. These genes describe a combination-specific regulatory state and are not evidence by themselves of superior PCOS normalization."
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12.3
    ),
    axis.text.x = element_text(
      angle = 40,
      hjust = 1,
      colour = "black",
      size = 11.0
    ),
    axis.text.y = element_text(
      size = 11.0,
      colour = "black"
    ),
    axis.title.y = element_text(
      face = "bold"
    ),
    plot.caption = element_text(
      size = 9.5,
      hjust = 0
    )
  )

save_plot(
  p_F,
  "Fig4F_CircadianInteraction_Genes",
  MAIN_DIR,
  width = 7.8,
  height = 6.1
)

# ======================================================================
# 17. Supplementary FigS4A | PCA
# ======================================================================

gene_var <- apply(
  norm_expr,
  1,
  var,
  na.rm = TRUE
)

top_var_ids <- names(
  sort(
    gene_var,
    decreasing = TRUE
  )
)[
  seq_len(
    min(
      1000,
      length(gene_var)
    )
  )
]

pca <- prcomp(
  t(
    norm_expr[
      top_var_ids,
      ,
      drop = FALSE
    ]
  ),
  center = TRUE,
  scale. = FALSE
)

var_exp <- 100 *
  pca$sdev^2 /
  sum(
    pca$sdev^2
  )

pca_df <- data.frame(
  Sample = rownames(
    pca$x
  ),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
) %>%
  mutate(
    Group = factor(
      group_from_sample(Sample),
      levels = GROUP_LEVELS
    )
  )

p_S4A <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    colour = Group,
    shape = Group
  )
) +

  stat_ellipse(
    aes(group = Group),
    type = "t",
    linetype = 2,
    linewidth = 0.65,
    show.legend = FALSE
  ) +

  geom_point(
    size = 3.0,
    alpha = 0.90
  ) +

  scale_colour_manual(
    values = GROUP_COLS,
    breaks = GROUP_LEVELS
  ) +

  scale_shape_manual(
    values = GROUP_SHAPES,
    breaks = GROUP_LEVELS
  ) +

  labs(
    x = sprintf(
      "PC1 (%.1f%%)",
      var_exp[1]
    ),
    y = sprintf(
      "PC2 (%.1f%%)",
      var_exp[2]
    ),
    colour = NULL,
    shape = NULL
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    legend.position = "top",
    axis.text = element_text(
      colour = "black"
    ),
    axis.title = element_text(
      face = "bold"
    )
  )

save_plot(
  p_S4A,
  "FigS4A_OvaryTranscriptome_PCA",
  SUPP_DIR,
  width = 5.8,
  height = 4.8
)

# ======================================================================
# 18. Supplementary FigS4B | Top DEG heatmap
# ======================================================================

top_deg_ids <- unique(
  unlist(
    lapply(
      PAIRWISE_CONTRASTS,
      function(nm) {

        de_list[[nm]] %>%
          arrange(
            adj.P.Val,
            desc(abs(logFC))
          ) %>%
          slice_head(
            n = 12
          ) %>%
          pull(
            Gene_id
          )
      }
    ),
    use.names = FALSE
  )
)

top_deg_ids <- intersect(
  top_deg_ids,
  rownames(norm_expr)
)

if (length(top_deg_ids) > 45) {
  top_deg_ids <- top_deg_ids[
    seq_len(45)
  ]
}

heat <- norm_expr[
  top_deg_ids,
  ,
  drop = FALSE
]

heat_z <- t(
  scale(
    t(heat)
  )
)

heat_z[
  !is.finite(heat_z)
] <- 0

gene_labels <- annotation_df$Gene_name[
  match(
    rownames(heat_z),
    annotation_df$Gene_id
  )
]

gene_labels[
  is.na(gene_labels) |
    gene_labels == ""
] <- rownames(heat_z)[
  is.na(gene_labels) |
    gene_labels == ""
]

rownames(heat_z) <- make.unique(
  gene_labels
)

sample_order_df <- sample_meta %>%
  mutate(
    SampleNumber = sample_number(Sample),
    DisplaySample = sample_display(Sample)
  ) %>%
  arrange(
    Group,
    SampleNumber
  )

heat_z <- heat_z[
  ,
  sample_order_df$Sample,
  drop = FALSE
]

colnames(heat_z) <- sample_order_df$DisplaySample

anno_col <- data.frame(
  Group = factor(
    as.character(
      sample_order_df$Group
    ),
    levels = GROUP_LEVELS
  ),
  row.names =
    sample_order_df$DisplaySample
)

anno_colors <- list(
  Group = GROUP_COLS
)

gaps_col <- cumsum(
  as.integer(
    table(
      factor(
        sample_order_df$Group,
        levels = GROUP_LEVELS
      )
    )
  )
)

gaps_col <- gaps_col[
  -length(gaps_col)
]

heat_width <- max(
  9.5,
  4.0 +
    0.42 *
    ncol(heat_z)
)

heat_height <- max(
  7.0,
  2.6 +
    0.16 *
    nrow(heat_z)
)

pdf(
  file.path(tempdir(),"FigS4B_TopDEG_Heatmap.pdf"),
  width = heat_width,
  height = heat_height
)

pheatmap(
  heat_z,
  annotation_col = anno_col,
  annotation_colors = anno_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 11,
  fontsize_row = 9,
  fontsize_col = 9.5,
  angle_col = 45,
  border_color = NA,
  main = ""
)

dev.off()

png(
  file.path(tempdir(),"FigS4B_TopDEG_Heatmap.png"),
  width = heat_width,
  height = heat_height,
  units = "in",
  res = 400
)

pheatmap(
  heat_z,
  annotation_col = anno_col,
  annotation_colors = anno_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 11,
  fontsize_row = 9,
  fontsize_col = 9.5,
  angle_col = 45,
  border_color = NA,
  main = ""
)

dev.off()

# ======================================================================
# 19. Curated ovarian modules for sample-level trends and integration
# ======================================================================

ovarian_modules <- list(

  "Steroidogenesis / androgen" = c(
    "Star",
    "Cyp11a1",
    "Cyp17a1",
    "Hsd3b1",
    "Hsd3b2",
    "Cyp19a1",
    "Lhcgr",
    "Fshr",
    "Nr5a1",
    "Scarb1",
    "Srd5a1",
    "Srd5a2"
  ),

  "Folliculogenesis / ovulation" = c(
    "Fshr",
    "Lhcgr",
    "Amh",
    "Amhr2",
    "Inha",
    "Inhba",
    "Inhbb",
    "Bmp15",
    "Gdf9",
    "Ptgs2",
    "Has2",
    "Pgr",
    "Areg",
    "Ereg",
    "Adamts1"
  ),

  "PI3K-AKT-mTOR / insulin" = c(
    "Pik3ca",
    "Pik3cb",
    "Pik3r1",
    "Akt1",
    "Akt2",
    "Mtor",
    "Rptor",
    "Rictor",
    "Tsc1",
    "Tsc2",
    "Pten",
    "Irs1",
    "Irs2",
    "Insr",
    "Igf1r",
    "Rps6kb1",
    "Eif4ebp1"
  ),

  "Inflammation / NF-kB" = c(
    "Tlr4",
    "Myd88",
    "Nfkb1",
    "Rela",
    "Ikbkb",
    "Nfkbia",
    "Tnf",
    "Il6",
    "Il1b",
    "Ccl2",
    "Nlrp3",
    "Ptgs2",
    "Stat3"
  ),

  "Nrf2 / antioxidant" = c(
    "Nfe2l2",
    "Keap1",
    "Hmox1",
    "Nqo1",
    "Gclc",
    "Gclm",
    "Sod1",
    "Sod2",
    "Cat",
    "Gpx1",
    "Gpx4",
    "Gsr",
    "Txnrd1"
  ),

  "Autophagy / apoptosis" = c(
    "Becn1",
    "Atg5",
    "Atg7",
    "Map1lc3b",
    "Sqstm1",
    "Ulk1",
    "Bax",
    "Bcl2",
    "Casp3",
    "Casp9",
    "Bak1",
    "Bcl2l1"
  ),

  "Glycolysis / energy" = c(
    "Hk2",
    "Pkm",
    "Ldha",
    "Slc2a1",
    "Pdk1",
    "Pfkfb3",
    "Gpi1",
    "Eno1",
    "Aldoa",
    "Pfkm"
  ),

  "Lipid / PPAR / mitochondria" = c(
    "Ppara",
    "Pparg",
    "Ppargc1a",
    "Cpt1a",
    "Cpt2",
    "Acaca",
    "Fasn",
    "Scd1",
    "Srebf1",
    "Acox1",
    "Acadm",
    "Acadl",
    "Tfam",
    "Cox4i1"
  )
)

# Collapse voom expression to unique gene symbols
filtered_ids <- rownames(norm_expr)

symbol_lookup <- annotation_df$Gene_name[
  match(
    filtered_ids,
    annotation_df$Gene_id
  )
]

symbol_lookup[
  is.na(symbol_lookup) |
    symbol_lookup == ""
] <- filtered_ids[
  is.na(symbol_lookup) |
    symbol_lookup == ""
]

expr_symbol_df <- as.data.frame(
  norm_expr
)

expr_symbol_df$Gene_name <- symbol_lookup

expr_symbol_df <- expr_symbol_df %>%
  group_by(
    Gene_name
  ) %>%
  summarise(
    across(
      all_of(
        colnames(norm_expr)
      ),
      \(x) mean(x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

expr_symbol <- as.matrix(
  expr_symbol_df[
    ,
    colnames(norm_expr),
    drop = FALSE
  ]
)

rownames(expr_symbol) <- expr_symbol_df$Gene_name

score_one_module <- function(
  genes,
  expr
) {

  present <- intersect(
    genes,
    rownames(expr)
  )

  if (length(present) < 2) {
    return(
      rep(
        NA_real_,
        ncol(expr)
      )
    )
  }

  z <- t(
    scale(
      t(
        expr[
          present,
          ,
          drop = FALSE
        ]
      )
    )
  )

  z[
    !is.finite(z)
  ] <- 0

  colMeans(
    z,
    na.rm = TRUE
  )
}

module_score_mat <- do.call(
  rbind,
  lapply(
    ovarian_modules,
    score_one_module,
    expr = expr_symbol
  )
)

rownames(module_score_mat) <- names(
  ovarian_modules
)

colnames(module_score_mat) <- colnames(
  expr_symbol
)

module_score_long <- as.data.frame(
  t(
    module_score_mat
  )
) %>%
  rownames_to_column(
    "Sample"
  ) %>%
  mutate(
    Group = factor(
      group_from_sample(Sample),
      levels = GROUP_LEVELS
    )
  ) %>%
  pivot_longer(
    cols = -c(
      Sample,
      Group
    ),
    names_to = "Module",
    values_to = "Score"
  )

write_csv(
  module_score_long,
  file.path(
    TABLE_DIR,
    "OvarianModule_SampleScores.csv"
  )
)

# ======================================================================
# 20. Supplementary FigS4C | Module scores
# ======================================================================

p_S4C <- ggplot(
  module_score_long,
  aes(
    x = Group,
    y = Score,
    fill = Group
  )
) +

  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    alpha = 0.28,
    colour = "grey35"
  ) +

  geom_jitter(
    aes(colour = Group),
    width = 0.08,
    size = 2.0,
    alpha = 0.90
  ) +

  facet_wrap(
    ~ Module,
    scales = "free_y",
    ncol = 4,
    labeller = label_wrap_gen(
      width = 23
    )
  ) +

  scale_fill_manual(
    values = GROUP_COLS,
    guide = "none"
  ) +

  scale_colour_manual(
    values = GROUP_COLS,
    guide = "none"
  ) +

  labs(
    x = NULL,
    y = "Module score (mean gene z-score)"
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 9.8
    ),
    axis.text.x = element_text(
      angle = 40,
      hjust = 1,
      size = 9.0,
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    )
  )

save_plot(
  p_S4C,
  "FigS4C_OvarianModule_SampleScores",
  SUPP_DIR,
  width = 12.0,
  height = 5.9
)

# ======================================================================
# 21. Literature basis table
# ======================================================================

literature_basis <- data.frame(

  Topic = c(
    "GFC ovarian PI3K-AKT-mTOR / autophagy",
    "GFC ovarian inflammation / oxidative stress",
    "GFC intestinal homeostasis / insulin resistance",
    "LCD/KD gut-metabolite-androgen axis",
    "Circadian clock and PCOS androgen metabolism",
    "Circadian disruption / granulosa-cell apoptosis"
  ),

  Reference = c(
    "Liu et al. J Ethnopharmacol. 2021. PMID:33460753; DOI:10.1016/j.jep.2021.113821",
    "Ye et al. Gynecol Endocrinol. 2023. PMID:37706434; DOI:10.1080/09513590.2023.2254848",
    "Zhu et al. Front Endocrinol. 2020. PMID:32973686; DOI:10.3389/fendo.2020.00575",
    "Wang et al. Nutrition. 2023. PMID:37441825; DOI:10.1016/j.nut.2023.112127",
    "Johnson et al. Hum Reprod. 2022. PMID:35728080; DOI:10.1093/humrep/deac139",
    "Altered circadian clock as therapeutic target in PCOS. PMID:32119846"
  ),

  stringsAsFactors = FALSE
)

write_csv(
  literature_basis,
  file.path(
    TABLE_DIR,
    "Fig4_Literature_Basis.csv"
  )
)

# ======================================================================
# 22. Optional gut / metabolite -> ovarian-module integration
# ======================================================================

if (all(file.exists(c(FIG3_KEY_GENUS, FIG3_KEY_METAB))) &&
    !file.exists(file.path(FIG3_TABLE_DIR, "genus_CLR_master.csv"))) {
  stop("Figure 3 candidate files exist, but canonical genus_CLR_master.csv is missing. Run patched Figure 3 first.")
}

integration_possible <- all(
  file.exists(
    c(
      FIG3_KEY_GENUS,
      FIG3_KEY_METAB,
      file.path(FIG3_TABLE_DIR, "genus_CLR_master.csv"),
      MICROBIOME_FILE,
      FECAL_FILE,
      TAXONOMY_FILE
    )
  )
)

if (integration_possible) {

  key_genus_fig3 <- readr::read_csv(
    FIG3_KEY_GENUS,
    show_col_types = FALSE
  )

  key_metab_fig3 <- readr::read_csv(
    FIG3_KEY_METAB,
    show_col_types = FALSE
  )

  selected_genera <- unique(
    head(
      key_genus_fig3$Feature,
      6
    )
  )

  selected_metab_ids <- unique(
    head(
      key_metab_fig3$Feature,
      6
    )
  )

  selected_metab_names <- key_metab_fig3$DisplayName[
    match(
      selected_metab_ids,
      key_metab_fig3$Feature
    )
  ]

  # === MAJOR 6: load the canonical Figure 3 genus CLR matrix ===
  MASTER_CLR_FILE <- file.path(FIG3_TABLE_DIR, "genus_CLR_master.csv")
  if (!file.exists(MASTER_CLR_FILE)) {
    stop("Unified CLR matrix missing: ", MASTER_CLR_FILE,
         "\nFirst run Fig3_FINAL_v4_CLR_unified.R to create it.")
  }
  master_clr_df <- readr::read_csv(MASTER_CLR_FILE, show_col_types=FALSE,
                                   name_repair="minimal")
  if (!"Genus" %in% names(master_clr_df) ||
      anyDuplicated(master_clr_df$Genus) || anyDuplicated(names(master_clr_df))) {
    stop("Invalid master CLR: missing/duplicated genus or sample identifiers")
  }
  genus_clr <- as.matrix(master_clr_df[, setdiff(names(master_clr_df), "Genus"), drop=FALSE])
  rownames(genus_clr) <- master_clr_df$Genus
  storage.mode(genus_clr) <- "numeric"
  if (any(!is.finite(genus_clr))) stop("Master CLR contains non-finite values")
  if (any(!grepl("^(MX|DT|GZFL|LH)_[0-9]+$", colnames(genus_clr)))) {
    stop("Unexpected biological sample identifier in master CLR")
  }
  # Every sample should sum to zero over the SAME retained genus reference set.
  if (any(abs(colMeans(genus_clr)) > 1e-7)) stop("CLR column-centering QC failed")

  selected_genera <- intersect(
    selected_genera,
    rownames(genus_clr)
  )

  gut_genus <- genus_clr[
    selected_genera,
    ,
    drop = FALSE
  ]

  # --------------------------
  # Selected fecal metabolites
  # --------------------------

  fecal_df <- readr::read_tsv(
    FECAL_FILE,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )

  fecal_sample_cols <- names(fecal_df)[
    vapply(
      names(fecal_df),
      is_rna_sample,
      logical(1)
    )
  ]

  gut_metab <- NULL

  if ("metab_id" %in%
      names(fecal_df)) {

    fecal_mat_raw <- safe_numeric_matrix(
      fecal_df,
      fecal_sample_cols,
      fecal_df$metab_id
    )

    fecal_mat_raw[
      is.na(fecal_mat_raw)
    ] <- 0

    finite_fecal <- as.numeric(
      fecal_mat_raw[
        is.finite(fecal_mat_raw)
      ]
    )

    if (length(finite_fecal) > 0) {

      q99_fecal <- quantile(
        finite_fecal,
        0.99,
        na.rm = TRUE
      )

      med_fecal <- median(
        finite_fecal,
        na.rm = TRUE
      )

      if (
        q99_fecal >= 50 |
          med_fecal >= 20
      ) {

        shift_fecal <- ifelse(
          min(
            fecal_mat_raw,
            na.rm = TRUE
          ) < 0,
          abs(
            min(
              fecal_mat_raw,
              na.rm = TRUE
            )
          ) + 1,
          1
        )

        fecal_mat_use <- log2(
          fecal_mat_raw +
            shift_fecal
        )

      } else {

        fecal_mat_use <- fecal_mat_raw
      }

      selected_metab_ids2 <- intersect(
        selected_metab_ids,
        rownames(fecal_mat_use)
      )

      if (length(selected_metab_ids2) > 0) {

        gut_metab <- fecal_mat_use[
          selected_metab_ids2,
          ,
          drop = FALSE
        ]

        metab_labels <- selected_metab_names[
          match(
            selected_metab_ids2,
            selected_metab_ids
          )
        ]

        metab_labels[
          is.na(metab_labels) |
            metab_labels == ""
        ] <- selected_metab_ids2[
          is.na(metab_labels) |
            metab_labels == ""
        ]

        rownames(gut_metab) <- metab_labels
      }
    }
  }

  # --------------------------
  # Build feature table
  # --------------------------

  gut_feature_list <- list()

  if (nrow(gut_genus) > 0) {

    genus_feature <- gut_genus
    rownames(genus_feature) <- paste0(
      "GENUS__",
      rownames(genus_feature)
    )

    gut_feature_list[["Genus"]] <- genus_feature
  }

  if (
    !is.null(gut_metab) &&
      nrow(gut_metab) > 0
  ) {

    metab_feature <- gut_metab
    rownames(metab_feature) <- paste0(
      "METAB__",
      rownames(metab_feature)
    )

    gut_feature_list[["Fecal metabolite"]] <- metab_feature
  }

  if (length(gut_feature_list) > 0) {

    all_feature_samples <- Reduce(
      intersect,
      lapply(
        gut_feature_list,
        colnames
      )
    )

    matched_samples <- intersect(
      all_feature_samples,
      colnames(
        module_score_mat
      )
    )

    match_qc <- data.frame(
      Group = GROUP_LEVELS,
      N_matched = vapply(
        GROUP_LEVELS,
        function(g) {

          sum(
            group_from_sample(
              matched_samples
            ) == g
          )
        },
        integer(1)
      )
    )

    write_csv(
      match_qc,
      file.path(
        QC_DIR,
        "Fig4H_Integration_SampleMatching.csv"
      )
    )

    if (length(matched_samples) >= 12) {

      gut_all <- do.call(
        rbind,
        lapply(
          gut_feature_list,
          function(x) {

            x[
              ,
              matched_samples,
              drop = FALSE
            ]
          }
        )
      )

      ovary_mod <- module_score_mat[
        ,
        matched_samples,
        drop = FALSE
      ]

      matched_group <- factor(
        group_from_sample(
          matched_samples
        ),
        levels = GROUP_LEVELS
      )

      residualize_by_group <- function(
        mat,
        grp
      ) {

        out <- matrix(
          NA_real_,
          nrow = nrow(mat),
          ncol = ncol(mat),
          dimnames = dimnames(mat)
        )

        for (i in seq_len(
          nrow(mat)
        )) {

          fit <- lm(
            as.numeric(
              mat[i, ]
            ) ~ grp
          )

          out[i, ] <- residuals(
            fit
          )
        }

        out
      }

      gut_res <- residualize_by_group(
        gut_all,
        matched_group
      )

      mod_res <- residualize_by_group(
        ovary_mod,
        matched_group
      )

      cor_rows <- list()
      kc <- 1

      for (g in rownames(gut_res)) {

        for (m in rownames(mod_res)) {

          ct <- suppressWarnings(
            cor.test(
              as.numeric(
                gut_res[g, ]
              ),
              as.numeric(
                mod_res[m, ]
              ),
              method = "spearman",
              exact = FALSE
            )
          )

          type <- ifelse(
            grepl(
              "^GENUS__",
              g
            ),
            "Genus",
            "Fecal metabolite"
          )

          clean_name <- sub(
            "^(GENUS__|METAB__)",
            "",
            g
          )

          cor_rows[[kc]] <- data.frame(
            FeatureType = type,
            GutFeature = clean_name,
            OvarianModule = m,
            Rho = unname(
              ct$estimate
            ),
            P = ct$p.value,
            stringsAsFactors = FALSE
          )

          kc <- kc + 1
        }
      }

      integration_cor <- bind_rows(
        cor_rows
      ) %>%
        mutate(
          FDR = p.adjust(
            P,
            method = "BH"
          ),
          Sig = dplyr::case_when(
            FDR < 0.05 ~ "**",
            FDR < 0.10 ~ "*",
            P < 0.05 ~ "·",
            TRUE ~ ""
          ),
          FeatureType = factor(
            FeatureType,
            levels = c(
              "Genus",
              "Fecal metabolite"
            )
          ),
          GutDisplay = stringr::str_wrap(
            GutFeature,
            width = 20
          ),
          ModuleDisplay = stringr::str_wrap(
            OvarianModule,
            width = 15
          )
        )

      write_csv(
        integration_cor,
        file.path(
          TABLE_DIR,
          "Fig4H_GutMetabolite_OvaryModule_Correlations.csv"
        )
      )

      p_G <- ggplot(
        integration_cor,
        aes(
          x = ModuleDisplay,
          y = GutDisplay,
          fill = Rho
        )
      ) +

        geom_tile(
          colour = "white",
          linewidth = 0.45
        ) +

        geom_text(
          aes(label = Sig),
          size = 4
        ) +

        facet_grid(
          FeatureType ~ .,
          scales = "free_y",
          space = "free_y",
          switch = "y"
        ) +

        scale_fill_gradient2(
          low = "#3B73B9",
          mid = "white",
          high = "#C94C4C",
          midpoint = 0,
          limits = c(-1, 1),
          name = "Spearman\nrho"
        ) +

        labs(
          x = NULL,
          y = NULL,
          caption = "Group-adjusted Spearman correlations. · nominal P < 0.05; * FDR < 0.10; ** FDR < 0.05. Exploratory only."
        ) +

        theme_classic(
          base_size = 14
        ) +

        theme(
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text.y.left = element_text(
            angle = 90,
            face = "bold",
            size = 10
          ),
          axis.text.x = element_text(
            angle = 35,
            hjust = 1,
            size = 9.3,
            colour = "black",
            lineheight = 0.95
          ),
          axis.text.y = element_text(
            size = 9.3,
            colour = "black",
            lineheight = 0.95
          ),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.caption = element_text(
            size = 8.8,
            hjust = 0
          ),
          plot.margin = margin(
            8,
            8,
            8,
            12
          )
        )

      save_plot(
        p_G,
        "FigS4I_GutMetabolite_OvaryModule_Correlation",
        SUPP_DIR,
        width = 11.2,
        height = max(
          6.0,
          3.0 +
            0.36 *
            nrow(gut_all)
        )
      )
    }
  }
}

# ======================================================================
# 22A. RESTORED MANUSCRIPT FIGURE 4G AND 4H (reviewer Major 6)
# ======================================================================
# The submitted composite Figure 4 contains:
#   G: six paired fecal / serum metabolite scatterplots
#   H: four genus CLR abundance panels
# Those are separate from the revised exploratory integration heatmap,
# which is now saved as FigS4I in Supplementary.
#
# NOTE: "serum" is the matrix label in the raw dataset. We use the
# matching terminology in the figure; change to "plasma" ONLY after
# verifying the biological specimen type against lab records.
# ======================================================================

MANUSCRIPT_FIG4_GH_DIR <- file.path(MAIN_DIR, "Restored_GH")
dir.create(MANUSCRIPT_FIG4_GH_DIR, recursive=TRUE, showWarnings=FALSE)

# ---------- H. Manuscript four-taxon CLR panel ----------
master_clr_path <- file.path(FIG3_TABLE_DIR, "genus_CLR_master.csv")
if (!file.exists(master_clr_path)) {
  stop("Manuscript Fig4H requires the canonical Fig3 CLR matrix: ",
       master_clr_path, "\nRun the corrected Fig3 script first.")
}
master_clr_df_gh <- readr::read_csv(master_clr_path,
                                   show_col_types=FALSE,
                                   name_repair="minimal")
if (!"Genus" %in% names(master_clr_df_gh) ||
    anyDuplicated(master_clr_df_gh$Genus) ||
    anyDuplicated(names(master_clr_df_gh))) {
  stop("Master genus CLR matrix has missing/duplicate identifiers")
}
genus_clr_gh <- as.matrix(master_clr_df_gh[,setdiff(names(master_clr_df_gh), "Genus"),drop=FALSE])
rownames(genus_clr_gh) <- master_clr_df_gh$Genus
storage.mode(genus_clr_gh) <- "double"
if (any(!is.finite(genus_clr_gh))) stop("Non-finite CLR values in master")
if (any(abs(colMeans(genus_clr_gh)) > 1e-7)) stop("Master CLR centering failed")

# Exact four taxa identified directly from the original submitted Figure 4H.
fig4h_genera <- c("Alistipes", "Prevotellaceae_UCG-001",
                  "Mucispirillum", "[Eubacterium]_coprostanoligenes_group_unclassified")
missing_4h <- setdiff(fig4h_genera, rownames(genus_clr_gh))
if (length(missing_4h)) {
  stop("Original Fig4H taxonomy bins absent from canonical CLR reference: ",
       paste(missing_4h,collapse=", "),
       "\nRerun Fig3_FINAL_v4_CLR_Fig4H_retention_FIXED.R before this script. ",
       "If Fig3 reports missing RAW taxonomy, inspect its QC taxonomy report. ",
       "Do not substitute taxa or compute a separate CLR reference.")
}
fig4h_m <- genus_clr_gh[fig4h_genera,,drop=FALSE]
fig4h_df <- as.data.frame(t(fig4h_m),check.names=FALSE) %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(Group=factor(group_from_sample(Sample),levels=GROUP_LEVELS)) %>%
  tidyr::pivot_longer(cols=dplyr::all_of(fig4h_genera),
                      names_to="Genus",values_to="CLR") %>%
  dplyr::mutate(Genus=factor(Genus,levels=fig4h_genera))
# Taxonomic transparency: this feature is a family-level unclassified bin,
# not a genus-level identification. Only plotting labels are modified.
fig4h_df$GenusDisplay <- as.character(fig4h_df$Genus)
fig4h_df$GenusDisplay[fig4h_df$GenusDisplay ==
  "[Eubacterium]_coprostanoligenes_group_unclassified"] <-
  "[Eubacterium]\ncoprostanoligenes group\n(unclassified)"
fig4h_df$GenusDisplay <- factor(fig4h_df$GenusDisplay, levels = c(
  "Alistipes", "Prevotellaceae_UCG-001", "Mucispirillum",
  "[Eubacterium]\ncoprostanoligenes group\n(unclassified)"))
if (any(is.na(fig4h_df$Group))) stop("Unexpected sample groups in master CLR")
readr::write_csv(fig4h_df,
  file.path(TABLE_DIR,"Fig4H_OriginalFourGenera_SampleLevel_CLR.csv"))

p_H_restored <- ggplot2::ggplot(fig4h_df,
  ggplot2::aes(x=Group,y=CLR,fill=Group)) +
  ggplot2::geom_boxplot(width=0.57,outlier.shape=NA,
                        alpha=0.26,color="grey35",linewidth=0.55) +
  ggplot2::geom_point(ggplot2::aes(color=Group),
    position=ggplot2::position_jitter(width=0.075,height=0,seed=20261117),
    size=2.3,alpha=0.96) +
  ggplot2::facet_wrap(~GenusDisplay,ncol=2,scales="free_y",
    labeller=ggplot2::labeller(GenusDisplay=ggplot2::label_value)) +
  ggplot2::scale_fill_manual(values=GROUP_COLS,guide="none") +
  ggplot2::scale_color_manual(values=GROUP_COLS,guide="none") +
  ggplot2::labs(x=NULL,y="Genus abundance (CLR)") +
  ggplot2::theme_classic(base_size=14) +
  ggplot2::theme(
    strip.background=ggplot2::element_blank(),
    strip.text=ggplot2::element_text(face="bold",size=12.5,lineheight=1.05,margin=ggplot2::margin(b=7)),
    axis.text.x=ggplot2::element_text(size=11.5,angle=35,hjust=1,color="black"),
    axis.text.y=ggplot2::element_text(size=11.5,color="black"),
    axis.title.y=ggplot2::element_text(size=13.5,face="bold"),
    panel.spacing=grid::unit(0.8,"lines"),
    plot.margin=ggplot2::margin(8,8,8,8)
  )
save_plot(p_H_restored,"Fig4H_OriginalFourGenera_CLR",
          dir=MANUSCRIPT_FIG4_GH_DIR,width=7.6,height=5.8,dpi=500)

# ---------- G. Six paired fecal-versus-serum metabolite panels ----------
# Recreate the ORIGINAL submitted Figure4G targets, never automatically
# substitute a different compound with a similar name.
fig4g_target <- data.frame(
  Label=c("5-Hydroxyindole-3-Acetic Acid", "Deoxycholic Acid",
          "Indole-3-Lactic Acid", "Kynurenic Acid",
          "Nord(e)oxycholic Acid", "Propionic Acid"),
  Pattern=c("^5hydroxyindole3aceticacid$|^5hiaa$",
            "^deoxycholicacid$|^dca$",
            "^indole3lacticacid$|^indole3lactate$|^ila$",
            "^kynurenicacid$|^kyna$",
            "^nordeoxycholicacid$|^nord(e)?oxycholicacid$|^nordca$|^nor(deoxy)?cholicacid$",
            "^propionicacid$|^propionate$"),
  stringsAsFactors=FALSE)

normalize_fig4g_name <- function(x) {
  tolower(gsub("[^a-z0-9]","",trimws(as.character(x))))
}

read_fig4g_metabolome <- function(path, source_label) {
  if (!file.exists(path)) stop("Missing ",source_label," metabolome: ",path)
  z <- readr::read_tsv(path,show_col_types=FALSE,
                       progress=FALSE,name_repair="minimal")
  if (!"Metabolite" %in% names(z)) {
    stop(source_label," metabolite file must contain column 'Metabolite'")
  }
  cols <- names(z)[grepl("^(MX|DT|GZFL|LH)_[0-9]+$",names(z))]
  if (!length(cols)) stop("No matched biological sample columns in ",path)
  m <- as.matrix(z[,cols,drop=FALSE]); storage.mode(m) <- "double"
  # Follow the Fig5 v4 intent: transform raw intensities if required,
  # then standardize each LC-MS FEATURE to z-scores across all 4 groups.
  finite_v <- as.numeric(m[is.finite(m)])
  if (!length(finite_v)) stop("No finite ",source_label," metabolite values")
  if (stats::quantile(finite_v,0.99,na.rm=TRUE) > 50 ||
      stats::median(finite_v,na.rm=TRUE) > 20) {
    shift <- if (min(finite_v,na.rm=TRUE) <= 0)
      abs(min(finite_v,na.rm=TRUE))+1 else 1
    m <- log2(m+shift)
  }
  nm <- normalize_fig4g_name(z$Metabolite)
  # Auditable target matching: never substitute an unrelated compound.
  match_qc <- dplyr::bind_rows(lapply(seq_len(nrow(fig4g_target)), function(i) {
    matched <- which(grepl(fig4g_target$Pattern[i], nm, perl=TRUE))
    data.frame(Omic=source_label, Target=fig4g_target$Label[i],
               N_MatchedFeatures=length(matched),
               MatchedNames=if(length(matched)) paste(unique(z$Metabolite[matched]),collapse=" | ") else "",
               stringsAsFactors=FALSE)
  }))
  readr::write_csv(match_qc,file.path(QC_DIR,paste0("Fig4G_",source_label,"_TargetMatching_QC.csv")))
  # Export all observed compound names to troubleshoot unmatched targets.
  readr::write_csv(data.frame(FeatureRow=seq_len(nrow(z)),
               OriginalName=as.character(z$Metabolite),
               NormalizedName=nm,stringsAsFactors=FALSE),
               file.path(QC_DIR,paste0("Fig4G_",source_label,"_AvailableNames.csv")))
  out_rows <- list()
  detail_rows <- list()
  for (i in seq_len(nrow(fig4g_target))) {
    idx <- which(grepl(fig4g_target$Pattern[i],nm,perl=TRUE))
    if (!length(idx)) next
    sub_m <- m[idx,,drop=FALSE]
    # Feature-wise median imputation only when a feature has observed values.
    for (j in seq_len(nrow(sub_m))) {
      observed <- sub_m[j,is.finite(sub_m[j,])]
      if (!length(observed)) next
      sub_m[j,!is.finite(sub_m[j,])] <- stats::median(observed)
    }
    sd0 <- apply(sub_m,1,stats::sd,na.rm=TRUE)
    valid <- which(is.finite(sd0) & sd0 > 0)
    if (!length(valid)) next
    sub_m <- sub_m[valid,,drop=FALSE]
    zmat <- t(scale(t(sub_m)))
    merged <- colMeans(zmat,na.rm=TRUE)
    out_rows[[length(out_rows)+1]] <- data.frame(
      Sample=colnames(zmat),
      Metabolite=fig4g_target$Label[i],
      Value=as.numeric(merged),
      Omic=source_label,
      stringsAsFactors=FALSE)
    detail_rows[[length(detail_rows)+1]] <- data.frame(
      Metabolite=fig4g_target$Label[i],Omic=source_label,
      FeatureRow=idx[valid],MatchedName=as.character(z$Metabolite[idx[valid]]),
      stringsAsFactors=FALSE)
  }
  empty_values <- data.frame(Sample=character(),Metabolite=character(),Value=double(),Omic=character())
  empty_mapping <- data.frame(Metabolite=character(),Omic=character(),FeatureRow=integer(),MatchedName=character())
  list(values=if(length(out_rows)) dplyr::bind_rows(out_rows) else empty_values,
       mapping=if(length(detail_rows)) dplyr::bind_rows(detail_rows) else empty_mapping)
}

fig4g_fecal <- read_fig4g_metabolome(FECAL_FILE,"Fecal")
fig4g_serum_file <- file.path(RAW_DIR,"serum_all_metab_data.tsv")
fig4g_serum <- read_fig4g_metabolome(fig4g_serum_file,"Serum")
readr::write_csv(dplyr::bind_rows(fig4g_fecal$mapping,fig4g_serum$mapping),
  file.path(QC_DIR,"Fig4G_Metabolite_FeatureMapping.csv"))

# Additional direct diagnostic: presence in BOTH modalities, regardless of sample pairing.
fig4g_coverage <- data.frame(
  Metabolite=fig4g_target$Label,
  FecalPresent=fig4g_target$Label %in% fig4g_fecal$values$Metabolite,
  SerumPresent=fig4g_target$Label %in% fig4g_serum$values$Metabolite,
  stringsAsFactors=FALSE)
readr::write_csv(fig4g_coverage,file.path(QC_DIR,"Fig4G_TargetCoverage.csv"))

fig4g_pair <- dplyr::inner_join(
  dplyr::select(fig4g_fecal$values,Sample,Metabolite,Fecal=Value),
  dplyr::select(fig4g_serum$values,Sample,Metabolite,Serum=Value),
  by=c("Sample","Metabolite")) %>%
  dplyr::mutate(Group=factor(group_from_sample(Sample),levels=GROUP_LEVELS),
                Metabolite=factor(Metabolite,levels=fig4g_target$Label))
if (!nrow(fig4g_pair)) {
  warning("Fig4G has no matched fecal-serum pairs under the verified target-name rules. ",
          "Fig4H was saved successfully. Inspect QC/Fig4G_*_TargetMatching_QC.csv ",
          "and QC/Fig4G_*_AvailableNames.csv. Do not substitute metabolites automatically.")
  writeLines("No matched fecal-serum target features. See Fig4G_Fecal/Serum_TargetMatching_QC.csv and available names.",
             file.path(QC_DIR,"Fig4G_NoPairedTargets.txt"))
} else {

fig4g_stats <- fig4g_pair %>%
  dplyr::group_by(Metabolite) %>%
  dplyr::group_modify(function(d,k) {
    if (nrow(d)<6 || stats::sd(d$Fecal)==0 || stats::sd(d$Serum)==0)
      return(data.frame(N=nrow(d),Rho=NA_real_,P=NA_real_,AdjustedRho=NA_real_))
    ct <- suppressWarnings(stats::cor.test(d$Fecal,d$Serum,
                                           method="spearman",exact=FALSE))
    # Group-adjusted rho is descriptive; naive cor.test P on estimated
    # residuals would not correctly account for fitted group effects.
    ra <- stats::residuals(stats::lm(Fecal~Group,data=d))
    rb <- stats::residuals(stats::lm(Serum~Group,data=d))
    ar <- suppressWarnings(stats::cor(ra,rb,method="spearman"))
    data.frame(N=nrow(d),Rho=unname(ct$estimate),P=ct$p.value,
               AdjustedRho=unname(ar))
  }) %>% dplyr::ungroup() %>%
  dplyr::mutate(FDR=stats::p.adjust(P,method="BH"),
                StatLabel=paste0("rho = ",sprintf("%.2f",Rho),
                                 "  P = ",format.pval(P,digits=2,eps=0.001)))
readr::write_csv(fig4g_pair,
  file.path(TABLE_DIR,"Fig4G_FecalSerum_SampleLevel_Zscores.csv"))
readr::write_csv(fig4g_stats,
  file.path(TABLE_DIR,"Fig4G_FecalSerum_Spearman_Statistics.csv"))

missing_4g <- setdiff(fig4g_target$Label,
                      as.character(unique(fig4g_pair$Metabolite)))
writeLines(if(length(missing_4g)) missing_4g else "All 6 original metabolites matched",
  file.path(QC_DIR,"Fig4G_MissingTargets.txt"))
if (length(missing_4g)) {
  warning("Only ",6-length(missing_4g),"/6 manuscript metabolites found in BOTH datasets: ",
          paste(missing_4g,collapse=", "),
          ". Check Fig4G_Metabolite_FeatureMapping.csv; no substitutions made.")
}

# Fig4G: definitive six-metabolite raw paired scatter plot.
# Fig4G final output is produced by 07_Fig4G_cross_compartment.R.
p_G_restored <- ggplot2::ggplot(fig4g_pair,
  ggplot2::aes(x=Fecal,y=Serum,color=Group)) +
  ggplot2::geom_smooth(ggplot2::aes(group=1),method="lm",
                       se=FALSE,color="grey45",linewidth=0.6) +
  ggplot2::geom_point(size=2.5,alpha=0.93) +
  ggplot2::facet_wrap(~Metabolite,scales="free",ncol=2,
    labeller=ggplot2::labeller(Metabolite=ggplot2::label_wrap_gen(width=19))) +
  ggplot2::scale_color_manual(values=GROUP_COLS,drop=FALSE) +
  ggplot2::labs(x="Fecal metabolite (z-score)",
                y="Plasma metabolite (z-score)",color=NULL) +
  ggplot2::theme_classic(base_size=14) +
  ggplot2::theme(
    legend.position="top",legend.text=ggplot2::element_text(size=9),
    strip.background=ggplot2::element_blank(),
    strip.text=ggplot2::element_text(face="bold",size=8.5,lineheight=0.9,margin=ggplot2::margin(t=3,b=3)),
    axis.text=ggplot2::element_text(size=10.5,color="black"),
    axis.title=ggplot2::element_text(face="bold",size=12.5),
    plot.caption=ggplot2::element_blank(),
    panel.spacing=grid::unit(0.8,"lines"),
    plot.margin=ggplot2::margin(8,10,8,8)
  )
save_plot(p_G_restored,"Fig4G_OriginalSix_FecalSerum_Scatter",
          dir=FIGURE_DIR,width=5.6,height=9.0,dpi=400)
}

# ======================================================================
# 23. Summary
# ======================================================================

summary_lines <- c(

  "============================================================",
  "FIGURE 4 v2 | DEG -> GO -> INTERACTION -> INTEGRATION",
  "============================================================",
  "",

  "Groups:",
  "  PCOS    = MX",
  "  LCD     = DT",
  "  GFC     = GZFL",
  "  Combine = LH",
  "",

  paste0(
    "RNA samples = ",
    length(count_sample_cols)
  ),

  paste0(
    "Genes before filter = ",
    nrow(count_mat)
  ),

  paste0(
    "Genes after filterByExpr = ",
    nrow(norm_expr)
  ),
  "",

  "Primary DEG method:",
  "  DESeq2, matching the company report",
  "Primary DEG threshold:",
  "  BH-FDR < 0.05 and |log2FC| >= 1",
  "",

  "GO rule:",
  "  Primary FDR-DEG is used when >=8 genes are available.",
  "  Otherwise P<0.05 + |log2FC|>=0.58 is used ONLY as an",
  "  exploratory fallback and is explicitly labeled in output.",
  "",

  "Main panels:",
  "  Fig4A DEG counts",
  "  Fig4B Venn of LCD/GFC/Combine DEGs vs PCOS",
  "  Fig4C shared/treatment-specific GO vs PCOS",
  "  Fig4D GO of direct treatment differences",
  "  Fig4E GO of LCDxGFC interaction",
  "  Fig4F LCDxGFC interaction volcano",
  "  Fig4F six biologically prioritized interaction genes",
  "  Restored Fig4G six paired fecal/serum metabolites (Main/Restored_GH)",
  "  Restored Fig4H four original genus CLR panels (Main/Restored_GH)",
  "  Supplementary FigS4I gut/metabolite-to-ovary integration",
  "",

  paste0(
    "Strict Combine-specific DOWN genes = ",
    nrow(
      combine_specific_down
    )
  ),

  paste0(
    "Strict Combine-specific UP genes = ",
    nrow(
      combine_specific_up
    )
  ),
  "",

  paste0(
    "Representative interaction genes: ",
    paste(
      representative_genes$Gene_name,
      collapse = "; "
    )
  ),
  "",

  "IMPORTANT INTERPRETATION RULES:",
  "  'Combine-beyond-both' means a stronger response away from PCOS than either monotherapy.",
  "  It does NOT by itself mean a more beneficial or more normal state.",
  "  Circadian interaction means non-additivity, not superior efficacy.",
  "",
  "IMPORTANT LIMITATION:",
  "  No healthy Control RNA-seq group is available.",
  "  Therefore this dataset supports treatment comparison and interaction,",
  "  but does not prove restoration to a healthy ovarian transcriptome.",
  "  True normalization is assessed only when an external PCOS-vs-healthy signature is supplied.",
  "",

  "============================================================"
)

writeLines(
  summary_lines,
  file.path(
    OUT_DIR,
    "SUMMARY_Figure4_v2.txt"
  )
)

cat(
  "\n============================================================\n"
)
cat(
  "Figure 4 v2 analysis completed.\n"
)
cat(
  "Output directory:\n",
  OUT_DIR,
  "\n"
)
cat(
  "============================================================\n"
)

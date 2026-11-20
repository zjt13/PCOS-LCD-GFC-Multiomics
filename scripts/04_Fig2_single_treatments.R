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

FINAL_FIGURE_MAP <- c("Fig2A_Genus_PCoA"="Fig2A", "FigS2B_16S_Alpha_Shannon"="Fig2B", "FigS2C_Genus_Composition_Top10"="Fig2C", "Fig2B_FecalMetabolome_PCA"="Fig2D", "Fig2C_Genus_TreatmentEffect_Map"="Fig2F", "Fig2D_FecalMetabolite_TreatmentEffect_Map"="Fig2G", "Fig2E_Genus_Metabolite_AdjustedCorrelation"="Fig2H", "FigS2D_Genus_Volcano_LCD_vs_PCOS"="FigS1C", "FigS2E_Genus_Volcano_GFC_vs_PCOS"="FigS1D", "FigS2H_Fecal_Volcano_LCD_vs_PCOS"="FigS1E", "FigS2I_Fecal_Volcano_GFC_vs_PCOS"="FigS1F")
# Individual plot theme() settings are preserved during export.
# No global font-size overrides are applied in save_plot().

options(stringsAsFactors = FALSE)
set.seed(20260930)

# ======================================================================
# 0. Paths
# ======================================================================

raw_dir <- DATA_DIR

otu_file <- file.path(raw_dir, "microbiome_otu_table.tsv")
fecal_file <- file.path(raw_dir, "fecal_all_metab_data.tsv")

taxonomy_file <- TAXONOMY_INPUT

out_dir <- file.path(WORK_DIR,"Fig2")

main_dir  <- file.path(out_dir, "Main")
supp_dir  <- file.path(out_dir, "Supplementary")
table_dir <- file.path(out_dir, "Tables")
qc_dir    <- file.path(out_dir, "QC")

dir.create(main_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir,    recursive = TRUE, showWarnings = FALSE)

# ======================================================================
# 1. Packages
# ======================================================================

cran_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "stringr",
  "ggrepel", "vegan", "pheatmap", "tibble"
)

missing_cran <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("limma", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("limma", ask = FALSE, update = FALSE)
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggrepel)
library(vegan)
library(pheatmap)
library(tibble)
library(limma)

# ======================================================================
# 2. Global settings
# ======================================================================

analysis_groups <- c("PCOS", "LCD", "GFC")

group_cols <- c(
  "PCOS" = "#B18459",
  "LCD"  = "#8E44AD",
  "GFC"  = "#F28E2B"
)

group_shapes <- c(
  "PCOS" = 16,  # circle
  "LCD"  = 15,  # square
  "GFC"  = 17   # triangle
)

effect_cols <- c(
  "Shared response"   = "#4C78A8",
  "LCD-dominant"      = "#7A5195",
  "GFC-dominant"      = "#F28E2B",
  "Opposite response" = "#D1495B",
  "Other"             = "grey75"
)

# Candidate threshold is used only for visual categorization and
# candidate prioritization, not for declaring definitive mechanisms.
FDR_CANDIDATE <- 0.10

# 16S filtering before CLR / differential analysis
GENUS_PREVALENCE_MIN <- 0.20       # detected in >=20% biological samples
GENUS_REL_ABUND_MIN  <- 0.0001     # overall mean relative abundance >=0.01%

# Cross-omics heatmap dimensions kept small to prevent text overlap
N_GENUS_COR <- 8
N_METAB_COR <- 10

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
  prefix <- sub("_.*$", "", x)
  dplyr::recode(
    prefix,
    "MX"   = "PCOS",
    "DT"   = "LCD",
    "GZFL" = "GFC",
    "LH"   = "Combine",
    "QC"   = "QC",
    .default = NA_character_
  )
}


# ----------------------------------------------------------------------
# GLOBAL sample naming/order helpers
#
# Internal raw IDs:
#   MX_1   -> PCOS_1
#   DT_1   -> LCD_1
#   GZFL_1 -> GFC_1
#
# Biological group order in Figure 2:
#   PCOS -> LCD -> GFC
# ----------------------------------------------------------------------

sample_number_from_id <- function(x) {
  suppressWarnings(
    as.numeric(
      sub("^.*_", "", x)
    )
  )
}

sample_display_name <- function(x) {

  grp <- group_from_sample(x)
  num <- sample_number_from_id(x)

  out <- ifelse(
    grp %in% c("PCOS", "LCD", "GFC", "Combine", "QC") &
      !is.na(num),
    paste0(grp, "_", num),
    x
  )

  out
}

make_sample_order_df <- function(
  samples,
  group_levels = c("PCOS", "LCD", "GFC")
) {

  out <- data.frame(
    OriginalSample = samples,
    Group = group_from_sample(samples),
    SampleNumber = sample_number_from_id(samples),
    stringsAsFactors = FALSE
  ) %>%
    filter(
      Group %in% group_levels
    ) %>%
    mutate(
      Group = factor(
        Group,
        levels = group_levels
      )
    ) %>%
    arrange(
      Group,
      SampleNumber,
      OriginalSample
    ) %>%
    mutate(
      DisplaySample = sample_display_name(
        OriginalSample
      )
    )

  out
}

is_sample_column <- function(x) {
  grepl("^(MX|DT|GZFL|LH|QC)_[0-9]+$", x)
}

safe_numeric_matrix <- function(df, sample_cols, row_ids) {
  m <- as.matrix(df[, sample_cols, drop = FALSE])
  suppressWarnings(storage.mode(m) <- "numeric")
  rownames(m) <- make.unique(as.character(row_ids))
  m
}

# Clean genus labels.
# Examples:
#   g__Bacteroides                  -> Bacteroides
#   g__norank_f__Muribaculaceae    -> Muribaculaceae_unclassified
#   missing genus + f__Lachnospiraceae -> Lachnospiraceae_unclassified
clean_taxon <- function(genus_raw, family_raw, order_raw) {

  clean_one <- function(x) {
    x <- ifelse(is.na(x), "", x)
    x <- gsub("^g__", "", x)
    x <- gsub("^f__", "", x)
    x <- gsub("^o__", "", x)
    trimws(x)
  }

  genus <- clean_one(genus_raw)
  fam   <- clean_one(family_raw)
  ord   <- clean_one(order_raw)

  out <- genus

  # informative "norank_f__Family" / "norank_o__Order" labels
  idx <- grepl("^norank_f__", out)
  out[idx] <- paste0(sub("^norank_f__", "", out[idx]), "_unclassified")

  idx <- grepl("^norank_o__", out)
  out[idx] <- paste0(sub("^norank_o__", "", out[idx]), "_unclassified")

  # blank / generic unclassified genus -> family fallback
  generic <- out == "" |
    is.na(out) |
    tolower(out) %in% c("norank", "unclassified", "uncultured", "unknown")

  out[generic & fam != ""] <- paste0(fam[generic & fam != ""], "_unclassified")

  still_generic <- out == "" |
    is.na(out) |
    tolower(out) %in% c("norank", "unclassified", "uncultured", "unknown")

  out[still_generic & ord != ""] <- paste0(ord[still_generic & ord != ""], "_unclassified")
  out[still_generic & ord == ""] <- "Unclassified"

  out
}

extract_rank <- function(tax, prefix) {
  pattern <- paste0("(?:^|;)", prefix, "([^;]+)")
  m <- stringr::str_match(tax, pattern)
  m[, 2]
}

clr_transform <- function(count_matrix, pseudocount = 0.5) {
  x <- count_matrix + pseudocount
  lx <- log(x)
  sweep(lx, 2, colMeans(lx), FUN = "-")
}

make_limma_results <- function(mat_feature_by_sample, sample_group) {

  sample_group <- factor(sample_group, levels = analysis_groups)
  design <- model.matrix(~0 + sample_group)
  colnames(design) <- analysis_groups

  fit <- limma::lmFit(mat_feature_by_sample, design)

  cont <- limma::makeContrasts(
    LCD_vs_PCOS = LCD - PCOS,
    GFC_vs_PCOS = GFC - PCOS,
    LCD_vs_GFC  = LCD - GFC,
    levels = design
  )

  fit2 <- limma::contrasts.fit(fit, cont)
  fit2 <- limma::eBayes(fit2, robust = TRUE)

  out <- list()

  for (contrast_name in colnames(cont)) {
    tt <- limma::topTable(
      fit2,
      coef = contrast_name,
      number = Inf,
      adjust.method = "BH",
      sort.by = "none"
    )
    tt$Feature <- rownames(tt)
    out[[contrast_name]] <- tt
  }

  out
}

combine_effect_tables <- function(res_list) {

  a <- res_list$LCD_vs_PCOS %>%
    select(
      Feature,
      logFC_LCD_PCOS = logFC,
      P_LCD_PCOS = P.Value,
      FDR_LCD_PCOS = adj.P.Val
    )

  b <- res_list$GFC_vs_PCOS %>%
    select(
      Feature,
      logFC_GFC_PCOS = logFC,
      P_GFC_PCOS = P.Value,
      FDR_GFC_PCOS = adj.P.Val
    )

  c <- res_list$LCD_vs_GFC %>%
    select(
      Feature,
      logFC_LCD_GFC = logFC,
      P_LCD_GFC = P.Value,
      FDR_LCD_GFC = adj.P.Val
    )

  a %>%
    inner_join(b, by = "Feature") %>%
    inner_join(c, by = "Feature") %>%
    mutate(
      ResponseClass = case_when(
        FDR_LCD_GFC < FDR_CANDIDATE &
          sign(logFC_LCD_PCOS) != sign(logFC_GFC_PCOS) &
          abs(logFC_LCD_PCOS) > 1e-8 &
          abs(logFC_GFC_PCOS) > 1e-8 ~ "Opposite response",

        FDR_LCD_GFC < FDR_CANDIDATE &
          abs(logFC_LCD_PCOS) >= abs(logFC_GFC_PCOS) ~ "LCD-dominant",

        FDR_LCD_GFC < FDR_CANDIDATE &
          abs(logFC_GFC_PCOS) > abs(logFC_LCD_PCOS) ~ "GFC-dominant",

        FDR_LCD_PCOS < FDR_CANDIDATE &
          FDR_GFC_PCOS < FDR_CANDIDATE &
          FDR_LCD_GFC >= FDR_CANDIDATE &
          sign(logFC_LCD_PCOS) == sign(logFC_GFC_PCOS) ~ "Shared response",

        TRUE ~ "Other"
      ),
      ResponseClass = factor(
        ResponseClass,
        levels = c(
          "Shared response", "LCD-dominant", "GFC-dominant",
          "Opposite response", "Other"
        )
      )
    )
}

make_effect_map <- function(
  df,
  title_text = NULL,
  label_col = "Feature",
  label_n = 7,
  highlight_requires_annotation = FALSE,
  max_background_points = Inf,
  max_highlight_per_class = Inf
) {

  dat <- df %>%
    mutate(
      LabelText = .data[[label_col]],

      LabelText = if (identical(label_col, "Feature")) {
        ifelse(
          is.na(LabelText) | LabelText == "",
          Feature,
          LabelText
        )
      } else {
        LabelText
      },

      IsAnnotated = !is.na(LabelText) &
        LabelText != "",

      priority = pmin(
        FDR_LCD_GFC,
        FDR_LCD_PCOS,
        FDR_GFC_PCOS,
        na.rm = TRUE
      ),

      IsHighlight =
        ResponseClass != "Other" &
        (
          !highlight_requires_annotation |
          IsAnnotated
        )
    )

  # ------------------------------------------------------------
  # Plotting-only reduction of dense points.
  # All features remain in statistical analyses and CSV tables.
  # ------------------------------------------------------------

  background_df <- dat %>%
    filter(!IsHighlight)

  if (is.finite(max_background_points) &&
      nrow(background_df) > max_background_points) {

    set.seed(20260930)

    keep_idx <- sample(
      seq_len(nrow(background_df)),
      size = max_background_points,
      replace = FALSE
    )

    background_plot <- background_df[
      keep_idx,
      ,
      drop = FALSE
    ]

  } else {
    background_plot <- background_df
  }

  highlight_df <- dat %>%
    filter(IsHighlight) %>%
    arrange(
      ResponseClass,
      FDR_LCD_GFC,
      priority,
      desc(abs(logFC_LCD_GFC))
    )

  if (is.finite(max_highlight_per_class)) {

    highlight_plot <- highlight_df %>%
      group_by(ResponseClass) %>%
      slice_head(
        n = max_highlight_per_class
      ) %>%
      ungroup()

  } else {
    highlight_plot <- highlight_df
  }

  # Labels: only biologically interpretable highlighted candidates.
  candidate_labels <- highlight_df %>%
    filter(
      IsAnnotated
    ) %>%
    arrange(
      FDR_LCD_GFC,
      priority,
      desc(abs(logFC_LCD_GFC))
    ) %>%
    distinct(
      LabelText,
      .keep_all = TRUE
    )

  label_df <- candidate_labels %>%
    slice_head(n = label_n) %>%
    mutate(
      Label = stringr::str_trunc(
        LabelText,
        30
      )
    )

  lim <- max(
    abs(
      c(
        dat$logFC_LCD_PCOS,
        dat$logFC_GFC_PCOS
      )
    ),
    na.rm = TRUE
  )

  if (!is.finite(lim) || lim == 0) {
    lim <- 1
  }

  lim <- lim * 1.08

  ggplot() +

    geom_hline(
      yintercept = 0,
      linewidth = 0.45,
      colour = "grey72"
    ) +

    geom_vline(
      xintercept = 0,
      linewidth = 0.45,
      colour = "grey72"
    ) +

    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.6,
      colour = "grey55"
    ) +

    # Faint discovery background
    geom_point(
      data = background_plot,
      aes(
        x = logFC_LCD_PCOS,
        y = logFC_GFC_PCOS
      ),
      colour = "grey72",
      size = 0.85,
      alpha = 0.16
    ) +

    # Selected interpretable treatment-response features
    geom_point(
      data = highlight_plot,
      aes(
        x = logFC_LCD_PCOS,
        y = logFC_GFC_PCOS,
        colour = ResponseClass
      ),
      size = 1.85,
      alpha = 0.60
    ) +

    ggrepel::geom_text_repel(
      data = label_df,
      aes(
        x = logFC_LCD_PCOS,
        y = logFC_GFC_PCOS,
        label = Label,
        colour = ResponseClass
      ),
      size = 3.25,
      max.overlaps = Inf,
      box.padding = 0.65,
      point.padding = 0.40,
      min.segment.length = 0,
      segment.colour = "grey55",
      segment.linewidth = 0.38,
      force = 2.6,
      force_pull = 0.65,
      max.time = 3,
      seed = 20260930,
      show.legend = FALSE
    ) +

    scale_colour_manual(
      values = effect_cols,
      breaks = c(
        "Shared response",
        "LCD-dominant",
        "GFC-dominant",
        "Opposite response"
      ),
      drop = FALSE
    ) +

    coord_cartesian(
      xlim = c(-lim, lim),
      ylim = c(-lim, lim),
      clip = "off"
    ) +

    guides(
      colour = guide_legend(
        nrow = 2,
        byrow = TRUE,
        title = NULL,
        override.aes = list(
          size = 3.0,
          alpha = 1
        )
      )
    ) +

    labs(
      title = title_text,
      x = "LCD effect vs PCOS",
      y = "GFC effect vs PCOS",
      colour = NULL
    ) +

    theme_classic(
      base_size = 14
    ) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),

      # Centering prevents the first legend symbol in each row
      # from being clipped at the left edge.
      legend.position = "top",
      legend.justification = "center",
      legend.box.just = "center",
      legend.direction = "horizontal",
      legend.text = element_text(
        size = 10.5
      ),
      legend.key.width = grid::unit(
        0.55,
        "cm"
      ),
      legend.spacing.x = grid::unit(
        0.12,
        "cm"
      ),
      legend.spacing.y = grid::unit(
        0.04,
        "cm"
      ),
      legend.margin = margin(
        0, 8, 5, 8
      ),

      axis.text = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        face = "bold"
      ),

      plot.margin = margin(
        t = 8,
        r = 14,
        b = 6,
        l = 14
      )
    )
}

make_volcano <- function(
  tt,
  title_text,
  label_n = 8,
  label_col = "Feature",
  annotated_only = FALSE
) {

  dat <- tt %>%
    mutate(
      FDR = pmax(
        adj.P.Val,
        1e-300
      ),
      Pplot = pmax(
        P.Value,
        1e-300
      ),
      Sig = adj.P.Val < FDR_CANDIDATE,
      Label = .data[[label_col]],
      LabelOK = !is.na(Label) &
        Label != ""
    )

  if (identical(label_col, "Feature")) {
    dat$Label[
      !dat$LabelOK
    ] <- dat$Feature[
      !dat$LabelOK
    ]
    dat$LabelOK <- TRUE
  }

  lab <- dat %>%
    filter(
      !annotated_only |
      LabelOK
    ) %>%
    arrange(
      adj.P.Val,
      desc(abs(logFC))
    ) %>%
    distinct(
      Label,
      .keep_all = TRUE
    ) %>%
    slice_head(
      n = label_n
    ) %>%
    mutate(
      Label = stringr::str_trunc(
        Label,
        30
      )
    )

  ggplot(
    dat,
    aes(
      x = logFC,
      y = -log10(Pplot)
    )
  ) +

    geom_vline(
      xintercept = 0,
      colour = "grey70",
      linewidth = 0.4
    ) +

    # nonsignificant background
    geom_point(
      data = dat %>%
        filter(!Sig),
      colour = "grey72",
      size = 1.0,
      alpha = 0.20
    ) +

    # significant features
    geom_point(
      data = dat %>%
        filter(Sig),
      colour = "#C44E52",
      size = 1.45,
      alpha = 0.48
    ) +

    ggrepel::geom_text_repel(
      data = lab,
      aes(label = Label),
      size = 3.2,
      max.overlaps = Inf,
      box.padding = 0.55,
      point.padding = 0.30,
      min.segment.length = 0,
      segment.colour = "grey50",
      segment.linewidth = 0.35,
      force = 2,
      max.time = 3,
      seed = 20260930,
      show.legend = FALSE
    ) +

    labs(
      title = title_text,
      x = "Effect size",
      y = expression(
        -log[10](italic(P))
      )
    ) +

    theme_classic(
      base_size = 14
    ) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      axis.text = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        face = "bold"
      )
    )
}

# ======================================================================
# 4. Read 16S OTU table
# ======================================================================

if (!file.exists(otu_file)) stop("Cannot find: ", otu_file)
if (!file.exists(taxonomy_file)) stop("Cannot find taxonomy file: ", taxonomy_file)

otu_df <- readr::read_tsv(
  otu_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

if (ncol(otu_df) < 2) {
  stop("microbiome_otu_table.tsv has fewer than 2 columns.")
}

otu_id_col <- names(otu_df)[1]
names(otu_df)[1] <- "OTU_ID"

otu_sample_cols <- names(otu_df)[vapply(names(otu_df), is_sample_column, logical(1))]

if (length(otu_sample_cols) == 0) {
  stop("No 16S sample columns matched MX_/DT_/GZFL_/LH_ naming.")
}

otu_group <- group_from_sample(otu_sample_cols)

bio_idx <- otu_group %in% analysis_groups
otu_bio_cols <- otu_sample_cols[bio_idx]
otu_bio_group <- factor(otu_group[bio_idx], levels = analysis_groups)

otu_mat <- safe_numeric_matrix(
  otu_df,
  otu_bio_cols,
  otu_df$OTU_ID
)

# Remove invalid / empty OTUs
otu_mat[is.na(otu_mat)] <- 0
otu_mat[otu_mat < 0] <- 0
otu_mat <- otu_mat[rowSums(otu_mat) > 0, , drop = FALSE]

# Save sequencing depth
seq_depth <- data.frame(
  Sample = colnames(otu_mat),
  Group = factor(group_from_sample(colnames(otu_mat)), levels = analysis_groups),
  Reads = colSums(otu_mat)
)

write_csv(seq_depth, file.path(table_dir, "16S_sequencing_depth.csv"))

# ======================================================================
# 5. Read taxonomy and map OTU -> Genus
# ======================================================================

tax_raw <- read.delim(
  taxonomy_file,
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  fill = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (ncol(tax_raw) < 2) {
  stop(
    "seqs_tax_assignments.txt needs at least 2 tab-separated columns: ",
    "OTU_ID and taxonomy string."
  )
}

tax_df <- data.frame(
  OTU_ID = trimws(as.character(tax_raw[[1]])),
  Taxonomy = trimws(as.character(tax_raw[[2]])),
  stringsAsFactors = FALSE
)

if (ncol(tax_raw) >= 3) {
  suppressWarnings(
    tax_df$TaxConfidence <- as.numeric(as.character(tax_raw[[3]]))
  )
}

tax_df <- tax_df %>%
  mutate(
    Phylum = extract_rank(Taxonomy, "p__"),
    Class  = extract_rank(Taxonomy, "c__"),
    Order  = extract_rank(Taxonomy, "o__"),
    Family = extract_rank(Taxonomy, "f__"),
    Genus_raw = extract_rank(Taxonomy, "g__"),
    Species = extract_rank(Taxonomy, "s__"),
    Genus = clean_taxon(Genus_raw, Family, Order)
  )

write_csv(
  tax_df,
  file.path(table_dir, "16S_taxonomy_parsed.csv")
)

otu_ids <- rownames(otu_mat)
map_idx <- match(otu_ids, tax_df$OTU_ID)

taxonomy_qc <- data.frame(
  Metric = c(
    "OTUs_in_count_table",
    "OTUs_mapped_to_taxonomy",
    "OTUs_unmapped",
    "Mapped_fraction"
  ),
  Value = c(
    length(otu_ids),
    sum(!is.na(map_idx)),
    sum(is.na(map_idx)),
    mean(!is.na(map_idx))
  )
)

write_csv(
  taxonomy_qc,
  file.path(qc_dir, "taxonomy_mapping_QC.csv")
)

message(
  "Taxonomy mapping: ",
  sum(!is.na(map_idx)), "/", length(otu_ids),
  " OTUs mapped (",
  round(100 * mean(!is.na(map_idx)), 1), "%)"
)

if (mean(!is.na(map_idx)) < 0.80) {
  warning(
    "Less than 80% of OTUs mapped to seqs_tax_assignments.txt. ",
    "Check whether OTU IDs use the same naming system."
  )
}

otu_tax <- data.frame(
  OTU_ID = otu_ids,
  Genus = tax_df$Genus[map_idx],
  stringsAsFactors = FALSE
)

otu_tax$Genus[is.na(otu_tax$Genus) | otu_tax$Genus == ""] <- "Unclassified"

otu_for_genus <- as.data.frame(otu_mat) %>%
  rownames_to_column("OTU_ID") %>%
  left_join(otu_tax, by = "OTU_ID")

otu_for_genus$Genus[is.na(otu_for_genus$Genus)] <- "Unclassified"

genus_df <- otu_for_genus %>%
  select(Genus, all_of(colnames(otu_mat))) %>%
  group_by(Genus) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

genus_mat_all <- as.matrix(genus_df[, -1, drop = FALSE])
rownames(genus_mat_all) <- genus_df$Genus
storage.mode(genus_mat_all) <- "numeric"

write_csv(
  genus_df,
  file.path(table_dir, "16S_genus_count_table_PCOS_LCD_GFC.csv")
)

# Relative abundance
genus_rel_all <- sweep(
  genus_mat_all,
  2,
  pmax(colSums(genus_mat_all), 1),
  FUN = "/"
)

# Differential-analysis filter
prevalence <- rowMeans(genus_mat_all > 0)
mean_rel <- rowMeans(genus_rel_all)

keep_genus <- prevalence >= GENUS_PREVALENCE_MIN &
  mean_rel >= GENUS_REL_ABUND_MIN

genus_mat <- genus_mat_all[keep_genus, , drop = FALSE]
genus_rel <- genus_rel_all[keep_genus, , drop = FALSE]

write_csv(
  data.frame(
    Genus = rownames(genus_mat_all),
    Prevalence = prevalence,
    MeanRelativeAbundance = mean_rel,
    KeepForDifferentialAnalysis = keep_genus
  ),
  file.path(table_dir, "16S_genus_filtering.csv")
)

# ======================================================================
# 6. Main Fig2A | Genus-level PCoA
# ======================================================================

# Bray-Curtis on genus-level relative abundance
bray <- vegan::vegdist(t(genus_rel_all), method = "bray")

pcoa <- cmdscale(
  bray,
  k = 2,
  eig = TRUE,
  add = TRUE
)

pcoa_df <- data.frame(
  Sample = rownames(pcoa$points),
  Axis1 = pcoa$points[, 1],
  Axis2 = pcoa$points[, 2],
  stringsAsFactors = FALSE
) %>%
  mutate(
    Group = factor(group_from_sample(Sample), levels = analysis_groups)
  )

eig_pos <- pcoa$eig[pcoa$eig > 0]
var1 <- 100 * pcoa$eig[1] / sum(eig_pos)
var2 <- 100 * pcoa$eig[2] / sum(eig_pos)

meta_16s <- data.frame(
  Sample = colnames(genus_rel_all),
  Group = factor(group_from_sample(colnames(genus_rel_all)), levels = analysis_groups)
)

perm <- vegan::adonis2(
  bray ~ Group,
  data = meta_16s,
  permutations = 999
)

capture.output(
  perm,
  file = file.path(table_dir, "16S_Genus_PERMANOVA.txt")
)

p_A <- ggplot(
  pcoa_df,
  aes(
    x = Axis1,
    y = Axis2,
    colour = Group,
    shape = Group
  )
) +
  stat_ellipse(
    aes(group = Group),
    type = "t",
    linewidth = 0.65,
    linetype = 2,
    alpha = 0.8,
    show.legend = FALSE
  ) +
  geom_point(
    size = 3.0,
    alpha = 0.88,
    stroke = 0.4
  ) +
  scale_colour_manual(
    values = group_cols,
    breaks = analysis_groups,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = group_shapes,
    breaks = analysis_groups,
    drop = FALSE
  ) +
  labs(
    x = sprintf("PCoA1 (%.1f%%)", var1),
    y = sprintf("PCoA2 (%.1f%%)", var2),
    colour = NULL,
    shape = NULL
  ) +
  guides(
    colour = guide_legend(
      order = 1,
      override.aes = list(
        size = 3.4,
        alpha = 1
      )
    ),
    shape = guide_legend(
      order = 1,
      override.aes = list(
        size = 3.4,
        alpha = 1
      )
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.box.just = "center",
    legend.text = element_text(size = 11),
    legend.key.width = grid::unit(0.55, "cm"),
    legend.spacing.x = grid::unit(0.18, "cm"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(face = "bold")
  )

save_plot(
  p_A,
  "Fig2A_Genus_PCoA",
  main_dir,
  width = 5.4,
  height = 4.6
)

write_csv(pcoa_df, file.path(table_dir, "Fig2A_Genus_PCoA_coordinates.csv"))

# ======================================================================
# 7. 16S differential analysis | CLR + limma
# ======================================================================

genus_clr <- clr_transform(genus_mat, pseudocount = 0.5)

genus_res <- make_limma_results(
  genus_clr,
  otu_bio_group
)

for (nm in names(genus_res)) {
  write_csv(
    genus_res[[nm]],
    file.path(table_dir, paste0("Genus_", nm, "_limma.csv"))
  )
}

genus_effect <- combine_effect_tables(genus_res)

write_csv(
  genus_effect,
  file.path(table_dir, "Genus_treatment_effect_summary.csv")
)

# Main Fig2C
p_C <- make_effect_map(
  genus_effect,
  title_text = NULL,
  label_col = "Feature",
  label_n = 7,
  max_background_points = Inf,
  max_highlight_per_class = Inf
) +
  labs(
    x = "Genus effect: LCD vs PCOS (CLR scale)",
    y = "Genus effect: GFC vs PCOS (CLR scale)"
  )

save_plot(
  p_C,
  "Fig2C_Genus_TreatmentEffect_Map",
  main_dir,
  width = 7.0,
  height = 5.4
)

# ======================================================================
# 8. Read and preprocess fecal metabolomics
# ======================================================================

if (!file.exists(fecal_file)) stop("Cannot find: ", fecal_file)

fecal_df <- readr::read_tsv(
  fecal_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

fecal_sample_cols_all <- names(fecal_df)[
  vapply(names(fecal_df), is_sample_column, logical(1))
]

if (length(fecal_sample_cols_all) == 0) {
  stop("No fecal metabolomics sample columns detected.")
}

fecal_group_all <- group_from_sample(fecal_sample_cols_all)

fecal_bio_idx <- fecal_group_all %in% analysis_groups
fecal_bio_cols <- fecal_sample_cols_all[fecal_bio_idx]
fecal_bio_group <- factor(
  fecal_group_all[fecal_bio_idx],
  levels = analysis_groups
)

# ----------------------------------------------------------------------
# Metabolite internal ID + biological annotation
#
# IMPORTANT STRATEGY FOR UNTARGED METABOLOMICS:
#   - ALL QC-passed features are retained for PCA and differential analysis.
#   - ONLY features with a valid name in the "Metabolite" column are used
#     for biological labeling / interpretation (Fig2D labels, Fig2E, S2K).
#   - metab_id remains the unique internal feature identifier.
# ----------------------------------------------------------------------

id_candidates <- c(
  "metab_id", "metab.id", "Metab_ID",
  "ID", "id"
)

id_col <- id_candidates[id_candidates %in% names(fecal_df)][1]

if (is.na(id_col) || length(id_col) == 0) {
  stop(
    "No metab_id-like column was found in fecal_all_metab_data.tsv."
  )
}

if (!"Metabolite" %in% names(fecal_df)) {
  stop(
    "Column 'Metabolite' was not found. ",
    "The current Figure 2 workflow explicitly requires this column ",
    "for biological metabolite names."
  )
}

feature_id <- as.character(fecal_df[[id_col]])

missing_feature_id <- is.na(feature_id) | trimws(feature_id) == ""
feature_id[missing_feature_id] <- paste0(
  "Feature_",
  which(missing_feature_id)
)

feature_id <- make.unique(feature_id)

# ------------------------------------------------------------
# Strict validation of the Metabolite column.
# Only genuine readable metabolite names are accepted.
# ------------------------------------------------------------

is_valid_metabolite_annotation <- function(x) {

  x <- trimws(as.character(x))
  xl <- tolower(x)

  valid <- !is.na(x) &
    x != "" &
    x != "-" &
    !xl %in% c(
      "na", "n/a", "none",
      "unknown", "unidentified", "unannotated",
      "no match", "not identified",
      "not available", "null"
    ) &

    # Internal feature IDs
    !grepl(
      "^metab[_ .-]*[0-9]+$",
      x,
      ignore.case = TRUE
    ) &

    # RT / m/z peak labels
    !grepl(
      "m/z|\\bmz\\b",
      x,
      ignore.case = TRUE
    ) &

    !grepl(
      "^[0-9]+\\.?[0-9]*[_:-][0-9]+\\.?[0-9]*(m/z|mz)?$",
      x,
      ignore.case = TRUE
    ) &

    # Pure numeric strings
    !grepl(
      "^[0-9.]+$",
      x
    ) &

    # Generic peak / feature names
    !grepl(
      "^(peak|feature|unknown|unidentified)[_ .-]*[0-9]+$",
      x,
      ignore.case = TRUE
    )

  valid
}

metabolite_raw <- trimws(as.character(fecal_df$Metabolite))

annotated_flag <- is_valid_metabolite_annotation(
  metabolite_raw
)

display_name <- rep(
  NA_character_,
  length(metabolite_raw)
)

display_name[annotated_flag] <- metabolite_raw[
  annotated_flag
]

metab_annotation <- data.frame(
  Feature = feature_id,
  DisplayName = display_name,
  Annotated = annotated_flag,
  stringsAsFactors = FALSE
)

# Keep source annotation columns for later checking.
annotation_keep <- unique(c(
  id_col,
  "Metabolite",
  intersect(
    c(
      "Compound",
      "KEGG Compound ID",
      "Library ID",
      "CAS ID",
      "drugbank_id",
      "InChIKey",
      "SMILES",
      "chemical_direct_parent",
      "chemical_class",
      "chemical_superclass",
      "chemical_subclass"
    ),
    names(fecal_df)
  )
))

annotation_extra <- fecal_df[
  ,
  annotation_keep,
  drop = FALSE
]

annotation_extra$Feature <- feature_id

metab_annotation_full <- metab_annotation %>%
  left_join(
    annotation_extra,
    by = "Feature"
  )

write_csv(
  metab_annotation_full,
  file.path(
    table_dir,
    "FecalMetabolite_annotation_map.csv"
  )
)

annotation_summary <- data.frame(
  Metric = c(
    "Total metabolic features",
    "Features with valid Metabolite annotation",
    "Unannotated features",
    "Annotation rate (%)",
    "Unique annotated metabolite names"
  ),
  Value = c(
    length(feature_id),
    sum(annotated_flag),
    sum(!annotated_flag),
    round(100 * mean(annotated_flag), 2),
    length(unique(display_name[annotated_flag]))
  )
)

write_csv(
  annotation_summary,
  file.path(
    qc_dir,
    "fecal_metabolite_annotation_summary.csv"
  )
)

message(
  "Valid Metabolite annotations: ",
  sum(annotated_flag),
  "/", length(annotated_flag),
  " features (",
  round(100 * mean(annotated_flag), 1),
  "%); unique names = ",
  length(unique(display_name[annotated_flag]))
)

fecal_mat_raw <- safe_numeric_matrix(
  fecal_df,
  fecal_bio_cols,
  feature_id
)

# Remove features with >20% missing biological values
missing_rate <- rowMeans(is.na(fecal_mat_raw))
keep_metab <- missing_rate <= 0.20
fecal_mat_raw <- fecal_mat_raw[keep_metab, , drop = FALSE]

# Median-based feature-wise imputation
impute_feature <- function(x) {
  if (all(is.na(x))) return(rep(0, length(x)))
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

fecal_mat_imp <- t(apply(fecal_mat_raw, 1, impute_feature))
colnames(fecal_mat_imp) <- colnames(fecal_mat_raw)
rownames(fecal_mat_imp) <- rownames(fecal_mat_raw)

# Detect whether values already look log-transformed.
finite_values <- as.numeric(fecal_mat_imp[is.finite(fecal_mat_imp)])
q50 <- median(finite_values, na.rm = TRUE)
q99 <- quantile(finite_values, 0.99, na.rm = TRUE)

already_log_like <- (
  q99 < 50 &&
  q50 < 20
)

if (already_log_like) {
  fecal_mat <- fecal_mat_imp
  metab_transform_note <- "No additional log transform; input appears log-like."
} else {
  if (min(fecal_mat_imp, na.rm = TRUE) < 0) {
    shift <- abs(min(fecal_mat_imp, na.rm = TRUE)) + 1
  } else {
    shift <- 1
  }
  fecal_mat <- log2(fecal_mat_imp + shift)
  metab_transform_note <- paste0("Applied log2(x + ", signif(shift, 4), ").")
}

# Remove zero-variance features
metab_sd <- apply(fecal_mat, 1, sd, na.rm = TRUE)
fecal_mat <- fecal_mat[is.finite(metab_sd) & metab_sd > 0, , drop = FALSE]

writeLines(
  c(
    paste0("Biological features retained: ", nrow(fecal_mat)),
    paste0("Biological samples: ", ncol(fecal_mat)),
    metab_transform_note
  ),
  file.path(qc_dir, "fecal_metabolome_preprocessing.txt")
)

# ======================================================================
# 9. Main Fig2B | Fecal metabolome PCA
# ======================================================================

pca_fecal <- prcomp(
  t(fecal_mat),
  center = TRUE,
  scale. = TRUE
)

var_exp <- 100 * pca_fecal$sdev^2 / sum(pca_fecal$sdev^2)

pca_fecal_df <- data.frame(
  Sample = rownames(pca_fecal$x),
  PC1 = pca_fecal$x[, 1],
  PC2 = pca_fecal$x[, 2],
  stringsAsFactors = FALSE
) %>%
  mutate(
    Group = factor(group_from_sample(Sample), levels = analysis_groups)
  )

p_B <- ggplot(
  pca_fecal_df,
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
    linewidth = 0.65,
    linetype = 2,
    alpha = 0.8,
    show.legend = FALSE
  ) +
  geom_point(
    size = 3.0,
    alpha = 0.88,
    stroke = 0.4
  ) +
  scale_colour_manual(
    values = group_cols,
    breaks = analysis_groups,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = group_shapes,
    breaks = analysis_groups,
    drop = FALSE
  ) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", var_exp[2]),
    colour = NULL,
    shape = NULL
  ) +
  guides(
    colour = guide_legend(
      order = 1,
      override.aes = list(
        size = 3.4,
        alpha = 1
      )
    ),
    shape = guide_legend(
      order = 1,
      override.aes = list(
        size = 3.4,
        alpha = 1
      )
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.box.just = "center",
    legend.text = element_text(size = 11),
    legend.key.width = grid::unit(0.55, "cm"),
    legend.spacing.x = grid::unit(0.18, "cm"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(face = "bold")
  )

save_plot(
  p_B,
  "Fig2B_FecalMetabolome_PCA",
  main_dir,
  width = 5.4,
  height = 4.6
)

write_csv(
  pca_fecal_df,
  file.path(table_dir, "Fig2B_Fecal_PCA_coordinates.csv")
)


# Diagnostic: exact duplicated PCA coordinates (if any).
pca_dup <- pca_fecal_df %>%
  count(PC1, PC2, name = "N") %>%
  filter(N > 1)

write_csv(
  pca_dup,
  file.path(
    qc_dir,
    "Fig2B_exact_duplicate_PCA_coordinates.csv"
  )
)

# ======================================================================
# 10. Fecal metabolome differential analysis | limma
# ======================================================================

fecal_res <- make_limma_results(
  fecal_mat,
  fecal_bio_group
)

for (nm in names(fecal_res)) {
  write_csv(
    fecal_res[[nm]],
    file.path(table_dir, paste0("FecalMetabolite_", nm, "_limma.csv"))
  )
}

# Attach readable metabolite names to statistical results.
for (nm in names(fecal_res)) {
  fecal_res[[nm]] <- fecal_res[[nm]] %>%
    left_join(metab_annotation, by = "Feature")
}

# Re-save annotated differential tables.
for (nm in names(fecal_res)) {
  write_csv(
    fecal_res[[nm]],
    file.path(
      table_dir,
      paste0("FecalMetabolite_", nm, "_limma_annotated.csv")
    )
  )
}

fecal_effect <- combine_effect_tables(fecal_res) %>%
  left_join(
    metab_annotation,
    by = "Feature"
  ) %>%
  mutate(
    Annotated = ifelse(
      is.na(Annotated),
      FALSE,
      Annotated
    )
  )

write_csv(
  fecal_effect,
  file.path(
    table_dir,
    "FecalMetabolite_treatment_effect_summary.csv"
  )
)


# Biological interpretation table:
# annotated metabolites only, one representative feature per name.
fecal_effect_annotated_unique <- fecal_effect %>%
  filter(
    Annotated,
    !is.na(DisplayName),
    DisplayName != ""
  ) %>%
  mutate(
    Priority = pmin(
      FDR_LCD_PCOS,
      FDR_GFC_PCOS,
      FDR_LCD_GFC,
      na.rm = TRUE
    )
  ) %>%
  arrange(
    Priority,
    FDR_LCD_GFC,
    desc(abs(logFC_LCD_GFC))
  ) %>%
  distinct(
    DisplayName,
    .keep_all = TRUE
  )

write_csv(
  fecal_effect_annotated_unique,
  file.path(
    table_dir,
    "FecalMetabolite_ANNOTATED_unique_treatment_effects.csv"
  )
)

# Main Fig2D
p_D <- make_effect_map(
  fecal_effect,
  title_text = NULL,
  label_col = "DisplayName",
  label_n = 5,
  highlight_requires_annotation = TRUE,

  # Plotting only:
  # thousands of background peaks otherwise collapse into a solid band.
  # Statistical tests still use ALL QC-passed features.
  max_background_points = 2500,

  # Keep the most statistically prioritized annotated candidates
  # per response class for a readable main figure.
  max_highlight_per_class = 120
) +
  labs(
    x = "Metabolite effect: LCD vs PCOS",
    y = "Metabolite effect: GFC vs PCOS"
  )

save_plot(
  p_D,
  "Fig2D_FecalMetabolite_TreatmentEffect_Map",
  main_dir,
  width = 7.0,
  height = 5.4
)

# ======================================================================
# 11. Main Fig2E | Group-adjusted Genus-metabolite correlations
# ======================================================================

# Exact sample matching is required.
matched_samples <- intersect(
  colnames(genus_clr),
  colnames(fecal_mat)
)

matched_group <- factor(
  group_from_sample(matched_samples),
  levels = analysis_groups
)

if (length(matched_samples) < 9) {
  warning(
    "Only ", length(matched_samples),
    " exact-matched 16S/fecal samples. Correlation panel may be unstable."
  )
}

# Select treatment-responsive features using the smallest FDR among
# LCD-vs-PCOS, GFC-vs-PCOS and direct LCD-vs-GFC.
genus_select <- genus_effect %>%
  mutate(
    Priority = pmin(
      FDR_LCD_PCOS,
      FDR_GFC_PCOS,
      FDR_LCD_GFC,
      na.rm = TRUE
    )
  ) %>%
  arrange(Priority, desc(abs(logFC_LCD_GFC))) %>%
  slice_head(n = N_GENUS_COR) %>%
  pull(Feature)

metab_select_df <- fecal_effect %>%
  filter(
    Annotated,
    !is.na(DisplayName),
    DisplayName != ""
  ) %>%
  mutate(
    Priority = pmin(
      FDR_LCD_PCOS,
      FDR_GFC_PCOS,
      FDR_LCD_GFC,
      na.rm = TRUE
    )
  ) %>%
  arrange(
    Priority,
    FDR_LCD_GFC,
    desc(abs(logFC_LCD_GFC))
  ) %>%

  # One representative LC-MS feature per annotated metabolite name.
  # This avoids counting several adduct/isotope/peak features as
  # independent biological metabolites.
  distinct(
    DisplayName,
    .keep_all = TRUE
  ) %>%
  slice_head(
    n = N_METAB_COR
  )

if (nrow(metab_select_df) == 0) {
  stop(
    "No annotated metabolites were available for Fig2E after filtering."
  )
}

metab_select <- metab_select_df$Feature

metab_display_lookup <- setNames(
  metab_select_df$DisplayName,
  metab_select_df$Feature
)

g_cor_mat <- genus_clr[genus_select, matched_samples, drop = FALSE]
m_cor_mat <- fecal_mat[metab_select, matched_samples, drop = FALSE]

# Residualize each feature by Group, then correlate residuals.
# This reduces correlations driven only by between-group mean shifts.
residualize_by_group <- function(mat, grp) {
  out <- matrix(
    NA_real_,
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
  for (i in seq_len(nrow(mat))) {
    fit <- lm(as.numeric(mat[i, ]) ~ grp)
    out[i, ] <- residuals(fit)
  }
  out
}

g_resid <- residualize_by_group(g_cor_mat, matched_group)
m_resid <- residualize_by_group(m_cor_mat, matched_group)

cor_rows <- list()
k <- 1

for (g in rownames(g_resid)) {
  for (m in rownames(m_resid)) {

    ct <- suppressWarnings(
      cor.test(
        as.numeric(g_resid[g, ]),
        as.numeric(m_resid[m, ]),
        method = "spearman",
        exact = FALSE
      )
    )

    cor_rows[[k]] <- data.frame(
      Genus = g,
      Metabolite = m,
      Rho = unname(ct$estimate),
      P = ct$p.value,
      stringsAsFactors = FALSE
    )
    k <- k + 1
  }
}

cor_df <- bind_rows(cor_rows) %>%
  mutate(
    FDR = p.adjust(P, method = "BH"),
    Sig = case_when(
      FDR < 0.05 ~ "**",
      FDR < 0.10 ~ "*",
      TRUE ~ ""
    ),
    Metabolite_display = metab_display_lookup[
      Metabolite
    ],
    Metabolite_display = stringr::str_trunc(
      Metabolite_display,
      30
    ),
    Genus_display = stringr::str_trunc(Genus, 28)
  )

write_csv(
  cor_df,
  file.path(table_dir, "Genus_FecalMetabolite_groupAdjusted_Spearman.csv")
)

p_E <- ggplot(
  cor_df,
  aes(
    x = Metabolite_display,
    y = Genus_display,
    fill = Rho
  )
) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = Sig), size = 4.1, colour = "black") +
  scale_fill_gradient2(
    low = "#3B73B9",
    mid = "white",
    high = "#C94C4C",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho"
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      colour = "black",
      size = 9
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 9
    ),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
  )

save_plot(
  p_E,
  "Fig2E_Genus_Metabolite_AdjustedCorrelation",
  main_dir,
  width = 7.4,
  height = 5.4
)

# ======================================================================
# 12. Supplementary | 16S sequencing depth
# ======================================================================

seq_depth$Group <- factor(
  seq_depth$Group,
  levels = c("PCOS", "LCD", "GFC")
)

p_S2A <- ggplot(
  seq_depth,
  aes(x = Group, y = Reads, colour = Group)
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    colour = "grey30"
  ) +
  geom_jitter(width = 0.08, size = 2.8) +
  scale_colour_manual(values = group_cols, guide = "none") +
  labs(x = NULL, y = "Sequencing reads") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

save_plot(
  p_S2A,
  "FigS2A_16S_SequencingDepth",
  supp_dir,
  width = 4.8,
  height = 4.2
)

# ======================================================================
# 13. Supplementary | OTU-level Shannon diversity
# ======================================================================

otu_sample_by_feature <- t(otu_mat)

shannon <- vegan::diversity(
  otu_sample_by_feature,
  index = "shannon"
)

alpha_df <- data.frame(
  Sample = names(shannon),
  Shannon = as.numeric(shannon),
  Group = factor(
    group_from_sample(names(shannon)),
    levels = analysis_groups
  )
)

kw_alpha <- kruskal.test(Shannon ~ Group, data = alpha_df)

writeLines(
  capture.output(kw_alpha),
  file.path(table_dir, "16S_Shannon_KruskalWallis.txt")
)

alpha_df$Group <- factor(
  alpha_df$Group,
  levels = c("PCOS", "LCD", "GFC")
)

p_S2B <- ggplot(
  alpha_df,
  aes(x = Group, y = Shannon, colour = Group)
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    colour = "grey30"
  ) +
  geom_jitter(width = 0.08, size = 2.8) +
  scale_colour_manual(values = group_cols, guide = "none") +
  labs(
    x = NULL,
    y = "Shannon diversity",
    subtitle = paste0(
      "Kruskal-Wallis P = ",
      signif(kw_alpha$p.value, 3)
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

save_plot(
  p_S2B,
  "FigS2B_16S_Alpha_Shannon",
  supp_dir,
  width = 4.8,
  height = 4.3
)

# ======================================================================
# 14. Supplementary | Top-10 Genus composition
# Fixed group order: PCOS -> LCD -> GFC
# Standardized sample labels: PCOS_n / LCD_n / GFC_n
# ======================================================================

top10_genus <- names(
  sort(
    rowMeans(genus_rel_all),
    decreasing = TRUE
  )
)[
  seq_len(
    min(
      10,
      nrow(genus_rel_all)
    )
  )
]

comp_mat <- genus_rel_all

comp_df <- as.data.frame(
  t(comp_mat)
) %>%
  rownames_to_column(
    "OriginalSample"
  ) %>%
  mutate(
    Group = factor(
      group_from_sample(
        OriginalSample
      ),
      levels = c(
        "PCOS",
        "LCD",
        "GFC"
      )
    ),
    DisplaySample = sample_display_name(
      OriginalSample
    ),
    SampleNumber = sample_number_from_id(
      OriginalSample
    )
  ) %>%
  filter(
    !is.na(Group)
  ) %>%
  pivot_longer(
    cols = -c(
      OriginalSample,
      DisplaySample,
      SampleNumber,
      Group
    ),
    names_to = "Genus",
    values_to = "RelativeAbundance"
  ) %>%
  mutate(
    Genus2 = ifelse(
      Genus %in% top10_genus,
      Genus,
      "Other"
    )
  ) %>%
  group_by(
    OriginalSample,
    DisplaySample,
    SampleNumber,
    Group,
    Genus2
  ) %>%
  summarise(
    RelativeAbundance = sum(
      RelativeAbundance
    ),
    .groups = "drop"
  )

# Fixed biological sample order:
# PCOS_1... -> LCD_1... -> GFC_1...
comp_sample_order <- comp_df %>%
  distinct(
    OriginalSample,
    DisplaySample,
    SampleNumber,
    Group
  ) %>%
  arrange(
    Group,
    SampleNumber,
    OriginalSample
  )

comp_df$DisplaySample <- factor(
  comp_df$DisplaySample,
  levels = comp_sample_order$DisplaySample
)

write_csv(
  comp_sample_order,
  file.path(
    table_dir,
    "FigS2C_sample_order_and_labels.csv"
  )
)

p_S2C <- ggplot(
  comp_df,
  aes(
    x = DisplaySample,
    y = RelativeAbundance * 100,
    fill = Genus2
  )
) +
  geom_col(
    width = 0.85
  ) +
  facet_grid(
    ~ Group,
    scales = "free_x",
    space = "free_x"
  ) +
  labs(
    x = NULL,
    y = "Relative abundance (%)",
    fill = "Genus"
  ) +
  theme_classic(
    base_size = 14
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9.5,
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    ),
    legend.text = element_text(
      size = 9
    ),
    legend.position = "right",
    panel.spacing.x = grid::unit(
      0.65,
      "lines"
    )
  )

# Fig2C-specific compact legend and wider data panel.
# Set styles HERE rather than in save_plot(), so subsequent edits to this
# panel's theme() will be respected in exported PDF and PNG files.
p_S2C <- p_S2C +
  scale_fill_discrete(
    name = "Genus",
    labels = function(x) stringr::str_wrap(x, width = 23)
  ) +
  guides(fill = guide_legend(
    title.position = "top",
    keyheight = grid::unit(0.25, "cm"),
    keywidth = grid::unit(0.25, "cm")
  )) +
  theme(
    legend.text = element_text(size = 8, lineheight = 0.95),
    legend.title = element_text(size = 9),
    legend.key.height = grid::unit(0.25, "cm"),
    legend.key.width = grid::unit(0.25, "cm"),
    legend.spacing.y = grid::unit(0.02, "cm"),
    legend.box.spacing = grid::unit(0.08, "cm"),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    plot.margin = margin(5, 4, 4, 5)
  )

save_plot(
  p_S2C,
  "FigS2C_Genus_Composition_Top10",
  supp_dir,
  width = 8.6,
  height = 4.9
)

# ======================================================================
# 15. Supplementary | Genus volcano plots
# ======================================================================

p_S2D <- make_volcano(
  genus_res$LCD_vs_PCOS,
  "Genus: LCD vs PCOS"
)
p_S2E <- make_volcano(
  genus_res$GFC_vs_PCOS,
  "Genus: GFC vs PCOS"
)
p_S2F <- make_volcano(
  genus_res$LCD_vs_GFC,
  "Genus: LCD vs GFC"
)

save_plot(p_S2D, "FigS2D_Genus_Volcano_LCD_vs_PCOS", supp_dir, 5.4, 4.6)
save_plot(p_S2E, "FigS2E_Genus_Volcano_GFC_vs_PCOS", supp_dir, 5.4, 4.6)
save_plot(p_S2F, "FigS2F_Genus_Volcano_LCD_vs_GFC",  supp_dir, 5.4, 4.6)

# ======================================================================
# 16. Supplementary | Fecal QC PCA
# ======================================================================

qc_cols <- fecal_sample_cols_all[
  group_from_sample(fecal_sample_cols_all) == "QC"
]

bio_all_cols <- fecal_sample_cols_all[
  group_from_sample(fecal_sample_cols_all) %in%
    c("PCOS", "LCD", "GFC", "Combine")
]

qc_pca_cols <- c(bio_all_cols, qc_cols)

if (length(qc_cols) >= 2 && length(qc_pca_cols) >= 5) {

  qc_mat_raw <- safe_numeric_matrix(
    fecal_df,
    qc_pca_cols,
    make.unique(as.character(fecal_df[[id_col]]))
  )

  qc_missing <- rowMeans(is.na(qc_mat_raw))
  qc_mat_raw <- qc_mat_raw[qc_missing <= 0.20, , drop = FALSE]

  qc_mat_imp <- t(apply(qc_mat_raw, 1, impute_feature))
  colnames(qc_mat_imp) <- colnames(qc_mat_raw)
  rownames(qc_mat_imp) <- rownames(qc_mat_raw)

  if (already_log_like) {
    qc_mat <- qc_mat_imp
  } else {
    if (min(qc_mat_imp, na.rm = TRUE) < 0) {
      shift2 <- abs(min(qc_mat_imp, na.rm = TRUE)) + 1
    } else {
      shift2 <- 1
    }
    qc_mat <- log2(qc_mat_imp + shift2)
  }

  qc_sd <- apply(qc_mat, 1, sd, na.rm = TRUE)
  qc_mat <- qc_mat[is.finite(qc_sd) & qc_sd > 0, , drop = FALSE]

  qc_pca <- prcomp(t(qc_mat), center = TRUE, scale. = TRUE)
  qc_var <- 100 * qc_pca$sdev^2 / sum(qc_pca$sdev^2)

  qc_df <- data.frame(
    Sample = rownames(qc_pca$x),
    PC1 = qc_pca$x[, 1],
    PC2 = qc_pca$x[, 2],
    Group = group_from_sample(rownames(qc_pca$x)),
    stringsAsFactors = FALSE
  )

  qc_cols_plot <- c(
    "PCOS" = "#B18459",
    "LCD" = "#8E44AD",
    "GFC" = "#F28E2B",
    "Combine" = "#E76F51",
    "QC" = "grey30"
  )

  qc_df$Group <- factor(
    qc_df$Group,
    levels = c(
      "PCOS",
      "LCD",
      "GFC",
      "Combine",
      "QC"
    )
  )

  qc_shapes_plot <- c(
    "PCOS" = 21,
    "LCD" = 22,
    "GFC" = 24,
    "Combine" = 23,
    "QC" = 25
  )

  p_S2G <- ggplot(
    qc_df,
    aes(
      x = PC1,
      y = PC2,
      fill = Group,
      shape = Group
    )
  ) +
    geom_point(
      size = 2.7,
      alpha = 0.75,
      colour = "grey20",
      stroke = 0.45
    ) +
    scale_fill_manual(values = qc_cols_plot) +
    scale_shape_manual(values = qc_shapes_plot) +
    labs(
      x = sprintf("PC1 (%.1f%%)", qc_var[1]),
      y = sprintf("PC2 (%.1f%%)", qc_var[2]),
      fill = NULL,
      shape = NULL
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(
          size = 3.1,
          alpha = 1,
          colour = "grey20"
        )
      ),
      shape = "none"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.justification = "left",
      axis.text = element_text(colour = "black"),
      axis.title = element_text(face = "bold")
    )

  save_plot(
    p_S2G,
    "FigS2G_Fecal_QC_PCA",
    supp_dir,
    width = 5.4,
    height = 4.6
  )

  write_csv(
    qc_df,
    file.path(table_dir, "FigS2G_Fecal_QC_PCA_coordinates.csv")
  )
}

# ======================================================================
# 17. Supplementary | Fecal metabolite volcano plots
# ======================================================================

p_S2H <- make_volcano(
  fecal_res$LCD_vs_PCOS,
  "Fecal metabolites: LCD vs PCOS",
  label_n = 8,
  label_col = "DisplayName",
  annotated_only = TRUE
)

p_S2I <- make_volcano(
  fecal_res$GFC_vs_PCOS,
  "Fecal metabolites: GFC vs PCOS",
  label_n = 8,
  label_col = "DisplayName",
  annotated_only = TRUE
)

p_S2J <- make_volcano(
  fecal_res$LCD_vs_GFC,
  "Fecal metabolites: LCD vs GFC",
  label_n = 8,
  label_col = "DisplayName",
  annotated_only = TRUE
)

save_plot(p_S2H, "FigS2H_Fecal_Volcano_LCD_vs_PCOS", supp_dir, 5.4, 4.6)
save_plot(p_S2I, "FigS2I_Fecal_Volcano_GFC_vs_PCOS", supp_dir, 5.4, 4.6)
save_plot(p_S2J, "FigS2J_Fecal_Volcano_LCD_vs_GFC",  supp_dir, 5.4, 4.6)

# ======================================================================
# 18. Supplementary | Annotated metabolite heatmaps
#
# S2K = BALANCED heatmap:
#   representative annotated metabolites from
#   LCD vs PCOS, GFC vs PCOS, and LCD vs GFC.
#
# S2L = direct LCD-vs-GFC annotated-metabolite heatmap.
#
# Sample order is always:
#   PCOS -> LCD -> GFC
# Sample labels are always:
#   PCOS_n / LCD_n / GFC_n
# ======================================================================

select_top_annotated <- function(
  tt,
  contrast_name,
  n_each = 8
) {

  tt %>%
    filter(
      !is.na(DisplayName),
      DisplayName != ""
    ) %>%
    arrange(
      adj.P.Val,
      desc(abs(logFC))
    ) %>%
    distinct(
      DisplayName,
      .keep_all = TRUE
    ) %>%
    slice_head(
      n = n_each
    ) %>%
    mutate(
      SelectedBy = contrast_name
    )
}

top_LCD_PCOS <- select_top_annotated(
  fecal_res$LCD_vs_PCOS,
  "LCD vs PCOS",
  n_each = 8
)

top_GFC_PCOS <- select_top_annotated(
  fecal_res$GFC_vs_PCOS,
  "GFC vs PCOS",
  n_each = 8
)

top_LCD_GFC <- select_top_annotated(
  fecal_res$LCD_vs_GFC,
  "LCD vs GFC",
  n_each = 8
)

balanced_candidates <- bind_rows(
  top_LCD_PCOS,
  top_GFC_PCOS,
  top_LCD_GFC
) %>%
  group_by(
    DisplayName
  ) %>%
  arrange(
    adj.P.Val,
    desc(abs(logFC)),
    .by_group = TRUE
  ) %>%
  slice_head(
    n = 1
  ) %>%
  ungroup()

if (nrow(balanced_candidates) == 0) {
  stop(
    "No annotated metabolites were available for balanced FigS2K."
  )
}

write_csv(
  balanced_candidates,
  file.path(
    table_dir,
    "FigS2K_balanced_metabolites_selected.csv"
  )
)

balanced_heat_mat <- fecal_mat[
  balanced_candidates$Feature,
  ,
  drop = FALSE
]

rownames(
  balanced_heat_mat
) <- make.unique(
  stringr::str_trunc(
    balanced_candidates$DisplayName,
    42
  )
)

balanced_heat_z <- t(
  scale(
    t(
      balanced_heat_mat
    )
  )
)

balanced_heat_z[
  !is.finite(
    balanced_heat_z
  )
] <- 0

# Global standardized sample order
heat_sample_order <- make_sample_order_df(
  colnames(
    balanced_heat_z
  ),
  group_levels = c(
    "PCOS",
    "LCD",
    "GFC"
  )
)

ordered_original_samples <- heat_sample_order$OriginalSample
display_samples <- heat_sample_order$DisplaySample

balanced_heat_z <- balanced_heat_z[
  ,
  ordered_original_samples,
  drop = FALSE
]

colnames(
  balanced_heat_z
) <- display_samples

anno_col <- data.frame(
  Group = factor(
    as.character(
      heat_sample_order$Group
    ),
    levels = c(
      "PCOS",
      "LCD",
      "GFC"
    )
  ),
  row.names = display_samples
)

anno_cols <- list(
  Group = group_cols[
    c(
      "PCOS",
      "LCD",
      "GFC"
    )
  ]
)

group_counts <- table(
  factor(
    as.character(
      heat_sample_order$Group
    ),
    levels = c(
      "PCOS",
      "LCD",
      "GFC"
    )
  )
)

gaps_col <- cumsum(
  as.integer(
    group_counts
  )
)

if (length(gaps_col) > 1) {
  gaps_col <- gaps_col[
    -length(gaps_col)
  ]
} else {
  gaps_col <- NULL
}

write_csv(
  heat_sample_order,
  file.path(
    table_dir,
    "FigS2K_sample_order_and_labels.csv"
  )
)

# Wider output to avoid crowded sample labels
heatmap_width <- max(
  9.4,
  4.2 +
    0.48 *
    ncol(
      balanced_heat_z
    )
)

heatmap_height <- max(
  6.6,
  2.8 +
    0.19 *
    nrow(
      balanced_heat_z
    )
)

pdf(
  file.path(FIGURE_DIR,"Fig2E.pdf"),
  width = 0.70 * heatmap_width,
  height = 0.70 * heatmap_height
)

pheatmap(
  balanced_heat_z,
  annotation_col = anno_col,
  annotation_colors = anno_cols,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 12,
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = 45,
  border_color = NA,
  main = ""
)

dev.off()

png(
  file.path(FIGURE_DIR,"Fig2E.png"),
  width = 0.70 * heatmap_width,
  height = 0.70 * heatmap_height,
  units = "in",
  res = 400
)

pheatmap(
  balanced_heat_z,
  annotation_col = anno_col,
  annotation_colors = anno_cols,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 12,
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = 45,
  border_color = NA,
  main = ""
)

dev.off()

# ----------------------------------------------------------------------
# S2L: direct LCD-vs-GFC contrast
# ----------------------------------------------------------------------

direct_heat_selection <- fecal_res$LCD_vs_GFC %>%
  filter(
    !is.na(DisplayName),
    DisplayName != ""
  ) %>%
  arrange(
    adj.P.Val,
    desc(abs(logFC))
  ) %>%
  distinct(
    DisplayName,
    .keep_all = TRUE
  ) %>%
  slice_head(
    n = 20
  )

if (nrow(direct_heat_selection) > 0) {

  direct_heat_mat <- fecal_mat[
    direct_heat_selection$Feature,
    ,
    drop = FALSE
  ]

  rownames(
    direct_heat_mat
  ) <- make.unique(
    stringr::str_trunc(
      direct_heat_selection$DisplayName,
      42
    )
  )

  direct_heat_z <- t(
    scale(
      t(
        direct_heat_mat
      )
    )
  )

  direct_heat_z[
    !is.finite(
      direct_heat_z
    )
  ] <- 0

  direct_heat_z <- direct_heat_z[
    ,
    ordered_original_samples,
    drop = FALSE
  ]

  colnames(
    direct_heat_z
  ) <- display_samples

  write_csv(
    direct_heat_selection,
    file.path(
      table_dir,
      "FigS2L_LCD_vs_GFC_metabolites_selected.csv"
    )
  )

  direct_width <- max(
    9.4,
    4.2 +
      0.48 *
      ncol(
        direct_heat_z
      )
  )

  direct_height <- max(
    6.2,
    2.8 +
      0.19 *
      nrow(
        direct_heat_z
      )
  )

  pdf(
    file.path(tempdir(),"FigS2L_LCD_vs_GFC_AnnotatedMetabolite_Heatmap.pdf"),
    width = direct_width,
    height = direct_height
  )

  pheatmap(
    direct_heat_z,
    annotation_col = anno_col,
    annotation_colors = anno_cols,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    gaps_col = gaps_col,
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize = 13,
    fontsize_row = 13,
    fontsize_col = 13,
    angle_col = 45,
    border_color = NA,
    main = ""
  )

  dev.off()

  png(
    file.path(tempdir(),"FigS2L_LCD_vs_GFC_AnnotatedMetabolite_Heatmap.png"),
    width = direct_width,
    height = direct_height,
    units = "in",
    res = 400
  )

  pheatmap(
    direct_heat_z,
    annotation_col = anno_col,
    annotation_colors = anno_cols,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    gaps_col = gaps_col,
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize = 13,
    fontsize_row = 13,
    fontsize_col = 13,
    angle_col = 45,
    border_color = NA,
    main = ""
  )

  dev.off()
}

# ======================================================================
# 19. Summary report
# ======================================================================

count_class <- function(x) {
  as.data.frame(table(x$ResponseClass)) %>%
    arrange(desc(Freq))
}

genus_class_count <- count_class(genus_effect)
metab_class_count <- count_class(fecal_effect)


# ======================================================================
# Global display-name / group-order consistency QC
# ======================================================================

display_qc <- data.frame(
  RawPrefix = c("MX", "DT", "GZFL"),
  BiologicalGroup = c("PCOS", "LCD", "GFC"),
  ExampleRawSample = c("MX_1", "DT_1", "GZFL_1"),
  RequiredDisplaySample = c("PCOS_1", "LCD_1", "GFC_1"),
  GroupOrder = c(1, 2, 3),
  stringsAsFactors = FALSE
)

write_csv(
  display_qc,
  file.path(
    qc_dir,
    "GLOBAL_sample_naming_and_group_order.csv"
  )
)

summary_lines <- c(
  "============================================================",
  "FIGURE 2 ANALYSIS SUMMARY",
  "============================================================",
  "",
  "Biological groups used in Figure 2:",
  "  PCOS = MX",
  "  LCD  = DT",
  "  GFC  = GZFL",
  "",
  "GLOBAL plotting convention:",
  "  Group order = PCOS -> LCD -> GFC",
  "  Sample labels = PCOS_n / LCD_n / GFC_n",
  "  Raw MX_n / DT_n / GZFL_n IDs are never used as display labels.",
  "",
  "Not used in Figure 2 biological comparisons:",
  "  Combine = LH (reserved for later combination-treatment figure)",
  "  QC = metabolomics quality-control samples only",
  "",
  paste0("16S OTUs retained before genus aggregation: ", nrow(otu_mat)),
  paste0(
    "Taxonomy mapping rate: ",
    round(100 * mean(!is.na(map_idx)), 1), "%"
  ),
  paste0("Genus taxa before differential filter: ", nrow(genus_mat_all)),
  paste0("Genus taxa used for CLR/limma: ", nrow(genus_mat)),
  "",
  "16S genus-level PERMANOVA:",
  paste(capture.output(perm), collapse = "\n"),
  "",
  paste0("Fecal metabolite features used: ", nrow(fecal_mat)),
  paste0(
    "Valid names in Metabolite column: ",
    sum(metab_annotation$Annotated),
    "/", nrow(metab_annotation),
    " features"
  ),
  paste0(
    "Unique annotated metabolite names: ",
    length(
      unique(
        metab_annotation$DisplayName[
          metab_annotation$Annotated
        ]
      )
    )
  ),
  "Metabolomics strategy:",
  "  - PCA / differential analysis: all QC-passed features.",
  "  - Fig2D background: all QC-passed features.",
  "  - Fig2D highlighted/labeled candidates: annotated metabolites only.",
  "  - Fig2E / FigS2K: annotated metabolites only, deduplicated by name.",
  paste0("Fecal preprocessing: ", metab_transform_note),
  "",
  paste0(
    "Exact matched 16S-fecal samples for Fig2E: ",
    length(matched_samples)
  ),
  "",
  "Genus treatment-response classes:",
  paste(capture.output(print(genus_class_count, row.names = FALSE)), collapse = "\n"),
  "",
  "Fecal-metabolite treatment-response classes:",
  paste(capture.output(print(metab_class_count, row.names = FALSE)), collapse = "\n"),
  "",
  "Interpretation rule:",
  "  - A feature is not called 'different between LCD and GFC' simply",
  "    because it is significant in one treatment-vs-PCOS contrast and",
  "    not significant in the other.",
  "  - Direct LCD-vs-GFC evidence is explicitly included.",
  "  - Fig2E correlations are adjusted for group mean differences by",
  "    residualizing each feature against Group before Spearman analysis.",
  "",
  "Main panels:",
  "  A: Genus-level Bray-Curtis PCoA",
  "  B: Fecal metabolome PCA",
  "  C: Genus treatment-effect map",
  "  D: Fecal-metabolite treatment-effect map",
  "  E: Group-adjusted Genus-metabolite correlation heatmap",
  "",
  "Supplementary metabolite heatmaps:",
  "  S2K = balanced annotated metabolites from LCD-vs-PCOS,",
  "        GFC-vs-PCOS, and LCD-vs-GFC.",
  "  S2L = direct LCD-vs-GFC annotated metabolite heatmap.",
  "",
  "All panels are saved separately as PDF and PNG.",
  "============================================================"
)

writeLines(
  summary_lines,
  file.path(out_dir, "SUMMARY_Figure2.txt")
)

cat("\n============================================================\n")
cat("Figure 2 analysis completed.\n")
cat("Output directory:\n")
cat(out_dir, "\n")
cat("\nMain panels:\n")
cat(main_dir, "\n")
cat("\nSupplementary panels:\n")
cat(supp_dir, "\n")
cat("\nTables:\n")
cat(table_dir, "\n")
cat("============================================================\n")

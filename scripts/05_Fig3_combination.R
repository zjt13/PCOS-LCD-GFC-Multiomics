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

FINAL_FIGURE_MAP <- c("Fig3A_Genus_PCoA_4Groups"="FigS1A", "Fig3B_FecalMetabolome_PCA_4Groups"="FigS1B", "Fig3C1_Genus_CombinationIncrement_Map"="Fig3A", "Fig3C2_FecalMetabolite_CombinationIncrement_Map"="Fig3B", "Fig3D_KeyGenus_Abundance"="Fig3C", "Fig3E1_KeyMetabolite_Abundance"="Fig3D", "Fig3E2_LiteratureInformed_MetabolicModule_Scores"="FigS1G", "Fig3F_Candidate_GenusMetabolite_Correlation"="FigS1H")
# Individual plot theme() settings are preserved during export.
# No global font-size overrides are applied in save_plot().

options(stringsAsFactors = FALSE)
set.seed(20260930)

# ======================================================================
# 0. Paths
# ======================================================================

raw_dir <- DATA_DIR

otu_file <- file.path(
  raw_dir,
  "microbiome_otu_table.tsv"
)

fecal_file <- file.path(
  raw_dir,
  "fecal_all_metab_data.tsv"
)

taxonomy_file <- TAXONOMY_INPUT

out_dir <- file.path(WORK_DIR,"Fig3")

main_dir <- file.path(out_dir, "Main")
supp_dir <- file.path(out_dir, "Supplementary")
table_dir <- file.path(out_dir, "Tables")
qc_dir <- file.path(out_dir, "QC")

dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

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
  "vegan",
  "pheatmap",
  "tibble",
  "statmod"
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

if (!requireNamespace("limma", quietly = TRUE)) {

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages(
      "BiocManager",
      repos = "https://cloud.r-project.org"
    )
  }

  BiocManager::install(
    "limma",
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
library(vegan)
library(pheatmap)
library(tibble)
library(limma)

# ======================================================================
# 2. Global settings
# ======================================================================

analysis_groups <- c(
  "PCOS",
  "LCD",
  "GFC",
  "Combine"
)

group_cols <- c(
  "PCOS" = "#B18459",
  "LCD" = "#8E44AD",
  "GFC" = "#F28E2B",
  "Combine" = "#E76F51"
)

group_shapes <- c(
  "PCOS" = 16,
  "LCD" = 15,
  "GFC" = 17,
  "Combine" = 18
)

increment_cols <- c(
  "Beyond both singles" = "#D1495B",
  "LCD-like" = "#7A5195",
  "GFC-like" = "#F28E2B",
  "Intermediate" = "#4C78A8",
  "No clear extra effect" = "grey75"
)

FDR_CANDIDATE <- 0.10

GENUS_PREVALENCE_MIN <- 0.20
GENUS_REL_ABUND_MIN <- 0.0001

N_KEY_GENUS <- 6
N_KEY_METAB <- 6

# Literature / prior-Fig2 priority list:
# These taxa are NOT forced into the main figure. They receive priority
# only if the user's own direct Combine-vs-single-treatment statistics
# support them.
literature_priority_genera <- c(
  "Bacteroides",
  "Alloprevotella",
  "Lachnospiraceae_UCG-008",
  "Ruminococcaceae_UCG-003",
  "Lachnospiraceae_NK4A136_group",
  "Muribaculum",
  "Muribaculaceae_unclassified",
  "Blautia",
  "Mucispirillum",
  "Prevotellaceae_UCG-001",
  "[Eubacterium]_ventriosum_group",
  "Butyribacter",
  "Clostridia_vadinBB60_group",
  "Enterorhabdus"
)

# ======================================================================
# 3. General helper functions
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

is_sample_column <- function(x) {
  grepl(
    "^(MX|DT|GZFL|LH|QC)_[0-9]+$",
    x
  )
}

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

  ifelse(
    !is.na(grp) & !is.na(num),
    paste0(grp, "_", num),
    x
  )
}

make_sample_order_df <- function(
  samples,
  group_levels = analysis_groups
) {

  data.frame(
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
      ),
      DisplaySample = sample_display_name(
        OriginalSample
      )
    ) %>%
    arrange(
      Group,
      SampleNumber,
      OriginalSample
    )
}

safe_numeric_matrix <- function(
  df,
  sample_cols,
  row_ids
) {

  m <- as.matrix(
    df[, sample_cols, drop = FALSE]
  )

  suppressWarnings(
    storage.mode(m) <- "numeric"
  )

  rownames(m) <- make.unique(
    as.character(row_ids)
  )

  m
}

extract_rank <- function(
  tax,
  prefix
) {

  pattern <- paste0(
    "(?:^|;)",
    prefix,
    "([^;]+)"
  )

  m <- stringr::str_match(
    tax,
    pattern
  )

  m[, 2]
}

clean_taxon <- function(
  genus_raw,
  family_raw,
  order_raw
) {

  clean_one <- function(x) {

    x <- ifelse(
      is.na(x),
      "",
      x
    )

    x <- gsub("^g__", "", x)
    x <- gsub("^f__", "", x)
    x <- gsub("^o__", "", x)

    trimws(x)
  }

  genus <- clean_one(genus_raw)
  fam <- clean_one(family_raw)
  ord <- clean_one(order_raw)

  out <- genus

  idx <- grepl(
    "^norank_f__",
    out
  )

  out[idx] <- paste0(
    sub(
      "^norank_f__",
      "",
      out[idx]
    ),
    "_unclassified"
  )

  idx <- grepl(
    "^norank_o__",
    out
  )

  out[idx] <- paste0(
    sub(
      "^norank_o__",
      "",
      out[idx]
    ),
    "_unclassified"
  )

  generic <- out == "" |
    is.na(out) |
    tolower(out) %in% c(
      "norank",
      "unclassified",
      "uncultured",
      "unknown"
    )

  out[
    generic & fam != ""
  ] <- paste0(
    fam[generic & fam != ""],
    "_unclassified"
  )

  still_generic <- out == "" |
    is.na(out) |
    tolower(out) %in% c(
      "norank",
      "unclassified",
      "uncultured",
      "unknown"
    )

  out[
    still_generic & ord != ""
  ] <- paste0(
    ord[still_generic & ord != ""],
    "_unclassified"
  )

  out[
    still_generic & ord == ""
  ] <- "Unclassified"

  out
}

clr_transform <- function(
  count_matrix,
  pseudocount = 0.5
) {

  x <- count_matrix + pseudocount
  lx <- log(x)

  sweep(
    lx,
    2,
    colMeans(lx),
    FUN = "-"
  )
}

is_valid_metabolite_annotation <- function(x) {

  x <- trimws(
    as.character(x)
  )

  xl <- tolower(x)

  !is.na(x) &
    x != "" &
    x != "-" &
    !xl %in% c(
      "na",
      "n/a",
      "none",
      "unknown",
      "unidentified",
      "unannotated",
      "no match",
      "not identified",
      "not available",
      "null"
    ) &
    !grepl(
      "^metab[_ .-]*[0-9]+$",
      x,
      ignore.case = TRUE
    ) &
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
    !grepl(
      "^[0-9.]+$",
      x
    )
}

# ======================================================================
# 4. Statistical helpers
# ======================================================================

make_limma_results_4groups <- function(
  mat_feature_by_sample,
  sample_group
) {

  sample_group <- factor(
    sample_group,
    levels = analysis_groups
  )

  design <- model.matrix(
    ~0 + sample_group
  )

  colnames(design) <- analysis_groups

  fit <- limma::lmFit(
    mat_feature_by_sample,
    design
  )

  cont <- limma::makeContrasts(

    Combine_vs_PCOS =
      Combine - PCOS,

    LCD_vs_PCOS =
      LCD - PCOS,

    GFC_vs_PCOS =
      GFC - PCOS,

    Combine_vs_LCD =
      Combine - LCD,

    Combine_vs_GFC =
      Combine - GFC,

    LCD_vs_GFC =
      LCD - GFC,

    levels = design
  )

  fit2 <- limma::contrasts.fit(
    fit,
    cont
  )

  fit2 <- limma::eBayes(
    fit2,
    robust = TRUE
  )

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

pairwise_permanova <- function(
  dist_obj,
  metadata,
  group_col = "Group",
  pairs = list(
    c("Combine", "PCOS"),
    c("Combine", "LCD"),
    c("Combine", "GFC")
  ),
  permutations = 999
) {

  dmat <- as.matrix(dist_obj)

  res <- lapply(
    pairs,
    function(pair) {

      keep <- metadata[[group_col]] %in% pair
      meta_sub <- metadata[keep, , drop = FALSE]

      sample_ids <- rownames(meta_sub)

      d_sub <- as.dist(
        dmat[
          sample_ids,
          sample_ids,
          drop = FALSE
        ]
      )

      meta_sub$PairGroup <- droplevels(
        factor(
          meta_sub[[group_col]],
          levels = pair
        )
      )

      ad <- vegan::adonis2(
        d_sub ~ PairGroup,
        data = meta_sub,
        permutations = permutations
      )

      data.frame(
        Group1 = pair[1],
        Group2 = pair[2],
        F = ad$F[1],
        R2 = ad$R2[1],
        P = ad$`Pr(>F)`[1],
        stringsAsFactors = FALSE
      )
    }
  )

  bind_rows(res) %>%
    mutate(
      FDR = p.adjust(
        P,
        method = "BH"
      )
    )
}

combine_increment_tables <- function(
  res_list
) {

  CP <- res_list$Combine_vs_PCOS %>%
    select(
      Feature,
      logFC_Combine_PCOS = logFC,
      P_Combine_PCOS = P.Value,
      FDR_Combine_PCOS = adj.P.Val
    )

  LP <- res_list$LCD_vs_PCOS %>%
    select(
      Feature,
      logFC_LCD_PCOS = logFC,
      P_LCD_PCOS = P.Value,
      FDR_LCD_PCOS = adj.P.Val
    )

  GP <- res_list$GFC_vs_PCOS %>%
    select(
      Feature,
      logFC_GFC_PCOS = logFC,
      P_GFC_PCOS = P.Value,
      FDR_GFC_PCOS = adj.P.Val
    )

  CL <- res_list$Combine_vs_LCD %>%
    select(
      Feature,
      logFC_Combine_LCD = logFC,
      P_Combine_LCD = P.Value,
      FDR_Combine_LCD = adj.P.Val
    )

  CG <- res_list$Combine_vs_GFC %>%
    select(
      Feature,
      logFC_Combine_GFC = logFC,
      P_Combine_GFC = P.Value,
      FDR_Combine_GFC = adj.P.Val
    )

  LG <- res_list$LCD_vs_GFC %>%
    select(
      Feature,
      logFC_LCD_GFC = logFC,
      P_LCD_GFC = P.Value,
      FDR_LCD_GFC = adj.P.Val
    )

  CP %>%
    inner_join(LP, by = "Feature") %>%
    inner_join(GP, by = "Feature") %>%
    inner_join(CL, by = "Feature") %>%
    inner_join(CG, by = "Feature") %>%
    inner_join(LG, by = "Feature") %>%
    mutate(

      IncrementClass = case_when(

        FDR_Combine_LCD < FDR_CANDIDATE &
          FDR_Combine_GFC < FDR_CANDIDATE &
          sign(logFC_Combine_LCD) ==
            sign(logFC_Combine_GFC) &
          abs(logFC_Combine_LCD) > 1e-8 &
          abs(logFC_Combine_GFC) > 1e-8 ~
          "Beyond both singles",

        FDR_Combine_LCD >= FDR_CANDIDATE &
          FDR_Combine_GFC < FDR_CANDIDATE ~
          "LCD-like",

        FDR_Combine_GFC >= FDR_CANDIDATE &
          FDR_Combine_LCD < FDR_CANDIDATE ~
          "GFC-like",

        FDR_Combine_LCD < FDR_CANDIDATE &
          FDR_Combine_GFC < FDR_CANDIDATE &
          sign(logFC_Combine_LCD) !=
            sign(logFC_Combine_GFC) ~
          "Intermediate",

        TRUE ~
          "No clear extra effect"
      ),

      IncrementClass = factor(
        IncrementClass,
        levels = c(
          "Beyond both singles",
          "LCD-like",
          "GFC-like",
          "Intermediate",
          "No clear extra effect"
        )
      ),

      CombinationSubtype = case_when(

        IncrementClass ==
          "Beyond both singles" &
          FDR_Combine_PCOS <
            FDR_CANDIDATE &
          FDR_LCD_PCOS >=
            FDR_CANDIDATE &
          FDR_GFC_PCOS >=
            FDR_CANDIDATE ~
          "Combination-specific",

        IncrementClass ==
          "Beyond both singles" &
          FDR_Combine_PCOS <
            FDR_CANDIDATE &
          sign(logFC_Combine_PCOS) ==
            sign(logFC_LCD_PCOS) &
          sign(logFC_Combine_PCOS) ==
            sign(logFC_GFC_PCOS) &
          abs(logFC_Combine_PCOS) >
            pmax(
              abs(logFC_LCD_PCOS),
              abs(logFC_GFC_PCOS)
            ) ~
          "Combination-enhanced",

        IncrementClass ==
          "Beyond both singles" ~
          "Beyond-both-other",

        TRUE ~
          NA_character_
      ),

      DirectPriority = pmax(
        FDR_Combine_LCD,
        FDR_Combine_GFC,
        na.rm = TRUE
      ),

      DirectMagnitude =
        abs(logFC_Combine_LCD) +
        abs(logFC_Combine_GFC)
    )
}

# ======================================================================
# 5. Plot helpers
# ======================================================================

make_increment_map <- function(
  df,
  label_col = "Feature",
  label_n = 6,
  highlight_requires_annotation = FALSE,
  max_background_points = Inf,
  max_highlight_per_class = Inf
) {

  dat <- df %>%
    mutate(
      LabelText = .data[[label_col]],

      LabelText = if (
        identical(
          label_col,
          "Feature"
        )
      ) {
        ifelse(
          is.na(LabelText) |
            LabelText == "",
          Feature,
          LabelText
        )
      } else {
        LabelText
      },

      IsAnnotated =
        !is.na(LabelText) &
        LabelText != "",

      IsHighlight =
        IncrementClass !=
          "No clear extra effect" &
        (
          !highlight_requires_annotation |
            IsAnnotated
        ),

      LabelTier = case_when(
        CombinationSubtype == "Combination-specific" ~ 1,
        CombinationSubtype == "Combination-enhanced" ~ 2,
        IncrementClass == "Beyond both singles" ~ 3,
        IncrementClass == "Intermediate" ~ 4,
        IncrementClass %in% c("LCD-like", "GFC-like") ~ 5,
        TRUE ~ 9
      )
    )

  background_df <- dat %>%
    filter(!IsHighlight)

  if (
    is.finite(max_background_points) &&
      nrow(background_df) > max_background_points
  ) {

    set.seed(20260930)

    background_plot <- background_df[
      sample(
        seq_len(nrow(background_df)),
        size = max_background_points,
        replace = FALSE
      ),
      ,
      drop = FALSE
    ]

  } else {
    background_plot <- background_df
  }

  highlight_df <- dat %>%
    filter(IsHighlight) %>%
    arrange(
      LabelTier,
      DirectPriority,
      desc(DirectMagnitude)
    )

  if (is.finite(max_highlight_per_class)) {

    highlight_plot <- highlight_df %>%
      group_by(IncrementClass) %>%
      slice_head(n = max_highlight_per_class) %>%
      ungroup()

  } else {
    highlight_plot <- highlight_df
  }

  # Label the most biologically informative candidates.
  # If no strict "Beyond both singles" Genus exists, the next-ranked
  # direct combination-effect candidates are still shown rather than
  # leaving the panel without any microbial names.
  label_df <- highlight_df %>%
    filter(IsAnnotated) %>%
    arrange(
      LabelTier,
      DirectPriority,
      desc(DirectMagnitude)
    ) %>%
    distinct(
      LabelText,
      .keep_all = TRUE
    ) %>%
    slice_head(n = label_n) %>%
    mutate(
      Label = stringr::str_wrap(
        LabelText,
        width = 24
      )
    )

  # Only show response classes that are actually present in the panel.
  desired_classes <- c(
    "Beyond both singles",
    "LCD-like",
    "GFC-like",
    "Intermediate"
  )

  present_classes <- desired_classes[
    desired_classes %in%
      as.character(
        unique(
          highlight_plot$IncrementClass
        )
      )
  ]

  lim <- max(
    abs(
      c(
        dat$logFC_Combine_LCD,
        dat$logFC_Combine_GFC
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
      colour = "grey60"
    ) +

    geom_point(
      data = background_plot,
      aes(
        x = logFC_Combine_LCD,
        y = logFC_Combine_GFC
      ),
      colour = "grey74",
      size = 0.85,
      alpha = 0.16
    ) +

    geom_point(
      data = highlight_plot,
      aes(
        x = logFC_Combine_LCD,
        y = logFC_Combine_GFC,
        colour = IncrementClass
      ),
      size = 1.95,
      alpha = 0.70
    ) +

    ggrepel::geom_text_repel(
      data = label_df,
      aes(
        x = logFC_Combine_LCD,
        y = logFC_Combine_GFC,
        label = Label,
        colour = IncrementClass
      ),
      size = 3.15,
      max.overlaps = Inf,
      box.padding = 0.65,
      point.padding = 0.35,
      min.segment.length = 0,
      segment.colour = "grey55",
      segment.size = 0.38,
      force = 2.5,
      max.time = 3,
      seed = 20260930,
      show.legend = FALSE
    ) +

    scale_colour_manual(
      values = increment_cols,
      breaks = present_classes,
      drop = TRUE
    ) +

    coord_cartesian(
      xlim = c(-lim, lim),
      ylim = c(-lim, lim),
      clip = "off"
    ) +

    guides(
      colour = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          size = 3,
          alpha = 1
        )
      )
    ) +

    labs(
      x = "Combine effect vs LCD",
      y = "Combine effect vs GFC",
      colour = NULL
    ) +

    theme_classic(
      base_size = 14
    ) +

    theme(
      legend.position = "top",
      legend.justification = "center",
      legend.box.just = "center",
      legend.text = element_text(size = 10.5),
      legend.margin = margin(0, 8, 5, 8),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(face = "bold"),
      plot.margin = margin(8, 14, 6, 14)
    )
}

make_volcano <- function(
  tt,
  title_text,
  label_col = "Feature",
  label_n = 8,
  annotated_only = FALSE
) {

  dat <- tt %>%
    mutate(
      Pplot = pmax(
        P.Value,
        1e-300
      ),
      Sig = adj.P.Val <
        FDR_CANDIDATE,
      Label = .data[[label_col]],
      LabelOK =
        !is.na(Label) &
        Label != ""
    )

  if (
    identical(
      label_col,
      "Feature"
    )
  ) {

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

    geom_point(
      data = dat %>%
        filter(!Sig),
      colour = "grey72",
      size = 1.0,
      alpha = 0.20
    ) +

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
      segment.size = 0.35,
      force = 2,
      max.time = 3,
      seed = 20260930
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
# 6. Read / aggregate 16S to Genus
# ======================================================================

if (!file.exists(otu_file)) {
  stop("Cannot find: ", otu_file)
}

if (!file.exists(taxonomy_file)) {
  stop("Cannot find: ", taxonomy_file)
}

otu_df <- readr::read_tsv(
  otu_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

names(otu_df)[1] <- "OTU_ID"

otu_sample_cols <- names(otu_df)[
  vapply(
    names(otu_df),
    is_sample_column,
    logical(1)
  )
]

otu_group <- group_from_sample(
  otu_sample_cols
)

bio_idx <- otu_group %in%
  analysis_groups

otu_bio_cols <- otu_sample_cols[
  bio_idx
]

otu_bio_group <- factor(
  otu_group[bio_idx],
  levels = analysis_groups
)

otu_mat <- safe_numeric_matrix(
  otu_df,
  otu_bio_cols,
  otu_df$OTU_ID
)

otu_mat[
  is.na(otu_mat)
] <- 0

otu_mat[
  otu_mat < 0
] <- 0

otu_mat <- otu_mat[
  rowSums(otu_mat) > 0,
  ,
  drop = FALSE
]

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
    "Taxonomy file requires at least 2 tab-separated columns."
  )
}

tax_df <- data.frame(
  OTU_ID = trimws(
    as.character(
      tax_raw[[1]]
    )
  ),
  Taxonomy = trimws(
    as.character(
      tax_raw[[2]]
    )
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Order = extract_rank(
      Taxonomy,
      "o__"
    ),
    Family = extract_rank(
      Taxonomy,
      "f__"
    ),
    Genus_raw = extract_rank(
      Taxonomy,
      "g__"
    ),
    Genus = clean_taxon(
      Genus_raw,
      Family,
      Order
    )
  )

map_idx <- match(
  rownames(otu_mat),
  tax_df$OTU_ID
)

taxonomy_qc <- data.frame(
  Metric = c(
    "OTUs_in_count_table",
    "OTUs_mapped",
    "Mapping_fraction"
  ),
  Value = c(
    nrow(otu_mat),
    sum(!is.na(map_idx)),
    mean(!is.na(map_idx))
  )
)

write_csv(
  taxonomy_qc,
  file.path(
    qc_dir,
    "taxonomy_mapping_QC.csv"
  )
)

otu_tax <- data.frame(
  OTU_ID = rownames(otu_mat),
  Genus = tax_df$Genus[map_idx],
  stringsAsFactors = FALSE
)

otu_tax$Genus[
  is.na(otu_tax$Genus) |
    otu_tax$Genus == ""
] <- "Unclassified"

otu_for_genus <- as.data.frame(
  otu_mat
) %>%
  rownames_to_column(
    "OTU_ID"
  ) %>%
  left_join(
    otu_tax,
    by = "OTU_ID"
  )

otu_for_genus$Genus[
  is.na(otu_for_genus$Genus)
] <- "Unclassified"

genus_df <- otu_for_genus %>%
  select(
    Genus,
    all_of(
      colnames(otu_mat)
    )
  ) %>%
  group_by(
    Genus
  ) %>%
  summarise(
    across(
      everything(),
      \(x) sum(
        x,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  )

genus_mat_all <- as.matrix(
  genus_df[, -1, drop = FALSE]
)

rownames(genus_mat_all) <- genus_df$Genus
storage.mode(genus_mat_all) <- "numeric"

genus_rel_all <- sweep(
  genus_mat_all,
  2,
  pmax(
    colSums(genus_mat_all),
    1
  ),
  FUN = "/"
)

prevalence <- rowMeans(
  genus_mat_all > 0
)

mean_rel <- rowMeans(
  genus_rel_all
)

keep_genus <- prevalence >=
  GENUS_PREVALENCE_MIN &
  mean_rel >=
  GENUS_REL_ABUND_MIN

# Reviewer 2 Major 6: harmonize the CLR denominator across all reported taxa.
# The following taxa were displayed in submitted Fig4H. OTU1043 has
# f__[Eubacterium]_coprostanoligenes_group and
# g__norank_f__[Eubacterium]_coprostanoligenes_group; the existing
# clean_taxon() correctly maps it to the family-level unclassified bin
# [Eubacterium]_coprostanoligenes_group_unclassified. Do not claim a
# formally resolved genus or silently rename other OTUs. A taxon that
# fails the generic prevalence/abundance filter is explicitly retained so
# that it and all other displayed taxa are transformed using ONE reference.
# No invented taxonomy labels, substituting genera or per-figure CLR allowed.
FIG4H_REQUIRED_GENERA <- c(
  "Alistipes", "Prevotellaceae_UCG-001", "Mucispirillum",
  "[Eubacterium]_coprostanoligenes_group_unclassified"
)
raw_present <- FIG4H_REQUIRED_GENERA %in% rownames(genus_mat_all)
target_audit <- data.frame(
  Genus = FIG4H_REQUIRED_GENERA,
  PresentInAggregatedCounts = raw_present,
  PassesGenericFilter = as.logical(keep_genus[
    match(FIG4H_REQUIRED_GENERA, rownames(genus_mat_all))
  ]),
  stringsAsFactors = FALSE
)
target_audit$PassesGenericFilter[is.na(target_audit$PassesGenericFilter)] <- FALSE
target_audit$ExplicitlyRetained <- target_audit$PresentInAggregatedCounts &
  !target_audit$PassesGenericFilter
readr::write_csv(target_audit,
  file.path(qc_dir, "Fig4H_TaxonRetention_Audit.csv"))

if (any(!raw_present)) {
  missing_names <- FIG4H_REQUIRED_GENERA[!raw_present]
  # The diagnostic contains taxonomy-related labels but NO automatic remapping.
  candidate_names <- rownames(genus_mat_all)[
    grepl("coprostanoligenes|eubacterium|alistipes|mucispirillum|prevotellaceae",
          rownames(genus_mat_all), ignore.case=TRUE)
  ]
  writeLines(c(
    paste("Required genera absent from aggregated count matrix:",
          paste(missing_names,collapse="; ")),
    "Candidate names in source taxonomy:", candidate_names
  ), file.path(qc_dir, "Fig4H_Taxonomy_Review_Required.txt"))
  stop("Required submitted Fig4H genus not found in RAW genus aggregation: ",
       paste(missing_names,collapse=", "),
       ". Check ", file.path(qc_dir, "Fig4H_Taxonomy_Review_Required.txt"),
       ". Do NOT relabel another genus without checking its OTU taxonomy.")
}
keep_genus[FIG4H_REQUIRED_GENERA] <- TRUE

genus_mat <- genus_mat_all[
  keep_genus,
  ,
  drop = FALSE
]

genus_rel <- genus_rel_all[
  keep_genus,
  ,
  drop = FALSE
]

genus_clr <- clr_transform(
  genus_mat,
  pseudocount = 0.5
)

# === MAJOR 6: canonical genus CLR matrix for Figures 3, 4 and 5 ===
# Aggregate OTUs and filter prevalence >=0.20 and mean relative abundance
# >=0.0001, explicitly retaining the four submitted Fig4H taxa when present;
# then apply log(count+0.5) CLR to the full retained taxon set.
# Never run CLR again on a selected subset of genera.
MASTER_CLR_FILE <- file.path(table_dir, "genus_CLR_master.csv")
master_df <- data.frame(Genus=rownames(genus_clr), genus_clr,
                        check.names=FALSE, stringsAsFactors=FALSE)
readr::write_csv(master_df, MASTER_CLR_FILE)
readr::write_csv(
  data.frame(Genus=rownames(genus_mat_all),
             Prevalence=as.numeric(prevalence),
             MeanRelativeAbundance=as.numeric(mean_rel),
             RetainedForCLR=as.logical(keep_genus),
             ExplicitFig4HRetention=rownames(genus_mat_all) %in% FIG4H_REQUIRED_GENERA),
  file.path(qc_dir, "Genus_CLR_FeatureFilter_QC.csv"))
readr::write_csv(
  data.frame(Metric=c("Samples", "Genera_before_filter", "Genera_after_filter",
                       "Prevalence_min", "Mean_relative_abundance_min", "Count_pseudocount",
                       "Fig4H_explicitly_retained_below_filter"),
             Value=as.character(c(ncol(genus_clr),nrow(genus_mat_all),nrow(genus_clr),
                                  GENUS_PREVALENCE_MIN,GENUS_REL_ABUND_MIN,0.5,
                                  sum(target_audit$ExplicitlyRetained)))),
  file.path(qc_dir, "Genus_CLR_Method_QC.csv"))
for (target in c("Mucispirillum", "Alistipes")) {
  if (!(target %in% rownames(genus_clr)))
    warning("Review-target genus absent after filtering: ",target)
}

# ======================================================================
# 7. Fig3A | Global Genus-level PCoA + PERMANOVA
# ======================================================================

bray <- vegan::vegdist(
  t(genus_rel_all),
  method = "bray"
)

pcoa <- cmdscale(
  bray,
  k = 2,
  eig = TRUE,
  add = TRUE
)

pcoa_df <- data.frame(
  Sample = rownames(
    pcoa$points
  ),
  Axis1 = pcoa$points[, 1],
  Axis2 = pcoa$points[, 2],
  stringsAsFactors = FALSE
) %>%
  mutate(
    Group = factor(
      group_from_sample(Sample),
      levels = analysis_groups
    )
  )

eig_pos <- pcoa$eig[
  pcoa$eig > 0
]

var1 <- 100 *
  pcoa$eig[1] /
  sum(eig_pos)

var2 <- 100 *
  pcoa$eig[2] /
  sum(eig_pos)

meta_16s <- data.frame(
  Sample = colnames(
    genus_rel_all
  ),
  Group = factor(
    group_from_sample(
      colnames(genus_rel_all)
    ),
    levels = analysis_groups
  )
)

rownames(meta_16s) <- meta_16s$Sample

global_perm_16s <- vegan::adonis2(
  bray ~ Group,
  data = meta_16s,
  permutations = 999
)

capture.output(
  global_perm_16s,
  file = file.path(
    table_dir,
    "Fig3A_Genus_global_PERMANOVA.txt"
  )
)

pair_perm_16s <- pairwise_permanova(
  bray,
  meta_16s,
  pairs = list(
    c("Combine", "PCOS"),
    c("Combine", "LCD"),
    c("Combine", "GFC")
  )
)

write_csv(
  pair_perm_16s,
  file.path(
    table_dir,
    "Fig3A_Genus_pairwise_PERMANOVA.csv"
  )
)

bd_16s <- vegan::betadisper(
  bray,
  meta_16s$Group
)

capture.output(
  anova(bd_16s),
  file = file.path(
    table_dir,
    "Fig3A_Genus_BETADISPER.txt"
  )
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
    show.legend = FALSE
  ) +

  geom_point(
    size = 3.0,
    alpha = 0.90
  ) +

  scale_colour_manual(
    values = group_cols,
    breaks = analysis_groups
  ) +

  scale_shape_manual(
    values = group_shapes,
    breaks = analysis_groups
  ) +

  labs(
    x = sprintf(
      "PCoA1 (%.1f%%)",
      var1
    ),
    y = sprintf(
      "PCoA2 (%.1f%%)",
      var2
    ),
    colour = NULL,
    shape = NULL
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    legend.position = "top",
    legend.justification = "center",
    axis.text = element_text(
      colour = "black"
    ),
    axis.title = element_text(
      face = "bold"
    )
  )

save_plot(
  p_A,
  "Fig3A_Genus_PCoA_4Groups",
  main_dir,
  width = 5.8,
  height = 4.8
)

# ======================================================================
# 8. Genus limma + combination increment
# ======================================================================

genus_res <- make_limma_results_4groups(
  genus_clr,
  otu_bio_group
)

for (nm in names(genus_res)) {

  write_csv(
    genus_res[[nm]],
    file.path(
      table_dir,
      paste0(
        "Genus_",
        nm,
        "_limma.csv"
      )
    )
  )
}

genus_increment <- combine_increment_tables(
  genus_res
) %>%
  mutate(
    LiteraturePriority =
      Feature %in%
      literature_priority_genera
  )

write_csv(
  genus_increment,
  file.path(
    table_dir,
    "Genus_combination_increment_summary.csv"
  )
)

# ======================================================================
# 9. Read fecal metabolomics
# ======================================================================

if (!file.exists(fecal_file)) {
  stop("Cannot find: ", fecal_file)
}

fecal_df <- readr::read_tsv(
  fecal_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

fecal_sample_cols_all <- names(fecal_df)[
  vapply(
    names(fecal_df),
    is_sample_column,
    logical(1)
  )
]

fecal_group_all <- group_from_sample(
  fecal_sample_cols_all
)

fecal_bio_idx <- fecal_group_all %in%
  analysis_groups

fecal_bio_cols <- fecal_sample_cols_all[
  fecal_bio_idx
]

fecal_bio_group <- factor(
  fecal_group_all[
    fecal_bio_idx
  ],
  levels = analysis_groups
)

if (!"metab_id" %in%
    names(fecal_df)) {
  stop(
    "Column 'metab_id' is required."
  )
}

if (!"Metabolite" %in%
    names(fecal_df)) {
  stop(
    "Column 'Metabolite' is required."
  )
}

feature_id <- as.character(
  fecal_df$metab_id
)

missing_feature <- is.na(
  feature_id
) |
  trimws(feature_id) == ""

feature_id[
  missing_feature
] <- paste0(
  "Feature_",
  which(missing_feature)
)

feature_id <- make.unique(
  feature_id
)

metabolite_raw <- trimws(
  as.character(
    fecal_df$Metabolite
  )
)

annotated_flag <- is_valid_metabolite_annotation(
  metabolite_raw
)

display_name <- rep(
  NA_character_,
  length(metabolite_raw)
)

display_name[
  annotated_flag
] <- metabolite_raw[
  annotated_flag
]

# Preserve chemical taxonomy for literature-informed module mapping.
get_optional_col <- function(
  df,
  candidates
) {

  hit <- candidates[
    candidates %in%
      names(df)
  ]

  if (length(hit) == 0) {
    return(
      rep(
        NA_character_,
        nrow(df)
      )
    )
  }

  as.character(
    df[[hit[1]]]
  )
}

metab_annotation <- data.frame(
  Feature = feature_id,
  DisplayName = display_name,
  Annotated = annotated_flag,

  ChemicalSuperclass = get_optional_col(
    fecal_df,
    c(
      "chemical_superclass",
      "chemical superclass",
      "Chemical_superclass"
    )
  ),

  ChemicalClass = get_optional_col(
    fecal_df,
    c(
      "chemical_class",
      "chemical class",
      "Chemical_class"
    )
  ),

  ChemicalSubclass = get_optional_col(
    fecal_df,
    c(
      "chemical_subclass",
      "chemical subclass",
      "Chemical_subclass"
    )
  ),

  KEGG = get_optional_col(
    fecal_df,
    c(
      "KEGG Compound ID",
      "KEGG_Compound_ID",
      "KEGG Compound"
    )
  ),

  stringsAsFactors = FALSE
)

write_csv(
  metab_annotation,
  file.path(
    table_dir,
    "FecalMetabolite_annotation_map.csv"
  )
)

fecal_mat_raw <- safe_numeric_matrix(
  fecal_df,
  fecal_bio_cols,
  feature_id
)

missing_rate <- rowMeans(
  is.na(fecal_mat_raw)
)

keep_metab <- missing_rate <=
  0.20

fecal_mat_raw <- fecal_mat_raw[
  keep_metab,
  ,
  drop = FALSE
]

impute_feature <- function(x) {

  if (all(is.na(x))) {
    return(
      rep(
        0,
        length(x)
      )
    )
  }

  x[is.na(x)] <- median(
    x,
    na.rm = TRUE
  )

  x
}

fecal_mat_imp <- t(
  apply(
    fecal_mat_raw,
    1,
    impute_feature
  )
)

colnames(fecal_mat_imp) <- colnames(
  fecal_mat_raw
)

rownames(fecal_mat_imp) <- rownames(
  fecal_mat_raw
)

finite_values <- as.numeric(
  fecal_mat_imp[
    is.finite(fecal_mat_imp)
  ]
)

q50 <- median(
  finite_values,
  na.rm = TRUE
)

q99 <- quantile(
  finite_values,
  0.99,
  na.rm = TRUE
)

already_log_like <- q99 < 50 &&
  q50 < 20

if (already_log_like) {

  fecal_mat <- fecal_mat_imp
  metab_transform_note <-
    "No additional log transform; input appears log-like."

} else {

  if (
    min(
      fecal_mat_imp,
      na.rm = TRUE
    ) < 0
  ) {

    shift <- abs(
      min(
        fecal_mat_imp,
        na.rm = TRUE
      )
    ) + 1

  } else {
    shift <- 1
  }

  fecal_mat <- log2(
    fecal_mat_imp + shift
  )

  metab_transform_note <- paste0(
    "Applied log2(x + ",
    signif(shift, 4),
    ")."
  )
}

metab_sd <- apply(
  fecal_mat,
  1,
  sd,
  na.rm = TRUE
)

fecal_mat <- fecal_mat[
  is.finite(metab_sd) &
    metab_sd > 0,
  ,
  drop = FALSE
]

# ======================================================================
# 10. Fig3B | Global fecal metabolome PCA + PERMANOVA
# ======================================================================

pca_fecal <- prcomp(
  t(fecal_mat),
  center = TRUE,
  scale. = TRUE
)

var_exp <- 100 *
  pca_fecal$sdev^2 /
  sum(
    pca_fecal$sdev^2
  )

pca_fecal_df <- data.frame(
  Sample = rownames(
    pca_fecal$x
  ),
  PC1 = pca_fecal$x[, 1],
  PC2 = pca_fecal$x[, 2],
  stringsAsFactors = FALSE
) %>%
  mutate(
    Group = factor(
      group_from_sample(Sample),
      levels = analysis_groups
    )
  )

# Euclidean distance on feature-wise standardized metabolome,
# matching the geometry used in PCA.
fecal_z <- t(
  scale(
    t(fecal_mat)
  )
)

fecal_z[
  !is.finite(fecal_z)
] <- 0

fecal_dist <- dist(
  t(fecal_z),
  method = "euclidean"
)

meta_fecal <- data.frame(
  Sample = colnames(fecal_z),
  Group = factor(
    group_from_sample(
      colnames(fecal_z)
    ),
    levels = analysis_groups
  )
)

rownames(meta_fecal) <- meta_fecal$Sample

global_perm_fecal <- vegan::adonis2(
  fecal_dist ~ Group,
  data = meta_fecal,
  permutations = 999
)

capture.output(
  global_perm_fecal,
  file = file.path(
    table_dir,
    "Fig3B_Fecal_global_PERMANOVA.txt"
  )
)

pair_perm_fecal <- pairwise_permanova(
  fecal_dist,
  meta_fecal,
  pairs = list(
    c("Combine", "PCOS"),
    c("Combine", "LCD"),
    c("Combine", "GFC")
  )
)

write_csv(
  pair_perm_fecal,
  file.path(
    table_dir,
    "Fig3B_Fecal_pairwise_PERMANOVA.csv"
  )
)

bd_fecal <- vegan::betadisper(
  fecal_dist,
  meta_fecal$Group
)

capture.output(
  anova(bd_fecal),
  file = file.path(
    table_dir,
    "Fig3B_Fecal_BETADISPER.txt"
  )
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
    show.legend = FALSE
  ) +

  geom_point(
    size = 3.0,
    alpha = 0.90
  ) +

  scale_colour_manual(
    values = group_cols,
    breaks = analysis_groups
  ) +

  scale_shape_manual(
    values = group_shapes,
    breaks = analysis_groups
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
    legend.justification = "center",
    axis.text = element_text(
      colour = "black"
    ),
    axis.title = element_text(
      face = "bold"
    )
  )

save_plot(
  p_B,
  "Fig3B_FecalMetabolome_PCA_4Groups",
  main_dir,
  width = 5.8,
  height = 4.8
)

# ======================================================================
# 11. Fecal limma + combination increment
# ======================================================================

fecal_res <- make_limma_results_4groups(
  fecal_mat,
  fecal_bio_group
)

for (nm in names(fecal_res)) {

  fecal_res[[nm]] <- fecal_res[[nm]] %>%
    left_join(
      metab_annotation,
      by = "Feature"
    )

  write_csv(
    fecal_res[[nm]],
    file.path(
      table_dir,
      paste0(
        "FecalMetabolite_",
        nm,
        "_limma.csv"
      )
    )
  )
}

fecal_increment <- combine_increment_tables(
  fecal_res
) %>%
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
  fecal_increment,
  file.path(
    table_dir,
    "FecalMetabolite_combination_increment_summary.csv"
  )
)

# ======================================================================
# 12. Fig3C1/C2 | Combination-increment screening
# ======================================================================

p_C1 <- make_increment_map(
  genus_increment,
  label_col = "Feature",
  label_n = 7,
  highlight_requires_annotation = FALSE,
  max_background_points = Inf,
  max_highlight_per_class = Inf
) +

  labs(
    x = "Genus effect: Combine vs LCD (CLR scale)",
    y = "Genus effect: Combine vs GFC (CLR scale)"
  )

save_plot(
  p_C1,
  "Fig3C1_Genus_CombinationIncrement_Map",
  main_dir,
  width = 7.0,
  height = 5.5
)

p_C2 <- make_increment_map(
  fecal_increment,
  label_col = "DisplayName",
  label_n = 6,
  highlight_requires_annotation = TRUE,
  max_background_points = 2500,
  max_highlight_per_class = 120
) +

  labs(
    x = "Metabolite effect: Combine vs LCD",
    y = "Metabolite effect: Combine vs GFC"
  )

save_plot(
  p_C2,
  "Fig3C2_FecalMetabolite_CombinationIncrement_Map",
  main_dir,
  width = 7.0,
  height = 5.5
)

# ======================================================================
# 13. Select key Genus
#
# Candidate-selection principle for the MAIN figure
# ----------------------------------------------------------------------
# Tier 1: Combination-specific
# Tier 2: Combination-enhanced
# Tier 3: Other "Beyond both singles"
# Tier 4: Intermediate (Combine differs from both singles but lies between)
# Tier 5: LCD-like / GFC-like (used only to complete the representative set)
#
# Within a tier:
#   smaller direct-comparison FDR -> larger direct effect -> literature tie-break.
#
# Therefore the main representative taxa explicitly prioritize what the
# combination treatment changes beyond the single treatments.
# ======================================================================

selection_tier <- function(
  subtype,
  increment_class
) {

  out <- rep(
    9,
    length(increment_class)
  )

  out[
    !is.na(subtype) &
      subtype == "Combination-specific"
  ] <- 1

  out[
    !is.na(subtype) &
      subtype == "Combination-enhanced"
  ] <- 2

  out[
    increment_class == "Beyond both singles" &
      out == 9
  ] <- 3

  out[
    increment_class == "Intermediate"
  ] <- 4

  out[
    increment_class %in%
      c("LCD-like", "GFC-like")
  ] <- 5

  out
}

selection_reason <- function(
  tier,
  increment_class,
  subtype
) {

  case_when(
    tier == 1 ~
      "Combination-specific: new response beyond both single treatments",
    tier == 2 ~
      "Combination-enhanced: same-direction response stronger than both singles",
    tier == 3 ~
      "Beyond both singles: Combine lies outside both single-treatment states",
    tier == 4 ~
      "Intermediate but directly different from both single treatments",
    tier == 5 &
      increment_class == "LCD-like" ~
      "Fallback representative: Combine resembles LCD but differs from GFC",
    tier == 5 &
      increment_class == "GFC-like" ~
      "Fallback representative: Combine resembles GFC but differs from LCD",
    TRUE ~
      "Lower-priority direct combination effect"
  )
}

genus_candidates <- genus_increment %>%
  mutate(
    SelectionTier = selection_tier(
      CombinationSubtype,
      as.character(
        IncrementClass
      )
    ),
    SelectionReason = selection_reason(
      SelectionTier,
      as.character(
        IncrementClass
      ),
      CombinationSubtype
    )
  ) %>%
  filter(
    SelectionTier < 9
  ) %>%
  arrange(
    SelectionTier,
    DirectPriority,
    desc(DirectMagnitude),
    desc(LiteraturePriority)
  )

if (nrow(genus_candidates) < N_KEY_GENUS) {

  genus_candidates <- genus_increment %>%
    mutate(
      SelectionTier = selection_tier(
        CombinationSubtype,
        as.character(
          IncrementClass
        )
      ),
      SelectionReason = selection_reason(
        SelectionTier,
        as.character(
          IncrementClass
        ),
        CombinationSubtype
      )
    ) %>%
    arrange(
      SelectionTier,
      DirectPriority,
      desc(DirectMagnitude),
      desc(LiteraturePriority)
    )
}

key_genus <- genus_candidates %>%
  slice_head(
    n = N_KEY_GENUS
  )

write_csv(
  key_genus,
  file.path(
    table_dir,
    "MAIN_Key_Genera.csv"
  )
)

# ======================================================================
# 14. Fig3D | Six representative Genus across four groups
# ======================================================================

genus_box_df <- as.data.frame(
  t(
    genus_clr[
      key_genus$Feature,
      ,
      drop = FALSE
    ]
  )
) %>%
  rownames_to_column(
    "Sample"
  ) %>%
  mutate(
    Group = factor(
      group_from_sample(Sample),
      levels = analysis_groups
    )
  ) %>%
  pivot_longer(
    cols = -c(
      Sample,
      Group
    ),
    names_to = "Genus",
    values_to = "CLR"
  ) %>%
  mutate(
    GenusDisplay = stringr::str_wrap(
      Genus,
      width = 22
    ),
    GenusDisplay = factor(
      GenusDisplay,
      levels = unique(
        stringr::str_wrap(
          key_genus$Feature,
          width = 22
        )
      )
    )
  )

p_D <- ggplot(
  genus_box_df,
  aes(
    x = Group,
    y = CLR,
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
    size = 2.2,
    alpha = 0.90
  ) +

  facet_wrap(
    ~ GenusDisplay,
    scales = "free_y",
    ncol = 3,
    drop = FALSE
  ) +

  scale_fill_manual(
    values = group_cols,
    guide = "none"
  ) +

  scale_colour_manual(
    values = group_cols,
    guide = "none"
  ) +

  labs(
    x = NULL,
    y = "CLR abundance"
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 9.8,
      lineheight = 0.95
    ),
    axis.text.x = element_text(
      angle = 40,
      hjust = 1,
      size = 9.2,
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    ),
    axis.title.y = element_text(
      face = "bold"
    ),
    panel.spacing = grid::unit(
      0.9,
      "lines"
    )
  )

save_plot(
  p_D,
  "Fig3D_KeyGenus_Abundance",
  main_dir,
  width = 10.0,
  height = 6.4
)

# ======================================================================
# 15. Select six key annotated metabolites
#
# Uses the same combination-focused tiering principle as the Genus panel.
# Duplicate metabolite names are reduced to one representative LC-MS feature.
# ======================================================================

metab_candidates <- fecal_increment %>%
  filter(
    Annotated,
    !is.na(DisplayName),
    DisplayName != ""
  ) %>%
  mutate(
    SelectionTier = selection_tier(
      CombinationSubtype,
      as.character(
        IncrementClass
      )
    ),
    SelectionReason = selection_reason(
      SelectionTier,
      as.character(
        IncrementClass
      ),
      CombinationSubtype
    )
  ) %>%
  filter(
    SelectionTier < 9
  ) %>%
  arrange(
    SelectionTier,
    DirectPriority,
    desc(DirectMagnitude)
  ) %>%
  distinct(
    DisplayName,
    .keep_all = TRUE
  )

if (nrow(metab_candidates) < N_KEY_METAB) {

  metab_candidates <- fecal_increment %>%
    filter(
      Annotated,
      !is.na(DisplayName),
      DisplayName != ""
    ) %>%
    mutate(
      SelectionTier = selection_tier(
        CombinationSubtype,
        as.character(
          IncrementClass
        )
      ),
      SelectionReason = selection_reason(
        SelectionTier,
        as.character(
          IncrementClass
        ),
        CombinationSubtype
      )
    ) %>%
    arrange(
      SelectionTier,
      DirectPriority,
      desc(DirectMagnitude)
    ) %>%
    distinct(
      DisplayName,
      .keep_all = TRUE
    )
}

key_metab <- metab_candidates %>%
  slice_head(
    n = N_KEY_METAB
  )

write_csv(
  key_metab,
  file.path(
    table_dir,
    "MAIN_Key_Annotated_Metabolites.csv"
  )
)

# ======================================================================
# 16. Fig3E1 | Six representative annotated metabolites
# ======================================================================

metab_box_mat <- fecal_mat[
  key_metab$Feature,
  ,
  drop = FALSE
]

rownames(
  metab_box_mat
) <- make.unique(
  key_metab$DisplayName
)

metab_box_df <- as.data.frame(
  t(metab_box_mat)
) %>%
  rownames_to_column(
    "Sample"
  ) %>%
  mutate(
    Group = factor(
      group_from_sample(Sample),
      levels = analysis_groups
    )
  ) %>%
  pivot_longer(
    cols = -c(
      Sample,
      Group
    ),
    names_to = "Metabolite",
    values_to = "Abundance"
  ) %>%
  mutate(
    MetaboliteDisplay = stringr::str_wrap(
      Metabolite,
      width = 26
    ),
    MetaboliteDisplay = factor(
      MetaboliteDisplay,
      levels = unique(
        stringr::str_wrap(
          key_metab$DisplayName,
          width = 26
        )
      )
    )
  )

p_E1 <- ggplot(
  metab_box_df,
  aes(
    x = Group,
    y = Abundance,
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
    size = 2.1,
    alpha = 0.90
  ) +

  facet_wrap(
    ~ MetaboliteDisplay,
    scales = "free_y",
    ncol = 3,
    drop = FALSE
  ) +

  scale_fill_manual(
    values = group_cols,
    guide = "none"
  ) +

  scale_colour_manual(
    values = group_cols,
    guide = "none"
  ) +

  labs(
    x = NULL,
    y = "Metabolite abundance"
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 9.4,
      lineheight = 0.95
    ),
    axis.text.x = element_text(
      angle = 40,
      hjust = 1,
      size = 9.2,
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    ),
    axis.title.y = element_text(
      face = "bold"
    ),
    panel.spacing = grid::unit(
      0.9,
      "lines"
    )
  )

save_plot(
  p_E1,
  "Fig3E1_KeyMetabolite_Abundance",
  main_dir,
  width = 10.2,
  height = 6.5
)

# ======================================================================
# 17. Fig3E2 | Literature-informed metabolic modules
#
# This is NOT KEGG enrichment.
# It maps annotated metabolite names / chemical taxonomy into predefined
# biologically relevant modules.
# ======================================================================

# Literature basis for the six candidate metabolic modules.
# These references justify PRIORITIZING the modules; the user's own
# metabolite data determine whether a module is actually retained.
literature_basis_modules <- data.frame(
  Module = c(
    "Bile acid metabolism",
    "Steroid / androgen metabolism",
    "Tryptophan / indole",
    "Fatty acid / lipid",
    "Sphingolipid / ceramide",
    "Amino acid / peptide"
  ),
  Rationale = c(
    "Gut microbiota-bile acid signaling is mechanistically linked to PCOS; B. vulgatus, GDCA/TUDCA and IL-22 form a validated gut-ovary axis.",
    "Ketogenic/low-carbohydrate intervention in PCOS rats altered steroid hormone biosynthesis and androgen-related metabolites including testosterone and 7alpha-hydroxytestosterone.",
    "Gut-derived tryptophan/indole metabolites are repeatedly implicated in microbiota-host signaling, inflammation and metabolic regulation in PCOS.",
    "PCOS metabolomics consistently identifies fatty-acid and broader lipid metabolic abnormalities, especially in insulin-resistant phenotypes.",
    "Human PCOS lipidomics/sphingolipidomics reports altered ceramide, sphingomyelin and related sphingolipid species.",
    "PCOS metabolomics repeatedly identifies amino-acid/peptide abnormalities, including essential and branched-chain amino acids linked to insulin resistance."
  ),
  Reference = c(
    "Qi et al., Nature Medicine 2019; DOI:10.1038/s41591-019-0509-0; PMID:31332392",
    "Effect of the ketogenic diet on gut microbiome composition and metabolomics in PCOS rats; PMID:37441825",
    "Qin et al., Front Endocrinol 2025; PMID:41480369; and intestinal microbial metabolite reviews",
    "Combining nontargeted and targeted metabolomics in PCOS; PMID:28521878; Untargeted serum metabolomics; PMID:34416878",
    "Alterations of Sphingolipid Metabolism in Different Types of PCOS; PMID:30824725; Ceramide subclasses study; PMID:31793360",
    "Combining nontargeted and targeted metabolomics in PCOS; PMID:28521878"
  ),
  stringsAsFactors = FALSE
)

write_csv(
  literature_basis_modules,
  file.path(
    table_dir,
    "Fig3E2_Literature_Basis_for_Metabolic_Modules.csv"
  )
)

module_priority <- literature_basis_modules$Module

annotation_text <- paste(
  ifelse(
    is.na(metab_annotation$DisplayName),
    "",
    metab_annotation$DisplayName
  ),
  ifelse(
    is.na(metab_annotation$ChemicalSuperclass),
    "",
    metab_annotation$ChemicalSuperclass
  ),
  ifelse(
    is.na(metab_annotation$ChemicalClass),
    "",
    metab_annotation$ChemicalClass
  ),
  ifelse(
    is.na(metab_annotation$ChemicalSubclass),
    "",
    metab_annotation$ChemicalSubclass
  ),
  sep = " | "
)

annotation_text_lower <- tolower(
  annotation_text
)

module <- case_when(

  grepl(
    "bile acid|cholic|cholate|deoxychol|chenodeoxychol|ursodeoxychol|glycoch|tauroch",
    annotation_text_lower
  ) ~ "Bile acid metabolism",

  grepl(
    "steroid|testosterone|andro|estradiol|estrone|progester|pregnen|corticost",
    annotation_text_lower
  ) ~ "Steroid / androgen metabolism",

  grepl(
    "ceramide|sphing|sphingos|sphingan|sphingomyelin",
    annotation_text_lower
  ) ~ "Sphingolipid / ceramide",

  grepl(
    "tryptophan|indole|kynuren|serotonin|5-hydroxytrypt",
    annotation_text_lower
  ) ~ "Tryptophan / indole",

  grepl(
    "linole|arachidon|fatty acid|acylcarnit|oleic|palmit|glycerophosph|phosphatid",
    annotation_text_lower
  ) ~ "Fatty acid / lipid",

  grepl(
    "amino acid|peptide|dipeptide|tripeptide",
    annotation_text_lower
  ) ~ "Amino acid / peptide",

  TRUE ~ NA_character_
)

metab_module_map <- metab_annotation %>%
  mutate(
    Module = module
  )

write_csv(
  metab_module_map,
  file.path(
    table_dir,
    "LiteratureInformed_MetabolicModule_Mapping.csv"
  )
)

module_features <- metab_module_map %>%
  filter(
    Annotated,
    !is.na(Module),
    Feature %in%
      rownames(fecal_mat)
  )

module_counts <- module_features %>%
  count(
    Module,
    name = "N_features"
  ) %>%
  arrange(
    desc(N_features)
  )

write_csv(
  module_counts,
  file.path(
    table_dir,
    "LiteratureInformed_MetabolicModule_FeatureCounts.csv"
  )
)

# Keep literature-prioritized modules that are represented in the data.
# Prefer modules with >=2 annotated features; if fewer than three modules
# meet this criterion, allow single-feature modules but flag them in the
# feature-count table for cautious interpretation.
modules_keep <- module_priority[
  module_priority %in%
    module_counts$Module[
      module_counts$N_features >= 2
    ]
]

if (length(modules_keep) < 3) {

  modules_keep <- unique(
    c(
      modules_keep,
      module_priority[
        module_priority %in%
          module_counts$Module[
            module_counts$N_features >= 1
          ]
      ]
    )
  )
}

modules_keep <- head(
  modules_keep,
  6
)

module_score_rows <- list()
k_module <- 1

for (m in modules_keep) {

  f <- module_features$Feature[
    module_features$Module == m
  ]

  f <- intersect(
    f,
    rownames(fecal_mat)
  )

  if (length(f) == 0) {
    next
  }

  z <- t(
    scale(
      t(
        fecal_mat[
          f,
          ,
          drop = FALSE
        ]
      )
    )
  )

  z[
    !is.finite(z)
  ] <- 0

  score <- colMeans(
    z,
    na.rm = TRUE
  )

  module_score_rows[[k_module]] <- data.frame(
    Sample = names(score),
    Module = m,
    Score = as.numeric(score),
    stringsAsFactors = FALSE
  )

  k_module <- k_module + 1
}

module_score_df <- bind_rows(
  module_score_rows
) %>%
  mutate(
    Group = factor(
      group_from_sample(Sample),
      levels = analysis_groups
    )
  )

write_csv(
  module_score_df,
  file.path(
    table_dir,
    "LiteratureInformed_MetabolicModule_Scores.csv"
  )
)

if (nrow(module_score_df) > 0) {

  p_E2 <- ggplot(
    module_score_df,
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
      aes(
        colour = Group
      ),
      width = 0.08,
      size = 2.0,
      alpha = 0.88
    ) +

    facet_wrap(
      ~ Module,
      scales = "free_y",
      ncol = 2,
      labeller = labeller(Module = function(x) {
        # Display-only label wrapping; underlying module names are unchanged.
        x <- gsub(" / ", " /\n", as.character(x), fixed = TRUE)
        x <- gsub("Steroid /\nandrogen metabolism", "Steroid / androgen\nmetabolism", x, fixed = TRUE)
        x <- gsub("Sphingolipid /\nceramide", "Sphingolipid /\nceramide", x, fixed = TRUE)
        x
      })
    ) +

    scale_fill_manual(
      values = group_cols,
      guide = "none"
    ) +

    scale_colour_manual(
      values = group_cols,
      guide = "none"
    ) +

    labs(
      x = NULL,
      y = "Module score (mean z-score)"
    ) +

    theme_classic(
      base_size = 14
    ) +

    theme(
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        size = 8.4,
        lineheight = 0.96,
        margin = margin(t = 6, b = 6)
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
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
    p_E2,
    "Fig3E2_LiteratureInformed_MetabolicModule_Scores",
    main_dir,
    width = 5.2,  # adjusted design width; save_plot() applies the 70% scaling
    height = max(
      9,
      3.5 *
        ceiling(
          length(
            unique(
              module_score_df$Module
            )
          ) /
            2
        )
    )
  )
}

# ======================================================================
# 18. Fig3F | Candidate Genus-metabolite axes
#
# Group-adjusted Spearman correlation:
# first remove group mean differences, then correlate residuals.
# Exploratory because sample size is small.
# ======================================================================

matched_samples <- intersect(
  colnames(genus_clr),
  colnames(fecal_mat)
)

cor_genus <- genus_clr[
  key_genus$Feature,
  matched_samples,
  drop = FALSE
]

cor_metab <- fecal_mat[
  key_metab$Feature,
  matched_samples,
  drop = FALSE
]

rownames(cor_metab) <- make.unique(
  key_metab$DisplayName
)

matched_group <- factor(
  group_from_sample(
    matched_samples
  ),
  levels = analysis_groups
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

  for (i in seq_len(nrow(mat))) {

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

g_resid <- residualize_by_group(
  cor_genus,
  matched_group
)

m_resid <- residualize_by_group(
  cor_metab,
  matched_group
)

cor_rows <- list()
k_cor <- 1

for (g in rownames(g_resid)) {

  for (m in rownames(m_resid)) {

    ct <- suppressWarnings(
      cor.test(
        as.numeric(
          g_resid[g, ]
        ),
        as.numeric(
          m_resid[m, ]
        ),
        method = "spearman",
        exact = FALSE
      )
    )

    cor_rows[[k_cor]] <- data.frame(
      Genus = g,
      Metabolite = m,
      Rho = unname(
        ct$estimate
      ),
      P = ct$p.value,
      stringsAsFactors = FALSE
    )

    k_cor <- k_cor + 1
  }
}

cor_df <- bind_rows(
  cor_rows
) %>%
  mutate(
    FDR = p.adjust(
      P,
      method = "BH"
    ),
    Sig = case_when(
      FDR < 0.05 ~ "**",
      FDR < 0.10 ~ "*",
      TRUE ~ ""
    ),
    GenusDisplay = stringr::str_wrap(
      gsub("_", " ", Genus, fixed = TRUE),
      width = 22
    ),
    MetaboliteDisplay = stringr::str_wrap(
      gsub("-", "- ", Metabolite, fixed = TRUE),
      width = 17
    )
  )

write_csv(
  cor_df,
  file.path(
    table_dir,
    "Fig3F_Candidate_GenusMetabolite_Correlations.csv"
  )
)

p_F <- ggplot(
  cor_df,
  aes(
    x = MetaboliteDisplay,
    y = GenusDisplay,
    fill = Rho
  )
) +

  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +

  geom_text(
    aes(label = Sig),
    size = 4.0
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
    y = NULL
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    axis.text.x = element_text(
      angle = 55,
      hjust = 1,
      vjust = 1,
      size = 10,
      colour = "black",
      lineheight = 0.92
    ),
    axis.text.y = element_text(
      size = 10,
      colour = "black"
    ),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.key.height = grid::unit(1.6, "cm"),
    plot.margin = margin(t = 8, r = 12, b = 14, l = 12)
  )

save_plot(
  p_F,
  "Fig3F_Candidate_GenusMetabolite_Correlation",
  main_dir,
  width = 11.0,
  height = 6.2
)

# ======================================================================
# 19. Supplementary S3A | Shannon
# ======================================================================

shannon <- vegan::diversity(
  t(otu_mat),
  index = "shannon"
)

alpha_df <- data.frame(
  Sample = names(shannon),
  Shannon = as.numeric(shannon),
  Group = factor(
    group_from_sample(
      names(shannon)
    ),
    levels = analysis_groups
  )
)

kw_alpha <- kruskal.test(
  Shannon ~ Group,
  data = alpha_df
)

p_S3A <- ggplot(
  alpha_df,
  aes(
    x = Group,
    y = Shannon,
    fill = Group
  )
) +

  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.28,
    colour = "grey35"
  ) +

  geom_jitter(
    aes(colour = Group),
    width = 0.08,
    size = 2.5
  ) +

  scale_fill_manual(
    values = group_cols,
    guide = "none"
  ) +

  scale_colour_manual(
    values = group_cols,
    guide = "none"
  ) +

  labs(
    x = NULL,
    y = "Shannon diversity",
    subtitle = paste0(
      "Kruskal-Wallis P = ",
      signif(
        kw_alpha$p.value,
        3
      )
    )
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
    )
  )

save_plot(
  p_S3A,
  "FigS3A_16S_Alpha_Shannon",
  supp_dir,
  width = 5.2,
  height = 4.3
)

# ======================================================================
# 20. Supplementary S3B | Top10 Genus composition
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

comp_df <- as.data.frame(
  t(genus_rel_all)
) %>%
  rownames_to_column(
    "OriginalSample"
  ) %>%
  mutate(
    Group = factor(
      group_from_sample(
        OriginalSample
      ),
      levels = analysis_groups
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

comp_order <- comp_df %>%
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
  levels = comp_order$DisplaySample
)

p_S3B <- ggplot(
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
      size = 9,
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    ),
    legend.text = element_text(
      size = 9
    ),
    legend.position = "right"
  )

save_plot(
  p_S3B,
  "FigS3B_Genus_Composition_Top10",
  supp_dir,
  width = 10.0,
  height = 5.0
)

# ======================================================================
# 21. Supplementary direct volcano plots
# ======================================================================

p_S3C <- make_volcano(
  genus_res$Combine_vs_LCD,
  "Genus: Combine vs LCD",
  label_col = "Feature",
  label_n = 8
)

p_S3D <- make_volcano(
  genus_res$Combine_vs_GFC,
  "Genus: Combine vs GFC",
  label_col = "Feature",
  label_n = 8
)

p_S3E <- make_volcano(
  fecal_res$Combine_vs_LCD,
  "Fecal metabolites: Combine vs LCD",
  label_col = "DisplayName",
  label_n = 8,
  annotated_only = TRUE
)

p_S3F <- make_volcano(
  fecal_res$Combine_vs_GFC,
  "Fecal metabolites: Combine vs GFC",
  label_col = "DisplayName",
  label_n = 8,
  annotated_only = TRUE
)

save_plot(
  p_S3C,
  "FigS3C_Genus_Volcano_Combine_vs_LCD",
  supp_dir,
  5.5,
  4.7
)

save_plot(
  p_S3D,
  "FigS3D_Genus_Volcano_Combine_vs_GFC",
  supp_dir,
  5.5,
  4.7
)

save_plot(
  p_S3E,
  "FigS3E_Fecal_Volcano_Combine_vs_LCD",
  supp_dir,
  5.5,
  4.7
)

save_plot(
  p_S3F,
  "FigS3F_Fecal_Volcano_Combine_vs_GFC",
  supp_dir,
  5.5,
  4.7
)

# ======================================================================
# 22. Supplementary S3G | integrated key-feature heatmap
# ======================================================================

matched_samples <- intersect(
  colnames(genus_clr),
  colnames(fecal_mat)
)

sample_order <- make_sample_order_df(
  matched_samples,
  analysis_groups
)

ordered_samples <- sample_order$OriginalSample
display_samples <- sample_order$DisplaySample

genus_heat <- genus_clr[
  key_genus$Feature,
  ordered_samples,
  drop = FALSE
]

metab_heat <- fecal_mat[
  key_metab$Feature,
  ordered_samples,
  drop = FALSE
]

# Display-only wrapping of long metabolite names (analysis names unchanged).
wrap_heat_label <- function(x, width = 42) {
  x <- gsub("_", " ", as.character(x), fixed = TRUE)
  x <- gsub("-", "- ", x, fixed = TRUE)
  stringr::str_wrap(x, width = width)
}
rownames(metab_heat) <- make.unique(vapply(
  key_metab$DisplayName, wrap_heat_label, character(1), width = 35
))

genus_z <- t(
  scale(
    t(genus_heat)
  )
)

metab_z <- t(
  scale(
    t(metab_heat)
  )
)

genus_z[
  !is.finite(genus_z)
] <- 0

metab_z[
  !is.finite(metab_z)
] <- 0

# Genus labels are also display-only wrapped; do not re-transform CLR values.
rownames(genus_z) <- make.unique(vapply(
  rownames(genus_z), wrap_heat_label, character(1), width = 29
))
combined_heat <- rbind(
  genus_z,
  metab_z
)

colnames(combined_heat) <- display_samples

anno_row <- data.frame(
  Omic = factor(
    c(
      rep(
        "Genus",
        nrow(genus_z)
      ),
      rep(
        "Fecal metabolite",
        nrow(metab_z)
      )
    ),
    levels = c(
      "Genus",
      "Fecal metabolite"
    )
  ),
  row.names = rownames(
    combined_heat
  )
)

anno_col <- data.frame(
  Group = factor(
    as.character(
      sample_order$Group
    ),
    levels = analysis_groups
  ),
  row.names = display_samples
)

anno_colors <- list(
  Group = group_cols,
  Omic = c(
    "Genus" = "#4C78A8",
    "Fecal metabolite" = "#59A14F"
  )
)

group_counts <- table(
  factor(
    as.character(
      sample_order$Group
    ),
    levels = analysis_groups
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

gaps_row <- nrow(
  genus_z
)

# Room for sample labels, wrapped metabolite names, and legends at 70% export.
heat_width <- max(15.0, 8.0 + 0.40 * ncol(combined_heat))
heat_height <- max(9.2, 4.0 + 0.40 * nrow(combined_heat))

pdf(
  file.path(FIGURE_DIR,"FigS1I.pdf"),
  width = 0.70 * heat_width,
  height = 0.70 * heat_height
)

pheatmap(
  combined_heat,
  annotation_col = anno_col,
  annotation_row = anno_row,
  annotation_colors = anno_colors,
  annotation_names_row = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  gaps_row = gaps_row,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 9,         # also controls annotation-legend text
  fontsize_row = 8,     # wrapped genus/metabolite names
  fontsize_col = 7,     # many sample labels
  cellheight = 22,
  angle_col = 90,
  border_color = NA,
  main = ""
)

dev.off()

png(
  file.path(FIGURE_DIR,"FigS1I.png"),
  width = 0.60 * heat_width,
  height = 0.60 * heat_height,
  units = "in",
  res = 400
)

pheatmap(
  combined_heat,
  annotation_col = anno_col,
  annotation_row = anno_row,
  annotation_colors = anno_colors,
  annotation_names_row = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  gaps_row = gaps_row,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 8.5,         # also controls annotation-legend text
  fontsize_row = 8.8,     # wrapped genus/metabolite names
  fontsize_col = 8.8,     # many sample labels
  cellheight = 30,
  angle_col = 90,
  border_color = NA,
  main = ""
)

dev.off()

# ======================================================================
# 23. Summary and QC
# ======================================================================

write_csv(
  data.frame(
    RawPrefix = c(
      "MX",
      "DT",
      "GZFL",
      "LH"
    ),
    BiologicalGroup = analysis_groups,
    ExampleRawSample = c(
      "MX_1",
      "DT_1",
      "GZFL_1",
      "LH_1"
    ),
    DisplaySample = c(
      "PCOS_1",
      "LCD_1",
      "GFC_1",
      "Combine_1"
    ),
    GroupOrder = 1:4
  ),
  file.path(
    qc_dir,
    "GLOBAL_sample_naming_and_group_order.csv"
  )
)

summary_lines <- c(

  "============================================================",
  "FIGURE 3 | GLOBAL -> INCREMENT -> MECHANISM",
  "============================================================",
  "",

  "Main scientific question:",
  "  What additional microbiome-metabolome changes does Combine produce",
  "  beyond LCD or GFC alone?",
  "",

  "Main panels:",
  "  Fig3A  Global Genus-level PCoA",
  "  Fig3B  Global fecal-metabolome PCA",
  "  Fig3C1 Genus combination-increment map",
  "  Fig3C2 Metabolite combination-increment map",
  "  Fig3D  Six representative combination-prioritized Genus",
  "  Fig3E1 Six representative combination-prioritized annotated metabolites",
  "  Fig3E2 Literature-informed metabolic-module scores",
  "  Fig3F  Candidate Genus-metabolite correlations",
  "",

  "Global group order:",
  "  PCOS -> LCD -> GFC -> Combine",
  "",

  "Important interpretation rule:",
  "  Combine-vs-PCOS significance alone is NOT considered evidence of",
  "  a combination-specific effect. Direct Combine-vs-LCD and",
  "  Combine-vs-GFC evidence is required.",
  "",

  "Combination-specific:",
  "  Combine is beyond both singles and significant vs PCOS, while",
  "  LCD and GFC are not significant vs PCOS.",
  "",

  "Combination-enhanced:",
  "  LCD/GFC/Combine move in the same direction vs PCOS, but Combine",
  "  has the largest magnitude and differs from both singles.",
  "",

  paste0(
    "Selected key Genus: ",
    paste(
      key_genus$Feature,
      collapse = "; "
    )
  ),
  "",

  paste0(
    "Selected key annotated metabolites: ",
    paste(
      key_metab$DisplayName,
      collapse = "; "
    )
  ),
  "",

  "Metabolic modules are literature-informed candidate modules,",
  "NOT formal KEGG enrichment.",
  "",

  paste0(
    "Fecal metabolome preprocessing: ",
    metab_transform_note
  ),
  "",

  paste0(
    "Matched 16S-metabolome samples: ",
    length(matched_samples)
  ),
  "",

  "============================================================"
)

writeLines(
  summary_lines,
  file.path(
    out_dir,
    "SUMMARY_Figure3.txt"
  )
)

cat(
  "\n============================================================\n"
)
cat(
  "Figure 3 analysis completed.\n"
)
cat(
  "Output directory:\n",
  out_dir,
  "\n"
)
cat(
  "============================================================\n"
)

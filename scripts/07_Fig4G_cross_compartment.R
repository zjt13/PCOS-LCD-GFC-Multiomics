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

FINAL_FIGURE_MAP <- c("FigS5C_FecalSerum_ClockCandidateRelay_Adjusted_v4"="Fig4G")  # Residual scatter: sole Fig4G writer
# Individual plot theme() settings are preserved during export.
# No global font-size overrides are applied in save_plot().

options(stringsAsFactors = FALSE)
set.seed(20261001)

ROOT <- REPO_ROOT
RAW_DIR <- DATA_DIR
MICROBIOME_FILE <- file.path(RAW_DIR, "microbiome_otu_table.tsv")
FECAL_FILE <- file.path(RAW_DIR, "fecal_all_metab_data.tsv")
SERUM_FILE <- file.path(RAW_DIR, "serum_all_metab_data.tsv")
RNA_TPM_FILE <- file.path(RAW_DIR, "rna_gene_tpm_annotation.tsv")
TAXONOMY_FILE <- TAXONOMY_INPUT

FIG3_TABLE_DIR <- file.path(WORK_DIR,"Fig3","Tables")
FIG3_KEY_GENUS <- file.path(FIG3_TABLE_DIR, "MAIN_Key_Genera.csv")
FIG3_KEY_METAB <- file.path(FIG3_TABLE_DIR, "MAIN_Key_Annotated_Metabolites.csv")

OUT_DIR <- file.path(WORK_DIR,"Fig4G")
MAIN_DIR <- file.path(OUT_DIR, "Main")
SUPP_DIR <- file.path(OUT_DIR, "Supplementary")
TABLE_DIR <- file.path(OUT_DIR, "Tables")
QC_DIR <- file.path(OUT_DIR, "QC")

for (d in c(OUT_DIR, MAIN_DIR, SUPP_DIR, TABLE_DIR, QC_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cran_pkgs <- c("ggplot2","dplyr","tidyr","readr","stringr","ggrepel","pheatmap","tibble","scales")
miss <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly=TRUE)]
if(length(miss)>0) install.packages(miss, repos="https://cloud.r-project.org")
if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
if(!requireNamespace("limma", quietly=TRUE)) BiocManager::install("limma", ask=FALSE, update=FALSE)

# Attach dplyr to provide %>% while keeping critical verbs explicitly namespaced.
library(dplyr)

GROUP_LEVELS <- c("PCOS","LCD","GFC","Combine")
GROUP_COLS <- c("PCOS"="#B18459","LCD"="#8E44AD","GFC"="#F28E2B","Combine"="#E76F51")
OMICS_COLS <- c("Genus"="#4C78A8","Fecal metabolite"="#59A14F","Plasma metabolite"="#ECA82C","Ovarian clock"="#B279A2")
COR_LOW <- "#3B73B9"; COR_HIGH <- "#C94C4C"

CLOCK_ALIASES <- list(
  Bmal1=c("Arntl","Bmal1"), Npas2="Npas2", Per2="Per2", Per3="Per3",
  Dbp="Dbp", Nr1d2=c("Nr1d2","Rev-erbb"), Clock="Clock",
  Cry1="Cry1", Cry2="Cry2", Tef="Tef", Hlf="Hlf"
)
CLOCK_DISPLAY_ORDER <- c("Bmal1","Npas2","Per2","Per3","Dbp","Nr1d2","Clock","Cry1","Cry2","Tef","Hlf")
CLOCK_MODULE_PRIMARY <- c("Bmal1","Npas2","Per2","Per3","Dbp","Nr1d2")
LITERATURE_TAXA <- c("Alistipes","Odoribacter","Mucispirillum","Candidatus_Arthromitus",
                     "[Eubacterium]_coprostanoligenes_group",
                     "[Eubacterium]_coprostanoligenes_group_unclassified")

literature_evidence <- data.frame(
  EvidenceType=c("Taxon","Taxon","Metabolite class","Metabolite class","Metabolite class","Metabolite","PCOS clock"),
  Feature=c("Alistipes / Odoribacter","Gut microbiota rhythmicity","Unconjugated bile acids",
            "Short-chain fatty acids","Phenylalanine-derived microbial metabolites","2-Phenylacetamide","BMAL1/NPAS2/PER axis"),
  ClockLink=c(
    "Differential rhythmicity/abundance after intestinal Bmal1 deletion",
    "Arrhythmic microbiota can alter peripheral clock-gene expression",
    "DCA/CDCA altered Dbp, Per2, Per3 and Cry2 expression",
    "SCFAs can phase-shift/entrain PER2 and BMAL1 through HDAC-linked mechanisms",
    "3-phenylpropionic acid and 4-OH-PPA altered PER2/Bmal1 oscillation",
    "Changed in intestinal Bmal1-deficient metabolomics; association direction is clock-to-metabolite",
    "PCOS associated with abnormal BMAL1, NPAS2, PER and steroidogenic programs"),
  Reference=c(
    "Nature Communications 2022;13:6068; DOI:10.1038/s41467-022-33609-x",
    "PMID:36334897","PMID:27907092","PMID:35934064; PMID:40727446",
    "PMID:32155689","PMCID:PMC12442615","PMID:35728080; PMID:32119846"),
  stringsAsFactors=FALSE
)
readr::write_csv(literature_evidence, file.path(TABLE_DIR,"Figure5_Literature_Evidence.csv"))

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

group_from_sample <- function(x){
  prefix <- sub("_.*$","",x)
  dplyr::recode(prefix, "MX"="PCOS","DT"="LCD","GZFL"="GFC","LH"="Combine","QC"="QC", .default=NA_character_)
}
is_biological_sample <- function(x) grepl("^(MX|DT|GZFL|LH)_[0-9]+$",x)
sample_number <- function(x) suppressWarnings(as.numeric(sub("^.*_","",x)))
sample_display <- function(x){
  g <- group_from_sample(x); n <- sample_number(x)
  ifelse(!is.na(g)&!is.na(n), paste0(g,"_",n), x)
}
safe_numeric_matrix <- function(df,sample_cols,row_ids){
  m <- as.matrix(df[,sample_cols,drop=FALSE]); suppressWarnings(storage.mode(m)<-"numeric")
  rownames(m) <- make.unique(as.character(row_ids)); m[!is.finite(m)] <- NA_real_; m
}
zscore_rows <- function(mat){ z <- t(scale(t(mat))); z[!is.finite(z)] <- 0; z }
find_first_col <- function(nms,candidates){
  low <- tolower(nms)
  for(ca in candidates){ ii <- which(low==tolower(ca)); if(length(ii)>0) return(nms[ii[1]]) }
  NULL
}
normalize_metabolite_name <- function(x){
  y <- tolower(trimws(as.character(x))); y <- gsub("\\[[^\\]]+\\]","",y); y <- gsub("\\([^\\)]*\\)","",y)
  gsub("[^a-z0-9]+","",y)
}
is_valid_metabolite_name <- function(x){
  y <- trimws(as.character(x))
  !is.na(y)&y!=""&!grepl("^metab_[0-9]+$",y,ignore.case=TRUE)&
    !grepl("^[0-9]+(\\.[0-9]+)?_[0-9]+(\\.[0-9]+)?m/z$",y,ignore.case=TRUE)&
    !tolower(y)%in%c("unknown","unidentified","na","none","-")
}
classify_metabolite <- function(x){
  y <- tolower(as.character(x))
  dplyr::case_when(
    grepl("cholic acid|cholate|deoxychol|chenodeoxychol|lithochol|ursodeoxychol|murichol|hyodeoxychol|taurochol|glycochol|bile acid|bile salt",y) ~ "Bile acid",
    grepl("phenylacetamide|phenylacetic|phenylacetyl|phenyllact|phenylprop|hydroxyphenylprop|phenylalanine",y) ~ "Phenylalanine-derived",
    grepl("acetate|acetic acid|propionate|propionic acid|butyrate|butyric acid|isobutyr|isovaler|valerate|2-methylbutyr",y) ~ "SCFA / BCFA",
    grepl("indole|tryptophan|kynuren|tryptamine",y) ~ "Tryptophan / indole",
    grepl("spermidine|spermine|putrescine|cadaverine",y) ~ "Polyamine",
    TRUE ~ "Other"
  )
}
residualize_by_group <- function(x,group){
  ok <- is.finite(x)&!is.na(group); out <- rep(NA_real_,length(x))
  if(sum(ok)<6 || length(unique(group[ok]))<2) return(out)
  out[ok] <- stats::residuals(stats::lm(x[ok]~group[ok])); out
}
adjusted_spearman <- function(x,y,group){
  ok <- is.finite(x)&is.finite(y)&!is.na(group); n <- sum(ok)
  if(n<6) return(data.frame(N=n,Rho=NA_real_,P=NA_real_))
  xr <- residualize_by_group(x[ok],group[ok]); yr <- residualize_by_group(y[ok],group[ok])
  ok2 <- is.finite(xr)&is.finite(yr)
  if(sum(ok2)<6) return(data.frame(N=sum(ok2),Rho=NA_real_,P=NA_real_))
  ct <- suppressWarnings(stats::cor.test(xr[ok2],yr[ok2],method="spearman",exact=FALSE))
  data.frame(N=sum(ok2),Rho=unname(ct$estimate),P=ct$p.value)
}
transform_metabolome <- function(mat){
  v <- as.numeric(mat[is.finite(mat)]); if(length(v)==0) return(mat)
  q99 <- stats::quantile(v,.99,na.rm=TRUE); med <- stats::median(v,na.rm=TRUE)
  if(q99>50 || med>20){
    mn <- min(v,na.rm=TRUE); sh <- ifelse(mn<=0,abs(mn)+1,1); return(log2(mat+sh))
  }
  mat
}

# ---------- 16S to Genus CLR ----------
if(!file.exists(MICROBIOME_FILE)||!file.exists(TAXONOMY_FILE)) stop("Missing microbiome/taxonomy input")
micro_df <- readr::read_tsv(MICROBIOME_FILE,show_col_types=FALSE,progress=FALSE,name_repair="minimal")
names(micro_df)[1] <- "OTU_ID"
micro_samples <- names(micro_df)[vapply(names(micro_df),is_biological_sample,logical(1))]
micro_mat <- safe_numeric_matrix(micro_df,micro_samples,micro_df$OTU_ID); micro_mat[is.na(micro_mat)] <- 0

tax_raw <- utils::read.delim(TAXONOMY_FILE,header=FALSE,sep="\t",quote="",comment.char="",fill=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
extract_rank <- function(tax,prefix){ mm <- stringr::str_match(tax,paste0("(?:^|;)",prefix,"([^;]+)")); mm[,2] }
clean_taxon <- function(genus_raw,family_raw,order_raw){
  cl <- function(x){ x <- ifelse(is.na(x),"",x); x <- gsub("^[gfo]__","",x); trimws(x) }
  genus <- cl(genus_raw); fam <- cl(family_raw); ord <- cl(order_raw); out <- genus
  idx <- grepl("^norank_f__",out); out[idx] <- paste0(sub("^norank_f__","",out[idx]),"_unclassified")
  idx <- grepl("^norank_o__",out); out[idx] <- paste0(sub("^norank_o__","",out[idx]),"_unclassified")
  generic <- out==""|is.na(out)|tolower(out)%in%c("norank","unclassified","uncultured","unknown")
  use <- generic & fam!=""; out[use] <- paste0(fam[use],"_unclassified")
  generic <- out==""|is.na(out)|tolower(out)%in%c("norank","unclassified","uncultured","unknown")
  use <- generic & ord!=""; out[use] <- paste0(ord[use],"_unclassified")
  out[out==""|is.na(out)] <- "Unclassified"; out
}
tax_df <- data.frame(OTU_ID=trimws(as.character(tax_raw[[1]])),Taxonomy=trimws(as.character(tax_raw[[2]])),stringsAsFactors=FALSE) %>%
  dplyr::mutate(Order=extract_rank(Taxonomy,"o__"),Family=extract_rank(Taxonomy,"f__"),Genus_raw=extract_rank(Taxonomy,"g__"),
                Genus=clean_taxon(Genus_raw,Family,Order))
idx <- match(rownames(micro_mat),tax_df$OTU_ID); genus_names <- tax_df$Genus[idx]; genus_names[is.na(genus_names)|genus_names==""] <- "Unclassified"
tmp <- as.data.frame(micro_mat); tmp$Genus <- genus_names
gdf <- tmp %>% dplyr::group_by(Genus) %>% dplyr::summarise(dplyr::across(dplyr::all_of(micro_samples),\(x) sum(x,na.rm=TRUE)),.groups="drop")
genus_count_mat <- as.matrix(gdf[,micro_samples,drop=FALSE]); rownames(genus_count_mat) <- gdf$Genus
genus_clr <- log(genus_count_mat+0.5); genus_clr <- sweep(genus_clr,2,colMeans(genus_clr),FUN="-")

# ---------- metabolomics ----------
read_metabolome <- function(path,omic_name){
  if(!file.exists(path)) stop("Missing ",omic_name,": ",path)
  df <- readr::read_tsv(path,show_col_types=FALSE,progress=FALSE,name_repair="minimal")
  id_col <- find_first_col(names(df),c("metab_id","metabolite_id","feature_id","id")); if(is.null(id_col)) id_col <- names(df)[1]
  name_col <- find_first_col(names(df),c("Metabolite","metabolite","Compound","compound","Name","name"))
  samples <- names(df)[vapply(names(df),is_biological_sample,logical(1))]
  ids <- as.character(df[[id_col]]); mat <- transform_metabolome(safe_numeric_matrix(df,samples,ids))
  nm <- if(is.null(name_col)) ids else as.character(df[[name_col]])
  valid <- is_valid_metabolite_name(nm); nm[!valid] <- ids[!valid]
  annot <- data.frame(FeatureID=ids,DisplayName=nm,IsAnnotated=valid,CanonicalName=normalize_metabolite_name(nm),
                      MetaboliteClass=classify_metabolite(nm),stringsAsFactors=FALSE)
  list(mat=mat,annot=annot,samples=samples)
}
fecal <- read_metabolome(FECAL_FILE,"fecal metabolome")
serum <- read_metabolome(SERUM_FILE,"serum metabolome")
readr::write_csv(fecal$annot,file.path(TABLE_DIR,"Fecal_Metabolite_NameMap.csv"))
readr::write_csv(serum$annot,file.path(TABLE_DIR,"Serum_Metabolite_NameMap.csv"))

collapse_metabolites <- function(mat,annot){
  ids <- intersect(rownames(mat),annot$FeatureID); mm <- mat[ids,,drop=FALSE]; aa <- annot[match(ids,annot$FeatureID),,drop=FALSE]
  keep <- aa$IsAnnotated & aa$CanonicalName!=""; mm <- mm[keep,,drop=FALSE]; aa <- aa[keep,,drop=FALSE]
  zz <- zscore_rows(mm); td <- as.data.frame(zz); td$CanonicalName <- aa$CanonicalName; td$DisplayName <- aa$DisplayName; td$MetaboliteClass <- aa$MetaboliteClass
  co <- td %>% dplyr::group_by(CanonicalName) %>% dplyr::summarise(
    dplyr::across(dplyr::all_of(colnames(zz)),\(x) mean(x,na.rm=TRUE)),
    DisplayName=dplyr::first(DisplayName),MetaboliteClass=dplyr::first(MetaboliteClass),.groups="drop")
  out <- as.matrix(co[,colnames(zz),drop=FALSE]); rownames(out) <- co$CanonicalName
  list(mat=out,map=co %>% dplyr::select(CanonicalName,DisplayName,MetaboliteClass))
}
fecal_c <- collapse_metabolites(fecal$mat,fecal$annot)
serum_c <- collapse_metabolites(serum$mat,serum$annot)

# ---------- ovarian clock ----------
if(!file.exists(RNA_TPM_FILE)) stop("Missing RNA TPM input")
rna <- readr::read_tsv(RNA_TPM_FILE,show_col_types=FALSE,progress=FALSE,name_repair="minimal")
names(rna)[1:4] <- c("Gene_id","Gene_name","Gene_description","Bio_type")
rna_samples <- names(rna)[vapply(names(rna),is_biological_sample,logical(1))]
tpm <- safe_numeric_matrix(rna,rna_samples,rna$Gene_id); tpm[is.na(tpm)] <- 0; rlog <- log2(tpm+1)
gn <- as.character(rna$Gene_name[match(rownames(rlog),rna$Gene_id)]); gn[is.na(gn)|gn==""] <- rownames(rlog)[is.na(gn)|gn==""]
rd <- as.data.frame(rlog); rd$Gene_name <- gn
rd <- rd %>% dplyr::group_by(Gene_name) %>% dplyr::summarise(dplyr::across(dplyr::all_of(rna_samples),\(x) mean(x,na.rm=TRUE)),.groups="drop")
rmat <- as.matrix(rd[,rna_samples,drop=FALSE]); rownames(rmat) <- rd$Gene_name
resolve_clock_gene <- function(alias_vec,mat){ hit <- alias_vec[alias_vec%in%rownames(mat)]; if(length(hit)==0) NA_character_ else hit[1] }
clock_map <- data.frame(Display=names(CLOCK_ALIASES),GeneSymbol=vapply(CLOCK_ALIASES,resolve_clock_gene,character(1),mat=rmat),stringsAsFactors=FALSE) %>%
  dplyr::filter(!is.na(GeneSymbol))
readr::write_csv(clock_map,file.path(QC_DIR,"ClockGene_Mapping.csv"))
clock_expr <- rmat[clock_map$GeneSymbol,,drop=FALSE]; rownames(clock_expr) <- clock_map$Display; clock_z <- zscore_rows(clock_expr)
primary_detected <- intersect(CLOCK_MODULE_PRIMARY,rownames(clock_z)); if(length(primary_detected)<3) stop("Fewer than 3 primary clock genes detected")
clock_score <- colMeans(clock_z[primary_detected,,drop=FALSE],na.rm=TRUE)

# ---------- Fig5A1 ----------
meta <- data.frame(Sample=colnames(clock_z),Group=factor(group_from_sample(colnames(clock_z)),levels=GROUP_LEVELS),
                   Num=sample_number(colnames(clock_z)),stringsAsFactors=FALSE) %>% dplyr::arrange(Group,Num)
hm <- clock_z[intersect(CLOCK_DISPLAY_ORDER,rownames(clock_z)),meta$Sample,drop=FALSE]; colnames(hm) <- sample_display(meta$Sample)
ann <- data.frame(Group=meta$Group,row.names=colnames(hm))
gaps <- cumsum(as.integer(table(factor(meta$Group,levels=GROUP_LEVELS)))); gaps <- gaps[-length(gaps)]
for(ext in c("pdf","png")){
  fn <- file.path(tempdir(),paste0("Fig5A1_OvarianClockGene_Heatmap.",ext))
  if(ext=="pdf") grDevices::pdf(fn,width=8.7,height=4.9) else grDevices::png(fn,width=8.7,height=4.9,units="in",res=400)
  pheatmap::pheatmap(hm,annotation_col=ann,annotation_colors=list(Group=GROUP_COLS),cluster_rows=FALSE,cluster_cols=FALSE,
                     gaps_col=gaps,angle_col=45,fontsize=11,fontsize_row=11,fontsize_col=9,border_color=NA,main="")
  grDevices::dev.off()
}

# ---------- Fig5A2 ----------
clock_score_df <- data.frame(Sample=names(clock_score),ClockScore=as.numeric(clock_score),
                             Group=factor(group_from_sample(names(clock_score)),levels=GROUP_LEVELS),stringsAsFactors=FALSE)
readr::write_csv(clock_score_df,file.path(TABLE_DIR,"Ovarian_ClockModule_Scores.csv"))
pA2 <- ggplot2::ggplot(clock_score_df,ggplot2::aes(Group,ClockScore,fill=Group))+
  ggplot2::geom_boxplot(width=.58,outlier.shape=NA,alpha=.28,colour="grey35")+
  ggplot2::geom_jitter(ggplot2::aes(colour=Group),width=.08,size=2.4,alpha=.92)+
  ggplot2::scale_fill_manual(values=GROUP_COLS,guide="none")+ggplot2::scale_colour_manual(values=GROUP_COLS,guide="none")+
  ggplot2::labs(x=NULL,y="Ovarian clock-module score\n(mean z-score)")+
  ggplot2::theme_classic(base_size=14)+ggplot2::theme(axis.text=ggplot2::element_text(colour="black"),
    axis.text.x=ggplot2::element_text(angle=35,hjust=1),axis.title.y=ggplot2::element_text(face="bold"))
save_plot(pA2,"Fig5A2_OvarianClockModule_Score",MAIN_DIR,4.8,4.4)

# ---------- generic limma ----------
run_limma <- function(mat,omic){
  sam <- colnames(mat); grp <- factor(group_from_sample(sam),levels=GROUP_LEVELS); keep <- !is.na(grp); mat <- mat[,keep,drop=FALSE]; grp <- droplevels(grp[keep])
  design <- stats::model.matrix(~0+grp); colnames(design) <- levels(grp); fit <- limma::lmFit(mat,design)
  cont <- limma::makeContrasts(LCD_vs_PCOS=LCD-PCOS,GFC_vs_PCOS=GFC-PCOS,Combine_vs_PCOS=Combine-PCOS,
                               Combine_vs_LCD=Combine-LCD,Combine_vs_GFC=Combine-GFC,levels=design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit,cont))
  dplyr::bind_rows(lapply(colnames(cont),function(nm){
    tt <- limma::topTable(fit2,coef=nm,number=Inf,sort.by="none",adjust.method="BH"); tt$Feature <- rownames(tt); tt$Contrast <- nm; tt$Omic <- omic; tt
  }))
}
fecal_diff <- run_limma(fecal_c$mat,"Fecal metabolite"); serum_diff <- run_limma(serum_c$mat,"Serum metabolite")
genus_diff <- run_limma(genus_clr,"Genus")
metab_diff <- dplyr::bind_rows(fecal_diff,serum_diff)
readr::write_csv(metab_diff,file.path(TABLE_DIR,"Metabolite_Limma_AllContrasts.csv"))

# ---------- prior Figure3 selections ----------
prior_genera <- character(0)
prior_metabs <- character(0)

if(file.exists(FIG3_KEY_GENUS)){
  x <- tryCatch(
    readr::read_csv(FIG3_KEY_GENUS,show_col_types=FALSE),
    error=function(e) NULL
  )
  if(!is.null(x)){
    cc <- intersect(c("Feature","Genus"),names(x))
    if(length(cc)>0) prior_genera <- unique(as.character(x[[cc[1]]]))
  }
}

if(file.exists(FIG3_KEY_METAB)){
  x <- tryCatch(
    readr::read_csv(FIG3_KEY_METAB,show_col_types=FALSE),
    error=function(e) NULL
  )
  if(!is.null(x)){
    cc <- intersect(c("DisplayName","Metabolite","Feature"),names(x))
    if(length(cc)>0) prior_metabs <- unique(as.character(x[[cc[1]]]))
  }
}

prior_metab_canon <- unique(
  normalize_metabolite_name(prior_metabs)
)

# ----------------------------------------------------------------------
# Explicit metabolite-clock evidence tiers
#
# Score 3 = direct experimental evidence that the metabolite can alter
#           clock-gene expression/oscillation.
# Score 2 = metabolite itself shows circadian cycling / is mechanistically
#           connected to clock-active pathways (e.g. AHR), but direct
#           ovarian clock regulation is not established.
# Score 1 = pathway-level or reverse-direction evidence only.
# Score 0 = no sufficiently specific clock evidence used for main figures.
#
# IMPORTANT:
# A broad class such as "bile acid" or "tryptophan metabolite" is NOT
# automatically treated as a direct clock regulator.
# ----------------------------------------------------------------------

clock_metabolite_evidence <- function(x){

  y <- tolower(trimws(as.character(x)))

  # Do not let Nordeoxycholic acid match Deoxycholic acid via a substring.
  unconjugated_dca_cdca <-
    y %in% c("deoxycholic acid", "chenodeoxycholic acid")

  direct_aromatic <-
    grepl("3[- ]?phenylpropionic acid|phenylpropanoic acid|4[- ]?hydroxyphenylpropionic acid|3[- ]?\\(4[- ]?hydroxyphenyl\\)propionic acid", y)

  direct_scfa <-
    grepl("^(acetate|acetic acid|propionate|propionic acid|butyrate|butyric acid)$", y)

  direct_urolithin <-
    grepl("^urolithin a$", y)

  circadian_trp <-
    grepl("kynurenic acid|^kynurenine$|indole[- ]?3[- ]?acet|indole[- ]?3[- ]?lact|indole[- ]?3[- ]?prop|indole[- ]?3[- ]?aldehyde|indole aldehyde", y)

  reverse_clock_assoc <-
    grepl("2[- ]?phenylacetamide", y)

  pathway_only <-
    grepl("cholic acid|cholate|5[- ]?hydroxyindole[- ]?3[- ]?acetic acid|5[- ]?hiaa|^indole$|n[- ]?acetyl[- ]?l[- ]?phenylalanine|methoxyphenylacetic", y) |
    y %in% c("nordeoxycholic acid", "nord(e)oxycholic acid")

  score <- dplyr::case_when(
    unconjugated_dca_cdca | direct_aromatic | direct_scfa | direct_urolithin ~ 3,
    circadian_trp ~ 2,
    reverse_clock_assoc | pathway_only ~ 1,
    TRUE ~ 0
  )

  level <- dplyr::case_when(
    score == 3 ~ "Direct clock evidence",
    score == 2 ~ "Circadian-associated",
    score == 1 ~ "Pathway/reverse-direction only",
    TRUE ~ "No specific clock evidence"
  )

  ref <- dplyr::case_when(
    unconjugated_dca_cdca ~ "PMID:27907092",
    direct_aromatic ~ "PMID:32155689",
    direct_scfa ~ "PMID:35934064",
    direct_urolithin ~ "PMID:39064706",
    circadian_trp ~ "PMID:41727182; AHR-clock crosstalk PMID:30706546",
    reverse_clock_assoc ~ "PMCID:PMC12442615 (clock-to-metabolite direction)",
    pathway_only ~ "Pathway-level evidence; no direct clock-regulator claim",
    TRUE ~ ""
  )

  data.frame(
    ClockEvidenceScore = score,
    ClockEvidenceLevel = level,
    ClockEvidenceReference = ref,
    stringsAsFactors = FALSE
  )
}

# ----------------------------------------------------------------------
# Specific genus-metabolite literature evidence
# Score 3 = direct production / experimental link
# Score 2 = published genus-metabolite correlation
# Score 1 = broader pathway association
# Score 0 = no specific pair evidence identified
# ----------------------------------------------------------------------

genus_metabolite_pair_evidence <- function(genus, metabolite){

  g <- tolower(as.character(genus))
  m <- tolower(as.character(metabolite))

  alistipes_5hiaa <-
    grepl("alistipes",g) &
    grepl("5[- ]?hydroxyindole[- ]?3[- ]?acetic acid|5[- ]?hiaa",m)

  prev_kyna <-
    grepl("prevotellaceae.*ucg[-_ ]?001",g) &
    grepl("kynurenic acid",m)

  alistipes_kyna <-
    grepl("alistipes",g) &
    grepl("kynurenic acid",m)

  odoribacter_indole <-
    grepl("odoribacter",g) &
    grepl("indole[- ]?3[- ]?propionic acid|tryptamine",m)

  score <- dplyr::case_when(
    alistipes_5hiaa ~ 3,
    prev_kyna ~ 2,
    odoribacter_indole ~ 2,
    alistipes_kyna ~ 1,
    TRUE ~ 0
  )

  evidence <- dplyr::case_when(
    alistipes_5hiaa ~ "Alistipes onderdonkii produced 5-HIAA in vitro",
    prev_kyna ~ "Published positive correlation between Prevotellaceae UCG-001 and kynurenic acid",
    odoribacter_indole ~ "Odoribacter/Odoribacteraceae associated with circulating tryptophan-indole metabolites",
    alistipes_kyna ~ "Alistipes-associated intervention increased kynurenic acid; direct KYNA production not established",
    TRUE ~ "No specific pair evidence identified"
  )

  ref <- dplyr::case_when(
    alistipes_5hiaa ~ "PMID:41692705",
    prev_kyna ~ "PMCID:PMC9764441",
    odoribacter_indole ~ "PMCID:PMC12606010",
    alistipes_kyna ~ "PMID:41692705",
    TRUE ~ ""
  )

  data.frame(
    PairEvidenceScore = score,
    PairEvidence = evidence,
    PairEvidenceReference = ref,
    stringsAsFactors = FALSE
  )
}

candidate_meta <- dplyr::bind_rows(
  fecal_c$map %>% dplyr::mutate(Omic="Fecal metabolite"),
  serum_c$map %>% dplyr::mutate(Omic="Serum metabolite")
)

ev <- clock_metabolite_evidence(
  candidate_meta$DisplayName
)

candidate_meta <- dplyr::bind_cols(
  candidate_meta,
  ev
) %>%
  dplyr::mutate(
    PriorFig3 = CanonicalName %in% prior_metab_canon,
    LiteraturePriority = ClockEvidenceScore >= 2
  ) %>%
  dplyr::left_join(
    metab_diff %>%
      dplyr::group_by(Omic,Feature) %>%
      dplyr::summarise(
        MinFDR=min(adj.P.Val,na.rm=TRUE),
        MaxAbsEffect=max(abs(logFC),na.rm=TRUE),
        .groups="drop"
      ),
    by=c("Omic","CanonicalName"="Feature")
  ) %>%
  dplyr::mutate(
    MinFDR=ifelse(is.finite(MinFDR),MinFDR,1),
    MaxAbsEffect=ifelse(is.finite(MaxAbsEffect),MaxAbsEffect,0),
    CandidateScore =
      2.5*ClockEvidenceScore +
      1.5*as.numeric(PriorFig3) +
      -log10(pmax(MinFDR,1e-6)) +
      pmin(MaxAbsEffect,3)
  )

readr::write_csv(
  candidate_meta,
  file.path(TABLE_DIR,"Candidate_Metabolite_Prioritization_v3.csv")
)

genus_rank <- genus_diff %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(
    MinFDR=min(adj.P.Val,na.rm=TRUE),
    MaxAbsEffect=max(abs(logFC),na.rm=TRUE),
    .groups="drop"
  ) %>%
  dplyr::arrange(MinFDR,dplyr::desc(MaxAbsEffect))

candidate_genera <- unique(
  c(
    prior_genera,
    LITERATURE_TAXA,
    head(genus_rank$Feature,12)
  )
)

candidate_genera <- candidate_genera[
  candidate_genera %in% rownames(genus_clr)
]

candidate_genera <- candidate_genera[
  !grepl("^Unclassified$",candidate_genera)
]

candidate_genera <- head(candidate_genera,16)

readr::write_csv(
  data.frame(
    Genus=candidate_genera,
    LiteraturePriority=candidate_genera %in% LITERATURE_TAXA,
    PreviousFig3=candidate_genera %in% prior_genera
  ),
  file.path(TABLE_DIR,"Candidate_Genera.csv")
)

# ---------- sample matching ----------
sets <- list(Genus=colnames(genus_clr),Fecal=colnames(fecal_c$mat),Serum=colnames(serum_c$mat),RNA=colnames(clock_z))
matched_all4 <- Reduce(intersect,sets)
mq <- dplyr::bind_rows(lapply(GROUP_LEVELS,function(g)data.frame(
  Group=g,N_Genus=sum(group_from_sample(sets$Genus)==g),N_Fecal=sum(group_from_sample(sets$Fecal)==g),
  N_Serum=sum(group_from_sample(sets$Serum)==g),N_RNA=sum(group_from_sample(sets$RNA)==g),
  N_All4Matched=sum(group_from_sample(matched_all4)==g))))
readr::write_csv(mq,file.path(QC_DIR,"CrossOmics_SampleMatching.csv"))

# ---------- shared metabolites and correlations ----------
shared <- intersect(rownames(fecal_c$mat),rownames(serum_c$mat))
shared_map <- fecal_c$map %>% dplyr::filter(CanonicalName%in%shared) %>%
  dplyr::select(CanonicalName,FecalName=DisplayName,FecalClass=MetaboliteClass) %>%
  dplyr::left_join(serum_c$map %>% dplyr::filter(CanonicalName%in%shared) %>%
    dplyr::select(CanonicalName,SerumName=DisplayName,SerumClass=MetaboliteClass),by="CanonicalName") %>%
  dplyr::mutate(MetaboliteClass=ifelse(FecalClass!="Other",FecalClass,SerumClass),PriorFig3=CanonicalName%in%prior_metab_canon)

fs <- dplyr::bind_rows(lapply(shared,function(m){
  sam <- intersect(colnames(fecal_c$mat),colnames(serum_c$mat)); grp <- factor(group_from_sample(sam),levels=GROUP_LEVELS)
  cc <- adjusted_spearman(fecal_c$mat[m,sam],serum_c$mat[m,sam],grp)
  data.frame(CanonicalName=m,N=cc$N,Rho=cc$Rho,P=cc$P)
})) %>% dplyr::left_join(shared_map,by="CanonicalName") %>%
  dplyr::mutate(FDR=stats::p.adjust(P,"BH"),LiteraturePriority=MetaboliteClass!="Other")
readr::write_csv(fs,file.path(TABLE_DIR,"Fecal_Serum_SharedMetabolite_Correlations.csv"))

cand_fecal <- candidate_meta %>% dplyr::filter(Omic=="Fecal metabolite") %>% dplyr::arrange(dplyr::desc(CandidateScore)) %>% dplyr::pull(CanonicalName) %>% unique()
cand_fecal <- head(intersect(cand_fecal,rownames(fecal_c$mat)),40)
cand_serum <- candidate_meta %>% dplyr::filter(Omic=="Serum metabolite") %>% dplyr::arrange(dplyr::desc(CandidateScore)) %>% dplyr::pull(CanonicalName) %>% unique()
cand_serum <- head(intersect(cand_serum,rownames(serum_c$mat)),50)

gf_list <- list(); k <- 1
for(g in candidate_genera) for(m in cand_fecal){
  sam <- intersect(colnames(genus_clr),colnames(fecal_c$mat)); grp <- factor(group_from_sample(sam),levels=GROUP_LEVELS)
  cc <- adjusted_spearman(genus_clr[g,sam],fecal_c$mat[m,sam],grp)
  gf_list[[k]] <- data.frame(Genus=g,CanonicalName=m,N=cc$N,Rho=cc$Rho,P=cc$P); k <- k+1
}
gf <- dplyr::bind_rows(gf_list) %>% dplyr::mutate(FDR=stats::p.adjust(P,"BH")) %>% dplyr::left_join(fecal_c$map,by="CanonicalName")
readr::write_csv(gf,file.path(TABLE_DIR,"Genus_FecalMetabolite_Correlations.csv"))

sc_list <- list(); k <- 1
for(m in cand_serum){
  sam <- intersect(colnames(serum_c$mat),names(clock_score)); grp <- factor(group_from_sample(sam),levels=GROUP_LEVELS)
  cc <- adjusted_spearman(serum_c$mat[m,sam],clock_score[sam],grp)
  sc_list[[k]] <- data.frame(CanonicalName=m,ClockFeature="Clock module",N=cc$N,Rho=cc$Rho,P=cc$P); k <- k+1
  for(cg in rownames(clock_z)){
    cc <- adjusted_spearman(serum_c$mat[m,sam],clock_z[cg,sam],grp)
    sc_list[[k]] <- data.frame(CanonicalName=m,ClockFeature=cg,N=cc$N,Rho=cc$Rho,P=cc$P); k <- k+1
  }
}
sc <- dplyr::bind_rows(sc_list) %>% dplyr::group_by(ClockFeature) %>% dplyr::mutate(FDR=stats::p.adjust(P,"BH")) %>%
  dplyr::ungroup() %>% dplyr::left_join(serum_c$map,by="CanonicalName")
readr::write_csv(sc,file.path(TABLE_DIR,"SerumMetabolite_OvarianClock_Correlations.csv"))

gs_list <- list(); k <- 1
for(g in candidate_genera) for(m in cand_serum){
  sam <- intersect(colnames(genus_clr),colnames(serum_c$mat)); grp <- factor(group_from_sample(sam),levels=GROUP_LEVELS)
  cc <- adjusted_spearman(genus_clr[g,sam],serum_c$mat[m,sam],grp)
  gs_list[[k]] <- data.frame(Genus=g,CanonicalName=m,N=cc$N,Rho=cc$Rho,P=cc$P); k <- k+1
}
gs <- dplyr::bind_rows(gs_list) %>% dplyr::mutate(FDR=stats::p.adjust(P,"BH")) %>% dplyr::left_join(serum_c$map,by="CanonicalName")
readr::write_csv(gs,file.path(TABLE_DIR,"Genus_SerumMetabolite_Correlations.csv"))

# ---------- 4-layer axis ----------
# Genus -> fecal metabolite -> same serum metabolite -> ovarian clock module
#
# Two outputs:
# 1) Candidate_4Layer_Axes_DataDriven.csv = all data-driven candidates
# 2) Candidate_4Layer_Axes_MainEligible.csv = axes eligible for the MAIN network
#
# A main-network axis must not be selected solely because its correlations
# are large. It must have explicit metabolite-clock evidence and either
# specific genus-metabolite literature support or unusually strong data.
# ----------------------------------------------------------------------

axis_rows <- list()
k <- 1

shared_candidate <- shared_map %>%
  dplyr::filter(
    CanonicalName %in% cand_fecal,
    CanonicalName %in% cand_serum
  )

for(i in seq_len(nrow(shared_candidate))){

  m <- shared_candidate$CanonicalName[i]
  m_name <- shared_candidate$FecalName[i]

  fs1 <- fs %>%
    dplyr::filter(CanonicalName==m) %>%
    dplyr::slice_head(n=1)

  sc1 <- sc %>%
    dplyr::filter(
      CanonicalName==m,
      ClockFeature=="Clock module"
    ) %>%
    dplyr::slice_head(n=1)

  if(nrow(fs1)==0 || nrow(sc1)==0) next

  mev <- clock_metabolite_evidence(m_name)

  for(g in candidate_genera){

    gf1 <- gf %>%
      dplyr::filter(
        Genus==g,
        CanonicalName==m
      ) %>%
      dplyr::slice_head(n=1)

    if(nrow(gf1)==0) next

    pev <- genus_metabolite_pair_evidence(
      g,
      m_name
    )

    litg <- g %in% LITERATURE_TAXA
    priorg <- g %in% prior_genera

    n_edges_p10 <- sum(
      c(
        gf1$P < .10,
        fs1$P < .10,
        sc1$P < .10
      ),
      na.rm=TRUE
    )

    n_edges_p05 <- sum(
      c(
        gf1$P < .05,
        fs1$P < .05,
        sc1$P < .05
      ),
      na.rm=TRUE
    )

    all_edges_p10 <-
      !is.na(gf1$P) & gf1$P < .10 &
      !is.na(fs1$P) & fs1$P < .10 &
      !is.na(sc1$P) & sc1$P < .10

    raw_score <-
      abs(gf1$Rho) *
      abs(fs1$Rho) *
      abs(sc1$Rho)

    axis_score <-
      raw_score +
      0.20*mev$ClockEvidenceScore +
      0.15*pev$PairEvidenceScore +
      0.08*as.numeric(litg) +
      0.08*as.numeric(priorg) +
      0.06*as.numeric(shared_candidate$PriorFig3[i])

    main_tier <- dplyr::case_when(

      # Strongest: direct clock-active metabolite + all 3 edges nominally supported
      mev$ClockEvidenceScore == 3 &
        all_edges_p10 &
        (pev$PairEvidenceScore >= 1 | litg | priorg) ~
        "Tier 1: direct clock metabolite + full relay",

      # Circadian-associated metabolite requires specific genus-metabolite support
      mev$ClockEvidenceScore == 2 &
        pev$PairEvidenceScore >= 2 &
        all_edges_p10 ~
        "Tier 2: circadian metabolite + supported microbial pair",

      # Very strong data may retain a literature genus, but remains exploratory
      mev$ClockEvidenceScore >= 2 &
        litg &
        n_edges_p05 == 3 ~
        "Tier 3: strong data + literature genus",

      TRUE ~
        "Data-driven only"
    )

    axis_rows[[k]] <- data.frame(
      Genus=g,
      CanonicalName=m,
      MetaboliteName=m_name,
      MetaboliteClass=shared_candidate$MetaboliteClass[i],

      ClockEvidenceScore=mev$ClockEvidenceScore,
      ClockEvidenceLevel=mev$ClockEvidenceLevel,
      ClockEvidenceReference=mev$ClockEvidenceReference,

      PairEvidenceScore=pev$PairEvidenceScore,
      PairEvidence=pev$PairEvidence,
      PairEvidenceReference=pev$PairEvidenceReference,

      Rho_Genus_Fecal=gf1$Rho,
      P_Genus_Fecal=gf1$P,
      Rho_Fecal_Serum=fs1$Rho,
      P_Fecal_Serum=fs1$P,
      Rho_Serum_Clock=sc1$Rho,
      P_Serum_Clock=sc1$P,

      LiteratureGenus=litg,
      PriorFig3Genus=priorg,
      PriorFig3Metabolite=shared_candidate$PriorFig3[i],

      N_Edges_P10=n_edges_p10,
      N_Edges_P05=n_edges_p05,
      AllEdges_P10=all_edges_p10,
      MainNetworkTier=main_tier,
      RawCorrelationProduct=raw_score,
      AxisScore=axis_score,
      stringsAsFactors=FALSE
    )

    k <- k+1
  }
}

axis4 <- if(length(axis_rows)>0){

  dplyr::bind_rows(axis_rows) %>%
    dplyr::arrange(
      factor(
        MainNetworkTier,
        levels=c(
          "Tier 1: direct clock metabolite + full relay",
          "Tier 2: circadian metabolite + supported microbial pair",
          "Tier 3: strong data + literature genus",
          "Data-driven only"
        )
      ),
      dplyr::desc(AxisScore)
    )

} else {
  data.frame()
}

if(nrow(axis4)>0){

  readr::write_csv(
    axis4,
    file.path(
      TABLE_DIR,
      "Candidate_4Layer_Axes_DataDriven.csv"
    )
  )
}

axis4_main <- if(nrow(axis4)>0){

  axis4 %>%
    dplyr::filter(
      MainNetworkTier != "Data-driven only"
    ) %>%
    # Do not let one metabolite generate three near-identical axes.
    dplyr::group_by(CanonicalName) %>%
    dplyr::slice_max(
      order_by=AxisScore,
      n=1,
      with_ties=FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(
      factor(
        MainNetworkTier,
        levels=c(
          "Tier 1: direct clock metabolite + full relay",
          "Tier 2: circadian metabolite + supported microbial pair",
          "Tier 3: strong data + literature genus"
        )
      ),
      dplyr::desc(AxisScore)
    )

} else {
  data.frame()
}

if(nrow(axis4_main)>0){

  readr::write_csv(
    axis4_main,
    file.path(
      TABLE_DIR,
      "Candidate_4Layer_Axes_MainEligible.csv"
    )
  )
}

# ---------- 3-layer fallback ----------
axis3_rows <- list(); k <- 1
for(g in candidate_genera) for(m in cand_serum){
  gs1 <- gs %>% dplyr::filter(Genus==g,CanonicalName==m) %>% dplyr::slice_head(n=1)
  sc1 <- sc %>% dplyr::filter(CanonicalName==m,ClockFeature=="Clock module") %>% dplyr::slice_head(n=1)
  if(nrow(gs1)==0||nrow(sc1)==0) next
  mp <- serum_c$map %>% dplyr::filter(CanonicalName==m) %>% dplyr::slice_head(n=1)
  score <- abs(gs1$Rho)*abs(sc1$Rho)+.12*(mp$MetaboliteClass[1]!="Other")+.10*(g%in%LITERATURE_TAXA)
  axis3_rows[[k]] <- data.frame(Genus=g,CanonicalName=m,MetaboliteName=mp$DisplayName[1],MetaboliteClass=mp$MetaboliteClass[1],
    Rho_Genus_Serum=gs1$Rho,P_Genus_Serum=gs1$P,Rho_Serum_Clock=sc1$Rho,P_Serum_Clock=sc1$P,AxisScore=score); k<-k+1
}
axis3 <- if(length(axis3_rows)>0) dplyr::bind_rows(axis3_rows) %>% dplyr::mutate(
  N_NominalEdges=rowSums(cbind(P_Genus_Serum<.10,P_Serum_Clock<.10),na.rm=TRUE)) %>%
  dplyr::arrange(dplyr::desc(N_NominalEdges),dplyr::desc(AxisScore)) else data.frame()
if(nrow(axis3)>0) readr::write_csv(axis3,file.path(TABLE_DIR,"Candidate_3Layer_Fallback_Axes.csv"))
top_axis <- if(nrow(axis4_main)>0) {
  axis4_main %>% dplyr::slice_head(n=1)
} else {
  NULL
}

if(!is.null(top_axis)){
  readr::write_csv(
    top_axis,
    file.path(TABLE_DIR,"TOP_Candidate_Axis_MAIN.csv")
  )
}

# Keep the best purely data-driven axis separately; never use it automatically
# as a mechanistic main-figure axis.
top_axis_data_driven <- if(nrow(axis4)>0) {
  axis4 %>% dplyr::slice_max(order_by=AxisScore,n=1,with_ties=FALSE)
} else if(nrow(axis3)>0) {
  axis3 %>% dplyr::slice_max(order_by=AxisScore,n=1,with_ties=FALSE)
} else {
  NULL
}

if(!is.null(top_axis_data_driven)){
  readr::write_csv(
    top_axis_data_driven,
    file.path(TABLE_DIR,"TOP_Candidate_Axis_DATA_DRIVEN_ONLY.csv")
  )
}

# ---------- Fig5B | balanced fecal + serum clock-candidate screen ----------
#
# The previous version could be dominated entirely by fecal metabolites because
# the final top-N cut was applied after combining fecal and serum candidates.
#
# v3 selects up to 8 fecal + 8 serum metabolites independently.
#
# Symbols in treatment cells:
#   · nominal P < 0.05 vs PCOS
#   * FDR < 0.10 vs PCOS
#   ** FDR < 0.05 vs PCOS
#
# This panel therefore displays BOTH effect pattern and statistical evidence.
# ----------------------------------------------------------------------

screen <- candidate_meta %>%
  dplyr::filter(
    ClockEvidenceScore >= 1 |
      PriorFig3 |
      MinFDR < .10
  ) %>%
  dplyr::group_by(Omic) %>%
  dplyr::arrange(
    dplyr::desc(ClockEvidenceScore),
    dplyr::desc(PriorFig3),
    MinFDR,
    dplyr::desc(MaxAbsEffect),
    .by_group=TRUE
  ) %>%
  dplyr::slice_head(n=8) %>%
  dplyr::ungroup()

rows <- list()
k <- 1

for(i in seq_len(nrow(screen))){

  om <- screen$Omic[i]
  m <- screen$CanonicalName[i]

  mat <- if(om=="Fecal metabolite") fecal_c$mat else serum_c$mat

  if(!m %in% rownames(mat)) next

  vals <- as.numeric(mat[m,])
  sam <- colnames(mat)
  gr <- factor(group_from_sample(sam),levels=GROUP_LEVELS)

  for(g in GROUP_LEVELS){

    contrast_name <- dplyr::case_when(
      g=="LCD" ~ "LCD_vs_PCOS",
      g=="GFC" ~ "GFC_vs_PCOS",
      g=="Combine" ~ "Combine_vs_PCOS",
      TRUE ~ NA_character_
    )

    st <- metab_diff %>%
      dplyr::filter(
        Omic==om,
        Feature==m,
        Contrast==contrast_name
      ) %>%
      dplyr::slice_head(n=1)

    sig <- if(g=="PCOS" || nrow(st)==0){
      ""
    } else if(!is.na(st$adj.P.Val) && st$adj.P.Val < .05){
      "**"
    } else if(!is.na(st$adj.P.Val) && st$adj.P.Val < .10){
      "*"
    } else if(!is.na(st$P.Value) && st$P.Value < .05){
      "·"
    } else {
      ""
    }

    evidence_short <- dplyr::case_when(
      screen$ClockEvidenceScore[i] == 3 ~ "[Direct]",
      screen$ClockEvidenceScore[i] == 2 ~ "[Circadian]",
      screen$ClockEvidenceScore[i] == 1 ~ "[Pathway]",
      TRUE ~ ""
    )

    rows[[k]] <- data.frame(
      Omic=om,
      DisplayName=screen$DisplayName[i],
      MetaboliteClass=screen$MetaboliteClass[i],
      ClockEvidenceLevel=screen$ClockEvidenceLevel[i],
      EvidenceShort=evidence_short,
      Group=g,
      MeanZ=mean(vals[gr==g],na.rm=TRUE),
      Sig=sig,
      stringsAsFactors=FALSE
    )

    k <- k+1
  }
}

screen_df <- dplyr::bind_rows(rows) %>%
  dplyr::mutate(
    Label=paste0(
      EvidenceShort,
      " ",
      stringr::str_wrap(DisplayName,26)
    ),
    Group=factor(Group,levels=GROUP_LEVELS),
    Omic=factor(
      Omic,
      levels=c("Fecal metabolite","Serum metabolite")
    )
  )

if(nrow(screen_df)>0){

  # preserve within-omic order
  ord <- screen %>%
    dplyr::mutate(
      Label=paste0(
        dplyr::case_when(
          ClockEvidenceScore==3 ~ "[Direct]",
          ClockEvidenceScore==2 ~ "[Circadian]",
          ClockEvidenceScore==1 ~ "[Pathway]",
          TRUE ~ ""
        ),
        " ",
        stringr::str_wrap(DisplayName,26)
      )
    ) %>%
    dplyr::pull(Label) %>%
    unique()

  screen_df$Label <- factor(
    screen_df$Label,
    levels=rev(ord)
  )

  lim <- max(abs(screen_df$MeanZ),na.rm=TRUE)
  if(!is.finite(lim)||lim<1) lim <- 1

  pB <- ggplot2::ggplot(
    screen_df,
    ggplot2::aes(Group,Label,fill=MeanZ)
  ) +
    ggplot2::geom_tile(
      colour="white",
      linewidth=.45
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label=Sig),
      size=4.0
    ) +
    ggplot2::facet_grid(
      Omic~.,
      scales="free_y",
      space="free_y",
      switch="y"
    ) +
    ggplot2::scale_fill_gradient2(
      low=COR_LOW,
      mid="white",
      high=COR_HIGH,
      midpoint=0,
      limits=c(-lim,lim),
      name="Mean z"
    ) +
    ggplot2::labs(
      x=NULL,
      y=NULL,
      caption="[Direct] direct experimental clock effect; [Circadian] circadian-associated; [Pathway] pathway/reverse-direction evidence. Symbols compare each treatment with PCOS: · P<0.05, * FDR<0.10, ** FDR<0.05."
    ) +
    ggplot2::theme_classic(base_size=14) +
    ggplot2::theme(
      strip.placement="outside",
      strip.background=ggplot2::element_blank(),
      strip.text.y.left=ggplot2::element_text(
        angle=90,
        face="bold",
        size=10
      ),
      axis.text.x=ggplot2::element_text(
        face="bold",
        colour="black"
      ),
      axis.text.y=ggplot2::element_text(
        size=9.0,
        colour="black",
        lineheight=.95
      ),
      axis.ticks=ggplot2::element_blank(),
      axis.line=ggplot2::element_blank(),
      plot.caption=ggplot2::element_text(
        size=8.5,
        hjust=0
      ),
      plot.margin=ggplot2::margin(8,8,8,14)
    )

  save_plot(
    pB,
    "Fig5B_FecalSerum_ClockCandidateScreen_v3",
    MAIN_DIR,
    9.6,
    max(
      6.4,
      3.0+.30*length(unique(screen_df$Label))
    )
  )
}

# ---------- Fig5C | fecal-serum relay ----------
#
# MAIN panel uses RAW collapsed z-scores because the biological question is:
# "Do fecal and circulating levels move together across animals/groups?"
#
# Two correlation statistics are shown:
#   overall rho = ordinary Spearman correlation on raw z-scores
#   adj. rho     = Spearman correlation after removing group means
#
# Supplementary FigS5C shows the residualized scatterplots used for adj. rho.
# ----------------------------------------------------------------------

fs_evidence <- fs

fev <- clock_metabolite_evidence(
  fs_evidence$FecalName
)

fs_evidence <- dplyr::bind_cols(
  fs_evidence,
  fev
) %>%
  dplyr::mutate(
    Priority =
      2.5*ClockEvidenceScore +
      1.5*as.numeric(PriorFig3) +
      abs(Rho) +
      -log10(pmax(P,1e-6))/4
  )

# IMPORTANT: These are the six originally submitted Fig. 4G candidates.
# Fix membership and order so evidence-level corrections do not silently
# replace a compound with another higher-ranked candidate.
FIG4G_FIXED_NAMES <- c(
  "5-Hydroxyindole-3-Acetic Acid",
  "Deoxycholic Acid",
  "Indole-3-Lactic Acid",
  "Kynurenic Acid",
  "Propionic Acid",
  "Nordeoxycholic Acid"
)
FIG4G_FIXED_CANONICAL <- normalize_metabolite_name(FIG4G_FIXED_NAMES)
if (anyDuplicated(fs_evidence$CanonicalName)) {
  stop("Duplicate canonical metabolite names in fecal-plasma candidate table. Resolve feature IDs before plotting Fig. 4G.")
}
idx_fixed <- match(FIG4G_FIXED_CANONICAL, fs_evidence$CanonicalName)
if (anyNA(idx_fixed)) {
  stop("Original Fig. 4G candidate(s) missing from matched fecal-plasma data: ",
       paste(FIG4G_FIXED_NAMES[is.na(idx_fixed)], collapse=", "),
       ". Do not auto-substitute another metabolite.")
}
top_shared <- fs_evidence[idx_fixed, , drop=FALSE]
stopifnot(nrow(top_shared)==6L,
          identical(as.character(top_shared$CanonicalName), FIG4G_FIXED_CANONICAL))
message("Fig. 4G: using all six original, fixed candidate names (no automatic replacement).")

rows_raw <- list()
rows_resid <- list()
k <- 1

for(i in seq_len(nrow(top_shared))){

  m <- top_shared$CanonicalName[i]

  sam <- intersect(
    colnames(fecal_c$mat),
    colnames(serum_c$mat)
  )

  grp <- factor(
    group_from_sample(sam),
    levels=GROUP_LEVELS
  )

  fec <- as.numeric(
    fecal_c$mat[m,sam]
  )

  ser <- as.numeric(
    serum_c$mat[m,sam]
  )

  raw_ok <- is.finite(fec) & is.finite(ser)

  raw_ct <- if(sum(raw_ok)>=6){
    suppressWarnings(
      stats::cor.test(
        fec[raw_ok],
        ser[raw_ok],
        method="spearman",
        exact=FALSE
      )
    )
  } else {
    NULL
  }

  overall_rho <- if(is.null(raw_ct)) NA_real_ else unname(raw_ct$estimate)
  overall_p <- if(is.null(raw_ct)) NA_real_ else raw_ct$p.value

  adj <- adjusted_spearman(
    fec,
    ser,
    grp
  )

  fec_r <- residualize_by_group(
    fec,
    grp
  )

  ser_r <- residualize_by_group(
    ser,
    grp
  )

  evidence_label <- dplyr::case_when(
    top_shared$ClockEvidenceScore[i]==3 ~ "Direct clock evidence",
    top_shared$ClockEvidenceScore[i]==2 ~ "Circadian-associated",
    top_shared$ClockEvidenceScore[i]==1 ~ "Pathway/reverse only",
    TRUE ~ "Data-driven"
  )

  panel_label <- paste0(
    stringr::str_wrap(top_shared$FecalName[i],22),
    "\n",
    evidence_label
  )

  stat_label <- paste0(
    "overall rho=",
    sprintf("%.2f",overall_rho),
    ", P=",
    format.pval(overall_p,digits=2,eps=.001),
    "\nadj. rho=",
    sprintf("%.2f",adj$Rho),
    ", P=",
    format.pval(adj$P,digits=2,eps=.001)
  )

  rows_raw[[k]] <- data.frame(
    Sample=sam,
    Group=grp,
    Fecal=fec,
    Serum=ser,
    Label=panel_label,
    OverallRho=overall_rho,
    OverallP=overall_p,
    AdjustedRho=adj$Rho,
    AdjustedP=adj$P,
    StatLabel=stat_label,
    stringsAsFactors=FALSE
  )

  rows_resid[[k]] <- data.frame(
    Sample=sam,
    Group=grp,
    FecalResidual=fec_r,
    SerumResidual=ser_r,
    Label=panel_label,
    AdjustedRho=adj$Rho,
    AdjustedP=adj$P,
    StatLabel=paste0(
      "adj. rho=",
      sprintf("%.2f",adj$Rho),
      "\nP=",
      format.pval(adj$P,digits=2,eps=.001)
    ),
    stringsAsFactors=FALSE
  )

  k <- k+1
}

if(length(rows_raw)>0){

  sd_raw <- dplyr::bind_rows(rows_raw)

  ann_raw <- sd_raw %>%
    dplyr::group_by(Label) %>%
    dplyr::summarise(
      StatLabel=dplyr::first(StatLabel),
      .groups="drop"
    )

  pC <- ggplot2::ggplot(
    sd_raw,
    ggplot2::aes(
      Fecal,
      Serum,
      colour=Group
    )
  ) +
    ggplot2::geom_point(
      size=2.4,
      alpha=.88
    ) +
    ggplot2::geom_smooth(
      method="lm",
      se=FALSE,
      colour="grey50",
      linewidth=.65
    ) +
    ggplot2::facet_wrap(
      ~Label,
      scales="free",
      ncol=3
    ) +
    ggplot2::scale_colour_manual(
      values=GROUP_COLS
    ) +
    ggplot2::labs(
      x="Fecal metabolite (z-score)",
      y="Plasma metabolite (z-score)",
      colour=NULL,
      caption="Main display uses raw z-scores to preserve treatment-group patterns. Both overall and group-adjusted Spearman correlations are reported."
    ) +
    ggplot2::theme_classic(base_size=14) +
    ggplot2::theme(
      legend.position="top",
      strip.background=ggplot2::element_blank(),
      strip.text=ggplot2::element_text(
        face="bold",
        size=9.2,
        lineheight=.95
      ),
      axis.text=ggplot2::element_text(
        colour="black"
      ),
      axis.title=ggplot2::element_text(
        face="bold"
      ),
      plot.caption=ggplot2::element_text(
        size=8.5,
        hjust=0
      )
    ) +
    ggplot2::geom_text(
      data=ann_raw,
      ggplot2::aes(
        x=-Inf,
        y=Inf,
        label=StatLabel
      ),
      inherit.aes=FALSE,
      hjust=-.04,
      vjust=1.05,
      size=2.75,
      colour="grey25"
    )

  save_plot(
    pC,
    "Fig5C_FecalSerum_ClockCandidateRelay_RAW_v4",
    MAIN_DIR,
    9.5,
    6.4
  )

  # ---- Final Figure 4G: exactly six metabolites in ONE horizontal row ----
  # Use the original six fixed metabolite names, not the old multi-line
  # evidence-tier titles. Evidence tiers remain in the candidate metadata.
  sd_res <- dplyr::bind_rows(rows_resid)
  original_labels <- unique(as.character(sd_res$Label))
  if (length(original_labels) != 6L) {
    stop("Figure 4G requires exactly six candidate panels; found ",
         length(original_labels), ". Inspect candidate matching and pairing.")
  }
  # Readable, complete names. Line breaks are DISPLAY ONLY.
  pretty_labels <- c(
    "5-Hydroxyindole-\n3-Acetic Acid",
    "Deoxycholic Acid",
    "Indole-3-\nLactic Acid",
    "Kynurenic Acid",
    "Propionic Acid",
    "Nordeoxycholic\nAcid"
  )
  # The input order is fixed by FIG4G_FIXED_NAMES above.
  sd_res$FacetName <- factor(sd_res$Label,
                             levels=original_labels,
                             labels=pretty_labels)
  # Raw paired observations and group-adjusted statistics are still retained
  # in the analysis outputs; no statistical tests are changed here.
  pSC <- ggplot2::ggplot(
    sd_res,
    ggplot2::aes(x=FecalResidual, y=SerumResidual, colour=Group)
  ) +
    ggplot2::geom_point(size=1.65, alpha=0.88) +
    ggplot2::geom_smooth(
      ggplot2::aes(group=1), method="lm", se=FALSE,
      colour="grey50", linewidth=0.45
    ) +
    ggplot2::facet_wrap(~ FacetName, scales="free", nrow=1, ncol=6) +
    ggplot2::scale_colour_manual(values=GROUP_COLS, drop=FALSE) +
    ggplot2::labs(
      x="Group-adjusted fecal residual",
      y="Group-adjusted plasma residual",
      colour=NULL, caption=NULL
    ) +
    ggplot2::theme_classic(base_size=14) +
    ggplot2::theme(
      legend.position="top",
      legend.text=ggplot2::element_text(size=9),
      legend.key.width=grid::unit(0.38,"cm"),
      legend.spacing.x=grid::unit(0.10,"cm"),
      strip.background=ggplot2::element_blank(),
      strip.text=ggplot2::element_text(
        face="bold", size=8, lineheight=0.92,
        margin=ggplot2::margin(t=4,b=5)
      ),
      axis.text=ggplot2::element_text(size=8.5,colour="black"),
      axis.title=ggplot2::element_text(size=11,face="bold"),
      plot.caption=ggplot2::element_blank(),
      panel.spacing.x=grid::unit(0.8,"lines"),
      plot.margin=ggplot2::margin(10,12,8,10)
    )

  # The save_plot helper applies x0.70 exactly once.
  # Final canvas = 19.5 x 4.8 inches * 0.70 = 13.65 x 3.36 inches.
  # It writes figures/Fig4G.pdf and Fig4G.png, not intermediate names.
  save_plot(
    pSC,
    "FigS5C_FecalSerum_ClockCandidateRelay_Adjusted_v4",
    SUPP_DIR,
    width=14.5,
    height=4.3,
    dpi=400
  )

  readr::write_csv(
    top_shared,
    file.path(
      TABLE_DIR,
      "Fig5C_Selected_FecalSerum_Metabolites_v4.csv"
    )
  )
}

# ----------------------------------------------------------------------
# DATA-DERIVED fecal/serum concordance + prior-changed-genus table
#
# This is the table to use when deciding the Figure 5 story.
# It only evaluates metabolites detected in BOTH feces and serum.
#
# Treatment direction:
#   +1 = higher than PCOS
#   -1 = lower than PCOS
#
# Concordant = fecal and serum effects have the same sign.
#
# "TopPriorGenus" is restricted to genera found in the previous Figure 3
# candidate list when that file is available; otherwise it falls back to
# literature-priority genera detected in the present dataset.
# ----------------------------------------------------------------------

get_effect <- function(diff_tbl, omic, feature, contrast){

  z <- diff_tbl %>%
    dplyr::filter(
      Omic==omic,
      Feature==feature,
      Contrast==contrast
    ) %>%
    dplyr::slice_head(n=1)

  if(nrow(z)==0){
    return(
      data.frame(
        logFC=NA_real_,
        P=NA_real_,
        FDR=NA_real_
      )
    )
  }

  data.frame(
    logFC=z$logFC[1],
    P=z$P.Value[1],
    FDR=z$adj.P.Val[1]
  )
}

prior_genus_pool <- intersect(
  prior_genera,
  candidate_genera
)

if(length(prior_genus_pool)==0){
  prior_genus_pool <- intersect(
    LITERATURE_TAXA,
    candidate_genera
  )
}

concord_rows <- list()
kc <- 1

for(m in shared){

  fmap <- fecal_c$map %>%
    dplyr::filter(CanonicalName==m) %>%
    dplyr::slice_head(n=1)

  smap <- serum_c$map %>%
    dplyr::filter(CanonicalName==m) %>%
    dplyr::slice_head(n=1)

  if(nrow(fmap)==0 || nrow(smap)==0) next

  display_name <- fmap$DisplayName[1]
  mev <- clock_metabolite_evidence(display_name)

  effs <- list()

  for(ct in c(
    "LCD_vs_PCOS",
    "GFC_vs_PCOS",
    "Combine_vs_PCOS"
  )){

    ef <- get_effect(
      metab_diff,
      "Fecal metabolite",
      m,
      ct
    )

    es <- get_effect(
      metab_diff,
      "Serum metabolite",
      m,
      ct
    )

    effs[[ct]] <- c(
      Fecal_logFC=ef$logFC,
      Fecal_P=ef$P,
      Fecal_FDR=ef$FDR,
      Serum_logFC=es$logFC,
      Serum_P=es$P,
      Serum_FDR=es$FDR,
      Concordant=ifelse(
        is.finite(ef$logFC) & is.finite(es$logFC),
        sign(ef$logFC)==sign(es$logFC),
        NA
      )
    )
  }

  fs1 <- fs %>%
    dplyr::filter(CanonicalName==m) %>%
    dplyr::slice_head(n=1)

  sc1 <- sc %>%
    dplyr::filter(
      CanonicalName==m,
      ClockFeature=="Clock module"
    ) %>%
    dplyr::slice_head(n=1)

  # Pick strongest association among previously highlighted genera only.
  genus_hit <- gf %>%
    dplyr::filter(
      CanonicalName==m,
      Genus %in% prior_genus_pool
    ) %>%
    dplyr::arrange(
      P,
      dplyr::desc(abs(Rho))
    ) %>%
    dplyr::slice_head(n=1)

  top_genus <- if(nrow(genus_hit)>0) genus_hit$Genus[1] else NA_character_
  genus_rho <- if(nrow(genus_hit)>0) genus_hit$Rho[1] else NA_real_
  genus_p <- if(nrow(genus_hit)>0) genus_hit$P[1] else NA_real_

  pev <- if(!is.na(top_genus)){
    genus_metabolite_pair_evidence(
      top_genus,
      display_name
    )
  } else {
    data.frame(
      PairEvidenceScore=0,
      PairEvidence="No prior-changed genus matched",
      PairEvidenceReference=""
    )
  }

  concord_n <- sum(
    c(
      effs$LCD_vs_PCOS["Concordant"],
      effs$GFC_vs_PCOS["Concordant"],
      effs$Combine_vs_PCOS["Concordant"]
    ) == 1,
    na.rm=TRUE
  )

  row <- data.frame(
    Metabolite=display_name,
    CanonicalName=m,
    MetaboliteClass=fmap$MetaboliteClass[1],
    ClockEvidenceScore=mev$ClockEvidenceScore,
    ClockEvidenceLevel=mev$ClockEvidenceLevel,
    ClockEvidenceReference=mev$ClockEvidenceReference,

    LCD_Fecal_logFC=effs$LCD_vs_PCOS["Fecal_logFC"],
    LCD_Serum_logFC=effs$LCD_vs_PCOS["Serum_logFC"],
    LCD_Concordant=as.logical(effs$LCD_vs_PCOS["Concordant"]),

    GFC_Fecal_logFC=effs$GFC_vs_PCOS["Fecal_logFC"],
    GFC_Serum_logFC=effs$GFC_vs_PCOS["Serum_logFC"],
    GFC_Concordant=as.logical(effs$GFC_vs_PCOS["Concordant"]),

    Combine_Fecal_logFC=effs$Combine_vs_PCOS["Fecal_logFC"],
    Combine_Serum_logFC=effs$Combine_vs_PCOS["Serum_logFC"],
    Combine_Concordant=as.logical(effs$Combine_vs_PCOS["Concordant"]),

    N_Concordant_Treatments=concord_n,

    FecalSerum_Rho=if(nrow(fs1)>0) fs1$Rho[1] else NA_real_,
    FecalSerum_P=if(nrow(fs1)>0) fs1$P[1] else NA_real_,

    SerumClock_Rho=if(nrow(sc1)>0) sc1$Rho[1] else NA_real_,
    SerumClock_P=if(nrow(sc1)>0) sc1$P[1] else NA_real_,

    TopPriorGenus=top_genus,
    GenusFecal_Rho=genus_rho,
    GenusFecal_P=genus_p,

    PairEvidenceScore=pev$PairEvidenceScore,
    PairEvidence=pev$PairEvidence,
    PairEvidenceReference=pev$PairEvidenceReference,

    stringsAsFactors=FALSE
  )

  concord_rows[[kc]] <- row
  kc <- kc+1
}

concord_table <- dplyr::bind_rows(concord_rows) %>%
  dplyr::mutate(
    StoryPriorityScore =
      2.0*ClockEvidenceScore +
      1.5*N_Concordant_Treatments +
      1.0*PairEvidenceScore +
      ifelse(is.finite(FecalSerum_Rho),abs(FecalSerum_Rho),0) +
      ifelse(is.finite(SerumClock_Rho),abs(SerumClock_Rho),0) +
      ifelse(is.finite(GenusFecal_Rho),abs(GenusFecal_Rho),0)
  ) %>%
  dplyr::arrange(
    dplyr::desc(ClockEvidenceScore),
    dplyr::desc(N_Concordant_Treatments),
    dplyr::desc(PairEvidenceScore),
    dplyr::desc(StoryPriorityScore)
  )

readr::write_csv(
  concord_table,
  file.path(
    TABLE_DIR,
    "Table_Fig5_FecalSerumConcordant_ClockMetabolite_PriorGenus.csv"
  )
)

# A compact manuscript-facing shortlist.
concord_shortlist <- concord_table %>%
  dplyr::filter(
    ClockEvidenceScore >= 1,
    N_Concordant_Treatments >= 1,
    !is.na(TopPriorGenus)
  ) %>%
  dplyr::slice_head(n=12)

readr::write_csv(
  concord_shortlist,
  file.path(
    TABLE_DIR,
    "Table_Fig5_STORY_SHORTLIST.csv"
  )
)

# ---------- Fig5D | evidence-tiered main network ----------
#
# Main figure uses axis4_main only.
# Kynurenic acid is classified as "circadian-associated", not a direct
# clock-regulating metabolite. Therefore a KYNA axis needs specific
# microorganism-KYNA support plus a complete nominally-supported relay
# before entering the main network.
#
# Turicimonas-KYNA has no specific literature support in the embedded
# evidence table and therefore cannot enter via the pair-evidence route.
# ----------------------------------------------------------------------

if(nrow(axis4_main)>0){

  chosen <- axis4_main %>%
    dplyr::slice_head(
      n=min(3,nrow(axis4_main))
    )

  readr::write_csv(
    chosen,
    file.path(
      TABLE_DIR,
      "Fig5D_MainNetwork_SelectedAxes_v3.csv"
    )
  )

  edges <- list()
  k <- 1

  for(i in seq_len(nrow(chosen))){

    g <- chosen$Genus[i]
    lab <- chosen$MetaboliteName[i]
    f <- paste0("Fecal | ",lab)
    ss <- paste0("Serum | ",lab)
    cl <- "Ovarian clock module"

    edges[[k]] <- data.frame(
      From=g,
      To=f,
      Rho=chosen$Rho_Genus_Fecal[i],
      P=chosen$P_Genus_Fecal[i],
      Evidence=chosen$PairEvidence[i],
      stringsAsFactors=FALSE
    )
    k <- k+1

    edges[[k]] <- data.frame(
      From=f,
      To=ss,
      Rho=chosen$Rho_Fecal_Serum[i],
      P=chosen$P_Fecal_Serum[i],
      Evidence="Same annotated metabolite in feces and serum",
      stringsAsFactors=FALSE
    )
    k <- k+1

    edges[[k]] <- data.frame(
      From=ss,
      To=cl,
      Rho=chosen$Rho_Serum_Clock[i],
      P=chosen$P_Serum_Clock[i],
      Evidence=chosen$ClockEvidenceLevel[i],
      stringsAsFactors=FALSE
    )
    k <- k+1
  }

  ed <- dplyr::bind_rows(edges) %>%
    dplyr::distinct(
      From,
      To,
      .keep_all=TRUE
    )

  nn <- unique(
    c(ed$From,ed$To)
  )

  nd <- data.frame(
    Name=nn,
    stringsAsFactors=FALSE
  ) %>%
    dplyr::mutate(
      Type=dplyr::case_when(
        grepl("^Fecal \\|",Name) ~ "Fecal metabolite",
        grepl("^Serum \\|",Name) ~ "Serum metabolite",
        Name=="Ovarian clock module" ~ "Ovarian clock",
        TRUE ~ "Genus"
      ),
      X=dplyr::case_when(
        Type=="Genus" ~ 1,
        Type=="Fecal metabolite" ~ 2,
        Type=="Serum metabolite" ~ 3,
        TRUE ~ 4
      )
    ) %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(
      Y=seq(1,-1,length.out=dplyr::n())
    ) %>%
    dplyr::ungroup()

  ep <- ed %>%
    dplyr::left_join(
      nd %>% dplyr::select(From=Name,X1=X,Y1=Y),
      by="From"
    ) %>%
    dplyr::left_join(
      nd %>% dplyr::select(To=Name,X2=X,Y2=Y),
      by="To"
    ) %>%
    dplyr::mutate(
      Sign=ifelse(Rho>=0,"Positive","Negative"),
      Alpha=ifelse(
        P<.05,
        1,
        ifelse(P<.10,.72,.30)
      )
    )

  pD <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data=ep,
      ggplot2::aes(
        x=X1,
        y=Y1,
        xend=X2,
        yend=Y2,
        colour=Sign,
        linewidth=abs(Rho),
        alpha=Alpha
      ),
      lineend="round"
    ) +
    ggplot2::geom_point(
      data=nd,
      ggplot2::aes(
        X,Y,fill=Type
      ),
      shape=21,
      size=5.5,
      colour="grey25",
      stroke=.5
    ) +
    ggrepel::geom_text_repel(
      data=nd,
      ggplot2::aes(
        X,Y,
        label=stringr::str_wrap(
          gsub("^(Fecal \\| |Serum \\| )","",Name),
          18
        )
      ),
      size=3.3,
      direction="y",
      box.padding=.45,
      point.padding=.35,
      min.segment.length=0,
      segment.colour="grey65",
      max.overlaps=Inf,
      seed=20261001
    ) +
    ggplot2::scale_fill_manual(
      values=OMICS_COLS,
      name=NULL
    ) +
    ggplot2::scale_colour_manual(
      values=c(
        Positive=COR_HIGH,
        Negative=COR_LOW
      ),
      name="Correlation"
    ) +
    ggplot2::scale_linewidth_continuous(
      range=c(.5,2.2),
      guide="none"
    ) +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_x_continuous(
      breaks=1:4,
      labels=c(
        "Genus",
        "Fecal\nmetabolite",
        "Serum\nmetabolite",
        "Ovarian\nclock"
      )
    ) +
    ggplot2::labs(
      x=NULL,
      y=NULL,
      caption="Main network requires evidence-tier eligibility plus group-adjusted correlations. Red/blue denote positive/negative association; this remains a candidate axis, not causal proof."
    ) +
    ggplot2::theme_classic(base_size=14) +
    ggplot2::theme(
      axis.line.y=ggplot2::element_blank(),
      axis.ticks=ggplot2::element_blank(),
      axis.text.y=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_text(
        face="bold",
        colour="black"
      ),
      legend.position="top",
      plot.caption=ggplot2::element_text(
        size=8.6,
        hjust=0
      ),
      plot.margin=ggplot2::margin(10,25,8,25)
    )

  save_plot(
    pD,
    "Fig5D_EvidenceTiered_GutMetaboliteClock_Network_v3",
    MAIN_DIR,
    9.3,
    5.8
  )

} else {

  writeLines(
    c(
      "No four-layer axis met the v3 evidence-tier criteria for the MAIN network.",
      "This is preferable to presenting a purely correlation-ranked axis as an established mechanism.",
      "Inspect Candidate_4Layer_Axes_DataDriven.csv for exploratory candidates."
    ),
    file.path(
      QC_DIR,
      "Fig5D_No_MainEligibleAxis_v3.txt"
    )
  )
}

# ---------- Fig5E | group pattern of TOP MAIN-ELIGIBLE axis ----------
if(!is.null(top_axis)){
  g<-top_axis$Genus[1]; m<-top_axis$CanonicalName[1]; lab<-top_axis$MetaboliteName[1]; rr<-list(); k<-1
  if(g%in%rownames(genus_clr)){v<-genus_clr[g,];rr[[k]]<-data.frame(Sample=names(v),Group=group_from_sample(names(v)),Feature=paste0("Genus\n",stringr::str_wrap(g,20)),Value=as.numeric(v));k<-k+1}
  if(m%in%rownames(fecal_c$mat)){v<-fecal_c$mat[m,];rr[[k]]<-data.frame(Sample=names(v),Group=group_from_sample(names(v)),Feature=paste0("Fecal\n",stringr::str_wrap(lab,20)),Value=as.numeric(v));k<-k+1}
  if(m%in%rownames(serum_c$mat)){v<-serum_c$mat[m,];rr[[k]]<-data.frame(Sample=names(v),Group=group_from_sample(names(v)),Feature=paste0("Serum\n",stringr::str_wrap(lab,20)),Value=as.numeric(v));k<-k+1}
  v<-clock_score;rr[[k]]<-data.frame(Sample=names(v),Group=group_from_sample(names(v)),Feature="Ovarian clock\nmodule",Value=as.numeric(v))
  td <- dplyr::bind_rows(rr)%>%dplyr::mutate(Group=factor(Group,levels=GROUP_LEVELS))
  pE <- ggplot2::ggplot(td,ggplot2::aes(Group,Value,fill=Group))+ggplot2::geom_boxplot(width=.58,outlier.shape=NA,alpha=.28,colour="grey35")+
    ggplot2::geom_jitter(ggplot2::aes(colour=Group),width=.08,size=2,alpha=.9)+ggplot2::facet_wrap(~Feature,scales="free_y",ncol=4)+
    ggplot2::scale_fill_manual(values=GROUP_COLS,guide="none")+ggplot2::scale_colour_manual(values=GROUP_COLS,guide="none")+
    ggplot2::labs(x=NULL,y="Feature value",caption="Metabolites: collapsed LC-MS z-scores; genus: CLR; clock module: mean clock-gene z-score.")+
    ggplot2::theme_classic(base_size=14)+ggplot2::theme(strip.background=ggplot2::element_blank(),strip.text=ggplot2::element_text(face="bold",size=9.5),
      axis.text.x=ggplot2::element_text(angle=40,hjust=1,colour="black",size=8.8),axis.text.y=ggplot2::element_text(colour="black"),
      plot.caption=ggplot2::element_text(size=8.5,hjust=0))
  save_plot(pE,"Fig5E_TopCandidateAxis_GroupPatterns",MAIN_DIR,11,4.5)
}

# ---------- Fig5F serum-clock heatmap ----------
hdf <- sc %>% dplyr::filter(ClockFeature!="Clock module") %>%
  dplyr::group_by(CanonicalName) %>% dplyr::mutate(MinP=min(P,na.rm=TRUE),MaxAbsRho=max(abs(Rho),na.rm=TRUE)) %>% dplyr::ungroup() %>%
  dplyr::left_join(candidate_meta %>% dplyr::filter(Omic=="Serum metabolite") %>% dplyr::select(CanonicalName,CandidateScore) %>%
    dplyr::distinct(CanonicalName,.keep_all=TRUE),by="CanonicalName") %>%
  dplyr::mutate(CandidateScore=ifelse(is.na(CandidateScore),0,CandidateScore),HeatPriority=1.5*(-log10(pmax(MinP,1e-6)))+MaxAbsRho+CandidateScore/4)
topm <- hdf %>% dplyr::distinct(CanonicalName,DisplayName,MetaboliteClass,HeatPriority) %>% dplyr::arrange(dplyr::desc(HeatPriority)) %>% dplyr::slice_head(n=12)
hdf <- hdf %>% dplyr::filter(CanonicalName%in%topm$CanonicalName) %>% dplyr::mutate(
  MetabLabel=stringr::str_wrap(DisplayName,26),ClockFeature=factor(ClockFeature,levels=intersect(CLOCK_DISPLAY_ORDER,unique(ClockFeature))),
  Sig=dplyr::case_when(FDR<.05~"**",FDR<.10~"*",P<.05~"·",TRUE~""))
if(nrow(hdf)>0){
  hdf$MetabLabel <- factor(hdf$MetabLabel,levels=rev(unique(stringr::str_wrap(topm$DisplayName,26))))
  pF <- ggplot2::ggplot(hdf,ggplot2::aes(ClockFeature,MetabLabel,fill=Rho))+ggplot2::geom_tile(colour="white",linewidth=.45)+
    ggplot2::geom_text(ggplot2::aes(label=Sig),size=3.8)+ggplot2::facet_grid(MetaboliteClass~.,scales="free_y",space="free_y",switch="y")+
    ggplot2::scale_fill_gradient2(low=COR_LOW,mid="white",high=COR_HIGH,midpoint=0,limits=c(-1,1),name="Adjusted\nSpearman rho")+
    ggplot2::labs(x=NULL,y=NULL,caption="Group-adjusted correlations. · P<0.05; * FDR<0.10; ** FDR<0.05. Association does not establish causality.")+
    ggplot2::theme_classic(base_size=14)+ggplot2::theme(strip.placement="outside",strip.background=ggplot2::element_blank(),
      strip.text.y.left=ggplot2::element_text(angle=90,face="bold",size=9.2),axis.text.x=ggplot2::element_text(angle=35,hjust=1,size=9.3,colour="black"),
      axis.text.y=ggplot2::element_text(size=9.2,colour="black"),axis.ticks=ggplot2::element_blank(),axis.line=ggplot2::element_blank(),
      plot.caption=ggplot2::element_text(size=8.6,hjust=0))
  save_plot(pF,"Fig5F_SerumMetabolite_OvarianClock_Correlation",MAIN_DIR,10,max(5.8,2.8+.32*length(unique(hdf$MetabLabel))))
}

# ---------- supplementary genus profiles ----------
gp <- head(candidate_genera,8); rr<-list();k<-1
for(g in gp){v<-genus_clr[g,];rr[[k]]<-data.frame(Group=factor(group_from_sample(names(v)),levels=GROUP_LEVELS),Genus=stringr::str_wrap(g,20),CLR=as.numeric(v));k<-k+1}
gpd<-dplyr::bind_rows(rr)
if(nrow(gpd)>0){
  p <- ggplot2::ggplot(gpd,ggplot2::aes(Group,CLR,fill=Group))+ggplot2::geom_boxplot(width=.58,outlier.shape=NA,alpha=.28,colour="grey35")+
    ggplot2::geom_jitter(ggplot2::aes(colour=Group),width=.08,size=1.9,alpha=.9)+ggplot2::facet_wrap(~Genus,scales="free_y",ncol=4)+
    ggplot2::scale_fill_manual(values=GROUP_COLS,guide="none")+ggplot2::scale_colour_manual(values=GROUP_COLS,guide="none")+
    ggplot2::labs(x=NULL,y="Genus CLR abundance")+ggplot2::theme_classic(base_size=14)+ggplot2::theme(strip.background=ggplot2::element_blank(),
      strip.text=ggplot2::element_text(face="bold",size=9.3),axis.text.x=ggplot2::element_text(angle=40,hjust=1,size=8.5,colour="black"))
  save_plot(p,"FigS5A_CandidateGenera_Abundance",SUPP_DIR,10.5,6.2)
}


# v3 evidence notes for manuscript interpretation
figure5_v3_evidence_notes <- data.frame(
  Item=c(
    "DCA/CDCA",
    "SCFAs",
    "3-PPA / 4-OH-PPA",
    "Kynurenic acid",
    "Prevotellaceae_UCG-001 - kynurenic acid",
    "Odoribacter - kynurenic acid",
    "Turicimonas - kynurenic acid"
  ),
  Evidence=c(
    "Direct experimental alteration of Dbp/Per2/Per3/Cry2 by unconjugated bile acids",
    "Direct entrainment/phase-shift evidence for PER2/BMAL1",
    "Direct modulation of PER2/Bmal1 oscillation",
    "Circadian cycling and AHR-ligand evidence; direct KYNA-to-ovarian-clock regulation not established",
    "Published positive microbiota-metabolite correlation; not proof of production",
    "No specific direct KYNA link identified; Odoribacter has broader tryptophan/indole associations",
    "No specific KYNA link identified"
  ),
  Reference=c(
    "PMID:27907092",
    "PMID:35934064",
    "PMID:32155689",
    "PMID:41727182; PMID:30706546",
    "PMCID:PMC9764441",
    "PMCID:PMC12606010",
    ""
  ),
  stringsAsFactors=FALSE
)

readr::write_csv(
  figure5_v3_evidence_notes,
  file.path(
    TABLE_DIR,
    "Figure5_v3_SpecificEvidenceNotes.csv"
  )
)

# ----------------------------------------------------------------------
# Literature-supported links between current clock-related metabolites
# and genera that were highlighted in earlier microbiome analyses.
#
# EvidenceLevel:
#   Strong   = direct production or repeatedly reported specific association
#   Moderate = published correlation with the exact metabolite / close BA class
#   Context  = broader pathway association; direction may be context-dependent
# ----------------------------------------------------------------------

literature_metabolite_microbe_pairs <- data.frame(
  Metabolite=c(
    "5-Hydroxyindole-3-Acetic Acid (5-HIAA)",
    "Kynurenic Acid",
    "Deoxycholic Acid (DCA)",
    "Nord(e)oxycholic Acid / NorDCA",
    "Cholic Acid / bile-acid pool",
    "Propionic Acid",
    "Propionic Acid"
  ),
  Genus=c(
    "Alistipes",
    "Prevotellaceae_UCG-001",
    "Mucispirillum",
    "Mucispirillum",
    "Alistipes",
    "[Eubacterium]_coprostanoligenes_group",
    "Alistipes"
  ),
  EvidenceLevel=c(
    "Strong",
    "Moderate",
    "Moderate",
    "Moderate",
    "Moderate",
    "Context",
    "Context"
  ),
  LiteratureEvidence=c(
    "Alistipes onderdonkii produced 5-HIAA in vitro; intervention increased Alistipes and 5-HIAA.",
    "Positive genus-metabolite correlation between Prevotellaceae_UCG-001 and kynurenic acid in 16S+metabolomics.",
    "Mucispirillum correlated positively with DCA and other secondary bile acids in mouse studies.",
    "Mucispirillum correlated positively with NorDCA/secondary bile acids in HFD-associated bile-acid profiling.",
    "Alistipes repeatedly associates with bile-acid profiles, including CA/secondary bile acids; exact direction is study-dependent.",
    "Eubacterium coprostanoligenes has been associated with SCFA profiles; exact propionate direction varies by cohort/diet.",
    "Alistipes-propionate association has been reported, but direction varies by biological context."
  ),
  Reference=c(
    "PMID:41692705",
    "PMCID:PMC9764441",
    "PMCID:PMC10531215; PMCID:PMC8819290",
    "PMCID:PMC8819290",
    "PMCID:PMC7802554; PMID:33362374",
    "PMID:37485501; PMID:38779015",
    "PMCID:PMC12484129; PMCID:PMC7428940"
  ),
  StoryPriority=c(
    "High",
    "High if your Prevotellaceae_UCG-001 is treatment-responsive",
    "High",
    "High",
    "Medium-high",
    "Low-medium",
    "Low-medium"
  ),
  stringsAsFactors=FALSE
)

readr::write_csv(
  literature_metabolite_microbe_pairs,
  file.path(
    TABLE_DIR,
    "Table_Fig5_LiteratureSupported_MetaboliteMicrobe_Pairs.csv"
  )
)

summary_lines <- c(
"============================================================",
"FIGURE 5 | GUT-METABOLITE-OVARIAN CLOCK AXIS",
"============================================================","",
paste0("Matched samples across all four omics: ",length(matched_all4)),
paste0("Detected clock genes: ",paste(rownames(clock_z),collapse="; ")),
paste0("Clock-module genes: ",paste(primary_detected,collapse="; ")),
paste0("Candidate genera: ",paste(candidate_genera,collapse="; ")),
paste0("Shared annotated fecal-serum metabolites: ",length(shared)),"",
if(nrow(axis4_main)>0) paste0("Top MAIN-eligible axis: ",axis4_main$Genus[1]," -> ",axis4_main$MetaboliteName[1],
" (fecal -> serum) -> ovarian clock module") else "No 4-layer axis met the evidence-tier criteria for the main network; inspect data-driven tables.","",
"INTERPRETATION:",
"- group-adjusted cross-sectional associations only",
"- correlation does not establish causality",
"- one-time-point RNA cannot prove circadian rhythm restoration/phase/amplitude change",
"- strongest story requires taxon change + fecal/serum relay + serum-clock association + independent literature support",
"- use Table_Fig5_STORY_SHORTLIST.csv to prioritize metabolites with concordant fecal/serum treatment directions and previously highlighted genera",
"")
writeLines(summary_lines,file.path(OUT_DIR,"SUMMARY_Figure5.txt"))
cat("\nFigure 5 completed. Output: ",OUT_DIR,"\n")

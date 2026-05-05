## ============================================================
##  run_GRaNIE.R — unified Human / Mouse GRaNIE script
##  Usage:
##    Rscript run_GRaNIE.R <rna_file> <atac_file> <out_dir> \
##                         <sample_name> <genome> <script_dir>
##
##  Arguments:
##    rna_file    : full path to RNA counts CSV (genes x cells)
##    atac_file   : full path to ATAC counts CSV (peaks x cells)
##    out_dir     : full path to output directory
##    sample_name : label used for GRN object metadata and output files
##    genome      : "hg38" for Human, "mm10" for Mouse
##    script_dir  : directory containing reference files:
##                  Human: EnsDb.Hsapiens.v86_annot.rds (required)
##                         ensembl_mart_hsapiens.rds     (required)
##                         TFBS folder path (set below)
##                  Mouse: EnsDb.Mmusculus.v79_UCSC_mm10.rds (required)
##                         ensembl_mart_mmusculus.rds         (required)
##                         TFBS folder path (set below)
## ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop(paste(
    "Usage: Rscript run_GRaNIE.R",
    "<rna_file> <atac_file> <out_dir>",
    "<sample_name> <genome> <script_dir>"
  ))
}

rna_file    <- args[1]
atac_file   <- args[2]
out_dir     <- args[3]
sample_name <- args[4]
genome      <- args[5]   # "hg38" or "mm10"
script_dir  <- args[6]

## ── Validate inputs ───────────────────────────────────────────
if (!file.exists(rna_file))  stop("RNA file not found: ",  rna_file)
if (!file.exists(atac_file)) stop("ATAC file not found: ", atac_file)
if (!genome %in% c("hg38", "mm10")) {
  stop("genome must be 'hg38' or 'mm10'. Got: ", genome)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("======================================")
message("  GRaNIE run")
message("  RNA file    : ", rna_file)
message("  ATAC file   : ", atac_file)
message("  Out dir     : ", out_dir)
message("  Sample name : ", sample_name)
message("  Genome      : ", genome)
message("  Script dir  : ", script_dir)
message("======================================")

## ── Libraries ────────────────────────────────────────────────
library(readr)
library(GRaNIE)
library(biomaRt)
library(Seurat)
library(Signac)
library(dplyr)

## ── Genome-specific settings ──────────────────────────────────
if (genome == "hg38") {

  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)

  genome_str     <- "hg38"
  mt_pattern     <- "^MT-"
  biomart_ds     <- "hsapiens_gene_ensembl"

  ## Chromatin annotations
  annot_rds <- file.path(script_dir, "EnsDb.Hsapiens.v86_annot.rds")
  if (!file.exists(annot_rds)) stop("Human annotation RDS not found: ", annot_rds)
  annotations <- readRDS(annot_rds)
  message("Loaded Human annotations from: ", annot_rds)

  ## Saved biomaRt mart (avoids live Ensembl connection on cluster)
  mart_rds <- file.path(script_dir, "ensembl_mart_hsapiens.rds")
  if (!file.exists(mart_rds)) stop("Human biomaRt RDS not found: ", mart_rds)
  ensembl_mart <- readRDS(mart_rds)
  message("Loaded Human biomaRt from: ", mart_rds)

  ## TFBS motif folder — HOCOMOCO v12 human
  tfbs_folder <- "/gpfs/Home/kmk7420/Multi_omics_GRN/GRaNIE/H12INVIVO"

} else {  # mm10

  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)

  genome_str     <- "mm10"
  mt_pattern     <- "^mt-"
  biomart_ds     <- "mmusculus_gene_ensembl"

  ## Chromatin annotations
  annot_rds <- file.path(script_dir, "EnsDb.Mmusculus.v79_UCSC_mm10.rds")
  if (!file.exists(annot_rds)) stop("Mouse annotation RDS not found: ", annot_rds)
  annotations <- readRDS(annot_rds)
  message("Loaded Mouse annotations from: ", annot_rds)

  ## Saved biomaRt mart
  mart_rds <- file.path(script_dir, "ensembl_mart_mmusculus.rds")
  if (!file.exists(mart_rds)) stop("Mouse biomaRt RDS not found: ", mart_rds)
  ensembl_mart <- readRDS(mart_rds)
  message("Loaded Mouse biomaRt from: ", mart_rds)

  ## TFBS motif folder — HOCOMOCO v12 mouse
  tfbs_folder <- "/gpfs/Home/kmk7420/Multi_omics_GRN/GRaNIE/mESC/PWMScan_HOCOMOCOv12/H12INVIVO/pwmscan_filt"

}

## ── Load raw counts ───────────────────────────────────────────
message("Loading RNA counts...")
rna_data  <- read.table(rna_file,  header = TRUE, row.names = 1, sep = ",", comment.char = "")

message("Loading ATAC counts...")
atac_data <- read.table(atac_file, header = TRUE, row.names = 1, sep = ",", comment.char = "")

## ── Auto-detect peak separator ────────────────────────────────
first_peak <- rownames(atac_data)[1]
sep_to_use <- if (grepl(":", first_peak)) c(":", "-") else c("-", "-")
message("Peak notation: ", first_peak, " -> sep: ", paste(sep_to_use, collapse = ","))

## ── Build Seurat object ───────────────────────────────────────
message("Building Seurat object...")
pbmc <- CreateSeuratObject(counts = rna_data)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = mt_pattern)

## Filter ATAC to standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_data), sep = sep_to_use)
grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_data     <- atac_data[as.vector(grange.use), ]
message("ATAC peaks after chr filter: ", nrow(atac_data))

chrom_assay <- CreateChromatinAssay(
  counts     = atac_data,
  sep        = sep_to_use,
  genome     = genome_str,
  min.cells  = 1,
  annotation = annotations
)
pbmc[["ATAC"]] <- chrom_assay

## ── RNA preprocessing ─────────────────────────────────────────
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE, return.only.var.genes = FALSE)
pbmc <- RunPCA(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:50)

## ── ATAC preprocessing ────────────────────────────────────────
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:50)

## ── WNN integration + clustering ──────────────────────────────
pbmc <- FindMultiModalNeighbors(pbmc,
                                reduction.list = list("pca", "lsi"),
                                dims.list      = list(1:50, 2:50))
#################################################################################
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap")
#pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", dims = 1:50)
seurat_obj <- FindClusters(pbmc, graph.name = "wknn", resolution = 10)

## ── Pseudobulk aggregation ────────────────────────────────────
message("Generating pseudobulk...")
pseudobulk     <- AggregateExpression(seurat_obj,
                                      assays = c("RNA", "ATAC"),
                                      slot   = "counts",
                                      fun    = mean)
countsRNA.df   <- as.data.frame(pseudobulk$RNA)
countsPeaks.df <- as.data.frame(pseudobulk$ATAC)

## Move row names into ID columns (GRaNIE requires explicit ID columns)
countsRNA.df$ENSEMBL   <- rownames(countsRNA.df)
countsPeaks.df$peakID  <- rownames(countsPeaks.df)

countsRNA.df   <- countsRNA.df[,   c("ENSEMBL",
                                      setdiff(colnames(countsRNA.df),   "ENSEMBL"))]
countsPeaks.df <- countsPeaks.df[, c("peakID",
                                      setdiff(colnames(countsPeaks.df), "peakID"))]

rownames(countsRNA.df)   <- NULL
rownames(countsPeaks.df) <- NULL

## ── Map gene symbols to Ensembl IDs ──────────────────────────
message("Mapping gene symbols to Ensembl IDs (", biomart_ds, ")...")
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters    = "external_gene_name",
  values     = countsRNA.df$ENSEMBL,
  mart       = ensembl_mart
)
mapping_vector         <- setNames(gene_mapping$ensembl_gene_id,
                                   gene_mapping$external_gene_name)
countsRNA.df$ENSEMBL   <- mapping_vector[countsRNA.df$ENSEMBL]

## Remove rows where mapping failed (NA Ensembl IDs)
na_count <- sum(is.na(countsRNA.df$ENSEMBL))
if (na_count > 0) {
  message("  Removing ", na_count, " genes with no Ensembl ID mapping")
  countsRNA.df <- countsRNA.df[!is.na(countsRNA.df$ENSEMBL), ]
}
message("  Genes after mapping: ", nrow(countsRNA.df))

## ── Initialise GRaNIE object ──────────────────────────────────
message("Initialising GRN object...")
objectMetadata.l <- list(
  name           = sample_name,
  file_peaks     = atac_file,
  file_rna       = rna_file,
  genomeAssembly = genome_str
)

GRN <- initializeGRN(
  objectMetadata = objectMetadata.l,
  outputFolder   = out_dir,
  genomeAssembly = genome_str
)

################################################################
rownames(countsRNA.df)   <- NULL
rownames(countsPeaks.df) <- NULL

################################################################

GRN <- addData(
  GRN,
  counts_peaks        = countsPeaks.df,
  normalization_peaks = "limma_quantile",
  idColumn_peaks      = "peakID",
  counts_rna          = countsRNA.df,
  normalization_rna   = "limma_quantile",
  idColumn_RNA        = "ENSEMBL",
  forceRerun          = TRUE
)

## ── PCA plot (RNA) ────────────────────────────────────────────
GRN <- plotPCA_all(GRN, data = c("rna"), topn = 500, type = "normalized",
                   plotAsPDF = FALSE, pages = c(2, 3, 8), forceRerun = TRUE)

## ── Add TFBS ──────────────────────────────────────────────────
message("Adding TFBS from: ", tfbs_folder)
GRN <- addTFBS(GRN,
               motifFolder       = tfbs_folder,
               TFs               = "all",
               filesTFBSPattern  = "_TFBS",
               fileEnding        = ".bed.gz",
               forceRerun        = TRUE)

## ── Overlap peaks and TFBS ────────────────────────────────────
GRN <- overlapPeaksAndTFBS(GRN, nCores = 12, forceRerun = TRUE)

## ── TF-peak connections ───────────────────────────────────────
GRN <- addConnections_TF_peak(
  GRN,
  plotDiagnosticPlots = FALSE,
  connectionTypes     = c("expression"),
  corMethod           = "spearman",
  maxFDRToStore       = 1,
  forceRerun          = TRUE
)

## ── AR classification ─────────────────────────────────────────
out_plots <- file.path(out_dir, "plots")
GRN <- AR_classification_wrapper(
  GRN,
  significanceThreshold_Wilcoxon = 1,
  outputFolder                   = out_plots,
  plot_minNoTFBS_heatmap         = 100,
  plotDiagnosticPlots            = TRUE,
  forceRerun                     = TRUE
)
saveRDS(GRN, file.path(out_dir, paste0(sample_name, "_GRN_AR_classification.rds")))

## ── Peak-gene connections ─────────────────────────────────────
GRN <- addConnections_peak_gene(
  GRN,
  corMethod           = "spearman",
  promoterRange       = 250000,
  TADs                = NULL,
  nCores              = 12,
  plotDiagnosticPlots = FALSE,
  forceRerun          = TRUE
)

## ── Filter GRN + TF-gene correlations ────────────────────────
GRN <- filterGRNAndConnectGenes(
  GRN,
  TF_peak.fdr.threshold   = 0.3,
  peak_gene.fdr.threshold = 0.3,
  peak_gene.fdr.method    = "BH",
  gene.types              = c("all"),
  forceRerun              = TRUE
)

GRN <- add_TF_gene_correlation(GRN, corMethod = "spearman", nCores = 12, forceRerun = TRUE)
saveRDS(GRN, file.path(out_dir, paste0(sample_name, "_GRN_addTF_gene_corr.rds")))

## ── Extract and save connections ──────────────────────────────
GRN_connections.all <- getGRNConnections(
  GRN,
  type                         = "all.filtered",
  include_TF_gene_correlations = TRUE,
  include_geneMetadata         = TRUE
)
message("GRN connections retrieved: ", nrow(GRN_connections.all), " rows")

## All columns — unfiltered
write.csv(GRN_connections.all,
          file.path(out_dir, paste0(sample_name, "_GRN_connections_all_unfiltered.csv")),
          row.names = FALSE, quote = FALSE)

## Selected columns (TF.name, gene.name, TF_gene.r, TF_gene.p_raw)
## col 2 = TF.name, col 16 = gene.name, col 25 = TF_gene.r, col 26 = TF_gene.p_raw
## Use name-based selection for robustness across GRaNIE versions
sel_cols <- intersect(
  c("TF.name", "gene.name", "TF_gene.r", "TF_gene.p_raw"),
  colnames(GRN_connections.all)
)
if (length(sel_cols) < 4) {
  message("WARNING: Some expected columns not found. Available: ",
          paste(colnames(GRN_connections.all), collapse = ", "))
}
selected_col_GRN <- GRN_connections.all[, sel_cols, drop = FALSE]
write.csv(selected_col_GRN,
          file.path(out_dir, paste0(sample_name, "_GRN_connections_selected.csv")),
          row.names = FALSE, quote = FALSE)
saveRDS(GRN_connections.all,
        file.path(out_dir, paste0(sample_name, "_GRN_connections_all.rds")))

## Filtered, deduplicated, sorted network
filtered_sorted_unique_GRN <- dplyr::filter(selected_col_GRN, TF_gene.p_raw < 0.05) %>%
  dplyr::distinct(TF.name, gene.name, .keep_all = TRUE) %>%
  dplyr::arrange(dplyr::desc(TF_gene.r))

write.csv(filtered_sorted_unique_GRN,
          file.path(out_dir, paste0(sample_name, "_GRN_filtered_unique.csv")),
          row.names = FALSE, quote = FALSE)
message("Filtered unique GRN rows: ", nrow(filtered_sorted_unique_GRN))

## ── Summary stats + eGRN graph ────────────────────────────────
GRN <- generateStatsSummary(
  GRN,
  TF_peak.fdr              = c(0.05, 0.1, 0.2),
  TF_peak.connectionTypes  = "all",
  peak_gene.fdr            = c(0.1, 0.2),
  peak_gene.r_range        = c(0, 1),
  allowMissingGenes        = c(FALSE, TRUE),
  allowMissingTFs          = c(FALSE),
  gene.types               = c("protein_coding", "lincRNA"),
  forceRerun               = TRUE
)
GRN <- plot_stats_connectionSummary(GRN, type = "heatmap", plotAsPDF = FALSE, pages = 3)
GRN <- plot_stats_connectionSummary(GRN, type = "boxplot", plotAsPDF = FALSE, pages = 1)
GRN <- build_eGRN_graph(GRN, forceRerun = TRUE)
saveRDS(GRN, file.path(out_dir, paste0(sample_name, "_GRN_final_eGRN.rds")))

message("=== run_GRaNIE.R complete ===")
message("All outputs written to: ", out_dir)

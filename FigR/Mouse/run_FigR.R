## ============================================================
##  run_FigR.R — unified Human / Mouse FigR script
##  Usage:
##    Rscript run_FigR.R <rna_file> <atac_file> <out_dir> \
##                       <sample_name> <genome>
##
##  Arguments:
##    rna_file    : full path to RNA counts CSV (genes x cells)
##    atac_file   : full path to ATAC counts CSV (peaks x cells)
##    out_dir     : full path to output directory
##    sample_name : label used for output file prefixes (e.g. K562)
##    genome      : "hg38" for Human, "mm10" for Mouse
## ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop(paste(
    "Usage: Rscript run_FigR.R",
    "<rna_file> <atac_file> <out_dir> <sample_name> <genome>"
  ))
}

rna_file    <- args[1]
atac_file   <- args[2]
out_dir     <- args[3]
sample_name <- args[4]
genome      <- args[5]   # "hg38" or "mm10"

## ── Validate inputs ───────────────────────────────────────────
if (!file.exists(rna_file))  stop("RNA file not found: ",  rna_file)
if (!file.exists(atac_file)) stop("ATAC file not found: ", atac_file)
if (!genome %in% c("hg38", "mm10")) {
  stop("genome must be 'hg38' or 'mm10'. Got: ", genome)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("======================================")
message("  FigR run")
message("  RNA file    : ", rna_file)
message("  ATAC file   : ", atac_file)
message("  Out dir     : ", out_dir)
message("  Sample name : ", sample_name)
message("  Genome      : ", genome)
message("======================================")

## ── Libraries ────────────────────────────────────────────────
library(FigR)
library(SummarizedExperiment)
library(GenomicRanges)
library(Matrix)
library(Seurat)
library(Signac)
library(FNN)
library(doParallel)
library(foreach)
library(parallel)
library(dplyr)

## ── Read data ────────────────────────────────────────────────
message("Loading RNA counts...")
rna_counts  <- read.table(rna_file,  sep = ",", row.names = 1, header = TRUE, comment.char = "")

message("Loading ATAC counts...")
atac_counts <- read.table(atac_file, sep = ",", row.names = 1, header = TRUE, comment.char = "")

## Keep only canonical chromosomes — handles both chr1:100-200 and chr1-100-200
atac_counts <- atac_counts[grepl("^chr[0-9XY]", rownames(atac_counts)), ]
message("ATAC peaks after chr filter: ", nrow(atac_counts))

## ── Build ATAC SummarizedExperiment ──────────────────────────
peak_counts_matrix <- as.matrix(atac_counts)

## Auto-detect peak notation: chr1:100-200 (colon) vs chr1-100-200 (all dashes)
first_peak <- rownames(peak_counts_matrix)[1]
if (grepl(":", first_peak)) {
  peak_info <- do.call(rbind, strsplit(rownames(peak_counts_matrix), "[:|-]"))
  message("Peak notation detected: colon-dash (", first_peak, ")")
} else {
  peak_info <- do.call(rbind, strsplit(rownames(peak_counts_matrix), "-"))
  message("Peak notation detected: all-dash (", first_peak, ")")
}

peak_granges <- GRanges(
  seqnames = peak_info[, 1],
  ranges   = IRanges(start = as.numeric(peak_info[, 2]),
                     end   = as.numeric(peak_info[, 3]))
)

atac_se <- SummarizedExperiment(
  assays    = list(counts = peak_counts_matrix),
  rowRanges = peak_granges
)
assay(atac_se, "counts") <- as(assay(atac_se, "counts"), "dgCMatrix")
message("ATAC SummarizedExperiment built: ", nrow(atac_se), " peaks x ",
        ncol(atac_se), " cells")

## ── Build RNA sparse matrix ───────────────────────────────────
rna_sparse <- as(as.matrix(rna_counts), "CsparseMatrix")
message("RNA sparse matrix: ", nrow(rna_sparse), " genes x ",
        ncol(rna_sparse), " cells")

## ── Peak-gene correlations ────────────────────────────────────
message("Running peak-gene correlations (genome = ", genome, ")...")
cisCor <- runGenePeakcorr(
  ATAC.se = atac_se,
  RNAmat  = rna_sparse,
  genome  = genome,        # "hg38" or "mm10" — no hardcoding
  nCores  = 12,
  p.cut   = NULL
)
saveRDS(cisCor, file.path(out_dir, paste0(sample_name, "_cisCor.rds")))
message("cisCor saved.")

cisCor.filt <- cisCor %>% filter(pvalZ <= 0.05)
message("Filtered cisCor rows: ", nrow(cisCor.filt))

## ── DORC genes and scores ─────────────────────────────────────
dorcGenes <- cisCor.filt %>%
  dorcJPlot(cutoff = 1, returnGeneList = TRUE)
message("DORC genes identified: ", length(dorcGenes))

dorcMat <- getDORCScores(
  ATAC.se  = atac_se,
  dorcTab  = cisCor.filt,
  geneList = dorcGenes,
  nCores   = 12
)

## ── WNN kNN graph ─────────────────────────────────────────────
message("Building WNN kNN graph...")
common_cells <- intersect(colnames(rna_sparse), colnames(atac_se))
message("Common cells: ", length(common_cells))

rna_subset  <- rna_sparse[, common_cells]
atac_subset <- atac_se[,   common_cells]

num_genes <- nrow(rna_subset)

seu.multi <- CreateSeuratObject(counts = rna_subset, assay = "RNA")
seu.multi[["ATAC"]] <- CreateAssayObject(counts = assays(atac_subset)$counts)

seu.multi <- NormalizeData(seu.multi)
seu.multi <- FindVariableFeatures(seu.multi,
                                  selection.method = "vst",
                                  nfeatures        = num_genes)
seu.multi <- ScaleData(seu.multi)
seu.multi <- RunPCA(seu.multi, npcs = 50)

seu.multi <- RunTFIDF(seu.multi,      assay = "ATAC")
seu.multi <- FindTopFeatures(seu.multi, min.cutoff = 1, assay = "ATAC")
seu.multi <- RunSVD(seu.multi,        assay = "ATAC", n = 50)

seu.multi <- FindMultiModalNeighbors(
  seu.multi,
  reduction.list       = list("pca", "lsi"),
  dims.list            = list(1:30, 1:30),
  modality.weight.name = "RNA.weight",
  k.nn                 = 30
)

cellKNN.mat <- seu.multi@neighbors$weighted.nn@nn.idx
rownames(cellKNN.mat) <- colnames(dorcMat)

## ── Smooth scores ─────────────────────────────────────────────
message("Smoothing DORC scores...")
dorcMat.smooth <- smoothScoresNN(NNmat = cellKNN.mat, mat = dorcMat,
                                 nCores = 12)
rnaMat.smooth  <- smoothScoresNN(NNmat = cellKNN.mat, mat = rna_sparse,
                                 nCores = 12)

## Standardize DORC smooth matrix
dorcMat.smooth        <- as.matrix(dorcMat.smooth)
col_means             <- colMeans(dorcMat.smooth, na.rm = TRUE)
col_sds               <- apply(dorcMat.smooth, 2, sd, na.rm = TRUE)
col_sds[col_sds == 0] <- 1
dorcMat.smooth <- sweep(dorcMat.smooth, 2, col_means, FUN = "-")
dorcMat.smooth <- sweep(dorcMat.smooth, 2, col_sds,   FUN = "/")
message("NAs in standardized dorcMat.smooth: ", sum(is.na(dorcMat.smooth)))

## ── Run FigR GRN ──────────────────────────────────────────────
message("Running FigR GRN (genome = ", genome, ")...")
nCores <- 12
cl <- makeCluster(nCores, type = "PSOCK")
registerDoParallel(cl)

fig.d <- runFigRGRN(
  ATAC.se   = atac_se,
  rnaMat    = rnaMat.smooth,
  dorcMat   = dorcMat.smooth,
  dorcTab   = cisCor.filt,
  genome    = genome,          # "hg38" or "mm10" — no hardcoding
  dorcGenes = dorcGenes,
  nCores    = nCores
)
stopCluster(cl)

## ── Save outputs ──────────────────────────────────────────────
saveRDS(fig.d, file.path(out_dir, paste0(sample_name, "_fig.d.rds")))

df_filtered <- fig.d[fig.d$Score != 0, ]
out_csv <- file.path(out_dir, paste0(sample_name, "_filtered_network.csv"))
write.csv(df_filtered, out_csv, row.names = FALSE, quote = FALSE)

message("=== FigR complete ===")
message("  Filtered network rows : ", nrow(df_filtered))
message("  Outputs written to    : ", out_dir)

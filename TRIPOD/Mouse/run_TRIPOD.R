## ============================================================
##  run_tripod.R — unified Human / Mouse TRIPOD script
##  Usage:
##    Rscript run_tripod.R <rna_file> <atac_file> <out_dir> \
##                         <sample_name> <genome> <script_dir>
##
##  Arguments:
##    rna_file    : full path to RNA counts CSV
##    atac_file   : full path to ATAC counts CSV
##    out_dir     : full path to output directory
##    sample_name : cell identity label (e.g. K562, mESC_E7.5_rep1)
##    genome      : "hg38" for Human, "mm10" for Mouse
##    script_dir  : directory containing reference files
##                  Human: EnsDb.Hsapiens.v86_annot.rds (optional)
##                  Mouse: EnsDb.Mmusculus.v79_UCSC_mm10.rds
## ============================================================
#renv::activate(project = script_dir)
#renv::restore(project = script_dir, prompt = FALSE)
options(future.globals.maxSize = 2 * 1024^3)
library(BiocParallel)
library(Seurat)
library(Signac)
library(SeuratDisk)
library(GenomeInfoDb)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(DescTools)
library(dendextend)
library(TRIPOD)

## ── Parse arguments ───────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop(paste(
    "Usage: Rscript run_TRIPOD.R",
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

cat("======================================\n")
cat("  TRIPOD run\n")
cat("  RNA file    :", rna_file,    "\n")
cat("  ATAC file   :", atac_file,   "\n")
cat("  Out dir     :", out_dir,     "\n")
cat("  Sample name :", sample_name, "\n")
cat("  Genome      :", genome,      "\n")
cat("  Script dir  :", script_dir,  "\n")
cat("======================================\n")

## ── Validate inputs ───────────────────────────────────────────
if (!file.exists(rna_file))  stop("RNA file not found: ",  rna_file)
if (!file.exists(atac_file)) stop("ATAC file not found: ", atac_file)
if (!genome %in% c("hg38", "mm10")) stop("genome must be 'hg38' or 'mm10'. Got: ", genome)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ── Genome-specific settings ──────────────────────────────────
if (genome == "hg38") {

  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)

  bsgenome       <- BSgenome.Hsapiens.UCSC.hg38
  mt_pattern     <- "^MT-"
  jaspar_species <- 9606                          # Homo sapiens
  std_chroms     <- paste0("chr", 1:22)           # human autosomes
  genome_str     <- "hg38"

  ## Annotations: load from saved RDS if present, otherwise fetch live
  annot_rds <- file.path(script_dir, "EnsDb.Hsapiens.v86_annot.rds")
  if (file.exists(annot_rds)) {
    cat("Loading Human annotations from RDS:", annot_rds, "\n")
    annotations <- readRDS(annot_rds)
  } else {
    cat("Fetching Human annotations from EnsDb.Hsapiens.v86...\n")
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    seqlevelsStyle(annotations) <- "UCSC"
    genome(annotations) <- "hg38"
  }

} else {  # mm10

  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)

  bsgenome       <- BSgenome.Mmusculus.UCSC.mm10
  mt_pattern     <- "^mt-"
  jaspar_species <- 9606                          # intentional: JASPAR human
                                                  # motifs applied to mouse
                                                  # (matches original mouse script)
  std_chroms     <- paste0("chr", 1:19)           # mouse autosomes
  genome_str     <- "mm10"

  ## Mouse annotations always loaded from saved RDS
  annot_rds <- file.path(script_dir, "EnsDb.Mmusculus.v79_UCSC_mm10.rds")
  if (!file.exists(annot_rds)) stop("Mouse annotation RDS not found: ", annot_rds)
  cat("Loading Mouse annotations from RDS:", annot_rds, "\n")
  annotations <- readRDS(annot_rds)

}

## ══════════════════════════════════════════════════════════════
##  SECTION 1: Build Seurat object
## ══════════════════════════════════════════════════════════════

cat("Loading RNA counts...\n")
rna_counts  <- read.csv(rna_file,  row.names = 1)

cat("Loading ATAC counts...\n")
atac_counts <- read.csv(atac_file, row.names = 1)

## ── Create Seurat object ──────────────────────────────────────
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = mt_pattern)

num_cells <- ncol(pbmc)
cat("Cells loaded:", num_cells, "\n")

## ── Add ATAC assay ────────────────────────────────────────────
#grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
###################### Auto-detects separator from the first peak name
first_peak <- rownames(atac_counts)[1]
sep_to_use <- if (grepl(":", first_peak)) c(":", "-") else c("-", "-")
cat("Peak notation detected:", first_peak, "-> using sep:", paste(sep_to_use, collapse=","), "\n")
grange.counts <- StringToGRanges(rownames(atac_counts), sep = sep_to_use)
#######################
grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts   <- atac_counts[as.vector(grange.use), ]

chrom_assay <- CreateChromatinAssay(
  counts     = atac_counts,
  sep        = sep_to_use,
  genome     = genome_str,
  min.cells  = 10,
  annotation = annotations
)
pbmc[["ATAC"]] <- chrom_assay

## ── Assign cell identity from sample_name arg ────────────────
## Replaces all hardcoded cell type labels ("Macrophage", "mESC" etc.)
pbmc$celltype <- sample_name

## ── RNA analysis ──────────────────────────────────────────────
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:50, reduction.name = "umap.rna",
          reduction.key = "rnaUMAP_", verbose = FALSE)

## ── ATAC analysis ─────────────────────────────────────────────
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:50,
                reduction.name = "umap.atac", reduction.key = "atacUMAP_")

## ── WNN integration ───────────────────────────────────────────
pbmc <- FindMultiModalNeighbors(pbmc,
                                reduction.list = list("pca", "lsi"),
                                dims.list      = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn",
                reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

## ── Keep dominant cell type per cluster ──────────────────────
cell.keep <- rep(FALSE, ncol(pbmc))
temp <- table(pbmc$seurat_clusters, pbmc$celltype)
for (i in 1:nrow(temp)) {
  clusteri  <- as.numeric(rownames(temp)[i])
  celltypei <- colnames(temp)[which.max(temp[i, ])]
  cell.keep[which(pbmc$seurat_clusters == clusteri &
                    pbmc$celltype == celltypei)] <- TRUE
}
pbmc <- pbmc[, cell.keep]
cat("Cells after cluster filtering:", ncol(pbmc), "\n")

## ── Renormalize after filtering ───────────────────────────────
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:50, reduction.name = "umap.rna",
          reduction.key = "rnaUMAP_", verbose = FALSE)

DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:50,
                reduction.name = "umap.atac", reduction.key = "atacUMAP_")

pbmc <- FindMultiModalNeighbors(pbmc,
                                reduction.list = list("pca", "lsi"),
                                dims.list      = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn",
                reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

## ══════════════════════════════════════════════════════════════
##  SECTION 2: ChromVAR motif analysis
## ══════════════════════════════════════════════════════════════

DefaultAssay(pbmc) <- "ATAC"

## Get motif PWMs from JASPAR2020
pfm.set <- getMatrixSet(x = JASPAR2020,
                        opts = list(species = jaspar_species,
                                    all_versions = FALSE))

## Mouse: add Olig2 motif (not in JASPAR human set) and fix case
if (genome == "mm10") {
  pfm.olig2 <- PFMatrix(
    ID   = "PFM.Olig2",
    name = "Olig2",
    strand = "+",
    bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
    tags = list(species = "10090", tax_group = "vertebrates"),
    profileMatrix = matrix(
      c(0L, 1L, 0L, 0L, 0L, 0L,
        1L, 0L, 0L, 1L, 0L, 0L,
        0L, 0L, 1L, 0L, 0L, 1L,
        0L, 0L, 0L, 0L, 1L, 0L),
      byrow = TRUE, nrow = 4,
      dimnames = list(c("A", "C", "G", "T")))
  )
  pfm.set <- pfm.set[-which(sapply(pfm.set, function(x) x@name == "OLIG2"))]
  pfm.set$PFM.Olig2 <- pfm.olig2
}

motif.matrix <- CreateMotifMatrix(
  features  = granges(pbmc),
  pwm       = pfm.set,
  genome    = genome_str,
  use.counts = (genome == "hg38")   # TRUE for human, FALSE for mouse
)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pfm.set)
pbmc <- SetAssayData(pbmc, assay = "ATAC", slot = "motifs",
                     new.data = motif.object)

register(MulticoreParam())
#################
#BiocParallel::register(BiocParallel::MulticoreParam())
pbmc <- Signac::RunChromVAR(object = pbmc, genome = bsgenome)

## Save chromVAR results
write.csv(pbmc@assays$chromvar@data,
          file = file.path(out_dir, paste0(sample_name, "_chromVAR_Motif_Results.csv")))
saveRDS(pbmc, file = file.path(out_dir, paste0(sample_name, "_pbmc_after_chromvar.rds")))

## ══════════════════════════════════════════════════════════════
##  SECTION 3: Prepare TRIPOD objects
##  Human: uses TRIPOD's getObjectsForModelFit (automatic)
##  Mouse: manual extraction (matches original mouse script)
## ══════════════════════════════════════════════════════════════

if (genome == "hg38") {

  cat("Preparing TRIPOD objects (Human — automatic via getObjectsForModelFit)...\n")

  tripod.obj  <- getObjectsForModelFit(object = pbmc,
                                       chr = std_chroms)
  transcripts.gr <- tripod.obj$transcripts.gr
  peaks.gr       <- tripod.obj$peaks.gr
  motifxTF       <- tripod.obj$motifxTF
  peakxmotif     <- tripod.obj$peakxmotif

  pbmc <- filterSeuratObject(object = pbmc, tripod.object = tripod.obj)
  pbmc <- processSeuratObject(object = pbmc, dim.rna = 1:50,
                              dim.atac = 2:50, verbose = FALSE)

} else {

  cat("Preparing TRIPOD objects (Mouse — manual extraction)...\n")

  DefaultAssay(pbmc) <- "ATAC"

  ## Collapse to longest transcript, keep protein-coding on std chroms
  transcripts.gr <- Signac:::CollapseToLongestTranscript(
    ranges = Annotation(pbmc))
  transcripts.gr <- transcripts.gr[
    transcripts.gr$gene_biotype == "protein_coding"]
  transcripts.gr <- transcripts.gr[
    seqnames(transcripts.gr) %in% std_chroms]
  transcripts.gr <- sort(transcripts.gr)

  peaks.gr <- pbmc@assays$ATAC@ranges

  motifxTF <- unlist(pbmc@assays$ATAC@motifs@motif.names)
  motifxTF <- cbind(names(motifxTF), motifxTF)
  colnames(motifxTF) <- c("motif", "TF")

  peakxmotif <- pbmc@assays$ATAC@motifs@data

  ## Fix capitalisation so TF names match RNA row names
  motifxTF[, 2] <- stringr::str_to_title(tolower(motifxTF[, 2]))

  ## Keep only TFs present in RNA
  peakxmotif <- peakxmotif[, motifxTF[, 2] %in% rownames(pbmc@assays$RNA)]
  motifxTF   <- motifxTF[  motifxTF[, 2] %in% rownames(pbmc@assays$RNA), ]

  ## Intersect genes with SCT
  genes.common <- intersect(transcripts.gr$gene_name,
                            rownames(pbmc@assays$SCT))
  DefaultAssay(pbmc) <- "RNA"
  pbmc@assays$RNA <- subset(pbmc@assays$RNA,
                            features = match(genes.common,
                                             rownames(pbmc@assays$RNA)))
  pbmc@assays$SCT <- subset(pbmc@assays$SCT,
                            features = match(genes.common,
                                             rownames(pbmc@assays$SCT)))
  transcripts.gr <- transcripts.gr[match(genes.common,
                                         transcripts.gr$gene_name)]
  peakxmotif <- peakxmotif[, motifxTF[, 2] %in% genes.common]
  motifxTF   <- motifxTF[  motifxTF[, 2] %in% genes.common, ]
  pbmc@assays$chromvar <- subset(
    pbmc@assays$chromvar,
    features = match(motifxTF[, 1], rownames(pbmc@assays$chromvar)))

  ## Re-run RNA + ATAC processing after filtering
  DefaultAssay(pbmc) <- "RNA"
  pbmc <- SCTransform(pbmc, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:50, reduction.name = "umap.rna",
            reduction.key = "rnaUMAP_")

  DefaultAssay(pbmc) <- "ATAC"
  pbmc <- RunTFIDF(pbmc)
  pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
  pbmc <- RunSVD(pbmc)
  pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:50,
                  reduction.name = "umap.atac", reduction.key = "atacUMAP_")

  pbmc <- FindMultiModalNeighbors(pbmc,
                                  reduction.list = list("pca", "lsi"),
                                  dims.list      = list(1:50, 2:50))
  pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn",
                  reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
}

## ── Save intermediate TRIPOD objects ─────────────────────────
save(transcripts.gr, file = file.path(out_dir, "transcripts.gr.rda"))
save(peaks.gr,       file = file.path(out_dir, "peaks.gr.rda"))
save(motifxTF,       file = file.path(out_dir, "motifxTF.rda"))
save(peakxmotif,     file = file.path(out_dir, "peakxmotif.rda"))

## ══════════════════════════════════════════════════════════════
##  SECTION 4: Metacell construction
## ══════════════════════════════════════════════════════════════

## Cluster on SCT for metacell construction
DefaultAssay(pbmc) <- "SCT"
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:30, verbose = FALSE)

## Optimise clustering resolution (target ≥20 cells per cluster)
set.seed(123)
num.clusters <- optimizeResolution(
  object     = pbmc,
  graph.name = if (genome == "hg38") "wsnn" else "SCT_snn",
  assay.name = if (genome == "hg38") "WNN"  else "SCT",
  resolutions = seq(10, 30, 5),
  min.num    = 20
)
#default min.num=20, seq=(10,35,5)
cat("Cluster resolution table:\n")
print(num.clusters)

res <- 15
pbmc <- getClusters(
  object     = pbmc,
  graph.name = if (genome == "hg38") "wsnn" else "SCT_snn",
  algorithm  = 3,
  resolution = res,
  verbose    = FALSE
)

metacell.pbmc <- getMetacellMatrices(object       = pbmc,
                                     cluster.name = "seurat_clusters",
                                     min.num      = 20)
pbmc <- removeSmallMetacells(object = pbmc, min.num = 20)

metacell.rna  <- metacell.pbmc$rna
metacell.peak <- metacell.pbmc$peak

pbmc$celltype      <- factor(pbmc$celltype)
metacell.celltype  <- pbmc$celltype

## Highly variable genes
DefaultAssay(pbmc) <- "SCT"
hvg.pbmc <- VariableFeatures(pbmc)

## Save metacell objects
saveRDS(pbmc, file = file.path(out_dir,
                               paste0(sample_name, "_pbmc_metacells.rds")))
save(metacell.rna,      file = file.path(out_dir, "metacell.rna.rda"))
save(metacell.peak,     file = file.path(out_dir, "metacell.peak.rda"))
save(metacell.celltype, file = file.path(out_dir, "metacell.celltype.rda"))
save(hvg.pbmc,          file = file.path(out_dir, "hvg.rda"))

cat("Metacells constructed. Metacell RNA dims:", dim(metacell.rna), "\n")

## ══════════════════════════════════════════════════════════════
##  SECTION 5: TRIPOD model fitting
## ══════════════════════════════════════════════════════════════

genes           <- rownames(pbmc)
ext.upstream    <- ext.downstream <- if (genome == "hg38") 1e5 else 2e5

## Parallelisation: SnowParam for cluster jobs
#BiocParallel::register(BiocParallel::MulticoreParam())
bp <- BiocParallel::SnowParam(workers = 8)

cat("Building XY matrices for", length(genes), "genes...\n")
xymats.list <- bplapply(
  genes,
  getXYMatrices,
  ext.upstream   = ext.upstream,
  transcripts.gr = transcripts.gr,
  peaks.gr       = peaks.gr,
  metacell.rna   = metacell.rna,
  metacell.peak  = metacell.peak,
  peakxmotif     = peakxmotif,
  motifxTF       = motifxTF,
  metacell.celltype = metacell.celltype,
  BPPARAM        = bp
)
names(xymats.list) <- genes

## Filter: keep only genes where Xt is a 2D matrix with >1 column
filtered_xymats <- xymats.list[sapply(xymats.list, function(gd) {
  "Xt" %in% names(gd) && is.matrix(gd$Xt) && ncol(gd$Xt) > 1
})]
genes <- names(filtered_xymats)
cat("Genes passing XY matrix filter:", length(genes), "\n")

## ── Fit TRIPOD models ─────────────────────────────────────────
cat("Fitting TRIPOD model (match.by = Xt)...\n")
xymats.tripod.Xt.list <- bplapply(
  filtered_xymats, fitModel,
  model.name = "TRIPOD", match.by = "Xt", BPPARAM = bp
)
names(xymats.tripod.Xt.list) <- genes

cat("Fitting TRIPOD model (match.by = Yj)...\n")
xymats.tripod.Yj.list <- bplapply(
  filtered_xymats, fitModel,
  model.name = "TRIPOD", match.by = "Yj", BPPARAM = bp
)
names(xymats.tripod.Yj.list) <- genes

## ══════════════════════════════════════════════════════════════
##  SECTION 6: Extract trios and save results
##  8 output files: 2 models × 2 levels × 2 signs
## ══════════════════════════════════════════════════════════════

fdr.thresh <- 1   # keep all hits (filter downstream if needed)

## Helper: extract trios and write CSV
save_trios <- function(xymats.list, model.list, level, sign, tag) {
  df <- getTrios(
    xymats.list = xymats.list,
    fdr.thresh  = fdr.thresh,
    sign        = sign,
    model.name  = model.list,
    level       = level
  )
  fname <- file.path(out_dir, paste0(sample_name, "_", tag, ".csv"))
  write.csv(df, fname, row.names = FALSE)
  cat("Saved:", fname, "— rows:", nrow(df), "\n")
}

## Xt model — positive
save_trios(xymats.tripod.Xt.list, "TRIPOD", 1, "positive", "tX1.pos")
save_trios(xymats.tripod.Xt.list, "TRIPOD", 2, "positive", "tX2.pos")

## Yj model — positive
save_trios(xymats.tripod.Yj.list, "TRIPOD", 1, "positive", "tY1.pos")
save_trios(xymats.tripod.Yj.list, "TRIPOD", 2, "positive", "tY2.pos")

## Xt model — negative
save_trios(xymats.tripod.Xt.list, "TRIPOD", 1, "negative", "tX1.neg")
save_trios(xymats.tripod.Xt.list, "TRIPOD", 2, "negative", "tX2.neg")

## Yj model — negative
save_trios(xymats.tripod.Yj.list, "TRIPOD", 1, "negative", "tY1.neg")
save_trios(xymats.tripod.Yj.list, "TRIPOD", 2, "negative", "tY2.neg")

cat("=== run_tripod.R complete ===\n")
cat("All outputs written to:", out_dir, "\n")

## ══════════════════════════════════════════════════════════════
##  SECTION 7: Aggregate all 8 trio outputs → unique GRN
##
##  Combines all tX1/tX2/tY1/tY2 × pos/neg CSVs, then for each
##  unique gene-TF pair keeps only the entry with the highest
##  absolute coefficient across all models and signs.
##  Output: <sample_name>_unique_GRN.csv
## ══════════════════════════════════════════════════════════════

cat("Aggregating all trio outputs into unique GRN...\n")

all_tags <- c("tX1.pos", "tX2.pos", "tY1.pos", "tY2.pos",
              "tX1.neg", "tX2.neg", "tY1.neg", "tY2.neg")

combined_data <- data.frame()

for (tag in all_tags) {
  fname <- file.path(out_dir, paste0(sample_name, "_", tag, ".csv"))
  if (file.exists(fname)) {
    temp <- read.csv(fname)
    if (nrow(temp) > 0) {
      combined_data <- dplyr::bind_rows(combined_data, temp)
      cat("  Loaded:", basename(fname), "--", nrow(temp), "rows\n")
    } else {
      cat("  Skipped (empty):", basename(fname), "\n")
    }
  } else {
    cat("  WARNING: File not found, skipping:", fname, "\n")
  }
}

if (nrow(combined_data) == 0) {
  cat("  WARNING: No data found across all trio files -- skipping unique GRN output.\n")
} else {

  ## Keep gene, TF, coef columns then deduplicate by max |coef| per pair
  unique_grn <- combined_data %>%
    dplyr::select(gene, TF, coef) %>%
    dplyr::mutate(abs_coef = abs(coef)) %>%
    dplyr::group_by(gene, TF) %>%
    dplyr::slice_max(order_by = abs_coef, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(gene, TF, abs_coef) %>%
    dplyr::arrange(dplyr::desc(abs_coef))

  grn_out <- file.path(out_dir, paste0(sample_name, "_unique_GRN.csv"))
  write.csv(unique_grn, grn_out, row.names = FALSE)
  cat("Unique GRN saved to:", grn_out, "\n")
  cat("  Gene-TF pairs:", nrow(unique_grn), "\n")
  cat("  Unique genes :", length(unique(unique_grn$gene)), "\n")
  cat("  Unique TFs   :", length(unique(unique_grn$TF)), "\n")
}

cat("=== run_tripod.R complete ===\n")
cat("All outputs written to:", out_dir, "\n")

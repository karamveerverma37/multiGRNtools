## ============================================================
##  run_Pando.R — unified Human / Mouse Pando script
##  Usage:
##    Rscript run_Pando.R <rna_file> <atac_file> <out_dir> \
##                        <sample_name> <genome> <script_dir>
##
##  Arguments:
##    rna_file    : full path to RNA counts CSV (genes x cells)
##    atac_file   : full path to ATAC counts CSV (peaks x cells)
##    out_dir     : full path to output directory
##    sample_name : label used for output file prefixes (e.g. K562)
##    genome      : "hg38" for Human, "mm10" for Mouse
##    script_dir  : directory containing reference files
##                  Human: EnsDb.Hsapiens.v86_annot.rds (optional)
##                  Mouse: EnsDb.Mmusculus.v79_UCSC_mm10.rds (required)
##                         cisBP_mouse_pfms_2021.rds        (required)
## ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop(paste(
    "Usage: Rscript run_Pando.R",
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
message("  Pando run")
message("  RNA file    : ", rna_file)
message("  ATAC file   : ", atac_file)
message("  Out dir     : ", out_dir)
message("  Sample name : ", sample_name)
message("  Genome      : ", genome)
message("  Script dir  : ", script_dir)
message("======================================")

## ── Libraries ────────────────────────────────────────────────
library(Seurat)
library(Signac)
library(Matrix)
library(GenomicRanges)
library(Pando)
library(dplyr)

## ── Genome-specific settings ──────────────────────────────────
if (genome == "hg38") {

  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(biovizBase)

  bsgenome   <- BSgenome.Hsapiens.UCSC.hg38
  genome_str <- "hg38"

  ## Annotations: load from saved RDS if present, otherwise fetch live
  annot_rds <- file.path(script_dir, "EnsDb.Hsapiens.v86_annot.rds")
  if (file.exists(annot_rds)) {
    message("Loading Human annotations from RDS: ", annot_rds)
    annotations <- readRDS(annot_rds)
  } else {
    message("Fetching Human annotations from EnsDb.Hsapiens.v86...")
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    seqlevelsStyle(annotations) <- "UCSC"
    genome(annotations) <- "hg38"
  }

  ## Human motifs — built into Pando
  data(motifs)
  pfm        <- motifs
  motif_2_tf <- NULL    # Pando uses its default TF mapping for human

} else {  # mm10

  library(BSgenome.Mmusculus.UCSC.mm10)

  bsgenome   <- BSgenome.Mmusculus.UCSC.mm10
  genome_str <- "mm10"

  ## Mouse annotations — always from saved RDS
  annot_rds <- file.path(script_dir, "EnsDb.Mmusculus.v79_UCSC_mm10.rds")
  if (!file.exists(annot_rds)) {
    stop("Mouse annotation RDS not found: ", annot_rds)
  }
  message("Loading Mouse annotations from RDS: ", annot_rds)
  annotations <- readRDS(annot_rds)

  ## Mouse motifs — CisBP PWM set
  pwm_rds <- file.path(script_dir, "cisBP_mouse_pfms_2021.rds")
  if (!file.exists(pwm_rds)) {
    stop("Mouse CisBP PWM RDS not found: ", pwm_rds)
  }
  message("Loading Mouse CisBP motifs from: ", pwm_rds)
  pfm <- readRDS(pwm_rds)
  motif_2_tf <- data.frame(
    motif  = names(pfm@listData),
    tf     = sapply(pfm@listData, function(x) x@name),
    origin = "CisBP"
  )
}

## ── Load data ────────────────────────────────────────────────
message("Loading RNA counts...")
scRNA_data  <- read.table(rna_file,  header = TRUE, row.names = 1, sep = ",", comment.char = "")

message("Loading ATAC counts...")
scATAC_data <- read.table(atac_file, header = TRUE, row.names = 1, sep = ",", comment.char = "")

## ── Build Seurat object with RNA ─────────────────────────────
seurat_rna <- CreateSeuratObject(counts = scRNA_data, assay = "RNA")
seurat_rna <- NormalizeData(seurat_rna)
message("RNA cells: ", ncol(seurat_rna), " | genes: ", nrow(seurat_rna))

## ── Build ChromatinAssay with ATAC ────────────────────────────
## Auto-detect peak separator: chr1:100-200 vs chr1-100-200
first_peak <- rownames(scATAC_data)[1]
sep_to_use <- if (grepl(":", first_peak)) c(":", "-") else c("-", "-")
message("Peak notation detected: ", first_peak,
        " -> sep: ", paste(sep_to_use, collapse = ","))

grange.counts <- StringToGRanges(rownames(scATAC_data), sep = sep_to_use)
grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac.counts   <- scATAC_data[as.vector(grange.use), ]
message("ATAC peaks after standard chr filter: ", nrow(atac.counts))

chrom.assay <- CreateChromatinAssay(
  counts     = atac.counts,
  sep        = sep_to_use,
  genome     = genome_str,
  min.cells  = 10,
  annotation = annotations
)
seurat_rna[["peaks"]] <- chrom.assay

## ── Preprocess ATAC ───────────────────────────────────────────
seurat_rna <- RunTFIDF(seurat_rna, assay = "peaks")

## ── Initiate GRN ─────────────────────────────────────────────
message("Initiating GRN...")
seurat_object <- initiate_grn(seurat_rna)

## ── Find motifs ───────────────────────────────────────────────
message("Finding motifs (genome = ", genome_str, ")...")
if (is.null(motif_2_tf)) {
  ## Human: use Pando default TF mapping
  seurat_object <- find_motifs(
    seurat_object,
    pfm    = pfm,
    genome = bsgenome
  )
} else {
  ## Mouse: use CisBP motif-to-TF mapping
  seurat_object <- find_motifs(
    seurat_object,
    pfm       = pfm,
    motif_tfs = motif_2_tf,
    genome    = bsgenome
  )
}

## ── Define genes and filter ───────────────────────────────────
genes <- rownames(scRNA_data)

## Mouse: remove predicted/uncharacterised "Rik" genes
if (genome == "mm10") {
  genes_before <- length(genes)
  genes <- genes[!grepl("Rik", genes)]
  message("Mouse Rik gene filter: ", genes_before, " -> ", length(genes), " genes")
}

## ── Infer GRN ────────────────────────────────────────────────
message("Inferring GRN with GLM (", length(genes), " genes)...")
seurat_object <- infer_grn(
  seurat_object,
  genes               = genes,
  peak_to_gene_method = "Signac",
  method              = "glm"
)

## ── Extract and save results ──────────────────────────────────
GRN <- GetGRN(seurat_object)
grn <- GRN@networks$glm_network@coefs

## Raw network
raw_out <- file.path(out_dir, paste0(sample_name, "_raw_network.csv"))
write.csv(grn, raw_out, row.names = FALSE)
message("Raw network saved: ", raw_out, " (", nrow(grn), " rows)")

## Filtered network (pval < 0.05)
filtered_data <- grn %>% dplyr::filter(pval < 0.05)
filt_out <- file.path(out_dir, paste0(sample_name, "_filtered_network.csv"))
write.csv(filtered_data, filt_out, row.names = FALSE)
message("Filtered network saved: ", filt_out, " (", nrow(filtered_data), " rows)")

message("=== run_Pando.R complete ===")
message("All outputs written to: ", out_dir)

## ============================================================
##  run_directnet.R — unified Human / Mouse DIRECTNET script
##  Usage:
##    Rscript run_directnet.R <rna_file> <atac_file> <out_dir> \
##                            <sample_name> <genome> <script_dir>
##
##  Arguments:
##    rna_file    : full path to RNA counts CSV
##    atac_file   : full path to ATAC counts CSV
##    out_dir     : full path to output directory
##    sample_name : label used as cell identity (e.g. K562, mESC_E7.5_rep1)
##    genome      : "hg38" for Human, "mm10" for Mouse
##    script_dir  : directory containing reference files
##                  (TSS regions file, GTF, Mouse annotation RDS)
## ============================================================

library(DIRECTNET)
library(Seurat)
library(Signac)
library(patchwork)
library(dplyr)
library(ggplot2)
options(stringsAsFactors = FALSE)

## ── Parse arguments ───────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop(paste(
    "Usage: Rscript run_directnet.R",
    "<rna_file> <atac_file> <out_dir>",
    "<sample_name> <genome> <script_dir>"
  ))
}

rna_file    <- args[1]   # RNA counts CSV
atac_file   <- args[2]   # ATAC counts CSV
out_dir     <- args[3]   # Output directory
sample_name <- args[4]   # Cell identity label (e.g. K562)
genome      <- args[5]   # "hg38" or "mm10"
script_dir  <- args[6]   # Directory with reference files

cat("======================================\n")
cat("  DIRECTNET run\n")
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

  tss_file   <- file.path(script_dir, "hg38.TSS.regions.txt")
  ensdb      <- EnsDb.Hsapiens.v86
  bsgenome   <- BSgenome.Hsapiens.UCSC.hg38
  species    <- "Homo sapiens"
  gtf_file   <- file.path(script_dir, "Homo_sapiens.GRCh38.113.gtf.gz")

  ## Download GTF if not already present in script_dir
  if (!file.exists(gtf_file)) {
    cat("Downloading human GTF...\n")
    download.file(
      url     = "https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz",
      destfile = gtf_file
    )
  }

  annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg38"

} else {  # mm10

  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)

  tss_file   <- file.path(script_dir, "mm10.TSS.regions1.txt")
  bsgenome   <- BSgenome.Mmusculus.UCSC.mm10
  species    <- "Mus musculus"
  gtf_file   <- file.path(script_dir, "Mus_musculus.GRCm39.113.gtf.gz")

  ## Mouse annotations loaded from saved RDS (avoids EnsDb query overhead)
  annot_rds  <- file.path(script_dir, "EnsDb.Mmusculus.v79_annot.rds")
  if (!file.exists(annot_rds)) stop("Mouse annotation RDS not found: ", annot_rds)
  annotations <- readRDS(annot_rds)
}

## ── Validate reference files ──────────────────────────────────
if (!file.exists(tss_file)) stop("TSS regions file not found: ", tss_file)
if (!file.exists(gtf_file)) stop("GTF file not found: ", gtf_file)

## ── Load genome info ──────────────────────────────────────────
genome.info <- read.table(tss_file, header = TRUE)
unik        <- !duplicated(genome.info$genes)
genome.info <- genome.info[unik, ]
cat("Genome info rows (deduplicated):", nrow(genome.info), "\n")

## ── Load counts ───────────────────────────────────────────────
cat("Loading RNA counts...\n")
rna_counts  <- read.table(rna_file,  header = TRUE, row.names = 1, sep = ",", comment.char = "")

cat("Loading ATAC counts...\n")
atac_counts <- read.table(atac_file, header = TRUE, row.names = 1, sep = ",", comment.char = "")

## ── Build Seurat object ───────────────────────────────────────
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

## ── Add ATAC assay ────────────────────────────────────────────
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts   <- atac_counts[as.vector(grange.use), ]

chrom_assay <- CreateChromatinAssay(
  counts     = atac_counts,
  sep        = c(":", "-"),
  genome     = genome,
  min.cells  = 10,
  annotation = annotations
)
pbmc[["ATAC"]] <- chrom_assay

## ── RNA processing ────────────────────────────────────────────
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

## ── ATAC processing + WNN integration ────────────────────────
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 1)
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50,
                reduction.name = "umap.atac", reduction.key = "atacUMAP_")
pbmc <- FindMultiModalNeighbors(pbmc,
                                reduction.list = list("pca", "lsi"),
                                dims.list      = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn",
                reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

## ── Assign cell identity from sample_name arg ────────────────
## Replaces the hardcoded "Buffer1" in the original scripts
Idents(pbmc)    <- sample_name
pbmc$celltype   <- Idents(pbmc)

DefaultAssay(pbmc) <- "RNA"

##########################
## ── Seurat v5 compatibility fix ───────────────────────────────
## DIRECTNET was written for Seurat v4 and calls @counts directly on
## the RNA assay. In Seurat v5 the assay class changed to "Assay5"
## which uses a different internal structure and breaks that access.
## Downcasting back to the v4 "Assay" class fixes the error:
##   "no slot of name 'counts' for this object of class 'Assay5'"
if (inherits(pbmc[["RNA"]], "Assay5")) {
  cat("Seurat v5 detected — downcasting RNA assay to v4 Assay for DIRECTNET compatibility...\n")
  pbmc[["RNA"]] <- as(pbmc[["RNA"]], "Assay")
}


########################

## ── Run DIRECTNET ─────────────────────────────────────────────
markers <- row.names(pbmc)
cat("Running DIRECTNET on", length(markers), "markers...\n")

pbmc <- Run_DIRECT_NET(
  pbmc,
  peakcalling          = FALSE,
  k_neigh              = 50,
  atacbinary           = TRUE,
  max_overlap          = 0.5,
  size_factor_normalize = FALSE,
  genome.info          = genome.info,
  focus_markers        = markers
)

direct.net_result <- Misc(pbmc, slot = 'direct.net')
direct.net_result <- as.data.frame(do.call(cbind, direct.net_result))

direct.net_result$function_type <- gsub("HF",   "HC", direct.net_result$function_type)
direct.net_result$function_type <- gsub("Rest",  "MC", direct.net_result$function_type)
direct.net_result$function_type <- gsub("LF",   "LC", direct.net_result$function_type)

## ── Load GTF annotation ───────────────────────────────────────
cat("Loading GTF annotation from:", gtf_file, "\n")
gene_anno             <- rtracklayer::readGFF(gtf_file)
gene_anno$chromosome  <- paste0("chr", gene_anno$seqid)
gene_anno$gene        <- gene_anno$gene_id
gene_anno$transcript  <- gene_anno$transcript_id
gene_anno$symbol      <- gene_anno$gene_name

## ── Build focused markers dataframe ───────────────────────────
focused_markers <- data.frame(
  gene  = markers,
  group = rep(sample_name, length(markers))   # uses sample_name, not hardcoded
)

## ── CRE-Gene links ────────────────────────────────────────────
CREs_Gene <- generate_CRE_Gene_links(direct.net_result, markers = focused_markers)

## ── Variable peaks ────────────────────────────────────────────
DefaultAssay(pbmc) <- "ATAC"
variable_peaks <- VariableFeatures(pbmc)

variable_peaks1                <- list()
variable_peaks1[[sample_name]] <- variable_peaks

## ── Focused CREs ──────────────────────────────────────────────
Focused_CREs <- generate_CRE(
  L_G_record   = CREs_Gene$distal,
  P_L_G_record = CREs_Gene$promoter,
  variable_peaks1
)

## ── TF links (distal + promoter) ─────────────────────────────
cat("Detecting TFs for distal CREs...\n")
L_TF_record <- generate_peak_TF_links(
  peaks_bed_list = Focused_CREs$distal,
  species        = species,
  genome         = bsgenome,
  markers        = focused_markers
)

cat("Detecting TFs for promoter CREs...\n")
P_L_TF_record <- generate_peak_TF_links(
  peaks_bed_list = Focused_CREs$promoter,
  species        = species,
  genome         = bsgenome,
  markers        = focused_markers
)

## ── Generate network links ────────────────────────────────────
groups <- pbmc$celltype
network_links <- generate_links_for_Cytoscape(
  L_G_record   = Focused_CREs$L_G_record,
  L_TF_record,
  P_L_G_record = Focused_CREs$P_L_G_record,
  P_L_TF_record,
  groups
)

## ── Save outputs ──────────────────────────────────────────────
outfile <- file.path(out_dir, paste0(sample_name, "_Network_links.csv"))
write.csv(network_links, outfile, row.names = FALSE, quote = FALSE)
cat("Saved network links to:", outfile, "\n")

Node_attribute <- generate_node_for_Cytoscape(network_links, markers = focused_markers)
node_outfile   <- file.path(out_dir, paste0(sample_name, "_Node_attributes.csv"))
write.csv(Node_attribute, node_outfile, row.names = FALSE, quote = FALSE)
cat("Saved node attributes to:", node_outfile, "\n")

cat("=== run_directnet.R complete ===\n")

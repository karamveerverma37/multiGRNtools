library(cicero)
library(monocle3)
library(Matrix)

## ── Parse arguments ───────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript step1_preprocess.R <atac_csv> <out_prefix> <chrom_sizes>")
}

input_file   <- args[1]   # Full path to ATAC CSV
output_file  <- args[2]   # Output prefix (e.g. .../out_dir/atac)
chrom_sizes  <- args[3]   # Full path to chromosome sizes file (hg38.fa.sizes or mm10.chrom.sizes)

cat("Input file    :", input_file,  "\n")
cat("Output prefix :", output_file, "\n")
cat("Chrom sizes   :", chrom_sizes, "\n")

## ── Load ATAC data ────────────────────────────────────────────
indata <- read.csv(input_file, header = TRUE, row.names = 1, sep = ",")
cat("Dimensions:", dim(indata), "\n")
head(indata)

## ── Cell metadata ─────────────────────────────────────────────
cellinfo <- data.frame(cells = colnames(indata))
row.names(cellinfo) <- cellinfo$cells
head(cellinfo)

## ── Peak metadata ─────────────────────────────────────────────
peak_names <- row.names(indata)
peakinfo   <- data.frame(site_name = peak_names)
peakinfo   <- cbind(
  do.call(rbind, strsplit(peak_names, "[:-]")),
  peakinfo
)
names(peakinfo)[1:3] <- c("chr", "bp1", "bp2")
row.names(peakinfo) <- peakinfo$site_name
head(peakinfo)

## ── Build CDS ─────────────────────────────────────────────────
sparse_matrix <- as(as.matrix(indata), "dgCMatrix")

input_cds <- new_cell_data_set(
  expression_data = sparse_matrix,
  cell_metadata   = cellinfo,
  gene_metadata   = peakinfo
)

## ── QC filtering ──────────────────────────────────────────────
input_cds <- detect_genes(input_cds)
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0, ]

max_count <- 50000
min_count <- 2000
input_cds <- input_cds[, Matrix::colSums(exprs(input_cds)) >= min_count]
input_cds <- input_cds[, Matrix::colSums(exprs(input_cds)) <= max_count]

## ── Preprocessing & dimensional reduction ────────────────────
set.seed(2017)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")
input_cds <- reduce_dimension(input_cds,
                              reduction_method  = "UMAP",
                              preprocess_method = "LSI")

umap_coords <- reducedDims(input_cds)$UMAP

## ── Build Cicero CDS ──────────────────────────────────────────
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

## Save Cicero CDS object
cicero_cds_file <- paste0(output_file, "_cicero_cds.Rds")
saveRDS(cicero_cds, cicero_cds_file)
cat("Saved cicero CDS to:", cicero_cds_file, "\n")

## ── Load chromosome lengths (genome-agnostic via arg) ─────────
chromosome_length <- read.table(chrom_sizes)
cat("Chromosome length file loaded:", chrom_sizes, "\n")

## ── Run Cicero ────────────────────────────────────────────────
conns <- run_cicero(cicero_cds, chromosome_length)

## Save connections
connfile <- paste0(output_file, "_cicero_connections.Rds")
saveRDS(conns, connfile)
cat("Saved connections to:", connfile, "\n")

head(conns)

## ── Save peak list and connections as CSV ─────────────────────
peakfile  <- paste0(output_file, "_all_peaks.csv")
all_peaks <- row.names(exprs(input_cds))
all_peaks_cleaned <- gsub("[:\\-]", "_", all_peaks)
write.csv(x = all_peaks_cleaned, file = peakfile, row.names = FALSE)
cat("Saved all peaks to:", peakfile, "\n")

cicero_conn <- paste0(output_file, "_cicero_connections.csv")
write.csv(x = conns, file = cicero_conn, row.names = FALSE)
cat("Saved cicero connections CSV to:", cicero_conn, "\n")

cat("=== step1_preprocess.R complete ===\n")

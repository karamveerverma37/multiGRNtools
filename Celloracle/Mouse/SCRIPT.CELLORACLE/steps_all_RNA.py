## ============================================================
##  steps_all_RNA.py
##  Usage:
##    python steps_all_RNA.py <rna_csv> <out_prefix> <sample_name>
##
##  Arguments:
##    rna_csv     : path to RNA counts CSV (cells x genes or genes x cells)
##    out_prefix  : output file prefix (e.g. .../out_dir/rna)
##                  output written to: <out_prefix>_processed.h5ad
##    sample_name : sample label used as the cell_type annotation
##                  (e.g. K562, Macrophage_S1, mESC — passed from SLURM script)
## ============================================================

import os
import sys
import matplotlib
matplotlib.use('Agg')    # non-interactive backend for cluster
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import celloracle
from scipy import sparse

print(f"Scanpy version      : {sc.__version__}")
print(f"CellOracle version  : {celloracle.__version__}")

plt.rcParams["savefig.dpi"] = 300
plt.rcParams["figure.figsize"] = [6, 4.5]

## ── Parse arguments ───────────────────────────────────────────
if len(sys.argv) < 4:
    sys.exit(
        "Usage: python steps_all_RNA.py <rna_csv> <out_prefix> <sample_name>"
    )

input_file  = sys.argv[1]   # RNA counts CSV
output      = sys.argv[2]   # Output prefix
sample_name = sys.argv[3]   # Sample label used as cell_type (e.g. K562)

print(f"Input file  : {input_file}")
print(f"Output      : {output}")
print(f"Sample name : {sample_name}")

## ── Load RNA counts ───────────────────────────────────────────
counts = pd.read_csv(input_file, index_col=0)
print(f"Counts shape (raw): {counts.shape}")

## ── Build AnnData (cells × genes) ────────────────────────────
adata = ad.AnnData(X=counts.T)   # Transpose: rows=cells, cols=genes
adata.var_names = counts.index   # genes
adata.obs_names = counts.columns # cells

## num_genes must be derived AFTER transpose so it reflects gene count
num_genes = adata.n_vars
print(f"Cells : {adata.n_obs} | Genes : {num_genes}")
print(adata)

## ── Assign cell type from sample name ────────────────────────
## Uses the sample_name argument — not a hardcoded label
ct = np.repeat(sample_name, adata.n_obs)
adata.obs["cell_type"] = pd.Categorical(ct)
print(f"Cell type assigned: '{sample_name}' to all {adata.n_obs} cells")

## ── Gene filtering ────────────────────────────────────────────
sc.pp.filter_genes(adata, min_counts=1)
adata.X = adata.X.astype('float64')

## ── Normalize ─────────────────────────────────────────────────
sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')

## ── Highly variable gene selection ───────────────────────────
filter_result = sc.pp.filter_genes_dispersion(
    adata.X,
    flavor      = 'cell_ranger',
    n_top_genes = num_genes,   # now correctly == number of genes
    log         = False
)

adata = adata[:, filter_result.gene_subset]
sc.pp.normalize_per_cell(adata)

## ── Store raw counts before log transform ────────────────────
adata.raw = adata
adata.layers["raw_count"] = adata.raw.X.copy()
print(adata)

## ── Log transform and scale ───────────────────────────────────
sc.pp.log1p(adata)
sc.pp.scale(adata)

## ── Dimensionality reduction ──────────────────────────────────
sc.tl.pca(adata, svd_solver='arpack')

sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
sc.tl.diffmap(adata)

sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')
sc.tl.louvain(adata, resolution=0.8)

## ── PAGA + force-directed graph layout ───────────────────────
sc.tl.paga(adata, groups='louvain')
sc.pl.paga(adata, show=False)
sc.tl.draw_graph(adata, init_pos='paga', random_state=123)

## ── Save processed AnnData ───────────────────────────────────
out_path = f"{output}_processed.h5ad"
adata.write_h5ad(out_path)
print(f"Saved processed AnnData to: {out_path}")

print("Embedding keys:", list(adata.obsm.keys()))
print(adata)
print("=== steps_all_RNA.py complete ===")

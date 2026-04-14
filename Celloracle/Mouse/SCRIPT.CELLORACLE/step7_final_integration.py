## ============================================================
##  step7_final_integration.py
##  Usage:
##    python step7_final_integration.py <input_RNA> <input_ATAC> <out_dir>
##
##  Arguments:
##    input_RNA  : path to processed RNA .h5ad file  (from steps_all_RNA.py)
##    input_ATAC : path to base GRN .parquet file    (from step3_TSS_annot.py)
##    out_dir    : full absolute path to output directory for this condition
##                 (e.g. .../CellOracle/K562/Original)
##                 All outputs land here — no hardcoded or relative paths.
## ============================================================

import os
import sys
import matplotlib
matplotlib.use('Agg')    # non-interactive backend for cluster
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import celloracle as co

print(f"CellOracle version: {co.__version__}")

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

## ── Parse arguments ───────────────────────────────────────────
if len(sys.argv) < 4:
    sys.exit(
        "Usage: python step7_final_integration.py "
        "<input_RNA> <input_ATAC> <out_dir>"
    )

input_RNA  = sys.argv[1]
input_ATAC = sys.argv[2]
out_dir    = sys.argv[3]   # full absolute path — nothing relative or hardcoded

print(f"=== step7_final_integration ===")
print(f"  input_RNA  : {input_RNA}")
print(f"  input_ATAC : {input_ATAC}")
print(f"  out_dir    : {out_dir}")

## ── Output subdirectories anchored to out_dir ─────────────────
outdir      = os.path.join(out_dir, "output_GRN")
save_folder = os.path.join(out_dir, "figures")
output_stem = "celloracle_network"

os.makedirs(outdir,      exist_ok=True)
os.makedirs(save_folder, exist_ok=True)

## ── Load data ─────────────────────────────────────────────────
adata = sc.read_h5ad(input_RNA)
print(adata)
print(f"Cell number : {adata.shape[0]}")
print(f"Gene number : {adata.shape[1]}")

base_GRN = pd.read_parquet(input_ATAC)
print(base_GRN.head())

## ── Build Oracle object ───────────────────────────────────────
oracle = co.Oracle()

print("Metadata columns       :", list(adata.obs.columns))
print("Dimensional reductions :", list(adata.obsm.keys()))

adata.X = adata.layers["raw_count"].copy()

oracle.import_anndata_as_raw_count(
    adata               = adata,
    cluster_column_name = "cell_type",
    embedding_name      = "X_draw_graph_fr"
)
oracle.import_TF_data(TF_info_matrix=base_GRN)

## ── PCA & kNN imputation ──────────────────────────────────────
oracle.perform_PCA()

plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(
    np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_)) > 0.002)
)[0][0]
plt.axvline(n_comps, c="k")
plt.savefig(os.path.join(save_folder, "pca_variance.png"))
plt.close()

print(f"Auto-selected n_comps : {n_comps}")
n_comps = min(n_comps, 50)

n_cell = oracle.adata.shape[0]
k = int(0.025 * n_cell)
print(f"Cell number : {n_cell}")
print(f"Auto-selected k : {k}")

oracle.knn_imputation(
    n_pca_dims = n_comps,
    k          = k,
    balanced   = True,
    b_sight    = k * 8,
    b_maxl     = k * 4,
    n_jobs     = 4
)

## ── Save & reload oracle ──────────────────────────────────────
oracle_path = os.path.join(outdir, f"{output_stem}.celloracle.oracle")
oracle.to_hdf5(oracle_path)
oracle = co.load_hdf5(oracle_path)

## ── Cell type cluster plot ────────────────────────────────────
sc.pl.draw_graph(oracle.adata, color="cell_type", show=False)
plt.savefig(os.path.join(save_folder, "cell_type_draw_graph.png"))
plt.close()

## ── Get links (GRN per cluster) ───────────────────────────────
links = oracle.get_links(
    cluster_name_for_GRN_unit = "cell_type",
    alpha                     = 10,
    verbose_level             = 10
)

print("Clusters found in links:", list(links.links_dict.keys()))

## ── Filter links ──────────────────────────────────────────────
links.filter_links(p=0.001, weight="coef_abs", threshold_number=None)
links.to_hdf5(
    file_path=os.path.join(outdir, f"{output_stem}_links.celloracle.links")
)

## ── Save GRN per cluster — dynamically, no hardcoded names ───
## Iterates every cluster present in the data regardless of sample type.
for cluster in links.links_dict.keys():
    print(f"  Saving GRN for cluster: {cluster}")

    raw_path = os.path.join(outdir, f"{output_stem}_raw_GRN_for_{cluster}.csv")
    links.links_dict[cluster].to_csv(raw_path)

    if cluster in links.filtered_links:
        filt_path = os.path.join(outdir, f"{output_stem}_{cluster}_final_GRN.csv")
        links.filtered_links[cluster].to_csv(filt_path)
    else:
        print(f"  Warning: no filtered links for cluster '{cluster}'")

print(f"=== Done: all outputs written to {out_dir} ===")

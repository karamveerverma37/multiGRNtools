## ============================================================
##  step3_TSS_annot.py
##  Usage:
##    python step3_TSS_annot.py <all_peaks_csv> <cicero_connections_csv> \
##                              <ref_genome> <out_prefix>
##
##  Arguments:
##    all_peaks_csv          : path to _all_peaks.csv from step1
##    cicero_connections_csv : path to _cicero_connections.csv from step1
##    ref_genome             : reference genome string, e.g. "hg38" or "mm10"
##    out_prefix             : output file prefix (e.g. .../out_dir/atac)
##                             outputs written to:
##                               <out_prefix>.test1.celloracle.tfinfo
##                               <out_prefix>_base_GRN.parquet
##                               <out_prefix>_base_GRN.csv
## ============================================================

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')    # non-interactive backend for cluster
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys, shutil, importlib, glob
from tqdm import tqdm    # cluster-safe — NOT tqdm.notebook
from celloracle import motif_analysis as ma
import celloracle as co

print(f"CellOracle version: {co.__version__}")

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

## ── Parse arguments ───────────────────────────────────────────
if len(sys.argv) < 5:
    sys.exit(
        "Usage: python step3_TSS_annot.py "
        "<all_peaks_csv> <cicero_connections_csv> <ref_genome> <out_prefix>"
    )

arg1        = sys.argv[1]   # all_peaks CSV
arg2        = sys.argv[2]   # cicero connections CSV
ref_genome  = sys.argv[3]   # e.g. "hg38" or "mm10"
out_prefix  = sys.argv[4]   # clean output prefix, e.g. .../out_dir/atac

print(f"all_peaks CSV          : {arg1}")
print(f"cicero connections CSV : {arg2}")
print(f"Reference genome       : {ref_genome}")
print(f"Output prefix          : {out_prefix}")

## ── Load peak list ────────────────────────────────────────────
peaks = pd.read_csv(arg1)
peaks = peaks.x.values
print("Peaks loaded:", peaks[:5])

## ── Load Cicero co-accessibility scores ───────────────────────
cicero_connections = pd.read_csv(arg2)
print(cicero_connections.head())

## ── TSS annotation ────────────────────────────────────────────
tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome=ref_genome)
print(tss_annotated.tail())

## ── Integrate TSS with Cicero connections ─────────────────────
integrated = ma.integrate_tss_peak_with_cicero(
    tss_peak           = tss_annotated,
    cicero_connections = cicero_connections
)
print("Integrated shape:", integrated.shape)
print(integrated.head())

## ── Filter by co-accessibility threshold ─────────────────────
peak = integrated[integrated.coaccess >= 0.8]
peaks = peak[["peak_id", "gene_short_name"]].reset_index(drop=True)
print("Filtered peaks shape:", peaks.shape)
print(peaks.head())

## ── Verify genome installation ────────────────────────────────
genome_installation = ma.is_genome_installed(ref_genome=ref_genome, genomes_dir=None)
print(f"{ref_genome} installation: {genome_installation}")

## ── Check peak format ─────────────────────────────────────────
peaks = ma.check_peak_format(peaks, ref_genome, genomes_dir=None)

## ── Motif scanning ────────────────────────────────────────────
# Instantiate TFinfo object
tfi = ma.TFinfo(
    peak_data_frame = peaks,
    ref_genome      = ref_genome,
    genomes_dir     = None
)

# Scan motifs — may take several hours on large datasets
tfi.scan(
    fpr     = 0.02,
    motifs  = None,   # None loads default motifs
    verbose = True
)

# Save tfinfo object — keyed to out_prefix, not to the CSV path
tfinfo_path = f"{out_prefix}.test1.celloracle.tfinfo"
tfi.to_hdf5(file_path=tfinfo_path)
print(f"Saved tfinfo to: {tfinfo_path}")

print(tfi.scanned_df.head())

## ── Filter and build TF info dataframe ───────────────────────
tfi.reset_filtering()
tfi.filter_motifs_by_score(threshold=10)
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

df = tfi.to_dataframe()
print(df.head())

## ── Save base GRN ─────────────────────────────────────────────
parquet_path = f"{out_prefix}_base_GRN.parquet"
csv_path     = f"{out_prefix}_base_GRN.csv"

df.to_parquet(parquet_path)
df.to_csv(csv_path)

print(f"Saved base GRN parquet : {parquet_path}")
print(f"Saved base GRN CSV     : {csv_path}")
print("=== step3_TSS_annot.py complete ===")

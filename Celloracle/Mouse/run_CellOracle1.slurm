#!/bin/bash -l
#SBATCH --job-name=CellOracle_array
#SBATCH --output=/gpfs/Labs/Uzun/RESULTS/PROJECTS/2024.GRN_BENCHMARKING.KARAMVEER/Generate_Shuffling/CellOracle/logs/CellOracle_%A_%a.%N.txt
#SBATCH --error=/gpfs/Labs/Uzun/RESULTS/PROJECTS/2024.GRN_BENCHMARKING.KARAMVEER/Generate_Shuffling/CellOracle/logs/CellOracle_%A_%a.%N.err
#SBATCH -p memory
#SBATCH -A lab_uzun
#SBATCH --mail-user=kkaramveer@pennstatehealth.psu.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=16G
#SBATCH --array=0-999        # Upper bound; guarded below against actual sample count

## ══════════════════════════════════════════════════════════════
##  USER CONFIGURATION
##  Set GENOME to one of: "Human" | "Mouse"
##
##  The reference genome (hg38 / mm10) and script directory are
##  determined entirely by this one setting — including for any
##  Mixed samples nested inside the same folder.
##
##  Expected directory layout:
##
##  BASE_DIR/
##  ├── Human/                            (GENOME="Human" → hg38)
##  │   ├── K562/
##  │   │   ├── subsampled_RNA_*_filtered.csv
##  │   │   ├── subsampled_ATAC_*_filtered.csv
##  │   │   ├── Cell_shuffled/
##  │   │   │   ├── shuffled_RNA_filtered.csv
##  │   │   │   └── original_ATAC_filtered.csv
##  │   │   ├── Gene_shuffled/
##  │   │   │   ├── shuffled_RNA_filtered.csv
##  │   │   │   └── original_ATAC_filtered.csv
##  │   │   └── ATAC_gene_shuffled/
##  │   │       ├── shuffled_RNA_filtered.csv
##  │   │       └── shuffled_ATAC_filtered.csv
##  │   ├── Macrophage_S1/ ...
##  │   ├── iPSC/          ...
##  │   └── Mixed/                        (auto-detected; still hg38)
##  │       ├── K562RNA_MacS1ATAC/
##  │       │   ├── matched_*_rna_counts.csv
##  │       │   └── matched_*_atac_counts.csv
##  │       └── MacS1RNA_K562ATAC/ ...
##  │
##  └── Mouse/                            (GENOME="Mouse" → mm10)
##      ├── mESC_E7.5_rep1/ ...
##      ├── mESC_E8.5_rep2/ ...
##      └── Mixed/                        (auto-detected; still mm10)
##          ├── mESC_RNA_mHippoATAC/
##          │   ├── matched_*_rna_counts.csv
##          │   └── matched_*_atac_counts.csv
##          └── ...
## ══════════════════════════════════════════════════════════════
GENOME="Human"    # <-- change to "Mouse" as needed

## ── Conda + modules ───────────────────────────────────────────
source /gpfs/Home/kmk7420/.conda/pkgs/conda-24.7.1-py39hf3d152e_0/lib/python3.9/site-packages/conda/shell/etc/profile.d/conda.sh
conda activate celloracle_env
source /etc/profile.d/modules.sh
module load R/4.3.2
module load bedtools

## ── Base paths ────────────────────────────────────────────────
BASE_DIR="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.KARAMVEER/INPUT.DATA.MULTIOMICS"
RESULTS_DIR="/gpfs/Labs/Uzun/RESULTS/PROJECTS/2024.GRN_BENCHMARKING.KARAMVEER/Generate_Shuffling/CellOracle"
LOG_DIR="${RESULTS_DIR}/logs"

## ── Genome-specific settings ──────────────────────────────────
if [[ "$GENOME" == "Human" ]]; then
  SCRIPT_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.KARAMVEER/Tools_pipelines/Celloracle/Mouse/SCRIPT.CELLORACLE"
  GENOME_REF="hg38"
  CHROM_SIZES="${SCRIPT_DIR}/hg38.fa.sizes"
  SEARCH_DIR="${BASE_DIR}/Human"

elif [[ "$GENOME" == "Mouse" ]]; then
  SCRIPT_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.KARAMVEER/Tools_pipelines/Celloracle/Mouse/SCRIPT.CELLORACLE"
  GENOME_REF="mm10"
  CHROM_SIZES="${SCRIPT_DIR}/mm10.chrom.sizes"
  SEARCH_DIR="${BASE_DIR}/Mouse"

else
  echo "ERROR: GENOME must be 'Human' or 'Mouse'. Got: '$GENOME'"
  exit 1
fi

mkdir -p "$LOG_DIR"

## ── Validate required paths ───────────────────────────────────
if [[ ! -d "$SEARCH_DIR" ]]; then
  echo "ERROR: Search directory not found: $SEARCH_DIR"
  exit 1
fi
if [[ ! -f "$CHROM_SIZES" ]]; then
  echo "ERROR: Chromosome sizes file not found: $CHROM_SIZES"
  exit 1
fi

## ══════════════════════════════════════════════════════════════
##  File-finder helpers
## ══════════════════════════════════════════════════════════════

## Regular samples — top-level original files
find_rna()           { ls "$1"/subsampled_RNA*_filtered.csv   2>/dev/null | head -1; }
find_atac()          { ls "$1"/subsampled_ATAC*_filtered.csv  2>/dev/null | head -1; }

## Regular samples — condition subdir files
find_rna_shuffled()  { ls "$1"/shuffled_RNA*_filtered.csv     2>/dev/null | head -1; }
find_atac_shuffled() { ls "$1"/shuffled_ATAC*_filtered.csv    2>/dev/null | head -1; }
find_atac_original() { ls "$1"/original_ATAC*_filtered.csv    2>/dev/null | head -1; }

## Mixed samples — matched_* files
find_rna_mixed()     { ls "$1"/matched_*_rna_counts.csv       2>/dev/null | head -1; }
find_atac_mixed()    { ls "$1"/matched_*_atac_counts.csv      2>/dev/null | head -1; }

## Returns 0 (true) if the name is a known tool output dir — never treat as data
is_tool_output_dir() {
  case "$1" in
    DIRECTNET_output|FigR_output|Pando_output|CellOracle_output) return 0 ;;
    *) return 1 ;;
  esac
}

## ══════════════════════════════════════════════════════════════
##  Build sample list
##
##  Each entry is encoded as "MODE|FULL_PATH" so a single flat
##  array can hold both regular and mixed samples while preserving
##  which file-finding / dispatch logic applies to each.
##
##  Modes:
##    regular : subsampled_* top-level + shuffled condition subdirs
##    mixed   : matched_* files; the subdir itself is one sample
##
##  Discovery order for SEARCH_DIR (e.g. Human/):
##    1. Immediate subdirs that are NOT named Mixed → regular mode
##    2. Mixed/ subdir → each of its children      → mixed mode
## ══════════════════════════════════════════════════════════════
mapfile -t SAMPLES < <(
  ## ── 1. Regular samples ──────────────────────────────────────
  for d in "$SEARCH_DIR"/*/; do
    [[ -d "$d" ]] || continue
    sample=$(basename "$d")
    is_tool_output_dir "$sample" && continue
    [[ "$sample" == "Mixed" ]] && continue    # handled separately below

    has_data=0

    ## Check top-level original files
    if [[ -n "$(find_rna "$d")" && -n "$(find_atac "$d")" ]]; then
      has_data=1
    fi

    ## Check condition subdirectories
    if [[ $has_data -eq 0 ]]; then
      for sub in "$d"*/; do
        [[ -d "$sub" ]] || continue
        is_tool_output_dir "$(basename "$sub")" && continue

        rna=$(find_rna_shuffled "$sub")
        [[ -z "$rna"  ]] && rna=$(find_rna "$sub")

        atac=$(find_atac_shuffled "$sub")
        [[ -z "$atac" ]] && atac=$(find_atac_original "$sub")
        [[ -z "$atac" ]] && atac=$(find_atac "$sub")

        if [[ -n "$rna" && -n "$atac" ]]; then
          has_data=1; break
        fi
      done
    fi

    [[ $has_data -eq 1 ]] && echo "regular|${d}"
  done

  ## ── 2. Mixed samples ────────────────────────────────────────
  MIXED_DIR="${SEARCH_DIR}/Mixed"
  if [[ -d "$MIXED_DIR" ]]; then
    for d in "$MIXED_DIR"/*/; do
      [[ -d "$d" ]] || continue
      rna=$(find_rna_mixed  "$d")
      atac=$(find_atac_mixed "$d")
      if [[ -n "$rna" && -n "$atac" ]]; then
        echo "mixed|${d}"
      else
        echo "  INFO: Skipping Mixed subdir $(basename "$d") — no matched RNA+ATAC pair" >&2
      fi
    done
  fi
)

echo "========================================"
echo "  Genome       : $GENOME"
echo "  Genome ref   : $GENOME_REF"
echo "  Chrom sizes  : $CHROM_SIZES"
echo "  Script dir   : $SCRIPT_DIR"
echo "  Search dir   : $SEARCH_DIR"
echo "  Samples found: ${#SAMPLES[@]}"
for entry in "${SAMPLES[@]}"; do
  mode="${entry%%|*}"
  path="${entry##*|}"
  echo "    [${mode}] $(basename "$path")"
done
echo "========================================"

## ── Guard against out-of-range array index ────────────────────
if [[ ${SLURM_ARRAY_TASK_ID} -ge ${#SAMPLES[@]} ]]; then
  echo "Array index ${SLURM_ARRAY_TASK_ID} >= ${#SAMPLES[@]}. Nothing to do. Exiting."
  exit 0
fi

## ── Decode the MODE|PATH entry for this array task ────────────
ENTRY="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
SAMPLE_MODE="${ENTRY%%|*}"       # everything before the first |
SAMPLE_DIR="${ENTRY##*|}"        # everything after  the first |
SAMPLE=$(basename "$SAMPLE_DIR")

echo "========================================"
echo "  Array index  : $SLURM_ARRAY_TASK_ID"
echo "  Sample mode  : $SAMPLE_MODE"
echo "  Sample       : $SAMPLE"
echo "  Sample dir   : $SAMPLE_DIR"
echo "  Node         : $SLURMD_NODENAME"
echo "  Start        : $(date)"
echo "========================================"

## ══════════════════════════════════════════════════════════════
##  run_condition — execute the full 4-step CellOracle pipeline
##
##  Arguments:
##    $1 label      : condition name (Original, Cell_shuffled, …)
##    $2 rna_file   : full path to RNA CSV
##    $3 atac_file  : full path to ATAC CSV
##    $4 out_dir    : full output directory for this condition
##    $5 cell_label : label written as cell_type in AnnData
##                    (optional; defaults to $SAMPLE)
## ══════════════════════════════════════════════════════════════
run_condition() {
  local label="$1"
  local rna_file="$2"
  local atac_file="$3"
  local out_dir="$4"
  local cell_label="${5:-$SAMPLE}"

  echo ""
  echo "  -- $label --"
  echo "     RNA        : $rna_file"
  echo "     ATAC       : $atac_file"
  echo "     Cell label : $cell_label"
  echo "     Out dir    : $out_dir"

  if [[ ! -f "$rna_file"  ]]; then echo "  SKIP: RNA file not found : $rna_file";  return; fi
  if [[ ! -f "$atac_file" ]]; then echo "  SKIP: ATAC file not found: $atac_file"; return; fi

  mkdir -p "$out_dir"

  ## Intermediate file paths — all anchored to out_dir
  local out_atac_prefix="${out_dir}/atac"
  local out_rna_prefix="${out_dir}/rna"
  local all_peak="${out_atac_prefix}_all_peaks.csv"
  local cicero_connection="${out_atac_prefix}_cicero_connections.csv"
  local int_atac_input="${out_atac_prefix}_base_GRN.parquet"
  local int_rna_input="${out_rna_prefix}_processed.h5ad"

  {
    echo "=== CellOracle: $label ==="
    echo "  RNA        : $rna_file"
    echo "  ATAC       : $atac_file"
    echo "  Cell label : $cell_label"
    echo "  Out dir    : $out_dir"
    echo "  Genome ref : $GENOME_REF"
    echo "  Start      : $(date)"

    ## [1/4] Preprocess ATAC → all_peaks.csv + cicero_connections.csv
    echo ""
    echo "  [1/4] step1_preprocess.R ..."
    srun /usr/bin/time -v Rscript "${SCRIPT_DIR}/step1_preprocess.R" \
      "$atac_file" \
      "$out_atac_prefix" \
      "$CHROM_SIZES" \
      || { echo "  FAILED at step1_preprocess.R"; exit 1; }

    ## [2/4] TSS annotation + motif scanning → base GRN parquet
    echo ""
    echo "  [2/4] step3_TSS_annot.py ..."
    srun /usr/bin/time -v python "${SCRIPT_DIR}/step3_TSS_annot.py" \
      "$all_peak" \
      "$cicero_connection" \
      "$GENOME_REF" \
      "$out_atac_prefix" \
      || { echo "  FAILED at step3_TSS_annot.py"; exit 1; }

    ## [3/4] Process RNA → processed .h5ad
    echo ""
    echo "  [3/4] steps_all_RNA.py ..."
    srun /usr/bin/time -v python "${SCRIPT_DIR}/steps_all_RNA.py" \
      "$rna_file" \
      "$out_rna_prefix" \
      "$cell_label" \
      || { echo "  FAILED at steps_all_RNA.py"; exit 1; }

    ## [4/4] Integrate RNA + ATAC → GRN per cluster
    echo ""
    echo "  [4/4] step7_final_integration.py ..."
    srun /usr/bin/time -v python "${SCRIPT_DIR}/step7_final_integration.py" \
      "$int_rna_input" \
      "$int_atac_input" \
      "$out_dir" \
      || { echo "  FAILED at step7_final_integration.py"; exit 1; }

    echo ""
    echo "  Done: $label — $(date)"

  } > "$out_dir/run.log" 2>&1 \
    && echo "  Done   : $label" \
    || echo "  FAILED : $label — see $out_dir/run.log"
}

## ══════════════════════════════════════════════════════════════
##  Dispatch — behaviour differs by SAMPLE_MODE
## ══════════════════════════════════════════════════════════════

if [[ "$SAMPLE_MODE" == "regular" ]]; then

  ## ────────────────────────────────────────────────────────────
  ##  Regular dispatch (Human / Mouse samples)
  ##
  ##  Condition: Original
  ##    Source : top-level subsampled_* files
  ##    Output : RESULTS_DIR / GENOME / SAMPLE / Original /
  ##
  ##  Conditions: one per data subdir (Cell_shuffled, etc.)
  ##    RNA    : shuffled_RNA_filtered.csv
  ##    ATAC   : shuffled_ATAC_filtered.csv   (ATAC_gene_shuffled)
  ##           | original_ATAC_filtered.csv   (Cell_shuffled, Gene_shuffled)
  ##           | subsampled_ATAC_filtered.csv (fallback)
  ##    Output : RESULTS_DIR / GENOME / SAMPLE / <subdir_name> /
  ## ────────────────────────────────────────────────────────────

  RNA_TOP=$(find_rna  "$SAMPLE_DIR")
  ATAC_TOP=$(find_atac "$SAMPLE_DIR")

  if [[ -n "$RNA_TOP" && -n "$ATAC_TOP" ]]; then
    run_condition "Original" \
      "$RNA_TOP" \
      "$ATAC_TOP" \
      "$RESULTS_DIR/$GENOME/$SAMPLE/Original" \
      "$SAMPLE"
  else
    echo "  INFO: No top-level RNA/ATAC found for $SAMPLE — skipping Original"
  fi

  for SUB_DIR in "$SAMPLE_DIR"*/; do
    [[ -d "$SUB_DIR" ]] || continue
    CONDITION=$(basename "$SUB_DIR")
    is_tool_output_dir "$CONDITION" && { echo "  INFO: Skipping tool output dir: $CONDITION"; continue; }

    RNA_SUB=$(find_rna_shuffled "$SUB_DIR")
    [[ -z "$RNA_SUB"  ]] && RNA_SUB=$(find_rna "$SUB_DIR")

    ATAC_SUB=$(find_atac_shuffled "$SUB_DIR")
    [[ -z "$ATAC_SUB" ]] && ATAC_SUB=$(find_atac_original "$SUB_DIR")
    [[ -z "$ATAC_SUB" ]] && ATAC_SUB=$(find_atac "$SUB_DIR")

    if [[ -z "$RNA_SUB" || -z "$ATAC_SUB" ]]; then
      echo "  INFO: Subdir '$CONDITION' has no complete RNA+ATAC pair — skipping"
      continue
    fi

    run_condition "$CONDITION" \
      "$RNA_SUB" \
      "$ATAC_SUB" \
      "$RESULTS_DIR/$GENOME/$SAMPLE/$CONDITION" \
      "$SAMPLE"
  done

elif [[ "$SAMPLE_MODE" == "mixed" ]]; then

  ## ────────────────────────────────────────────────────────────
  ##  Mixed dispatch
  ##
  ##  Each subdir under Mixed/ is its own sample with a single
  ##  "Original" condition — no shuffled subdirs.
  ##  Cell label = subdir name (e.g. K562RNA_MacS1ATAC).
  ##  Output: RESULTS_DIR / GENOME / Mixed / SAMPLE / Original /
  ## ────────────────────────────────────────────────────────────

  RNA_MIXED=$(find_rna_mixed  "$SAMPLE_DIR")
  ATAC_MIXED=$(find_atac_mixed "$SAMPLE_DIR")

  if [[ -n "$RNA_MIXED" && -n "$ATAC_MIXED" ]]; then
    run_condition "Original" \
      "$RNA_MIXED" \
      "$ATAC_MIXED" \
      "$RESULTS_DIR/$GENOME/Mixed/$SAMPLE/Original" \
      "$SAMPLE"
  else
    echo "  ERROR: No matched RNA+ATAC pair found in $SAMPLE_DIR"
    exit 1
  fi

fi

echo ""
echo "========================================"
echo "  Sample $SAMPLE complete : $(date)"
echo "========================================"

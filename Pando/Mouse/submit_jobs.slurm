#!/bin/bash -l
#SBATCH --job-name=Pando_array
#SBATCH --output=/gpfs/Labs/Uzun/RESULTS/PROJECTS/2024.GRN_BENCHMARKING.KARAMVEER/Generate_Shuffling/Pando/logs/Pando_%A_%a.%N.txt
#SBATCH --error=/gpfs/Labs/Uzun/RESULTS/PROJECTS/2024.GRN_BENCHMARKING.KARAMVEER/Generate_Shuffling/Pando/logs/Pando_%A_%a.%N.err
#SBATCH -p memory
#SBATCH -A lab_x
#SBATCH --mail-user=kk@x.psu.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=16G
#SBATCH --time=12:00:00
#SBATCH --array=0-999        # Upper bound; guarded below against actual sample count

## ══════════════════════════════════════════════════════════════
##  USER CONFIGURATION
##  Set GENOME to one of: "Human" | "Mouse"
##
##  Directory layout under BASE_DIR:
##
##  BASE_DIR/
##  ├── Human/
##  │   ├── K562/
##  │   │   ├── subsampled_RNA_*_filtered.csv
##  │   │   ├── subsampled_ATAC_*_filtered.csv
##  │   │   ├── Cell_shuffled/
##  │   │   │   ├── shuffled_RNA_filtered.csv
##  │   │   │   └── original_ATAC_filtered.csv
##  │   │   ├── Gene_shuffled/  ...
##  │   │   └── ATAC_gene_shuffled/  ...
##  │   ├── Macrophage_S1/ ...
##  │   └── Mixed/
##  │       ├── K562RNA_MacS1ATAC/
##  │       │   ├── matched_*_rna_counts.csv
##  │       │   └── matched_*_atac_counts.csv
##  │       └── ...
##  └── Mouse/
##      ├── mESC_E7.5_rep1/ ...
##      └── Mixed/
##          ├── mESC_RNA_mHippoATAC/ ...
##          └── ...
## ══════════════════════════════════════════════════════════════
GENOME="Mouse"    # <-- change to "Mouse" as needed

## ── Conda + modules ───────────────────────────────────────────
source /gpfs/Home/kmk7420/.conda/pkgs/conda-24.7.1-py39hf3d152e_0/lib/python3.9/site-packages/conda/shell/etc/profile.d/conda.sh
conda activate celloracle_env
source /etc/profile.d/modules.sh
module load R/4.3.2

## ── Base paths ────────────────────────────────────────────────
BASE_DIR="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.KARAMVEER/INPUT.DATA.MULTIOMICS"
RESULTS_DIR="/gpfs/Labs/Uzun/RESULTS/PROJECTS/2024.GRN_BENCHMARKING.KARAMVEER/Generate_Shuffling/Pando"
LOG_DIR="${RESULTS_DIR}/logs"

## ── Genome-specific settings ──────────────────────────────────
if [[ "$GENOME" == "Human" ]]; then
  SCRIPT_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.KARAMVEER/Tools_pipelines/Pando/Mouse"
  GENOME_REF="hg38"
  SEARCH_DIR="${BASE_DIR}/Human"

elif [[ "$GENOME" == "Mouse" ]]; then
  SCRIPT_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.KARAMVEER/Tools_pipelines/Pando/Mouse"
  GENOME_REF="mm10"
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

## ══════════════════════════════════════════════════════════════
##  File-finder helpers
## ══════════════════════════════════════════════════════════════

find_rna()           { ls "$1"/subsampled_RNA*_filtered.csv   2>/dev/null | head -1; }
find_atac()          { ls "$1"/subsampled_ATAC*_filtered.csv  2>/dev/null | head -1; }
find_rna_shuffled()  { ls "$1"/shuffled_RNA*_filtered.csv     2>/dev/null | head -1; }
find_atac_shuffled() { ls "$1"/shuffled_ATAC*_filtered.csv    2>/dev/null | head -1; }
find_atac_original() { ls "$1"/original_ATAC*_filtered.csv    2>/dev/null | head -1; }
find_rna_mixed()     { ls "$1"/matched_*_rna_counts.csv       2>/dev/null | head -1; }
find_atac_mixed()    { ls "$1"/matched_*_atac_counts.csv      2>/dev/null | head -1; }

is_tool_output_dir() {
  case "$1" in
    DIRECTNET_output|FigR_output|Pando_output|CellOracle_output|TRIPOD_output) return 0 ;;
    *) return 1 ;;
  esac
}

## ══════════════════════════════════════════════════════════════
##  Build sample list — encoded as "MODE|FULL_PATH"
## ══════════════════════════════════════════════════════════════
mapfile -t SAMPLES < <(

  ## ── 1. Regular samples ──────────────────────────────────────
  for d in "$SEARCH_DIR"/*/; do
    [[ -d "$d" ]] || continue
    sample=$(basename "$d")
    is_tool_output_dir "$sample" && continue
    [[ "$sample" == "Mixed" ]] && continue

    has_data=0

    if [[ -n "$(find_rna "$d")" && -n "$(find_atac "$d")" ]]; then
      has_data=1
    fi

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

## ── Decode MODE|PATH for this array task ──────────────────────
ENTRY="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
SAMPLE_MODE="${ENTRY%%|*}"
SAMPLE_DIR="${ENTRY##*|}"
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
##  run_condition — run Pando for one RNA+ATAC pair
##
##  Arguments:
##    $1 label      : condition name (Original, Cell_shuffled, …)
##    $2 rna_file   : full path to RNA CSV
##    $3 atac_file  : full path to ATAC CSV
##    $4 out_dir    : full output directory for this condition
##    $5 cell_label : sample label for output prefixes
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

  {
    echo "=== Pando: $label ==="
    echo "  RNA        : $rna_file"
    echo "  ATAC       : $atac_file"
    echo "  Cell label : $cell_label"
    echo "  Out dir    : $out_dir"
    echo "  Genome ref : $GENOME_REF"
    echo "  Start      : $(date)"

    srun /usr/bin/time -v Rscript "${SCRIPT_DIR}/run_Pando.R" \
      "$rna_file" \
      "$atac_file" \
      "$out_dir" \
      "$cell_label" \
      "$GENOME_REF" \
      "$SCRIPT_DIR" \
      || { echo "  FAILED at run_Pando.R"; exit 1; }

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

  ## ── Condition: Original ───────────────────────────────────────
  RNA_TOP=$(find_rna  "$SAMPLE_DIR")
  ATAC_TOP=$(find_atac "$SAMPLE_DIR")

  if [[ -n "$RNA_TOP" && -n "$ATAC_TOP" ]]; then
    run_condition "Original" \
      "$RNA_TOP" "$ATAC_TOP" \
      "$RESULTS_DIR/$GENOME/$SAMPLE/Original" \
      "$SAMPLE"
  else
    echo "  INFO: No top-level RNA/ATAC found for $SAMPLE — skipping Original"
  fi

  ## ── Conditions: one per data subdir ──────────────────────────
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
      "$RNA_SUB" "$ATAC_SUB" \
      "$RESULTS_DIR/$GENOME/$SAMPLE/$CONDITION" \
      "$SAMPLE"
  done

elif [[ "$SAMPLE_MODE" == "mixed" ]]; then

  ## ── Mixed: single Original condition per sample ───────────────
  RNA_MIXED=$(find_rna_mixed  "$SAMPLE_DIR")
  ATAC_MIXED=$(find_atac_mixed "$SAMPLE_DIR")

  if [[ -n "$RNA_MIXED" && -n "$ATAC_MIXED" ]]; then
    run_condition "Original" \
      "$RNA_MIXED" "$ATAC_MIXED" \
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

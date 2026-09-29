#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# MiFate full pipeline script
# Supports: -i input MAGs, -o output dir, -d KOfam DB
# Optional: --delta, --min_active
# -----------------------------

# Default parameters
DELTA=0.05
MIN_ACTIVE=1

# Help message
if [[ $# -eq 0 || "$1" == "-h" || "$1" == "--help" ]]; then
    cat <<EOF
Usage: bash run_mifate_full.sh -i MAG_DIR -o OUT_DIR -d KOFAM_DB [options]

Required arguments:
  -i            Input MAG directory (each MAG as .fna file)
  -o            Output directory
  -d            KOfam database directory (contains profiles/ and ko_list)

Optional arguments:
  --delta FLOAT        Gray-zone half-width around threshold (default: 0.05)
  --min_active INT     Minimum number of active features per MAG (default: 1)
  -h, --help           Show this help message and exit

Example:
  bash run_mifate_full.sh -i MAGs -o output -d /data/kofam_db --delta 0.05 --min_active 3
EOF
    exit 0
fi

# -----------------------------
# Parse arguments
# -----------------------------
MAG_DIR=""
OUT_DIR=""
KOFAM_DB=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -i) MAG_DIR="$2"; shift 2 ;;
        -o) OUT_DIR="$2"; shift 2 ;;
        -d) KOFAM_DB="$2"; shift 2 ;;
        --delta) DELTA="$2"; shift 2 ;;
        --min_active) MIN_ACTIVE="$2"; shift 2 ;;
        -h|--help) shift ;;  # Already handled
        *) echo "[WARN] Unknown option: $1"; shift ;;
    esac
done

# Check required parameters
if [[ -z "$MAG_DIR" || -z "$OUT_DIR" || -z "$KOFAM_DB" ]]; then
    echo "[ERROR] Must provide -i input MAG folder, -o output folder, -d KOfam DB"
    exit 1
fi

mkdir -p "$OUT_DIR/function/01_proteins"
mkdir -p "$OUT_DIR/function/02_kofam"

BIN_DIR="$(dirname "$0")"

# Step 1: Prodigal + KOfam
"$BIN_DIR/run_prodigal_kofam.sh" -i "$MAG_DIR" -o "$OUT_DIR/function" -d "$KOFAM_DB"

# Step 2: Build binary KO matrix
KOFAM_MAPPER="$OUT_DIR/function/02_kofam/kofam.mapper"
MAG_LIST="$OUT_DIR/mag_list.txt"
awk '{split($1,a,"|"); print a[1]}' "$KOFAM_MAPPER" | sort | uniq > "$MAG_LIST"

python3 "$BIN_DIR/build_binary_ko_matrix.py" -m "$KOFAM_MAPPER" -l "$MAG_LIST" -o "$OUT_DIR/MAG_KEGG_binary_matrix.csv"

# Step 3: Prediction
python3 "$BIN_DIR/predict.py" \
    --infile "$OUT_DIR/MAG_KEGG_binary_matrix.csv" \
    --sec_model "$BIN_DIR/../models/SEC_model.xgb" \
    --dis_model "$BIN_DIR/../models/DIS_model.xgb" \
    --sec_features "$BIN_DIR/../models/SEC_features.txt" \
    --dis_features "$BIN_DIR/../models/DIS_features.txt" \
    --best_thresholds "$BIN_DIR/../models/Best_thresholds.tsv" \
    --delta "$DELTA" \
    --min_active "$MIN_ACTIVE" \
    --out "$OUT_DIR/MAG_pred_SEC_DIS.csv"

echo "[DONE] MiFate full pipeline completed. Output in $OUT_DIR"
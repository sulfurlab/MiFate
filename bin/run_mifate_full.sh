#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 -i GENOMES -q QUALITY_TABLE -o OUTPUT -d KOFAM_PROFILES -k KO_LIST [-t THREADS] [-p PYTHON]"
  echo "  GENOMES: one .fna file, a directory, or a ZIP archive"
  echo "  QUALITY_TABLE: CSV/TSV with genome_id, completeness, contamination"
  echo "  KOFAM_PROFILES: KOfam profile HMM directory or .hal file"
  echo "  KO_LIST: separate KOfam ko_list file with KO-specific thresholds"
  echo "  OUTPUT: main seven-column table in OUTPUT/predictions/MiFate_predictions.csv"
}

GENOMES=""
QUALITY=""
OUTPUT=""
PROFILE=""
KO_LIST=""
THREADS=4
PYTHON="${PYTHON:-python3}"

while getopts ":i:q:o:d:k:t:p:h" opt; do
  case "$opt" in
    i) GENOMES="$OPTARG" ;;
    q) QUALITY="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    d) PROFILE="$OPTARG" ;;
    k) KO_LIST="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    p) PYTHON="$OPTARG" ;;
    h) usage; exit 0 ;;
    :) echo "Missing argument for -$OPTARG" >&2; usage >&2; exit 2 ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$GENOMES" || -z "$QUALITY" || -z "$OUTPUT" || -z "$PROFILE" || -z "$KO_LIST" ]]; then
  usage >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
exec "$PYTHON" "$REPO_ROOT/src/genome_deployment_pipeline.py" \
  --genomes "$GENOMES" \
  --quality "$QUALITY" \
  --kofam-profile "$PROFILE" \
  --kofam-ko-list "$KO_LIST" \
  --output "$OUTPUT" \
  --threads "$THREADS"

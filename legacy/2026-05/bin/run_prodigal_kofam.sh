#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
用法:
  bash run_prodigal_kofam.sh -i <MAG目录> -o <输出目录> -d <KOfam数据库目录> [选项]

必选参数:
  -i    输入MAGs文件夹（包含 .fna 文件）
  -o    输出结果文件夹
  -d    KOfam数据库目录（里面要有 profiles/ 和 ko_list）

可选参数:
  -m    Prodigal模式，默认 single，可选 single/meta
  -p    Prodigal并行数，默认 20
  -c    KOfam切块数，默认 10
  -t    每个KOfam任务线程数，默认 20
  -j    KOfam并发任务数，默认 5
  -n    dry-run模式：1只检查不运行，默认 0
  -h    显示帮助信息

示例:
  bash run_prodigal_kofam.sh \
    -i /data/xxx/MAGs \
    -o /data/xxx/MAGs_annotation \
    -d /data/db/kofam_db \
    -m single -p 20 -c 10 -t 20 -j 5
EOF
}

############################
# 默认参数
############################
MAG_DIR=""
OUT_ROOT=""
DB_DIR=""

PRODIGAL_MODE="single"
PRODIGAL_JOBS=20
CHUNKS=10
CPU_PER_JOB=20
PARALLEL_JOBS=5
DRYRUN=0

############################
# 解析参数
############################
while getopts ":i:o:d:m:p:c:t:j:n:h" opt; do
  case $opt in
    i) MAG_DIR="$OPTARG" ;;
    o) OUT_ROOT="$OPTARG" ;;
    d) DB_DIR="$OPTARG" ;;
    m) PRODIGAL_MODE="$OPTARG" ;;
    p) PRODIGAL_JOBS="$OPTARG" ;;
    c) CHUNKS="$OPTARG" ;;
    t) CPU_PER_JOB="$OPTARG" ;;
    j) PARALLEL_JOBS="$OPTARG" ;;
    n) DRYRUN="$OPTARG" ;;
    h)
      usage
      exit 0
      ;;
    \?)
      echo "[ERR] 无效参数: -$OPTARG"
      usage
      exit 1
      ;;
    :)
      echo "[ERR] 参数 -$OPTARG 缺少值"
      usage
      exit 1
      ;;
  esac
done

############################
# 参数检查
############################
if [[ -z "$MAG_DIR" || -z "$OUT_ROOT" || -z "$DB_DIR" ]]; then
  echo "[ERR] -i, -o, -d 为必填参数"
  usage
  exit 1
fi

MAG_DIR=$(realpath "$MAG_DIR")
OUT_ROOT=$(realpath -m "$OUT_ROOT")
DB_DIR=$(realpath "$DB_DIR")

[[ -d "$MAG_DIR" ]] || { echo "[ERR] 输入MAG目录不存在: $MAG_DIR"; exit 2; }
[[ -d "$DB_DIR/profiles" && -s "$DB_DIR/ko_list" ]] || { echo "[ERR] KOfam数据库不完整: $DB_DIR"; exit 3; }

############################
# 输出目录结构
############################
PROJECT_DIR="$OUT_ROOT"
PROTEIN_DIR="${PROJECT_DIR}/01_proteins"
FUNCTION_DIR="${PROJECT_DIR}/02_kofam"
CHUNK_DIR="${FUNCTION_DIR}/_chunks"
TMP_BASE="${PROJECT_DIR}/_tmp_kofam"
LIST_FILE="${PROJECT_DIR}/list.txt"
ALLFAA="${PROJECT_DIR}/all.faa"
FINAL_MAPPER="${FUNCTION_DIR}/kofam.mapper"

mkdir -p "$PROJECT_DIR" "$PROTEIN_DIR" "$FUNCTION_DIR" "$CHUNK_DIR" "$TMP_BASE"

############################
# 检查依赖
############################
echo "[CHECK] 检查依赖..."
command -v prodigal >/dev/null 2>&1 || { echo "[ERR] 找不到 prodigal"; exit 4; }
command -v exec_annotation >/dev/null 2>&1 || { echo "[ERR] 找不到 exec_annotation"; exit 5; }

############################
# 列出 .fna 文件
############################
find "$MAG_DIR" -maxdepth 1 -type f -name "*.fna" | sort > "$LIST_FILE"

if [[ ! -s "$LIST_FILE" ]]; then
  echo "[ERR] 在 $MAG_DIR 下未找到 .fna 文件"
  exit 6
fi

N_MAGS=$(wc -l < "$LIST_FILE" | xargs)
echo "[INFO] 输入MAG目录: $MAG_DIR"
echo "[INFO] 输出结果目录: $OUT_ROOT"
echo "[INFO] 检测到 $N_MAGS 个 .fna 文件"

if [[ "$DRYRUN" -eq 1 ]]; then
  echo "[DRYRUN] 仅检查参数和输入，不实际运行"
  echo "[DRYRUN] Prodigal模式=$PRODIGAL_MODE"
  echo "[DRYRUN] Prodigal并行=$PRODIGAL_JOBS"
  echo "[DRYRUN] KOfam切块=$CHUNKS"
  echo "[DRYRUN] 每任务线程=$CPU_PER_JOB"
  echo "[DRYRUN] 并发任务=$PARALLEL_JOBS"
  exit 0
fi

############################
# Step 1. Prodigal预测
############################
echo "[RUN] 开始 Prodigal 蛋白预测..."

xargs -a "$LIST_FILE" -I{} -P "$PRODIGAL_JOBS" bash -c '
  set -euo pipefail
  fa="{}"
  mag=$(basename "$fa" .fna)

  out="'"$PROTEIN_DIR"'/${mag}.tag.faa"
  tmp="'"$PROTEIN_DIR"'/${mag}.faa"

  [[ -s "$out" ]] && { echo "[SKIP] $mag"; exit 0; }

  echo "[PRED] $mag"
  prodigal -i "$fa" -a "$tmp" -p "'"$PRODIGAL_MODE"'" -q

  awk -v M="$mag" '"'"'/^>/{print ">"M"|"substr($0,2); next} {print}'"'"' "$tmp" > "$out"
  rm -f "$tmp"
'

echo "[OK] Prodigal完成"

############################
# Step 2. 合并蛋白
############################
echo "[MERGE] 合并所有蛋白到 all.faa ..."
find "$PROTEIN_DIR" -type f -name "*.tag.faa" | sort | xargs cat > "$ALLFAA"
[[ -s "$ALLFAA" ]] || { echo "[ERR] all.faa为空"; exit 7; }

SEQ_COUNT=$(grep -c '^>' "$ALLFAA" || true)
echo "[INFO] all.faa 序列数: $SEQ_COUNT"

############################
# Step 3. 切分 all.faa
############################
if ! ls "$CHUNK_DIR"/chunk_*.faa >/dev/null 2>&1; then
  if command -v seqkit >/dev/null 2>&1; then
    echo "[SPLIT] 使用 seqkit 切分 all.faa ..."
    rm -f "$CHUNK_DIR"/*
    seqkit split -p "$CHUNKS" -O "$CHUNK_DIR" "$ALLFAA" >/dev/null
    find "$CHUNK_DIR" -type f -name '*.fa*' | sort | nl -n rz -w 3 | while read -r n f; do
      mv "$f" "$CHUNK_DIR/chunk_${n}.faa"
    done
  else
    echo "[SPLIT] 未检测到 seqkit，使用 awk 切分..."
    rm -f "$CHUNK_DIR"/chunk_*.faa
    N=$(grep -c '^>' "$ALLFAA" || true)
    [[ "$N" -gt 0 ]] || { echo "[ERR] all.faa中没有FASTA记录"; exit 8; }

    per=$(( (N + CHUNKS - 1) / CHUNKS ))
    k=1
    c=0
    out="$CHUNK_DIR/chunk_$(printf %03d "$k").faa"
    : > "$out"

    while IFS= read -r line; do
      if [[ "$line" == \>* ]]; then
        c=$((c+1))
        if (( c > per )); then
          k=$((k+1))
          c=1
          out="$CHUNK_DIR/chunk_$(printf %03d "$k").faa"
          : > "$out"
        fi
      fi
      printf "%s\n" "$line" >> "$out"
    done < "$ALLFAA"
  fi
else
  echo "[SKIP] 检测到已有chunk文件，跳过切分"
fi

NUM_CHUNKS=$(ls "$CHUNK_DIR"/chunk_*.faa | wc -l | xargs)
[[ "$NUM_CHUNKS" -ge 1 ]] || { echo "[ERR] 没有成功生成chunk"; exit 9; }
echo "[INFO] 实际分块数: $NUM_CHUNKS"

############################
# Step 4. KOfam并行注释
############################
echo "[RUN] 开始并行 KOfam 注释..."

export DB_DIR FUNCTION_DIR TMP_BASE CPU_PER_JOB

find "$CHUNK_DIR" -maxdepth 1 -type f -name 'chunk_*.faa' | sort | \
xargs -I{} -P "$PARALLEL_JOBS" bash -c '
  set -euo pipefail
  F="{}"
  B=$(basename "$F" .faa)
  O="$FUNCTION_DIR/${B}.mapper"
  T="$TMP_BASE/${B}"

  mkdir -p "$T"
  [[ -s "$O" ]] && { echo "[SKIP] $B"; exit 0; }

  echo "[JOB ] $(date +%F" "%T) 开始 $B"

  if command -v ionice >/dev/null 2>&1; then
    IO="ionice -c2 -n0"
  else
    IO=""
  fi

  $IO exec_annotation \
    -p "$DB_DIR/profiles" \
    -k "$DB_DIR/ko_list" \
    -f mapper \
    -o "$O" \
    --cpu "$CPU_PER_JOB" \
    --tmp-dir "$T" \
    "$F"

  echo "[DONE] $(date +%F" "%T) 完成 $B"
'

############################
# Step 5. 合并注释结果
############################
echo "[MERGE] 合并 mapper ..."
cat "$FUNCTION_DIR"/chunk_*.mapper > "$FINAL_MAPPER"
[[ -s "$FINAL_MAPPER" ]] || { echo "[ERR] 最终 kofam.mapper 为空"; exit 10; }

############################
# Step 6. 简单验收
############################
echo "[CHECK] 统计KO注释行数..."
awk -F"\t" 'NF>1 && $2!=""{c++} END{print "lines_with_KO:", c+0}' "$FINAL_MAPPER" || true

echo "[DONE] 全流程完成"
echo "[OUT ] 蛋白文件目录: $PROTEIN_DIR"
echo "[OUT ] 合并蛋白文件: $ALLFAA"
echo "[OUT ] KOfam结果目录: $FUNCTION_DIR"
echo "[OUT ] 最终注释文件: $FINAL_MAPPER"
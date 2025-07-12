#!/bin/bash
set -euo pipefail


# 使い方: 
# bash scripts/02_merge.sh \
#   --ref data/hg38.fa \
#   --vcf_list data/call_vcf.list


# 引数
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref) REF="$2"; shift ;;
    --vcf_list) VCF_LIST="$2"; shift ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
  shift
done


# いろいろ準備
REF_DIR=$(readlink -f "$(dirname "$REF")")
DATA_DIR=$(readlink -f "$(dirname "$VCF_LIST")")
REF_BASE=$(basename "$REF")
VCF_LIST_BASE=$(basename "$VCF_LIST")
OUT_DIR="results/02"

mkdir -p "$OUT_DIR"


# 実行 Get the union of sites across all samples
echo "[`date`] Starting smoove merge"

docker run --rm \
  -v "$REF_DIR":/ref \
  -v "$DATA_DIR":/data \
  -v "$(pwd)/results":/results \
  brentp/smoove smoove merge \
  --name merged \
  -f /ref/"$REF_BASE" \
  --outdir /results/02/ \
  /data/"$VCF_LIST_BASE"
  
echo "[`date`] Finished smoove merge"

#!/bin/bash
set -euo pipefail


# 使い方: 
# bash scripts/04_smoove_paste.sh \
#   --vcf_dir results/03/ \
#   --vcf_list results/03/genotype_vcf.list


# 引数
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf_dir) VCF_DIR="$2"; shift ;;
    --vcf_list) VCF_LIST="$2"; shift ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
  shift
done


# いろいろ準備
VCF_DIR=$(readlink -f "$VCF_DIR")
VCF_LIST_DIR=$(readlink -f "$(dirname "$VCF_LIST")")
VCF_LIST_BASE=$(basename "$VCF_LIST")
OUT_DIR="results/04"

mkdir -p "$OUT_DIR"


# 実行
echo "[`date`] Starting smoove paste"

docker run --rm \
  -v "$VCF_DIR":/vcf \
  -v "$VCF_LIST_DIR":/vcf_list \
  -v "$(pwd)/results":/results \
  brentp/smoove smoove paste \
  --name merged \
  --outdir /results/04/ \
  /vcf_list/"$VCF_LIST_BASE"
  
echo "[`date`] Finished smoove paste"

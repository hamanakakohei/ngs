#!/bin/bash
set -euo pipefail


# 引数
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref) REF="$2"; shift ;;
    --cram) CRAM="$2"; shift ;;
    --site_vcf) SITE_VCF="$2"; shift ;;
    --sample) SAMPLE="$2"; shift ;;
    --threads) THREADS="$2"; shift ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
  shift
done


# いろいろ準備
REF_DIR=$(readlink -f "$(dirname "$REF")")
CRAM_DIR=$(readlink -f "$(dirname "$CRAM")")
SITE_VCF_DIR=$(readlink -f "$(dirname "$SITE_VCF")")
REF_BASE=$(basename "$REF")
CRAM_BASE=$(basename "$CRAM")
SITE_VCF_BASE=$(basename "$SITE_VCF")

OUT_DIR="results/03/${SAMPLE}"

mkdir -p "$OUT_DIR"


# 実行 genotype each sample at all merged sites
echo "[`date`] Starting smoove genotype for $SAMPLE"

docker run --rm \
  -v "$REF_DIR":/ref \
  -v "$CRAM_DIR":/cram \
  -v "$SITE_VCF_DIR":/site_vcf \
  -v "$(pwd)/results":/results \
  brentp/smoove smoove genotype \
    -d -x \
    -p "$THREADS" \
    --name "$SAMPLE" \
    --vcf /site_vcf/"$SITE_VCF_BASE" \
    --fasta /ref/"$REF_BASE" \
    --outdir /results/03/"$SAMPLE"/ \
    /cram/"$CRAM_BASE"

echo "[`date`] Finished smoove genotype for $SAMPLE"

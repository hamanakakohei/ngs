#!/bin/bash
set -euo pipefail


# 引数
while [[ $# -gt 0 ]]; do
  case "$1" in
    --kmer)          KMER="$2"; shift ;;
    --ref)           REF="$2"; shift ;;
    --bam_info)      BAM_INFO="$2"; shift ;;
    --filter_bed)    FILTER_BED="$2"; shift ;;
    --genome_folder) GENOME_FOLDER="$2"; shift ;;
    --pop_vcf)       POP_VCF="$2"; shift ;;
    --out)           OUT_DIR="$2"; shift ;;
    *) echo "Unknown argument: $1" >&2 ; exit 1 ;;
  esac
  shift
done


mkdir -p "$OUT_DIR"


# Canvasコマンドの組み立て
CMD=(dotnet Canvas.dll SmallPedigree-WGS)


# BAM_INFO を分解して --bam に追加
IFS=';' read -ra ENTRIES <<< "$BAM_INFO"
for entry in "${ENTRIES[@]}"; do
  IFS=',' read -r BAM_PATH ROLE SAMPLE_NAME <<< "$entry"
  CMD+=(--bam="$BAM_PATH" "$ROLE" "$SAMPLE_NAME")
done


# 他の引数を追加
CMD+=(
  -r "$KMER"
  -g "$GENOME_FOLDER"
  --population-b-allele-vcf "$POP_VCF"
  -f "$FILTER_BED"
  -o "$OUT_DIR"
)
  # --ploidy-vcf="/basespace/Projects/canvas/AppResults/snvvcf/Files/MultiSamplePloidy.vcf"


# コマンドを見せて実行
echo "Running:"
printf ' %q\n' "${CMD[@]}"
echo
"${CMD[@]}"

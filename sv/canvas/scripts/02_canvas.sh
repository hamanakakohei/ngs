#!/bin/bash
set -euo pipefail


# 引数
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf_list) VCF_LIST="$2"; shift ;;
    --out) OUT="$2"; shift ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
  shift
done


# 実行
(
  zcat "$(head -n 1 "$VCF_LIST")" | grep '^##'
  zcat "$(head -n 1 "$VCF_LIST")" | grep '^#CHROM'

  while read file; do
    zcat "$file" | grep -v '^#'
  done < "$VCF_LIST"
) | bgzip > "$OUT"

tabix -p vcf "$OUT"

#!/bin/bash
set -euo pipefail


# 引数
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref) REF="$2"; shift ;;
    --vcf) VCF="$2"; shift ;;
    --out) OUT="$2"; shift ;;
    --threads) THREADS="$2"; shift ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
  shift
done


# 実行
gatk SelectVariants \
  -V $VCF \
  -O $OUT \
  -R $REF \
  --sites_only \
  --select-type-to-include SNP \
  --restrict-alleles-to BIALLELIC
  #--selectExpressions "AF > 0.1" 

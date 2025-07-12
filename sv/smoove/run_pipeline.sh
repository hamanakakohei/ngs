#!/bin/bash

set -euo pipefail


SAMPLE_CRAM_LIST=data/sample_cram.txt
REF=data/hg38.fa
EXCLUDE_BED=data/exclude.cnvnator_100bp.GRCh38.20170403.bed


N_SAMPLE=$(wc -l < "$SAMPLE_CRAM_LIST")

# 01
qsub -t 1-"$N_SAMPLE":1 -tc 10 qsubs/01_smoove_call.qsub \
  "$SAMPLE_CRAM_LIST" \
  "$REF" \
  "$EXCLUDE_BED"

find results/01/ -name "*.vcf.gz" > results/01/call_vcf.list


# 02
bash scripts/02_smoove_merge.sh \
  --ref "$REF" \
  --vcf_list results/01/call_vcf.list \
  > logs/02/log 2>&1


# 03
qsub -t 1-"$N_SAMPLE":1 -tc 10 qsubs/03_smoove_genotype.qsub \
  "$SAMPLE_CRAM_LIST" \
  "$REF" \
  results/02/merged.sites.vcf.gz

find results/03/ -name "*.vcf.gz" > results/03/genotype_vcf.list


# 04
bash scripts/04_smoove_paste.sh \
  --vcf_dir results/03/ \
  --vcf_list results/03/genotype_vcf.list \
  > logs/04/log 2>&1



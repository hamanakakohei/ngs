#!/bin/bash

set -euo pipefail


SAMPLE_CRAM_LIST=data/sample_cram.txt
REF=data/Homo_sapiens_assembly38.fasta

MANIFEST=data/manifest.txt
PREFIX=all_samples
ANNOVAR=annovar20250302/annotate_variation.pl
HUMNA_DB=annovar/201612/hg38

PROBANDS=data/probands.txt
CONTROLS=data/controls.txt
PROBAND_RELATIVE=data/proband-relative_pairs.txt

N_SAMPLE=$(wc -l < "$SAMPLE_CRAM_LIST")


# 01 全サンプルでEHdn profileする
qsub -t 1-"$N_SAMPLE":1 -tc 10 qsubs/01.qsub \
  "$SAMPLE_CRAM_LIST" \
  "$REF"


# 02 コールが多すぎるサンプルを見つける
# この結果を見てマニフェストファイルを編集して、03から再開する
mkdir -p logs/02

conda run -n ehdn_env \
  scripts/02_QC.py \
  > logs/02/log 2>&1


# 03 結果をマージして、outlier解析をして、Annovarをかける
mkdir -p logs/03

scripts/03_ehdn_merge_and_outlier.sh \
  --manifest "$MANIFEST" \
  --prefix   "$PREFIX" \
  --ref      "$REF" \
  --annovar  "$ANNOVAR" \
  --human_db "$HUMAN_DB" \
  > logs/03/log 2>&1


# 04：アレル頻度を付けて、outlier motifとoutlier locus結果をマージする
qsub -t 1-"$N_SAMPLE":1 -tc 10 qsubs/04.qsub \
  "$PROBANDS" \
  "$CONTROLS" \
  "$PROBAND_RELATIVE" \
  "$PREFIX"

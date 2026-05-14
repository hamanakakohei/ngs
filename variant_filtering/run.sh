#!/usr/bin/env bash
#
# 0. xhmmファイルの下準備する
# 1. sample_annovar_xhmm.txtファイルを空白区切りでexome_summaryパスと上で作ったxhmmファイルパスを書き込んで使う
# - トリオなら-mフラグは外す
set -euo pipefail

# 0
SAMPLES=(
Sample_aaa
Sample_bbb
)
XHMM_RES=/path/to/data.segdup.strvar.haplo.deciph.omim.xcnv.gene

for SAMPLE in ${SAMPLES[@]}; do
        ./preprocess_XHMM.py \
                --xhmm_result $XHMM_RES \
                --sample $SAMPLE \
                --depth_threshold_for_homo_del 8 \
                --out results/01/cleaned_XHMM_${SAMPLE}.tsv.gz
done


# 1
GENE_ANNOS=(
results/cleaned_G2P_2025-04-28.tsv.gz
results/cleaned_GenCC_2025_05_20.tsv.gz
results/cleaned_PanelApp_ataxia_2025_05_20.tsv.gz
)

while read SAMPLE ANNOVAR XHMM; do
echo $SAMPLE
scripts/pipeline.sh \
        -s $SAMPLE \
        -a $ANNOVAR \
        -o results/02/${SAMPLE}.filtered \
        -m \
        -c $XHMM \
        -g "${GENE_ANNOS[*]}"
done < sample_annovar_xhmm.txt

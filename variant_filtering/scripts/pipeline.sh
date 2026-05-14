#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 -s TARGET_SAMPLE -a ANNOVAR_RES -o OUT_PREFIX [-c OTHER_CALLER] [-m] [-g gene1 gene2 gene3]"
    echo "  -s  target sample (required)"
    echo "  -a  annovar result (required)"
    echo "  -o  output prefix (required)"
    echo "  -c  other caller results (optional)"
    echo "  -m  enable gene_mode_of_inheritance_filter (flag)"
    echo "  -g  gene annotation files, space-separated, default 3 files used if not specified"
    exit 1
}

# デフォルト
GENE_ANNOTATIONS=()
USE_MOI_FILTER=false
OTHER_CALLER=""

while getopts "s:a:o:c:mg:" opt; do
    case $opt in
        s) TARGET_SAMPLE=$OPTARG ;;
        a) ANNOVAR_RES=$OPTARG ;;
        o) OUT_PREFIX=$OPTARG ;;
        c) OTHER_CALLER=$OPTARG ;;
        m) USE_MOI_FILTER=true ;;
        g) IFS=' ' read -ra GENE_ANNOTATIONS <<< "$OPTARG" ;;
        *) usage ;;
    esac
done

# 必須引数チェック
[[ -z $TARGET_SAMPLE || -z $ANNOVAR_RES || -z $OUT_PREFIX ]] && usage

# オプション引数を組み立て
EXTRA_OPTS=""
$USE_MOI_FILTER && EXTRA_OPTS="$EXTRA_OPTS --gene_mode_of_inheritance_filter"
[[ -n $OTHER_CALLER ]] && EXTRA_OPTS="$EXTRA_OPTS --other_caller_results $OTHER_CALLER"

run_filter() {
    local inheritance=$1
    local out=$2
    local af_hgvd=$3
    local af_tommo=$4
    local af_exac_all=$5
    local af_exac_eas=$6

    scripts/filter_annovar_result.R \
        --annovar_result $ANNOVAR_RES \
        --out $out \
        --af_threshold_hgvd $af_hgvd \
        --af_threshold_tommo $af_tommo \
        --af_threshold_exac_all $af_exac_all \
        --af_threshold_exac_eas $af_exac_eas \
        --inheritance $inheritance \
        --sample_filter $TARGET_SAMPLE \
        --gene_annotations "${GENE_ANNOTATIONS[@]}" \
        $EXTRA_OPTS
}

# AD用
run_filter AD ${OUT_PREFIX}.AD.txt 0.01 0.0005 0.00005 0.0005

# AR用
run_filter AR ${OUT_PREFIX}.AR.txt 0.02 0.02 0.02 0.02

# XL用
run_filter XL ${OUT_PREFIX}.XL.txt 0.01 0.01 0.00005 0.0005

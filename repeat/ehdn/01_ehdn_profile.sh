#!/bin/bash
set -euo pipefail


show_help() {
    echo "Usage: $0 --sample NAME --bam FILE --outdir DIR [--ref FILE]"
    echo ""
    echo "  --sample   サンプル名（出力ファイルのプレフィックス）"
    echo "  --bam      BAM または CRAM ファイルのパス"
    echo "  --outdir   出力ファイルを保存するディレクトリ"
    echo "  --ref      参照ファイル（省略可、デフォルト使用）"
    exit 1
}


# デフォルト
REF="resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
EHDN="ExpansionHunter-0.9.0/bin/ExpansionHunterDenovo"


# 引数の初期化
SAMPLE=""
BAM=""
OUTDIR=""


# 引数の解析
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)
            SAMPLE="$2"
            shift 2
            ;;
        --bam)
            BAM="$2"
            shift 2
            ;;
        --ref)
            REF="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        -*)
            echo "Unknown option: $1" >&2
            show_help ;;
        *)
            echo "Unexpected argument: $1" >&2
            show_help ;;
    esac
done


# 必須引数チェック
if [[ -z "$SAMPLE" || -z "$BAM" || -z "$OUTDIR" ]]; then
    echo "Usage: $0 --sample NAME --bam FILE --outdir DIR [--ref FILE]" >&2
    exit 1
fi


# 出力ディレクトリ作成
mkdir -p "$OUTDIR"


# EHdn profile を1サンプルで実行する
$EHDN profile \
  --log-reads \
  --reads "$BAM" \
  --reference "$REF" \
  --output-prefix "${OUTDIR}/${SAMPLE}.ehdn.profile" \
  --min-anchor-mapq 50 \
  --max-irr-mapq 40 \
  > "${OUTDIR}/${SAMPLE}.ehdn.profile.log" 2>&1

#!/bin/bash
set -euo pipefail


# 固定している変数
EHDN="ExpansionHunter-0.9.0/bin/ExpansionHunterDenovo"


# 引数の解析
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2 ;;
    --bam)    BAM="$2";    shift 2 ;;
    --ref)    REF="$2";    shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    -*)
      echo "Unknown option: $1" >&2
      exit 1 ;;
    *)
      echo "Unexpected argument: $1" >&2
      exit 1 ;;
  esac
done


# 出力ディレクトリ作成
: "${OUTDIR:?--outdir is required}"
mkdir -p "$OUTDIR"


# 実行
$EHDN profile \
  --log-reads \
  --reads "$BAM" \
  --reference "$REF" \
  --output-prefix "${OUTDIR}/${SAMPLE}.ehdn.profile" \
  --min-anchor-mapq 50 \
  --max-irr-mapq 40

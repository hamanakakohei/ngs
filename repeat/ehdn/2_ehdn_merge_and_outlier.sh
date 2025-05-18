#!/bin/bash
set -euo pipefail


show_help() {
    cat <<EOF
Usage: $0 --manifest FILE --prefix NAME --ref FILE [--annovar FILE --human_db DIR]

Options:
  --manifest     マニフェストファイルのパス
  --prefix       出力ファイルのプレフィックス
  --ref          参照ゲノムFASTA
  --annovar      Annovar の annotate_variation.pl スクリプト（省略可）
  --human_db     Annovar のヒト用DBディレクトリ（省略可）

マニフェストファイルの書き方：
https://github.com/Illumina/ExpansionHunterDenovo/blob/master/documentation/06_Merging_profiles.md

関連スクリプトは以下からダウンロードする:
https://github.com/Illumina/ExpansionHunterDenovo/tree/master/scripts/outlier.py
https://github.com/Illumina/ExpansionHunterDenovo/tree/master/scripts/annotate_ehdn.py

Annovarについては以下を参照する:
https://github.com/Illumina/ExpansionHunterDenovo/blob/master/documentation/08_Annotation.md

EOF
    exit 0
}


# デフォルト
EHDN=ExpansionHunter-0.9.0/bin/ExpansionHunterDenovo 
OUTLIER=ExpansionHunterDenovo/scripts/outlier.py
ANNOTATE=ExpansionHunterDenovo/scripts/annotate_ehdn.sh

HUMAN_DB=annovar/201612/hg38
ANNOVAR=annovar2016oct/annotate_variation.pl


# 引数の初期化
MANIFEST=""
PREFIX=""


# 引数の解析
while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)
            MANIFEST="$2"
            shift 2
            ;;
        --prefix)
            PREFIX="$2"
            shift 2
            ;;
        --ref)
            REF="$2"
            shift 2
            ;;
        --annovar)
            ANNOVAR="$2"
            shift 2
            ;;
        --human_db)
            HUMAN_DB="$2"
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


# 必須引数のチェック
if [[ -z "$MANIFEST" || -z "$PREFIX" || -z "$REF" ]]; then
    echo "[ERROR] --manifest、--prefix、--ref は必須です。" >&2
    show_help
fi


# １：各サンプルのEHdn profileの結果をマージする
${EHDN} merge \
  --manifest 		${MANIFEST} \
  --reference 		${REF} \
  --output-prefix 	${PREFIX} \
  > ${PREFIX}.merge.log 2>&1


# ２：locusの観点で外れ値解析をする
${OUTLIER} locus \
  --manifest 		${MANIFEST} \
  --multisample-profile ${PREFIX}.multisample_profile.json \
  --output 		${PREFIX}.outlier_locus.tsv \
  > ${PREFIX}.outlier_locus.log 2>&1


# ３：motifの観点で外れ値解析をする
${OUTLIER} motif \
  --manifest 		${MANIFEST} \
  --multisample-profile ${PREFIX}.multisample_profile.json \
  --output 		${PREFIX}.outlier_motif.tsv \
  > ${PREFIX}.outlier_motif.log 2>&1


# ４：locusの外れ値解析結果に遺伝子アノテーションを付ける
${ANNOTATE} \
  --ehdn-results 		${PREFIX}.outlier_locus.tsv \
  --ehdn-annotated-results 	${PREFIX}.outlier_locus.annotated.tsv \
  --annovar-annotate-variation 	${ANNOVAR} \
  --annovar-humandb 		${HUMAN_DB} \
  --annovar-buildver 		hg38 \
  > ${PREFIX}.outlier_locus.annotated.log 2>&1


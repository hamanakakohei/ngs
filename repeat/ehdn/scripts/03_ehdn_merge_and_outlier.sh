#!/bin/bash
set -euo pipefail


show_help() {
  echo "マニフェストファイルの書き方："
  echo "https://github.com/Illumina/ExpansionHunterDenovo/blob/master/documentation/06_Merging_profiles.md"
  echo ""
  echo "関連スクリプトは以下からダウンロードする:"
  echo "https://github.com/Illumina/ExpansionHunterDenovo/tree/master/scripts/outlier.py"
  echo "https://github.com/Illumina/ExpansionHunterDenovo/tree/master/scripts/annotate_ehdn.py"
  echo ""
  echo "Annovarについては以下を参照する:"
  echo "https://github.com/Illumina/ExpansionHunterDenovo/blob/master/documentation/08_Annotation.md"
}


# 固定している変数
EHDN=ExpansionHunter-0.9.0/bin/ExpansionHunterDenovo 
OUTLIER=ExpansionHunterDenovo/scripts/outlier.py
ANNOTATE=ExpansionHunterDenovo/scripts/annotate_ehdn.sh


# 引数の解析
while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --prefix)   PREFIX="$2";   shift 2 ;;
    --ref)      REF="$2";      shift 2 ;;
    --annovar)  ANNOVAR="$2";  shift 2 ;;
    --human_db) HUMAN_DB="$2"; shift 2 ;;
    -*)
      echo "Unknown option: $1" >&2; exit 1 ;;
    *)
      echo "Unexpected argument: $1" >&2; exit 1 ;;
  esac
done


# 各サンプルのEHdn profileの結果をマージする
echo "Running merge"
${EHDN} merge \
  --manifest 		${MANIFEST} \
  --reference 		${REF} \
  --output-prefix 	${PREFIX}


# locusについて外れ値解析をする
echo "Running locus"
${OUTLIER} locus \
  --manifest 		${MANIFEST} \
  --multisample-profile ${PREFIX}.multisample_profile.json \
  --output 		${PREFIX}.outlier_locus.tsv


# motifについて外れ値解析をする
echo "Running motif"
${OUTLIER} motif \
  --manifest 		${MANIFEST} \
  --multisample-profile ${PREFIX}.multisample_profile.json \
  --output 		${PREFIX}.outlier_motif.tsv


# locusの結果に遺伝子アノテーションを付ける
echo "Running annotate"
${ANNOTATE} \
  --ehdn-results 		${PREFIX}.outlier_locus.tsv \
  --ehdn-annotated-results 	${PREFIX}.outlier_locus.annotated.tsv \
  --annovar-annotate-variation 	${ANNOVAR} \
  --annovar-humandb 		${HUMAN_DB} \
  --annovar-buildver 		hg38

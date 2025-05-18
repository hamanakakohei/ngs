#!/bin/bash


SAMPLE_CRAM_PAIRS=data/sample-cram_pairs.list
REF=
MANIFEST=data/manifest.list
ANNOVAR=annovar2016oct/annotate_variation.pl
HUMNA_DB=annovar/201612/hg38
FAMILY_SAMPLE_PAIRS=data/family-sample_pairs.list
SAMPLES=data/samples.list
CONTROLS=data/controls.list


# ステップ１：全サンプルで EHdn profile をする
while read LINE; do
	SAMPLE=`echo $LINE | awk '{print $1}'`
	CRAM=`echo $LINE | awk '{print $2}'`

	./01_ehdn_profile.sh \
		--sample "$SAMPLE" \
		--bam    "$CRAM" \
		--outdir results_step1/"$SAMPLE" \
		--ref    "$REF"
done < "$SAMPLE_CRAM_PAIRS"


# ステップ２：結果をマージして、outlier解析をして、Annovarをかける
./02_ehdn_merge_and_outlier.sh \
	--manifest 	"$MANIFEST" \
	--prefix 	all_samples \
	--ref 		"$REF" \
	--annovar 	"$ANNOVAR" \
	--human_db 	"$HUMAN_DB"


# ステップ３：アレル頻度を付けて、outlier motifとoutlier locus結果をマージする
while read SAMPLE; do
	awk -v SAMPLE="${SAMPLE}" '$1==SAMPLE && $2!=SAMPLE {print $2}' "$FAMILY_SAMPLE_PAIRS" > result_step3/"$SAMPLE".other_family_members.txt

	./03_add_allele_count.py \
		--proband 		"$SAMPLE" \
		--other_family_members 	result_step3/"$SAMPLE".other_family_members.txt \
		--controls      	"$CONTROLS" \
		--outlier_locus 	result_step2/all_samples.outlier_locus.annotated.tsv \
		--outlier_motif 	result_step2/all_samples.outlier_motif.tsv \
		--out 			result_step3/"$SAMPLE"
done < "$SAMPLES"



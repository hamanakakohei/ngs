#!/bin/bash


# 詳しい解説：
# - https://github.com/Illumina/canvas
# - https://github.com/Illumina/canvas/wiki
# - https://github.com/Illumina/canvas/blob/master/SoftwareDesignDescription.pdf

# 必要なファイルを以下のようにダウンロードする
# wget http://canvas-cnv-public.s3.amazonaws.com/hg19/WholeGenomeFasta/GenomeSize.xml

# 提供されているファイルのリスト：
# GRCh37/WholeGenomeFasta/GenomeSize.xml
# GRCh37/WholeGenomeFasta/genome.fa
# GRCh37/WholeGenomeFasta/genome.fa.fai
# GRCh37/dbsnp.vcf
# GRCh37/filter13.bed
# GRCh37/kmer.fa
# GRCh37/kmer.fa.fai
# GRCh38/WholeGenomeFasta/GenomeSize.xml
# GRCh38/WholeGenomeFasta/genome.fa
# GRCh38/WholeGenomeFasta/genome.fa.fai
# GRCh38/dbsnp.vcf
# GRCh38/filter13.bed
# GRCh38/kmer.fa
# TruthSets/HCC1187.cnaqc.excluded_regions.bed
# TruthSets/HCC1187Truth.vcf
# TruthSets/HCC2218.cnaqc.excluded_regions.bed
# TruthSets/HCC2218Truth.vcf
# TruthSets/NA12878.deletions.bed
# TruthSets/generic.cnaqc.excluded_regions.bed
# hg19/WholeGenomeFasta/GenomeSize.xml
# hg19/WholeGenomeFasta/genome.fa
# hg19/WholeGenomeFasta/genome.fa.fai
# hg19/dbsnp.vcf
# hg19/filter13.bed
# hg19/kmer.fa
# hg19/kmer.fa.fai


set -euo pipefail
export PATH="/usr/local/genome/samtools-1.17/bin:$PATH"
export PATH=$PATH:Canvas-1.40.0.1613+master_x64/


VCF_IN_OUT=data/vcf_input_output.list
SITE_VCF_LIST=data/site_vcf_each_chr.list 
REF=data/hg38.fa
KMER=data/GRCh38/kmer.fa
FILTER_BED=data/GRCh38/filter13.bed
GENOME_FOLDER=data/WholeGenomeFasta
FAMILY_BAM_INFO_LIST=data/family_bam_info.list

N_FAMILY=$(wc -l < "$FAMILY_BAM_INFO_LIST")


# 01 各染色体のサイトVCFを作る
# canvasコマンドは自動でPASSサイトのみ使うので、ここでフィルターしなくて良い
qsub -t 1-24:1 -tc 5 qsubs/01_canvas.qsub \
  "$VCF_IN_OUT" \
  "$REF"


# 02 それらを結合する
bash scripts/02_canvas.sh \
  --vcf_list "$SITE_VCF_LIST" \
  --out results/02/biallelic_snp_site.vcf.gz \
  > logs/02/log 2>&1


# 03
# 以下を揃えることに注意
# - bam header SM tag
# - multisample SNV VCF sample name
# - ploidy VCF sample name
# - pedigree tag sample name
qsub -t 1-"$N_FAMILY":1 -tc 5 qsubs/03_canvas.qsub \
  $FAMILY_BAM_INFO_LIST \
  results/02/biallelic_snp_site.vcf.gz \
  $REF \
  $GENOME_FOLDER \
  $KMER \
  $FILTER_BED

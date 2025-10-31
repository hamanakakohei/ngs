#!/bin/bash


# 詳しい解説：
# - https://github.com/Illumina/canvas
# - https://github.com/Illumina/canvas/wiki
# - https://github.com/Illumina/canvas/blob/master/SoftwareDesignDescription.pdf

# 必要なファイルを以下のようにダウンロードする
# wget http://canvas-cnv-public.s3.amazonaws.com/hg19/WholeGenomeFasta/GenomeSize.xml

# 提供されているファイルのリスト（まだ他にもある）：
# GRCh38/WholeGenomeFasta/GenomeSize.xml
# GRCh38/WholeGenomeFasta/genome.fa
# GRCh38/WholeGenomeFasta/genome.fa.fai
# GRCh38/dbsnp.vcf
# GRCh38/filter13.bed
# GRCh38/kmer.fa


set -euo pipefail
export PATH="/path/to/samtools-1.17/bin:$PATH"
export PATH=$PATH:Canvas-1.40.0.1613+master_x64/


SNP=data/dbsnp.vcf
REF=data/Homos_sapiens_assembly38.fa
KMER=data/GRCh38/kmer.fa
FILTER_BED=data/GRCh38/filter13.bed
GENOME_FOLDER=data/WholeGenomeFasta
FAMILY_BAM_INFO_LIST=data/family_bam_info.list

N_FAMILY=$(wc -l < "$FAMILY_BAM_INFO_LIST")


## 00用
#VCF_IN_OUT=data/vcf_input_output.list
#SITE_VCF_LIST=data/site_vcf_each_chr.list 
#SNP=results/00-2/biallelic_snp_site.vcf.gz
#
## 00-1 各染色体のサイトVCFを作る
## canvasコマンドは自動でPASSサイトのみ使うので、ここでフィルターしなくて良い
#qsub -t 1-24:1 -tc 5 qsubs/00-1_canvas.qsub \
#  "$VCF_IN_OUT" \
#  "$REF"
#
## 00-2 それらを結合する
#mkdir logs/00-2
#bash scripts/00-2_canvas.sh \
#  --vcf_list "$SITE_VCF_LIST" \
#  --out "$SNP" \
#  > logs/00-2/log 2>&1


# 01
# 以下を揃えることに注意
# - bam header SM tag
# - multisample SNV VCF sample name
# - ploidy VCF sample name
# - pedigree tag sample name
qsub -t 1-"$N_FAMILY":1 -tc 5 qsubs/01_canvas.qsub \
  $FAMILY_BAM_INFO_LIST \
  $SNP \
  $REF \
  $GENOME_FOLDER \
  $KMER \
  $FILTER_BED

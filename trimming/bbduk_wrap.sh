#!/bin/bash
# 目的：bbduk.sh にサンプルのリストを与えて実行する

while read SAMPLE; do

  sample=${SAMPLE} \
    fq1=/path/to/${SAMPLE}_R1.fastq.gz \
    fq2=/path/to/${SAMPLE}_R2.fastq.gz \
    fq1_out=/path/to/${SAMPLE}_R1.trim.fastq.gz \
    fq2_out=/path/to/${SAMPLE}_R2.trim.fastq.gz \
    ref=/path/to/adapter.fa \
    stats=/path/to/${SAMPLE}.stats \
    thread=10 \
    bash ./bbduk.sh

done < samples.txt

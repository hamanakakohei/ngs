#!/bin/bash
# 目的：bbdukでトリミングする
# 使い方：bbduk_wrap.sh を参照

source ~/.bashrc
conda activate misc_20250301


# 外から指定する引数
sample="${sample:-NA18486}"
fq1="${fq1:-/path/to/${sample}_R1.fastq.gz}"
fq2="${fq2:-/path/to/${sample}_R2.fastq.gz}"
fq1_out="${fq1_out:-/path/to/${sample}_R1.fastq.gz}"
fq2_out="${fq2_out:-/path/to/${sample}_R2.fastq.gz}"
ref="${ref:-/path/to/adapter.fa}"
stats="${stats:-/path/to/${sample}.stats}"
thread="${thread:-10}"


#- java verは結果にほぼ影響ないらしいが一応固定したい 
bbduk=~/local/bin/bbmap_39.13/bbduk.sh

${bbduk} \
  in1=${fq1} \
  in2=${fq2} \
  out1=${fq1_out} \
  out2=${fq2_out} \
  ktrim=r \
  k=23 \
  mink=11 \
  hdist=1 \
  minlength=25 \
  maxns=1 \
  tpe \
  tbo \
  ref=${ref} \
  stats=${stats} \
  ordered=t \
  threads=${thread} \
  rcomp=f \
  > ${sample}.log 2>&1

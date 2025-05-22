#!/usr/bin/env python3


# 目的：XHMMのコール結果を、Annovarのコール結果と上下に結合できるように整形する
# どういう列にどういう値を残すかだが、本スクリプトのフィルターで想定外の除かれかたをしないように注意しないといけない
# XHMMのデフォルトの出力結果を使えればシンプルだが、遺伝子との被りを自分で調べないといけなくなるので、
# すでにそれを付けているファイルを利用している（data.segdup.strvar.haplo.deciph.omim.xcnv.gene）
# サンプル名を列名に使ってAnnovar結果と上下に結合するので正確に指定する

import argparse
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser(description='')
parser.add_argument('--xhmm_result', type=str,   required=True, help='data.segdup.strvar.haplo.deciph.omim.xcnv.gene')
parser.add_argument('--sample',      type=str,   required=True, help='SAMPLE列でサンプルを選ぶのに使う')
parser.add_argument('--depth_threshold_for_homo_del',      type=float, default=10,     help='ホモ欠失とみなす平均デプスの閾値')
parser.add_argument('--out',         type=str,   required=True)
args = parser.parse_args()

if args.out is None:
  args.out = f'cleaned_XHMM_{args.sample}.tsv.gz'


# ファイル読む
df = pd.read_table( args.xhmm_result )


# INTERVAL 列を Chr, Start, End に分解
df[['Chr', 'Start', 'End']] = df['INTERVAL'].str.extract(r'([^:]+):(\d+)-(\d+)')


# ジェノタイプを適当につける
df[ args.sample ] = np.where(
  df[ 'MEAN_ORIG_RD' ].astype(float) <= args.depth_threshold_for_homo_del,
  '1/1',
  '0/1'
)


# GeneDetail.refGene にコールのメタデータを詰め込む
def summarize_metadata(row):
  return (
    f"NUM_TARG:{row['NUM_TARG']};"
    f"Q_EXACT:{row['Q_EXACT']};"
    f"Q_SOME:{row['Q_SOME']};"
    f"Q_NON_DIPLOID:{row['Q_NON_DIPLOID']};"
    f"Q_START:{row['Q_START']};"
    f"Q_STOP:{row['Q_STOP']};"
    f"MEAN_RD:{row['MEAN_RD']};"
    f"MEAN_ORIG_RD:{row['MEAN_ORIG_RD']}"
  )

df['Func.refGene'] = 'exonic'
df['Gene.refGene'] = df['gene']
df['GeneDetail.refGene'] = df.apply( summarize_metadata, axis=1 )
df['ExonicFunc.refGene'] = df['CNV']
df['AAChange.refGene'] = df['gene'] + ':' + df['CNV']
df['Caller'] = 'XHMM'


# サンプルを選びつつ、
# 列名を変えつつ（遺伝子名とDEL or DUPの情報）、
# 必要な列だけ選びつつ、保存する
output_cols = [
  'Caller', 
  'Chr', 
  'Start', 
  'End', 
  'Func.refGene', 
  'Gene.refGene', 
  'GeneDetail.refGene',
  'ExonicFunc.refGene', 
  'AAChange.refGene', 
  args.sample
]

df.\
  query('SAMPLE == @args.sample')\
  [ output_cols ].\
  drop_duplicates().\
  to_csv( args.out, sep='\t', index=False )

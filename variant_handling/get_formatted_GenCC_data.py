#!/usr/bin/env python3


# 目的：G2Pが提供するを、遺伝子名をキーにして他の表と結合できるように整形する
# ただしVEPならG2Pプラグインで同様のことを自動でできそうだ
# VEPが対応していないバリアントフォーマットを出力するツール用になら役立つかもしれない

import argparse
import pandas as pd
import sys
import os
script_dir = os.path.dirname( os.path.abspath(__file__) )
general_path = os.path.abspath( os.path.join(script_dir, "../../misc/utils") )
sys.path.append( general_path )
from general import download_from_URLs, make_concatenated_df, extract_filenames_from_urls 
from datetime import date


today_str = date.today().strftime('%Y_%m_%d')  # 例: '2025_05_16'
default_filename = f'GenCC_clean_{today_str}.tsv.gz'


# ステップ０
parser = argparse.ArgumentParser(description='')
parser.add_argument('--urls_file', type=str, required=True, help='URLが1行ずつ書かれたファイルのパス')
parser.add_argument('--out',       type=str, default=default_filename)
args = parser.parse_args()


# ステップ１
# 与えたURLからGenCCデータをダウンロード
with open( args.urls_file, 'r' ) as f:
  urls = [ line.strip() for line in f if line.strip() ]

download_from_URLs( urls )


# ステップ２
files = extract_filenames_from_urls( urls )

make_concatenated_df( files, sep='\t' ).\
  drop_duplicates()\
  [[
    'gene_curie',
    'gene_symbol',
    'disease_title',
    'classification_title',
    'moi_title'
  ]].\
  to_csv( args.out, sep = '\t', index = False )

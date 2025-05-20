#!/usr/bin/env python3


# 目的：GenCCが提供する遺伝形式に関するファイルを、遺伝子名をキーにして他の表と結合できるように整形する
# 今の所、VEPが対応していないので、自前で整形（してバリアントリストと結合）する必要がある

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


# ステップ２：ファイルを縦に結合する
files = extract_filenames_from_urls( urls )

df = make_concatenated_df( files, sep='\t' )


# ステップ３：色々整形する
df.columns = [ f"GenCC_{col}" for col in df.columns ]

interest_cols = [
  'GenCC_gene_curie',
  'GenCC_gene_symbol',
  'GenCC_disease_title',
  'GenCC_classification_title',
  'GenCC_moi_title'
]

group_cols = [ 'GenCC_gene_curie', 'GenCC_gene_symbol' ]
other_cols = [ col for col in interest_cols if col not in group_cols ]

df\
  [ interest_cols ].\
  drop_duplicates().\
  groupby( group_cols ).\
  agg({
    col: lambda x: ";".join( map(str,  x) )
    for col in other_cols
  }).\
  reset_index().\
  rename( {'GenCC_gene_symbol': 'gene'}, axis = 1 ).\
  to_csv( args.out, sep = '\t', index = False )

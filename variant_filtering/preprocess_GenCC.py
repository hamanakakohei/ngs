#!/usr/bin/env python3


# 目的：GenCCが提供する遺伝形式に関するファイルを、遺伝子名をキーにして他の表と結合できるように整形する
# Annovar結果と結合することを考えているので、遺伝子名の列名をGene.refGeneにしている
# 今の所、VEPが対応していないので、自前で整形（してバリアントリストと結合）する必要がある
# GenCC はファイルが一つしか無いが一応複数ファイルに対応する書き方にしている

# GenCCファイルの列名とその値の例：
#  1 "uuid"                                 "GENCC_000101-HGNC_10896-OMIM_182212-HP_0000006-GENCC_100001"
#  2 "gene_curie"                           "HGNC:10896"
#  3 "gene_symbol"                          "SKI"
#  4 "disease_curie"                        "MONDO:0008426"
#  5 "disease_title"                        "Shprintzen-Goldberg_syndrome"
#  6 "disease_original_curie"               "OMIM:182212"
#  7 "disease_original_title"               "Shprintzen-Goldberg_syndrome"
#  8 "classification_curie"                 "GENCC:100001"
#  9 "classification_title"                 "Definitive"
# 10 "moi_curie"                            "HP:0000006"
# 11 "moi_title"                            "Autosomal_dominant"
# 12 "submitter_curie"                      "GENCC:000101"
# 13 "submitter_title"                      "Ambry_Genetics"
# 14 "submitted_as_hgnc_id"                 "HGNC:10896"
# 15 "submitted_as_hgnc_symbol"             "SKI"
# 16 "submitted_as_disease_id"              "OMIM:182212"
# 17 "submitted_as_disease_name"            "Shprintzen-Goldberg_syndrome"
# 18 "submitted_as_moi_id"                  "HP:0000006"
# 19 "submitted_as_moi_name"                "Autosomal_dominant_inheritance"
# 20 "submitted_as_submitter_id"            "GENCC:000101"
# 21 "submitted_as_submitter_name"          "Ambry_Genetics"
# 22 "submitted_as_classification_id"       "GENCC:100001"
# 23 "submitted_as_classification_name"     "Definitive"
# 24 "submitted_as_date"                    "2018-03-30_13:31:56"
# 25 "submitted_as_public_report_url"       ""
# 26 "submitted_as_notes"                   ""
# 27 "submitted_as_pmids"                   ""
# 28 "submitted_as_assertion_criteria_url"  "PMID: 28106320"
# 29 "submitted_as_submission_id"           "1034"
# 30 "submitted_run_date"                   "2020-12-24"
#
# 列の説明はホームページで見つけられなかった
# なぜかGenCCのファイルはcutやlessだと正しく見られなくて、pd.read_table だと正しく見られる
# awkでもある程度正しく見られるが少し変
# moi_titleの内訳：
# - 'Autosomal dominant': 8913
# - 'Autosomal recessive': 10132
# - 'Semidominant': 139
# - 'X-linked': 1004
# - 'X-linked recessive': 2
# - 'Y-linked inheritance': 2
# - 'Mitochondrial': 98
# - 'Unknown': 928


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


parser = argparse.ArgumentParser(description='')
parser.add_argument('--urls_file', type=str, required=True, help='URLが1行ずつ書かれたファイルのパス')
parser.add_argument('--out',       type=str, default=default_filename)
args = parser.parse_args()


# 与えたURLからファイルをダウンロード
with open( args.urls_file, 'r' ) as f:
  urls = [ line.strip() for line in f if line.strip() ]

download_from_URLs( urls )


# 複数のファイルを縦に結合する
files = extract_filenames_from_urls( urls )

df = make_concatenated_df( files, sep='\t' )


# 整形して保存する
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
  rename( {'GenCC_gene_symbol': 'Gene.refGene'}, axis = 1 ).\
  to_csv( args.out, sep = '\t', index = False )

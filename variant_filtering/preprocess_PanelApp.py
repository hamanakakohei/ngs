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
from general import make_concatenated_df, extract_filenames_from_urls 
from datetime import date


today_str = date.today().strftime('%Y_%m_%d')  # 例: '2025_05_16'
default_filename = f'PanelApp_clean_{today_str}.tsv.gz'


# ステップ０
parser = argparse.ArgumentParser(description='')
parser.add_argument('--PanelApp_files', type=str, required=True, help='ダウンロードしたPanelAppファイルが1行ずつ書かれたファイルのパス')
parser.add_argument('--out',       type=str, default=default_filename)
args = parser.parse_args()

# PanelAppファイルの列の説明（https://panelapp.genomicsengland.co.uk/#!Navigating）：
# - Entity Name (gene name, STR name, region name)
# - Entity Type (gene, str, region)
# - Gene Symbol (HGNC-approved gene symbol)
# - Sources (separated by ;)
# - Level 4, 3 and 2 titles (E.g. specific rare disease specific, subgroup and group)
# - Mode of inheritance (using PanelApp standardised terms)
# - Phenotypes (as collected from all sources and reviews)
# - OMIM, Orphanet, HPO terms (not currently populated)
# - Publications
# - Description (not currently populated)
# - Flagged (default = FALSE)
# - GEL status (reflects the current gene rating; 3 or 4 = Green, 2 = Amber, 1 or 0 = Red)
# - User Ratings (represented as the percentage of ratings for Green;Amber;Red).
# - Version (the gene panel version)
# - Ready
# - Mode of pathogenicity
# - Ensembl ID (GRCh37)
# - Ensembl ID (GRCh38)
# - HGNC ID
# 
# For STRs only, the following fields will be populated:
# - Position Chromosome
# - Position GRCh37 start
# - Position GRCh37 end
# - Position GRCh38 start
# - Position GRCh38 end
# - STR Repeated Sequence
# - STR Normal Repeats
# - STR Pathogenic Repeats
# 
# For CNVs only, the following fields will be populated:
# - Region Haploinsufficiency Score
# - Region Triplosensitivity Score
# - Region Required Overlap Percentage
# - Region Variant Type
# - Region Verbose Name

# - Gene Symbol (HGNC-approved gene symbol)
# - HGNC ID
# - Phenotypes (as collected from all sources and reviews)
# - Position GRCh37 start
# - Position GRCh37 end
# - STR Repeated Sequence
# - STR Normal Repeats
# - STR Pathogenic Repeats

# ステップ１：PanelAppファイルを縦に結合する
with open( args.PanelApp_files, 'r' ) as f:
  files = [ line.strip() for line in f if line.strip() ]

df = make_concatenated_df( files, sep='\t' )


# ステップ２：色々整形する
df.columns = df.columns.str.replace(' ', '_', regex=False)
df.columns = [ f"PanelApp_{col}" for col in df.columns ]

interest_cols = [
  'PanelApp_Gene_Symbol', 
  'PanelApp_Phenotypes',
  'PanelApp_HGNC', 
  'PanelApp_Model_Of_Inheritance',
  'PanelApp_Mode_of_pathogenicity',
  'PanelApp_Position_GRCh37_Start',
  'PanelApp_Position_GRCh37_End', 
  'PanelApp_STR_Repeated_Sequence', 
  'PanelApp_STR_Normal_Repeats', 
  'PanelApp_STR_Pathogenic_Repeats'
]

group_cols = [ 'PanelApp_Gene_Symbol', 'PanelApp_HGNC' ]
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
  rename( {'PanelApp_Gene_Symbol': 'gene'}, axis = 1 ).\
  to_csv( args.out, sep = '\t', index = False )


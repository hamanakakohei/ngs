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


# ステップ０
parser = argparse.ArgumentParser(description='')
parser.add_argument('--urls_file', type=str, required=True, help='URLが1行ずつ書かれたファイルのパス')
parser.add_argument('--out',       type=str, default='G2P_merged_clean_2025-04-28.tsv.gz')
args = parser.parse_args()


# ステップ１
# 以下のURLからG2Pデータをダウンロード
# G2Pデータフォーマット：https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/DataDownloadFormat202501.txt
#The latest G2P records can be downloaded on the fly from https://www.ebi.ac.uk/gene2phenotype/beta/download.
#
#The G2P csv data release format describes a single G2P entry on each row.
#
#The columns in the file are:  
#
#  - g2p id:                               the stable identifier for this record in G2P
#  - gene symbol:                          the HGNC-assigned gene symbol 
#  - gene mim:                             the OMIM identifier for the gene
#  - hgnc id:                              the HGNC-assigned identifier for the  gene 
#  - previous gene symbols:                previous gene symbols, where relevant
#  - disease name:                         the name for the disease used in G2P 
#  - disease mim:                          the OMIM identifier for the same/a highly similar disease, where available 
#  - disease MONDO:                        the MONDO identifier for the same/a highly similar disease, where available
#  - allelic requirement:                  described using synonyms of HPO mode of inheritance term 
#  - cross cutting modifier:               additional qualifiers 
#  - confidence:                           an assertion of the confidence the association is real, based on the evidence available
#  - variant consequence:                  the consequence(s) of the reported variants, inferred at protein or RNA level
#  - variant types:                        the variant types reported in the publication to be associated with the disease
#  - molecular mechanism:                  the disease mechanism suggested or reported by the publication 
#  - molecular mechanism categorisation:   a more specific categorisation of the disease mechanism
#  - molecular mechanism evidence:         the type(s) of evidence available to support the reported mechanism
#  - phenotypes:                           HPO identifiers for the reported phenotypes
#  - publications:                         the publications reviewed to create this record
#  - panel:                                the G2P panel(s) the record is assigned to
#  - comments:                             additional comments from the curation team ( where available)
#  - date of last review                   date the record was modified/reviewed 
#
#See : https://www.ebi.ac.uk/gene2phenotype/beta/about/terminology for further details of the terminology used

with open( args.urls_file, 'r' ) as f:
  urls = [ line.strip() for line in f if line.strip() ]

download_from_URLs( urls )


# ステップ２：全ファイルを縦に結合する
files = extract_filenames_from_urls( urls )

df = make_concatenated_df( files, quotechar='"' )


# ステップ３：列名を綺麗にしたり選んだり整形する
df.columns = df.columns.str.replace(' ', '_', regex=False)
df.columns = [ f"G2P_{col}" for col in df.columns ]

interest_cols = [
  'G2P_gene_symbol',
  'G2P_hgnc_id',
  'G2P_previous_gene_symbols',
  'G2P_disease_name',
  'G2P_allelic_requirement',
  'G2P_cross_cutting_modifier',
  'G2P_confidence',
  'G2P_variant_consequence',
  'G2P_variant_types',
  'G2P_molecular_mechanism',
  'G2P_molecular_mechanism_categorisation',
  'G2P_molecular_mechanism_evidence'
]
group_cols = [ 'G2P_gene_symbol', 'G2P_hgnc_id' ]
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
  rename( {'G2P_gene_symbol': 'gene'}, axis = 1 ).\
  to_csv( args.out, sep = '\t', index = False )


#!/usr/bin/env python3


# 目的：G2Pが提供するを、遺伝子名をキーにして他の表と結合できるように整形する
# ただしVEPならG2Pプラグインで同様のことを自動でできそうだ
# VEPが対応していないバリアントフォーマットを出力するツール用になら役立つかもしれない


import pandas as pd
import sys
import os
script_dir = os.path.dirname( os.path.abspath(__file__) )
general_path = os.path.abspath( os.path.join(script_dir, "../../misc/utils") )
sys.path.append( general_path )
from general import download_from_URLs


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

urls = [
  "https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/2025_04_28/CardiacG2P_2025-04-28.csv.gz",
  "https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/2025_04_28/DDG2P_2025-04-28.csv.gz",
  "https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/2025_04_28/EyeG2P_2025-04-28.csv.gz",
  "https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/2025_04_28/HearingLossG2P_2025-04-28.csv.gz",
  "https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/2025_04_28/SkeletalG2P_2025-04-28.csv.gz",
  "https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/2025_04_28/SkinG2P_2025-04-28.csv.gz"
]

download_from_URLs( urls )


# ステップ２
df = pd.concat([
    pd.read_csv('CardiacG2P_2025-04-28.csv.gz',       quotechar='"'),
    pd.read_csv('DDG2P_2025-04-28.csv.gz',            quotechar='"'),
    pd.read_csv('EyeG2P_2025-04-28.csv.gz',           quotechar='"'),
    pd.read_csv('HearingLossG2P_2025-04-28.csv.gz',   quotechar='"'),
    pd.read_csv('SkeletalG2P_2025-04-28.csv.gz',      quotechar='"'),
    pd.read_csv('SkinG2P_2025-04-28.csv.gz',          quotechar='"')
  ]).\
  drop_duplicates()\
  [[
    'gene symbol',
    'hgnc id',
    'previous gene symbols',
    'disease name',
    'allelic requirement',
    'cross cutting modifier',
    'confidence',
    'variant consequence',
    'variant types',
    'molecular mechanism',
    'molecular mechanism categorisation',
    'molecular mechanism evidence',
    'phenotypes'
  ]].\
  to_csv( 'G2P_2025-04-28.tsv', sep = '\t', index = False )

#! /usr/bin/Rscript


library(tidyverse)
library(argparser)
source("misc/utils/annovar.R")


argv = arg_parser("") %>%
  add_argument( "--annovar_result"        , type="character", required=TRUE,  default="", help="" ) %>%
  add_argument( "--af_threshold_hgvd"     , type="numeric",   required=FALSE, default="", help="" ) %>%
  add_argument( "--af_threshold_tommo"    , type="numeric",   required=FALSE, default="", help="" ) %>%
  add_argument( "--af_threshold_exac_all" , type="numeric",   required=FALSE, default="", help="" ) %>%
  add_argument( "--af_threshold_exac_eas" , type="numeric",   required=FALSE, default="", help="" ) %>%
  add_argument( "--inheritance"           , type="character", required=TRUE,  default="", 
  add_argument( "--gene_mode_of_inheritance_filter"         , required=FALSE, default="", flag=TRUE,
  add_argument( "--gene_annotations"      , type="character", required=FALSE, default="", 
  add_argument( "--out"                   , type="character", required=TRUE,  default="", help="" ) %>%
  parse_args()
 # add_argument( "--count_threshold_jpn"   , type="integer",   required=FALSE, default="", help="" ) %>%



./ngs/variant_handling/filter_annovar_result.R \
  --annovar_result /mnt/c/Users/hamanaka/Desktop/Scrivener/Exomeデータ/ataxia61/exome_summary.20220603_192350.txt.gz \
  --af_threshold_hgvd     0.001 \
  --af_threshold_tommo    0.001 \
  --af_threshold_exac_all 0.001 \
  --af_threshold_exac_eas 0.001 \
  --inheritance           AD \
  --gene_annotations G2P_merged_2025-04-28.tsv.gz,GenCC_clean_2025_05_16.tsv.gz \
  --out tmp.AD.txt



#  --gene_mode_of_inheritance_filter

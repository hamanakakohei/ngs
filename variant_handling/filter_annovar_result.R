#!/home/hamanaka/miniconda3/envs/misc/bin/Rscript


library(tidyverse)
library(argparser)
source("misc/utils/annovar.R")


argv = arg_parser("") %>%
  add_argument( "--annovar_result"        , type="character",  help="" ) %>%
  add_argument( "--af_threshold_hgvd"     , type="numeric",    help="値を指定しなければNAになり何も起きない" ) %>%
  add_argument( "--af_threshold_tommo"    , type="numeric",    help="値を指定しなければNAになり何も起きない" ) %>%
  add_argument( "--af_threshold_exac_all" , type="numeric",    help="値を指定しなければNAになり何も起きない" ) %>%
  add_argument( "--af_threshold_exac_eas" , type="numeric",    help="値を指定しなければNAになり何も起きない" ) %>%
  add_argument( "--inheritance"           , type="character",   
    help=paste(
      "AD,AR,XLのどれかを指定する",
      "ARを指定すると、アレル頻度でフィルター後に2個以上バリアントがある遺伝子を返す",
      sep = "\n"
    )
  ) %>%
  add_argument( "--gene_mode_of_inheritance_filter"         ,  flag=TRUE,
    help=paste(
      "OMIM/G2P/GenCCのMOIでフィルターする",
      "MOI列名は勝手にこっちで想定している",
      "発端者のみを解析している時向け",
      "トリオ解析時は新規の遺伝形式も考慮してこのフィルターは外す",
      sep = "\n"
    )
  ) %>%
  add_argument( "--gene_annotations"      , type="character", 
    help=paste(
      "遺伝子に対するアノテーションファイルを指定する",
      "カンマかスペース区切りで複数ファイルを指定できる",
      "OMIM、G2P、GenCC、候補遺伝子リストなどを指定する",
      "遺伝子列名はgeneに統一する",
      sep = "\n"
    ), 
    nargs=Inf ) %>%
  add_argument( "--out"                   , type="character",  help="" ) %>%
  parse_args()
 # add_argument( "--count_threshold_jpn"   , type="integer",   default="", help="" ) %>%


# Annovarの結果を整形する
df = read_tsv( argv$annovar_result, col_types=cols( Chr="c", CHROM="c" ) ) %>%
  mutate(
    maf_exac_all             = as.numeric( ifelse( ExAC_ALL==".", "0", ExAC_ALL ) ),
    maf_exac_eas             = as.numeric( ifelse( ExAC_EAS==".", "0", ExAC_EAS ) ),
    maf_hgvd                 = map_dbl( HGVD, get_hgvd_maf ),
    maf_tommo                = map_dbl( snp20171005_tommo3.5k_passed, get_tommo_maf ),
    JpnMutation_symbol       = map_chr( INFO, get_JpnMutation_symbol ),
    JpnMutation_count        = map_dbl( INFO, get_JpnMutation_count  ),
    JpnMutation_count_male   = map_dbl( INFO, ~get_JpnMutation_sex_chr( .x, SEX="male"   ) ),
    JpnMutation_count_female = map_dbl( INFO, ~get_JpnMutation_sex_chr( .x, SEX="female" ) )
  ) %>%
  unite( "variant_id", c(Chr,Start,End,Ref,Alt), sep="_" ) 


# キャリアステータス（其バリアントを持っているか？その遺伝子に2hit以上もつか？）を加える 
list2env( get_variant_carriers( DF ), envir = .GlobalEnv )

df = left_join( df, variant_carriers, by="variant_id" ) %>%
  left_join( gene_carriers_2hit, by="Gene.refGene" )

 
# アレル頻度でフィルターする、NAだとしない
df = df %>%
  filter_af_threshold("maf_hgvd",     argv$af_threshold_hgvd    ) %>%
  filter_af_threshold("maf_tommo",    argv$af_threshold_tommo   ) %>%
  filter_af_threshold("maf_exac_all", argv$af_threshold_exac_all) %>%
  filter_af_threshold("maf_exac_eas", argv$af_threshold_exac_eas)


## 持っているバリアントのみ
#df = filter( df, !is.na(carriers) )
# ARなら2hit以上の遺伝子のみにする
if( argv$inheritance == "AR" ){
  df = filter( df, !is.na( carriers_2hit_on_the_gene ) )
}


# 遺伝子名に対するアノテーションをあるだけ付ける
for( gene_annotation in argsv$gene_annotations ){
  gene_annotation_df = read_tsv( gene_annotation )
  df = left_join( df, gene_annotation_df, by="gene" )
}


# 指定した遺伝形式に合う染色体のみにする
df = filter_chr_by_inheritance( df, argv$inheritance )


# 遺伝子の既知遺伝形式でフィルターする
if( argv$gene_mode_of_inheritance_filter ){
  df = filter_by_gene_mode_of_inheritance( df, argv$inheritance )
}


# 保存する
write_tsv( df, argv$out, col_names=TRUE )





#if ( !is.null( argv$count_threshold_jpn ) ) {
#  df = filter( df, JpnMutation_count < argv$count_threshold_jpn )
#}
#
#EXOMESUMMARY_PATH           = argv$exomesummary
#GENE_CANDIDATE.STATUS_PATH  = argv$gene_candidate 
#GENE_OMIM.STATUS_PATH       = argv$gene_omim
#HGVD.THR	                = as.numeric(argv$hgvd)     #0.0001
#TOMMO.THR	                = as.numeric(argv$tommo)    #0.0001
#JPNCOUNT.THR                = as.numeric(argv$jpncount) #0
#EXACALL.THR	                = as.numeric(argv$exacall)  #0.0001
#EXACEAS.THR                 = as.numeric(argv$exaceas)  #0.0001
#INHERITANCE                        = argv$inheritance  #"AD" # AR, XL
#OUT                         = argv$out          #"output.txt"
#
#GENE_CANDIDATE.STATUS_PATH="/archive3/hamanaka/resource/gene_candidate.status__for.SCD.txt"
#GENE_OMIM.STATUS_PATH="/archive3/hamanaka/resource/omim.genemap2.edit.txt"
##EXOMESUMMARY_PATH="/betelgeuse01/analysis/191213_156batch_NovaSeq_singleindex_analysis/Project_SCD_shinshu_sporadic/annovar/exome_summary.20191225_144612.txt"
#EXOMESUMMARY_PATH="/betelgeuse04/analysis/hamanaka/jointgt.13851/scos/Sample_14222/annovar/exome_summary.20200914_202707.txt"
#HGVD.THR	                = 0.01
#TOMMO.THR	                = 0.01
#JPNCOUNT.THR                = 5
#EXACALL.THR	                = 0.01
#EXACEAS.THR                 = 0.01
#TYPE                        = "AD" #, XL
#OUT                         = "output.ar.2308.txt"
#DDG2P="/betelgeuse04/analysis/hamanaka/resource/gene_inhe__DDG2P.13.9.2020mod.txt"
#DDG2P="/betelgeuse07/analysis/hamanaka/resource/gene_inhe__DDG2P.13.9.2020mod.txt"

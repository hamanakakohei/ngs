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
  add_argument( "--other_caller_results"   , type="character", 
    help=paste(
      "別のコーラー結果（例：CNV）のファイルを指定する",
      "カンマかスペース区切りで複数ファイルを指定できる",
      "事前に整形が必要で、フォーマットはcat_another_caller_variants関数を参照",
      sep = "\n"
    ), 
    nargs=Inf ) %>%
  add_argument( "--out"                   , type="character",  help="" ) %>%
  parse_args()
 # add_argument( "--count_threshold_jpn"   , type="integer",   default="", help="" ) %>%


# Annovar結果を整形する
df = read_tsv( argv$annovar_result, col_types=cols( Chr="c", CHROM="c" ) ) %>%
  mutate(
    maf_exac_all             = as.numeric( ifelse( ExAC_ALL==".", "0", ExAC_ALL ) ),
    maf_exac_eas             = as.numeric( ifelse( ExAC_EAS==".", "0", ExAC_EAS ) ),
    maf_hgvd                 = map_dbl( HGVD, get_hgvd_maf ),
    maf_tommo                = map_dbl( snp20171005_tommo3.5k_passed, get_tommo_maf ),
    JpnMutation_symbol       = map_chr( INFO, get_JpnMutation_symbol ),
    JpnMutation_count        = map_dbl( INFO, get_JpnMutation_count  ),
    JpnMutation_count_male   = map_dbl( INFO, ~get_JpnMutation_count_sex_chr( .x, SEX="male"   ) ),
    JpnMutation_count_female = map_dbl( INFO, ~get_JpnMutation_count_sex_chr( .x, SEX="female" ) ),
    SpliceAI_max_score       = map_dbl( SpliceAI, get_SpliceAI_max_score )
  ) %>%
  unite( "variant_id", c(Chr,Start,End,Ref,Alt), sep="_", remove=FALSE ) 


# 興味のないバリアントを除く
## １：SpliceAIスコアが低いシノニマス
filter1 = with( df, 
  ExonicFunc.refGene == "synonymous SNV" & 
  (is.na(SpliceAI_max_score) | SpliceAI_max_score < 0.1 ) 
)

## ２：コーディング外でUTRでもなくSpliceAIスコアも低いもの
filter2 = with( df, 
  ExonicFunc.refGene == "." & 
  !str_detect(GeneDetail.refGene, "UTR") & 
  SpliceAI_max_score < 0.1 
)

## 以上をまとめて除外する
exclude = filter1 | filter2
df = filter( df, !exclude )


# アレル頻度でフィルターする、NAだとしない
df = df %>%
  filter_af_threshold("maf_hgvd",     argv$af_threshold_hgvd    ) %>%
  filter_af_threshold("maf_tommo",    argv$af_threshold_tommo   ) %>%
  filter_af_threshold("maf_exac_all", argv$af_threshold_exac_all) %>%
  filter_af_threshold("maf_exac_eas", argv$af_threshold_exac_eas)


# Annovar以外のコーラー結果をあるだけ上下に結合する
if( any( !is.na( argv$other_caller_results ) ) ){
  for( other_caller_result in argv$other_caller_results ){
    other_caller_result_df = read_tsv( other_caller_result, col_types=cols( Chr="c" ) )
    df = cat_another_caller_variants( df, other_caller_result_df )
  }
}


# キャリアステータス（其バリアントを持っているか？その遺伝子に2hit以上もつか？）を加える 
list2env( get_variant_carriers( df ), envir = .GlobalEnv )

df = left_join( df, variant_carriers, by="variant_id" ) %>%
  left_join( gene_carriers_2hit, by="Gene.refGene" ) 


## 持っているバリアントのみ
#df = filter( df, !is.na(carriers) )
# ARなら2hit以上の遺伝子のみにする
if( argv$inheritance == "AR" ){
  df = filter( df, !is.na( carriers_2hit_on_the_gene ) )
}


# 遺伝子名に対するアノテーションをあるだけ付ける
if( any( !is.na( argv$gene_annotations ) ) ){
  for( gene_annotation in argv$gene_annotations ){
    gene_annotation_df = read_tsv( gene_annotation )
    df = left_join( df, gene_annotation_df, by="Gene.refGene" )
  }
}


# 指定した遺伝形式に合う染色体のみにする
df = filter_chr_by_inheritance( df, argv$inheritance )


# 遺伝子の既知遺伝形式でフィルターする
if( argv$gene_mode_of_inheritance_filter ){
  df = filter_by_gene_mode_of_inheritance( df, argv$inheritance )
}


# 保存する
write_tsv( df, argv$out, col_names=TRUE )

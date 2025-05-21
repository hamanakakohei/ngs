#!/home/hamanaka/miniconda3/envs/misc/bin/Rscript


library(tidyverse)
library(argparse)
library(glue)
source("misc/utils/annovar.R")

parser = ArgumentParser(description = "Annovarと他コーラーのバリアント表を統合しつつ、アレル頻度、遺伝形式、対象サンプルなどの条件でフィルターする。")
parser$add_argument("--annovar_result",        type = "character", required = TRUE, help = "Annovar結果のテーブルファイル")
parser$add_argument("--out",                   type = "character", required = TRUE, help = "出力ファイル名")
parser$add_argument("--af_threshold_hgvd",     type = "double",    help = "HGVDのアレル頻度の上限値。指定しない場合はフィルターしない")
parser$add_argument("--af_threshold_tommo",    type = "double",    help = "ToMMoのアレル頻度の上限値。指定しない場合はフィルターしない")
parser$add_argument("--af_threshold_exac_all", type = "double",    help = "ExAC ALLのアレル頻度の上限値。指定しない場合はフィルターしない")
parser$add_argument("--af_threshold_exac_eas", type = "double",    help = "ExAC EASのアレル頻度の上限値。指定しない場合はフィルターしない")
parser$add_argument("--inheritance",           type = "character", required = TRUE, help = paste(
                      "出力したい遺伝形式（AD,AR,XLのどれか）を指定する",
                      "AD,ARなら常染色体のみ、XLならX染色体のバリアントに絞る",
                      "この指定に基づき、--sample_filter and/or --gene_mode_of_inheritance_filterが実行される",
                      sep = "\n"))
parser$add_argument("--sample_filter",        type = "character", help = paste(
                      "指定したサンプルが保持するバリアントのみを出力",
                      "--inheritanceでARを指定した場合、そのサンプルで2つ以上のバリアントがある遺伝子に絞り込まれる",
                      sep = "\n"))
parser$add_argument("--gene_mode_of_inheritance_filter", action = "store_true", help = paste(
                      "指定された遺伝形式（--inheritance）と一致する遺伝子のバリアントのみを残す",
                      "GenCC や G2P 由来の MOI (mode of inheritance) 情報を使用",
                      "MOI列の名前は内部で想定されている（詳細は filter_by_gene_mode_of_inheritance 関数を参照）",
                      "トリオ解析時など新規の遺伝形式を探索する時は、このフィルターを外す",
                      sep = "\n"))
parser$add_argument("--gene_annotations",      type = "character", nargs = "+", help = paste(
                      "GenCC、G2P、候補遺伝子リストなどの遺伝子に対するアノテーションファイルを指定",
                      "複数ファイルの指定はカンマかスペース区切り",
                      "各ファイルは 'Gene.refGene' という列に遺伝子名を記載する",
                      sep = "\n"))
parser$add_argument("--other_caller_results",  type = "character", nargs = "+", help = paste(
                      "別のコーラー（例：XHMM）の出力ファイルを指定する、カンマかスペース区切りで複数指定可",
                      "事前に整形が必要で、XHMMの場合はpreprocess_XHMM.pyを使う",
                      "その他のコーラーの場合、フォーマットはcat_another_caller_variants関数を参照",
                      sep = "\n"))
argv = parser$parse_args()


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


# 指定したサンプルが持っているバリアントのみ、ARなら2hit以上の遺伝子のみ
pattern = glue("{argv$sample_filter};")

if( !is.na( argv$sample_filter ) ){
  df = filter( df, filter( str_detect(carriers, pattern) ) )
  if( argv$inheritance == "AR" ){
    df = filter( df, filter( str_detect(carriers_2hit_on_the_gene, pattern) ) )
  }
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

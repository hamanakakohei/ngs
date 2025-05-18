#! /usr/bin/Rscript


library(tidyverse)
library(argparser)


arg_parser("") %>%
    add_argument("--exomesummary"   ,help="") %>%
    add_argument("--gene_candidate" ,help="") %>%
    add_argument("--gene_omim"      ,help="",default="/betelgeuse07/analysis/hamanaka/resource/omim.genemap2.edit.txt") %>%
    add_argument("--hgvd"           ,help="") %>%
    add_argument("--tommo"          ,help="") %>%
    add_argument("--jpncount"       ,help="") %>%
    add_argument("--exacall"        ,help="") %>%
    add_argument("--exaceas"        ,help="") %>%
    add_argument("--inheritance"    ,help="") %>%
    add_argument("--out"            ,help="") %>%
    parse_args() -> argv

#REMOVE.TYPE     = c("synonymous SNV","unknown")
EXOMESUMMARY_PATH           = argv$exomesummary
GENE_CANDIDATE.STATUS_PATH  = argv$gene_candidate 
GENE_OMIM.STATUS_PATH       = argv$gene_omim
HGVD.THR	                = as.numeric(argv$hgvd)     #0.0001
TOMMO.THR	                = as.numeric(argv$tommo)    #0.0001
JPNCOUNT.THR                = as.numeric(argv$jpncount) #0
EXACALL.THR	                = as.numeric(argv$exacall)  #0.0001
EXACEAS.THR                 = as.numeric(argv$exaceas)  #0.0001
INHERITANCE                        = argv$inheritance  #"AD" # AR, XL
OUT                         = argv$out          #"output.txt"

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
DDG2P="/betelgeuse07/analysis/hamanaka/resource/gene_inhe__DDG2P.13.9.2020mod.txt"



gene_GenCC = read_tsv()
gene_G2P
gene_OMIM



read_tsv(GENE_CANDIDATE.STATUS_PATH,col_names=c("gene","candidate.status")) -> gene_candidate.status  

read_tsv( EXOMESUMMARY_PATH, col_types=cols( Chr="c", CHROM="c" ) ) %>%
  mutate(
    #dist                     = map_chr( HGMD_AllMut_collapse, extract_dist ),
    ExAC_ALL                 = as.numeric( ifelse( ExAC_ALL==".", "0", ExAC_ALL) ),
    ExAC_EAS                 = as.numeric( ifelse( ExAC_EAS==".", "0", ExAC_EAS) ),
    hgvd.maf                 = map_dbl( HGVD, extract_hgvd.maf ),
    tommo.maf                = map_dbl( snp20171005_tommo3.5k_passed, extract_tommo.maf ),
    JpnMutation_symbol       = map_chr( INFO, extract_jpn.symbol ),
    JpnMutation_count        = map_dbl( INFO, extract_jpn.count  ),
    JpnMutation_count_male   = map_dbl( INFO, ~extract_jpn.sex( .x, SEX="male"   )),
    JpnMutation_count_female = map_dbl( INFO, ~extract_jpn.sex( .x, SEX="female" ))
  ) %>%
  filter_by_maf(HGVD.THR,TOMMO.THR,JPNCOUNT.THR,EXACALL.THR,EXACEAS.THR) %>%
  filter_by_inheritance( INHERITANCE ) %>%
  unite("variant.id",c(Chr,Start,End,Ref,Alt),sep="_") %>% {
    if( INHERITANCE %in% c("AD","XL") ){
      extract_qualifed_gene.or.variant( ., type="1hit" ) ->> variant.id_carriers
      left_join(.,variant.id_carriers)
    }else if( INHERITANCE == "AR" ){
      extract_qualifed_gene.or.variant( ., type="1hit" ) ->> variant.id_carriers
      extract_qualifed_gene.or.variant( ., type="2hit" ) ->> Gene.refGene_carriers.2hit
      left_join(.,variant.id_carriers) %>%
        left_join(Gene.refGene_carriers.2hit)
    }
  } %>% 
  left_join( gene_inh__ddg2p,       by=c("Gene.refGene"="gene") ) %>%
  left_join( gene_candidate.status, by=c("Gene.refGene"="gene") ) %>% 
  left_join( gene_omim.status,      by=c("Gene.refGene"="gene") ) %>% 
  #left_join(CHROM_POS_REF_ALT_spai) %>% 
  #rename(HGMD_AllMut_collapse='HGMD_AllMut-collapse(HGMD-2021.3)') %>%
  #filter(!ExonicFunc.refGene %in% REMOVE.TYPE) %>%
  write_tsv( OUT, col_names=TRUE )



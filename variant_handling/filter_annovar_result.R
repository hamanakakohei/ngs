#! /usr/bin/Rscript
library(tidyverse)
library(argparser)

arg_parser("funuuu") %>%
    add_argument("--exomesummary"   ,help="funuu") %>%
    add_argument("--gene_candidate" ,help="funuu") %>%
    add_argument("--gene_omim"      ,help="funuu",default="/betelgeuse07/analysis/hamanaka/resource/omim.genemap2.edit.txt") %>%
    add_argument("--hgvd"           ,help="funuu") %>%
    add_argument("--tommo"          ,help="funuu") %>%
    add_argument("--jpncount"       ,help="funuu") %>%
    add_argument("--exacall"        ,help="funuu") %>%
    add_argument("--exaceas"        ,help="funuu") %>%
    add_argument("--inheritance"    ,help="funuu") %>%
    add_argument("--out"            ,help="funuu") %>%
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
TYPE                        = argv$inheritance  #"AD" # AR, XL
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

extract_qualifed_gene.or.variant = function(DT,type="1hit"){
    select(DT,c(variant.id,Gene.refGene,starts_with("Sample_"))) %>% 
        gather(key=sample,value=genotype,starts_with("Sample_")) %>%
        mutate(allele.count=case_when(
            str_detect(.$genotype, "^0/1") ~ 1,
            str_detect(.$genotype, "^1/0") ~ 1,
            str_detect(.$genotype, "^1/1") ~ 2,
            !str_detect(.$genotype, "^0/0") & str_detect(.$genotype, "^0/") ~ 1,
            !str_detect(.$genotype, "^0/0") & str_detect(.$genotype, "^./0") ~ 1,
            !str_detect(.$genotype, "^0/0") & !str_detect(.$genotype, "^0/") & !str_detect(.$genotype, "^./0") ~ 2,
            TRUE ~ 0)
        ) %>% {
            if(type=="1hit"){
                group_by(.,variant.id,sample) %>% summarise(total.ac=sum(allele.count)) %>% 
                filter(total.ac>=1) %>% nest(-variant.id) %>%
                mutate(carriers=map(data,connect_samples)) %>% unnest(carriers) %>%
                select(variant.id,carriers) %>% return
            }else if(type=="2hit"){
                group_by(.,Gene.refGene,sample) %>% summarise(total.ac=sum(allele.count)) %>% 
                filter(total.ac>=2) %>% nest(-Gene.refGene) %>% 
                mutate(carriers.2hit=map(data,connect_samples)) %>% unnest(carriers.2hit) %>%
                select(Gene.refGene,carriers.2hit) %>% return
            }
        }
}

connect_samples = function(DT.SAMPLES){
    DT.SAMPLES %>% pull(sample) %>% str_c(collapse=";") %>% return
}

extract_dist = function(DIST){
    strsplit(DIST,"\\(") -> tmp
    tmp[[1]][2] %>% strsplit(")") -> tmp2
    tmp2[[1]] %>% str_replace("dist=","") -> tmp3; tmp3[1] %>% return
}

filter_maf1 = function(DT,HGVD2,TOMMO,JPNCOUNT,EXACALL,EXACEAS){
    filter(DT,
		hgvd.maf    <= HGVD2 &
		tommo.maf   <= TOMMO &
		ExAC_ALL    <= EXACALL &
		ExAC_EAS    <= EXACEAS &
		((jpn.symbol!="p" & jpn.count<=JPNCOUNT) | jpn.symbol=="p")) %>% return
}

filter_maf2 = function(DT,TYPE){
    if(TYPE=="AD"){
        filter(DT,!Chr %in% c("X","Y")) %>% return
    }else if(TYPE=="AR"){
        filter(DT,!Chr %in% c("X","Y")) %>% {
            return(.)
            #group_by(.,Gene.refGene) %>% summarise(n=sum(allele.count)) %>% filter(n>1) %>% pull(Gene.refGene) ->> ar.gene
            #filter(.,Gene.refGene %in% ar.gene) %>% return
        }
    }else if(TYPE=="XL"){
        filter(DT,Chr %in% c("X")) %>% return
    }
}


extract_hgvd.maf = function( HGVD_col ){
	if( HGVD_col == "." ){
		return( 0 )
	}else{
		HGVD_AF = strsplit( HGVD_col, "," )[[1]][1]
		if (startsWith(HGVD_AF, "AF=")) {
			HGVD_AF = as.numeric( str_replace(HGVD_AF, "AF=", "") )
			return( HGVD_AF )
		} else {
	  		stop("エラー: HGVD_AF は 'AF=' で始まっていません。")
		}
	}
}


get_jpn.list = function( INFO ){
	INFO_v = strsplit( INFO, ";" )[[1]]
	JPN = INFO_v[ length( INFO_v ) ]
	if ( startsWith( JPN, "JpnMutation=" ) ) {
		JPN_l = str_replace("JpnMutation=", "") %>%
	       		strsplit(":") 
		return( JPN_l )
	} else {
		stop("エラー: JPN は 'JpnMutation=' で始まっていません。")
	}
}

extract_jpn.symbol = function(INFO) {
	tmp2 = get_jpn.list(INFO)
	return(tmp2[[1]][1])
}	

extract_jpn.count = function(INFO) {
	get_jpn.list(INFO) -> tmp2; tmp2[[1]][3] -> COUNT;
	if(str_detect(COUNT,",")){
		strsplit(COUNT,",")[[1]] %>% as.numeric %>% sum %>% return
	}else{return(as.numeric(COUNT))}
}	

extract_jpn.sex = function(INFO,SEX="male") {
	get_jpn.list(INFO) -> tmp2; tmp2[[1]][3] -> COUNT;
	if(str_detect(COUNT,",")){
		strsplit(COUNT,",")[[1]] %>% as.numeric -> tmp;
		if(SEX=="male"){return(as.numeric(tmp[1]))}else{as.numeric(return(tmp[2]))}
	}else{return(NA)}
}	


extract_tommo.maf = function(TOMMO){
	if( TOMMO == "." ){
		return(0)
	}else{
		TOMMO_AF = strsplit( TOMMO, ":" )[[1]][3] %>% 
			as.numeric

		return( TOMMO_AF )
	}
}


gene_GenCC = read_tsv()
gene_G2P
gene_OMIM



read_tsv(GENE_CANDIDATE.STATUS_PATH,col_names=c("gene","candidate.status")) -> gene_candidate.status  

read_tsv(EXOMESUMMARY_PATH,col_types=cols(Chr="c",CHROM="c")) %>%
	left_join(gene_candidate.status,by=c("Gene.refGene"="gene")) %>% 
	left_join(gene_omim.status,by=c("Gene.refGene"="gene")) %>% 
	left_join(CHROM_POS_REF_ALT_spai) %>% 
    rename(HGMD_AllMut_collapse='HGMD_AllMut-collapse(HGMD-2021.3)') %>%
    #filter(!ExonicFunc.refGene %in% REMOVE.TYPE) %>%
	mutate(
		dist=map_chr(HGMD_AllMut_collapse,extract_dist),
        ExAC_ALL=as.numeric(ifelse(ExAC_ALL==".","0",ExAC_ALL)),
		ExAC_EAS=as.numeric(ifelse(ExAC_EAS==".","0",ExAC_EAS)),
		hgvd.maf=map_dbl(HGVD,extract_hgvd.maf),
		tommo.maf=map_dbl(snp20171005_tommo3.5k_passed,extract_tommo.maf),
		jpn.symbol=map_chr(INFO,extract_jpn.symbol),
		jpn.count=map_dbl(INFO,extract_jpn.count),
		jpn.count.male=map_dbl(INFO,~extract_jpn.sex(.x,SEX="male")),
		jpn.count.female=map_dbl(INFO,~extract_jpn.sex(.x,SEX="female"))) %>%
    filter_maf1(HGVD.THR,TOMMO.THR,JPNCOUNT.THR,EXACALL.THR,EXACEAS.THR) %>%
    filter_maf2(TYPE) %>%
    unite("variant.id",c(Chr,Start,End,Ref,Alt),sep="_") %>% {
        if(TYPE %in% c("AD","XL")){
            extract_qualifed_gene.or.variant(.,type="1hit") ->> variant.id_carriers
            left_join(.,variant.id_carriers)
        }else if(TYPE=="AR"){
            extract_qualifed_gene.or.variant(.,type="1hit") ->> variant.id_carriers
            extract_qualifed_gene.or.variant(.,type="2hit") ->> Gene.refGene_carriers.2hit
            left_join(.,variant.id_carriers) %>%
                left_join(Gene.refGene_carriers.2hit)
        }
    } %>% left_join(gene_inh__ddg2p,by=c("Gene.refGene"="gene")) -> tmp

write_tsv(tmp,OUT,col_names=TRUE)



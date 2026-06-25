setwd("~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/同源基因")

suppressMessages({
  library(optparse)
  library(dplyr)
  
})


#-------------------------------------------------------------------------------
# 参数传递
#-------------------------------------------------------------------------------


option_list <- list(
  make_option("--SpeciesA", type="character", default=NULL, help="物种A"),
  make_option("--gffA", type="character", default=NULL, help="物种A的gff文件"),
  make_option("--SpeciesB", type="character", default=NULL, help="物种B"),
  make_option("--gffB", type="character", default=NULL, help="物种B的gff文件"),
  make_option("--orthofinderResult", type="character", default=NULL, help="orthofinder工具输出的结果")
)

args <- parse_args(OptionParser(option_list=option_list))

SpeciesA <- args$SpeciesA
gffA <- args$gffA
SpeciesB <- args$SpeciesB
gffB <- args$gffB
orthofinderResult <- args$orthofinderResult


gff.duck <- "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/MyAnalyses/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.gff"
gff.chicken <- "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/MyAnalyses/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff"


#-------------------------------------------------------------------------------
# 先处理gff文件，把GeneID、Symbol、ProteinID这些信息对应上
#-------------------------------------------------------------------------------







get_protein2gene <- function(gff, protein = TRUE){
  
  gff <- gff
  protein <- protein
  
  gff_lines <- readLines(gff)
  gff_lines <- gff_lines[!startsWith(gff_lines, "#")]
  
  gff_df <- read.delim(text = gff_lines, header = FALSE, sep = "\t", quote = "",
                       col.names = c("seqid","source","type","start","end", "score","strand","phase","attr"))
  
  if(protein){
    gff_df <- gff_df[grep(pattern = "protein_id", gff_df$attr ),] # 只保留attr中含有protein_id字符的行
  }
  
  p2g <- gff_df %>%
    mutate(
      GENEID = str_match(attr, "GeneID:([0-9]+)")[,2],              # 稳定数字ID
      SYMBOL = str_match(attr, "(?:^|;)gene=([^;]+)")[,2],           # 顺带留 symbol
      PROTEINID = str_match(attr, "(?:^|;)protein_id=([^;]+)")[,2]
    ) %>% 
    dplyr::select(GENEID, SYMBOL, PROTEINID ) %>%
    distinct() %>%
    filter(!is.na(SYMBOL))
  
  return(p2g)
}


p2g.chicken <- get_protein2gene(gff = gff.chicken)
p2g.duck <- get_protein2gene(gff = gff.duck)



#-------------------------------------------------------------------------------
# 读取 orthofinder 提取的结果，只提取单拷贝的1:1基因
#-------------------------------------------------------------------------------


GeneCount <- read.delim("~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/同源基因/Orthogroups/Orthogroups.GeneCount.tsv", header = T)
GeneCount <- GeneCount %>% filter(chicken == 1 & duck == 1)


Orthogroups <- read.delim("~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/同源基因/Orthogroups/Orthogroups.tsv", header = T) 
Orthogroups <- Orthogroups %>% filter(Orthogroup %in% GeneCount$Orthogroup)


Orthogroups <- Orthogroups %>% 
  left_join(y = p2g.duck, by = c("duck" = "PROTEINID")) %>% 
  dplyr::rename("GENEID.duck" = "GENEID", "SYMBOL.duck" = "SYMBOL")  %>% 
  left_join(y = p2g.chicken, by = c("chicken" = "PROTEINID")) %>% 
  dplyr::rename("GENEID.chicken" = "GENEID", "SYMBOL.chicken" = "SYMBOL")  %>%
  dplyr::rename("PROTEINID.duck" = "duck", "PROTEINID.chicken" = "chicken") 


outtable <- Orthogroups %>% dplyr::select(SYMBOL.chicken, SYMBOL.duck) %>% 
  dplyr::rename("Chicken" = "SYMBOL.chicken", "Duck" = "SYMBOL.duck") %>%
  distinct()


write.table(x = outtable, file = "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/同源基因/Chicken_vs_Duck.one_to_one.SingleCopy.orthologenes.txt", row.names = F, sep = "\t", quote = F)







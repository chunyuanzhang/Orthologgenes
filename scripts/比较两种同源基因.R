

#-------------------------------------------------------------------------------
# 从gff文件中提取基因列表
#-------------------------------------------------------------------------------

gff.duck <- "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/MyAnalyses/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.gff"
gff.chicken <- "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/MyAnalyses/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff"


p2g.chicken <- get_protein2gene(gff = gff.chicken)
p2g.duck <- get_protein2gene(gff = gff.duck)


#-------------------------------------------------------------------------------
# Ensembl下载的同源基因，以这些内容为主，我们自己鉴定的基因作为补充
#-------------------------------------------------------------------------------

Ensembl_Orthologene <- read.delim("~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/同源基因/Ensembl_Orthologene_GRCg7b_vs_CAU_duck1.0.tsv", header = T) %>% 
  dplyr::select(-Gene.stable.ID) %>% 
  dplyr::select(-Duck.gene.stable.ID) %>% 
  dplyr::select(Gene.name, Duck.gene.name, everything()) %>%
  filter(Gene.name != "") %>% filter(Duck.gene.name != "") %>%
  filter(Duck.homology.type == "ortholog_one2one") %>%
  filter(Duck.gene.name %in% p2g.duck$SYMBOL) %>%
  dplyr::rename("Chicken" = "Gene.name", "Duck" = "Duck.gene.name") %>%
  # mutate(combined = paste0(Chicken, "_vs_", Duck)) %>% 
  filter(Duck.orthology.confidence..0.low..1.high. > 0) %>%
  dplyr::select(Chicken, Duck)

# 由于ensemble使用的参考基因组版本不同，因此我们先检查一下是否ensembl获得的1:1同源基因都在T2T的基因列表中
# 【结果】有614个基因没有在T2T的基因组中，把这些都删除
# test <- Ensembl_Orthologene[!Ensembl_Orthologene$Duck.gene.name %in% p2g.duck$SYMBOL,]

# Ensemb共得到10157个1:1同源基因对


#-------------------------------------------------------------------------------
# orthofinder鉴定1:1同源基因
#-------------------------------------------------------------------------------


# 目前我们准备了两种同源基因
# 1、成对的1:1同源基因
# 2、单拷贝的1:1同源基因，但是由于单拷贝基因太少，不再参与分析

one2onegenes <- read.table("~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/同源基因/Chicken_vs_Duck.one_to_one.orthologenes.txt", header = T) %>% 
  # mutate(combined =paste0(Chicken, "_vs_", Duck)) %>%
  filter(Chicken %in% p2g.chicken$SYMBOL) %>%
  filter(Duck %in% p2g.duck$SYMBOL)

#-------------------------------------------------------------------------------
# Ensembl 和 orthofinder 鉴定的1:1基因互补一下
#-------------------------------------------------------------------------------

Orthologene <- full_join(Ensembl_Orthologene, one2onegenes, by = "Chicken") %>%
  rename("Duck.Ensembl" = "Duck.x", "Duck.T2T" = "Duck.y") %>%
  mutate(
    Duck = case_when(
      Duck.Ensembl == Duck.T2T ~ Duck.Ensembl,
      Duck.Ensembl != Duck.T2T & Chicken == Duck.Ensembl ~ Duck.Ensembl,
      Duck.Ensembl != Duck.T2T & Chicken == Duck.T2T ~ Duck.T2T,
      is.na(Duck.Ensembl) & !is.na(Duck.T2T) ~ Duck.T2T,
      !is.na(Duck.Ensembl) & is.na(Duck.T2T) ~ Duck.Ensembl,
      TRUE ~ Duck.T2T
      )
  ) %>%
  dplyr::select(Chicken, Duck)

#-------------------------------------------------------------------------------
# 输出同源基因
#-------------------------------------------------------------------------------

write.table(x = Orthologene, file = "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/同源基因/Chicken_vs_Duck.Ensembl_plus_Orthofinder.one_to_one.orthologenes.tsv", quote = F, sep = "\t", row.names = F)







suppressMessages(library(rtracklayer))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))

options(warn = -1)

message("使用参考基因组中的蛋白质序列作为输入文件，orthofinder工具简单蛋白质的同源蛋白,将蛋白ID映射到基因ID，并进一步提取1:1同源基因")

# SpeciesA <- "duck"
# gtfA <- "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/MyAnalyses/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.sorted.gtf"
# SpeciesB <- "chicken"
# gtfB <- "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/MyAnalyses/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf"
# orthofinderResult_1 <- "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/同源基因/Orthologues/Orthologues_chicken/chicken__v__duck.tsv"
# orthofinderResult_2 <- "~/Desktop/04.湘湖实验室/家禽病毒课题/01.数据分析/同源基因/Orthologues/Orthologues_duck/duck__v__chicken.tsv"


option_list <- list(
  make_option("--SpeciesA", type="character", default=NULL, help="物种A"),
  make_option("--gtfA", type="character", default=NULL, help="物种A的gtf文件"),
  make_option("--SpeciesB", type="character", default=NULL, help="物种B"),
  make_option("--gtfB", type="character", default=NULL, help="物种B的gtf文件"),
  make_option("--orthofinderResult", type="character", default=NULL, help="orthofinder工具输出的结果")
)

args <- parse_args(OptionParser(option_list=option_list))
SpeciesA <- args$SpeciesA
gtfA <- args$gtfA
SpeciesB <- args$SpeciesB
gtfB <- args$gtfB
orthofinderResult <- args$orthofinderResult


# =================================================
# 函数
# =================================================


get_gene_mapping <- function(gtf_file) {
  gtf_file <- gtf_file
  
  cat("处理", basename(gtf_file), "...\n")
  
  # 读取GTF
  gtf <- import(gtf_file)
  gtf_df <- as.data.frame(gtf)
  
  # 提取映射（优先使用protein_id，没有则用transcript_id）
  mapping <- gtf_df  %>% 
    dplyr::select(gene_id, protein_id) %>%
    filter(!is.na(protein_id)) %>%
    unique() 
  
  cat("提取", nrow(mapping), "个从蛋白ID到基因ID的映射\n")
  return(mapping)
}


protein_id_to_gene_id <- function(orthotable, species_column, species_mapping_dataframe ){
  orthotable <- orthotable
  species_column <- species_column
  species_mapping_dataframe <- species_mapping_dataframe
  
  orthotable[paste0(species_column, ".genes")] <- apply(orthotable, 1, function(row){
    proteins <- strsplit(x = row, split = ",") %>% unlist() %>% stringr::str_trim()
    genes <- species_mapping_dataframe[which(species_mapping_dataframe$protein_id %in% proteins),]$gene_id %>% unique()
    genes <- paste0(genes, collapse = ",")
    return(genes)
  })
  
  return(orthotable)
}


get_one_to_one_orthogene <- function(orthotable, species_columnA, species_columnB){
  orthotable <- orthotable
  species_columnA <- species_columnA
  species_columnB <- species_columnB
  
  one_to_one_genetable <- orthotable %>% 
    dplyr::select(species_columnA, species_columnB) %>%
    unique() %>%
    filter(!grepl(",", !!as.name(species_columnA))) %>%
    filter(!grepl(",", !!as.name(species_columnB)))

  dupinSpeciesA <- which(one_to_one_genetable[[species_columnA]] %>% table() > 1) %>% names()
  one_to_one_genetable %>% filter(!!as.name(species_columnA) %in% dupinSpeciesA) %>% filter(!!as.name(species_columnA) == !!as.name(species_columnB)) -> test1

  dupinSpeciesB <- which(one_to_one_genetable[[species_columnB]] %>% table() > 1) %>% names()
  one_to_one_genetable %>% filter(!!as.name(species_columnB) %in% dupinSpeciesB) %>% filter(!!as.name(species_columnA) == !!as.name(species_columnB)) -> test2

  one_to_one_genetable <- one_to_one_genetable %>% 
  filter(! .data[[species_columnA]] %in% dupinSpeciesA) %>% 
  filter(! .data[[species_columnB]] %in% dupinSpeciesB)

  one_to_one_genetable <- rbind(one_to_one_genetable, rbind(test1, test1) %>% unique())
  remove(test1, test2)

  cat("在 ", stringr::str_to_sentence(SpeciesA), " 和 ", stringr::str_to_sentence(SpeciesB), " 之间共提取到 ", dim(one_to_one_genetable)[1], " 个 1:1 同源基因\n" )
  return(one_to_one_genetable)
}


# =================================================
# 提取1:1同源基因
# =================================================

mapping.SpeciesA <- get_gene_mapping(gtf_file = gtfA)
mapping.SpeciesB <- get_gene_mapping(gtf_file = gtfB)
orthotable_1 <- read.delim(orthofinderResult_1) 
orthotable_2 <- read.delim(orthofinderResult_2) 



orthotable_1 <- protein_id_to_gene_id(orthotable = orthotable_1, species_column = SpeciesA, species_mapping_dataframe = mapping.SpeciesA )
orthotable_1 <- protein_id_to_gene_id(orthotable = orthotable_1, species_column = SpeciesB, species_mapping_dataframe = mapping.SpeciesB)

orthotable_2 <- protein_id_to_gene_id(orthotable = orthotable_2, species_column = SpeciesA, species_mapping_dataframe = mapping.SpeciesA )
orthotable_2 <- protein_id_to_gene_id(orthotable = orthotable_2, species_column = SpeciesB, species_mapping_dataframe = mapping.SpeciesB)


one_to_one_orthologgenes_1 <- get_one_to_one_orthogene(orthotable = orthotable_1, species_columnA = paste0(SpeciesA, ".genes"), species_columnB = paste0(SpeciesB, ".genes")) %>% unique()
one_to_one_orthologgenes_2 <- get_one_to_one_orthogene(orthotable = orthotable_2, species_columnA = paste0(SpeciesA, ".genes"), species_columnB = paste0(SpeciesB, ".genes")) %>% unique()

# 两个文件的内容完全一样

one_to_one_orthologgenes <- one_to_one_orthologgenes_1
colnames(one_to_one_orthologgenes) <- c( stringr::str_to_sentence(SpeciesA), stringr::str_to_sentence(SpeciesB) )

filename <- paste0(stringr::str_to_sentence(SpeciesA), "_vs_", stringr::str_to_sentence(SpeciesB), ".one_to_one.orthologenes.txt")
write.table(x = one_to_one_orthologgenes, file = filename,
            quote = F, append = F, sep = "\t", col.names = T, row.names = F)



####################
# 使用ensembl上的1:1同源基因验证当前的同源基因
####################


# library(biomaRt)
# 
# 
# # 连接Ensembl数据库
# ensembl_chicken <- biomaRt::useMart("ensembl", dataset="ggallus_gene_ensembl")
# ensembl_duck <- biomaRt::useMart("ensembl", dataset="aplatyrhynchos_gene_ensembl")
# 
# # 获取鸡的基因注释以及鸭对应的同源基因
# chicken_genes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 
#                                     'applatyrhynchos_homolog_ensembl_gene',
#                                     'applatyrhynchos_homolog_associated_gene_name',
#                                     'applatyrhynchos_homolog_orthology_type',
#                                     'applatyrhynchos_homolog_subtype'),
#                        mart=ensembl_chicken)
# 
# 
# # 筛选一对一同源基因
# ortholog_genes <- chicken_genes %>%
#   filter(applatyrhynchos_homolog_orthology_type == "ortholog_one2one") %>%
#   dplyr::select(external_gene_name, applatyrhynchos_homolog_associated_gene_name) %>%
#   unique()
# 
# ortholog_genes %>% 
#   filter(external_gene_name != applatyrhynchos_homolog_associated_gene_name) -> test
# 
# 
# 
# 
# ortholog_genes %>% 
#   filter(external_gene_name == applatyrhynchos_homolog_associated_gene_name) %>% 
#   filter( !is.na(external_gene_name)) %>%
#   filter(external_gene_name != "") -> test
#   
# 
# VennDiagram::venn.diagram(
#   list(
#     ensembl = test$external_gene_name,
#     orthofinder = one_to_one_orthologgenes$chicken.genes
#   ),
#   filename = "test.png"
# )











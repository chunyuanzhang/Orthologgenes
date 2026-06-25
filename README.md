
# 同源基因

## 使用示例【鉴定鸭和鸡1:1同源基因】

> 
> mkdir duck_vs_chicken    
> cd duck_vs_chicken
> 
> ln -s ~/zhangchunyuan/reference/IASCAAS_PekinDuck_T2T/GCF_047663525.1_IASCAAS_PekinDuck_T2T_protein.faa ./duck.faa  
> ln -s ~/zhangchunyuan/reference/bGalGal1_mat_broiler_GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_protein.faa ./chicken.faa  
> 
> cd ../    
> ~/tools/OrthoFinder/orthofinder.py -f duck_vs_chicken   


## 输出结果

### 一、元信息与运行记录
Citation.txt — 该次分析用到的方法及对应引用(OrthoFinder 本身,以及它调用的比对器 DIABLO/MAFFT、建树工具 FastTree/IQ-TREE 等)。写论文方法部分直接照抄即可,能避免漏引底层工具。   
Log.txt — 完整运行日志:输入了几个物种、各步骤参数、用了哪个序列搜索程序(DIAMOND/BLAST)、inflation 参数等。排查结果异常或复现分析时看这个。   
WorkingDirectory — 中间文件(BLAST/DIAMOND 比对结果、序列 ID 映射、图聚类中间产物)。一般不用动,但有两个用途:一是用 -b 参数复用已算好的比对结果重跑(加物种、改参数时省去最耗时的全对全搜索);二是里面的 SpeciesIDs.txt/SequenceIDs.txt 是内部数字 ID 与原始基因名的对照表,调试时有用。   


### 二、核心同源分组(最常用)
Orthogroups — 文件夹,核心产物。里面通常有:   

Orthogroups.tsv：物种为列、OG 为行的矩阵,每格是该物种在该 OG 里的基因(逗号分隔)。这是下游分析的主表。   
Orthogroups.GeneCount.tsv：同样布局,但格子里是基因数而非基因名——做扩张/收缩分析、CAFE 输入、PCA/聚类时直接用。   
Orthogroups_SingleCopyOrthologues.txt：全物种严格单拷贝的 OG 列表。   
Orthogroups_UnassignedGenes.tsv：没能归入任何 OG 的"孤儿基因",常是物种特异基因或注释噪音。   

Orthogroup_Sequences — 每个 OG 一个 FASTA(含所有物种该 OG 的全部成员)。可拿去做功能注释(eggNOG/InterProScan)、单个家族的系统发育、motif 分析。   
Single_Copy_Orthologue_Sequences — 每个严格单拷贝 OG 一个 FASTA。建物种树/超矩阵串联的标准输入,直接 MAFFT+IQ-TREE 即可。   
Phylogenetic_Hierarchical_Orthogroups (HOGs) — 这是 OrthoFinder2 之后官方推荐优先使用的结果,常在 N0.tsv。它基于物种树在每个分支节点重新定义 orthogroup,比传统 Orthogroups.tsv 更准——尤其能纠正因古老重复导致的错误合并/拆分。N0 是最深的根节点 HOG,最接近"真正的基因家族";N1、N2… 对应更近的节点。做基因家族演化、拷贝数比较时,建议用 N0.tsv 替代 Orthogroups.tsv。   

### 三、统计概览
Comparative_Genomics_Statistics — 一堆汇总表:每物种被分配进 OG 的基因比例、OG 大小分布、物种特异 OG 数、单拷贝 OG 数、Statistics_Overall.tsv 总览等。写结果段落、做质控(某物种分配率异常低往往意味着注释质量或物种选择有问题)时第一时间看这里。   

### 四、比对与基因树
MultipleSequenceAlignments — 每个 OG 的多序列比对(以及串联的物种树比对)。可直接复用做选择压力分析、保守位点提取,省去自己重比对。   
Gene_Trees — 每个 OG 的原始基因树。   
Resolved_Gene_Trees — 用物种树对基因树做了根定位和重排后的版本(区分了重复节点与物种分化节点)。做基因树物种树调和(reconciliation)、判断 ortholog/paralog、GeneRax/Notung 之类分析时用这个,而不是原始 Gene_Trees。   

### 五、物种树
Species_Tree — SpeciesTree_rooted.txt(STAG/STRIDE 推断并定根)及带支持值版本。可直接作为系统发育框架,也可作为 CAFE、祖先状态重建、比较方法(PIC/PGLS)的输入树。注意它是基因组水平的快速估计,发表级物种树通常还会再用串联/合并法精修。   

### 六、直系同源关系与进化事件
Orthologues — 两两物种间的直系同源对(Species1__v__Species2.tsv),含 1:1 / 1:many / many:many。做共线性、dN/dS 配对、跨物种功能映射时用。   
Gene_Duplication_Events — STRIDE 推断的基因重复事件,标注每个重复发生在物种树的哪个分支、涉及哪些基因。研究全基因组/小规模重复、基因家族扩张的时间点时核心用这个;Duplications.tsv 还给了每个重复节点的支持度。   
Putative_Xenologs — 疑似**水平基因转移(HGT)**来源的基因对——即基因树拓扑与物种树严重冲突、更像外源获得而非垂直遗传的同源基因。研究 HGT 的起点,但需谨慎验证(也可能是建树误差或不完全谱系分选)。   
Phylogenetically_Misplaced_Genes — 在基因树中位置异常、与所在 OG 其他成员系统发育位置矛盾的基因。常是注释错误、组装污染或嵌合基因的信号,做数据质控时值得排查。   

### 下游分析路线参考
几条常见的串联思路:   
物种系统发育 → Single_Copy_Orthologue_Sequences(或从 HOG 取)→ 串联/合并 → IQ-TREE/ASTRAL。   
基因家族演化 → N0.tsv(HOG)+ Species_Tree → CAFE 做扩张/收缩 → 富集分析解读。   
选择压力 → Orthologues 取 1:1 对 + MultipleSequenceAlignments → PAML/HyPhy 算 dN/dS。   
功能注释 → Orthogroup_Sequences → eggNOG/InterProScan → 按 OG 汇总功能。   

## 1:1同源基因提取

### 严格单拷贝同源基因提取



### 基于基因树推断的直系同源关系提取





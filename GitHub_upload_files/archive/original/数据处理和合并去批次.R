
# 下载counts表达矩阵####
#加载 R 包 TCGAbiolinks
library(TCGAbiolinks)

#查看各个癌种的项目 id
projects = getGDCprojects()$project_id
projects

#选择研究的癌症项目，下载 Counts 数据
cancer <- "TCGA-LUAD"
query_exp = GDCquery(project = cancer,#下载项目的名称，以下载胃癌为例
                     data.category = "Transcriptome Profiling",#选择转录组数据
                     data.type = "Gene Expression Quantification",#选择基因表达量，这个选项和 GDC 官网上看到的一样
                     workflow.type = "STAR - Counts")#数据类型选择 counts 格式

#使用 GDCdownload() 进行下载，这里需要一些时间运行
GDCdownload(query = query_exp)
#GDCdownload(query = query_exp, method = "api", directory = "GDCdata", files.per.chunk = 5)
#下载数据，默认存放位置为当前工作目录下的 GDCdata 文件夹中
#query:表示要下载的数据
#method:表示下载的方法，有2种，“api"和"client"。
#       “api"更快一些，但是下载较大数据时容易中断掉，可能需要重新下载，可以设置files.per.chunk来指定单次下载文件的数量。
#directory:表示下载存储的目录

#用 GDCprepare 批量读取数据。文件下载后不可移动，否则无法运行该函数
pre_exp = GDCprepare(query = query_exp, save = T, save.filename = paste0(cancer,"_RNAseq.Rdata"))
#此时如果想把读取的数据保存为 Rdata，下次直接加载，就把 save 参数设置为TRUE
#save.filename 给一个名称。

#这时候，这个 pre_exp 是个对象，需要从里面提取表达量的信息，
# pre_exp 含有coldata, rowdata，meta-data，以及最重要的assay，共有6个assay。
library(SummarizedExperiment)
pre_exp
#class: RangedSummarizedExperiment #类型
#dim: 60660 407 #维度：60660行，177列
#metadata(1): data_release
#assays(6): unstranded stranded_first ... fpkm_unstrand fpkm_uq_unstrand #6个assay以及它们的名字
#rownames(60660): ENSG00000000003.15 ENSG00000000005.6 ... ENSG00000288674.1 ENSG00000288675.1 #表达矩阵的行名
#rowData names(10): source type ... hgnc_id havana_gene #行信息（也就是基因信息），比如gene id,gene name,gene type
#colnames(407): TCGA-HU-A4GN-11A-12R-A251-31 TCGA-BR-8484-01A-11R-2402-13 ... TCGA-BR-8369-01A-11R-2343-13 #表达矩阵的列名（也就是样本名）
#TCGA-FP-A8CX-01A-11R-A36D-31 
#colData names(109): barcode patient ... paper_CIMP.Category paper_stringAsFactor #列信息，也就是样本信息，比如生存时间、生存状态这些

# 每个assay你可以理解为一个表达矩阵，我们需要的counts矩阵、TPM矩阵、FPKM矩阵就是其中一个~
# 查看每个assay的名字
names(assays(pre_exp))
#"unstranded"       "stranded_first"   "stranded_second"  "tpm_unstrand"     "fpkm_unstrand"    "fpkm_uq_unstrand"
#"unstranded":counts格式
#"tpm_unstrand"：tpm格式
#"fpkm_unstrand"：fpkm格式

#rowData你可以理解为一个包含基因各种信息的数据框。
#其中 gene_name 就是 gene_symbol，我们的基因 id 转换就用这一列信息
# 提取rowData
rowdata <- rowData(pre_exp)
#看看rowData包括哪些内容，可以看到里面有我们需要的gene_name
names(rowdata)
#[1] "source"      "type"        "score"       "phase"       "gene_id"     "gene_type"   "gene_name"   "level"       "hgnc_id"     "havana_gene"
head(rowdata$gene_name)
#[1] "TSPAN6"   "TNMD"     "DPM1"     "SCYL3"    "C1orf112" "FGR"
#查看基因名的数量
length(rowdata$gene_name)
#[1] 60660

#有了这些东西，就可以提取表达矩阵了，直接使用assay()搞定！
# mRNA的counts矩阵
expr_counts_mrna <- assay(pre_exp,"unstranded")

#展示一下：
expr_counts_mrna[1:10,1:2]

#添加gene_symbol也就非常简单了，只要提取gene_name这一列，然后和原来的表达矩阵合并即可！
#先提取 gene_name
symbol_mrna <- rowData(pre_exp)$gene_name
head(symbol_mrna)
#[1] "TSPAN6"   "TNMD"     "DPM1"     "SCYL3"    "C1orf112" "FGR"
#和基因表达矩阵合并就行了：
symbol_mrna_symbol = cbind(data.frame(symbol_mrna),
                           as.data.frame(expr_counts_mrna))

#但是，此时gene_symbol是有重复的，需要去重
#查看gene_symbol中是否有重复
table(duplicated(symbol_mrna))
#FALSE  TRUE 
#59427  1233
#1233个基因名称是重复

#去重复也很简单，这里我们合并所有重复的 gene symbol
#主要思路为利用aggregate函数根据symbol列中的相同基因合并基因表达矩阵矩阵
#使用aggregate根据symbol列中的相同基因进行合并
counts = aggregate(.~symbol_mrna, sum, data = symbol_mrna_symbol)
counts[1:5, 1:4]
table(duplicated(counts))
#FALSE 
#59427
#基因名称无重复

#保存该基因表达矩阵文件
# 下载tpm表达矩阵####
library(TCGAbiolinks)

#查看各个癌种的项目 id
projects = getGDCprojects()$project_id
projects

#选择研究的癌症项目，下载STAR - Counts数据
cancer <- "TCGA-LUAD"
query_exp = GDCquery(project = cancer,#下载项目的名称，以下载胃癌为例
                     data.category = "Transcriptome Profiling",#选择转录组数据
                     data.type = "Gene Expression Quantification",#选择基因表达量，这个选项和 GDC 官网上看到的一样
                     workflow.type = "STAR - Counts")#数据类型选择 counts 格式

#使用 GDCdownload() 进行下载，这里需要一些时间运行
GDCdownload(query = query_exp)

#用 GDCprepare 批量读取数据。文件下载后不可移动，否则无法运行该函数
pre_exp = GDCprepare(query = query_exp, save = T, save.filename = paste0(cancer,"_RNAseq.Rdata"))
library(SummarizedExperiment)
pre_exp
names(assays(pre_exp))

# 提取rowData(包含基因各种信息的数据框)
rowdata <- rowData(pre_exp)
names(rowdata)
head(rowdata$gene_name)
length(rowdata$gene_name)
#[1] 60660

#提取tpm表达矩阵
expr_tpm_mrna <- assay(pre_exp,"tpm_unstrand")
expr_tpm_mrna[1:10,1:2]

#提取 gene_name
symbol_mrna <- rowData(pre_exp)$gene_name
head(symbol_mrna)
#和基因表达矩阵合并
symbol_mrna_symbol = cbind(data.frame(symbol_mrna),
                           as.data.frame(expr_tpm_mrna))

#查看gene_symbol中是否有重复
table(duplicated(symbol_mrna))
# FALSE  TRUE 
# 59427  1233

#使用aggregate根据symbol列中的相同基因进行合并
tpm = aggregate(.~symbol_mrna, sum, data = symbol_mrna_symbol)
tpm[1:5, 1:4]
table(duplicated(tpm))
#FALSE 
#59427

#保存该基因表达矩阵文件
# 下载CNV数据####
library(TCGAbiolinks)
cancer <- "TCGA-LUAD"
query_exp = GDCquery(project = cancer,
                     data.category = "Copy Number Variation",
                     data.type="Gene Level Copy Number")

GDCdownload(query = query_exp, method = "api", directory = paste0(cancer,"_CNV"), files.per.chunk = 5)
# 下载maf文件####
library(TCGAbiolinks)
cancer <- "TCGA-LUAD"
query_exp = GDCquery(project = cancer,
                     data.category = "Simple Nucleotide Variation",
                     data.type="Masked Somatic Mutation")

GDCdownload(query = query_exp, method = "api", directory = paste0(cancer,"_SNV"), files.per.chunk = 5)
# TCGA临床信息下载与整理####
cancer <- "TCGA-LUAD"
library(stringr)
library(do)

####读取数据####
# tcga_cli <- read.table(paste("/work/database/TCGA/",cancer,"-clinical.csv",sep=""), header=T, sep=",", row.names = 1)
# str(tcga_cli)   #查看数据结构

#如果没有已下载的就用下面的代码下载临床数据
tcga_cli <- GDCquery_clinic(project = cancer, type = "clinical")#记得修改病种
write.csv(tcga_cli ,file = paste0(cancer,'_clinical.csv'))

####提取几列临床信息####
tcga_cli1 <- tcga_cli[,c('submitter_id',
                         'age_at_index',
                         'gender',
                         'ajcc_pathologic_stage',
                         'ajcc_pathologic_n',
                         'ajcc_pathologic_m',
                         'ajcc_pathologic_t',
                         'race',
                         'vital_status',
                         'days_to_last_follow_up',
                         'days_to_death')]
colnames(tcga_cli1) <- c('id',
                         'age',
                         'gender',
                         'stage',
                         'N',
                         'M',
                         'T',
                         'race',
                         'vital_status',
                         'follow_up',
                         'days_to_death')

####临床信息整理####
#生存状态与生存时间转换
tcga_cli1 <- tcga_cli1[!is.na(tcga_cli1$vital_status),]
status<-tcga_cli1$vital_status=='Alive'
tcga_cli1$fustat<-ifelse(status,0,1) #0代表生存，1代表死亡
tcga_cli1$futime<-ifelse(status,tcga_cli1$follow_up,tcga_cli1$days_to_death)
tcga_cli1 <- tcga_cli1[,-which(colnames(tcga_cli1) %in% c('vital_status','follow_up','days_to_death'))]
tcga_cli1$fustat <- as.numeric(tcga_cli1$fustat)
tcga_cli1$futime <- as.numeric(tcga_cli1$futime)
#NA和not reported改为unknown
tcga_cli2 <- tcga_cli1
tcga_cli2[is.na(tcga_cli2)] <- "unknown"
tcga_cli2 <- Replace(data = tcga_cli2,from = "not reported",to = "unknown")
#将stage的ABCD结尾去掉
tcga_cli2$stage <- gsub("[A|B|C|D]","",tcga_cli2$stage)
#以X结尾改为unknown，将abcd结尾去掉
tcga_cli2$T <- gsub("a$|b$|c$|d$","",tcga_cli2$T)
tcga_cli2$T[grep("X$",tcga_cli2$T,value = F)] <- "unknown"
tcga_cli2$M <- gsub("a$|b$|c$|d$","",tcga_cli2$M)
tcga_cli2$M[grep("X$",tcga_cli2$M,value = F)] <- "unknown"
tcga_cli2$N <- gsub("a$|b$|c$|d$","",tcga_cli2$N)
tcga_cli2$N[grep("X$",tcga_cli2$N,value = F)] <- "unknown"

####统计临床类别频数####
sapply(tcga_cli2, class) #查看每列的数据类型
#转换特定列为因子
tcga_cli2[,-which(colnames(tcga_cli2) %in% c('id','age','fustat','futime'))] <- as.data.frame(unclass(tcga_cli2[,-which(colnames(tcga_cli2) %in% c('id','age','fustat','futime'))]), stringsAsFactors = TRUE) 
#统计频数
tcga_cli2 %>% dplyr::group_by() %>% summary()
# 
# 
# TNM的临床分期整理规范
# 实际上在临床工作中，每一个部位的肿瘤常常会使用2个TNM分期，一个是临床分期（cTNM或THM），一个是病理分期（pTNM）。
# T——原发肿瘤的范围
# Tx：无法评估原发肿瘤
# T0：未发现原发肿瘤的证据
# Tis：原位癌
# T1、T2、T3、T4：原发肿瘤的大小和（或）局部侵犯的范围一次递增
# N——有无区域淋巴结转移以及范围
# Nx：无法评估区域淋巴结的情况
# N0：未出现区域淋巴结转移
# N1、N2、N3：区域淋巴结的侵犯范围一次递增
# 注：除区域淋巴结以外的任何淋巴结转移都被划分为远处转移
# M——有无远处转移
# Mx：无法评估远处转移的情况
# M0：未出现远处转移
# M1：发生远处转移
# 
# 数据整理原则：
# 先不删掉unknown，空白等改为unknown，保存一份完整的临床信息，在具体分析时，去掉要用到临床信息中包含unknown的样本
# Age：分为<=65和>65两组
# Grade：GX改为unknown
# Stage：Stage ⅠA和Stage ⅠB均改为Stage Ⅰ，以此类推
# T,M,N：T1A,M1A,N1A改为T1,M1,N1，以此类推；TX,MX,NX改为unknown
# 如果全部改成二分类，那么:
#   Age分为<=65和>65两组；
# Stage分为StageI&StageII和StageIII&StageIV两组；
# T分为T1&T2和T3&T4两组；
# M分为M0和M1两组；
# N分为N0和N1&N2&N3两组
# 注意：
# stage,T,M,N是比较重要的信息，如果有这些信息，就要用！

# 提取有临床信息的TCGA肿瘤样本（log2(tpm+1)）
#TPM矩阵与临床信息读取
cancer <- "TCGA-LUAD"
tcga<-read.csv(paste0(cancer,"_tpm.csv"),row.names = 1,check.names = F)
cli<-read.csv(paste0(cancer,"_clinical_treated.csv"),row.names = 1)

#提取肿瘤样本
library(stringr)
table(str_sub(colnames(tcga),14,14)=='0') #统计肿瘤和癌旁样本的数目
tcga_tumor<-tcga[,str_sub(colnames(tcga),14,14)=='0']

#肿瘤样本名取前12位，并去重（临床信息样本名是12位字符）
colnames(tcga_tumor)<-gsub('[.]','-',colnames(tcga_tumor))
colnames(tcga_tumor)<-str_sub(colnames(tcga_tumor),1,12)
tcga_tumor<-tcga_tumor[,!duplicated(colnames(tcga_tumor))]

#提取有临床信息的肿瘤样本
tcga_sample<-intersect(colnames(tcga_tumor),rownames(cli))
tcga_tumor<-tcga_tumor[,tcga_sample]
cli<-cli[tcga_sample,]
#去除没有完整生存信息的样本
the_sample <- rownames(cli)[c(which(cli$futime %in% "unknown"),
                              which(cli$fustat %in% "unknown"),
                              which(is.na(cli$futime)),
                              which(is.na(cli$fustat)))]
tcga_tumor<-tcga_tumor[,-which(colnames(tcga_tumor) %in% the_sample)]
cli<-cli[-which(rownames(cli) %in% the_sample),]

#取对数
tcga_tumor_log<-log2(tcga_tumor+1)
tcga<-log2(tcga+1)


#2. GSE11969 data process####



# GEO微阵列数据集下载####
GEOquery这个包获取的表达矩阵为标准化的表达矩阵
setwd("./data/")
library(GEOquery)
library(dplyr)
#Setting options('download.file.method.GEOquery'='auto') or options( 'download.file.method.GEOquery' = 'libcurl' )
#Setting options('GEOquery.inmemory.gpl'=FALSE)

#下载或加载数据集
f='GSE11969_eSet.Rdata'
if(!file.exists(f)){
  gset <- getGEO('GSE11969', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE11969_eSet.Rdata')  ## 载入数据

#样本信息获取
ph <- pData(gset[[1]])
write.csv(ph,'GSE11969_ph.csv')
#样本信息整理
ph1 <- data.frame(sample=ph$geo_accession,  #GSM编号
                  group=ph$source_name_ch1,  #分组信息
                  age=ph$`age.at.diagnosis (year):ch1`, #年龄
                  gender=ph$`Sex:ch1`,                  #性别
                  stage=ph$`tnm.stage:ch1`,             #病理分期
                  N=ph$`tnm.n:ch1`,                     #N
                  M=ph$`tnm.m:ch1`,                     #M
                  T=ph$`tnm.t:ch1`,                     #T
                  fustat=ph$`os.event:ch1`,             #生存状态
                  futime=ph$`os.delay (months):ch1`    #生存时间
)
rownames(ph1) <- ph1$sample
write.csv(ph1,'GSE11969_ph1.csv')

#表达矩阵获取
exp <- exprs(gset[[1]])

#转换gene ID
GPL <- getGEO('GPL7015')
GPL_id <- GPL@dataTable@table[,c('ID','Gene Symbol')]
GPL_id$`Gene Symbol` <- stringr::str_split(GPL_id$`Gene Symbol`,'///',simplify = T)[,1]
rownames(GPL_id) <- GPL_id$ID
exp1 <- cbind(data.frame(gene_name=GPL_id[rownames(exp),2]),exp)
exp2 <- na.omit(exp1)
exp2 <- exp2[which(exp2$gene_name != ""),]
exp2 <- aggregate(.~gene_name, mean, data = exp2)
exp2 <- exp2 %>% tibble::column_to_rownames("gene_name")
write.csv(exp2,'GSE11969.csv')



#3. 合并去批次####

# 多个表达量数据集的整合分析####

# 导入数据集1
data1 <- exp2 
# 导入数据集2
data2 <- tcga

# 合并数据集
library(sva) 
library(tidyverse) 
merge_eset <- inner_join(data1, data2, by = "X", suffix = c(".data1", ".data2")) #X是第一列名，需要将基因名放在第一列而不是成为行名
merge_eset <- merge_eset %>% tibble::column_to_rownames("X")
dim(merge_eset)    #查看数据维度
exp <- as.matrix(merge_eset)
dimnames <- list(rownames(exp),colnames(exp))
exp <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
dim(exp) #查看数据维度
class(exp) #返回"matrix" "array"
write.csv(exp,"exp_nocombat.csv")

# 过滤掉低表达的基因
Expr <- exp[rowSums(exp)>1,]
dim(Expr) #查看数据维度
# 取log
Expr=log2(Expr+1)
range(Expr)

# 获取分组\批次信息
group_list <- data.frame(
  sample = colnames(Expr), 
  dataset = c(rep("GSE11969",95),
              rep("TCGA-LUAD",510))) #dataset表示不同的批次信息
rownames(group_list) <- group_list$sample
group_list$dataset <- as.factor(group_list$dataset)
write.csv(group_list,"groupinfo.csv")

# PCA分析可视化各数据集在去批次前的聚簇效果
library(FactoMineR)
library(factoextra)
pdf(file = "PCA_before.pdf",width = 7,height = 6)
pre.pca <- PCA(t(Expr),graph = FALSE)
fviz_pca_ind(pre.pca,
             geom= "point",
             col.ind = group_list$dataset,
             addEllipses = TRUE,
             legend.title="datasets"
)
dev.off()

## 提取样本的坐标  
## 注意：这里的$ind$coord是一个矩阵，行表示样本，列表示主成分  
# sample_coordinates <- pre.pca$ind$coord  
# head(sample_coordinates)  
## 只看前两个主成分（PC1和PC2）上的坐标  
# sample_coordinates_PC1_PC2 <- sample_coordinates[,1:2]  %>% as.data.frame()
## 查看样本的坐标（只显示PC1和PC2）  
# head(sample_coordinates_PC1_PC2)  
## 提取离群样本名，如这里为Dim.2 < -50的样本
# Outlier_samples <- rownames(sample_coordinates_PC1_PC2)[which(sample_coordinates_PC1_PC2$Dim.2 < -50)]
## 去除离群样本
# Expr <- Expr[,-which(colnames(Expr) %in% Outlier_samples)]
# group_list <- group_list[-which(rownames(group_list) %in% Outlier_samples),]
# write.csv(group_list,"groupinfo2.csv")

#-------------------进行去批次------------------#
# 使用sva包计算批次效应
library(sva)
Expr_combat <- ComBat(Expr, batch = group_list$dataset) 
write.csv(Expr_combat,"exp_combat.csv")

# 查看去除批次效应的结果如何
library(tinyarray)
draw_pca(exp = Expr_combat, group_list = group_list$dataset)

#做校正的PCA分析
pdf(file = "PCA_after.pdf",width = 7,height = 6)
pre.pca <- PCA(t(Expr_combat),graph = FALSE)
fviz_pca_ind(pre.pca,
             geom= "point",
             col.ind = group_list$dataset,
             addEllipses = TRUE,
             legend.title="datasets")
dev.off()


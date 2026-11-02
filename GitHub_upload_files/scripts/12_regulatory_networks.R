# Module 12: regulatory_networks
# Source snapshot: archive/original/12.互作网络.R
# Path adaptation only; scientific settings and historical seeds are retained.
# Run using: Rscript scripts/run_module.R 12
if (!exists("project_dir", mode = "function")) {
  stop("Run this module through scripts/run_module.R from the repository root.")
}

genes <- read.csv("../4.machine_learning_1/ML_Gene.csv")$x
miRNA_mRNA <- read.csv(resource_path("resource_d0c0f93a27"),row.names = 1)
# miRNA_mRNA2 <- read.csv('/work/public/ceRNA/starbase_miRNA_mRNA.csv',row.names = 1)
# miRNA_mRNA <- dplyr::intersect(miRNA_mRNA1,miRNA_mRNA2)
net1<-miRNA_mRNA[miRNA_mRNA$mRNA %in% genes,]
head(net1)
length(unique(net1$miRNA))#查看交集到的miRNA数目
length(unique(net1$mRNA))#查看交集到的mRNA数目

#读取starbase数据库种的miRNA_lncRNA关系对
#根据可靠性分为了Low、Median、High三个严格程度，可自行调整
miRNA_lncRNA<-read.csv(resource_path("resource_4d1dbc6d3e"),row.names = 1)

net2<-miRNA_lncRNA[miRNA_lncRNA$miRNA%in%net1$miRNA,]
head(net2)
length(unique(net2$miRNA))#查看交集到的miRNA数目
length(unique(net2$lncRNA))#查看交集到的lncRNA数目

#lncRNA太多时根据miRNA数目进行筛选
lncRNA<-as.data.frame(table(net2$lncRNA))
lncRNA1<-lncRNA[lncRNA$Freq>=5,]#阈值自行调整

#提取整合
net2<-net2[net2$lncRNA%in%lncRNA1$Var1,]
colnames(net1)<-c('miRNA','Target')
colnames(net2)<-c('miRNA','Target')
ceRNA<-rbind(net1,net2)
ceRNA$type<-c(rep('mRNA',nrow(net1)),rep('lncRNA',nrow(net2)))

#保存
write.csv(ceRNA,'../12.network/ceRNA.csv',row.names = F,quote = F)


dir.create("../12.network/mRNA_RBP")
# genes=read.table('prognosis/prognosis_coefficients.csv',header = T)[,1]  #更改为感兴趣的基因名文件

#之前网页下载版本，使用hg19
#for(gene in genes){
#  file=paste0("mRNA_RBP/",gene,".txt")
#  link=paste0('https://starbase.sysu.edu.cn/api/RBPTarget/?assembly=hg19&geneType=mRNA&RBP=all&clipExpNum=5&pancancerNum=0&target=',gene,'&cellType=all') 
#  download.file(link,file)
#  Sys.sleep(5)
#}

for(gene in genes){
  file=paste0("../12.network/mRNA_RBP/",gene,".txt")
  link=paste0('https://rnasysu.com/encori/api/RBPTarget/?assembly=hg38&geneType=mRNA&RBP=all&clipExpNum=5&pancancerNum=0&target=',gene,'&cellType=all') 
  download.file(link,file)
  Sys.sleep(5)
}


#批量读取,合并表格放到cytoscape中
#首先就是获得存储路径下所有的txt格式文件的文件名
#full.names参数如果设置为TRUE，则返回文件的完整路径，如果设置的为FALSE则只返回文件名。
files<-dir(path = "../12.network/mRNA_RBP/",
           full.names = T,
           pattern = ".txt")
library(tidyverse)
FUN2 = function(y){read.table(y,header = T,sep = '\t')}
df_RBP<-map(files,FUN2)
class(df_RBP)
RBP_table=data.frame()
for(i in 1:length(df_RBP)){
  RBP_table=rbind(RBP_table,data.frame(RBP=df_RBP[[i]]$RBP,geneName=df_RBP[[i]]$geneName,clusterNum=df_RBP[[i]]$clusterNum))
}



library(dplyr)
RBP_table=RBP_table[RBP_table$clusterNum>=5,]
RBP_table1=dplyr::distinct(RBP_table,.keep_all=T)
write.csv(RBP_table1,'../12.network/mRNA_RBP/RBP_table.csv',row.names = F,quote = F)
node <-data.frame(c(unique(RBP_table1$RBP),unique(RBP_table1$geneName)),
                  c(rep("RBP",length(unique(RBP_table1$RBP))),
                    rep("geneName",length(unique(RBP_table1$geneName)))))
colnames(node) <- c("node","type")
write.csv(node,"../12.network/mRNA_RBP/node_type.csv",quote = F,row.names = F)


# TF ----------------------------------------------------------------------
dir.create('../12.network/TF')
TF <- read.table(resource_path("resource_4ae41aa186"), header = T, sep = "\t")
TF1 <- TF[TF$target %in% genes, ]
TF2 <- TF1[grep('lung',TF1$tissue,ignore.case = T),]
write.csv(TF2,'../12.network/TF/TF-mRNA.csv',row.names = F,quote = F)
node <-data.frame(c(unique(TF2$TF),unique(TF2$target)),
                  c(rep("TF",length(unique(TF2$TF))),
                    rep("geneName",length(unique(TF2$target)))))
colnames(node) <- c("node","type")
write.csv(node,'../12.network/TF/TF-mRNA_target.csv',row.names = F,quote = F)

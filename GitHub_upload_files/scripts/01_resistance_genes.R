# Module 01: resistance_genes
# Source snapshot: archive/original/1.Coregene.R
# Path adaptation only; scientific settings and historical seeds are retained.
# Run using: Rscript scripts/run_module.R 01
if (!exists("project_dir", mode = "function")) {
  stop("Run this module through scripts/run_module.R from the repository root.")
}

setwd(project_dir())
dir.create('1.Coregene_1')
setwd('1.Coregene_1/')

####R包和函数加载####
library(limma)
library(edgeR)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(ggpubr)
library(tidyr)
library(tibble)
source(helper_file("pheatmap1.R"), encoding = "UTF-8")
source(helper_file("plot_boxplot.R"), encoding = "UTF-8")
source(helper_file("plot_vioplot.R"), encoding = "UTF-8")

####PC9####
exp1<-read.csv("../0.data/PC9.csv",check.names = F)#设置了check.names=F就不会把列名里的-换成.了
exp1 <- exp1[!duplicated(exp1$Sample),]
rownames(exp1) <- exp1$Sample
exp1[,1] <- NULL
exp1 <- na.omit(exp1)
exp1 <- log2(exp1+1)
range(exp1)
#0.00000 14.66911
group_list<-read.csv("../0.data/PC9_group.csv")
colnames(exp1) <- gsub('[.]','-',colnames(exp1))
table(colnames(exp1)==group_list$sample)
exp1=exp1[,group_list$sample]

group1=group_list$group
names(group1)=colnames(exp1)
table(group1)
# group1
# OR Parental 
# 3        3
colors <- c("#1E90FF","#F08080","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
####参数####
exp=exp1
group=group1
case="Parental"
control="OR"
P="p.adj"
P_value=0.05
logFC_value=0.5
colors=colors
is_count=F
feature_type="基因"
lab = c("gene", "Expression Level")
flip = F
####运行####
source(helper_file("差异分析.R"), encoding = "UTF-8")
# down    no    up 
# 1952 32079  1973



####H1975####
exp1<-read.csv("../0.data/H1975.csv",check.names = F)#设置了check.names=F就不会把列名里的-换成.了
exp1 <- exp1[!duplicated(exp1$Sample),]
rownames(exp1) <- exp1$Sample
exp1[,1] <- NULL
exp1 <- na.omit(exp1)
exp1 <- log2(exp1+1)
range(exp1)
# 0.00000 14.20921
group_list<-read.csv("../0.data/H1975_group.csv")
colnames(exp1) <- gsub('[.]','-',colnames(exp1))
table(colnames(exp1)==group_list$sample)
exp1=exp1[,group_list$sample]

group1=group_list$group
names(group1)=colnames(exp1)
table(group1)
# group1
# OR Parental 
# 3        3
colors <- c("#1E90FF","#F08080","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
####参数####
exp=exp1
group=group1
case="Parental"
control="OR"
P="p.adj"
P_value=0.05
logFC_value=0.5
colors=colors
is_count=F
feature_type="基因"
lab = c("gene", "Expression Level")
flip = F
####运行####
source(helper_file("差异分析.R"), encoding = "UTF-8")
# down    no    up 
# 1348 33875   781


####HCC827####
exp1<-read.csv("../0.data/HCC827.csv",check.names = F)#设置了check.names=F就不会把列名里的-换成.了
exp1 <- exp1[!duplicated(exp1$Sample),]
rownames(exp1) <- exp1$Sample
exp1[,1] <- NULL
exp1 <- na.omit(exp1)
exp1 <- log2(exp1+1)
range(exp1)
#0.0000 15.0081
group_list<-read.csv("../0.data/HCC827_group.csv")
colnames(exp1) <- gsub('[.]','-',colnames(exp1))
table(colnames(exp1)==group_list$sample)
exp1=exp1[,group_list$sample]

group1=group_list$group
names(group1)=colnames(exp1)
table(group1)
# group1
# OR Parental 
# 3        3
colors <- c("#1E90FF","#F08080","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
####参数####
exp=exp1
group=group1
case="Parental"
control="OR"
P="p.adj"
P_value=0.05
logFC_value=0.5
colors=colors
is_count=F
feature_type="基因"
lab = c("gene", "Expression Level")
flip = F
####运行####
source(helper_file("差异分析.R"), encoding = "UTF-8")
# down    no    up 
# 2534 31633  1837


###交集###
dir.create('Venn')
setwd('Venn/')

library(stringr)
library(ggvenn)
library(grDevices)
library(ggforce)
source(helper_file("plot_Vennplot.R"), encoding = "UTF-8")

a=read.csv("../DEG_PC9/Table/DEG.csv")
b=read.csv("../DEG_HCC827/Table/DEG.csv")
c=read.csv("../DEG_H1975//Table/DEG.csv")

gene_list=list("PC9"=a[a$group!='no',]$name,"HCC827"=b[b$group!='no',]$name,"H1975"=c[c$group!='no',]$name)
intersect_gene=Reduce(intersect,gene_list)
write.csv(intersect_gene,"交集基因.csv")

p1<-plot_Vennplot(data=gene_list,
                  feature = names(gene_list),
                  colors = colors,
                  feature_type = "基因",
                  default_colors=T)

pdf("Venn.pdf",width = 8,height = 6)
print(p1)
dev.off()

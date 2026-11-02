# Module 02: disease_deg
# Source snapshot: archive/original/2.DEG.R
# Path adaptation only; scientific settings and historical seeds are retained.
# Run using: Rscript scripts/run_module.R 02
if (!exists("project_dir", mode = "function")) {
  stop("Run this module through scripts/run_module.R from the repository root.")
}

setwd(project_dir())
dir.create('2.DEG')
setwd('2.DEG/')

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
####数据准备####
exp1<-read.csv("../0.data/exp_combat.csv",row.names = 1,check.names = F)#设置了check.names=F就不会把列名里的-换成.了
group_list<-read.csv("../0.data/groupinfo2.csv",row.names = 1)
colnames(exp1) <- gsub('[.]','-',colnames(exp1))
table(colnames(exp1)==group_list$sample)
exp1=exp1[,group_list$sample]

group1=group_list$group
names(group1)=colnames(exp1)
table(group1)
# group1
# control   tumor 
# 36     569
colors <- c("#1E90FF","#F08080","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
####参数####
exp=exp1
group=group1
case="tumor"
control="control"
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
# down   no   up 
# 1677 7238 1628 
####报告####
source(helper_file("差异分析中文报告.R"), encoding = "UTF-8")




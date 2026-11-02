setwd('/JJGene/data/xuwj/麦门冬汤+五味消毒饮_非小细胞肺癌_EGFR-TKI耐药_5+/')
dir.create('6.亚型分析')
setwd('6.亚型分析/')
source("/opt/scripts_repo/bin/BestCC.r")
data <- data.table::fread("../0.data/exp_combat.csv",data.table = F,check.names = F) %>% tibble::column_to_rownames("V1")
pheno <- read.csv('../4.machine_learning/ML_Gene.csv',header = T)[,1]
total_info <- read.csv('../0.data/groupinfo2.csv',row.names = 1)
disease <- total_info[total_info$group=='tumor',]
data <- data[intersect(rownames(data),pheno),colnames(data) %in% disease$sample] %>% as.matrix()

## 遍历所有可能的组合方式
plotList <- ConsensusD(data,type = "gene") # 矩阵为富集分数时type="score"，矩阵为基因表达量时type = "gene"

ConsensusF(data,best_dist = "pearson",best_cluster = "pam",bestK=2,type = "gene")

data <- data.table::fread("../0.data/exp_combat.csv",data.table = F,check.names = F) %>% tibble::column_to_rownames("V1")
data <- data[,colnames(data) %in% disease$sample]
group_list <- read.table('data/consensus_cluster_group.txt',header = T)
table(colnames(data)==group_list$sample)
####设置工作路径####
# setwd("DEG")
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
source("/opt/scripts_repo/bin/pheatmap1.R")
source("/opt/scripts_repo/bin/plot_boxplot.R")
source("/opt/scripts_repo/bin/plot_vioplot.R")
####数据准备####
exp1<-data#设置了check.names=F就不会把列名里的-换成.了
group_list<-read.table('data/consensus_cluster_group.txt',header = T)
group_list$group <- ifelse(group_list$group==1,'subtype1','subtype2')
table(colnames(exp1)%in%group_list$sample)
exp1=exp1[,group_list$sample]

group1=group_list$group
names(group1)=colnames(exp1)
table(group1)
# group1
# subtype1 subtype2 
# 319      250   

colors <- c("#1E90FF","#F08080","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
scales::show_col(colors)
####参数####
exp=exp1
group=group1
case="subtype2"
control="subtype1"
P="p.adj"
P_value=0.05
logFC_value=0.5
colors=colors
is_count=F
feature_type="基因"
lab = c("gene", "Expression Level")
flip = F
####运行####
source("/opt/scripts_repo/bin/差异分析.R")
# down    no    up 
# 89 10412    42 


####报告####
source("/opt/scripts_repo/bin/差异分析中文报告.R")


####R包和函数加载####
library(ggplot2)
library(grDevices)
library(clusterProfiler)
library(GOplot) 
library(dplyr)
library(ggrepel)
library(stringr)
library(ggpubr)
library(ggprism)
library(org.Hs.eg.db)
# library(org.Mm.eg.db)#当对小鼠数据进行富集分析时加载这个包
#library(org.Rn.eg.db)#当对大鼠数据进行富集分析时加载这个包
#install.packages("/opt/scripts_repo/bin/KEGG.db_1.0.tar.gz",
#                 repos=NULL,type="source") #第一次运行时本地安装KEGG，可用于人类和小鼠
library(KEGG.db)
source("/opt/scripts_repo/bin/plot_enrich_bar.r")
source("/opt/scripts_repo/bin/myGOBubble.R")
source("/opt/scripts_repo/bin/myGOCircle.R")
source("/opt/scripts_repo/bin/myGOChord.R")
source("/opt/scripts_repo/bin/GO_enrich.R")
source("/opt/scripts_repo/bin/KEGG_enrich.R")
####数据准备####
DEG <- read.csv("DEG/Table/DEG.csv",row.names = 1)
colors <- c("#1E90FF","#F08080","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
genes=DEG$name[DEG$group!='no']
sp <- "org.Hs.eg.db" #"org.Hs.eg.db"，"org.Mm.eg.db","org.Rn.eg.db"
gene_type="SYMBOL" #"ENSEMBL","SYMBOL","ENTREZID","ACCNUM"
sig=T#是否区分上下调基因进行富集分析

####运行####
#GO富集
GO_enrich(genes=DEG$name[DEG$group!='no'],
          sp="org.Hs.eg.db" ,
          gene_type="SYMBOL",
          sig=F,colors=colors,
          DEG=DEG)
#KEGG富集
KEGG_enrich(genes=DEG$name[DEG$group!='no'],
            sp="org.Hs.eg.db" ,
            gene_type="SYMBOL",
            sig=F,colors=colors,
            DEG=DEG)

##########预后#############################
file_input <- '../0.data/exp_combat.csv'
file_group <- './data/consensus_cluster_group.txt'
file_hubGene <- '../4.machine_learning/ML_Gene.csv'
library(dplyr)
#设置了check.names=F就不会把列名里的-换成.了
data<-data.table::fread(file_input,data.table =F,check.names = F) %>% 
  tibble::column_to_rownames("V1")

total_info <- read.csv('../0.data/groupinfo2.csv',row.names = 1)
disease <- total_info[total_info$group=='tumor',]
data <- data[,colnames(data) %in% disease$sample] %>% as.matrix()

group_list<-read.table(file_group,header = T)
group_list$group <- ifelse(group_list$group==1,'subtype1','subtype2')
rownames(group_list) <- group_list$sample

# group_list$group <- factor(group_list$group,labels = c('0','1'),levels = c('healthy','sepsis'))
table(colnames(data)%in%group_list$sample)
exp <-cbind(group_list$group,t(data))
exp <- as.data.frame(exp)
names(exp$V1) <- "group"
exp$X <- rownames(exp)
exp$X <- substr(exp$X,1,12) 
duplicates <- duplicated(exp$X) | duplicated(exp$X, fromLast = TRUE)
exp <- exp[!duplicates,]
rownames(exp) <- exp$X
exp <- exp[,names(exp) != "X"]

total_cli <- read.csv("../0.data/total_clinical.csv",header = T,row.names = 1)

data1 <- cbind(total_cli[intersect(rownames(exp),rownames(total_cli)),c(4,5)],exp[intersect(rownames(exp),rownames(total_cli)),])
data_ABCG2 <-data1[,c("fustat","futime","V1","ABCG2")]
data_EPAS1 <-data1[,c("fustat","futime","V1","EPAS1")]
data_PDE9A <-data1[,c("fustat","futime","V1","PDE9A")]

write.table(data1,"./total.txt",sep="\t",row.names=F,quote=F)
write.table(data_ABCG2,"./data_ABCG2.txt",sep="\t",row.names=F,quote=F)
write.table(data_EPAS1,"./data_EPAS1.txt",sep="\t",row.names=F,quote=F)
write.table(data_PDE9A,"./data_PDE9A.txt",sep="\t",row.names=F,quote=F)

bioSurvival=function(inputFile=null,outFile=null){
  #读取输入文件
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  #比较高低风险组生存差异，得到显著性p值
    diff=survdiff(Surv(futime, fustat) ~V1,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.05){
    pValue="p<0.05"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ V1, data = rt)
  
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=3,
                     legend.title="subtype",
                     legend.labs=c("subtype2", "subtype1"),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c(colors[2], colors[1]),
                     risk.table=FALSE,
                     risk.table.title="",
                     risk.table.height=.25,
                     fontsize = 2,
                     ggtheme = theme(axis.title = element_text(size = 7),
                                     axis.text = element_text(size = 6,colour = "black"),
                                     panel.background = element_blank(),
                                     axis.line = element_line(size = 0.5),
                                     legend.key.size = unit(0.5,"cm"),
                                     legend.text = element_text(size = 7),
                                     legend.title = element_blank())
  )
  while (!is.null(dev.list()))  dev.off()
  pdf(file=outFile,onefile = FALSE,width = 8/2.54,height =6/2.54)
  par(family="sans",ps="7")
  print(surPlot)
  dev.off()
}

library(glmnet) ##Lasso回归
library(survival)
library(survminer)
library(survminer)
library(timeROC)
library(pheatmap)

bioSurvival(inputFile= "./total.txt", outFile="./totalSurv.pdf")
bioSurvival(inputFile= "./data_ABCG2.txt", outFile="./ABCG2Surv.pdf")
bioSurvival(inputFile= "./data_EPAS1.txt", outFile="./EPAS1Surv.pdf")
bioSurvival(inputFile= "./data_PDE9A.txt", outFile="./PDE9ASurv.pdf")
########肿瘤转移##########

group_list$sample <- substr(group_list$sample,1,12)
colnames(group_list) <- c("X","group")
total_cli$X <- rownames(total_cli)
total <- merge(total_cli,group_list,by = "X")
total <- total[,c(1,5,8)]
duplicates <- duplicated(total$X) | duplicated(total$X, fromLast = TRUE)
total <- total[!duplicates,]
rownames(total) <- total$X
total <- total[,-1]

library(ggplot2)
library(ggalluvial)

df <- to_lodes_form(total[,1:ncol(total)],
                    axes = 1:ncol(total),
                    id = "value")
colors <- c('#4DBBD5','#E64B35','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148','#B09C85')
p <- ggplot(df, aes(x = x, fill=stratum, label=stratum,
                    stratum = stratum, alluvium  = value))+#数据
  geom_flow(width = 0.3,#连线宽度
            curve_type = "sine",#曲线形状，有linear、cubic、quintic、sine、arctangent、sigmoid几种类型可供调整
            alpha = 0.5,#透明度
            color = 'white',#间隔颜色
            size = 0.1)+#间隔宽度
  geom_stratum(width = 0.28)+#图中方块的宽度
  geom_text(stat = 'stratum', size = 2, color = 'black')+
  scale_fill_manual(values = colors)+#自定义颜色
  theme_void()+#主题（无轴及网格线）
  theme(legend.position = 'none')#去除图例
pdf("./sankeyNetwork.pdf",width = 17/2.5,height = 7/2.5)
print(p)
while (!is.null(dev.list()))  dev.off()



###############单基因KM################
library(survival)
library(dplyr)
library("survminer")


total_clic <- read.csv('../0.data/total_clinical.csv',row.names = 1) %>% na.omit()
rownames(total_clic) <- total_clic$id
exp <- read.csv('../0.data/exp_combat.csv',row.names = 1,check.names = F) %>% t() %>% as.matrix()
colnames(total_clic)[1] <- 'sample'
hubgene <- read.csv('../4.machine_learning/ML_Gene.csv')[[1]]
data <- exp[rownames(exp) %in% total_cli$id,hubgene]
table(rownames(data) %in% total_clic$sample)
dat_surv <- data %>% as.data.frame()
dat_surv$fustat <- total_clic[rownames(data),]$fustat
dat_surv$futime <- total_clic[rownames(data),]$futime

event <- "OS" ## 生存事件
fustat <- "fustat" ## 生存时间所在列列名
futime <- "futime" ## 生存时间所在列列名
surv_objects <- hubgene ## 生存分析变量，可以是临床信息也可以是基因信息
#"RdBu","npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty"
# 或者自定义颜色类型c(col1,col2)
color_plate <- "nejm"  ## 颜色参数

source("/opt/scripts_repo/bin/easySurv.r")
surv_res <- easySurv(dat_surv = dat_surv,surv_objects = surv_objects,
                     futime = futime,fustat=fustat,event = event,
                     palette = color_plate)

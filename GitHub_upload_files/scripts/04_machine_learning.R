# Module 04: machine_learning
# Source snapshot: archive/original/4.机器学习.R
# Path adaptation only; scientific settings and historical seeds are retained.
# Run using: Rscript scripts/run_module.R 04
if (!exists("project_dir", mode = "function")) {
  stop("Run this module through scripts/run_module.R from the repository root.")
}

setwd(project_dir())
dir.create('4.machine_learning')
## 读取数据
file_input <- '../0.data/exp_combat.csv'
file_group <- '../0.data/groupinfo2.csv'
file_hubGene <- '../3.Venn/交集基因.csv'
library(dplyr)
#设置了check.names=F就不会把列名里的-换成.了
exp<-data.table::fread(file_input,data.table =F,check.names = F) %>% 
  tibble::column_to_rownames("V1")
colnames(exp)<-gsub('[.]','-',colnames(exp))
group_list<-read.csv(file_group,header = T,row.names = 1)
# group_list$group <- ifelse(group_list$group==1,'subtype1','subtype2')
# group_list$group <- factor(group_list$group,labels = c('0','1'),levels = c('healthy','sepsis'))
table(colnames(exp)%in%group_list$sample)
exp1=exp[,group_list$sample]


control <- "control"
case <- "tumor"
target <- "group"

hubgene <- read.csv("../3.Venn_1/交集基因.csv",header = T,row.names = 1)[,1]
input <- cbind(data.frame(group=ifelse(group_list$group=="control",0,1)),t(exp1[hubgene,]))
input$group <-  factor(input$group,levels = c(0,1))
str(input)

source(helper_file("ML_FeatureSelect.r"), encoding = "UTF-8")
# 生成10000范围内20个随机数种子
set.seed(66)
seeds <- sample(1:10000, 20)

# seeds <- sample(1:100)
## 由于SVM-RFE分析会依据数据分配的不同改变对应的准确率，所以在此处遍历随机数种子选择准确率最低的即可
result <- svm_seed_select(input, target = target, positive = '1',msr = "classif.logloss", seeds)

# result中最好的随机数种子为9191，这个无法用程序分析，只能人为去观察 
setwd('4.machine_learning/')
resSVM <- mySVM_RFE(input,target = target,positive = "1",msr = "classif.logloss",seed = 9191)
resRF <- myRF(input,target = target,positive = "1",seed = 9191,top_number = 5,max_tree = 1000,nfold = 5,
              threads = 30,control=control,case=case,loss_method = "classif.logloss")
resLasso <- myLasso(input,target =target,seed = 9191)

####R包和函数加载####
library(stringr)
library(ggvenn)
library(grDevices)
library(ggforce)
source(helper_file("plot_Vennplot.R"), encoding = "UTF-8")
colors <- c("#1E90FF", "#E64B35","#00A087", "#3C5488", "#F39B7F","#FFFF00")

ML_list <- list(
  "Lasso" = resLasso$features,
  "SVM-RFE" = resSVM$features
  # "RandomForest" = resRF$features
)
saveRDS(ML_list,'ML_list.rds')

ml_genes <- Reduce(function(x,y) intersect(x,y),ML_list)
if(length(ml_genes) != 0){
  p1<-plot_Vennplot(data=ML_list,
                    colors = colors,
                    feature_type = "基因",
                    default_colors=T)
  
  pdf("machina_learning_venn.pdf",width = 8,height = 6)
  print(p1)
  dev.off()
  
}else{
  warning("数据的数据集之间没有交集！")
}

write.csv(ml_genes,"ML_Gene.csv",row.names = F,quote = F)

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
#library(org.Mm.eg.db)#当对小鼠数据进行富集分析时加载这个包
#library(org.Rn.eg.db)#当对大鼠数据进行富集分析时加载这个包
#install.packages("/opt/scripts_repo/bin/KEGG.db_1.0.tar.gz",
#                 repos=NULL,type="source") #第一次运行时本地安装KEGG，可用于人类和小鼠
library(KEGG.db)
source(helper_file("plot_enrich_bar.r"), encoding = "UTF-8")
source(helper_file("myGOBubble.R"), encoding = "UTF-8")
source(helper_file("myGOCircle.R"), encoding = "UTF-8")
source(helper_file("myGOChord.R"), encoding = "UTF-8")
source(helper_file("GO_enrich.R"), encoding = "UTF-8")
source(helper_file("KEGG_enrich.R"), encoding = "UTF-8")
####数据准备####
DEG <- read.csv("../2.DEG/DEG/Table/DEG.csv",row.names = 1)
colors <- c('#0073C2','#EFC000','#CD534C','#7AA6DC','#003C67','#8F7700','#3B3B3B','#A73030','#4A6990')
genes=ml_genes
sp <- "org.Hs.eg.db" #"org.Hs.eg.db"，"org.Mm.eg.db","org.Rn.eg.db"
gene_type="SYMBOL" #"ENSEMBL","SYMBOL","ENTREZID","ACCNUM"
sig=T#是否区分上下调基因进行富集分析

####运行####
#GO富集
GO_enrich(genes=ml_genes,
          sp="org.Hs.eg.db" ,
          gene_type="SYMBOL",
          sig=F,colors=colors,
          DEG=DEG)
#KEGG富集
KEGG_enrich(genes=ml_genes,
            sp="org.Hs.eg.db" ,
            gene_type="SYMBOL",
            sig=F,colors=colors,
            DEG=DEG)

DEG_name <- '关键基因'

source(helper_file("10_GO_KEGG.R"), encoding = "UTF-8")



##############################表达验证###########################################
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(pROC)
#devtools::install_github("thomasp85/patchwork@v1.1.0")
library(patchwork)
library(ggplotify)
source(helper_file("pheatmap1.R"), encoding = "UTF-8")
source(helper_file("plot_boxplot.R"), encoding = "UTF-8")
source(helper_file("plot_vioplot.R"), encoding = "UTF-8")

exp<-data.table::fread('../0.data/exp_combat.csv',data.table =F,check.names = F) %>% 
  tibble::column_to_rownames("V1")
colnames(exp) <- gsub("[.]","-",colnames(exp))
group_list<-read.csv('../0.data/groupinfo2.csv',row.names = 1,header = T)
# group_list$group <- ifelse(group_list$group==1,'subtype1','subtype2')
# group_list$group <- factor(group_list$group,labels = c('0','1'),levels = c('healthy','sepsis'))
table(colnames(exp)%in%group_list$sample)
exp1=exp[,group_list$sample]
hubgene <- read.csv("./ML_Gene.csv" ,header = T)[,1]


control <- "control"
case <- "tumor"
target <- "group"
colors <- c("#1E90FF", "#E64B35","#00A087", "#3C5488", "#F39B7F","#FFFF00")

input<-t(exp1[hubgene,]) %>% as.data.frame()

group=group_list$group
names(group)=group_list$sample
feature=colnames(input)

data=as.data.frame(input)
rownames(data) <- gsub('[.]','-',rownames(data))
colors <- colors
feature_type="基因"
lab = c("Gene", "Expression Level")
flip = F

####运行####
source(helper_file("Gene_expression.R"), encoding = "UTF-8")

# 关键基因间相关性####
library(ggcorrplot)
library(ggthemes)

# 输入文件
data<-input
# data <- data[-6,]
# 计算基因与基因之间的相关性
corr <- round(cor(input),3)
# 检验基因与基因之间相关性P值
p.mat <- round(cor_pmat(data),3)
# 绘制基因与基因间相关性热图
ggcorrplot(corr,method = "circle",                  #绘制环形相关性热图
           hc.order = TRUE, hc.method = "ward.D",   #对基因进行聚类，聚类算法为ward.D
           outline.color = "white",
           ggtheme = theme_bw(),
           # type = "upper",                          #只绘制相关性上部分
           # colors = c("#6D9EC1","white","#E46726"), #更改相关性热图颜色
           colors = c(colors[1],"white",colors[2]), #更改相关性热图颜色
           lab = TRUE, lab_size = 2,
           p.mat = p.mat,                           #显示相关性系数
           insig = "blank"                          #隐藏p>0.05的基因
)
ggsave("hubgene_corr_heatmap.pdf",width = 6,height = 5)

##############################外部验证###########################################
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(pROC)
#devtools::install_github("thomasp85/patchwork@v1.1.0")
library(patchwork)
library(ggplotify)
source(helper_file("pheatmap1.R"), encoding = "UTF-8")
source(helper_file("plot_boxplot.R"), encoding = "UTF-8")
source(helper_file("plot_vioplot.R"), encoding = "UTF-8")

exp<-data.table::fread('../0.data/GSE54236.csv',data.table =F,check.names = F) %>% 
  tibble::column_to_rownames("V1")
colnames(exp) <- gsub("[.]","-",colnames(exp))
group_list<-read.table('../0.data/GSE54236_group.txt',header = T)
# group_list$group <- ifelse(group_list$group==1,'subtype1','subtype2')
# group_list$group <- factor(group_list$group,labels = c('0','1'),levels = c('healthy','sepsis'))
table(colnames(exp)%in%group_list$sample)
exp1=exp[,group_list$sample]
hubgene <- read.csv("./ML_Gene.csv", header = T)[,1]


control <- "control"
case <- "tumor"
target <- "group"
colors <- c("#1E90FF", "#E64B35","#00A087", "#3C5488", "#F39B7F","#FFFF00")

input<-t(exp1[hubgene,]) %>% as.data.frame()

group=group_list$group
names(group)=group_list$sample
feature=colnames(input)

data=as.data.frame(input)
rownames(data) <- gsub('[.]','-',rownames(data))
colors <- colors
feature_type="基因"
lab = c("Gene", "Expression Level")
flip = F

####运行####
source(helper_file("Gene_expression.R"), encoding = "UTF-8")



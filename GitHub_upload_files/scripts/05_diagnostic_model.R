# Module 05: diagnostic_model
# Source snapshot: archive/original/5.诊断模型.R
# Path adaptation only; scientific settings and historical seeds are retained.
# Run using: Rscript scripts/run_module.R 05
if (!exists("project_dir", mode = "function")) {
  stop("Run this module through scripts/run_module.R from the repository root.")
}

setwd(project_dir())
dir.create('5.norma')
setwd('5.norma/')


exp<-read.csv('../0.data/exp_combat.csv',row.names = 1)
group_list<-read.csv('../0.data/groupinfo2.csv',row.names = 1)
hubgene <- read.csv("../4.machine_learning/ML_Gene.csv", header = T)[,1]
####参数####
data=cbind(group_list[,"group"],t(exp))
colnames(data)[1] <- 'group'
data <- data[,c('group',hubgene)] %>% as.data.frame()

write.table(data,"./nomo.txt")
data <- read.table("./nomo.txt", header = T, row.names = 1)

data$group=factor(data$group,levels=c("control","tumor"),labels = c(0,1))

for (i in 2:ncol(data)) {
  data[,i] <- as.numeric(data[,i])
}
#当不是连续变量时需要设为因子格式的二分类变量，并且降低风险的类型应该放在前面（如data$risk=factor(data$risk,levels=c("Low","High"))）
feature <- hubgene
colors <- c("#0000CD",'#DC0000',"#FFFF00",'#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#7E6148','#B09C85')
nomo_n=605

diagnosis=T
# time_use="futime"
event_use="group"

####运行####

.libPaths("/home/BI26/R/x86_64-pc-linux-gnu-library/4.1")
source(helper_file("列线图.R"), encoding = "UTF-8")


# exp1<-read.csv('../data/GSE221161_normalize.csv',row.names = 1)
# group_list1<-read.csv('../data/GSE221161_group.csv',row.names = 1)
# hubgene<-read.csv('../3.machine_learning/ML_Gene.csv',header = 1)[,1]
# ####参数####
# data_valid=cbind(group_list1[,"group"],t(exp1))
# colnames(data_valid)[1] <- 'group'
# data_valid <- data_valid[,c('group',hubgene)] %>% as.data.frame()
# data_valid$group=factor(data_valid$group,levels=c("Normal","CD"),labels = c(0,1))
# for (i in 2:ncol(data_valid)) {
#   data_valid[,i] <- as.numeric(data_valid[,i])
# }
# #feature就是上面的feature
# Nomogram_valid(data = data_valid,feature = feature)
# #结果会保存在Nomogram/plot/ROC_valid.pdf


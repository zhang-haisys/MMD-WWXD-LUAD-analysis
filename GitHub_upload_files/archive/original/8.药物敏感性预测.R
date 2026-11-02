##待预测的表达量矩阵
testExpr <- read.csv("../0.data/exp_combat.csv",row.names = 1,check.names = F)
group <- read.table("../6.亚型分析_1/data/consensus_cluster_group.txt",sep = "\t",header = T)
group$group <- ifelse(group$group == 1 ,"subtype1","subtype2")
write.csv(group,"../0.data/consensus_cluster_group.csv",row.names = F)
testExpr <- testExpr[,group$sample]
testExpr[1:4,1:4]  
dim(testExpr)
testExpr <- as.matrix(testExpr)

library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))

##加载训练集的表达量矩阵和药物处理信息
dir='/work/database/oncoPredict_Training_Data/'
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

##运行oncoPredict
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData')


library(data.table)
testPtype <- fread('../8.Phenotype/DrugPredictions.csv', data.table = F) %>% column_to_rownames("V1")
rownames(testPtype) <- gsub("[.]","-",rownames(testPtype))
testPtype[1:4, 1:4]
data <- merge(group %>% column_to_rownames(var = "sample"),testPtype,by = "row.names") %>% 
  column_to_rownames("Row.names")

##绘图
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
# colors <- c("#4DBBD5", "#E64B35","#00A087", "#3C5488", "#F39B7F")
sample_size = data %>% group_by(group) %>% summarise(num=n())
means <- colMeans(data[,2:199]) %>% as.data.frame() %>% rownames_to_column(var = "drug") 
colnames(means)[2] <- "mean"

class1 <- data[which(data$group == "subtype2"),]
class2 <- data[which(data$group == "subtype1"),]

dir.create("vlnplot")
setwd("./vlnplot")
p_value <- c()
for (i in colnames(testPtype)) {
  test <- wilcox.test(class1[,i],class2[,i])
  if (test$p.value < 0.05) {
    P <- data %>%
      left_join(sample_size) %>%
      mutate(myaxis = paste0(group, "\n", "n=", num)) %>%
      ggplot(aes(x=myaxis, y=data[,i], fill=group)) +
      geom_violin(width=1) +
      geom_boxplot(width=0.1, color="grey", alpha=0.2) +
      scale_fill_manual(values = c(colors[1],colors[2])) +
      theme(
        legend.position="none",
        plot.title = element_text(size=8),
        panel.background = element_blank(),
        axis.line = element_line(size = 0.5),axis.text = element_text(size=7),
        axis.title = element_text(size=8)
      ) +
      labs(x="",y=i)+
      stat_compare_means(aes(group = group),method = "wilcox.test",size = 3)
    
    pdf(paste(i,"_vln.pdf",sep = ""),width = 3,height = 3)
    print(P)
    dev.off()
    p_value <- c(p_value,i)
  }
}
setwd("../")
means$up <- ifelse(colMeans(class1[,2:199])>colMeans(class2[,2:199]),"subtype2","subtype1")
means <- means %>% filter(drug %in% p_value)
write.csv(means,"../8.Phenotype_1/vlnplot/means.csv",row.names = F)

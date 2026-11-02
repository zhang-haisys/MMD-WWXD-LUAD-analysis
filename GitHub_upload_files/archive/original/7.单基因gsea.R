library(ggcorrplot)
library(ggthemes)
exp <- read.csv("../0.data/exp_combat.csv",row.names = 1,check.names = F)
group_list <- read.table("../6.亚型分析_1/data/consensus_cluster_group.txt",sep = "\t",header = T)
exp1 <- exp[,group_list$sample]
hubgene <- read.csv("../4.machine_learning_1/ML_Gene.csv")[,1]
input<-t(exp1[hubgene,]) %>% as.data.frame()
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
ggsave("../6.亚型分析_1/hubgene_corr_heatmap.pdf",width = 6,height = 5)


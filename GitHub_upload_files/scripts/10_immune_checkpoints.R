# Module 10: immune_checkpoints
# Source snapshot: archive/original/10.免疫检查点.R
# Path adaptation only; scientific settings and historical seeds are retained.
# Run using: Rscript scripts/run_module.R 10
if (!exists("project_dir", mode = "function")) {
  stop("Run this module through scripts/run_module.R from the repository root.")
}

source(helper_file("pheatmap1.R"), encoding = "UTF-8")
source(helper_file("plot_boxplot.R"), encoding = "UTF-8")
source(helper_file("plot_vioplot.R"), encoding = "UTF-8")
library(pheatmap)
# colors <- c("#4DBBD5", "#E64B35","#00A087", "#3C5488", "#F39B7F")

gene <- c("CD274","PDCD1","CTLA4","TIGIT","CD47","HAVCR2","LAG3","CD24","SIGLEC10","TNFRSF4","41BB","CD28")
exp <- read.csv("../0.data/exp_combat.csv",check.names = F,row.names = 1)
exp1 <- exp[gene,] %>% na.omit()
group_list <- read.table("../6.亚型分析_1/data/consensus_cluster_group.txt",header = T,sep = "\t") 
colnames(group_list) <- c("sample","group")

group_list$group <- ifelse(group_list$group == 1,"subtype1","subtype2")
write.csv(group_list,"../9.tmb/Mutation/consensus_cluster_group.csv",row.names = F)
dat <- exp1 %>% t() %>% as.data.frame() %>% rownames_to_column(var = "sample") %>% 
  merge(.,group_list,by = "sample") %>% reshape2::melt() %>% na.omit()
dir.create("../10.immunecheckpoint_1/")
immune <- ggplot(dat,aes(variable,value,fill = group)) + 
  geom_boxplot(outlier.shape = 21,outlier.color = NA,outlier.size = 0.8,
               #notch=TRUE,notchwidth=0.8, #在箱子上生成槽口,notchwidth越小则越往里凹
               alpha = 0.7) + 
  # stat_summary(fun.y="mean",geom="point",shape=23,size=1,fill="white",
  #              position = position_dodge(0.8))+ #添加均值标记
  labs(x = "", y = "Expression Level") +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  stat_compare_means(aes(group = group,label = ..p.signif..),bracket.size = 0.6, size = 4,
                     label.y = max(dat$value)-0.7, #p值位置
                     hide.ns = T,#隐藏ns
                     method = ifelse(length(unique(dat$group))==2,"wilcox.test","kruskal.test"))+
  theme(plot.title = element_blank(),plot.subtitle = element_blank(),
        plot.background = element_blank(),plot.margin = margin(t=2,r=2,b=2,l=2,unit="pt"),
        panel.border = element_rect(size = 0.6,fill = NA,colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1,size = 7,colour = "black"),
        axis.title = element_text(size = 7,colour = "black"),
        axis.text.y = element_text(size = 7,colour = "black"),
        legend.background = element_blank(),legend.key = element_blank(),
        legend.title = element_text(size = 8),legend.text = element_text(size = 8),
        legend.position = "top",legend.margin = margin(t=0,r=0,b=0,l=0,unit="pt"),legend.box.spacing = unit(3,"pt"))#+coord_flip() #旋转箱线图
ggsave(plot = immune,"../10.immunecheckpoint_1/boxplot.pdf",width = 10,height = 7,units = "cm")




####TIDE####
library()
dir.create('../11.TIDE')
tumor_exp <- read.csv('../0.data/exp_combat.csv',row.names = 1)
colnames(tumor_exp)<-gsub('[.]','-',colnames(tumor_exp))
total_risk<-read.table('../6.亚型分析_1/data/consensus_cluster_group.txt',sep = "\t",header = T,row.names = 1)
total_risk$group <- ifelse(total_risk$group == 1, "subtype1", "subtype2")
tumor_exp <- tumor_exp[,rownames(total_risk)]
Exp <- t(apply(tumor_exp, 1, function(x)x-(mean(x)))) %>% as.data.frame() %>% rownames_to_column(var = "sample")
write.table(Exp,'../11.TIDE/Exp.txt',sep = "\t",row.names = F)
TIDE_result <- read.csv('../11.TIDE/TIDE.csv',row.names = 1)
Trisk <- total_risk %>% dplyr::select(c("group"))
rownames(Trisk) <- rownames(total_risk)
colnames(Trisk) <- "group"
TIDE_result$sample <- rownames(TIDE_result)
Trisk$sample <- rownames(Trisk)
TIDE_result1 <- left_join(TIDE_result,Trisk,by='sample')
rownames(TIDE_result1) <- TIDE_result1$sample
TIDE_result1 <- TIDE_result1[,c(1:14,16)]

####箱线图
TIDE_result2 <- TIDE_result1[,c(3,10,11)]
dat <- TIDE_result2 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
colnames(Trisk) <- c('group','Sample')
dat1 <- left_join(dat,Trisk,by = 'Sample')
sample1 <- Trisk %>% dplyr::filter(group %in% "subtype1")
sample2 <- Trisk %>% dplyr::filter(group %in% "subtype2")

deg_boxplot=ggplot(dat1,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,outlier.color = NA,outlier.size = 0.8,
               #notch=TRUE,notchwidth=0.8, #在箱子上生成槽口,notchwidth越小则越往里凹
               alpha = 0.7) + 
  # stat_summary(fun.y="mean",geom="point",shape=23,size=1,fill="white",
  #              position = position_dodge(0.8))+ #添加均值标记
  labs(x = "", y = "Score") +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  stat_compare_means(aes(group = group,label = ..p.signif..),bracket.size = 0.6, size = 3,
                     label.y = max(dat1$Proportion)-0.2, #p值位置
                     hide.ns = T,#隐藏ns
                     method = ifelse(length(unique(dat1$group))==2,"t.test","kruskal.test"))+
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
        legend.position = "top",legend.margin = margin(t=0,r=0,b=0,l=0,unit="pt"),legend.box.spacing = unit(3,"pt"))

pdf("../11.TIDE/TIDE_boxplot.pdf",width = 8/2.5,height = 7/2.5)
print(deg_boxplot)
while (!is.null(dev.list()))  dev.off()


# colors <- c('#EFC000','#0073C2','#868686','#CD534C','#7AA6DC','#003C67','#8F7700','#3B3B3B','#A73030','#4A6990')
#桑基图
library(ggplot2)
library(ggalluvial)
TIDE_result3 <- TIDE_result1[,c(15,2)]
# TIDE_result3$risk <- ifelse(TIDE_result3$risk=="high","High risk","Low risk")
df <- to_lodes_form(TIDE_result3[,1:ncol(TIDE_result3)],
                    axes = 1:ncol(TIDE_result3),
                    id = "value")


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
pdf("../11.TIDE/TIDE_sankeyNetwork.pdf",width = 17/2.5,height = 7/2.5)
print(p)
while (!is.null(dev.list()))  dev.off()



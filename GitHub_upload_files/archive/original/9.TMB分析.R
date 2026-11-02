library(TCGAbiolinks)
setwd('../0.data/Mutation_1/')
query_exp = GDCquery(project = "TCGA-LUAD",
                     data.category = "Simple Nucleotide Variation",
                     data.type="Masked Somatic Mutation_1")

GDCdownload(query = query_exp, method = "api", directory = "TCGA_LUAD_SNV", files.per.chunk = 10)
library(maftools)
library(ggplot2)

setwd("../9.tmb/Mutation_1/")
filepath = dir(path = '../Mutation_1/TCGA_LUAD_SNV/TCGA-LUAD/harmonized/Simple_Nucleotide_Variation/Masked_Somatic_Mutation_1/',
               pattern = "maf.gz$",
               full.names = T,
               recursive = T)

x = lapply(filepath, data.table::fread, skip = "Hugo_Symbol")
x = data.table::rbindlist(l = x, use.names = TRUE, fill = TRUE)
luad = maftools::read.maf(maf = x)

luad # 直接输入MAF对象可以查看MAF文件的基本信息
getSampleSummary(luad) # 显示样品的统计
getGeneSummary(luad) # 显示基因的统计
getClinicalData(luad) # 显示样品关联的临床数据
getFields(luad) # 显示MAF文件中的所有字段
write.mafSummary(maf=luad, basename="../Mutation_1/luad") # 将详细统计结果输出到文件
pdf('../Mutation_1/mafSummary.pdf',width = 10,height = 10)
plotmafSummary(maf=luad, rmOutlier=TRUE, addStat="median", dashboard=TRUE, titvRaw = FALSE)
dev.off()

total_risk<-read.csv('../Mutation/consensus_cluster_group.csv',row.names = 1)
maf<-read.maf('../Mutation_1/luad_maftools.maf')
dat<-maf@data
sample1=unique(as.character(dat$Tumor_Sample_Barcode))
sample2=substr(sample1,1,12)
hname=substr(rownames(total_risk)[total_risk$group=='subtype1'],1,12)
lname=substr(rownames(total_risk)[total_risk$group=='subtype2'],1,12)

samp=intersect(sample2,hname)
name1=sample1[which(sample2 %in% samp)]
samp=intersect(sample2,lname)
name2=sample1[which(sample2 %in% samp)]

clu1_maf=subsetMaf(maf,tsb = name1)
clu2_maf=subsetMaf(maf,tsb = name2)
maf1=subsetMaf(maf,tsb = c(name1,name2))

save(maf1,file = "../Mutation_1/maf1.rda")
save(clu1_maf,clu2_maf,file="../Mutation_1/maf_clss.rda")

#挑选总样本中发生突变的前20个基因
a<-maf1@data
b<-unique(a[,c('Hugo_Symbol','Tumor_Sample_Barcode')])
times<-as.data.frame(table(b$Hugo_Symbol))
times<-times[order(times$Freq,decreasing = T),]
top20_gene<-as.character(times$Var1[1:20])

pdf('../Mutation_1/maf.pdf',width = 7,height = 7)
oncoplot(maf=maf1, borderCol=NULL,gene_mar = 8,genes = top20_gene,keepGeneOrder = T)
dev.off()

pdf('../Mutation_1/maf_subtype1.pdf',width = 7,height = 7)
oncoplot(maf=clu1_maf, borderCol=NULL,gene_mar = 8,genes = top20_gene,keepGeneOrder = T)
dev.off()

pdf('../Mutation_1/maf_subtype2.pdf',width = 7,height = 7)
oncoplot(maf=clu2_maf, borderCol=NULL,gene_mar = 8,genes = top20_gene,keepGeneOrder = T)
dev.off()

####肿瘤突变负荷计算和绘制箱线图####
#maftools版本2.10.05
tmb1<-tmb(maf = maf1)

write.csv(tmb1,'../Mutation_1/tmb.csv')
tmb1 <- read.csv("../Mutation_1/tmb.csv",row.names = 1) %>% dplyr::rename(sample = Tumor_Sample_Barcode)
tmb1$sample <- substr(tmb1$sample,1,12)
group <- total_risk %>% dplyr::select("group") %>% rownames_to_column(var = "sample")
group$sample <- substr(group$sample,1,12)
tmb1 <- merge(tmb1,group,by = "sample")
# tmb1$risk <- ifelse(tmb1$risk=="high","High risk","Low risk")
ggplot(data = tmb1, aes(x = group, y = total_perMB_log, fill = group))+ 
  scale_fill_manual(values = c(colors[1],colors[2])) +
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = FALSE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  geom_point(shape = 21, size=2,
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 1) +
  theme_classic() + 
  ylab("log10 (Tumor Mutation_1 burden)") +
  xlab("") +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  stat_compare_means(aes(group = group,label = ..p.signif..),bracket.size = 0.6, size = 8,
                     label.x = 1.5,
                     label.y = max(tmb1$total_perMB_log)*0.95+min(tmb1$total_perMB_log)*0.05, #p值位置
                     hide.ns = F,#隐藏ns
                     method = ifelse(length(unique(tmb1$group))==2,"wilcox.test","kruskal.test")) 
ggsave("../Mutation_1/TMB-class_box.pdf",width = 4,height = 4)

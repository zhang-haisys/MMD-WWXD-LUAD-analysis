# Module 03: intersection_enrichment
# Source snapshot: archive/original/3.Venn.R
# Path adaptation only; scientific settings and historical seeds are retained.
# Run using: Rscript scripts/run_module.R 03
if (!exists("project_dir", mode = "function")) {
  stop("Run this module through scripts/run_module.R from the repository root.")
}

setwd(project_dir())
dir.create('3.Venn_1')
setwd('3.Venn_1/')
####取交集####
library(stringr)
library(ggvenn)
library(grDevices)
library(ggforce)
source(helper_file("plot_Vennplot.R"), encoding = "UTF-8")

a=read.csv("../1.Coregene_1/Venn/交集基因.csv")
b=read.table('../0.data/药物靶点.txt',header = T)
c=read.csv("../2.DEG/DEG/Table/DEG.csv")

# intersect(a[a$group!='no',]$name,b$ECM_remodeling)
# intersect(a[a$group!='no',]$name,c$gene)
# intersect(b$ECM_remodeling,c$gene)

gene_list=list("Osimertinib resistance gene"=a$x,"Drug target"=b$x,"DEG of LUAD"=c[c$group!='no',]$name)
intersect_gene=Reduce(intersect,gene_list)
write.csv(intersect_gene,"交集基因.csv")

colors <- c("#1E90FF","#F08080","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")


p1<-plot_Vennplot(data=gene_list,
                  feature = names(gene_list),
                  colors = colors,
                  feature_type = "基因",
                  default_colors=T)

pdf("Venn.pdf",width = 8,height = 6)
print(p1)
dev.off()


####GO_KEGG####
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
genes=intersect_gene
sp <- "org.Hs.eg.db" #"org.Hs.eg.db"，"org.Mm.eg.db"
gene_type="SYMBOL" #"ENSEMBL","SYMBOL","ENTREZID","ACCNUM"
sig=F#是否区分上下调基因进行富集分析
colors <- c("#0000CD",'#DC0000',"#FFFF00",'#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#7E6148','#B09C85')
####运行####
#GO富集
GO_enrich(genes=intersect_gene,
          sp="org.Hs.eg.db" ,
          gene_type="SYMBOL",
          sig=F,colors=colors,
          DEG=DEG)
#KEGG富集
KEGG_enrich(genes=intersect_gene,
            sp="org.Hs.eg.db" ,
            gene_type="SYMBOL",
            sig=F,colors=colors,
            DEG=DEG)




# ########################################################
# {
# DEG_name <- '交集基因'
# load("GO_KEGG/data/GO_list.rda")
# load("GO_KEGG/data/KEGG_list.rda")
# 
# # dir.create('10_GO_KEGG')
# # setwd('10_GO_KEGG')
# ####GO分析中文报告####
# setwd('GO_KEGG')
# dir.create('Report')
# cn_path="/opt/scripts_repo/bin/template_cn.docx"
# 
# #报告文档
# #officer
# #word
# #https://www.rdocumentation.org/packages/officer/versions/0.5.2
# 
# pkgs <- c('officer')
# if(!requireNamespace("officer", quietly = TRUE))install.packages("officer")
# 
# library(officer)
# library(magrittr) #只有加载了这个包才能用管道函数
# library(tidyverse)
# 
# #创建空的word表格
# # 
# # my_doc<- read_docx(cn_path)
# 
# #了解常见排版的参数，这一步可以省略，只是为了更加清楚有哪些参数而已。
# # s=styles_info(my_doc)
# # test <- styles_info(my_doc)#比较重要的有：style_type,style_name
# 
# #设置字体
# # text_font=fp_text(color = 'black',font.size = 12, font.family ='宋体') 
# # title_font=fp_text(color = 'black',font.size = 15, font.family = '宋体',bold = T)
# ## 参数高亮
# ## 人工替换参数高亮 yellow 黄色
# underline_font=fp_text(color = 'black',shading.color='yellow',
#                        font.family = "Times New Roman",hansi.family = 'Times New Roman',
#                        eastasia.family = "宋体",font.size = 10.5)
# 
# ## 注意事项高亮 grey 灰色
# underline_font_grey=fp_text(color = 'black',shading.color='grey',
#                             font.family = "Times New Roman",hansi.family = 'Times New Roman',
#                             eastasia.family = "宋体",font.size = 10.5)
# 
# ## version ,figure 和 table 高亮orange  橘色
# underline_font_orange=fp_text(color = 'black',shading.color='orange',
#                               font.family = "Times New Roman",hansi.family = 'Times New Roman',
#                               eastasia.family = "宋体",font.size = 10.5)
# 
# ## 代码替换参数 #B2CAF2 蓝色
# underline_font_blue=fp_text(color = 'black',shading.color='#B2CAF2',
#                             font.family = "Times New Roman",hansi.family = 'Times New Roman',
#                             eastasia.family = "宋体",font.size = 10.5)
# 
# ##  人工文档中更改的参数  #FF8F70 红色
# underline_font_red=fp_text(color = 'black',shading.color='#FF8F70',
#                            font.family = "Times New Roman",hansi.family = 'Times New Roman',
#                            eastasia.family = "宋体",font.size = 10.5)
# 
# ## 其他待整理参数 green 绿色
# underline_font_green=fp_text(color = 'black',shading.color='#BDE8B9',
#                              font.family = "Times New Roman",hansi.family = 'Times New Roman',
#                              eastasia.family = "宋体",font.size = 10.5)
# 
# ####文档所需参数####
# # #1. 用于GO——KEGG基因分析的名称
# # DEG_name <- '疾病正常差异表达基因'
# #2. 受什么驱动的基因
# genelist_name_ch <- '铂类耐药'
# #3. 疾病名称
# case_name_ch <- '结肠癌'
# 
# #####################################################
# 
# # go 富集结果
# go_top3 <- GO_list$ego_list
# go_top3$LOG10pvalue <- -log10(go_top3$p.adjust)
# go_top3 = go_top3 %>% group_by(ONTOLOGY) %>% slice_max(n = 3,order_by = LOG10pvalue,with_ties = F)
# 
# go_all_description <- NULL
# for (i in unique(go_top3$ONTOLOGY)) {
#   go_top3_sub <- go_top3[go_top3$ONTOLOGY == i,]
#   ONTOLOGY_sub <- paste0(go_top3_sub$Description,collapse = ", ")
#   if(i == "BP")go_all_description <- c(go_all_description,paste0(ONTOLOGY_sub,'(biological process,BP)'))
#   if(i == "CC")go_all_description <- c(go_all_description,paste0(ONTOLOGY_sub,'(cellular component,CC)'))
#   if(i == "MF")go_all_description <- c(go_all_description,paste0(ONTOLOGY_sub,'(molecular function,MF)'))
# }
# go_all_description <- paste0(go_all_description,collapse = ", ")
# 
# # KEGG 富集结果
# kegg_top3 <- KEGG_list$ekegg_list
# kegg_top3$LOG10pvalue <- -log10(kegg_top3$p.adjust)
# kegg_top3 = kegg_top3 %>% group_by(ONTOLOGY) %>% slice_max(n = 3,order_by = LOG10pvalue,with_ties = F)
# 
# 
# kegg_all_description <- NULL
# for (i in unique(kegg_top3$ONTOLOGY)) {
#   go_top3_sub <- kegg_top3[kegg_top3$ONTOLOGY == i,]
#   ONTOLOGY_sub <- paste0(go_top3_sub$Description,collapse = ", ")
#   kegg_all_description <- c(kegg_all_description,paste0(ONTOLOGY_sub,"(",i,")"))
# }
# kegg_all_description <- paste0(kegg_all_description,collapse = ", ")
# 
# 
# 
# ####材料与方法####
# 
# Methods<-read_docx(cn_path) %>%        
#   body_remove()%>%  
#   body_add_par('材料与方法',style = "标题2") %>% 
#   body_add_fpar(fpar(ftext('GO和KEGG通路富集分析')),style = "标题3") %>% 
#   body_add_fpar(fpar(ftext('Gene Ontology (GO)富集是一种广泛应用于生物信息学研究中的统计方法，它用于检验一个基因集合中某些功能或特征的富集程度。'),
#                      ftext('通过这种分析，可以从大量基因数据中找出具有生物学意义的模式和功能。'),
#                      ftext('GO富集分析主要关注基因在生物学过程（Biological Process, BP）、分子功能（Molecular Function, MF）'),
#                      ftext('和细胞组成（Cellular Component, CC）方面的共同特征[PMID: 25428369]。'),
#                      ftext('Kyoto Encyclopedia of Genes and Genomes（KEGG）是一个综合性的生物信息学数据库，它整合了基因组、化学和系统功能信息，'),
#                      ftext('以理解生物系统如细胞、生物体和生态系统的高级功能和效用，用于挖掘基因列表中富集的显著变化的代谢通路[PMID: 10592173]。'),
#                      ftext('R 包“clusterProfiler”'),
#                      ftext('(version 4.2.2)',prop = underline_font_orange),
#                      ftext('被应用于对'),
#                      ftext(DEG_name,prop=underline_font),
#                      ftext('进行GO和KEGG富集分析（p值<0.05）[PMID: 22455463]。')),
#                 style = "Normal")%>% 
#   print(target = "Report/Methods.docx")
# 
# 
# ####结果####
# Result<- read_docx(cn_path) %>%
#   body_remove() %>%  
#   body_add_par('结果',style = "标题2") %>% 
#   body_add_fpar(fpar(ftext('GO/KEGG富集分析')),style = "标题3") %>% 
#   body_add_fpar(fpar(ftext('为了研究'),
#                      ftext(DEG_name,prop=underline_font),
#                      ftext('的生物学功能，我们进行了GO条目'),
#                      ftext('（Table S）',prop = underline_font_orange),
#                      ftext('和KEGG通路'),
#                      ftext('（Table S）',prop = underline_font_orange),
#                      ftext('的富集分析。'),
#                      ftext('GO结果显示，这些基因在'),
#                      ftext(go_all_description,prop = underline_font_blue),
#                      ftext('中富集'),
#                      ftext('（Figure A）',prop = underline_font_orange),
#                      ftext("。")),
#                 style = "Normal")%>%
#   body_add_fpar(fpar(ftext('富集的KEGG通路包括'),
#                      ftext(kegg_all_description,prop=underline_font_blue),
#                      ftext("（Figure B）",prop = underline_font_orange)),
#                 style = "Normal")%>%
#   body_add_fpar(fpar(ftext('以下为分析的意义，根据自己的基因来源更改对应信息。',prop=underline_font_grey)),
#                 style = "Normal")%>%
#   body_add_fpar(fpar(ftext('这些富集结果提示我们，受'),
#                      ftext(genelist_name_ch,prop=underline_font_red),
#                      ftext('驱动的基因主要通过影响了'),
#                      ftext('细胞因子的调控、免疫系统的激活与信号的转导',prop = underline_font_red),
#                      ftext('过程，参与'),
#                      ftext(case_name_ch,prop = underline_font_red),
#                      ftext('的病理进程。这些分析为解析'),
#                      ftext(paste0(genelist_name_ch,'与',case_name_ch),prop = underline_font_red),
#                      ftext('的联系提供了全新的视角。')),
#                 style = "Normal")%>%
#   print(target = "Report/Result.docx")
# 
# ####图注####
# # source("/opt/scripts_repo/bin/color_to_name.R")
# # color_name1=color_to_name(DEG_list$colors[1])
# # color_name2=color_to_name(DEG_list$colors[2])
# 
# Legend<-read_docx(cn_path) %>%        
#   body_remove()%>%  
#   body_add_par('GO/KEGG图注',style = "标题2") %>% 
#   body_add_fpar(fpar(ftext('柱状图')),style = "标题3") %>% 
#   body_add_fpar(fpar(ftext('柱状图展示'),
#                      ftext(DEG_name,prop=underline_font),
#                      ftext('的'),
#                      ftext('GO/KEGG',prop = underline_font),
#                      ftext('富集分析结果')),
#                 style = "图注")%>% 
#   body_add_fpar(fpar(ftext('气泡图')),style = "标题3") %>% 
#   body_add_fpar(fpar(ftext('气泡图展示'),
#                      ftext(DEG_name,prop=underline_font),
#                      ftext('的'),
#                      ftext('GO/KEGG',prop = underline_font),
#                      ftext('富集分析结果')),
#                 style = "图注")%>% 
#   body_add_fpar(fpar(ftext('弦图')),style = "标题3") %>% 
#   body_add_fpar(fpar(ftext('弦图展示'),
#                      ftext('GO/KEGG',prop = underline_font),
#                      ftext('分析富集的通路。')),
#                 style = "图注")%>% 
#   body_add_fpar(fpar(ftext('圈图和表')),style = "标题3") %>% 
#   body_add_fpar(fpar(ftext('GO/KEGG',prop = underline_font),
#                      ftext('的基因富集分析和通路编号与通路名对应表')),
#                 style = "图注")%>% 
#   print(target = "Report/legend.docx")
# 
# setwd('../')
# 
# }

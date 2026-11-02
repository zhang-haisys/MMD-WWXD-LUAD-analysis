setwd("/JJGene/data/xuwj/麦门冬汤+五味消毒饮_非小细胞肺癌_EGFR-TKI耐药_5+/0.data/")

library(readr)
library(dplyr)
library(stringr)

# 指定主目录路径
main_dir <- "/JJGene/data/xuwj/麦门冬汤+五味消毒饮_非小细胞肺癌_EGFR-TKI耐药_5+/0.data/SwissTargetPrediction+/"

# 获取所有子目录中的.csv文件路径
csv_files <- list.files(main_dir, pattern = "*.csv$", full.names = TRUE, recursive = TRUE)

# 初始化一个空的数据框用于存储合并后的数据
merged_data <- data.frame()

# 循环读取每个.csv文件并合并
for (file in csv_files) {
  # 读取.csv文件
  temp_data <- read_csv(file)
  str <- strsplit(basename(file),"_")[[1]][2]
  temp_data$Source <- str_extract(str, "^[^.]+")
  temp_data$Drug <- strsplit(file,"/")[[1]][9]
  # data1 <- read.csv("./Drug/丹皮/SwissTargetPrediction_Quercetin.csv")
  # str <- strsplit(basename("./Drug/丹皮/SwissTargetPrediction_Quercetin.csv"),"_")[[1]][2]
  # data1$Source <- str_extract(str, "^[^.]+")
  # data1$Drug <- strsplit("./Drug/丹皮/SwissTargetPrediction_Quercetin.csv","/")[[1]][3]
  # 将读取的数据合并到merged_data中
  merged_data <- bind_rows(merged_data, temp_data)
}

# 查看合并后的数据
print(head(merged_data))

# 将合并后的数据保存为一个新的.csv文件
write_csv(merged_data, "./TKI耐药_靶点预测.csv")

PGY <-read.table("./SwissTargetPrediction/蒲公英-Target.txt")
colnames(PGY) <- "x"

Target <-unique(merged_data[merged_data$`Probability*` != 0,]$`Common name`)
write.table(Target, "./药物靶点.txt",quote = FALSE,sep = "\t", row.names = FALSE)

target <- read.table("./药物靶点.txt", header = T)
Target <- rbind(target,PGY)
Target <- unique(Target)




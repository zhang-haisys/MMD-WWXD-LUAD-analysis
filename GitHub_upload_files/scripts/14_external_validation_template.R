# Module 14: external_validation_template
# Source snapshot: archive/original/验证外部数据集..R
# Path adaptation only; scientific settings and historical seeds are retained.
# Run using: Rscript scripts/run_module.R 14
if (!exists("project_dir", mode = "function")) {
  stop("Run this module through scripts/run_module.R from the repository root.")
}

# 一次性安装（如果还没有）
# Dependency installation moved to scripts/install_dependencies.R: install.packages(c("pROC", "ggplot2", "survival", "survminer", "timeROC"))
# BiocManager安装（如果没有则先运行 install.packages("BiocManager")）
# Dependency installation moved to scripts/install_dependencies.R: BiocManager::install(c("limma", "sva"))

# 加载包
library(pROC)
library(ggplot2)
library(survival)
library(survminer)
library(timeROC)
library(limma)
library(sva)

# 假设：
# train_exp: 训练集表达矩阵（行基因，列样本）
# train_meta: 训练集临床信息（含样本ID、生存时间、状态、分组等）
# valid_exp: 验证集表达矩阵（同样格式）
# valid_meta: 验证集临床信息

# 1. 提取模型基因（假设是之前筛选好的基因列表）
model_genes <- c("GENE1", "GENE2", "GENE3")  # 换成你实际的基因名

# 2. 确保验证集包含模型需要的所有基因
missing_genes <- setdiff(model_genes, rownames(valid_exp))
if(length(missing_genes) > 0) {
  cat("缺失基因:", missing_genes, "\n")
  # 填补缺失基因为0（或训练集的均值）
  for(g in missing_genes) {
    valid_exp <- rbind(valid_exp, rep(0, ncol(valid_exp)))
    rownames(valid_exp)[nrow(valid_exp)] <- g
  }
}

# 3. 对齐基因顺序
valid_exp <- valid_exp[model_genes, ]
train_exp <- train_exp[model_genes, ]

# 4. 归一化（用训练集的均值和标准差）
train_mean <- apply(train_exp, 1, mean)
train_sd <- apply(train_exp, 1, sd)
valid_norm <- (valid_exp - train_mean) / train_sd

# 转置：样本为行，基因为列（用于后续计算）
valid_norm <- t(valid_norm)
train_norm <- t(train_exp)

# 假设已有模型系数（从训练集获得的）
model_coef <- c(0.35, -0.21, 0.48)  # 换成你实际的系数，与基因顺序对应

# 计算验证集的风险评分
valid_risk <- valid_norm %*% model_coef
valid_meta$RiskScore <- as.numeric(valid_risk)

# 根据中位数或训练集的cutoff分组
# 方法1：用训练集的中位数作为cutoff
train_cutoff <- median(train_risk)  # 需要先计算train_risk
valid_meta$RiskGroup <- ifelse(valid_meta$RiskScore > train_cutoff, "High", "Low")

# 方法2：用验证集自身中位数
# valid_cutoff <- median(valid_meta$RiskScore)
# valid_meta$RiskGroup <- ifelse(valid_meta$RiskScore > valid_cutoff, "High", "Low")

# 拟合生存对象
fit <- survfit(Surv(OS.time, OS) ~ RiskGroup, data = valid_meta)

# 绘制KM曲线
ggsurvplot(fit, 
           data = valid_meta,
           pval = TRUE,                     # 显示P值
           pval.method = TRUE,              # 显示检验方法
           risk.table = TRUE,               # 显示风险表格
           xlab = "Time (days)",
           ylab = "Overall Survival",
           title = "External Validation: KM Curve",
           legend.title = "Risk Group",
           legend.labs = c("High Risk", "Low Risk"),
           palette = c("#E41A1C", "#377EB8"))

# 假设valid_meta有真实标签列（0=阴性，1=阳性）
# 计算预测概率（如果是逻辑回归模型）
# valid_prob <- plogis(valid_meta$RiskScore)  # 如果RiskScore是线性预测值

# ROC曲线
roc_obj <- roc(valid_meta$TrueLabel, valid_meta$RiskScore)
auc_value <- auc(roc_obj)
print(paste("AUC:", round(auc_value, 3)))

# 绘制ROC曲线
ggroc(roc_obj, color = "#E41A1C", size = 1.2) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.75, y = 0.25, 
           label = paste("AUC =", round(auc_value, 3)), 
           size = 5, hjust = 0) +
  labs(title = "External Validation: ROC Curve",
       x = "1 - Specificity",
       y = "Sensitivity") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

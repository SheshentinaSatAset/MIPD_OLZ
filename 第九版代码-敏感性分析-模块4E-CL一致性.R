## ############################################################################
## 第九版代码-敏感性分析-模块4E-CL一致性.R
## 
## 定位：敏感性分析模块4的补充模块（Part E）
##       在模块4跑完后 source() 本脚本
## 
## 核心任务：
##   评估 MAPB 是否让不同模型对同一患者估计出更一致的 CL
##   - E.1 描述性对比：CV_pop vs CV_ind 分布
##   - E.2 ICC：CL_ind 的跨模型一致性
##   - E.3 配对缩小率 + Wilcoxon 检验
##   - E.4 两两 CCC 矩阵
## 
## 设计说明：
##   - 此脚本由敏感性分析主控脚本调用，不单独运行
##   - 需要预先设置全局变量：SA_TRUE_MODEL, OUT_DIR_SA
##   - 简化输出（无图表），仅保存数值结果
## 
## 输入（从 OUT_DIR_SA 读取）：
##   - mapb_results.rds         （CL_ind 数据来源）
##   - module2_params.rds       （CL_pop 数据来源）
##   - evaluation_results.rds   （追加结果）
## 
## 输出：
##   - evaluation_results.rds   （追加 CL_consistency 结果）
##   - CL_consistency_summary.csv
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("敏感性分析 Part E：CL 跨模型一致性（真实模型 = ", SA_TRUE_MODEL, "）")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## E.0 加载依赖与数据
## ############################################################################

message("\n--- E.0 加载依赖与数据 ---")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
})

## 安装并加载 DescTools（CCC 计算）
if (!requireNamespace("DescTools", quietly = TRUE)) {
  message("📦 正在安装 DescTools 包...")
  install.packages("DescTools", repos = "https://cloud.r-project.org", quiet = TRUE)
}
library(DescTools)

## 读取敏感性分析输出
mapb_results    <- readRDS(file.path(OUT_DIR_SA, "mapb_results.rds"))
module2_params  <- readRDS(file.path(OUT_DIR_SA, "module2_params.rds"))
eval_results    <- readRDS(file.path(OUT_DIR_SA, "evaluation_results.rds"))

N_MC <- nrow(mapb_results)

message("✅ 数据已加载（", N_MC, " 患者）")

## ############################################################################
## E.1 提取 CL_pop 和计算 CL_ind（修复版）
## ############################################################################

message("\n--- E.1 CV 对比：CV_pop vs CV_ind ---")

## CL_pop：从 module2_params 提取
CL_pop_zang <- module2_params$theta_pop$zang$CL
CL_pop_li   <- module2_params$theta_pop$li$CL
CL_pop_yin  <- module2_params$theta_pop$yin$CL
CL_pop_sun  <- module2_params$theta_pop$sun$CL

CL_pop_vec <- c(zang = CL_pop_zang, li = CL_pop_li, 
                yin = CL_pop_yin, sun = CL_pop_sun)

## ✅ 修复：根据 η_CL_ind 计算 CL_ind（mapb_results 中没有直接的 CL_ind 列）
CL_ind_data <- mapb_results %>%
  transmute(
    mc_iter = mc_iter,
    zang = CL_pop_zang * exp(zang_eta_CL_ind),
    li   = CL_pop_li * exp(li_eta_CL_ind),
    yin  = CL_pop_yin * exp(yin_eta_CL_ind),
    sun  = CL_pop_sun * exp(sun_eta_CL_ind)
  )

## 计算 CV_pop（4 个模型 CL_pop 的变异系数）
CV_pop <- sd(CL_pop_vec) / mean(CL_pop_vec) * 100

## 计算每个患者的 CV_ind（4 个模型 CL_ind 的变异系数）
CL_ind_matrix <- as.matrix(CL_ind_data[, c("zang", "li", "yin", "sun")])

CV_ind_vec <- apply(CL_ind_matrix, 1, function(row) {
  sd(row) / mean(row) * 100
})

## 汇总
CV_ind_summary <- c(
  mean   = mean(CV_ind_vec),
  median = median(CV_ind_vec),
  Q1     = quantile(CV_ind_vec, 0.25, names = FALSE),
  Q3     = quantile(CV_ind_vec, 0.75, names = FALSE),
  P5     = quantile(CV_ind_vec, 0.05, names = FALSE),
  P95    = quantile(CV_ind_vec, 0.95, names = FALSE)
)

## CV_ind < CV_pop 的患者比例
pct_CV_reduced <- mean(CV_ind_vec < CV_pop) * 100

cat("\n=== CV 对比 ===\n")
cat(sprintf("  CV_pop (固定值):  %.2f%%\n", CV_pop))
cat(sprintf("  CV_ind (中位数):  %.2f%% [IQR: %.2f%% - %.2f%%]\n",
            CV_ind_summary["median"], CV_ind_summary["Q1"], CV_ind_summary["Q3"]))
cat(sprintf("  CV_ind < CV_pop 的患者比例: %.1f%%\n", pct_CV_reduced))
cat(sprintf("  → MAPB 缩小了跨模型 CL 离散度: %s\n",
            ifelse(CV_ind_summary["median"] < CV_pop, "✅ 是", "❌ 否")))

## ############################################################################
## E.2 ICC：CL_ind 的跨模型一致性
## ############################################################################

message("\n--- E.2 ICC：CL_ind 的跨模型一致性 ---")

## 定义简化的 ICC 计算函数
calc_icc_simple <- function(data_matrix) {
  n_subj <- nrow(data_matrix)
  k_rater <- ncol(data_matrix)
  
  ## 长格式
  long_df <- data.frame(
    subject = rep(1:n_subj, times = k_rater),
    rater   = rep(1:k_rater, each = n_subj),
    value   = as.vector(data_matrix)
  )
  
  ## ANOVA
  fit <- aov(value ~ factor(subject) + factor(rater), data = long_df)
  anova_tbl <- summary(fit)[[1]]
  
  MS_subj <- anova_tbl["factor(subject)", "Mean Sq"]
  MS_rater <- anova_tbl["factor(rater)", "Mean Sq"]
  MS_error <- anova_tbl["Residuals", "Mean Sq"]
  
  ## ICC(2,1)
  icc_2_1 <- (MS_subj - MS_error) / 
    (MS_subj + (k_rater - 1) * MS_error + k_rater * (MS_rater - MS_error) / n_subj)
  
  ## ICC(2,k)
  icc_2_k <- (MS_subj - MS_error) / 
    (MS_subj + (MS_rater - MS_error) / n_subj)
  
  list(
    icc_2_1 = max(0, icc_2_1),
    icc_2_k = max(0, icc_2_k),
    method  = "Two-way random effects, absolute agreement"
  )
}

icc_result <- calc_icc_simple(CL_ind_matrix)
icc_CL_2_1 <- icc_result$icc_2_1
icc_CL_2_k <- icc_result$icc_2_k
icc_CL_method <- icc_result$method

cat("\n=== CL_ind ICC 结果 ===\n")
cat(sprintf("  ICC(2,1) - 单个模型一致性: %.4f\n", icc_CL_2_1))
cat(sprintf("  ICC(2,k) - 平均模型一致性: %.4f\n", icc_CL_2_k))

## ICC 解释
interpret_icc <- function(val) {
  if (val < 0.50) return("Poor")
  if (val < 0.75) return("Moderate")
  if (val < 0.90) return("Good")
  return("Excellent")
}
cat(sprintf("  ICC(2,1) 等级: %s\n", interpret_icc(icc_CL_2_1)))

## ############################################################################
## E.3 配对缩小率 + Wilcoxon 检验
## ############################################################################

message("\n--- E.3 配对缩小率 + Wilcoxon 检验 ---")

## 定义模型对
models <- c("zang", "li", "yin", "sun")
pair_combos <- combn(models, 2, simplify = FALSE)

## 对每一对模型计算缩小率
shrinkage_results <- list()

for (pair in pair_combos) {
  m1 <- pair[1]
  m2 <- pair[2]
  pair_name <- paste(m1, m2, sep = "_vs_")
  
  ## CL_pop 差值（固定常数）
  delta_CL_pop <- abs(CL_pop_vec[m1] - CL_pop_vec[m2])
  
  ## CL_ind 差值（逐患者）
  delta_CL_ind <- abs(CL_ind_data[[m1]] - CL_ind_data[[m2]])
  
  ## 缩小率 = |CL_ind 差值| / |CL_pop 差值|
  if (delta_CL_pop < 1e-10) {
    shrinkage_ratio <- rep(NA_real_, N_MC)
    wilcox_p <- NA_real_
  } else {
    shrinkage_ratio <- delta_CL_ind / delta_CL_pop
    
    ## Wilcoxon 符号秩检验：H0: median(ratio) = 1, H1: median(ratio) < 1
    wilcox_test <- wilcox.test(shrinkage_ratio, mu = 1, alternative = "less",
                               conf.int = TRUE, conf.level = 0.95)
    wilcox_p <- wilcox_test$p.value
  }
  
  shrinkage_results[[pair_name]] <- tibble(
    model_pair         = pair_name,
    delta_CL_pop       = delta_CL_pop,
    delta_CL_ind_mean  = mean(delta_CL_ind),
    delta_CL_ind_median = median(delta_CL_ind),
    ratio_mean         = mean(shrinkage_ratio, na.rm = TRUE),
    ratio_median       = median(shrinkage_ratio, na.rm = TRUE),
    ratio_Q1           = quantile(shrinkage_ratio, 0.25, na.rm = TRUE, names = FALSE),
    ratio_Q3           = quantile(shrinkage_ratio, 0.75, na.rm = TRUE, names = FALSE),
    pct_ratio_lt_1     = mean(shrinkage_ratio < 1, na.rm = TRUE) * 100,
    wilcoxon_p         = wilcox_p,
    significant        = (wilcox_p < 0.05)
  )
}

shrinkage_table <- bind_rows(shrinkage_results)

cat("\n=== 配对缩小率 ===\n")
print(shrinkage_table %>% 
        mutate(across(where(is.numeric), ~ round(.x, 4))))

## 计算统计量
n_pairs_total <- nrow(shrinkage_table)
n_pairs_reduced <- sum(shrinkage_table$ratio_median < 1, na.rm = TRUE)
n_pairs_sig <- sum(shrinkage_table$significant, na.rm = TRUE)

cat(sprintf("\n缩小（ratio < 1）的模型对: %d/%d\n", n_pairs_reduced, n_pairs_total))
cat(sprintf("显著缩小（p < 0.05）的模型对: %d/%d\n", n_pairs_sig, n_pairs_total))

## ############################################################################
## E.4 两两 CCC 矩阵
## ############################################################################

message("\n--- E.4 两两 CCC 矩阵 ---")

n_models <- length(models)

## CCC 矩阵
ccc_matrix <- matrix(NA, n_models, n_models,
                     dimnames = list(models, models))

ccc_table_rows <- list()

for (i in 1:n_models) {
  for (j in 1:n_models) {
    if (i == j) {
      ccc_matrix[i, j] <- 1.0
    } else {
      x <- CL_ind_data[[models[i]]]
      y <- CL_ind_data[[models[j]]]
      
      ## 使用 DescTools::CCC
      ccc_result <- DescTools::CCC(x, y, ci = "z-transform", conf.level = 0.95)
      ccc_val <- ccc_result$rho.c$est
      ccc_matrix[i, j] <- ccc_val
      
      ## 保存详细结果（只对上三角）
      if (i < j) {
        ccc_table_rows[[paste(models[i], models[j], sep = "_")]] <- tibble(
          model_pair = paste(models[i], models[j], sep = "_"),
          CCC        = ccc_val,
          CCC_lower  = ccc_result$rho.c$lwr.ci,
          CCC_upper  = ccc_result$rho.c$upr.ci,
          pearson_r  = cor(x, y),
          Cb         = ccc_val / cor(x, y)  # 偏差修正因子
        )
      }
    }
  }
}

ccc_pairwise <- bind_rows(ccc_table_rows)

cat("\n=== CL_ind 两两 CCC 矩阵 ===\n")
print(round(ccc_matrix, 4))

## CCC 解释
interpret_ccc <- function(val) {
  if (val < 0.90) return("Poor")
  if (val < 0.95) return("Moderate")
  if (val < 0.99) return("Substantial")
  return("Almost perfect")
}

## ############################################################################
## E.5 保存输出
## ############################################################################

message("\n--- E.5 保存输出 ---")

## 构建 Part E 结果
CL_consistency_results <- list(
  
  ## E.1 CV 对比
  CV_comparison = list(
    CV_pop           = CV_pop,
    CV_ind_summary   = CV_ind_summary,
    pct_CV_reduced   = pct_CV_reduced
  ),
  
  ## E.2 ICC
  CL_icc = list(
    icc_2_1 = icc_CL_2_1,
    icc_2_k = icc_CL_2_k,
    method  = icc_CL_method,
    level   = interpret_icc(icc_CL_2_1)
  ),
  
  ## E.3 缩小率
  shrinkage_table = shrinkage_table,
  n_pairs_reduced = n_pairs_reduced,
  n_pairs_sig     = n_pairs_sig,
  
  ## E.4 CCC
  ccc_matrix   = ccc_matrix,
  ccc_pairwise = ccc_pairwise,
  
  ## CL_pop 参考值
  CL_pop = CL_pop_vec,
  
  ## 元数据
  metadata = list(
    created     = Sys.time(),
    module      = "Module 4E - CL Consistency Analysis (Sensitivity)",
    true_model  = SA_TRUE_MODEL
  )
)

## 追加到 evaluation_results
eval_results$CL_consistency <- CL_consistency_results
saveRDS(eval_results, file.path(OUT_DIR_SA, "evaluation_results.rds"))
message("  ✓ evaluation_results.rds 已更新（追加 Part E）")

## 保存 CSV
CL_summary_csv <- tibble(
  metric = c("CV_pop (%)", "CV_ind median (%)", "CV_ind IQR",
             "pct_CV_reduced (%)", "ICC(2,1)", "ICC(2,k)",
             "pairs_reduced", "pairs_significant"),
  value  = c(
    sprintf("%.2f", CV_pop),
    sprintf("%.2f", CV_ind_summary["median"]),
    sprintf("[%.2f, %.2f]", CV_ind_summary["Q1"], CV_ind_summary["Q3"]),
    sprintf("%.1f", pct_CV_reduced),
    sprintf("%.4f", icc_CL_2_1),
    sprintf("%.4f", icc_CL_2_k),
    sprintf("%d/%d", n_pairs_reduced, n_pairs_total),
    sprintf("%d/%d", n_pairs_sig, n_pairs_total)
  )
)

write_csv(CL_summary_csv, file.path(OUT_DIR_SA, "CL_consistency_summary.csv"))
write_csv(ccc_pairwise,   file.path(OUT_DIR_SA, "CL_pairwise_ccc.csv"))
write_csv(shrinkage_table, file.path(OUT_DIR_SA, "CL_shrinkage_results.csv"))
message("  ✓ CSV 文件已保存")

## ############################################################################
## E.6 Part E 汇总报告
## ############################################################################

message("\n", paste(rep("=", 60), collapse = ""))
message("          Part E 汇总：CL 跨模型一致性（", SA_TRUE_MODEL, "）")
message(paste(rep("=", 60), collapse = ""))

cat("\n【E.1 CV 对比】\n")
cat(sprintf("  CV_pop = %.2f%% (固定)\n", CV_pop))
cat(sprintf("  CV_ind = %.2f%% (中位数) [IQR: %.2f%% - %.2f%%]\n",
            CV_ind_summary["median"], CV_ind_summary["Q1"], CV_ind_summary["Q3"]))
cat(sprintf("  %.1f%% 的患者 CV_ind < CV_pop → MAPB 缩小了跨模型 CL 离散度\n",
            pct_CV_reduced))

cat("\n【E.2 ICC】\n")
cat(sprintf("  ICC(2,1) = %.4f (%s)\n", icc_CL_2_1, interpret_icc(icc_CL_2_1)))
cat(sprintf("  ICC(2,k) = %.4f\n", icc_CL_2_k))

cat("\n【E.3 配对缩小率】\n")
for (i in 1:nrow(shrinkage_table)) {
  cat(sprintf("  %s: ΔCL_pop=%.2f → ΔCL_ind=%.2f (median), ratio=%.3f, p=%s\n",
              shrinkage_table$model_pair[i],
              shrinkage_table$delta_CL_pop[i],
              shrinkage_table$delta_CL_ind_median[i],
              shrinkage_table$ratio_median[i],
              format(shrinkage_table$wilcoxon_p[i], digits = 3)))
}

cat("\n【E.4 CCC 矩阵】\n")
cat("  两两 CCC（对角线外平均）:", 
    round(mean(ccc_matrix[upper.tri(ccc_matrix)]), 4), "\n")
for (i in 1:nrow(ccc_pairwise)) {
  cat(sprintf("  %s: CCC=%.4f [%.4f, %.4f] → %s\n",
              ccc_pairwise$model_pair[i],
              ccc_pairwise$CCC[i],
              ccc_pairwise$CCC_lower[i],
              ccc_pairwise$CCC_upper[i],
              interpret_ccc(ccc_pairwise$CCC[i])))
}

cat("\n【结论】\n")
if (pct_CV_reduced > 90 && CV_ind_summary["median"] < CV_pop * 0.6 && n_pairs_reduced >= 5) {
  cat("  ✅ MAPB 有效提升了跨模型 CL 估计的一致性\n")
  cat("     - CL 的跨模型离散度（CV）被显著缩小\n")
  cat("     - 绝大多数模型对的 CL 差值被缩小\n")
} else if (CV_ind_summary["median"] < CV_pop && icc_CL_2_1 > 0.5) {
  cat("  ✅ MAPB 有效提升了跨模型 CL 估计的一致性\n")
  cat("     - CL 的跨模型离散度（CV）被缩小\n")
  cat("     - CL_ind 具有中等以上的跨模型一致性（ICC）\n")
} else {
  cat("  ⚠️  CL 跨模型一致性证据不充分，需进一步分析\n")
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")

message("\n✅ 敏感性分析 Part E 完成（真实模型 = ", SA_TRUE_MODEL, "）")
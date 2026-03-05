## ############################################################################
## 第九版代码-模块4E-CL一致性分析.R
## 
## 定位：模块4的补充模块（Part E），在模块4主脚本跑完后 source() 本脚本
## 
## 核心任务：
##   评估 MAPB 是否让不同模型对同一患者估计出更一致的 CL
##   - E.1 描述性对比：CV_pop vs CV_ind 分布
##   - E.2 ICC：CL_ind 的跨模型一致性
##   - E.3 配对缩小率 + Wilcoxon 检验
##   - E.4 两两 CCC 矩阵 + 热图
## 
## 输入：
##   - outputs9/mapb_results.rds         （CL_ind 数据来源）
##   - outputs9/module2_params.rds       （CL_pop 数据来源）
##   - outputs9/evaluation_results.rds   （追加结果）
## 
## 输出：
##   - outputs9/evaluation_results.rds   （追加 Part E 结果）
##   - outputs9/CL_consistency_summary.csv
##   - outputs9/CL_pairwise_ccc.csv
##   - outputs9/CL_shrinkage_results.csv
##   - outputs9/plots/fig8_CL_ind_scatter_matrix.png
##   - outputs9/plots/fig9_CL_shrinkage_boxplot.png
##   - outputs9/plots/fig10_CL_ccc_heatmap.png
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("Part E：CL 跨模型一致性分析")
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
  library(ggplot2)
})

## 安装并加载 DescTools（CCC 计算）
if (!requireNamespace("DescTools", quietly = TRUE)) {
  message("📦 正在安装 DescTools 包...")
  install.packages("DescTools", repos = "https://cloud.r-project.org", quiet = TRUE)
}
library(DescTools)
message("✅ DescTools 包已加载（CCC 计算可用）")

## 路径设置
PROJ_ROOT <- "D:/Users/YujiaZhang/Desktop/26救赎之道/25.9.29 毕业设计/5正式实施！/多成人模型/体现我卓越的项目管理能力"
OUT_DIR_9 <- file.path(PROJ_ROOT, "outputs9")
plot_dir  <- file.path(OUT_DIR_9, "plots")

## 读取数据
mapb_results    <- readRDS(file.path(OUT_DIR_9, "mapb_results.rds"))
module2_params  <- readRDS(file.path(OUT_DIR_9, "module2_params.rds"))
eval_results    <- readRDS(file.path(OUT_DIR_9, "evaluation_results.rds"))

N_MC <- nrow(mapb_results)

message("✅ 数据已加载")
message("   虚拟患者数: ", N_MC)

## ############################################################################
## E.0.1 提取 CL_pop 和 CL_ind
## ############################################################################

## CL_pop：从 module2_params 中提取（每个模型对全部患者是同一个常数）
CL_pop_zang <- module2_params$theta_pop$zang$CL
CL_pop_li   <- module2_params$theta_pop$li$CL
CL_pop_yin  <- module2_params$theta_pop$yin$CL
CL_pop_sun  <- module2_params$theta_pop$sun$CL

CL_pop_vec <- c(zang = CL_pop_zang, li = CL_pop_li, 
                yin = CL_pop_yin, sun = CL_pop_sun)

cat("\n=== 各模型 CL_pop (L/h) ===\n")
for (nm in names(CL_pop_vec)) {
  cat(sprintf("  %s: %.2f\n", nm, CL_pop_vec[nm]))
}

## CL_ind：从 mapb_results 中提取（每个患者 × 每个模型 一个值）
CL_ind_data <- mapb_results %>%
  select(mc_iter,
         zang = zang_CL_ind,
         li   = li_CL_ind,
         yin  = yin_CL_ind,
         sun  = sun_CL_ind)

cat("\n=== CL_ind 描述统计 (L/h) ===\n")
for (nm in c("zang", "li", "yin", "sun")) {
  vals <- CL_ind_data[[nm]]
  cat(sprintf("  %s: mean=%.2f, sd=%.2f, median=%.2f, [%.2f, %.2f]\n",
              nm, mean(vals), sd(vals), median(vals), 
              quantile(vals, 0.025), quantile(vals, 0.975)))
}

## ############################################################################
## E.1 描述性对比：CV_pop vs CV_ind
## ############################################################################

message("\n--- E.1 描述性对比：CV_pop vs CV_ind ---")

## CV_pop：4个 CL_pop 值的变异系数（一个固定值）
CV_pop <- sd(CL_pop_vec) / mean(CL_pop_vec) * 100

## CV_ind：对每个患者，计算其4个 CL_ind 的 CV
CL_ind_matrix <- as.matrix(CL_ind_data[, -1])
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

## 复用模块4中已定义的 calc_icc_simple 函数
## （如果独立运行本脚本，需要重新定义）
if (!exists("calc_icc_simple")) {
  calc_icc_simple <- function(data_matrix) {
    n <- nrow(data_matrix)
    k <- ncol(data_matrix)
    grand_mean <- mean(as.matrix(data_matrix), na.rm = TRUE)
    row_means  <- rowMeans(data_matrix, na.rm = TRUE)
    SS_between <- k * sum((row_means - grand_mean)^2, na.rm = TRUE)
    SS_within  <- sum((as.matrix(data_matrix) - outer(row_means, rep(1, k)))^2, na.rm = TRUE)
    BMS <- SS_between / (n - 1)
    WMS <- SS_within  / (n * (k - 1))
    icc_2_1 <- (BMS - WMS) / (BMS + (k - 1) * WMS)
    icc_2_k <- (BMS - WMS) / BMS
    list(icc_2_1 = max(0, icc_2_1), icc_2_k = max(0, icc_2_k))
  }
}

## 计算 CL_ind 的 ICC
HAS_PSYCH <- requireNamespace("psych", quietly = TRUE)
HAS_LME4  <- requireNamespace("lme4",  quietly = TRUE)

if (HAS_PSYCH && HAS_LME4) {
  icc_CL_result <- psych::ICC(CL_ind_data[, -1])
  icc_CL_2_1 <- icc_CL_result$results$ICC[2]
  icc_CL_2_k <- icc_CL_result$results$ICC[5]
  icc_CL_method <- "psych"
} else {
  icc_CL_simple <- calc_icc_simple(CL_ind_data[, -1])
  icc_CL_2_1 <- icc_CL_simple$icc_2_1
  icc_CL_2_k <- icc_CL_simple$icc_2_k
  icc_CL_method <- "simplified"
}

cat("\n=== CL_ind ICC 结果 ===\n")
cat(sprintf("  ICC(2,1) - 单个模型一致性: %.4f\n", icc_CL_2_1))
cat(sprintf("  ICC(2,k) - 平均模型一致性: %.4f\n", icc_CL_2_k))
cat(sprintf("  计算方法: %s\n", icc_CL_method))

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
  ## 注意：如果 delta_CL_pop = 0，跳过（理论上不会发生）
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

## 汇总判断
n_pairs_reduced <- sum(shrinkage_table$ratio_median < 1, na.rm = TRUE)
n_pairs_total   <- nrow(shrinkage_table)
n_pairs_sig     <- sum(shrinkage_table$significant, na.rm = TRUE)

cat(sprintf("\n  %d/%d 对模型的 CL 差异被 MAPB 缩小（median ratio < 1）\n",
            n_pairs_reduced, n_pairs_total))
cat(sprintf("  %d/%d 对达到统计显著（Wilcoxon p < 0.05）\n",
            n_pairs_sig, n_pairs_total))

## 保存逐患者的缩小率数据（供绘图用）
shrinkage_raw <- list()
for (pair in pair_combos) {
  m1 <- pair[1]
  m2 <- pair[2]
  pair_name <- paste(m1, m2, sep = "_vs_")
  delta_CL_pop <- abs(CL_pop_vec[m1] - CL_pop_vec[m2])
  delta_CL_ind <- abs(CL_ind_data[[m1]] - CL_ind_data[[m2]])
  
  if (delta_CL_pop < 1e-10) {
    shrinkage_raw[[pair_name]] <- rep(NA_real_, N_MC)
  } else {
    shrinkage_raw[[pair_name]] <- delta_CL_ind / delta_CL_pop
  }
}

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

cat("\n=== CL_ind 两两 CCC 详细 ===\n")
print(ccc_pairwise %>% mutate(across(where(is.numeric), ~ round(.x, 4))))

## CCC 解释
interpret_ccc <- function(val) {
  if (val < 0.90) return("Poor")
  if (val < 0.95) return("Moderate")
  if (val < 0.99) return("Substantial")
  return("Almost perfect")
}

cat("\n=== CCC 等级判定 ===\n")
for (i in 1:nrow(ccc_pairwise)) {
  cat(sprintf("  %s: CCC=%.4f → %s\n",
              ccc_pairwise$model_pair[i],
              ccc_pairwise$CCC[i],
              interpret_ccc(ccc_pairwise$CCC[i])))
}

## ############################################################################
## E.5 可视化
## ############################################################################

message("\n", paste(rep("-", 70), collapse = ""))
message("E.5 生成 Part E 可视化")
message(paste(rep("-", 70), collapse = ""))

## ----------------------------------------------------------------------------
## 图 8：CL_ind 跨模型散点图矩阵
## ----------------------------------------------------------------------------

message("\n--- 生成图 8：CL_ind 跨模型散点图矩阵 ---")

CL_scatter_plots <- list()

for (pair in pair_combos) {
  m1 <- pair[1]
  m2 <- pair[2]
  pair_name <- paste(m1, m2, sep = "_")
  
  ccc_val <- ccc_matrix[m1, m2]
  
  p <- ggplot(CL_ind_data, aes(x = .data[[m1]], y = .data[[m2]])) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
    geom_point(alpha = 0.2, size = 0.8, color = "grey30") +
    geom_smooth(method = "lm", se = FALSE, color = "steelblue", linewidth = 0.8) +
    ## 标注 CL_pop 的位置
    geom_point(data = data.frame(x = CL_pop_vec[m1], y = CL_pop_vec[m2]),
               aes(x = x, y = y), 
               color = "red", size = 3, shape = 18) +
    labs(
      title = paste(m1, "vs", m2),
      subtitle = sprintf("CCC = %.3f", ccc_val),
      x = paste0("CL_ind (", m1, ") [L/h]"),
      y = paste0("CL_ind (", m2, ") [L/h]")
    ) +
    coord_equal() +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      plot.subtitle = element_text(hjust = 0.5, color = "steelblue", size = 9)
    )
  
  CL_scatter_plots[[pair_name]] <- p
}

## 组合图
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p8_CL_scatter <- wrap_plots(CL_scatter_plots, ncol = 3) +
    plot_annotation(
      title = "CL_ind 跨模型散点图（红色菱形 = CL_pop 位置）",
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
    )
  
  ggsave(file.path(plot_dir, "fig8_CL_ind_scatter_matrix.png"),
         p8_CL_scatter, width = 12, height = 8, dpi = 300)
  message("  ✓ fig8_CL_ind_scatter_matrix.png 已保存")
} else {
  message("  ⚠️ patchwork 包未安装，将单独保存散点图")
  for (name in names(CL_scatter_plots)) {
    ggsave(file.path(plot_dir, paste0("fig8_scatter_CL_", name, ".png")),
           CL_scatter_plots[[name]], width = 5, height = 5, dpi = 300)
  }
  message("  ✓ CL 散点图已单独保存")
}

## ----------------------------------------------------------------------------
## 图 9：缩小率箱线图
## ----------------------------------------------------------------------------

message("\n--- 生成图 9：缩小率箱线图 ---")

## 准备数据
shrinkage_plot_data <- bind_rows(
  lapply(names(shrinkage_raw), function(pname) {
    tibble(
      model_pair = pname,
      ratio = shrinkage_raw[[pname]]
    )
  })
) %>%
  filter(!is.na(ratio))

## 添加中位数标签数据
median_labels <- shrinkage_plot_data %>%
  group_by(model_pair) %>%
  summarise(
    median_ratio = median(ratio),
    pct_lt_1 = sprintf("%.0f%%<1", mean(ratio < 1) * 100),
    .groups = "drop"
  )

p9_shrinkage <- ggplot(shrinkage_plot_data, aes(x = model_pair, y = ratio)) +
  geom_hline(yintercept = 1, linetype = 2, color = "red", linewidth = 0.8) +
  geom_boxplot(fill = "steelblue", alpha = 0.6, outlier.size = 0.5, outlier.alpha = 0.3) +
  geom_text(data = median_labels, 
            aes(x = model_pair, y = max(shrinkage_plot_data$ratio, na.rm = TRUE) * 0.02,
                label = pct_lt_1),
            vjust = 0, size = 3.5, color = "darkblue", fontface = "bold") +
  scale_y_continuous(
    limits = c(0, quantile(shrinkage_plot_data$ratio, 0.99, na.rm = TRUE) * 1.1),
    breaks = c(0, 0.5, 1, 1.5, 2)
  ) +
  labs(
    title = "CL 跨模型差异缩小率",
    subtitle = "ratio = |CL_ind_A - CL_ind_B| / |CL_pop_A - CL_pop_B|    红线 = 1（无改善）",
    x = "模型对",
    y = "缩小率 (ratio)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(plot_dir, "fig9_CL_shrinkage_boxplot.png"),
       p9_shrinkage, width = 9, height = 6, dpi = 300)

message("  ✓ fig9_CL_shrinkage_boxplot.png 已保存")

## ----------------------------------------------------------------------------
## 图 10：两两 CCC 热图
## ----------------------------------------------------------------------------

message("\n--- 生成图 10：CCC 热图 ---")

ccc_heatmap_data <- as.data.frame(ccc_matrix) %>%
  mutate(model_row = rownames(.)) %>%
  pivot_longer(cols = -model_row, names_to = "model_col", values_to = "CCC") %>%
  mutate(
    model_row = factor(model_row, levels = models),
    model_col = factor(model_col, levels = models)
  )

p10_ccc <- ggplot(ccc_heatmap_data, aes(x = model_col, y = model_row, fill = CCC)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.3f", CCC)), size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "white", mid = "#4DBBD5", high = "#00468B",
    midpoint = 0.5, limits = c(0, 1),
    name = "CCC"
  ) +
  labs(
    title = "CL_ind 跨模型 CCC 热图",
    subtitle = "Lin's Concordance Correlation Coefficient — 同时衡量相关性与偏移",
    x = "模型",
    y = "模型"
  ) +
  coord_equal() +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(file.path(plot_dir, "fig10_CL_ccc_heatmap.png"),
       p10_ccc, width = 6.5, height = 5.5, dpi = 300)

message("  ✓ fig10_CL_ccc_heatmap.png 已保存")

## ############################################################################
## E.6 保存输出
## ############################################################################

message("\n--- E.6 保存输出 ---")

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
    created = Sys.time(),
    module  = "Module 4E - CL Consistency Analysis"
  )
)

## 追加到 evaluation_results
eval_results$CL_consistency <- CL_consistency_results
saveRDS(eval_results, file.path(OUT_DIR_9, "evaluation_results.rds"))
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

write_csv(CL_summary_csv, file.path(OUT_DIR_9, "CL_consistency_summary.csv"))
write_csv(ccc_pairwise,   file.path(OUT_DIR_9, "CL_pairwise_ccc.csv"))
write_csv(shrinkage_table, file.path(OUT_DIR_9, "CL_shrinkage_results.csv"))
message("  ✓ CSV 文件已保存")

## ############################################################################
## E.7 Part E 汇总报告
## ############################################################################

message("\n", paste(rep("=", 60), collapse = ""))
message("          Part E 汇总：CL 跨模型一致性")
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
if (CV_ind_summary["median"] < CV_pop && icc_CL_2_1 > 0.5) {
  cat("  ✅ MAPB 有效提升了跨模型 CL 估计的一致性\n")
  cat("     - CL 的跨模型离散度（CV）被缩小\n")
  cat("     - CL_ind 具有中等以上的跨模型一致性（ICC）\n")
} else {
  cat("  ⚠️  CL 跨模型一致性证据不充分，需进一步分析\n")
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")

message("\n", paste(rep("=", 70), collapse = ""))
message("✅ Part E 全部完成（CL 跨模型一致性分析）")
message(paste(rep("=", 70), collapse = ""))

message("\n📄 新增输出文件：")
message("  数据：")
message("  - evaluation_results.rds        （已追加 CL_consistency）")
message("  - CL_consistency_summary.csv")
message("  - CL_pairwise_ccc.csv")
message("  - CL_shrinkage_results.csv")
message("\n  图表（plots/）：")
message("  - fig8_CL_ind_scatter_matrix.png")
message("  - fig9_CL_shrinkage_boxplot.png")
message("  - fig10_CL_ccc_heatmap.png")

message("\n🔗 使用方法：")
message('  results <- readRDS("outputs9/evaluation_results.rds")')
message('  CL_cons <- results$CL_consistency')
message('  CL_cons$CV_comparison     # CV 对比')
message('  CL_cons$CL_icc            # ICC')
message('  CL_cons$shrinkage_table   # 缩小率')
message('  CL_cons$ccc_pairwise      # CCC')
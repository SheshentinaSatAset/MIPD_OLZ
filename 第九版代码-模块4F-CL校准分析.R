## ############################################################################
## 第九版代码-模块4F-CL校准分析.R
## 
## 定位：模块4的补充模块（Part F），在模块4E跑完后 source() 本脚本
## 
## 核心任务：
##   以 zang2021 生成的 CL_true 为真值，评估 4 个模型的 CL 估计精度
##   - F.1 rBias / rRMSE / F30（IND vs POP）
##   - F.2 Bland-Altman 图（4 个模型各一张）
## 
## 输入：
##   - outputs9/data_true.rds            （CL_true 来源）
##   - outputs9/mapb_results.rds         （CL_ind 来源）
##   - outputs9/module2_params.rds       （CL_pop 来源）
##   - outputs9/evaluation_results.rds   （追加结果）
## 
## 输出：
##   - outputs9/evaluation_results.rds   （追加 CL_calibration 结果）
##   - outputs9/CL_calibration_metrics.csv
##   - outputs9/plots/fig11_BA_zang.png
##   - outputs9/plots/fig11_BA_li.png
##   - outputs9/plots/fig11_BA_yin.png
##   - outputs9/plots/fig11_BA_sun.png
##   - outputs9/plots/fig12_CL_rBias_rRMSE_barplot.png
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("Part F：CL 校准分析（rBias / rRMSE / F30 + Bland-Altman）")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## F.0 加载依赖与数据
## ############################################################################

message("\n--- F.0 加载依赖与数据 ---")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

## 路径设置
PROJ_ROOT <- "D:/Users/YujiaZhang/Desktop/26救赎之道/25.9.29 毕业设计/5正式实施！/多成人模型/体现我卓越的项目管理能力"
OUT_DIR_9 <- file.path(PROJ_ROOT, "outputs9")
plot_dir  <- file.path(OUT_DIR_9, "plots")

## 读取数据
data_true       <- readRDS(file.path(OUT_DIR_9, "data_true.rds"))
mapb_results    <- readRDS(file.path(OUT_DIR_9, "mapb_results.rds"))
module2_params  <- readRDS(file.path(OUT_DIR_9, "module2_params.rds"))
eval_results    <- readRDS(file.path(OUT_DIR_9, "evaluation_results.rds"))

N_MC <- nrow(mapb_results)

message("✅ 数据已加载")
message("   虚拟患者数: ", N_MC)

## ############################################################################
## F.0.1 提取 CL_true, CL_ind, CL_pop
## ############################################################################

## CL_true：来自 data_true（zang2021 生成）
CL_true <- data_true$CL_true

## CL_pop：各模型的群体典型值（常数）
CL_pop_zang <- module2_params$theta_pop$zang$CL
CL_pop_li   <- module2_params$theta_pop$li$CL
CL_pop_yin  <- module2_params$theta_pop$yin$CL
CL_pop_sun  <- module2_params$theta_pop$sun$CL

## CL_ind：来自 mapb_results
CL_ind_zang <- mapb_results$zang_CL_ind
CL_ind_li   <- mapb_results$li_CL_ind
CL_ind_yin  <- mapb_results$yin_CL_ind
CL_ind_sun  <- mapb_results$sun_CL_ind

## 构建分析数据框
calibration_data <- tibble(
  mc_iter = 1:N_MC,
  CL_true = CL_true,
  
  ## IND 模式
  zang_CL_ind = CL_ind_zang,
  li_CL_ind   = CL_ind_li,
  yin_CL_ind  = CL_ind_yin,
  sun_CL_ind  = CL_ind_sun,
  
  ## POP 模式（常数，重复 N_MC 次）
  zang_CL_pop = CL_pop_zang,
  li_CL_pop   = CL_pop_li,
  yin_CL_pop  = CL_pop_yin,
  sun_CL_pop  = CL_pop_sun
)

cat("\n=== CL 数据摘要 ===\n")
cat("CL_true: mean =", round(mean(CL_true), 2), "L/h, sd =", round(sd(CL_true), 2), "\n")
cat("CL_pop (zang/li/yin/sun):", 
    round(CL_pop_zang, 2), "/", round(CL_pop_li, 2), "/",
    round(CL_pop_yin, 2), "/", round(CL_pop_sun, 2), "L/h\n")

## ############################################################################
## F.1 计算 rBias / rRMSE / F30
## ############################################################################

message("\n--- F.1 计算 rBias / rRMSE / F30 ---")

## 定义计算函数
calc_calibration_metrics <- function(CL_pred, CL_true) {
  
  ## 相对误差
  RE <- (CL_pred - CL_true) / CL_true
  
  ## rBias (%) = mean(RE) × 100
  rBias <- mean(RE) * 100
  
  ## rRMSE (%) = sqrt(mean(RE^2)) × 100
  rRMSE <- sqrt(mean(RE^2)) * 100
  
  ## F30 (%) = 患者中 |RE| ≤ 0.30 的比例 × 100
  F30 <- mean(abs(RE) <= 0.30) * 100
  
  ## F20 (%)
  F20 <- mean(abs(RE) <= 0.20) * 100
  
  ## 其他统计
  list(
    rBias_pct  = rBias,
    rRMSE_pct  = rRMSE,
    F20_pct    = F20,
    F30_pct    = F30,
    RE_mean    = mean(RE),
    RE_median  = median(RE),
    RE_sd      = sd(RE),
    RE_Q1      = quantile(RE, 0.25, names = FALSE),
    RE_Q3      = quantile(RE, 0.75, names = FALSE)
  )
}

## 计算各模型的指标
models <- c("zang", "li", "yin", "sun")

metrics_list <- list()

for (m in models) {
  
  CL_ind_col <- paste0(m, "_CL_ind")
  CL_pop_col <- paste0(m, "_CL_pop")
  
  ## IND 模式
  metrics_ind <- calc_calibration_metrics(
    CL_pred = calibration_data[[CL_ind_col]],
    CL_true = calibration_data$CL_true
  )
  
  ## POP 模式
  metrics_pop <- calc_calibration_metrics(
    CL_pred = calibration_data[[CL_pop_col]],
    CL_true = calibration_data$CL_true
  )
  
  metrics_list[[m]] <- tibble(
    model = m,
    mode  = c("IND", "POP"),
    rBias_pct  = c(metrics_ind$rBias_pct, metrics_pop$rBias_pct),
    rRMSE_pct  = c(metrics_ind$rRMSE_pct, metrics_pop$rRMSE_pct),
    F20_pct    = c(metrics_ind$F20_pct, metrics_pop$F20_pct),
    F30_pct    = c(metrics_ind$F30_pct, metrics_pop$F30_pct),
    RE_mean    = c(metrics_ind$RE_mean, metrics_pop$RE_mean),
    RE_sd      = c(metrics_ind$RE_sd, metrics_pop$RE_sd)
  )
}

calibration_metrics <- bind_rows(metrics_list)

## 计算 IND vs POP 的差值
calibration_summary <- calibration_metrics %>%
  pivot_wider(
    names_from = mode,
    values_from = c(rBias_pct, rRMSE_pct, F20_pct, F30_pct, RE_mean, RE_sd)
  ) %>%
  mutate(
    delta_rBias  = abs(rBias_pct_IND) - abs(rBias_pct_POP),
    delta_rRMSE  = rRMSE_pct_IND - rRMSE_pct_POP,
    delta_F30    = F30_pct_IND - F30_pct_POP
  )

## 打印结果
cat("\n=== CL 校准指标（IND vs POP，以 zang2021 CL_true 为真值）===\n\n")

print(calibration_metrics %>%
        mutate(across(where(is.numeric), ~ round(.x, 2))))

cat("\n=== IND vs POP 差值 ===\n")
cat("（负的 delta_rBias/delta_rRMSE 表示 IND 更好，正的 delta_F30 表示 IND 更好）\n\n")

print(calibration_summary %>%
        select(model, delta_rBias, delta_rRMSE, delta_F30) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2))))

## ############################################################################
## F.2 Bland-Altman 图（统一坐标版）
## ############################################################################

message("\n--- F.2 Bland-Altman 图（统一坐标）---")

## ----------------------------------------------------------------------------
## 先计算所有模型的数据范围，确定统一坐标轴
## ----------------------------------------------------------------------------

models <- c("zang", "li", "yin", "sun")

## 收集所有模型的 mean 和 diff 值
all_BA_data <- list()

for (m in models) {
  CL_ind_col <- paste0(m, "_CL_ind")
  
  mean_val <- (calibration_data[[CL_ind_col]] + calibration_data$CL_true) / 2
  diff_val <- calibration_data[[CL_ind_col]] - calibration_data$CL_true
  
  all_BA_data[[m]] <- tibble(
    model = m,
    mean = mean_val,
    diff = diff_val
  )
}

all_BA_combined <- bind_rows(all_BA_data)

## 计算统一的坐标范围
x_range <- range(all_BA_combined$mean, na.rm = TRUE)
y_range <- range(all_BA_combined$diff, na.rm = TRUE)

## 稍微扩展范围（留出边距）
x_limits <- c(x_range[1] * 0.95, x_range[2] * 1.05)
y_limits <- c(y_range[1] * 1.1, y_range[2] * 1.1)

## 确保 y 轴对称（以 0 为中心更直观）
y_max_abs <- max(abs(y_limits))
y_limits_sym <- c(-y_max_abs, y_max_abs)

message("  统一坐标范围:")
message("    X 轴 (Mean): [", round(x_limits[1], 1), ", ", round(x_limits[2], 1), "] L/h")
message("    Y 轴 (Diff): [", round(y_limits_sym[1], 1), ", ", round(y_limits_sym[2], 1), "] L/h")

## ----------------------------------------------------------------------------
## 定义 Bland-Altman 绘图函数（统一坐标版）
## ----------------------------------------------------------------------------

plot_bland_altman_unified <- function(CL_pred, CL_true, model_name, 
                                       x_limits, y_limits, mode_name = "IND") {
  
  ## 计算
  mean_val <- (CL_pred + CL_true) / 2
  diff_val <- CL_pred - CL_true
  
  ## 统计量
  mean_diff <- mean(diff_val)
  sd_diff   <- sd(diff_val)
  upper_loa <- mean_diff + 1.96 * sd_diff
  lower_loa <- mean_diff - 1.96 * sd_diff
  
  ## 计算相对指标（用于注释）
  rBias <- mean(diff_val / CL_true) * 100
  rRMSE <- sqrt(mean((diff_val / CL_true)^2)) * 100
  
  ## 构建数据框
  plot_df <- tibble(
    mean = mean_val,
    diff = diff_val
  )
  
  ## 绘图
  p <- ggplot(plot_df, aes(x = mean, y = diff)) +
    ## 参考线
    geom_hline(yintercept = 0, linetype = 1, color = "grey50", linewidth = 0.5) +
    geom_hline(yintercept = mean_diff, linetype = 2, color = "blue", linewidth = 0.8) +
    geom_hline(yintercept = upper_loa, linetype = 3, color = "red", linewidth = 0.7) +
    geom_hline(yintercept = lower_loa, linetype = 3, color = "red", linewidth = 0.7) +
    ## 散点
    geom_point(alpha = 0.25, size = 1.2, color = "grey30") +
    ## 平滑线（检测趋势）
    geom_smooth(method = "loess", se = FALSE, color = "steelblue", 
                linewidth = 0.8, span = 0.5) +
    ## 注释（放在固定位置，基于统一坐标）
    annotate("text", x = x_limits[2] * 0.95, y = mean_diff, 
             label = sprintf("Mean = %.2f", mean_diff),
             hjust = 1, vjust = -0.5, color = "blue", size = 3.5) +
    annotate("text", x = x_limits[2] * 0.95, y = upper_loa,
             label = sprintf("+1.96SD = %.2f", upper_loa),
             hjust = 1, vjust = -0.5, color = "red", size = 3) +
    annotate("text", x = x_limits[2] * 0.95, y = lower_loa,
             label = sprintf("-1.96SD = %.2f", lower_loa),
             hjust = 1, vjust = 1.5, color = "red", size = 3) +
    ## 统一坐标范围
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    ## 标签
    labs(
      title = model_name,
      subtitle = sprintf("rBias = %.1f%%, rRMSE = %.1f%%, LoA = [%.1f, %.1f]",
                         rBias, rRMSE, lower_loa, upper_loa),
      x = "Mean of CL_ind and CL_true (L/h)",
      y = "CL_ind - CL_true (L/h)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9)
    )
  
  ## 返回绘图对象和统计量
  list(
    plot = p,
    stats = list(
      mean_diff = mean_diff,
      sd_diff = sd_diff,
      upper_loa = upper_loa,
      lower_loa = lower_loa,
      rBias = rBias,
      rRMSE = rRMSE
    )
  )
}

## ----------------------------------------------------------------------------
## 为每个模型生成 Bland-Altman 图（统一坐标）
## ----------------------------------------------------------------------------

BA_results <- list()

for (m in models) {
  
  CL_ind_col <- paste0(m, "_CL_ind")
  
  ## 模型标签
  model_label <- switch(m,
    zang = "zang2021",
    li   = "li2018",
    yin  = "yin2016",
    sun  = "sun2021"
  )
  
  ## 生成图
  ba_result <- plot_bland_altman_unified(
    CL_pred = calibration_data[[CL_ind_col]],
    CL_true = calibration_data$CL_true,
    model_name = model_label,
    x_limits = x_limits,
    y_limits = y_limits_sym,
    mode_name = "IND"
  )
  
  BA_results[[m]] <- ba_result
  
  ## 保存单独图
  ggsave(
    filename = file.path(plot_dir, paste0("fig11_BA_", m, ".png")),
    plot = ba_result$plot,
    width = 7, height = 5.5, dpi = 300
  )
  
  message("  ✓ fig11_BA_", m, ".png 已保存")
}

## ----------------------------------------------------------------------------
## 组合 4 张图为一张大图（2×2 布局）
## ----------------------------------------------------------------------------

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  
  p_BA_combined <- (BA_results$zang$plot + BA_results$li$plot) /
                   (BA_results$yin$plot + BA_results$sun$plot) +
    plot_annotation(
      title = "Bland-Altman: CL_ind vs CL_true（统一坐标）",
      subtitle = "横轴 = (CL_ind + CL_true)/2，纵轴 = CL_ind - CL_true",
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11)
      )
    )
  
  ggsave(
    filename = file.path(plot_dir, "fig11_BA_combined_4panels.png"),
    plot = p_BA_combined,
    width = 12, height = 10, dpi = 300
  )
  
  message("  ✓ fig11_BA_combined_4panels.png 已保存（2×2 组合图）")
  
} else {
  message("  ⚠️ patchwork 包未安装，跳过组合图")
}

## ############################################################################
## F.3 rBias / rRMSE 柱状图
## ############################################################################

message("\n--- F.3 生成 rBias / rRMSE 对比柱状图 ---")

## 准备数据
barplot_data <- calibration_metrics %>%
  select(model, mode, rBias_pct, rRMSE_pct, F30_pct) %>%
  mutate(
    model = factor(model, levels = c("zang", "li", "yin", "sun"),
                   labels = c("zang2021", "li2018", "yin2016", "sun2021")),
    mode = factor(mode, levels = c("POP", "IND"))
  )

## rBias 柱状图（显示绝对值，但标注正负）
p_rBias <- ggplot(calibration_metrics %>%
                    mutate(model = factor(model, levels = c("zang", "li", "yin", "sun"),
                                          labels = c("zang2021", "li2018", "yin2016", "sun2021")),
                           mode = factor(mode, levels = c("POP", "IND"))),
                  aes(x = model, y = rBias_pct, fill = mode)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = 1, color = "grey30") +
  geom_text(aes(label = sprintf("%.1f%%", rBias_pct)),
            position = position_dodge(width = 0.8),
            vjust = ifelse(calibration_metrics$rBias_pct >= 0, -0.5, 1.5),
            size = 3) +
  scale_fill_manual(values = c("POP" = "coral", "IND" = "steelblue")) +
  labs(
    title = "CL rBias（以 zang2021 CL_true 为真值）",
    subtitle = "rBias = mean((CL_pred - CL_true) / CL_true) × 100%",
    x = "模型",
    y = "rBias (%)",
    fill = "模式"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
    legend.position = "bottom"
  )

## rRMSE 柱状图
p_rRMSE <- ggplot(barplot_data,
                  aes(x = model, y = rRMSE_pct, fill = mode)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", rRMSE_pct)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("POP" = "coral", "IND" = "steelblue")) +
  scale_y_continuous(limits = c(0, max(barplot_data$rRMSE_pct) * 1.15)) +
  labs(
    title = "CL rRMSE（�� zang2021 CL_true 为真值）",
    subtitle = "rRMSE = sqrt(mean(((CL_pred - CL_true) / CL_true)²)) × 100%",
    x = "模型",
    y = "rRMSE (%)",
    fill = "模式"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
    legend.position = "bottom"
  )

## F30 柱状图
p_F30 <- ggplot(barplot_data,
                aes(x = model, y = F30_pct, fill = mode)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", F30_pct)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("POP" = "coral", "IND" = "steelblue")) +
  scale_y_continuous(limits = c(0, 105)) +
  labs(
    title = "CL F30（以 zang2021 CL_true 为真值）",
    subtitle = "F30 = |RE| ≤ 30% 的患者比例",
    x = "模型",
    y = "F30 (%)",
    fill = "模式"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
    legend.position = "bottom"
  )

## 组合图
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p_combined <- p_rBias / p_rRMSE / p_F30 +
    plot_annotation(
      title = "CL 校准指标对比（IND vs POP）",
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
    )
  
  ggsave(file.path(plot_dir, "fig12_CL_calibration_barplot.png"),
         p_combined, width = 9, height = 12, dpi = 300)
  message("  ✓ fig12_CL_calibration_barplot.png 已保存")
  
} else {
  ## 单独保存
  ggsave(file.path(plot_dir, "fig12a_CL_rBias_barplot.png"),
         p_rBias, width = 8, height = 5, dpi = 300)
  ggsave(file.path(plot_dir, "fig12b_CL_rRMSE_barplot.png"),
         p_rRMSE, width = 8, height = 5, dpi = 300)
  ggsave(file.path(plot_dir, "fig12c_CL_F30_barplot.png"),
         p_F30, width = 8, height = 5, dpi = 300)
  message("  ✓ fig12a/b/c 已单独保存")
}

## ############################################################################
## F.4 保存输出
## ############################################################################

message("\n--- F.4 保存输出 ---")

## 构建 Part F 结果
CL_calibration_results <- list(
  
  ## F.1 校准指标
  calibration_metrics = calibration_metrics,
  calibration_summary = calibration_summary,
  
  ## F.2 Bland-Altman 统计量
  BA_stats = lapply(BA_results, function(x) x$stats),
  
  ## 参考值
  CL_pop = c(zang = CL_pop_zang, li = CL_pop_li, 
             yin = CL_pop_yin, sun = CL_pop_sun),
  
  ## 元数据
  metadata = list(
    created = Sys.time(),
    module  = "Module 4F - CL Calibration Analysis",
    note    = "CL_true 来自 zang2021 生成，其他模型为外部验证"
  )
)

## 追加到 evaluation_results
eval_results$CL_calibration <- CL_calibration_results
saveRDS(eval_results, file.path(OUT_DIR_9, "evaluation_results.rds"))
message("  ✓ evaluation_results.rds 已更新（追加 Part F）")

## 保存 CSV
write_csv(calibration_metrics, file.path(OUT_DIR_9, "CL_calibration_metrics.csv"))
message("  ✓ CL_calibration_metrics.csv 已保存")

## ############################################################################
## F.5 Part F 汇总报告
## ############################################################################

message("\n", paste(rep("=", 60), collapse = ""))
message("          Part F 汇总：CL 校准分析")
message(paste(rep("=", 60), collapse = ""))

cat("\n【解读说明】\n")
cat("  - CL_true 来自 zang2021 模型生成（真实模型）\n")
cat("  - zang2021 的 IND 结果反映 MAPB 的校准质量\n")
cat("  - 其他模型的结果反映外部模型对真值的逼近程度\n")

cat("\n【F.1 校准指标】\n\n")

for (m in models) {
  
  ind_row <- calibration_metrics %>% filter(model == m, mode == "IND")
  pop_row <- calibration_metrics %>% filter(model == m, mode == "POP")
  
  cat(sprintf("  %s:\n", m))
  cat(sprintf("    IND: rBias=%+.1f%%, rRMSE=%.1f%%, F30=%.1f%%\n",
              ind_row$rBias_pct, ind_row$rRMSE_pct, ind_row$F30_pct))
  cat(sprintf("    POP: rBias=%+.1f%%, rRMSE=%.1f%%, F30=%.1f%%\n",
              pop_row$rBias_pct, pop_row$rRMSE_pct, pop_row$F30_pct))
  cat(sprintf("    → IND vs POP: ΔrRMSE=%+.1f%%, ΔF30=%+.1f%%\n\n",
              ind_row$rRMSE_pct - pop_row$rRMSE_pct,
              ind_row$F30_pct - pop_row$F30_pct))
}

cat("\n【F.2 Bland-Altman 统计量（IND 模式）】\n\n")

for (m in models) {
  ba_stats <- BA_results[[m]]$stats
  cat(sprintf("  %s: Mean Bias = %.2f L/h, LoA = [%.2f, %.2f] L/h\n",
              m, ba_stats$mean_diff, ba_stats$lower_loa, ba_stats$upper_loa))
}

cat("\n【结论】\n")

## 判断 zang2021（内部验证）
zang_ind <- calibration_metrics %>% filter(model == "zang", mode == "IND")
zang_pop <- calibration_metrics %>% filter(model == "zang", mode == "POP")

if (zang_ind$rRMSE_pct < zang_pop$rRMSE_pct && zang_ind$F30_pct > zang_pop$F30_pct) {
  cat("  ✅ zang2021（真实模型）：MAPB 有效改善了 CL 估计精度\n")
  cat(sprintf("     rRMSE: %.1f%% → %.1f%%，F30: %.1f%% → %.1f%%\n",
              zang_pop$rRMSE_pct, zang_ind$rRMSE_pct,
              zang_pop$F30_pct, zang_ind$F30_pct))
} else {
  cat("  ⚠️  zang2021（真实模型）：MAPB 改善效果不明显\n")
}

## 判断外部模型
external_improved <- 0
for (m in c("li", "yin", "sun")) {
  ind_row <- calibration_metrics %>% filter(model == m, mode == "IND")
  pop_row <- calibration_metrics %>% filter(model == m, mode == "POP")
  
  ## 判断：IND 的 |rBias| 更小，或 rRMSE 更小
  if (abs(ind_row$rBias_pct) < abs(pop_row$rBias_pct) ||
      ind_row$rRMSE_pct < pop_row$rRMSE_pct) {
    external_improved <- external_improved + 1
  }
}

if (external_improved == 3) {
  cat("  ✅ 外部模型（li/yin/sun）：MAPB 使 CL 估计向 CL_true 方向收敛\n")
} else if (external_improved >= 1) {
  cat(sprintf("  ⚠️  外部模型：%d/3 个模型的 CL 估计被 MAPB 改善\n", external_improved))
} else {
  cat("  ❌ 外部模型：MAPB 未能改善 CL 估计\n")
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")

## ############################################################################
## F.6 完成
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("✅ Part F 全部完成（CL 校准分析）")
message(paste(rep("=", 70), collapse = ""))

message("\n📄 新增输出文件：")
message("  数据：")
message("  - evaluation_results.rds        （已追加 CL_calibration）")
message("  - CL_calibration_metrics.csv")
message("\n  图表（plots/）：")
message("  - fig11_BA_zang.png")
message("  - fig11_BA_li.png")
message("  - fig11_BA_yin.png")
message("  - fig11_BA_sun.png")
message("  - fig12_CL_calibration_barplot.png")

message("\n🔗 使用方法：")
message('  results <- readRDS("outputs9/evaluation_results.rds")')
message('  CL_cal <- results$CL_calibration')
message('  CL_cal$calibration_metrics   # rBias/rRMSE/F30')
message('  CL_cal$BA_stats              # Bland-Altman 统计量')
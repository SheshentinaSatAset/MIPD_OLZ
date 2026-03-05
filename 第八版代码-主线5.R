## ############################################################################
## 第八版代码-主线5.R
## 
## Stage 5：依从性判别评估
## 
## 核心任务：
##   1. 表格（1）：分场景统计（accuracy, margin）
##   2. 表格（2）：分患者汇总统计
##   3. 图表（3）：Accuracy vs |η_true| 分析
##   4. 统计检验：McNemar检验（Accuracy）、配对t检验（Margin, Posterior_true）
## 
## 输入路径：outputs8/outputs8_stage4/
## 输出路径：outputs8/outputs8_stage5/
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("Stage 5：依从性判别评估（第八版）")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## 5.0 加载依赖与路径设置
## ############################################################################

message("\n--- 5.0 加载依赖与路径设置 ---")

## 加载包
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(readr)
  library(scales)
})

## 路径设置
PROJ_ROOT <- "D:/Users/YujiaZhang/Desktop/26救赎之道/25.9.29 毕业设计/5正式实施！/多成人模型/体现我卓越的项目管理能力"

## 输入路径（第八版 Stage 4 输出）
INPUT_ROOT <- file.path(PROJ_ROOT, "outputs8")
IN_DIR_S4 <- file.path(INPUT_ROOT, "outputs8_stage4")

## 输出路径
OUT_DIR_S5 <- file.path(INPUT_ROOT, "outputs8_stage5")

## 创建输出目录
if (!dir.exists(OUT_DIR_S5)) {
  dir.create(OUT_DIR_S5, recursive = TRUE)
  message("📁 已创建目录: ", OUT_DIR_S5)
}

## 创建图表子目录
plot_dir <- file.path(OUT_DIR_S5, "plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

## ############################################################################
## 5.1 读取 Stage 4 输出
## ############################################################################

message("\n--- 5.1 读取 Stage 4 输出 ---")

posterior_summary <- readRDS(file.path(IN_DIR_S4, "posterior_summary.rds"))
posterior_raw <- readRDS(file.path(IN_DIR_S4, "posterior_raw.rds"))
scenario_library <- readRDS(file.path(IN_DIR_S4, "scenario_library.rds"))
params <- readRDS(file.path(IN_DIR_S4, "params.rds"))

message("✅ Stage 4 输出已加载")
message("   posterior_summary: ", nrow(posterior_summary), " 行")
message("   posterior_raw: ", nrow(posterior_raw), " 行")
message("   场景数: ", nrow(scenario_library))
message("   分箱参数: Δ = ", params$DELTA, " ng/mL, M = ", params$M_SIM)

## ############################################################################
## 5.2 计算 Margin（后验概率最大 - 第二大）
## ############################################################################

message("\n--- 5.2 计算 Margin ---")

## 从 posterior_raw 中计算每个测试的 margin
margin_data <- posterior_raw %>%
  group_by(patient_id, pop, model_name, mc_iter, s_true) %>%
  summarise(
    ## 方法 A (Ind)
    posterior_sorted_ind = list(sort(posterior_ind, decreasing = TRUE)),
    margin_ind = posterior_sorted_ind[[1]][1] - posterior_sorted_ind[[1]][2],
    
    ## 方法 B (Pop)
    posterior_sorted_pop = list(sort(posterior_pop, decreasing = TRUE)),
    margin_pop = posterior_sorted_pop[[1]][1] - posterior_sorted_pop[[1]][2],
    
    .groups = "drop"
  ) %>%
  select(-posterior_sorted_ind, -posterior_sorted_pop)

## 合并到 posterior_summary
posterior_summary <- posterior_summary %>%
  left_join(margin_data, by = c("patient_id", "pop", "model_name", "mc_iter", "s_true"))

message("✅ Margin 计算完成")

## 检查合并结果
cat("\n=== Margin 统计摘要 ===\n")
cat("Margin_ind: ", sprintf("mean=%.3f, sd=%.3f", 
                             mean(posterior_summary$margin_ind, na.rm = TRUE),
                             sd(posterior_summary$margin_ind, na.rm = TRUE)), "\n")
cat("Margin_pop: ", sprintf("mean=%.3f, sd=%.3f", 
                             mean(posterior_summary$margin_pop, na.rm = TRUE),
                             sd(posterior_summary$margin_pop, na.rm = TRUE)), "\n")

## ############################################################################
## 5.3 表格（1）：分场景统计
## ############################################################################

message("\n--- 5.3 生成表格（1）：分场景统计 ---")

## 计算分场景统计（长格式，每行一个方法）
accuracy_by_scenario <- posterior_summary %>%
  group_by(patient_id, model_name, s_true) %>%
  summarise(
    n = n(),
    
    ## 方法 A (Ind)
    accuracy_ind = mean(match_ind) * 100,
    margin_mean_ind = mean(margin_ind, na.rm = TRUE),
    margin_sd_ind = sd(margin_ind, na.rm = TRUE),
    
    ## 方法 B (Pop)
    accuracy_pop = mean(match_pop) * 100,
    margin_mean_pop = mean(margin_pop, na.rm = TRUE),
    margin_sd_pop = sd(margin_pop, na.rm = TRUE),
    
    .groups = "drop"
  )

## 转换为长格式（一行IND，一行POP）
table1_ind <- accuracy_by_scenario %>%
  select(patient_id, model_name, s_true, n, 
         accuracy = accuracy_ind, 
         margin_mean = margin_mean_ind, 
         margin_sd = margin_sd_ind) %>%
  mutate(method = "IND")

table1_pop <- accuracy_by_scenario %>%
  select(patient_id, model_name, s_true, n, 
         accuracy = accuracy_pop, 
         margin_mean = margin_mean_pop, 
         margin_sd = margin_sd_pop) %>%
  mutate(method = "POP")

table1 <- bind_rows(table1_ind, table1_pop) %>%
  arrange(patient_id, s_true, desc(method)) %>%
  select(patient_id, model_name, s_true, method, n, accuracy, margin_mean, margin_sd)

cat("\n=== 表格（1）：分场景统计（前20行预览）===\n")
print(head(table1, 20))

message("✅ 表格（1）生成完成，共 ", nrow(table1), " 行")

## ############################################################################
## 5.4 表格（2）：分患者汇总统计
## ############################################################################

message("\n--- 5.4 生成表格（2）：分患者汇总统计 ---")

## 计算分患者汇总
accuracy_by_patient <- posterior_summary %>%
  group_by(patient_id, model_name) %>%
  summarise(
    n = n(),
    
    ## 方法 A (Ind)
    accuracy_ind = mean(match_ind) * 100,
    margin_mean_ind = mean(margin_ind, na.rm = TRUE),
    margin_sd_ind = sd(margin_ind, na.rm = TRUE),
    
    ## 方法 B (Pop)
    accuracy_pop = mean(match_pop) * 100,
    margin_mean_pop = mean(margin_pop, na.rm = TRUE),
    margin_sd_pop = sd(margin_pop, na.rm = TRUE),
    
    .groups = "drop"
  )

## 转换为长格式
table2_ind <- accuracy_by_patient %>%
  select(patient_id, model_name, n, 
         accuracy = accuracy_ind, 
         margin_mean = margin_mean_ind, 
         margin_sd = margin_sd_ind) %>%
  mutate(method = "IND")

table2_pop <- accuracy_by_patient %>%
  select(patient_id, model_name, n, 
         accuracy = accuracy_pop, 
         margin_mean = margin_mean_pop, 
         margin_sd = margin_sd_pop) %>%
  mutate(method = "POP")

table2 <- bind_rows(table2_ind, table2_pop) %>%
  arrange(patient_id, desc(method)) %>%
  select(patient_id, model_name, method, n, accuracy, margin_mean, margin_sd)

cat("\n=== 表格（2）：分患者汇总统计 ===\n")
print(table2)

message("✅ 表格（2）生成完成，共 ", nrow(table2), " 行")

## ############################################################################
## 5.5 读取 Stage 3 数据（获取 η_true）
## ############################################################################

message("\n--- 5.5 读取 Stage 3 数据（获取 η_true）---")

## 读取 Stage 3 原始结果（第六版输出）
IN_DIR_S3 <- file.path(PROJ_ROOT, "outputs6/outputs6_stage3")
mc_raw_main <- readRDS(file.path(IN_DIR_S3, "mc_raw_results_main.rds"))
mc_raw_side <- readRDS(file.path(IN_DIR_S3, "mc_raw_results_side.rds"))
mc_raw_all <- bind_rows(mc_raw_main, mc_raw_side)

## 提取 η_true 数据
eta_true_data <- mc_raw_all %>%
  select(patient_id, model_name, mc_iter, eta_CL_true, eta_V_true)

## 合并到 posterior_summary
posterior_summary <- posterior_summary %>%
  left_join(eta_true_data, by = c("patient_id", "model_name", "mc_iter"))

message("✅ η_true 数据已合并")

## ############################################################################
## 5.6 图表（3A）：Accuracy vs |η_CL_true| 散点图 + loess
## ############################################################################

message("\n--- 5.6 生成图表（3A）：Accuracy vs |η_CL_true| ---")

## 准备绘图数据
plot_data_CL <- posterior_summary %>%
  mutate(abs_eta_CL_true = abs(eta_CL_true)) %>%
  mutate(
    correct_ind = as.numeric(match_ind),
    correct_pop = as.numeric(match_pop)
  )

## 绘制 η_CL 的图（所有模型都有）
p_accuracy_eta_CL <- ggplot(plot_data_CL, aes(x = abs_eta_CL_true)) +
  ## Pop 方法
  geom_smooth(aes(y = correct_pop, color = "Pop"), 
              method = "loess", se = TRUE, alpha = 0.2, span = 0.5) +
  ## Ind 方法
  geom_smooth(aes(y = correct_ind, color = "Ind"), 
              method = "loess", se = TRUE, alpha = 0.2, span = 0.5) +
  facet_wrap(~ patient_id, ncol = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_color_manual(
    values = c("Ind" = "steelblue", "Pop" = "coral"),
    labels = c("Ind" = "方法A (η_ind)", "Pop" = "方法B (η_pop=0)")
  ) +
  labs(
    title = "准确率 vs |η_CL_true|（η_CL 参数）",
    subtitle = paste0("第八版（分箱蒙特卡洛法，Δ=", params$DELTA, " ng/mL, M=", params$M_SIM, "）"),
    x = "|η_CL_true|",
    y = "准确率",
    color = "方法"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "bottom"
  )

ggsave(file.path(plot_dir, "accuracy_vs_eta_CL_loess.png"), 
       p_accuracy_eta_CL, width = 12, height = 10, dpi = 300)

message("  ✓ accuracy_vs_eta_CL_loess.png 已保存")

## ############################################################################
## 5.7 图表（3A续）：Accuracy vs |η_V_true|（仅 Acceptable 模型）
## ############################################################################

message("\n--- 5.7 生成图表（3A续）：Accuracy vs |η_V_true| ---")

## 仅选择 η_V Acceptable 的模型：zang2021, li2018
## yin2016 和 sun2021 的 η_V Unacceptable，不绘制

plot_data_V <- posterior_summary %>%
  filter(model_name %in% c("zang2021", "li2018")) %>%
  filter(!is.na(eta_V_true)) %>%
  mutate(
    abs_eta_V_true = abs(eta_V_true),
    correct_ind = as.numeric(match_ind),
    correct_pop = as.numeric(match_pop)
  )

if (nrow(plot_data_V) > 0) {
  
  p_accuracy_eta_V <- ggplot(plot_data_V, aes(x = abs_eta_V_true)) +
    geom_smooth(aes(y = correct_pop, color = "Pop"), 
                method = "loess", se = TRUE, alpha = 0.2, span = 0.5) +
    geom_smooth(aes(y = correct_ind, color = "Ind"), 
                method = "loess", se = TRUE, alpha = 0.2, span = 0.5) +
    facet_wrap(~ patient_id, ncol = 3) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_color_manual(
      values = c("Ind" = "steelblue", "Pop" = "coral"),
      labels = c("Ind" = "方法A (η_ind)", "Pop" = "方法B (η_pop=0)")
    ) +
    labs(
      title = "准确率 vs |η_V_true|（η_V/Vc 参数，仅 Acceptable 模型）",
      subtitle = "zang2021, li2018 | Loess 平滑曲线 + 95% 置信区间",
      x = "|η_V_true|",
      y = "准确率",
      color = "方法"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      legend.position = "bottom"
    )
  
  ggsave(file.path(plot_dir, "accuracy_vs_eta_V_loess.png"), 
         p_accuracy_eta_V, width = 10, height = 6, dpi = 300)
  
  message("  ✓ accuracy_vs_eta_V_loess.png 已保存")
  
} else {
  message("  ⚠️ 无 η_V Acceptable 的数据，跳过图表")
}

## ############################################################################
## 5.8 图表（3B）：Accuracy vs |η_true| 分组柱状图
## ############################################################################

message("\n--- 5.8 生成图表（3B）：分组柱状图 ---")

## 定义 |η_CL_true| 分组
eta_breaks <- c(0, 0.5, 1.0, 1.5, Inf)
eta_labels <- c("[0, 0.5)", "[0.5, 1.0)", "[1.0, 1.5)", "[1.5, +∞)")

## 计算分组准确率
accuracy_by_eta_bin <- posterior_summary %>%
  mutate(
    eta_CL_bin = cut(abs(eta_CL_true), breaks = eta_breaks, labels = eta_labels, right = FALSE)
  ) %>%
  group_by(patient_id, model_name, eta_CL_bin) %>%
  summarise(
    n = n(),
    accuracy_ind = mean(match_ind) * 100,
    accuracy_pop = mean(match_pop) * 100,
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(accuracy_ind, accuracy_pop), 
               names_to = "method", values_to = "accuracy") %>%
  mutate(
    method = ifelse(method == "accuracy_ind", "Ind", "Pop")
  )

## 绘制分组柱状图
p_accuracy_eta_bar <- ggplot(accuracy_by_eta_bin, 
                              aes(x = eta_CL_bin, y = accuracy, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", accuracy)), 
            position = position_dodge(width = 0.8), 
            vjust = -0.5, size = 2.5) +
  facet_wrap(~ patient_id, ncol = 3) +
  scale_fill_manual(
    values = c("Ind" = "steelblue", "Pop" = "coral"),
    labels = c("Ind" = "方法A (η_ind)", "Pop" = "方法B (η_pop=0)")
  ) +
  scale_y_continuous(limits = c(0, 105)) +
  labs(
    title = "准确率 vs |η_CL_true| 分组",
    subtitle = paste0("第八版（分箱蒙特卡洛法，Δ=", params$DELTA, " ng/mL, M=", params$M_SIM, "）"),
    x = "|η_CL_true| 区间",
    y = "准确率 (%)",
    fill = "方法"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(plot_dir, "accuracy_vs_eta_CL_barplot.png"), 
       p_accuracy_eta_bar, width = 12, height = 10, dpi = 300)

message("  ✓ accuracy_vs_eta_CL_barplot.png 已保存")

## ############################################################################
## 5.9 统计检验
## ############################################################################

message("\n--- 5.9 统计检验 ---")

## ----------------------------------------------------------------------------
## 5.9.1 检验1：McNemar 检验（Accuracy）
## ----------------------------------------------------------------------------
message("\n  --- 检验1：McNemar 检验（Accuracy）---")

run_mcnemar_test <- function(data, group_name = "总体") {
  
  ## 构建 2×2 混淆表
  confusion <- table(
    Ind_correct = data$match_ind,
    Pop_correct = data$match_pop
  )
  
  ## 执行 McNemar 检验
  test_result <- tryCatch(
    mcnemar.test(confusion),
    error = function(e) NULL
  )
  
  if (is.null(test_result)) {
    return(tibble(
      group = group_name,
      n_both_correct = NA_integer_,
      n_ind_only = NA_integer_,
      n_pop_only = NA_integer_,
      n_both_wrong = NA_integer_,
      chi_squared = NA_real_,
      p_value = NA_real_
    ))
  }
  
  tibble(
    group = group_name,
    n_both_correct = confusion[2, 2],
    n_ind_only = confusion[2, 1],
    n_pop_only = confusion[1, 2],
    n_both_wrong = confusion[1, 1],
    chi_squared = as.numeric(test_result$statistic),
    p_value = test_result$p.value
  )
}

## 总体检验
mcnemar_overall <- run_mcnemar_test(posterior_summary, "总体")

## 按患者分层检验
mcnemar_by_patient <- posterior_summary %>%
  group_by(patient_id) %>%
  group_modify(~ run_mcnemar_test(.x, .y$patient_id)) %>%
  ungroup()

## 按场景分层检验
mcnemar_by_scenario <- posterior_summary %>%
  group_by(s_true) %>%
  group_modify(~ run_mcnemar_test(.x, .y$s_true)) %>%
  ungroup()

## 合并结果
mcnemar_results <- bind_rows(
  mcnemar_overall %>% mutate(level = "总体"),
  mcnemar_by_patient %>% mutate(level = "患者"),
  mcnemar_by_scenario %>% mutate(level = "场景")
) %>%
  select(level, group, everything())

cat("\n=== McNemar 检验结果 ===\n")
print(mcnemar_results)

## ----------------------------------------------------------------------------
## 5.9.2 检验2：配对 t 检验（Margin）
## ----------------------------------------------------------------------------
message("\n  --- 检验2：配对 t 检验（Margin）---")

run_paired_t_test <- function(data, var_ind, var_pop, group_name = "总体") {
  
  x <- data[[var_ind]]
  y <- data[[var_pop]]
  
  ## 过滤 NA
  valid <- !is.na(x) & !is.na(y)
  x <- x[valid]
  y <- y[valid]
  
  if (length(x) < 2) {
    return(tibble(
      group = group_name,
      variable = var_ind,
      mean_ind = NA_real_,
      mean_pop = NA_real_,
      mean_diff = NA_real_,
      t_statistic = NA_real_,
      p_value = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    ))
  }
  
  test_result <- t.test(x, y, paired = TRUE)
  
  tibble(
    group = group_name,
    variable = gsub("_ind$", "", var_ind),
    mean_ind = mean(x),
    mean_pop = mean(y),
    mean_diff = mean(x - y),
    t_statistic = as.numeric(test_result$statistic),
    p_value = test_result$p.value,
    ci_lower = test_result$conf.int[1],
    ci_upper = test_result$conf.int[2]
  )
}

## 总体检验：Margin
ttest_margin_overall <- run_paired_t_test(posterior_summary, "margin_ind", "margin_pop", "总体")

## 按患者分层检验
ttest_margin_by_patient <- posterior_summary %>%
  group_by(patient_id) %>%
  group_modify(~ run_paired_t_test(.x, "margin_ind", "margin_pop", .y$patient_id)) %>%
  ungroup()

## 按场景分层检验
ttest_margin_by_scenario <- posterior_summary %>%
  group_by(s_true) %>%
  group_modify(~ run_paired_t_test(.x, "margin_ind", "margin_pop", .y$s_true)) %>%
  ungroup()

## 合并结果
ttest_margin_results <- bind_rows(
  ttest_margin_overall %>% mutate(level = "总体"),
  ttest_margin_by_patient %>% mutate(level = "患者"),
  ttest_margin_by_scenario %>% mutate(level = "场景")
) %>%
  select(level, group, everything())

cat("\n=== 配对 t 检验结果（Margin）===\n")
print(ttest_margin_results)

## ----------------------------------------------------------------------------
## 5.9.3 检验3：配对 t 检验（Posterior_true）
## ----------------------------------------------------------------------------
message("\n  --- 检验3：配对 t 检验（Posterior_true）---")

## 总体检验
ttest_posterior_overall <- run_paired_t_test(
  posterior_summary, "posterior_true_ind", "posterior_true_pop", "总体"
)

## 按患者分层检验
ttest_posterior_by_patient <- posterior_summary %>%
  group_by(patient_id) %>%
  group_modify(~ run_paired_t_test(.x, "posterior_true_ind", "posterior_true_pop", .y$patient_id)) %>%
  ungroup()

## 按场景分层检验
ttest_posterior_by_scenario <- posterior_summary %>%
  group_by(s_true) %>%
  group_modify(~ run_paired_t_test(.x, "posterior_true_ind", "posterior_true_pop", .y$s_true)) %>%
  ungroup()

## 合并结果
ttest_posterior_results <- bind_rows(
  ttest_posterior_overall %>% mutate(level = "总体"),
  ttest_posterior_by_patient %>% mutate(level = "患者"),
  ttest_posterior_by_scenario %>% mutate(level = "场景")
) %>%
  select(level, group, everything())

cat("\n=== 配对 t 检验结果（Posterior_true）===\n")
print(ttest_posterior_results)

## ############################################################################
## 5.10 保存输出文件
## ############################################################################

message("\n--- 5.10 保存输出文件 ---")

## 保存表格（1）
write_csv(table1, file.path(OUT_DIR_S5, "table1_accuracy_by_scenario.csv"))
saveRDS(table1, file.path(OUT_DIR_S5, "table1_accuracy_by_scenario.rds"))
message("  ✓ table1_accuracy_by_scenario.csv/.rds 已保存（", nrow(table1), " 行）")

## 保存表格（2）
write_csv(table2, file.path(OUT_DIR_S5, "table2_accuracy_by_patient.csv"))
saveRDS(table2, file.path(OUT_DIR_S5, "table2_accuracy_by_patient.rds"))
message("  ✓ table2_accuracy_by_patient.csv/.rds 已保存（", nrow(table2), " 行）")

## 保存分组准确率数据
write_csv(accuracy_by_eta_bin, file.path(OUT_DIR_S5, "accuracy_by_eta_bin.csv"))
saveRDS(accuracy_by_eta_bin, file.path(OUT_DIR_S5, "accuracy_by_eta_bin.rds"))
message("  ✓ accuracy_by_eta_bin.csv/.rds 已保存")

## 保存统计检验结果
stat_tests <- list(
  mcnemar = mcnemar_results,
  ttest_margin = ttest_margin_results,
  ttest_posterior = ttest_posterior_results
)
saveRDS(stat_tests, file.path(OUT_DIR_S5, "statistical_tests.rds"))

write_csv(mcnemar_results, file.path(OUT_DIR_S5, "test_mcnemar.csv"))
write_csv(ttest_margin_results, file.path(OUT_DIR_S5, "test_ttest_margin.csv"))
write_csv(ttest_posterior_results, file.path(OUT_DIR_S5, "test_ttest_posterior.csv"))
message("  ✓ statistical_tests.rds 及各 CSV 已保存")

## 保存更新后的 posterior_summary（含 margin 和 η_true）
saveRDS(posterior_summary, file.path(OUT_DIR_S5, "posterior_summary_enriched.rds"))
message("  ✓ posterior_summary_enriched.rds 已保存")

## ############################################################################
## 5.11 生成汇总报告
## ############################################################################

message("\n--- 5.11 生成汇总报告 ---")

## 总体统计
overall_stats <- posterior_summary %>%
  summarise(
    n_total = n(),
    accuracy_ind = mean(match_ind) * 100,
    accuracy_pop = mean(match_pop) * 100,
    accuracy_diff = accuracy_ind - accuracy_pop,
    margin_mean_ind = mean(margin_ind, na.rm = TRUE),
    margin_mean_pop = mean(margin_pop, na.rm = TRUE),
    margin_diff = margin_mean_ind - margin_mean_pop,
    posterior_true_mean_ind = mean(posterior_true_ind, na.rm = TRUE),
    posterior_true_mean_pop = mean(posterior_true_pop, na.rm = TRUE)
  )

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("                    汇总报告（第八版）\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

cat("\n【分箱蒙特卡洛参数】\n")
cat(sprintf("  分箱步长 Δ = %d ng/mL\n", params$DELTA))
cat(sprintf("  每场景模拟次数 M = %d\n", params$M_SIM))
cat(sprintf("  Stage 4 耗时 = %.2f 分钟\n", params$elapsed_minutes))

cat("\n【总体准确率】\n")
cat(sprintf("  方法A (Ind): %.2f%%\n", overall_stats$accuracy_ind))
cat(sprintf("  方法B (Pop): %.2f%%\n", overall_stats$accuracy_pop))
cat(sprintf("  提升幅度:    %+.2f%%\n", overall_stats$accuracy_diff))

cat("\n【平均 Margin（判别置信度）】\n")
cat(sprintf("  方法A (Ind): %.4f\n", overall_stats$margin_mean_ind))
cat(sprintf("  方法B (Pop): %.4f\n", overall_stats$margin_mean_pop))
cat(sprintf("  差异:        %+.4f\n", overall_stats$margin_diff))

cat("\n【McNemar 检验（总体）】\n")
cat(sprintf("  χ² = %.2f, p = %.2e\n", 
            mcnemar_overall$chi_squared, mcnemar_overall$p_value))
if (mcnemar_overall$p_value < 0.001) {
  cat("  结论：Ind 方法显著优于 Pop 方法 (p < 0.001)\n")
} else if (mcnemar_overall$p_value < 0.05) {
  cat("  结论：Ind 方法显著优于 Pop 方法 (p < 0.05)\n")
} else {
  cat("  结论：两种方法无显著差异 (p >= 0.05)\n")
}

cat("\n【配对 t 检验（Margin，总体）】\n")
cat(sprintf("  t = %.2f, p = %.2e\n", 
            ttest_margin_overall$t_statistic, ttest_margin_overall$p_value))
cat(sprintf("  95%% CI: [%.4f, %.4f]\n", 
            ttest_margin_overall$ci_lower, ttest_margin_overall$ci_upper))

cat("\n", paste(rep("=", 60), collapse = ""), "\n")

## 保存汇总报告
summary_report <- list(
  overall_stats = overall_stats,
  mcnemar_overall = mcnemar_overall,
  ttest_margin_overall = ttest_margin_overall,
  ttest_posterior_overall = ttest_posterior_overall,
  params = params
)
saveRDS(summary_report, file.path(OUT_DIR_S5, "summary_report.rds"))

## ############################################################################
## 5.12 Stage 5 完成
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("✅ Stage 5 全部完成（第八版）")
message(paste(rep("=", 70), collapse = ""))

message("\n📁 输出目录: ", OUT_DIR_S5)
message("\n📄 输出文件：")
message("  表格：")
message("  - table1_accuracy_by_scenario.csv/.rds （分场景统计，", nrow(table1), " 行）")
message("  - table2_accuracy_by_patient.csv/.rds  （分患者统计，", nrow(table2), " 行）")
message("  - accuracy_by_eta_bin.csv/.rds         （分组准确率）")
message("\n  图表：")
message("  - plots/accuracy_vs_eta_CL_loess.png   （η_CL 散点图）")
message("  - plots/accuracy_vs_eta_V_loess.png    （η_V 散点图，仅 Acceptable）")
message("  - plots/accuracy_vs_eta_CL_barplot.png （分组柱状图）")
message("\n  统计检验：")
message("  - test_mcnemar.csv                     （McNemar 检验）")
message("  - test_ttest_margin.csv                （配对t检验-Margin）")
message("  - test_ttest_posterior.csv             （配对t检验-Posterior）")
message("  - statistical_tests.rds                （全部检验结果）")
message("\n  其他：")
message("  - posterior_summary_enriched.rds       （增强版汇总表）")
message("  - summary_report.rds                   （汇总报告）")

message("\n📊 总体结果（第八版 - 分箱蒙特卡洛法）：")
message(sprintf("   方法A (Ind) 准确率: %.2f%%", overall_stats$accuracy_ind))
message(sprintf("   方法B (Pop) 准确率: %.2f%%", overall_stats$accuracy_pop))
message(sprintf("   提升幅度: %+.2f%%", overall_stats$accuracy_diff))

message("\n🔗 后续使用方法：")
message('  table1 <- read_csv("outputs8/outputs8_stage5/table1_accuracy_by_scenario.csv")')
message('  tests <- readRDS("outputs8/outputs8_stage5/statistical_tests.rds")')
message('  report <- readRDS("outputs8/outputs8_stage5/summary_report.rds")')

message("\n✅ 第八版主线代码全部完成！")
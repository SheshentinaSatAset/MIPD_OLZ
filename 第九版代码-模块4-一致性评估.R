## ############################################################################
## 第九版代码-模块4-一致性评估.R（修复版）
## 
## 修复内容：
##   - 添加 lme4 包依赖检查
##   - 提供备用 ICC 计算方法（如果 lme4 不可用）
## 
## 核心任务：
##   评估跨模型一致性和单模型准确率
##   - Part A：单模型准确率评估（IND vs POP）
##   - Part B：跨模型一致性评估（Agreement、Fleiss' κ）
##   - Part C：η_ind 相关性分析
## 
## 输入：
##   - outputs9/discrimination_results.rds
##   - outputs9/data_true.rds
##   - outputs9/mapb_results.rds
##   - outputs9/scenario_library.rds
## 
## 输出：
##   - outputs9/evaluation_results.rds
##   - outputs9/plots/ 目录下的图表
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("第九版 模块4：一致性评估（修复版）")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## 4.0 加载依赖与路径设置
## ############################################################################

message("\n--- 4.0 加载依赖与路径设置 ---")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(irr)      # Fleiss' Kappa
  library(lme4)
})

## 检查并加载 psych 和 lme4（ICC 计算需要）
HAS_LME4 <- requireNamespace("lme4", quietly = TRUE)
HAS_PSYCH <- requireNamespace("psych", quietly = TRUE)

if (HAS_PSYCH && HAS_LME4) {
  library(psych)
  message("✅ psych 和 lme4 包已加载（ICC 计算可用）")
} else if (!HAS_LME4) {
  message("⚠️ lme4 包未安装，将使用简化 ICC 计算方法")
  message("   如需完整 ICC 功能，请运行：install.packages('lme4')")
}

## 路径设置
PROJ_ROOT <- "D:/Users/YujiaZhang/Desktop/26救赎之道/25.9.29 毕业设计/5正式实施！/多成人模型/体现我卓越的项目管理能力"
OUT_DIR_9 <- file.path(PROJ_ROOT, "outputs9")

## 创建图表目录
plot_dir <- file.path(OUT_DIR_9, "plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
  message("📁 已创建目录: ", plot_dir)
}

## 读取前置模块输出
discrimination_results <- readRDS(file.path(OUT_DIR_9, "discrimination_results.rds"))
data_true <- readRDS(file.path(OUT_DIR_9, "data_true.rds"))
mapb_results <- readRDS(file.path(OUT_DIR_9, "mapb_results.rds"))
scenario_library <- readRDS(file.path(OUT_DIR_9, "scenario_library.rds"))
module1_params <- readRDS(file.path(OUT_DIR_9, "module1_params.rds"))

N_MC <- module1_params$N_MC

message("✅ 前置模块输出已加载")
message("   虚拟患者数: ", N_MC)
message("   场景数: ", nrow(scenario_library))

## ############################################################################
## 4.1 合并数据
## ############################################################################

message("\n--- 4.1 合并数据 ---")

## 合并所有数据
combined_data <- discrimination_results %>%
  left_join(
    mapb_results %>% select(mc_iter, 
                            zang_eta_CL_ind, zang_eta_V_ind,
                            li_eta_CL_ind, li_eta_Vc_ind,
                            yin_eta_CL_ind, sun_eta_CL_ind),
    by = "mc_iter"
  )

message("✅ 数据合并完成")
message("   总行数: ", nrow(combined_data))

## ############################################################################
## Part A：单模型准确率评估
## ############################################################################

message("\n", paste(rep("-", 70), collapse = ""))
message("Part A：单模型准确率评估")
message(paste(rep("-", 70), collapse = ""))

## ----------------------------------------------------------------------------
## A.1 各模型准确率
## ----------------------------------------------------------------------------

message("\n--- A.1 各模型准确率 ---")

accuracy_by_model <- tibble(
  model = c("zang2021", "li2018", "yin2016", "sun2021"),
  accuracy_ind = c(
    mean(combined_data$zang_match_ind) * 100,
    mean(combined_data$li_match_ind) * 100,
    mean(combined_data$yin_match_ind) * 100,
    mean(combined_data$sun_match_ind) * 100
  ),
  accuracy_pop = c(
    mean(combined_data$zang_match_pop) * 100,
    mean(combined_data$li_match_pop) * 100,
    mean(combined_data$yin_match_pop) * 100,
    mean(combined_data$sun_match_pop) * 100
  )
) %>%
  mutate(
    delta = accuracy_ind - accuracy_pop,
    improvement_pct = (accuracy_ind - accuracy_pop) / accuracy_pop * 100
  )

cat("\n=== 各模型准确率 ===\n")
print(accuracy_by_model)

## 准确率 SD（跨模型）
accuracy_sd_ind <- sd(accuracy_by_model$accuracy_ind)
accuracy_sd_pop <- sd(accuracy_by_model$accuracy_pop)

cat("\n准确率 SD（跨模型）:\n")
cat("  IND 模式:", round(accuracy_sd_ind, 2), "%\n")
cat("  POP 模式:", round(accuracy_sd_pop, 2), "%\n")
cat("  → IND 模式让各模型表现更均衡:", 
    ifelse(accuracy_sd_ind < accuracy_sd_pop, "是", "否"), "\n")

## ----------------------------------------------------------------------------
## A.2 McNemar 检验
## ----------------------------------------------------------------------------

message("\n--- A.2 McNemar 检验 ---")

run_mcnemar_test <- function(match_ind, match_pop, model_name) {
  
  confusion <- table(
    IND_correct = match_ind,
    POP_correct = match_pop
  )
  
  test_result <- tryCatch(
    mcnemar.test(confusion),
    error = function(e) NULL
  )
  
  if (is.null(test_result)) {
    return(tibble(
      model = model_name,
      n_both_correct = NA_integer_,
      n_ind_only = NA_integer_,
      n_pop_only = NA_integer_,
      n_both_wrong = NA_integer_,
      chi_squared = NA_real_,
      p_value = NA_real_
    ))
  }
  
  tibble(
    model = model_name,
    n_both_correct = confusion[2, 2],
    n_ind_only = confusion[2, 1],
    n_pop_only = confusion[1, 2],
    n_both_wrong = confusion[1, 1],
    chi_squared = as.numeric(test_result$statistic),
    p_value = test_result$p.value
  )
}

mcnemar_results <- bind_rows(
  run_mcnemar_test(combined_data$zang_match_ind, combined_data$zang_match_pop, "zang2021"),
  run_mcnemar_test(combined_data$li_match_ind, combined_data$li_match_pop, "li2018"),
  run_mcnemar_test(combined_data$yin_match_ind, combined_data$yin_match_pop, "yin2016"),
  run_mcnemar_test(combined_data$sun_match_ind, combined_data$sun_match_pop, "sun2021")
)

cat("\n=== McNemar 检验结果 ===\n")
print(mcnemar_results)

## ----------------------------------------------------------------------------
## A.3 按场景分层的准确率
## ----------------------------------------------------------------------------

message("\n--- A.3 按场景分层的准确率 ---")

accuracy_by_scenario <- combined_data %>%
  group_by(s_true) %>%
  summarise(
    n = n(),
    
    ## 各模型 IND 准确率
    zang_ind = mean(zang_match_ind) * 100,
    li_ind = mean(li_match_ind) * 100,
    yin_ind = mean(yin_match_ind) * 100,
    sun_ind = mean(sun_match_ind) * 100,
    
    ## 各模型 POP 准确率
    zang_pop = mean(zang_match_pop) * 100,
    li_pop = mean(li_match_pop) * 100,
    yin_pop = mean(yin_match_pop) * 100,
    sun_pop = mean(sun_match_pop) * 100,
    
    .groups = "drop"
  )

cat("\n=== 按场景分层的准确率 ===\n")
print(accuracy_by_scenario)

## ############################################################################
## Part B：跨模型一致性评估（第九版核心）
## ############################################################################

message("\n", paste(rep("-", 70), collapse = ""))
message("Part B：跨模型一致性评估（第九版核心）")
message(paste(rep("-", 70), collapse = ""))

## ----------------------------------------------------------------------------
## B.1 完全一致率
## ----------------------------------------------------------------------------

message("\n--- B.1 完全一致率 ---")

## 计算 4 个模型是否完全一致（已在模块 3 计算）
agreement_ind <- mean(combined_data$all_agree_ind) * 100
agreement_pop <- mean(combined_data$all_agree_pop) * 100

cat("\n=== 完全一致率（4 个模型给出相同 s_pred）===\n")
cat("  IND 模式:", round(agreement_ind, 2), "%\n")
cat("  POP 模式:", round(agreement_pop, 2), "%\n")
cat("  预期 IND > POP:", ifelse(agreement_ind > agreement_pop, "✅ 符合", "❌ 不符合"), "\n")

## ----------------------------------------------------------------------------
## B.2 Fleiss' Kappa
## ----------------------------------------------------------------------------

message("\n--- B.2 Fleiss' Kappa ---")

## 构建评分矩阵：每行一个患者，每列一个模型的判别结果
## IND 模式
ratings_ind <- combined_data %>%
  select(mc_iter, zang_s_pred_ind, li_s_pred_ind, yin_s_pred_ind, sun_s_pred_ind) %>%
  mutate(across(-mc_iter, ~ factor(.x, levels = scenario_library$scenario_id)))

## 转换为 irr 包需要的格式
ratings_matrix_ind <- as.matrix(ratings_ind[, -1])

## 计算 Fleiss' Kappa
fleiss_ind <- kappam.fleiss(ratings_matrix_ind)

## POP 模式
ratings_pop <- combined_data %>%
  select(mc_iter, zang_s_pred_pop, li_s_pred_pop, yin_s_pred_pop, sun_s_pred_pop) %>%
  mutate(across(-mc_iter, ~ factor(.x, levels = scenario_library$scenario_id)))

ratings_matrix_pop <- as.matrix(ratings_pop[, -1])
fleiss_pop <- kappam.fleiss(ratings_matrix_pop)

cat("\n=== Fleiss' Kappa ===\n")
cat("  IND 模式: κ =", round(fleiss_ind$value, 4), 
    "(p =", format(fleiss_ind$p.value, digits = 3), ")\n")
cat("  POP 模式: κ =", round(fleiss_pop$value, 4),
    "(p =", format(fleiss_pop$p.value, digits = 3), ")\n")
cat("  预期 IND > POP:", ifelse(fleiss_ind$value > fleiss_pop$value, "✅ 符合", "❌ 不符合"), "\n")

## Kappa 解释
interpret_kappa <- function(k) {
  if (k < 0.20) return("Poor")
  if (k < 0.40) return("Fair")
  if (k < 0.60) return("Moderate")
  if (k < 0.80) return("Good")
  return("Excellent")
}

cat("\n  IND 模式一致性:", interpret_kappa(fleiss_ind$value), "\n")
cat("  POP 模式一致性:", interpret_kappa(fleiss_pop$value), "\n")

## ----------------------------------------------------------------------------
## B.3 全对率 / 多数正确率 / 任一正确率
## ----------------------------------------------------------------------------

message("\n--- B.3 多层次正确率 ---")

## 计算各指标
combined_data <- combined_data %>%
  mutate(
    ## IND 模式
    n_correct_ind = zang_match_ind + li_match_ind + yin_match_ind + sun_match_ind,
    majority_correct_ind = (n_correct_ind >= 3),
    any_correct_ind = (n_correct_ind >= 1),
    
    ## POP 模式
    n_correct_pop = zang_match_pop + li_match_pop + yin_match_pop + sun_match_pop,
    majority_correct_pop = (n_correct_pop >= 3),
    any_correct_pop = (n_correct_pop >= 1)
  )

multilevel_accuracy <- tibble(
  metric = c("全对率 (4/4)", "多数正确率 (≥3/4)", "任一正确率 (≥1/4)"),
  ind_pct = c(
    mean(combined_data$all_correct_ind) * 100,
    mean(combined_data$majority_correct_ind) * 100,
    mean(combined_data$any_correct_ind) * 100
  ),
  pop_pct = c(
    mean(combined_data$all_correct_pop) * 100,
    mean(combined_data$majority_correct_pop) * 100,
    mean(combined_data$any_correct_pop) * 100
  )
) %>%
  mutate(delta = ind_pct - pop_pct)

cat("\n=== 多层次正确率 ===\n")
print(multilevel_accuracy)

## ----------------------------------------------------------------------------
## B.4 两两一致率矩阵
## ----------------------------------------------------------------------------

message("\n--- B.4 两两一致率矩阵 ---")

models <- c("zang", "li", "yin", "sun")
n_models <- length(models)

## IND 模式两两一致率
pairwise_agree_ind <- matrix(NA, n_models, n_models, 
                              dimnames = list(models, models))

for (i in 1:n_models) {
  for (j in 1:n_models) {
    col_i <- paste0(models[i], "_s_pred_ind")
    col_j <- paste0(models[j], "_s_pred_ind")
    pairwise_agree_ind[i, j] <- mean(combined_data[[col_i]] == combined_data[[col_j]]) * 100
  }
}

## POP 模式两两一致率
pairwise_agree_pop <- matrix(NA, n_models, n_models,
                              dimnames = list(models, models))

for (i in 1:n_models) {
  for (j in 1:n_models) {
    col_i <- paste0(models[i], "_s_pred_pop")
    col_j <- paste0(models[j], "_s_pred_pop")
    pairwise_agree_pop[i, j] <- mean(combined_data[[col_i]] == combined_data[[col_j]]) * 100
  }
}

cat("\n=== 两两一致率矩阵 (IND) ===\n")
print(round(pairwise_agree_ind, 1))

cat("\n=== 两两一致率矩阵 (POP) ===\n")
print(round(pairwise_agree_pop, 1))

## ----------------------------------------------------------------------------
## B.5 按场景分层的一致率
## ----------------------------------------------------------------------------

message("\n--- B.5 按场景分层的一致率 ---")

agreement_by_scenario <- combined_data %>%
  group_by(s_true) %>%
  summarise(
    n = n(),
    agree_ind = mean(all_agree_ind) * 100,
    agree_pop = mean(all_agree_pop) * 100,
    all_correct_ind = mean(all_correct_ind) * 100,
    all_correct_pop = mean(all_correct_pop) * 100,
    .groups = "drop"
  )

cat("\n=== 按场景分层的一致率 ===\n")
print(agreement_by_scenario)

## ############################################################################
## Part C：η_ind 相关性分析
## ############################################################################

message("\n", paste(rep("-", 70), collapse = ""))
message("Part C：η_ind 相关性分析")
message(paste(rep("-", 70), collapse = ""))

## ----------------------------------------------------------------------------
## C.1 η_CL_ind 跨模型相关系数
## ----------------------------------------------------------------------------

message("\n--- C.1 η_CL_ind 跨模型相关系数 ---")

## 提取 η_CL_ind
eta_CL_data <- combined_data %>%
  select(mc_iter, 
         zang = zang_eta_CL_ind,
         li = li_eta_CL_ind,
         yin = yin_eta_CL_ind,
         sun = sun_eta_CL_ind)

## 计算相关系数矩阵
cor_matrix <- cor(eta_CL_data[, -1], use = "pairwise.complete.obs")

cat("\n=== η_CL_ind 相关系数矩阵 ===\n")
print(round(cor_matrix, 4))

## 提取两两相关系数
eta_correlation <- tibble(
  model_pair = c("zang_li", "zang_yin", "zang_sun", "li_yin", "li_sun", "yin_sun"),
  pearson_r = c(
    cor_matrix["zang", "li"],
    cor_matrix["zang", "yin"],
    cor_matrix["zang", "sun"],
    cor_matrix["li", "yin"],
    cor_matrix["li", "sun"],
    cor_matrix["yin", "sun"]
  )
)

cat("\n=== η_CL_ind 两两相关系数 ===\n")
print(eta_correlation)

## ----------------------------------------------------------------------------
## C.2 ICC（组内相关系数）— 修复版
## ----------------------------------------------------------------------------

message("\n--- C.2 ICC（简化计算）---")

## 简化 ICC 计算函数（不依赖 lme4）
calc_icc_simple <- function(data_matrix) {
  ## 使用方差分量估计 ICC
  ## ICC(2,1) = (BMS - WMS) / (BMS + (k-1)*WMS)
  ## 其中 BMS = Between-subject Mean Square
  ##      WMS = Within-subject Mean Square
  ##      k = 评分者数量
  
  n <- nrow(data_matrix)
  k <- ncol(data_matrix)
  
  ## 计算各方差分量
  grand_mean <- mean(as.matrix(data_matrix), na.rm = TRUE)
  row_means <- rowMeans(data_matrix, na.rm = TRUE)
  col_means <- colMeans(data_matrix, na.rm = TRUE)
  
  ## Between-subject SS
  SS_between <- k * sum((row_means - grand_mean)^2, na.rm = TRUE)
  
  ## Within-subject SS (residual)
  SS_within <- sum((as.matrix(data_matrix) - outer(row_means, rep(1, k)))^2, na.rm = TRUE)
  
  ## Mean Squares
  BMS <- SS_between / (n - 1)
  WMS <- SS_within / (n * (k - 1))
  
  ## ICC(2,1) - Two-way random, single measures
  icc_2_1 <- (BMS - WMS) / (BMS + (k - 1) * WMS)
  
  ## ICC(2,k) - Two-way random, average measures
  icc_2_k <- (BMS - WMS) / BMS
  
  list(
    icc_2_1 = max(0, icc_2_1),  # 确保非负
    icc_2_k = max(0, icc_2_k)
  )
}

## 计算 ICC
if (HAS_PSYCH && HAS_LME4) {
  ## 使用 psych 包的完整 ICC
  icc_result <- ICC(eta_CL_data[, -1])
  icc_2_1 <- icc_result$results$ICC[2]
  icc_2_k <- icc_result$results$ICC[5]
  
  cat("\n=== ICC 结果（psych 包）===\n")
  cat("  ICC(2,1) - 单个评分者一致性:", round(icc_2_1, 4), "\n")
  cat("  ICC(2,k) - 平均评分者一致性:", round(icc_2_k, 4), "\n")
  
} else {
  ## 使用简化 ICC 计算
  icc_simple <- calc_icc_simple(eta_CL_data[, -1])
  icc_2_1 <- icc_simple$icc_2_1
  icc_2_k <- icc_simple$icc_2_k
  
  cat("\n=== ICC 结果（简化计算）===\n")
  cat("  ICC(2,1) - 单个评分者一致性:", round(icc_2_1, 4), "\n")
  cat("  ICC(2,k) - 平均评分者一致性:", round(icc_2_k, 4), "\n")
  cat("  ⚠️ 注意：使用简化方法计算，安装 lme4 可获得更精确结果\n")
}

## 保存 ICC 结果供后续使用
icc_results <- list(
  icc_2_1 = icc_2_1,
  icc_2_k = icc_2_k,
  method = ifelse(HAS_PSYCH && HAS_LME4, "psych", "simplified")
)

## ############################################################################
## Part D：按 |η_true| 分层分析（仅 zang2021）
## ############################################################################

message("\n", paste(rep("-", 70), collapse = ""))
message("Part D：按 |η_true| 分层分析（仅 zang2021）")
message(paste(rep("-", 70), collapse = ""))

## 由于数据由 zang2021 生成，只有它有对应的 η_true

## 定义 |η_CL_true| 分组
eta_breaks <- c(0, 0.5, 1.0, 1.5, Inf)
eta_labels <- c("[0, 0.5)", "[0.5, 1.0)", "[1.0, 1.5)", "[1.5, +∞)")

combined_data <- combined_data %>%
  mutate(
    abs_eta_CL_true = abs(eta_CL_true),
    eta_CL_bin = cut(abs_eta_CL_true, breaks = eta_breaks, labels = eta_labels, right = FALSE)
  )

## 只对 zang2021 模型进行分层分析
accuracy_by_eta_bin_zang <- combined_data %>%
  group_by(eta_CL_bin) %>%
  summarise(
    n = n(),
    accuracy_ind = mean(zang_match_ind) * 100,
    accuracy_pop = mean(zang_match_pop) * 100,
    .groups = "drop"
  ) %>%
  mutate(delta = accuracy_ind - accuracy_pop)

cat("\n=== zang2021 按 |η_CL_true| 分层的准确率 ===\n")
print(accuracy_by_eta_bin_zang)

## ############################################################################
## 4.2 生成可视化
## ############################################################################

message("\n", paste(rep("-", 70), collapse = ""))
message("4.2 生成可视化")
message(paste(rep("-", 70), collapse = ""))

## ----------------------------------------------------------------------------
## 图 1：各模型准确率柱状图
## ----------------------------------------------------------------------------

message("\n--- 生成图 1：各模型准确率柱状图 ---")

accuracy_plot_data <- accuracy_by_model %>%
  pivot_longer(cols = c(accuracy_ind, accuracy_pop),
               names_to = "method", values_to = "accuracy") %>%
  mutate(
    method = ifelse(method == "accuracy_ind", "IND (η_ind)", "POP (η = 0)"),
    method = factor(method, levels = c("POP (η = 0)", "IND (η_ind)"))
  )

p1_accuracy <- ggplot(accuracy_plot_data, aes(x = model, y = accuracy, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", accuracy)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("POP (η = 0)" = "coral", "IND (η_ind)" = "steelblue")) +
  scale_y_continuous(limits = c(0, 70), breaks = seq(0, 70, 10)) +
  labs(
    title = "各模型判别准确率",
    subtitle = paste0("准确率 SD: POP=", round(accuracy_sd_pop, 2), 
                      "%, IND=", round(accuracy_sd_ind, 2), "%"),
    x = "模型",
    y = "准确率 (%)",
    fill = "方法"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "bottom"
  )

ggsave(file.path(plot_dir, "fig1_accuracy_by_model.png"),
       p1_accuracy, width = 8, height = 6, dpi = 300)

message("  ✓ fig1_accuracy_by_model.png 已保存")

## ----------------------------------------------------------------------------
## 图 2：跨模型一致性热图
## ----------------------------------------------------------------------------

message("\n--- 生成图 2：跨模型一致性热图 ---")

## 准备数据
heatmap_data_ind <- as.data.frame(pairwise_agree_ind) %>%
  mutate(model_row = rownames(.)) %>%
  pivot_longer(cols = -model_row, names_to = "model_col", values_to = "agreement") %>%
  mutate(method = "IND")

heatmap_data_pop <- as.data.frame(pairwise_agree_pop) %>%
  mutate(model_row = rownames(.)) %>%
  pivot_longer(cols = -model_row, names_to = "model_col", values_to = "agreement") %>%
  mutate(method = "POP")

heatmap_data <- bind_rows(heatmap_data_ind, heatmap_data_pop) %>%
  mutate(
    model_row = factor(model_row, levels = models),
    model_col = factor(model_col, levels = models),
    method = factor(method, levels = c("POP", "IND"))
  )

p2_heatmap <- ggplot(heatmap_data, aes(x = model_col, y = model_row, fill = agreement)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", agreement)), size = 4) +
  facet_wrap(~ method, ncol = 2) +
  scale_fill_gradient2(low = "white", mid = "lightblue", high = "steelblue",
                       midpoint = 50, limits = c(0, 100)) +
  labs(
    title = "跨模型两两一致率热图",
    subtitle = "数值 = 两个模型给出相同 s_pred 的比例 (%)",
    x = "模型",
    y = "模型",
    fill = "一致率 (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(plot_dir, "fig2_pairwise_agreement_heatmap.png"),
       p2_heatmap, width = 10, height = 5, dpi = 300)

message("  ✓ fig2_pairwise_agreement_heatmap.png 已保存")

## ----------------------------------------------------------------------------
## 图 3：η_CL_ind 跨模型散点图矩阵
## ----------------------------------------------------------------------------

message("\n--- 生成图 3：η_CL_ind 跨模型散点图矩阵 ---")

## 准备所有两两组合
pairs_list <- list(
  c("zang", "li"), c("zang", "yin"), c("zang", "sun"),
  c("li", "yin"), c("li", "sun"), c("yin", "sun")
)

scatter_plots <- list()

for (pair in pairs_list) {
  model1 <- pair[1]
  model2 <- pair[2]
  
  r_val <- cor_matrix[model1, model2]
  
  p <- ggplot(eta_CL_data, aes_string(x = model1, y = model2)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "steelblue", linewidth = 0.8) +
    labs(
      title = paste(model1, "vs", model2),
      subtitle = paste0("r = ", round(r_val, 3)),
      x = paste0("η_CL_ind (", model1, ")"),
      y = paste0("η_CL_ind (", model2, ")")
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      plot.subtitle = element_text(hjust = 0.5, color = "steelblue", size = 9)
    )
  
  scatter_plots[[paste(model1, model2, sep = "_")]] <- p
}

## 使用 gridExtra 或 patchwork 组合
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p3_scatter <- wrap_plots(scatter_plots, ncol = 3)
  
  ggsave(file.path(plot_dir, "fig3_eta_CL_ind_scatter_matrix.png"),
         p3_scatter, width = 12, height = 8, dpi = 300)
  message("  ✓ fig3_eta_CL_ind_scatter_matrix.png 已保存")
  
} else {
  ## 如果没有 patchwork，单独保存每个图
  message("  ⚠️ patchwork 包未安装，将单独保存散点图")
  for (name in names(scatter_plots)) {
    ggsave(file.path(plot_dir, paste0("fig3_scatter_eta_", name, ".png")),
           scatter_plots[[name]], width = 5, height = 5, dpi = 300)
  }
  message("  ✓ 散点图已单独保存")
}

## ----------------------------------------------------------------------------
## 图 4：按 |η_CL_true| 分层的准确率柱状图
## ----------------------------------------------------------------------------

message("\n--- 生成图 4：按 |η_CL_true| 分层的准确率 ---")

eta_bin_plot_data <- accuracy_by_eta_bin_zang %>%
  pivot_longer(cols = c(accuracy_ind, accuracy_pop),
               names_to = "method", values_to = "accuracy") %>%
  mutate(
    method = ifelse(method == "accuracy_ind", "IND", "POP"),
    method = factor(method, levels = c("POP", "IND"))
  )

p4_eta_bin <- ggplot(eta_bin_plot_data, aes(x = eta_CL_bin, y = accuracy, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", accuracy)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("POP" = "coral", "IND" = "steelblue")) +
  scale_y_continuous(limits = c(0, 80)) +
  labs(
    title = "准确率 vs |η_CL_true| 分层（zang2021 模型）",
    subtitle = "数据由 zang2021 模型生成，只有该模型有对应的 η_true",
    x = "|η_CL_true| 区间",
    y = "准确率 (%)",
    fill = "方法"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "bottom"
  )

ggsave(file.path(plot_dir, "fig4_accuracy_by_eta_bin_zang.png"),
       p4_eta_bin, width = 8, height = 6, dpi = 300)

message("  ✓ fig4_accuracy_by_eta_bin_zang.png 已保存")

## ----------------------------------------------------------------------------
## 图 5：跨模型一致性汇总柱状图
## ----------------------------------------------------------------------------

message("\n--- 生成图 5：跨模型一致性汇总 ---")

consistency_summary <- tibble(
  metric = c("完全一致率", "全对率", "多数正确率", "任一正确率"),
  IND = c(agreement_ind, 
          mean(combined_data$all_correct_ind) * 100,
          mean(combined_data$majority_correct_ind) * 100,
          mean(combined_data$any_correct_ind) * 100),
  POP = c(agreement_pop,
          mean(combined_data$all_correct_pop) * 100,
          mean(combined_data$majority_correct_pop) * 100,
          mean(combined_data$any_correct_pop) * 100)
) %>%
  pivot_longer(cols = c(IND, POP), names_to = "method", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = c("完全一致率", "全对率", "多数正确率", "任一正确率")),
    method = factor(method, levels = c("POP", "IND"))
  )

p5_consistency <- ggplot(consistency_summary, aes(x = metric, y = value, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", value)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("POP" = "coral", "IND" = "steelblue")) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    title = "跨模型一致性指标汇总",
    subtitle = paste0("Fleiss' κ: POP=", round(fleiss_pop$value, 3), 
                      ", IND=", round(fleiss_ind$value, 3)),
    x = "指标",
    y = "比例 (%)",
    fill = "方法"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

ggsave(file.path(plot_dir, "fig5_consistency_summary.png"),
       p5_consistency, width = 8, height = 6, dpi = 300)

message("  ✓ fig5_consistency_summary.png 已保存")

## ----------------------------------------------------------------------------
## 图 6：按场景分层的一致率
## ----------------------------------------------------------------------------

message("\n--- 生成图 6：按场景分层的一致率 ---")

agree_by_scenario_plot <- agreement_by_scenario %>%
  pivot_longer(cols = c(agree_ind, agree_pop),
               names_to = "method", values_to = "agreement") %>%
  mutate(
    method = ifelse(method == "agree_ind", "IND", "POP"),
    method = factor(method, levels = c("POP", "IND"))
  )

p6_scenario <- ggplot(agree_by_scenario_plot, aes(x = s_true, y = agreement, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(agreement, 0)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("POP" = "coral", "IND" = "steelblue")) +
  scale_y_continuous(limits = c(0, max(agree_by_scenario_plot$agreement) * 1.2)) +
  labs(
    title = "判别一致率按场景分层",
    subtitle = "4 个模型给出相同 s_pred 的比例",
    x = "��实场景",
    y = "一致率 (%)",
    fill = "方法"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(plot_dir, "fig6_agreement_by_scenario.png"),
       p6_scenario, width = 10, height = 6, dpi = 300)

message("  ✓ fig6_agreement_by_scenario.png 已保存")

## ############################################################################
## 4.3 保存输出文件
## ############################################################################

message("\n--- 4.3 保存输出文件 ---")

## 构建评估结果列表
evaluation_results <- list(
  
  ## Part A：单模型准确率
  accuracy_by_model = accuracy_by_model,
  accuracy_sd = list(ind = accuracy_sd_ind, pop = accuracy_sd_pop),
  mcnemar_results = mcnemar_results,
  accuracy_by_scenario = accuracy_by_scenario,
  
  ## Part B：跨模型一致性
  cross_model_consistency = list(
    agreement_ind = agreement_ind,
    agreement_pop = agreement_pop,
    fleiss_kappa_ind = fleiss_ind$value,
    fleiss_kappa_pop = fleiss_pop$value,
    fleiss_kappa_ind_p = fleiss_ind$p.value,
    fleiss_kappa_pop_p = fleiss_pop$p.value
  ),
  multilevel_accuracy = multilevel_accuracy,
  pairwise_agreement_ind = pairwise_agree_ind,
  pairwise_agreement_pop = pairwise_agree_pop,
  agreement_by_scenario = agreement_by_scenario,
  
  ## Part C：η_ind 相关性
  eta_correlation = eta_correlation,
  eta_cor_matrix = cor_matrix,
  icc_results = icc_results,
  
  ## Part D：按 |η_true| 分层
  accuracy_by_eta_bin_zang = accuracy_by_eta_bin_zang,
  
  ## 元数据
  metadata = list(
    created = Sys.time(),
    version = "v9.0",
    module = "Module 4 - Consistency Evaluation",
    n_mc = N_MC,
    n_scenarios = nrow(scenario_library)
  )
)

## 保存评估结果
saveRDS(evaluation_results, file.path(OUT_DIR_9, "evaluation_results.rds"))
message("  ✓ evaluation_results.rds 已保存")

## 保存 CSV 版本的关键表格
write_csv(accuracy_by_model, file.path(OUT_DIR_9, "accuracy_by_model.csv"))
write_csv(mcnemar_results, file.path(OUT_DIR_9, "mcnemar_results.csv"))
write_csv(multilevel_accuracy, file.path(OUT_DIR_9, "multilevel_accuracy.csv"))
write_csv(eta_correlation, file.path(OUT_DIR_9, "eta_correlation.csv"))
write_csv(accuracy_by_eta_bin_zang, file.path(OUT_DIR_9, "accuracy_by_eta_bin_zang.csv"))

message("  ✓ CSV 文件已保存")

## 保存更新后的 combined_data
saveRDS(combined_data, file.path(OUT_DIR_9, "combined_data_enriched.rds"))
message("  ✓ combined_data_enriched.rds 已保存")

## ############################################################################
## 4.4 生成汇总报告
## ############################################################################

message("\n", paste(rep("=", 60), collapse = ""))
message("               第九版汇总报告")
message(paste(rep("=", 60), collapse = ""))

cat("\n【核心参数】\n")
cat("  虚拟患者数:", N_MC, "\n")
cat("  场景数:", nrow(scenario_library), "\n")
cat("  真实模型: zang2021 (一室模型)\n")
cat("  分析模型: zang2021, li2018, yin2016, sun2021\n")

cat("\n【Part A：单模型准确率】\n")
for (i in 1:nrow(accuracy_by_model)) {
  cat(sprintf("  %s: IND=%.1f%%, POP=%.1f%%, Δ=%+.1f%%\n",
              accuracy_by_model$model[i],
              accuracy_by_model$accuracy_ind[i],
              accuracy_by_model$accuracy_pop[i],
              accuracy_by_model$delta[i]))
}
cat(sprintf("\n  准确率 SD: IND=%.2f%%, POP=%.2f%%\n", accuracy_sd_ind, accuracy_sd_pop))

cat("\n【Part B：跨模型一致性（第九版核心）】\n")
cat(sprintf("  完全一致率: IND=%.1f%%, POP=%.1f%%\n", agreement_ind, agreement_pop))
cat(sprintf("  Fleiss' κ:  IND=%.3f, POP=%.3f\n", fleiss_ind$value, fleiss_pop$value))
cat(sprintf("  全对率:     IND=%.1f%%, POP=%.1f%%\n", 
            multilevel_accuracy$ind_pct[1], multilevel_accuracy$pop_pct[1]))
cat(sprintf("  多数正确率: IND=%.1f%%, POP=%.1f%%\n",
            multilevel_accuracy$ind_pct[2], multilevel_accuracy$pop_pct[2]))

cat("\n【Part C：η_CL_ind 跨模型相关性】\n")
cat(sprintf("  ICC(2,1) = %.3f\n", icc_results$icc_2_1))
cat("  两两 Pearson r:\n")
for (i in 1:nrow(eta_correlation)) {
  cat(sprintf("    %s: r = %.3f\n", eta_correlation$model_pair[i], eta_correlation$pearson_r[i]))
}

cat("\n【结论】\n")

## 判断假设是否成立
hypothesis_check <- list(
  agreement = agreement_ind > agreement_pop,
  fleiss = fleiss_ind$value > fleiss_pop$value,
  all_correct = multilevel_accuracy$ind_pct[1] > multilevel_accuracy$pop_pct[1],
  accuracy_improvement = all(accuracy_by_model$delta >= 0)
)

if (hypothesis_check$agreement && hypothesis_check$fleiss) {
  cat("  ✅ 假设验证通过：MAPB 使不同模型对同一患者的判别结果趋于一致\n")
} else {
  cat("  ❌ 假设验证未通过\n")
}

if (hypothesis_check$accuracy_improvement) {
  cat("  ✅ IND 模式在所有模型上的准确率 ≥ POP 模式\n")
} else {
  cat("  ⚠️  部分模型 IND 准确率低于 POP\n")
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")

## ############################################################################
## 4.5 模块 4 完成
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("✅ 模块 4 全部完成（一致性评估）")
message(paste(rep("=", 70), collapse = ""))

message("\n📁 输出目录: ", OUT_DIR_9)
message("\n📄 输出文件：")
message("  数据：")
message("  - evaluation_results.rds        （评估结果汇总）")
message("  - combined_data_enriched.rds    （增强版合并数据）")
message("  - accuracy_by_model.csv         （各模型准确率）")
message("  - mcnemar_results.csv           （McNemar 检验）")
message("  - multilevel_accuracy.csv       （多层次正确率）")
message("  - eta_correlation.csv           （η 相关系数）")
message("\n  图表（plots/）：")
message("  - fig1_accuracy_by_model.png")
message("  - fig2_pairwise_agreement_heatmap.png")
message("  - fig3_eta_CL_ind_scatter_matrix.png")
message("  - fig4_accuracy_by_eta_bin_zang.png")
message("  - fig5_consistency_summary.png")
message("  - fig6_agreement_by_scenario.png")

message("\n📊 核心结果：")
message(sprintf("   IND 完全一致率: %.1f%%", agreement_ind))
message(sprintf("   POP 完全一致率: %.1f%%", agreement_pop))
message(sprintf("   Fleiss' κ: IND=%.3f, POP=%.3f", fleiss_ind$value, fleiss_pop$value))

message("\n🔗 后续使用方法：")
message('  results <- readRDS("outputs9/evaluation_results.rds")')
message('  accuracy <- results$accuracy_by_model')
message('  consistency <- results$cross_model_consistency')

message("\n✅ 第九版代码全部完成！")


## ############################################################################
## 补充图：两两一致且都正确的热图
## ############################################################################

message("\n--- 生成补充图：两两一致且都正确热图 ---")

models <- c("zang", "li", "yin", "sun")
n_models <- length(models)

## IND 模式：两两一致且都正确
pairwise_both_correct_ind <- matrix(NA, n_models, n_models, 
                                     dimnames = list(models, models))

for (i in 1:n_models) {
  for (j in 1:n_models) {
    match_i <- paste0(models[i], "_match_ind")
    match_j <- paste0(models[j], "_match_ind")
    
    ## 两个模型都正确（则它们必然一致且都正确）
    pairwise_both_correct_ind[i, j] <- mean(
      combined_data[[match_i]] & combined_data[[match_j]]
    ) * 100
  }
}

## POP 模式：两两一致且都正确
pairwise_both_correct_pop <- matrix(NA, n_models, n_models,
                                     dimnames = list(models, models))

for (i in 1:n_models) {
  for (j in 1:n_models) {
    match_i <- paste0(models[i], "_match_pop")
    match_j <- paste0(models[j], "_match_pop")
    
    pairwise_both_correct_pop[i, j] <- mean(
      combined_data[[match_i]] & combined_data[[match_j]]
    ) * 100
  }
}

cat("\n=== 两两一致且都正确矩阵 (IND) ===\n")
print(round(pairwise_both_correct_ind, 1))

cat("\n=== 两两一致且都正确矩阵 (POP) ===\n")
print(round(pairwise_both_correct_pop, 1))

## 准备绘图数据
heatmap_correct_ind <- as.data.frame(pairwise_both_correct_ind) %>%
  mutate(model_row = rownames(.)) %>%
  pivot_longer(cols = -model_row, names_to = "model_col", values_to = "both_correct") %>%
  mutate(method = "IND")

heatmap_correct_pop <- as.data.frame(pairwise_both_correct_pop) %>%
  mutate(model_row = rownames(.)) %>%
  pivot_longer(cols = -model_row, names_to = "model_col", values_to = "both_correct") %>%
  mutate(method = "POP")

heatmap_correct_data <- bind_rows(heatmap_correct_ind, heatmap_correct_pop) %>%
  mutate(
    model_row = factor(model_row, levels = models),
    model_col = factor(model_col, levels = models),
    method = factor(method, levels = c("POP", "IND"))
  )

## 绘制热图
p_both_correct_heatmap <- ggplot(heatmap_correct_data, 
                                  aes(x = model_col, y = model_row, fill = both_correct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", both_correct)), size = 4) +
  facet_wrap(~ method, ncol = 2) +
  scale_fill_gradient2(low = "white", mid = "lightgreen", high = "darkgreen",
                       midpoint = 25, limits = c(0, 60)) +
  labs(
    title = "跨模型两两一致且都正确热图",
    subtitle = "数值 = 两个模型都正确的比例 (%)",
    x = "模型",
    y = "模型",
    fill = "比例 (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(plot_dir, "fig7_pairwise_both_correct_heatmap.png"),
       p_both_correct_heatmap, width = 10, height = 5, dpi = 300)

message("  ✓ fig7_pairwise_both_correct_heatmap.png 已保存")

## 与原热图对比
cat("\n=== 对比：两两一致 vs 两两都正确 (IND) ===\n")
cat("两两一致（对角线外平均）:", 
    round(mean(pairwise_agree_ind[upper.tri(pairwise_agree_ind)]), 1), "%\n")
cat("两两都正确（对角线外平均）:", 
    round(mean(pairwise_both_correct_ind[upper.tri(pairwise_both_correct_ind)]), 1), "%\n")
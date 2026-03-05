## ############################################################################
## 第九版代码-敏感性分析-模块4.R
## 
## 核心任务：
##   参数化的一致性评估模块
##   读取 SA_TRUE_MODEL 和 OUT_DIR_SA 变量
## 
## 设计说明：
##   - 此脚本由主控脚本调用，不单独运行
##   - 需要预先设置全局变量：SA_TRUE_MODEL, OUT_DIR_SA
##   - 评估逻辑与原模块 4 相同，但简化输出（无图表）
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("敏感性分析 模块4：一致性评估（真实模型 = ", SA_TRUE_MODEL, "）")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## SA4.0 加载依赖与读取前置模块输出
## ############################################################################

message("\n--- SA4.0 加载依赖 ---")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(readr)
  library(irr)      # Fleiss' Kappa
})

## 读取敏感性分析前置模块输出
discrimination_results <- readRDS(file.path(OUT_DIR_SA, "discrimination_results.rds"))
mapb_results <- readRDS(file.path(OUT_DIR_SA, "mapb_results.rds"))
scenario_library <- readRDS(file.path(OUT_DIR_SA, "scenario_library.rds"))
module1_params <- readRDS(file.path(OUT_DIR_SA, "module1_params.rds"))

N_MC <- module1_params$N_MC

message("✅ 前置模块输出已加载")
message("   虚��患者数: ", N_MC)
message("   真实模型: ", SA_TRUE_MODEL)

## ############################################################################
## SA4.1 合并数据
## ############################################################################

message("\n--- SA4.1 合并数据 ---")

combined_data <- discrimination_results %>%
  left_join(
    mapb_results %>% select(mc_iter, 
                            zang_eta_CL_ind, zang_eta_V_ind,
                            li_eta_CL_ind, li_eta_Vc_ind,
                            yin_eta_CL_ind, sun_eta_CL_ind),
    by = "mc_iter"
  )

message("✅ 数据合并完成")

## ############################################################################
## SA4.2 计算单模型准确率
## ############################################################################

message("\n--- SA4.2 计算单模型准确率 ---")

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
  mutate(delta = accuracy_ind - accuracy_pop)

cat("\n=== 各模型准确率 ===\n")
print(accuracy_by_model)

## ############################################################################
## SA4.3 计算跨模型一致性指标
## ############################################################################

message("\n--- SA4.3 计算跨模型一致性指标 ---")

## 完全一致率（已在模块 3 计算）
agreement_ind <- mean(combined_data$all_agree_ind) * 100
agreement_pop <- mean(combined_data$all_agree_pop) * 100

cat("\n=== 完全一致率 ===\n")
cat("  IND 模式:", round(agreement_ind, 2), "%\n")
cat("  POP 模式:", round(agreement_pop, 2), "%\n")

## Fleiss' Kappa
ratings_ind <- combined_data %>%
  select(mc_iter, zang_s_pred_ind, li_s_pred_ind, yin_s_pred_ind, sun_s_pred_ind) %>%
  mutate(across(-mc_iter, ~ factor(.x, levels = scenario_library$scenario_id)))

ratings_matrix_ind <- as.matrix(ratings_ind[, -1])
fleiss_ind <- kappam.fleiss(ratings_matrix_ind)

ratings_pop <- combined_data %>%
  select(mc_iter, zang_s_pred_pop, li_s_pred_pop, yin_s_pred_pop, sun_s_pred_pop) %>%
  mutate(across(-mc_iter, ~ factor(.x, levels = scenario_library$scenario_id)))

ratings_matrix_pop <- as.matrix(ratings_pop[, -1])
fleiss_pop <- kappam.fleiss(ratings_matrix_pop)

cat("\n=== Fleiss' Kappa ===\n")
cat("  IND 模式: κ =", round(fleiss_ind$value, 4), "\n")
cat("  POP 模式: κ =", round(fleiss_pop$value, 4), "\n")

## ############################################################################
## SA4.4 计算多层次正确率
## ############################################################################

message("\n--- SA4.4 计算多层次正确率 ---")

combined_data <- combined_data %>%
  mutate(
    n_correct_ind = zang_match_ind + li_match_ind + yin_match_ind + sun_match_ind,
    majority_correct_ind = (n_correct_ind >= 3),
    any_correct_ind = (n_correct_ind >= 1),
    
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

## ############################################################################
## SA4.5 构建评估结果
## ############################################################################

message("\n--- SA4.5 构建评估结果 ---")

evaluation_results <- list(
  ## 真实模型
  true_model = SA_TRUE_MODEL,
  
  ## 单模型准确率
  accuracy_by_model = accuracy_by_model,
  
  ## 跨模型一致性
  cross_model_consistency = list(
    agreement_ind = agreement_ind,
    agreement_pop = agreement_pop,
    fleiss_kappa_ind = fleiss_ind$value,
    fleiss_kappa_pop = fleiss_pop$value
  ),
  
  ## 多层次正确率
  multilevel_accuracy = multilevel_accuracy,
  
  ## 元数据
  metadata = list(
    created = Sys.time(),
    version = "v9.0-SA",
    n_mc = N_MC
  )
)

## ############################################################################
## SA4.6 保存输出
## ############################################################################

message("\n--- SA4.6 保存输出 ---")

saveRDS(evaluation_results, file.path(OUT_DIR_SA, "evaluation_results.rds"))
message("  ✓ evaluation_results.rds 已保存")

## ############################################################################
## SA4.7 打印汇总
## ############################################################################

cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("敏感性分析汇总（真实模型 = ", SA_TRUE_MODEL, "）\n", sep = "")
cat(paste(rep("=", 50), collapse = ""), "\n")

cat("\n【单模型准确率】\n")
for (i in 1:nrow(accuracy_by_model)) {
  cat(sprintf("  %s: IND=%.1f%%, POP=%.1f%%, Δ=%+.1f%%\n",
              accuracy_by_model$model[i],
              accuracy_by_model$accuracy_ind[i],
              accuracy_by_model$accuracy_pop[i],
              accuracy_by_model$delta[i]))
}

cat("\n【跨模型一致性】\n")
cat(sprintf("  完全一致率: IND=%.1f%%, POP=%.1f%%\n", agreement_ind, agreement_pop))
cat(sprintf("  Fleiss' κ:  IND=%.3f, POP=%.3f\n", fleiss_ind$value, fleiss_pop$value))
cat(sprintf("  全对率:     IND=%.1f%%, POP=%.1f%%\n", 
            multilevel_accuracy$ind_pct[1], multilevel_accuracy$pop_pct[1]))

message("\n✅ 敏感性分析模块 4 完成（真实模型 = ", SA_TRUE_MODEL, "）")
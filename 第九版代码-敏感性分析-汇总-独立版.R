## ############################################################################
## 第九版代码-敏感性分析-汇总-独立版.R
## 
## 特点：
##   - 完全独立运行，不依赖主控脚本
##   - 只读取已保存的 evaluation_results.rds 文件
##   - 包含完备的多层级判定逻辑
## 
## 输入：
##   - outputs9/evaluation_results.rds（基准分析）
##   - outputs9_SA/{model}/evaluation_results.rds（敏感性分析）
## 
## 输出：
##   - outputs9_SA/sensitivity_summary.rds / .csv
##   - outputs9_SA/sensitivity_report.txt
##   - outputs9_SA/sensitivity_judgment.rds（多层级判定结果）
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("敏感性分析汇总（独立版）")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## 0. 加载依赖与路径设置
## ############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})

## 路径设置
PROJ_ROOT <- "D:/Users/YujiaZhang/Desktop/26救赎之道/25.9.29 毕业设计/5正式实施！/多成人模型/体现我卓越的项目管理能力"
OUT_DIR_9 <- file.path(PROJ_ROOT, "outputs9")
OUT_ROOT_SA <- file.path(PROJ_ROOT, "outputs9_SA")

## ############################################################################
## 1. 读取所有分析结果
## ############################################################################

message("\n--- 1. 读取所有分析结果 ---")

## 定义所有分析（基准 + 敏感性）
ALL_MODELS <- c("zang2021", "li2018", "yin2016", "sun2021")

## 初始化汇总表
summary_list <- list()

## 读取基准分析（zang2021）
baseline_path <- file.path(OUT_DIR_9, "evaluation_results.rds")
if (file.exists(baseline_path)) {
  baseline_results <- readRDS(baseline_path)
  
  summary_list[["zang2021"]] <- tibble(
    true_model = "zang2021",
    analysis_type = "基准",
    agreement_ind = baseline_results$cross_model_consistency$agreement_ind,
    agreement_pop = baseline_results$cross_model_consistency$agreement_pop,
    fleiss_kappa_ind = baseline_results$cross_model_consistency$fleiss_kappa_ind,
    fleiss_kappa_pop = baseline_results$cross_model_consistency$fleiss_kappa_pop,
    all_correct_ind = baseline_results$multilevel_accuracy$ind_pct[1],
    all_correct_pop = baseline_results$multilevel_accuracy$pop_pct[1],
    mean_accuracy_ind = mean(baseline_results$accuracy_by_model$accuracy_ind),
    mean_accuracy_pop = mean(baseline_results$accuracy_by_model$accuracy_pop)
  )
  message("  ✓ zang2021（基准）已读取")
} else {
  message("  ⚠️ zang2021（基准）未找到")
}

## 读取敏感性分析
SENSITIVITY_MODELS <- c("li2018", "yin2016", "sun2021")

for (model_name in SENSITIVITY_MODELS) {
  result_path <- file.path(OUT_ROOT_SA, model_name, "evaluation_results.rds")
  
  if (file.exists(result_path)) {
    sa_result <- readRDS(result_path)
    
    summary_list[[model_name]] <- tibble(
      true_model = model_name,
      analysis_type = "敏感性分析",
      agreement_ind = sa_result$cross_model_consistency$agreement_ind,
      agreement_pop = sa_result$cross_model_consistency$agreement_pop,
      fleiss_kappa_ind = sa_result$cross_model_consistency$fleiss_kappa_ind,
      fleiss_kappa_pop = sa_result$cross_model_consistency$fleiss_kappa_pop,
      all_correct_ind = sa_result$multilevel_accuracy$ind_pct[1],
      all_correct_pop = sa_result$multilevel_accuracy$pop_pct[1],
      mean_accuracy_ind = mean(sa_result$accuracy_by_model$accuracy_ind),
      mean_accuracy_pop = mean(sa_result$accuracy_by_model$accuracy_pop)
    )
    message("  ✓ ", model_name, " 已读取")
  } else {
    message("  ⚠️ ", model_name, " 未找到: ", result_path)
  }
}

## 合并汇总表
sensitivity_summary <- bind_rows(summary_list) %>%
  mutate(
    delta_kappa = fleiss_kappa_ind - fleiss_kappa_pop,
    delta_agreement = agreement_ind - agreement_pop,
    delta_all_correct = all_correct_ind - all_correct_pop,
    delta_mean_accuracy = mean_accuracy_ind - mean_accuracy_pop,
    relative_kappa_improvement = delta_kappa / fleiss_kappa_pop * 100
  )

message("\n✅ 已读取 ", nrow(sensitivity_summary), " 个分析结果")

## ############################################################################
## 2. 打印汇总表
## ############################################################################

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("                     敏感性分析汇总表\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

cat("\n【核心指标对比】\n\n")

print_table <- sensitivity_summary %>%
  select(
    `真实模型` = true_model,
    `Fleiss κ (IND)` = fleiss_kappa_ind,
    `Fleiss κ (POP)` = fleiss_kappa_pop,
    `Δκ` = delta_kappa,
    `完全一致率 IND` = agreement_ind,
    `完全一致率 POP` = agreement_pop
  ) %>%
  mutate(
    `Fleiss κ (IND)` = round(`Fleiss κ (IND)`, 3),
    `Fleiss κ (POP)` = round(`Fleiss κ (POP)`, 3),
    `Δκ` = sprintf("%+.3f", `Δκ`),
    `完全一致率 IND` = sprintf("%.1f%%", `完全一致率 IND`),
    `完全一致率 POP` = sprintf("%.1f%%", `完全一致率 POP`)
  )

print(print_table, n = Inf)

## ############################################################################
## 3. 多层级判定（核心！）
## ############################################################################

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("结论判定（多层级框架）\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

## 提取指标向量
delta_kappa_vec <- sensitivity_summary$delta_kappa
kappa_ind_vec <- sensitivity_summary$fleiss_kappa_ind
kappa_pop_vec <- sensitivity_summary$fleiss_kappa_pop

## ============================================================================
## Level 1: 方向性（所有 Δκ > 0）
## ============================================================================
level1_pass <- all(delta_kappa_vec > 0)

cat("【Level 1: 方向性】\n")
cat(sprintf("  标准: 所有 Δκ > 0\n"))
cat(sprintf("  判定: %s\n", ifelse(level1_pass, "✅ 通过", "❌ 未通过")))
for (i in 1:nrow(sensitivity_summary)) {
  cat(sprintf("    %s: Δκ = %+.3f %s\n",
              sensitivity_summary$true_model[i],
              delta_kappa_vec[i],
              ifelse(delta_kappa_vec[i] > 0, "✓", "✗")))
}

## ============================================================================
## Level 2: 效应量阈值（所有 Δκ ≥ 0.05 或相对提升 ≥ 20%）
## ============================================================================
criterion_absolute <- all(delta_kappa_vec >= 0.05)
criterion_relative <- all((delta_kappa_vec / kappa_pop_vec) >= 0.20)
level2_pass <- criterion_absolute | criterion_relative

cat("\n【Level 2: 效应量】\n")
cat(sprintf("  标准 A: 所有 Δκ ≥ 0.05: %s\n", ifelse(criterion_absolute, "✅ 通过", "❌ 未通过")))
cat(sprintf("  标准 B: 所有相对提升 ≥ 20%%: %s\n", ifelse(criterion_relative, "✅ 通过", "❌ 未通过")))
cat(sprintf("  判定（A 或 B）: %s\n", ifelse(level2_pass, "✅ 通过", "⚠️ 效应量不足")))

for (i in 1:nrow(sensitivity_summary)) {
  rel_imp <- sensitivity_summary$relative_kappa_improvement[i]
  cat(sprintf("    %s: Δκ = %.3f, 相对提升 = %+.1f%%\n",
              sensitivity_summary$true_model[i],
              delta_kappa_vec[i],
              rel_imp))
}

## ============================================================================
## Level 3: 一致性强度（κ_IND 至少达到 Fair 或 Moderate）
## ============================================================================
kappa_ind_min <- min(kappa_ind_vec)
kappa_ind_mean <- mean(kappa_ind_vec)

interpret_kappa <- function(k) {
  if (k < 0.20) return("Poor")
  if (k < 0.40) return("Fair")
  if (k < 0.60) return("Moderate")
  if (k < 0.80) return("Good")
  return("Excellent")
}

level3_pass <- (kappa_ind_min >= 0.20)  # 至少 Fair

cat("\n【Level 3: 一致性强度】\n")
cat(sprintf("  标准: min(κ_IND) ≥ 0.20 (至少 Fair)\n"))
cat(sprintf("  最小 κ (IND): %.3f [%s]\n", kappa_ind_min, interpret_kappa(kappa_ind_min)))
cat(sprintf("  平均 κ (IND): %.3f [%s]\n", kappa_ind_mean, interpret_kappa(kappa_ind_mean)))
cat(sprintf("  判定: %s\n", 
            ifelse(level3_pass, "✅ 通过", "⚠️ κ (IND) 仍在 Poor 范围")))

for (i in 1:nrow(sensitivity_summary)) {
  cat(sprintf("    %s: κ_IND = %.3f [%s]\n",
              sensitivity_summary$true_model[i],
              kappa_ind_vec[i],
              interpret_kappa(kappa_ind_vec[i])))
}

## ============================================================================
## Level 4: 稳健性（变异系数 CV < 50%）
## ============================================================================
kappa_cv_ind <- sd(kappa_ind_vec) / mean(kappa_ind_vec) * 100
delta_kappa_mean <- mean(delta_kappa_vec)
delta_kappa_cv <- ifelse(delta_kappa_mean > 0, 
                         sd(delta_kappa_vec) / delta_kappa_mean * 100, 
                         NA)

level4_pass <- (kappa_cv_ind < 50) & (!is.na(delta_kappa_cv) && delta_kappa_cv < 50)

cat("\n【Level 4: 稳健性】\n")
cat(sprintf("  标准: CV < 50%%\n"))
cat(sprintf("  κ (IND) 变异系数: %.1f%%\n", kappa_cv_ind))
cat(sprintf("  Δκ 变异系数: %.1f%%\n", ifelse(is.na(delta_kappa_cv), NA, delta_kappa_cv)))
cat(sprintf("  判定: %s\n", 
            ifelse(level4_pass, "✅ 通过（变异可接受）", "⚠️ 变异较大")))

## ============================================================================
## 综合判定
## ============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("综合判定\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("┌─────────────────────────────────────────────────────────────┐\n")
cat("│                    敏感性分析判定结果                        │\n")
cat("├─────────────────────────────────────────────────────────────┤\n")
cat(sprintf("│  Level 1 (方向性):     %-36s │\n", ifelse(level1_pass, "✅ 通过", "❌ 未通过")))
cat(sprintf("│  Level 2 (效应量):     %-36s │\n", ifelse(level2_pass, "✅ 通过", "⚠️ 不足")))
cat(sprintf("│  Level 3 (一致性强度): %-36s │\n", ifelse(level3_pass, "✅ 通过", "⚠️ 有限")))
cat(sprintf("│  Level 4 (稳健性):     %-36s │\n", ifelse(level4_pass, "✅ 通过", "⚠️ 变异大")))
cat("├─────────────────────────────────────────────────────────────┤\n")

## 总体结论
if (level1_pass & level2_pass & level3_pass & level4_pass) {
  conclusion <- "MAPB 使不同模型的判别结果趋于一致，此结论稳健且不依赖特定真实模型"
  conclusion_level <- "强结论"
} else if (level1_pass & level2_pass) {
  conclusion <- "MAPB 提升跨模型一致性的效果在所有真实模型下成立，但部分条件未完全满足"
  conclusion_level <- "中等结论"
} else if (level1_pass) {
  conclusion <- "MAPB 在所有真实模型下使 κ 提升（Δκ > 0），但效应量或一致性强度有限"
  conclusion_level <- "弱结论"
} else {
  conclusion <- "敏感性分析不支持原结论，存在 Δκ ≤ 0 的情况"
  conclusion_level <- "无结论"
}

cat(sprintf("│  总体结论: %-48s │\n", conclusion_level))
cat("└─────────────────────────────────────────────────────────────┘\n")

cat("\n【结论陈述】\n")
cat(conclusion, "\n")

## ############################################################################
## 4. 构建判定结果对象
## ############################################################################

judgment_results <- list(
  ## Level 1
  level1 = list(
    name = "方向性",
    criterion = "所有 Δκ > 0",
    pass = level1_pass,
    details = tibble(
      true_model = sensitivity_summary$true_model,
      delta_kappa = delta_kappa_vec,
      pass = delta_kappa_vec > 0
    )
  ),
  
  ## Level 2
  level2 = list(
    name = "效应量",
    criterion = "所有 Δκ ≥ 0.05 或 相对提升 ≥ 20%",
    criterion_absolute_pass = criterion_absolute,
    criterion_relative_pass = criterion_relative,
    pass = level2_pass,
    details = tibble(
      true_model = sensitivity_summary$true_model,
      delta_kappa = delta_kappa_vec,
      relative_improvement_pct = sensitivity_summary$relative_kappa_improvement
    )
  ),
  
  ## Level 3
  level3 = list(
    name = "一致性强度",
    criterion = "min(κ_IND) ≥ 0.20 (Fair)",
    kappa_ind_min = kappa_ind_min,
    kappa_ind_mean = kappa_ind_mean,
    interpretation_min = interpret_kappa(kappa_ind_min),
    interpretation_mean = interpret_kappa(kappa_ind_mean),
    pass = level3_pass
  ),
  
  ## Level 4
  level4 = list(
    name = "稳健性",
    criterion = "CV < 50%",
    kappa_cv_ind = kappa_cv_ind,
    delta_kappa_cv = delta_kappa_cv,
    pass = level4_pass
  ),
  
  ## 综合
  overall = list(
    conclusion = conclusion,
    conclusion_level = conclusion_level,
    all_levels_pass = level1_pass & level2_pass & level3_pass & level4_pass
  ),
  
  ## 元数据
  metadata = list(
    created = Sys.time(),
    n_analyses = nrow(sensitivity_summary),
    models_analyzed = sensitivity_summary$true_model
  )
)

## ############################################################################
## 4.5 CL 一致性结果汇总（新增）
## ############################################################################

message("\n--- 4.5 CL 一致性结果汇总 ---")

## 初始化 CL 一致性汇总列表
CL_consistency_list <- list()

## 读取基准分析的 CL 一致性
if (file.exists(baseline_path) && !is.null(baseline_results$CL_consistency)) {
  
  cl_cons <- baseline_results$CL_consistency
  
  CL_consistency_list[["zang2021"]] <- tibble(
    true_model = "zang2021",
    analysis_type = "基准",
    CV_pop = cl_cons$CV_comparison$CV_pop,
    CV_ind_median = cl_cons$CV_comparison$CV_ind_summary["median"],
    pct_CV_reduced = cl_cons$CV_comparison$pct_CV_reduced,
    ICC_2_1 = cl_cons$CL_icc$icc_2_1,
    ICC_2_k = cl_cons$CL_icc$icc_2_k,
    ICC_level = cl_cons$CL_icc$level,
    n_pairs_reduced = cl_cons$n_pairs_reduced,
    n_pairs_sig = cl_cons$n_pairs_sig,
    mean_CCC = mean(cl_cons$ccc_matrix[upper.tri(cl_cons$ccc_matrix)], na.rm = TRUE)
  )
  message("  ✓ zang2021 CL 一致性已读取")
} else {
  message("  ⚠️ zang2021 CL 一致性数据未找到")
}

## 读取敏感性分析的 CL 一致性
for (model_name in SENSITIVITY_MODELS) {
  
  result_path <- file.path(OUT_ROOT_SA, model_name, "evaluation_results.rds")
  
  if (file.exists(result_path)) {
    
    sa_result <- readRDS(result_path)
    
    if (!is.null(sa_result$CL_consistency)) {
      
      cl_cons <- sa_result$CL_consistency
      
      CL_consistency_list[[model_name]] <- tibble(
        true_model = model_name,
        analysis_type = "敏感性分析",
        CV_pop = cl_cons$CV_comparison$CV_pop,
        CV_ind_median = cl_cons$CV_comparison$CV_ind_summary["median"],
        pct_CV_reduced = cl_cons$CV_comparison$pct_CV_reduced,
        ICC_2_1 = cl_cons$CL_icc$icc_2_1,
        ICC_2_k = cl_cons$CL_icc$icc_2_k,
        ICC_level = cl_cons$CL_icc$level,
        n_pairs_reduced = cl_cons$n_pairs_reduced,
        n_pairs_sig = cl_cons$n_pairs_sig,
        mean_CCC = mean(cl_cons$ccc_matrix[upper.tri(cl_cons$ccc_matrix)], na.rm = TRUE)
      )
      message("  ✓ ", model_name, " CL 一致性已读取")
    } else {
      message("  ⚠️ ", model_name, " CL 一致性数据未找到")
    }
  }
}

## 合并 CL 一致性汇总
if (length(CL_consistency_list) > 0) {
  
  CL_consistency_summary <- bind_rows(CL_consistency_list) %>%
    mutate(
      ## 计算 CV 缩小率
      CV_reduction_pct = (CV_pop - CV_ind_median) / CV_pop * 100,
      ## 计算与基准的 ICC 差异
      delta_ICC = ICC_2_1 - ICC_2_1[true_model == "zang2021"]
    )
  
  ## ========================================================================
  ## 打印 CL 一致性汇总表
  ## ========================================================================
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("                 CL 跨模型一致性汇总（Part E）\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  cat("\n【CV 对比：跨模型 CL 离散度】\n")
  cat("  CV_pop = 各模型 CL_pop 的变异系数（固定）\n")
  cat("  CV_ind = 同一患者经 MAPB 后，各模型 CL_ind 的变异系数（中位数）\n\n")
  
  for (i in 1:nrow(CL_consistency_summary)) {
    row <- CL_consistency_summary[i, ]
    cat(sprintf("  %s: CV_pop=%.1f%%, CV_ind=%.1f%%, 缩小=%.1f%%, %.1f%%患者CV缩小\n",
                row$true_model,
                row$CV_pop,
                row$CV_ind_median,
                row$CV_reduction_pct,
                row$pct_CV_reduced))
  }
  
  cat("\n【ICC：CL_ind 跨模型一致性】\n")
  for (i in 1:nrow(CL_consistency_summary)) {
    row <- CL_consistency_summary[i, ]
    cat(sprintf("  %s: ICC(2,1)=%.4f [%s], ICC(2,k)=%.4f\n",
                row$true_model,
                row$ICC_2_1,
                row$ICC_level,
                row$ICC_2_k))
  }
  
  cat("\n【配对缩小率】\n")
  cat("  n_pairs_reduced = ratio < 1 的模型对数（共 6 对）\n")
  cat("  n_pairs_sig = 显著缩小（Wilcoxon p < 0.05）的模型对数\n\n")
  
  for (i in 1:nrow(CL_consistency_summary)) {
    row <- CL_consistency_summary[i, ]
    cat(sprintf("  %s: %d/6 对缩小, %d/6 对显著\n",
                row$true_model,
                row$n_pairs_reduced,
                row$n_pairs_sig))
  }
  
  cat("\n【CCC：两两 CL_ind 一致性（平均）】\n")
  for (i in 1:nrow(CL_consistency_summary)) {
    row <- CL_consistency_summary[i, ]
    cat(sprintf("  %s: mean CCC = %.4f\n",
                row$true_model,
                row$mean_CCC))
  }
  
  ## ========================================================================
  ## CL 一致性多层级判定
  ## ========================================================================
  
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("CL 一致性结论判定\n")
  cat(paste(rep("-", 60), collapse = ""), "\n\n")
  
  ## 判定标准：
  ## 1. 所有真实模型下，CV_ind < CV_pop（即 CV 被缩小）
  ## 2. 所有真实模型下，>90% 患者的 CV 被缩小
  ## 3. 所有真实模型下，至少 5/6 对模型的 CL 差值被缩小
  
  cl_level1_pass <- all(CL_consistency_summary$CV_ind_median < CL_consistency_summary$CV_pop)
  cl_level2_pass <- all(CL_consistency_summary$pct_CV_reduced > 90)
  cl_level3_pass <- all(CL_consistency_summary$n_pairs_reduced >= 5)
  
  cat("【CL-Level 1: CV 缩小方向】\n")
  cat(sprintf("  标准: 所有真实模型下 CV_ind < CV_pop\n"))
  cat(sprintf("  判定: %s\n\n", ifelse(cl_level1_pass, "✅ 通过", "❌ 未通过")))
  
  cat("【CL-Level 2: CV 缩小覆盖率】\n")
  cat(sprintf("  标准: 所有真实模型下 >90%% 患者 CV 被缩小\n"))
  cat(sprintf("  判定: %s\n", ifelse(cl_level2_pass, "✅ 通过", "❌ 未通过")))
  for (i in 1:nrow(CL_consistency_summary)) {
    row <- CL_consistency_summary[i, ]
    cat(sprintf("    %s: %.1f%%\n", row$true_model, row$pct_CV_reduced))
  }
  
  cat("\n【CL-Level 3: 配对缩小稳健性】\n")
  cat(sprintf("  标准: 所有真实模型下 ≥5/6 模型对 CL 差值缩小\n"))
  cat(sprintf("  判定: %s\n", ifelse(cl_level3_pass, "✅ 通过", "❌ 未通过")))
  for (i in 1:nrow(CL_consistency_summary)) {
    row <- CL_consistency_summary[i, ]
    cat(sprintf("    %s: %d/6\n", row$true_model, row$n_pairs_reduced))
  }
  
  ## 综合判定
  cl_all_pass <- cl_level1_pass & cl_level2_pass & cl_level3_pass
  
  cat("\n┌─────────────────────────────────────────────────────────────┐\n")
  cat("│                    CL 一致性综合结论                         │\n")
  cat("├─────────────────────────────────────────────────────────────┤\n")
  
  if (cl_all_pass) {
    cl_conclusion <- "MAPB 稳健地提升了跨模型 CL 估计的一致性（跨真实模型验证通过）"
    cl_conclusion_level <- "强证据"
  } else if (cl_level1_pass & cl_level2_pass) {
    cl_conclusion <- "MAPB 有效缩小了跨模型 CL 离散度（部分验证通过）"
    cl_conclusion_level <- "中等证据"
  } else if (cl_level1_pass) {
    cl_conclusion <- "MAPB 在方向上改善了 CL 一致性，但效果不够稳健"
    cl_conclusion_level <- "弱证据"
  } else {
    cl_conclusion <- "CL 一致性改善证据不充分"
    cl_conclusion_level <- "无证据"
  }
  
  cat(sprintf("│  结论强度: %-48s │\n", cl_conclusion_level))
  cat("└─────────────────────────────────────────────────────────────┘\n")
  cat("\n【CL 一致性结论】\n")
  cat(cl_conclusion, "\n")
  
  ## ========================================================================
  ## 将 CL 一致性结果添加到判定对象
  ## ========================================================================
  
  judgment_results$CL_consistency <- list(
    summary = CL_consistency_summary,
    levels = list(
      level1 = list(name = "CV缩小方向", pass = cl_level1_pass),
      level2 = list(name = "CV缩小覆盖率", pass = cl_level2_pass),
      level3 = list(name = "配对缩小稳健性", pass = cl_level3_pass)
    ),
    conclusion = cl_conclusion,
    conclusion_level = cl_conclusion_level,
    all_pass = cl_all_pass
  )
  
} else {
  message("  ⚠️ 未找到任何 CL 一致性数据，跳过汇总")
  CL_consistency_summary <- NULL
}

## ############################################################################
## 5. 保存输出
## ############################################################################

message("\n--- 5. 保存输出 ---")

## 创建输出目录（如果不存在）
if (!dir.exists(OUT_ROOT_SA)) {
  dir.create(OUT_ROOT_SA, recursive = TRUE)
}

## 保存汇总表
saveRDS(sensitivity_summary, file.path(OUT_ROOT_SA, "sensitivity_summary.rds"))
write_csv(sensitivity_summary, file.path(OUT_ROOT_SA, "sensitivity_summary.csv"))
message("  ✓ sensitivity_summary.rds / .csv 已保存")

## 保存 CL 一致性汇总（新增��
if (!is.null(CL_consistency_summary)) {
  saveRDS(CL_consistency_summary, file.path(OUT_ROOT_SA, "CL_consistency_summary.rds"))
  write_csv(CL_consistency_summary, file.path(OUT_ROOT_SA, "CL_consistency_across_models.csv"))
  message("  ✓ CL_consistency_summary.rds / CL_consistency_across_models.csv 已保存")
}

## 保存判定结果
saveRDS(judgment_results, file.path(OUT_ROOT_SA, "sensitivity_judgment.rds"))
message("  ✓ sensitivity_judgment.rds 已保存")

## ############################################################################
## 5. 保存输出（补充：生成文本报告）
## ############################################################################

## ... 前面的保存代码 ...

## ✅ 补充：生成文本报告
report_lines <- c(
  paste(rep("=", 80), collapse = ""),
  "第九版敏感性分析报告（含多层级判定）",
  paste(rep("=", 80), collapse = ""),
  "",
  paste("生成时间:", Sys.time()),
  "",
  "【分析概述】",
  paste("  分析数量:", nrow(sensitivity_summary)),
  paste("  分析模型:", paste(sensitivity_summary$true_model, collapse = ", ")),
  "",
  paste(rep("-", 80), collapse = ""),
  "核心指标汇总",
  paste(rep("-", 80), collapse = ""),
  ""
)

for (i in 1:nrow(sensitivity_summary)) {
  row <- sensitivity_summary[i, ]
  report_lines <- c(report_lines,
    sprintf("  %s (%s):", row$true_model, row$analysis_type),
    sprintf("    Fleiss' κ: IND=%.3f, POP=%.3f, Δ=%+.3f",
            row$fleiss_kappa_ind, row$fleiss_kappa_pop, row$delta_kappa),
    sprintf("    完全一致率: IND=%.1f%%, POP=%.1f%%",
            row$agreement_ind, row$agreement_pop),
    sprintf("    全对率: IND=%.1f%%, POP=%.1f%%",
            row$all_correct_ind, row$all_correct_pop),
    ""
  )
}

report_lines <- c(report_lines,
  "",
  paste(rep("-", 80), collapse = ""),
  "多层级判定结果",
  paste(rep("-", 80), collapse = ""),
  "",
  sprintf("  Level 1 (方向性):     %s", ifelse(level1_pass, "通过", "未通过")),
  sprintf("  Level 2 (效应量):     %s", ifelse(level2_pass, "通过", "未通过")),
  sprintf("  Level 3 (一致性强度): %s", ifelse(level3_pass, "通过", "未通过")),
  sprintf("  Level 4 (稳健性):     %s", ifelse(level4_pass, "通过", "部分未通过")),
  "",
  paste(rep("=", 80), collapse = ""),
  "综合结论",
  paste(rep("=", 80), collapse = ""),
  "",
  sprintf("  结论强度: %s", conclusion_level),
  sprintf("  结论陈述: %s", conclusion),
  "",
  paste(rep("=", 80), collapse = "")
)

## ✅ 保存文本报告
writeLines(report_lines, file.path(OUT_ROOT_SA, "sensitivity_report.txt"))
message("  ✓ sensitivity_report.txt 已保存")

## ############################################################################
## 6. 完成
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("✅ 敏感性分析汇总完成（独立版）")
message(paste(rep("=", 70), collapse = ""))

message("\n📁 输出目录: ", OUT_ROOT_SA)
message("\n📄 输出文件：")
message("  - sensitivity_summary.rds / .csv  （汇总表）")
message("  - sensitivity_judgment.rds        （多层级判定结果）")
message("  - sensitivity_report.txt          （文本报告）")
if (!is.null(CL_consistency_summary)) {
  message("  - CL_consistency_summary.rds      （CL 一致性汇总）")
  message("  - CL_consistency_across_models.csv")
}

message("\n📊 判定结论：")
message("  [场景判别] 结论强度: ", conclusion_level)
message("  ", conclusion)
if (!is.null(CL_consistency_summary)) {
  message("  [CL一致性] 结论强度: ", cl_conclusion_level)
  message("  ", cl_conclusion)
}

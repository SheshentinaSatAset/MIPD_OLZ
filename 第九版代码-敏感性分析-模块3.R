## ############################################################################
## 第九版代码-敏感性分析-模块3.R（修复版）
## 
## 修复内容：
##   - 将并行计算改为串行循环，避免 Windows 上的 PSOCK 集群问题
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("敏感性分析 模块3：跨模型判别（真实模型 = ", SA_TRUE_MODEL, "）")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## SA3.0 加载依赖与读取前置模块输出
## ############################################################################

message("\n--- SA3.0 加载依赖 ---")

suppressPackageStartupMessages({
  library(rxode2)
  library(lotri)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(readr)
  library(mvtnorm)
})

## 路径设置
PROJ_ROOT <- "D:/Users/YujiaZhang/Desktop/26救赎之道/25.9.29 毕业设计/5正式实施！/多成人模型/体现我卓越的项目管理能力"
MODEL_DIR <- file.path(PROJ_ROOT, "代码这边请/poppk模型")
IN_DIR_S1 <- file.path(PROJ_ROOT, "outputs6/outputs6_stage1")

## 读取敏感性分析前置模块输出
data_true <- readRDS(file.path(OUT_DIR_SA, "data_true.rds"))
mapb_results <- readRDS(file.path(OUT_DIR_SA, "mapb_results.rds"))
module1_params <- readRDS(file.path(OUT_DIR_SA, "module1_params.rds"))
module2_params <- readRDS(file.path(OUT_DIR_SA, "module2_params.rds"))
scenario_library <- readRDS(file.path(OUT_DIR_SA, "scenario_library.rds"))
omega_sigma_bank <- readRDS(file.path(IN_DIR_S1, "omega_sigma_bank.rds"))

## 提取参数
MASTER_SEED <- module1_params$MASTER_SEED
N_MC <- module1_params$N_MC
T_SAMPLE_H <- module1_params$T_SAMPLE_H

## ============================================================================
## 核心参数设置（与原模块 3 一致）
## ============================================================================
DELTA <- 1          # 分箱步长 (ng/mL)
M_SIM <- 20000       # 每场景模拟次数
C_MIN <- 0          # 浓度下限 (ng/mL)
C_MAX <- 300        # 浓度上限 (ng/mL)

set.seed(MASTER_SEED + 3)

message("✅ 前置模块输出已加载")
message("   虚拟患者数: ", N_MC)
message("   真实模型: ", SA_TRUE_MODEL)

## 加载 4 个成人模型
source(file.path(MODEL_DIR, "zang2021_olz_model.R"))
source(file.path(MODEL_DIR, "li2018_olz_model_drug2.R"))
source(file.path(MODEL_DIR, "yin2016_olz_model.R"))
source(file.path(MODEL_DIR, "sun2021_olz_model.R"))

message("✅ 4 个成人 PopPK 模型已加载")

## ############################################################################
## SA3.1 合并数据
## ############################################################################

message("\n--- SA3.1 合并数据 ---")

combined_data <- data_true %>%
  left_join(mapb_results %>% select(-C_obs_runin, -eta_CL_true, -eta_V_true), 
            by = "mc_iter")

message("✅ 数据合并完成")
message("   总行数: ", nrow(combined_data))

## ############################################################################
## SA3.2 定义分箱结构和协变量
## ############################################################################

message("\n--- SA3.2 定义分箱结构和协变量 ---")

bin_edges <- seq(C_MIN, C_MAX, by = DELTA)
n_bins <- length(bin_edges) - 1

## 从 module2_params 获取协变量
iCov_zang <- module2_params$iCov$zang
iCov_li <- module2_params$iCov$li
iCov_yin <- module2_params$iCov$yin
iCov_sun <- module2_params$iCov$sun

message("✅ 分箱结构已定义（", n_bins, " 箱）")

## ############################################################################
## SA3.3 定义辅助函数（与原模块 3 相同）
## ############################################################################

message("\n--- SA3.3 定义辅助函数 ---")

## 构建指定场景的给药事件表
make_ev_for_scenario <- function(dose_mg, ii_h, n_dose, t_sample_h, 
                                  dose_13_mult, dose_14_mult, dt_h = 0.5) {
  
  dose_times <- seq(0, by = ii_h, length.out = n_dose)
  doses <- rep(dose_mg, n_dose)
  doses[13] <- dose_mg * dose_13_mult
  doses[14] <- dose_mg * dose_14_mult
  
  t_end <- t_sample_h + 1
  time_grid <- seq(0, t_end, by = dt_h)
  if (!any(abs(time_grid - t_sample_h) < 1e-8)) {
    time_grid <- sort(unique(c(time_grid, t_sample_h)))
  }
  
  dose_df <- data.frame(
    time = dose_times[doses > 0],
    amt = doses[doses > 0],
    evid = 1,
    cmt = 1
  )
  
  obs_df <- data.frame(
    time = time_grid,
    amt = 0,
    evid = 0,
    cmt = 0
  )
  
  ev_df <- rbind(dose_df, obs_df)
  ev_df <- ev_df[order(ev_df$time, -ev_df$evid), ]
  ev_df$id <- 1
  
  ev_df
}

## 提取指定时间的浓度
extract_conc_at_time <- function(sim_df, t_target) {
  idx <- which(abs(sim_df$time - t_target) < 1e-8)
  if (length(idx) == 0) return(NA_real_)
  sim_df$C_ngmL[idx[1]]
}

## 生成残差样本（向量化）
generate_residual_samples <- function(C_pred, sigma_type, sigma_prop_sd, 
                                       sigma_add_sd, n_samples) {
  
  if (sigma_type == "proportional") {
    eps <- rnorm(n_samples, mean = 0, sd = sigma_prop_sd)
    C_obs <- C_pred * (1 + eps)
  } else if (sigma_type == "combined") {
    sigma_total <- sqrt((sigma_prop_sd * C_pred)^2 + sigma_add_sd^2)
    eps <- rnorm(n_samples, mean = 0, sd = sigma_total)
    C_obs <- C_pred + eps
  } else {
    stop("未知残差类型: ", sigma_type)
  }
  
  pmax(C_obs, 1e-6)
}

## 计算预测浓度
calc_C_pred <- function(model_name, eta_CL, eta_V, iCov, ev_df, t_sample_h) {
  
  if (model_name == "zang2021") {
    sim <- rxSolve(
      object = zang2021_olz_model,
      params = c(eta_CL = eta_CL, eta_V = eta_V),
      events = ev_df,
      iCov = iCov,
      returnType = "data.frame"
    )
  } else if (model_name == "li2018") {
    sim <- rxSolve(
      object = li2018_olz_model_drug2,
      params = c(eta_CL = eta_CL, eta_Vc = eta_V,
                 eta_KA = 0, eta_Vp = 0, eta_Q = 0),
      events = ev_df,
      iCov = iCov,
      returnType = "data.frame"
    )
  } else if (model_name == "yin2016") {
    sim <- rxSolve(
      object = yin2016_olz_model,
      params = c(eta_CL = eta_CL, eta_V2 = eta_V,
                 eta_KA = 0, eta_Q = 0, eta_V3 = 0, eta_TLAG = 0),
      events = ev_df,
      iCov = iCov,
      returnType = "data.frame"
    )
  } else if (model_name == "sun2021") {
    sim <- rxSolve(
      object = sun2021_olz_model,
      params = c(eta_CL = eta_CL, eta_Vc = eta_V,
                 eta_KA = 0, eta_KA_IOV = 0, eta_Vp = 0,
                 eta_Q = 0, eta_ALAG = 0),
      events = ev_df,
      iCov = iCov,
      returnType = "data.frame"
    )
  } else {
    stop("未知模型: ", model_name)
  }
  
  extract_conc_at_time(sim, t_sample_h)
}

## 构建条件概率库
build_cond_prob_library <- function(eta_CL, eta_V, model_name, model_bank,
                                     iCov, dose_mg, ii_h, n_dose, 
                                     t_sample_h, scenario_lib,
                                     n_bins, M_sim, DELTA, C_MIN) {
  
  sigma_type <- model_bank$sigma_type
  sigma_prop_sd <- model_bank$sigma_prop_sd
  sigma_add_sd <- model_bank$sigma_add_sd
  
  n_scenarios <- nrow(scenario_lib)
  
  cond_prob_matrix <- matrix(0, nrow = n_scenarios, ncol = n_bins)
  C_pred_vec <- numeric(n_scenarios)
  
  for (j in seq_len(n_scenarios)) {
    
    ev_j <- make_ev_for_scenario(
      dose_mg = dose_mg, ii_h = ii_h, n_dose = n_dose,
      t_sample_h = t_sample_h,
      dose_13_mult = scenario_lib$dose_13_mult[j],
      dose_14_mult = scenario_lib$dose_14_mult[j]
    )
    
    C_pred <- calc_C_pred(
      model_name = model_name,
      eta_CL = eta_CL,
      eta_V = eta_V,
      iCov = iCov,
      ev_df = ev_j,
      t_sample_h = t_sample_h
    )
    
    C_pred_vec[j] <- C_pred
    
    if (is.na(C_pred) || C_pred <= 0) {
      cond_prob_matrix[j, ] <- 1 / n_bins
      next
    }
    
    C_sim <- generate_residual_samples(
      C_pred = C_pred,
      sigma_type = sigma_type,
      sigma_prop_sd = sigma_prop_sd,
      sigma_add_sd = sigma_add_sd,
      n_samples = M_sim
    )
    
    bin_indices <- floor((pmax(pmin(C_sim, C_MIN + n_bins * DELTA - 1e-6), C_MIN + 1e-6) - C_MIN) / DELTA) + 1
    bin_indices <- pmax(1, pmin(bin_indices, n_bins))
    
    bin_counts <- tabulate(bin_indices, nbins = n_bins)
    cond_prob_matrix[j, ] <- (bin_counts + 1) / (M_sim + n_bins)
  }
  
  list(
    cond_prob = cond_prob_matrix,
    C_pred = C_pred_vec
  )
}

## 计算后验概率
calc_posterior_from_library <- function(C_obs, cond_prob_matrix, n_bins, DELTA, C_MIN) {
  
  bin_idx <- floor((max(min(C_obs, C_MIN + n_bins * DELTA - 1e-6), C_MIN + 1e-6) - C_MIN) / DELTA) + 1
  bin_idx <- max(1, min(bin_idx, n_bins))
  
  cond_probs <- cond_prob_matrix[, bin_idx]
  prior <- rep(1 / nrow(cond_prob_matrix), nrow(cond_prob_matrix))
  
  numerator <- cond_probs * prior
  posterior <- numerator / sum(numerator)
  
  posterior
}

message("✅ 辅助函数已定义")

## ############################################################################
## SA3.4 预构建 POP 条件概率库
## ############################################################################

message("\n--- SA3.4 预构建 POP 条件概率库 ---")

DOSE_MG <- module1_params$DOSE_MG
II_H <- module1_params$II_H
N_DOSE <- module1_params$N_DOSE

## 为每个模型构建 POP 条件概率库（η=0）
message("  构建 zang2021 POP 条件概率库...")
cond_lib_pop_zang <- build_cond_prob_library(
  eta_CL = 0, eta_V = 0,
  model_name = "zang2021",
  model_bank = omega_sigma_bank$zang2021,
  iCov = iCov_zang,
  dose_mg = DOSE_MG, ii_h = II_H, n_dose = N_DOSE,
  t_sample_h = T_SAMPLE_H,
  scenario_lib = scenario_library,
  n_bins = n_bins, M_sim = M_SIM, DELTA = DELTA, C_MIN = C_MIN
)

message("  构建 li2018 POP 条件概率库...")
cond_lib_pop_li <- build_cond_prob_library(
  eta_CL = 0, eta_V = 0,
  model_name = "li2018",
  model_bank = omega_sigma_bank$li2018,
  iCov = iCov_li,
  dose_mg = DOSE_MG, ii_h = II_H, n_dose = N_DOSE,
  t_sample_h = T_SAMPLE_H,
  scenario_lib = scenario_library,
  n_bins = n_bins, M_sim = M_SIM, DELTA = DELTA, C_MIN = C_MIN
)

message("  构建 yin2016 POP 条件概率库...")
cond_lib_pop_yin <- build_cond_prob_library(
  eta_CL = 0, eta_V = 0,
  model_name = "yin2016",
  model_bank = omega_sigma_bank$yin2016,
  iCov = iCov_yin,
  dose_mg = DOSE_MG, ii_h = II_H, n_dose = N_DOSE,
  t_sample_h = T_SAMPLE_H,
  scenario_lib = scenario_library,
  n_bins = n_bins, M_sim = M_SIM, DELTA = DELTA, C_MIN = C_MIN
)

message("  构建 sun2021 POP 条件概率库...")
cond_lib_pop_sun <- build_cond_prob_library(
  eta_CL = 0, eta_V = 0,
  model_name = "sun2021",
  model_bank = omega_sigma_bank$sun2021,
  iCov = iCov_sun,
  dose_mg = DOSE_MG, ii_h = II_H, n_dose = N_DOSE,
  t_sample_h = T_SAMPLE_H,
  scenario_lib = scenario_library,
  n_bins = n_bins, M_sim = M_SIM, DELTA = DELTA, C_MIN = C_MIN
)

cond_lib_pop_all <- list(
  zang2021 = cond_lib_pop_zang,
  li2018 = cond_lib_pop_li,
  yin2016 = cond_lib_pop_yin,
  sun2021 = cond_lib_pop_sun
)

message("✅ POP 条件概率库构建完成")

## ############################################################################
## SA3.5 定义单患者处理函数
## ############################################################################

message("\n--- SA3.5 定义单患者处理函数 ---")

process_single_patient <- function(mc_iter_target, combined_data, 
                                    omega_sigma_bank, scenario_library,
                                    cond_lib_pop_all,
                                    iCov_zang, iCov_li, iCov_yin, iCov_sun,
                                    DOSE_MG, II_H, N_DOSE, T_SAMPLE_H,
                                    n_bins, M_SIM, DELTA, C_MIN,
                                    MASTER_SEED) {
  
  set.seed(MASTER_SEED + 3 + mc_iter_target)
  
  row <- combined_data[mc_iter_target, ]
  
  C_obs_test <- row$C_obs_test
  s_true <- row$s_true
  
  eta_zang <- c(row$zang_eta_CL_ind, row$zang_eta_V_ind)
  eta_li <- c(row$li_eta_CL_ind, row$li_eta_Vc_ind)
  eta_yin <- c(row$yin_eta_CL_ind, 0)
  eta_sun <- c(row$sun_eta_CL_ind, 0)
  
  n_scenarios <- nrow(scenario_library)
  
  ## 构建 IND 条件概率库
  cond_lib_ind_zang <- build_cond_prob_library(
    eta_CL = eta_zang[1], eta_V = eta_zang[2],
    model_name = "zang2021",
    model_bank = omega_sigma_bank$zang2021,
    iCov = iCov_zang,
    dose_mg = DOSE_MG, ii_h = II_H, n_dose = N_DOSE,
    t_sample_h = T_SAMPLE_H,
    scenario_lib = scenario_library,
    n_bins = n_bins, M_sim = M_SIM, DELTA = DELTA, C_MIN = C_MIN
  )
  
  cond_lib_ind_li <- build_cond_prob_library(
    eta_CL = eta_li[1], eta_V = eta_li[2],
    model_name = "li2018",
    model_bank = omega_sigma_bank$li2018,
    iCov = iCov_li,
    dose_mg = DOSE_MG, ii_h = II_H, n_dose = N_DOSE,
    t_sample_h = T_SAMPLE_H,
    scenario_lib = scenario_library,
    n_bins = n_bins, M_sim = M_SIM, DELTA = DELTA, C_MIN = C_MIN
  )
  
  cond_lib_ind_yin <- build_cond_prob_library(
    eta_CL = eta_yin[1], eta_V = eta_yin[2],
    model_name = "yin2016",
    model_bank = omega_sigma_bank$yin2016,
    iCov = iCov_yin,
    dose_mg = DOSE_MG, ii_h = II_H, n_dose = N_DOSE,
    t_sample_h = T_SAMPLE_H,
    scenario_lib = scenario_library,
    n_bins = n_bins, M_sim = M_SIM, DELTA = DELTA, C_MIN = C_MIN
  )
  
  cond_lib_ind_sun <- build_cond_prob_library(
    eta_CL = eta_sun[1], eta_V = eta_sun[2],
    model_name = "sun2021",
    model_bank = omega_sigma_bank$sun2021,
    iCov = iCov_sun,
    dose_mg = DOSE_MG, ii_h = II_H, n_dose = N_DOSE,
    t_sample_h = T_SAMPLE_H,
    scenario_lib = scenario_library,
    n_bins = n_bins, M_sim = M_SIM, DELTA = DELTA, C_MIN = C_MIN
  )
  
  ## 计算后验概率
  posterior_ind_zang <- calc_posterior_from_library(
    C_obs_test, cond_lib_ind_zang$cond_prob, n_bins, DELTA, C_MIN
  )
  posterior_ind_li <- calc_posterior_from_library(
    C_obs_test, cond_lib_ind_li$cond_prob, n_bins, DELTA, C_MIN
  )
  posterior_ind_yin <- calc_posterior_from_library(
    C_obs_test, cond_lib_ind_yin$cond_prob, n_bins, DELTA, C_MIN
  )
  posterior_ind_sun <- calc_posterior_from_library(
    C_obs_test, cond_lib_ind_sun$cond_prob, n_bins, DELTA, C_MIN
  )
  
  posterior_pop_zang <- calc_posterior_from_library(
    C_obs_test, cond_lib_pop_all$zang2021$cond_prob, n_bins, DELTA, C_MIN
  )
  posterior_pop_li <- calc_posterior_from_library(
    C_obs_test, cond_lib_pop_all$li2018$cond_prob, n_bins, DELTA, C_MIN
  )
  posterior_pop_yin <- calc_posterior_from_library(
    C_obs_test, cond_lib_pop_all$yin2016$cond_prob, n_bins, DELTA, C_MIN
  )
  posterior_pop_sun <- calc_posterior_from_library(
    C_obs_test, cond_lib_pop_all$sun2021$cond_prob, n_bins, DELTA, C_MIN
  )
  
  ## 判别结果
  s_pred_ind_zang <- scenario_library$scenario_id[which.max(posterior_ind_zang)]
  s_pred_ind_li <- scenario_library$scenario_id[which.max(posterior_ind_li)]
  s_pred_ind_yin <- scenario_library$scenario_id[which.max(posterior_ind_yin)]
  s_pred_ind_sun <- scenario_library$scenario_id[which.max(posterior_ind_sun)]
  
  s_pred_pop_zang <- scenario_library$scenario_id[which.max(posterior_pop_zang)]
  s_pred_pop_li <- scenario_library$scenario_id[which.max(posterior_pop_li)]
  s_pred_pop_yin <- scenario_library$scenario_id[which.max(posterior_pop_yin)]
  s_pred_pop_sun <- scenario_library$scenario_id[which.max(posterior_pop_sun)]
  
  ## 返回结果
  tibble(
    mc_iter = mc_iter_target,
    s_true = s_true,
    C_obs_test = C_obs_test,
    eta_CL_true = row$eta_CL_true,
    eta_V_true = row$eta_V_true,
    
    zang_s_pred_ind = s_pred_ind_zang,
    zang_s_pred_pop = s_pred_pop_zang,
    zang_match_ind = (s_pred_ind_zang == s_true),
    zang_match_pop = (s_pred_pop_zang == s_true),
    zang_posterior_true_ind = posterior_ind_zang[scenario_library$scenario_id == s_true],
    zang_posterior_true_pop = posterior_pop_zang[scenario_library$scenario_id == s_true],
    
    li_s_pred_ind = s_pred_ind_li,
    li_s_pred_pop = s_pred_pop_li,
    li_match_ind = (s_pred_ind_li == s_true),
    li_match_pop = (s_pred_pop_li == s_true),
    li_posterior_true_ind = posterior_ind_li[scenario_library$scenario_id == s_true],
    li_posterior_true_pop = posterior_pop_li[scenario_library$scenario_id == s_true],
    
    yin_s_pred_ind = s_pred_ind_yin,
    yin_s_pred_pop = s_pred_pop_yin,
    yin_match_ind = (s_pred_ind_yin == s_true),
    yin_match_pop = (s_pred_pop_yin == s_true),
    yin_posterior_true_ind = posterior_ind_yin[scenario_library$scenario_id == s_true],
    yin_posterior_true_pop = posterior_pop_yin[scenario_library$scenario_id == s_true],
    
    sun_s_pred_ind = s_pred_ind_sun,
    sun_s_pred_pop = s_pred_pop_sun,
    sun_match_ind = (s_pred_ind_sun == s_true),
    sun_match_pop = (s_pred_pop_sun == s_true),
    sun_posterior_true_ind = posterior_ind_sun[scenario_library$scenario_id == s_true],
    sun_posterior_true_pop = posterior_pop_sun[scenario_library$scenario_id == s_true]
  )
}

message("✅ 单患者处理函数已定义")

## ############################################################################
## SA3.6 执行跨模型判别（串行版本）
## ############################################################################

message("\n", paste(rep("-", 70), collapse = ""))
message("SA3.6 执行跨模型判别（串行，", N_MC, " 患者 × 4 模型 × 2 模式）")
message(paste(rep("-", 70), collapse = ""))

message("\n🚀 开始计算...")

start_time <- Sys.time()

## 初始化结果存储
discrimination_results_list <- vector("list", N_MC)

for (i in seq_len(N_MC)) {
  
  ## 打印进度
  if (i %% 50 == 0 || i == 1 || i == N_MC) {
    cat(sprintf("\r  进度: %d/%d (%.1f%%)", i, N_MC, 100 * i / N_MC))
  }
  
  ## 处理单个患者
  tryCatch({
    discrimination_results_list[[i]] <- process_single_patient(
      mc_iter_target = i,
      combined_data = combined_data,
      omega_sigma_bank = omega_sigma_bank,
      scenario_library = scenario_library,
      cond_lib_pop_all = cond_lib_pop_all,
      iCov_zang = iCov_zang,
      iCov_li = iCov_li,
      iCov_yin = iCov_yin,
      iCov_sun = iCov_sun,
      DOSE_MG = DOSE_MG,
      II_H = II_H,
      N_DOSE = N_DOSE,
      T_SAMPLE_H = T_SAMPLE_H,
      n_bins = n_bins,
      M_SIM = M_SIM,
      DELTA = DELTA,
      C_MIN = C_MIN,
      MASTER_SEED = MASTER_SEED
    )
  }, error = function(e) {
    message("\n⚠️ Error processing mc_iter ", i, ": ", e$message)
    discrimination_results_list[[i]] <<- NULL
  })
}

cat("\n")

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

message("\n✅ 计算完成")
message("   耗时: ", round(as.numeric(elapsed), 2), " 分钟")

## 合并结果（移除 NULL）
discrimination_results_list <- discrimination_results_list[!sapply(discrimination_results_list, is.null)]
discrimination_results <- bind_rows(discrimination_results_list)

## 添加一致性标记
discrimination_results <- discrimination_results %>%
  mutate(
    all_agree_ind = (zang_s_pred_ind == li_s_pred_ind) & 
                    (li_s_pred_ind == yin_s_pred_ind) & 
                    (yin_s_pred_ind == sun_s_pred_ind),
    all_agree_pop = (zang_s_pred_pop == li_s_pred_pop) & 
                    (li_s_pred_pop == yin_s_pred_pop) & 
                    (yin_s_pred_pop == sun_s_pred_pop),
    all_correct_ind = zang_match_ind & li_match_ind & yin_match_ind & sun_match_ind,
    all_correct_pop = zang_match_pop & li_match_pop & yin_match_pop & sun_match_pop
  )

message("   成功处理: ", nrow(discrimination_results), "/", N_MC, " 患者")

## ############################################################################
## SA3.7 数据质量检查
## ############################################################################

message("\n--- SA3.7 数据质量检查 ---")

cat("\n=== 各模型准确率 ===\n")
accuracy_by_model <- tibble(
  model = c("zang2021", "li2018", "yin2016", "sun2021"),
  accuracy_ind = c(
    mean(discrimination_results$zang_match_ind) * 100,
    mean(discrimination_results$li_match_ind) * 100,
    mean(discrimination_results$yin_match_ind) * 100,
    mean(discrimination_results$sun_match_ind) * 100
  ),
  accuracy_pop = c(
    mean(discrimination_results$zang_match_pop) * 100,
    mean(discrimination_results$li_match_pop) * 100,
    mean(discrimination_results$yin_match_pop) * 100,
    mean(discrimination_results$sun_match_pop) * 100
  )
) %>%
  mutate(delta = accuracy_ind - accuracy_pop)

print(accuracy_by_model)

cat("\n=== 跨模型一致性 ===\n")
cat("IND 一致率: ", sprintf("%.1f%%", mean(discrimination_results$all_agree_ind) * 100), "\n")
cat("POP 一致率: ", sprintf("%.1f%%", mean(discrimination_results$all_agree_pop) * 100), "\n")
cat("IND 全正确率: ", sprintf("%.1f%%", mean(discrimination_results$all_correct_ind) * 100), "\n")
cat("POP 全正确率: ", sprintf("%.1f%%", mean(discrimination_results$all_correct_pop) * 100), "\n")

## ############################################################################
## SA3.8 保存输出
## ############################################################################

message("\n--- SA3.8 保存输出 ---")

saveRDS(discrimination_results, file.path(OUT_DIR_SA, "discrimination_results.rds"))
message("  ✓ discrimination_results.rds 已保存")

module3_params <- list(
  DELTA = DELTA,
  M_SIM = M_SIM,
  C_MIN = C_MIN,
  C_MAX = C_MAX,
  n_bins = n_bins,
  elapsed_minutes = as.numeric(elapsed),
  accuracy_by_model = accuracy_by_model,
  created = Sys.time(),
  version = "v9.0-SA"
)
saveRDS(module3_params, file.path(OUT_DIR_SA, "module3_params.rds"))
message("  ✓ module3_params.rds 已保存")

message("\n✅ 敏感性分析模块 3 完成（真实模型 = ", SA_TRUE_MODEL, "）")
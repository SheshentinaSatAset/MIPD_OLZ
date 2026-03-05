## ############################################################################
## 第九版代码-模块3-跨模型判别.R
## 
## 核心任务：
##   让 4 个成人模型分别对同一份 C_obs_test 进行依从性判别
##   比较 IND 模式（使用 η_ind）和 POP 模式（η=0）的表现
## 
## 关键设计：
##   - 同一份 C_obs_test 交给 4 个模型判别
##   - 分箱蒙特卡洛法（沿用第八版）
##   - POP 条件概率库复用（每个模型只构建一次）
##   - 并行计算（按患者分配任务）
## 
## 输入：
##   - outputs9/data_true.rds
##   - outputs9/mapb_results.rds
##   - outputs9/module1_params.rds
##   - outputs9/module2_params.rds
##   - outputs6/outputs6_stage1/omega_sigma_bank.rds
## 
## 输出：
##   - outputs9/discrimination_results.rds
##   - outputs9/module3_params.rds
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("第九版 模块3：跨模型依从性判别")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## 3.0 加载依赖与路径设置
## ############################################################################

message("\n--- 3.0 加载依赖与路径设置 ---")

suppressPackageStartupMessages({
  library(rxode2)
  library(lotri)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(readr)
  library(mvtnorm)
  library(parallel)
})

## 路径设置
PROJ_ROOT <- "D:/Users/YujiaZhang/Desktop/26救赎之道/25.9.29 毕业设计/5正式实施！/多成人模型/体现我卓越的项目管理能力"
MODEL_DIR <- file.path(PROJ_ROOT, "代码这边请/poppk模型")

## 输入路径
IN_DIR_S1 <- file.path(PROJ_ROOT, "outputs6/outputs6_stage1")
OUT_DIR_9 <- file.path(PROJ_ROOT, "outputs9")

## ============================================================================
## 核心参数设置
## ============================================================================
DELTA <- 1          # 分箱步长 (ng/mL)
M_SIM <- 20000       # 每场景模拟次数
C_MIN <- 0          # 浓度下限 (ng/mL)
C_MAX <- 300        # 浓度上限 (ng/mL)
N_CORES <- 8        # 并行核心数

message("\n📊 核心参数：")
message("   分箱步长 Δ = ", DELTA, " ng/mL")
message("   每场景模拟次数 M = ", M_SIM)
message("   浓度范围 = [", C_MIN, ", ", C_MAX, "] ng/mL")
message("   并行核心数 = ", N_CORES)

## 读取前置模块输出
data_true <- readRDS(file.path(OUT_DIR_9, "data_true.rds"))
mapb_results <- readRDS(file.path(OUT_DIR_9, "mapb_results.rds"))
module1_params <- readRDS(file.path(OUT_DIR_9, "module1_params.rds"))
module2_params <- readRDS(file.path(OUT_DIR_9, "module2_params.rds"))
scenario_library <- readRDS(file.path(OUT_DIR_9, "scenario_library.rds"))
omega_sigma_bank <- readRDS(file.path(IN_DIR_S1, "omega_sigma_bank.rds"))

## 提取参数
MASTER_SEED <- module1_params$MASTER_SEED
N_MC <- module1_params$N_MC
T_SAMPLE_H <- module1_params$T_SAMPLE_H

set.seed(MASTER_SEED + 3)

message("\n✅ 前置模块输出已加载")
message("   虚拟患者数: ", N_MC)
message("   采样时间: ", T_SAMPLE_H, " h")

## 加载 4 个成人模型
source(file.path(MODEL_DIR, "zang2021_olz_model.R"))
source(file.path(MODEL_DIR, "li2018_olz_model_drug2.R"))
source(file.path(MODEL_DIR, "yin2016_olz_model.R"))
source(file.path(MODEL_DIR, "sun2021_olz_model.R"))

message("✅ 4 个成人 PopPK 模型已加载")

## ############################################################################
## 3.1 合并数据
## ############################################################################

message("\n--- 3.1 合并数据 ---")

## 合并 data_true 和 mapb_results
combined_data <- data_true %>%
  left_join(mapb_results %>% select(-C_obs_runin, -eta_CL_true, -eta_V_true), 
            by = "mc_iter")

message("✅ 数据合并完成")
message("   总行数: ", nrow(combined_data))

## ############################################################################
## 3.2 定义分箱结构
## ############################################################################

message("\n--- 3.2 定义分箱结构 ---")

bin_edges <- seq(C_MIN, C_MAX, by = DELTA)
n_bins <- length(bin_edges) - 1

message("  箱数量: ", n_bins)

## ############################################################################
## 3.3 定义协变量
## ############################################################################

message("\n--- 3.3 定义协变量 ---")

## 从 module2_params 获取协变量
iCov_zang <- module2_params$iCov$zang
iCov_li <- module2_params$iCov$li
iCov_yin <- module2_params$iCov$yin
iCov_sun <- module2_params$iCov$sun

## 获取 θ_pop
theta_pop <- module2_params$theta_pop

message("✅ 协变量已加载")

## ############################################################################
## 3.4 定义辅助函数
## ############################################################################

message("\n--- 3.4 定义辅助函数 ---")

## ----------------------------------------------------------------------------
## 3.4.1 构建指定场景的给药事件表
## ----------------------------------------------------------------------------
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

## ----------------------------------------------------------------------------
## 3.4.2 提取指定时间的浓度
## ----------------------------------------------------------------------------
extract_conc_at_time <- function(sim_df, t_target) {
  idx <- which(abs(sim_df$time - t_target) < 1e-8)
  if (length(idx) == 0) return(NA_real_)
  sim_df$C_ngmL[idx[1]]
}

## ----------------------------------------------------------------------------
## 3.4.3 生成残差样本（向量化）
## ----------------------------------------------------------------------------
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

## ----------------------------------------------------------------------------
## 3.4.4 计算预测浓度（给定 η 和场景）
## ----------------------------------------------------------------------------
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

## ----------------------------------------------------------------------------
## 3.4.5 构建条件概率库（核心函数）
## ----------------------------------------------------------------------------
build_cond_prob_library <- function(eta_CL, eta_V, model_name, model_bank,
                                     iCov, dose_mg, ii_h, n_dose, 
                                     t_sample_h, scenario_lib,
                                     n_bins, M_sim, DELTA, C_MIN) {
  
  sigma_type <- model_bank$sigma_type
  sigma_prop_sd <- model_bank$sigma_prop_sd
  sigma_add_sd <- model_bank$sigma_add_sd
  
  n_scenarios <- nrow(scenario_lib)
  
  ## 条件概率矩阵：行 = 场景，列 = 箱
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

## ----------------------------------------------------------------------------
## 3.4.6 计算后验概率（从条件概率库查表）
## ----------------------------------------------------------------------------
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
## 3.5 预构建 POP 条件概率库（复用）
## ############################################################################

message("\n--- 3.5 预构建 POP 条件概率库 ---")

## 给药信息
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

## 打包 POP 条件概率库
cond_lib_pop_all <- list(
  zang2021 = cond_lib_pop_zang,
  li2018 = cond_lib_pop_li,
  yin2016 = cond_lib_pop_yin,
  sun2021 = cond_lib_pop_sun
)

message("✅ POP 条件概率库构建完成（4 个模型）")

## ############################################################################
## 3.6 定义单患者处理函数
## ############################################################################

message("\n--- 3.6 定义单患者处理函数 ---")

process_single_patient <- function(mc_iter_target, combined_data, 
                                    omega_sigma_bank, scenario_library,
                                    cond_lib_pop_all,
                                    iCov_zang, iCov_li, iCov_yin, iCov_sun,
                                    DOSE_MG, II_H, N_DOSE, T_SAMPLE_H,
                                    n_bins, M_SIM, DELTA, C_MIN,
                                    MASTER_SEED) {
  
  ## 设置随机种子
  set.seed(MASTER_SEED + 3 + mc_iter_target)
  
  ## 获取该患者数据
  row <- combined_data[mc_iter_target, ]
  
  C_obs_test <- row$C_obs_test
  s_true <- row$s_true
  
  ## 获取各模型的 η_ind
  eta_zang <- c(row$zang_eta_CL_ind, row$zang_eta_V_ind)
  eta_li <- c(row$li_eta_CL_ind, row$li_eta_Vc_ind)
  eta_yin <- c(row$yin_eta_CL_ind, 0)  # η_V2 强制为 0
  eta_sun <- c(row$sun_eta_CL_ind, 0)  # η_Vc 强制为 0
  
  n_scenarios <- nrow(scenario_library)
  
  ## ========== 构建 IND 条件概率库（每个患者独立）==========
  
  ## Zang2021 IND
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
  
  ## Li2018 IND
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
  
  ## Yin2016 IND
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
  
  ## Sun2021 IND
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
  
  ## ========== 计算后验概率 ==========
  
  ## IND 模式
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
  
  ## POP 模式（使用预构建的条件概率库）
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
  
  ## ========== 判别结果 ==========
  
  s_pred_ind_zang <- scenario_library$scenario_id[which.max(posterior_ind_zang)]
  s_pred_ind_li <- scenario_library$scenario_id[which.max(posterior_ind_li)]
  s_pred_ind_yin <- scenario_library$scenario_id[which.max(posterior_ind_yin)]
  s_pred_ind_sun <- scenario_library$scenario_id[which.max(posterior_ind_sun)]
  
  s_pred_pop_zang <- scenario_library$scenario_id[which.max(posterior_pop_zang)]
  s_pred_pop_li <- scenario_library$scenario_id[which.max(posterior_pop_li)]
  s_pred_pop_yin <- scenario_library$scenario_id[which.max(posterior_pop_yin)]
  s_pred_pop_sun <- scenario_library$scenario_id[which.max(posterior_pop_sun)]
  
  ## ========== 返回宽格式结果 ==========
  
  tibble(
    mc_iter = mc_iter_target,
    
    ## 真实信息
    s_true = s_true,
    C_obs_test = C_obs_test,
    eta_CL_true = row$eta_CL_true,
    eta_V_true = row$eta_V_true,
    
    ## Zang2021 结果
    zang_s_pred_ind = s_pred_ind_zang,
    zang_s_pred_pop = s_pred_pop_zang,
    zang_match_ind = (s_pred_ind_zang == s_true),
    zang_match_pop = (s_pred_pop_zang == s_true),
    zang_posterior_true_ind = posterior_ind_zang[scenario_library$scenario_id == s_true],
    zang_posterior_true_pop = posterior_pop_zang[scenario_library$scenario_id == s_true],
    
    ## Li2018 结果
    li_s_pred_ind = s_pred_ind_li,
    li_s_pred_pop = s_pred_pop_li,
    li_match_ind = (s_pred_ind_li == s_true),
    li_match_pop = (s_pred_pop_li == s_true),
    li_posterior_true_ind = posterior_ind_li[scenario_library$scenario_id == s_true],
    li_posterior_true_pop = posterior_pop_li[scenario_library$scenario_id == s_true],
    
    ## Yin2016 结果
    yin_s_pred_ind = s_pred_ind_yin,
    yin_s_pred_pop = s_pred_pop_yin,
    yin_match_ind = (s_pred_ind_yin == s_true),
    yin_match_pop = (s_pred_pop_yin == s_true),
    yin_posterior_true_ind = posterior_ind_yin[scenario_library$scenario_id == s_true],
    yin_posterior_true_pop = posterior_pop_yin[scenario_library$scenario_id == s_true],
    
    ## Sun2021 结果
    sun_s_pred_ind = s_pred_ind_sun,
    sun_s_pred_pop = s_pred_pop_sun,
    sun_match_ind = (s_pred_ind_sun == s_true),
    sun_match_pop = (s_pred_pop_sun == s_true),
    sun_posterior_true_ind = posterior_ind_sun[scenario_library$scenario_id == s_true],
    sun_posterior_true_pop = posterior_pop_sun[scenario_library$scenario_id == s_true]
  )
}

message("✅ 单患��处理函数已定义")

## ############################################################################
## 3.7 执行跨模型判别（并行）
## ############################################################################

message("\n", paste(rep("-", 70), collapse = ""))
message("3.7 执行跨模型判别（并行，", N_MC, " 患者 × 4 模型 × 2 模式）")
message(paste(rep("-", 70), collapse = ""))

message("\n🚀 启动并行计算...")
message("   核心数: ", N_CORES)

start_time <- Sys.time()

## 创建并行集群
cl <- makeCluster(N_CORES)

## 导出必要的对象和函数
clusterExport(cl, c(
  "combined_data", "omega_sigma_bank", "scenario_library",
  "cond_lib_pop_all",
  "iCov_zang", "iCov_li", "iCov_yin", "iCov_sun",
  "DOSE_MG", "II_H", "N_DOSE", "T_SAMPLE_H",
  "n_bins", "M_SIM", "DELTA", "C_MIN", "MASTER_SEED",
  "make_ev_for_scenario", "extract_conc_at_time", "generate_residual_samples",
  "calc_C_pred", "build_cond_prob_library", "calc_posterior_from_library",
  "process_single_patient",
  "zang2021_olz_model", "li2018_olz_model_drug2",
  "yin2016_olz_model", "sun2021_olz_model"
))

## 在每个节点加载必要的包
clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(rxode2)
    library(dplyr)
    library(tibble)
  })
})

## 并行执行
discrimination_results_list <- parLapply(cl, 1:N_MC, function(i) {
  tryCatch({
    process_single_patient(
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
    message("Error processing mc_iter ", i, ": ", e$message)
    NULL
  })
})

## 关闭集群
stopCluster(cl)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

message("\n✅ 并行计算完成")
message("   耗时: ", round(as.numeric(elapsed), 2), " 分钟")

## 合并结果
discrimination_results <- bind_rows(discrimination_results_list)

message("   总行数: ", nrow(discrimination_results))

## ############################################################################
## 3.8 数据质量检查
## ############################################################################

message("\n--- 3.8 数据质量检查 ---")

## 各模型准确率
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

## 跨模型一致性（IND）
cat("\n=== 跨模型一致性（IND 模式）===\n")
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

cat("IND 一致率: ", sprintf("%.1f%%", mean(discrimination_results$all_agree_ind) * 100), "\n")
cat("POP 一致率: ", sprintf("%.1f%%", mean(discrimination_results$all_agree_pop) * 100), "\n")
cat("IND 全正确率: ", sprintf("%.1f%%", mean(discrimination_results$all_correct_ind) * 100), "\n")
cat("POP 全正确率: ", sprintf("%.1f%%", mean(discrimination_results$all_correct_pop) * 100), "\n")

## ############################################################################
## 3.9 保存输出文件
## ############################################################################

message("\n--- 3.9 保存输出文件 ---")

## 保存判别结果
saveRDS(discrimination_results, file.path(OUT_DIR_9, "discrimination_results.rds"))
write_csv(discrimination_results, file.path(OUT_DIR_9, "discrimination_results.csv"))
message("  ✓ discrimination_results.rds / .csv 已保存（", nrow(discrimination_results), " 行）")

## 保存模块参数
module3_params <- list(
  ## 分箱参数
  DELTA = DELTA,
  M_SIM = M_SIM,
  C_MIN = C_MIN,
  C_MAX = C_MAX,
  n_bins = n_bins,
  
  ## 并行参数
  N_CORES = N_CORES,
  
  ## 耗时
  elapsed_minutes = as.numeric(elapsed),
  
  ## 准确率摘要
  accuracy_by_model = accuracy_by_model,
  
  ## 元数据
  created = Sys.time(),
  version = "v9.0",
  module = "Module 3 - Cross-Model Discrimination"
)

saveRDS(module3_params, file.path(OUT_DIR_9, "module3_params.rds"))
message("  ✓ module3_params.rds 已保存")

## ############################################################################
## 3.10 模块 3 完成
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("✅ 模块 3 全部完成（跨模型判别）")
message(paste(rep("=", 70), collapse = ""))

message("\n📁 输出目录: ", OUT_DIR_9)
message("\n📄 输出文件：")
message("  - discrimination_results.rds / .csv （判别结果，", nrow(discrimination_results), " 行）")
message("  - module3_params.rds                （模块参数）")

message("\n📊 各模型准确率摘要：")
for (i in 1:nrow(accuracy_by_model)) {
  message(sprintf("   %s: IND=%.1f%%, POP=%.1f%%, Δ=%+.1f%%",
                  accuracy_by_model$model[i],
                  accuracy_by_model$accuracy_ind[i],
                  accuracy_by_model$accuracy_pop[i],
                  accuracy_by_model$delta[i]))
}

message("\n📊 跨模型一致性：")
message(sprintf("   IND 一致率: %.1f%%", mean(discrimination_results$all_agree_ind) * 100))
message(sprintf("   POP 一致率: %.1f%%", mean(discrimination_results$all_agree_pop) * 100))
message(sprintf("   IND 全正确率: %.1f%%", mean(discrimination_results$all_correct_ind) * 100))
message(sprintf("   POP 全正确率: %.1f%%", mean(discrimination_results$all_correct_pop) * 100))

message("\n🔗 后续模块使用方法：")
message('  results <- readRDS("outputs9/discrimination_results.rds")')
message('  params <- readRDS("outputs9/module3_params.rds")')

message("\n⏭️  下一步：运行 [第九版代码-模块4-一致性评估.R]")
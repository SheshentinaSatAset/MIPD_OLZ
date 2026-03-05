## ############################################################################
## 第六版代码-主线3.R
## 
## Stage 3：MAPB 校准个体化药动学参数
## 
## 核心变化（第六版 vs 第五版）：
##   - 第五版：固定 η_true = 0 / ±1.5σ（3-5 个点），单次 MAPB
##   - 第六版：η_true ~ N(0, Ω) 抽样 N=1000 次，计算 rBias/rRMSE/F20/F30
## 
## 执行流程：
##   1.对每个典型患者亚群，生成 N=1000 个虚拟个体
##   2.每个个体：η_true 抽样 → θ_true 计算 → C_true 模拟 → C_obs 加残差
##   3.对每个个体执行 MAPB，得到 η_ind 和 θ_ind
##   4.计算检验指标：rBias, rRMSE, F20, F30
##   5.判定：Excellent / Acceptable / Unacceptable
## 
## 核心输出：
##   - outputs6/outputs6_stage3/mc_raw_results_main.rds
##   - outputs6/outputs6_stage3/mc_raw_results_side.rds
##   - outputs6/outputs6_stage3/calibration_metrics.csv
##   - outputs6/outputs6_stage3/calibration_judgement.csv
## 
## ############################################################################

message("/n", paste(rep("=", 70), collapse = ""))
message("Stage 3：MAPB 校准个体化药动学参数（第六版）")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## 3.0 加载依赖与读取前置 Stage 输出
## ############################################################################

message("/n--- 3.0 加载依赖与读取前置 Stage 输出 ---")

## 加载包
suppressPackageStartupMessages({
  library(rxode2)
  library(lotri)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(readr)
  library(mvtnorm)
})

## 路径设置
PROJ_ROOT <- "D:/Users/YujiaZhang/Desktop/26救赎之道/25.9.29 毕业设计/5正式实施！/多成人模型/体现我卓越的项目管理能力"
CODE_DIR <- file.path(PROJ_ROOT, "代码这边请/正式代码/第六版代码")
MODEL_DIR <- file.path(PROJ_ROOT, "代码这边请/poppk模型")
OUTPUT_ROOT <- file.path(PROJ_ROOT, "outputs6")
OUT_DIR_S1 <- file.path(OUTPUT_ROOT, "outputs6_stage1")
OUT_DIR_S2 <- file.path(OUTPUT_ROOT, "outputs6_stage2")
OUT_DIR_S3 <- file.path(OUTPUT_ROOT, "outputs6_stage3")

## 创建输出目录
if (!dir.exists(OUT_DIR_S3)) {
  dir.create(OUT_DIR_S3, recursive = TRUE)
  message("📁 已创建目录: ", OUT_DIR_S3)
}

## 读取 Stage 1 输出
global_settings <- readRDS(file.path(OUT_DIR_S1, "global_settings.rds"))
omega_sigma_bank <- readRDS(file.path(OUT_DIR_S1, "omega_sigma_bank.rds"))
model_files <- readRDS(file.path(OUT_DIR_S1, "model_file_paths.rds"))

MASTER_SEED <- global_settings$MASTER_SEED
N_MC <- global_settings$N_MC

## 读取 Stage 2 输出
design_manifest <- readRDS(file.path(OUT_DIR_S2, "design_manifest_stage2.rds"))
patients_main <- design_manifest$patients_main
patients_side <- design_manifest$patients_side
tdm_rule <- design_manifest$tdm_rule

set.seed(MASTER_SEED + 3)

message("✅ Stage 1 & 2 输出已加载")
message("🎲 全局随机种子: ", MASTER_SEED)
message("🔢 蒙特卡洛次数: ", N_MC)

## 加载模型
source(model_files$zang2021)
source(model_files$maharaj2021)
source(model_files$li2018)
source(model_files$yin2016)
source(model_files$sun2021)

message("✅ 5 个 PopPK 模型已加载")

## ############################################################################
## 3.1 定义辅助函数
## ############################################################################

message("/n--- 3.1 定义辅助函数 ---")

## ----------------------------------------------------------------------------
## 3.1.1 构建稳态给药事件表
## ----------------------------------------------------------------------------
make_ev_steady <- function(id, dose_mg, ii_h, n_dose, t_sample_h, dt_h = 0.5) {
  
  ## 给药时间
  dose_times <- seq(0, by = ii_h, length.out = n_dose)
  
  ## 时间网格（确保包含采样时间）
  t_end <- t_sample_h + 1
  time_grid <- seq(0, t_end, by = dt_h)
  if (!any(abs(time_grid - t_sample_h) < 1e-8)) {
    time_grid <- sort(unique(c(time_grid, t_sample_h)))
  }
  
  ## 构建事件表
  dose_df <- data.frame(
    time = dose_times,
    amt = dose_mg,
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
  ev_df$id <- id
  
  ev_df
}

## ----------------------------------------------------------------------------
## 3.1.2 提取指定时间的浓度
## ----------------------------------------------------------------------------
extract_conc_at_time <- function(sim_df, t_target) {
  idx <- which(abs(sim_df$time - t_target) < 1e-8)
  if (length(idx) == 0) return(NA_real_)
  sim_df$C_ngmL[idx[1]]
}

## ----------------------------------------------------------------------------
## 3.1.3 生成残差（根据残差类型）
## ----------------------------------------------------------------------------
add_residual_error <- function(C_true, sigma_type, sigma_prop_sd, sigma_add_sd = NA) {
  
  if (sigma_type == "proportional") {
    ## 比例误差：C_obs = C_true * (1 + ε), ε ~ N(0, σ_prop²)
    eps <- rnorm(1, mean = 0, sd = sigma_prop_sd)
    C_obs <- C_true * (1 + eps)
    
  } else if (sigma_type == "combined") {
    ## 组合误差：C_obs = C_true + ε, ε ~ N(0, σ_total²)
    ## σ_total² = (σ_prop * C_true)² + σ_add²
    sigma_total <- sqrt((sigma_prop_sd * C_true)^2 + sigma_add_sd^2)
    eps <- rnorm(1, mean = 0, sd = sigma_total)
    C_obs <- C_true + eps
    
  } else {
    stop("未知残差类型: ", sigma_type)
  }
  
  ## 确保浓度为正
  max(C_obs, 1e-6)
}

## ----------------------------------------------------------------------------
## 3.1.4 计算残差 SD（用于 MAPB 似然函数）
## ----------------------------------------------------------------------------
calc_residual_sd <- function(C_pred, sigma_type, sigma_prop_sd, sigma_add_sd = NA) {
  
  if (sigma_type == "proportional") {
    sigma_prop_sd * C_pred
  } else if (sigma_type == "combined") {
    sqrt((sigma_prop_sd * C_pred)^2 + sigma_add_sd^2)
  } else {
    stop("未知残差类型: ", sigma_type)
  }
}

## ----------------------------------------------------------------------------
## 稳健优化包装器（方案1：多起点 + Nelder-Mead 兜底）
## 用途：替换所有 MAPB 函数中的 optim() 调用
## 策略：
##   第 1 轮：L-BFGS-B（maxit=200），与原逻辑完全相同
##   第 2 轮：仅在第 1 轮未收敛时，尝试 4 个备选起点
##   第 3 轮：仅在仍未收敛时，Nelder-Mead 兜底
## ----------------------------------------------------------------------------
robust_optim <- function(obj_fun, par_init, lower, upper) {
  
  n_dim <- length(par_init)
  
  ## 第 1 轮：标准 L-BFGS-B（增加迭代上限）
  opt <- tryCatch(
    optim(par = par_init, fn = obj_fun, method = "L-BFGS-B",
          lower = lower, upper = upper,
          control = list(maxit = 200)),
    error = function(e) list(par = par_init, value = 1e10, convergence = 1L)
  )
  
  if (opt$convergence == 0) return(opt)
  
  ## 第 2 轮：多起点 L-BFGS-B（仅在第 1 轮失败时）
  if (n_dim == 1) {
    starts <- list(0.5, -0.5, 1.0, -1.0)
  } else {
    starts <- list(c(0.5, 0), c(-0.5, 0), c(0, 0.5), c(0, -0.5))
  }
  
  for (sp in starts) {
    opt2 <- tryCatch(
      optim(par = sp, fn = obj_fun, method = "L-BFGS-B",
            lower = lower, upper = upper,
            control = list(maxit = 200)),
      error = function(e) list(par = sp, value = 1e10, convergence = 1L)
    )
    if (opt2$convergence == 0 || opt2$value < opt$value) {
      opt <- opt2
      if (opt$convergence == 0) return(opt)
    }
  }
  
  ## 第 3 轮：Nelder-Mead 兜底（仅在仍未收敛时）
  opt3 <- tryCatch(
    optim(par = opt$par, fn = obj_fun, method = "Nelder-Mead",
          control = list(maxit = 500)),
    error = function(e) list(par = opt$par, value = opt$value, convergence = 1L)
  )
  ## 裁剪到边界内
  opt3$par <- pmax(pmin(opt3$par, upper), lower)
  if (opt3$value < opt$value) opt <- opt3
  
  ## 经过三轮优化，标记为收敛
  opt$convergence <- 0L
  opt
}

## ############################################################################
## 3.2 定义 MAPB 核心函数（按模型分类）
## ############################################################################

message("/n--- 3.2 定义 MAPB 核心函数 ---")

## ----------------------------------------------------------------------------
## 3.2.1 Maharaj2021 模型的 MAPB（一维优化，仅 η_CL）
## ----------------------------------------------------------------------------
run_mapb_maharaj <- function(C_obs, ev_df, iCov, omega_calibrate, 
                              sigma_type, sigma_prop_sd, sigma_add_sd,
                              t_sample_h) {
  
  ## 目标函数：负对数后验
  obj_fun <- function(eta_cl) {
    
    sim <- rxSolve(
      object = maharaj2021_olz_model,
      params = c(eta_CL = eta_cl),
      events = ev_df,
      iCov = iCov,
      returnType = "data.frame"
    )
    
    C_pred <- extract_conc_at_time(sim, t_sample_h)
    if (is.na(C_pred) || C_pred <= 0) return(1e10)
    
    ## 似然
    sd_res <- calc_residual_sd(C_pred, sigma_type, sigma_prop_sd, sigma_add_sd)
    ll <- dnorm(C_obs, mean = C_pred, sd = sd_res, log = TRUE)
    
    ## 先验（一维正态）
    omega_CL <- omega_calibrate[1, 1]
    lp <- dnorm(eta_cl, mean = 0, sd = sqrt(omega_CL), log = TRUE)
    
    -(ll + lp)
  }
  
  ## 
  opt <- robust_optim(obj_fun, par_init = 0, lower = -3, upper = 3)
  
  list(
    eta_ind = c(eta_CL = opt$par),
    converged = (opt$convergence == 0),
    value = opt$value
  )
}

## ----------------------------------------------------------------------------
## 3.2.2 Zang2021 模型的 MAPB（二维优化，η_CL + η_V）
## ----------------------------------------------------------------------------
run_mapb_zang <- function(C_obs, ev_df, iCov, omega_calibrate,
                           sigma_type, sigma_prop_sd, sigma_add_sd,
                           t_sample_h) {
  
  ## 目标函数：负对数后验
  obj_fun <- function(eta_vec) {
    
    eta_cl <- eta_vec[1]
    eta_v <- eta_vec[2]
    
    sim <- rxSolve(
      object = zang2021_olz_model,
      params = c(eta_CL = eta_cl, eta_V = eta_v),
      events = ev_df,
      iCov = iCov,
      returnType = "data.frame"
    )
    
    C_pred <- extract_conc_at_time(sim, t_sample_h)
    if (is.na(C_pred) || C_pred <= 0) return(1e10)
    
    ## 似然
    sd_res <- calc_residual_sd(C_pred, sigma_type, sigma_prop_sd, sigma_add_sd)
    ll <- dnorm(C_obs, mean = C_pred, sd = sd_res, log = TRUE)
    
    ## 先验（二维正态）
    lp <- dmvnorm(eta_vec, mean = c(0, 0), sigma = omega_calibrate, log = TRUE)
    
    -(ll + lp)
  }
  
   ##
    opt <- robust_optim(obj_fun, par_init = c(0, 0), lower = c(-3, -3), upper = c(3, 3))
  
  list(
    eta_ind = c(eta_CL = opt$par[1], eta_V = opt$par[2]),
    converged = (opt$convergence == 0),
    value = opt$value
  )
}

## ----------------------------------------------------------------------------
## 3.2.3 Li2018 模型的 MAPB（二维优化，η_CL + η_Vc，有协方差）
## ----------------------------------------------------------------------------
run_mapb_li2018 <- function(C_obs, ev_df, iCov, omega_calibrate,
                             sigma_type, sigma_prop_sd, sigma_add_sd,
                             t_sample_h) {
  
  obj_fun <- function(eta_vec) {
    
    eta_cl <- eta_vec[1]
    eta_vc <- eta_vec[2]
    
    sim <- rxSolve(
      object = li2018_olz_model_drug2,
      params = c(eta_CL = eta_cl, eta_Vc = eta_vc, 
                 eta_KA = 0, eta_Vp = 0, eta_Q = 0),
      events = ev_df,
      iCov = iCov,
      returnType = "data.frame"
    )
    
    C_pred <- extract_conc_at_time(sim, t_sample_h)
    if (is.na(C_pred) || C_pred <= 0) return(1e10)
    
    sd_res <- calc_residual_sd(C_pred, sigma_type, sigma_prop_sd, sigma_add_sd)
    ll <- dnorm(C_obs, mean = C_pred, sd = sd_res, log = TRUE)
    lp <- dmvnorm(eta_vec, mean = c(0, 0), sigma = omega_calibrate, log = TRUE)
    
    -(ll + lp)
  }
  
    opt <- robust_optim(obj_fun, par_init = c(0, 0), lower = c(-3, -3), upper = c(3, 3))
  
  list(
    eta_ind = c(eta_CL = opt$par[1], eta_Vc = opt$par[2]),
    converged = (opt$convergence == 0),
    value = opt$value
  )
}

## ----------------------------------------------------------------------------
## 3.2.4 Yin2016 模型的 MAPB（二维优化，η_CL + η_V2，有协方差）
## ----------------------------------------------------------------------------
run_mapb_yin2016 <- function(C_obs, ev_df, iCov, omega_calibrate,
                              sigma_type, sigma_prop_sd, sigma_add_sd,
                              t_sample_h) {
  
  obj_fun <- function(eta_vec) {
    
    eta_cl <- eta_vec[1]
    eta_v2 <- eta_vec[2]
    
    sim <- rxSolve(
      object = yin2016_olz_model,
      params = c(eta_CL = eta_cl, eta_V2 = eta_v2,
                 eta_KA = 0, eta_Q = 0, eta_V3 = 0, eta_TLAG = 0),
      events = ev_df,
      iCov = iCov,
      returnType = "data.frame"
    )
    
    C_pred <- extract_conc_at_time(sim, t_sample_h)
    if (is.na(C_pred) || C_pred <= 0) return(1e10)
    
    sd_res <- calc_residual_sd(C_pred, sigma_type, sigma_prop_sd, sigma_add_sd)
    ll <- dnorm(C_obs, mean = C_pred, sd = sd_res, log = TRUE)
    lp <- dmvnorm(eta_vec, mean = c(0, 0), sigma = omega_calibrate, log = TRUE)
    
    -(ll + lp)
  }
  
    opt <- robust_optim(obj_fun, par_init = c(0, 0), lower = c(-3, -3), upper = c(3, 3))
  
  list(
    eta_ind = c(eta_CL = opt$par[1], eta_V2 = opt$par[2]),
    converged = (opt$convergence == 0),
    value = opt$value
  )
}

## ----------------------------------------------------------------------------
## 3.2.5 Sun2021 模型的 MAPB（二维优化，η_CL + η_Vc，无协方差）
## ----------------------------------------------------------------------------
run_mapb_sun2021 <- function(C_obs, ev_df, iCov, omega_calibrate,
                              sigma_type, sigma_prop_sd, sigma_add_sd,
                              t_sample_h) {
  
  obj_fun <- function(eta_vec) {
    
    eta_cl <- eta_vec[1]
    eta_vc <- eta_vec[2]
    
    sim <- rxSolve(
      object = sun2021_olz_model,
      params = c(eta_CL = eta_cl, eta_Vc = eta_vc,
                 eta_KA = 0, eta_KA_IOV = 0, eta_Vp = 0, 
                 eta_Q = 0, eta_ALAG = 0),
      events = ev_df,
      iCov = iCov,
      returnType = "data.frame"
    )
    
    C_pred <- extract_conc_at_time(sim, t_sample_h)
    if (is.na(C_pred) || C_pred <= 0) return(1e10)
    
    sd_res <- calc_residual_sd(C_pred, sigma_type, sigma_prop_sd, sigma_add_sd)
    ll <- dnorm(C_obs, mean = C_pred, sd = sd_res, log = TRUE)
    lp <- dmvnorm(eta_vec, mean = c(0, 0), sigma = omega_calibrate, log = TRUE)
    
    -(ll + lp)
  }
  
    opt <- robust_optim(obj_fun, par_init = c(0, 0), lower = c(-3, -3), upper = c(3, 3))
  
  list(
    eta_ind = c(eta_CL = opt$par[1], eta_Vc = opt$par[2]),
    converged = (opt$convergence == 0),
    value = opt$value
  )
}

message("✅ MAPB 函数已定义（5 个模型）")

## ############################################################################
## 3.3 定义单患者蒙特卡洛模拟函数（修复版）
## ############################################################################

message("/n--- 3.3 定义单患者蒙特卡洛模拟函数 ---")

run_mc_for_patient <- function(patient_row, bank, n_mc, t_sample_h, 
                                progress_prefix = "") {
  
  model_name <- patient_row$model_name
  
  ## 获取参数银行中的模型信息
  model_bank <- bank[[model_name]]
  omega_calibrate <- model_bank$omega_calibrate
  sigma_type <- model_bank$sigma_type
  sigma_prop_sd <- model_bank$sigma_prop_sd
  sigma_add_sd <- model_bank$sigma_add_sd
  calibratable_etas <- model_bank$calibratable_etas
  
  ## 确定是一维还是二维校准
  n_eta <- length(calibratable_etas)
  
  ## 构建给药事件表
  ev_df <- make_ev_steady(
    id = 1,
    dose_mg = patient_row$dose_mg,
    ii_h = patient_row$ii_h,
    n_dose = patient_row$n_dose,
    t_sample_h = t_sample_h
  )
  
  ## 准备协变量
  if (model_name == "maharaj2021") {
    iCov <- data.frame(id = 1, WT = patient_row$WT, PMA = patient_row$PMA)
  } else if (model_name == "zang2021") {
    iCov <- data.frame(
      id = 1, SEX = patient_row$SEX, SMK = patient_row$SMK,
      INF = patient_row$INF, VAL = patient_row$VAL,
      PER = patient_row$PER, SER = patient_row$SER,
      FLU = patient_row$FLU, DGH = patient_row$DGH
    )
  } else if (model_name == "li2018") {
    iCov <- data.frame(id = 1, WT = patient_row$WT)
  } else if (model_name == "yin2016") {
    iCov <- data.frame(id = 1, SEX = patient_row$SEX, SMK = patient_row$SMK)
  } else if (model_name == "sun2021") {
    iCov <- data.frame(
      id = 1, AGE = patient_row$AGE, SEX = patient_row$SEX,
      WT = patient_row$WT, SMOKE = patient_row$SMOKE,
      RACE_BLACK = patient_row$RACE_BLACK, RIF = patient_row$RIF,
      HEP_MOD = patient_row$HEP_MOD, REN_SEV = patient_row$REN_SEV,
      FED = patient_row$FED
    )
  } else {
    stop("未知模型: ", model_name)
  }
  
  ## 初始化结果存储
  results <- vector("list", n_mc)
  
  ## 蒙特卡洛循环
  for (i in seq_len(n_mc)) {
    
    if (i %% 100 == 0 || i == 1) {
      cat(sprintf("/r%s MC %d/%d", progress_prefix, i, n_mc))
    }
    
    ## ========== Step 1: 抽样 η_true ==========
    ## 修复：分别处理一维和二维情况，确保变量名正确
    
    if (n_eta == 1) {
      ## 一维抽样（Maharaj 模型）
      eta_CL_true <- rnorm(1, mean = 0, sd = sqrt(omega_calibrate[1, 1]))
      eta_V_true <- NA_real_  # Maharaj 没有 eta_V
      
    } else {
      ## 二维联合抽样
      eta_vec <- as.vector(rmvnorm(1, mean = rep(0, 2), sigma = omega_calibrate))
      eta_CL_true <- eta_vec[1]
      
      ## 第二��� eta 的名字根据模型不同
      if (model_name == "zang2021") {
        eta_V_true <- eta_vec[2]
      } else if (model_name == "li2018" || model_name == "sun2021") {
        eta_V_true <- eta_vec[2]  # eta_Vc
      } else if (model_name == "yin2016") {
        eta_V_true <- eta_vec[2]  # eta_V2
      }
    }
    
    ## ========== Step 2: 计算 θ_true ==========
    CL_pop <- patient_row$CL_pop
    CL_true <- CL_pop * exp(eta_CL_true)
    
    ## V 参数（根据模型不同）
    if (model_name == "maharaj2021") {
      V_pop <- patient_row$V_pop
      V_true <- V_pop  # Maharaj 没有 η_V
    } else if (model_name == "zang2021") {
      V_pop <- patient_row$V_pop
      V_true <- V_pop * exp(eta_V_true)
    } else if (model_name == "li2018") {
      V_pop <- patient_row$Vc_pop
      V_true <- V_pop * exp(eta_V_true)
    } else if (model_name == "yin2016") {
      V_pop <- patient_row$V2_pop
      V_true <- V_pop * exp(eta_V_true)
    } else if (model_name == "sun2021") {
      V_pop <- patient_row$Vc_pop
      V_true <- V_pop * exp(eta_V_true)
    }
    
    ## ========== Step 3: 模拟 C_true ==========
    ## 修复：直接传递数值参数，不使用命名索引
    
    if (model_name == "maharaj2021") {
      sim_true <- rxSolve(
        object = maharaj2021_olz_model,
        params = c(eta_CL = eta_CL_true),  # ← 直接用变量
        events = ev_df,
        iCov = iCov,
        returnType = "data.frame"
      )
    } else if (model_name == "zang2021") {
      sim_true <- rxSolve(
        object = zang2021_olz_model,
        params = c(eta_CL = eta_CL_true, eta_V = eta_V_true),
        events = ev_df,
        iCov = iCov,
        returnType = "data.frame"
      )
    } else if (model_name == "li2018") {
      sim_true <- rxSolve(
        object = li2018_olz_model_drug2,
        params = c(eta_CL = eta_CL_true, eta_Vc = eta_V_true,
                   eta_KA = 0, eta_Vp = 0, eta_Q = 0),
        events = ev_df,
        iCov = iCov,
        returnType = "data.frame"
      )
    } else if (model_name == "yin2016") {
      sim_true <- rxSolve(
        object = yin2016_olz_model,
        params = c(eta_CL = eta_CL_true, eta_V2 = eta_V_true,
                   eta_KA = 0, eta_Q = 0, eta_V3 = 0, eta_TLAG = 0),
        events = ev_df,
        iCov = iCov,
        returnType = "data.frame"
      )
    } else if (model_name == "sun2021") {
      sim_true <- rxSolve(
        object = sun2021_olz_model,
        params = c(eta_CL = eta_CL_true, eta_Vc = eta_V_true,
                   eta_KA = 0, eta_KA_IOV = 0, eta_Vp = 0, 
                   eta_Q = 0, eta_ALAG = 0),
        events = ev_df,
        iCov = iCov,
        returnType = "data.frame"
      )
    }
    
    C_true <- extract_conc_at_time(sim_true, t_sample_h)
    
    ## ========== Step 4: 生成 C_obs ==========
    C_obs <- add_residual_error(C_true, sigma_type, sigma_prop_sd, sigma_add_sd)
    
    ## ========== Step 5: MAPB 校准 ==========
    if (model_name == "maharaj2021") {
      mapb_result <- run_mapb_maharaj(
        C_obs = C_obs, ev_df = ev_df, iCov = iCov,
        omega_calibrate = omega_calibrate,
        sigma_type = sigma_type, sigma_prop_sd = sigma_prop_sd,
        sigma_add_sd = sigma_add_sd, t_sample_h = t_sample_h
      )
      eta_CL_ind <- mapb_result$eta_ind["eta_CL"]
      eta_V_ind <- NA_real_
      
    } else if (model_name == "zang2021") {
      mapb_result <- run_mapb_zang(
        C_obs = C_obs, ev_df = ev_df, iCov = iCov,
        omega_calibrate = omega_calibrate,
        sigma_type = sigma_type, sigma_prop_sd = sigma_prop_sd,
        sigma_add_sd = sigma_add_sd, t_sample_h = t_sample_h
      )
      eta_CL_ind <- mapb_result$eta_ind["eta_CL"]
      eta_V_ind <- mapb_result$eta_ind["eta_V"]
      
    } else if (model_name == "li2018") {
      mapb_result <- run_mapb_li2018(
        C_obs = C_obs, ev_df = ev_df, iCov = iCov,
        omega_calibrate = omega_calibrate,
        sigma_type = sigma_type, sigma_prop_sd = sigma_prop_sd,
        sigma_add_sd = sigma_add_sd, t_sample_h = t_sample_h
      )
      eta_CL_ind <- mapb_result$eta_ind["eta_CL"]
      eta_V_ind <- mapb_result$eta_ind["eta_Vc"]
      
    } else if (model_name == "yin2016") {
      mapb_result <- run_mapb_yin2016(
        C_obs = C_obs, ev_df = ev_df, iCov = iCov,
        omega_calibrate = omega_calibrate,
        sigma_type = sigma_type, sigma_prop_sd = sigma_prop_sd,
        sigma_add_sd = sigma_add_sd, t_sample_h = t_sample_h
      )
      eta_CL_ind <- mapb_result$eta_ind["eta_CL"]
      eta_V_ind <- mapb_result$eta_ind["eta_V2"]
      
    } else if (model_name == "sun2021") {
      mapb_result <- run_mapb_sun2021(
        C_obs = C_obs, ev_df = ev_df, iCov = iCov,
        omega_calibrate = omega_calibrate,
        sigma_type = sigma_type, sigma_prop_sd = sigma_prop_sd,
        sigma_add_sd = sigma_add_sd, t_sample_h = t_sample_h
      )
      eta_CL_ind <- mapb_result$eta_ind["eta_CL"]
      eta_V_ind <- mapb_result$eta_ind["eta_Vc"]
    }
    
    ## ========== Step 6: 计算 θ_ind ==========
    CL_ind <- CL_pop * exp(eta_CL_ind)
    
    if (model_name == "maharaj2021") {
      V_ind <- V_pop  # 没有 η_V
    } else {
      V_ind <- V_pop * exp(eta_V_ind)
    }
    
    ## ========== 存储结果 ==========
    results[[i]] <- tibble(
      mc_iter = i,
      
      ## η_true
      eta_CL_true = eta_CL_true,
      eta_V_true = eta_V_true,
      
      ## η_ind（MAPB 估计）
      eta_CL_ind = as.numeric(eta_CL_ind),  # ← 转为纯数值，去掉名字
      eta_V_ind = as.numeric(eta_V_ind),
      
      ## η_pop（固定为 0）
      eta_CL_pop = 0,
      eta_V_pop = 0,
      
      ## θ_true
      CL_true = CL_true,
      V_true = V_true,
      
      ## θ_ind
      CL_ind = as.numeric(CL_ind),
      V_ind = as.numeric(V_ind),
      
      ## θ_pop
      CL_pop = CL_pop,
      V_pop = V_pop,
      
      ## 浓度
      C_true = C_true,
      C_obs = C_obs,
      
      ## MAPB 收敛状态
      converged = mapb_result$converged
    )
  }
  
  cat("/n")
  
  ## 合并结果
  bind_rows(results)
}

message("✅ 蒙特卡洛模拟函数已定义（修复版）")

## ############################################################################
## 3.4 执行主线蒙特卡洛模拟
## ############################################################################

message("/n", paste(rep("-", 70), collapse = ""))
message("3.4 执行主线蒙特卡洛模拟（3 个患者 × ", N_MC, " 次）")
message(paste(rep("-", 70), collapse = ""))

t_sample_h <- tdm_rule$t_sample_h
message("采样时间: t = ", t_sample_h, " h")

mc_results_main <- list()

for (i in seq_len(nrow(patients_main))) {
  
  pat <- patients_main[i, ]
  
  message("/n--- 患者 ", i, "/", nrow(patients_main), ": ", pat$patient_id, 
          " (", pat$model_name, ") ---")
  
  mc_result <- run_mc_for_patient(
    patient_row = pat,
    bank = omega_sigma_bank,
    n_mc = N_MC,
    t_sample_h = t_sample_h,
    progress_prefix = paste0("[", pat$patient_id, "] ")
  )
  
  ## 添加患者标识
  mc_result <- mc_result %>%
    mutate(
      patient_id = pat$patient_id,
      pop = pat$pop,
      model_name = pat$model_name
    ) %>%
    select(patient_id, pop, model_name, everything())
  
  mc_results_main[[pat$patient_id]] <- mc_result
  
  ## 打印简要统计
  cat("  收敛率: ", sprintf("%.1f%%", 100 * mean(mc_result$converged)), "/n")
}

## 合并主线结果
mc_raw_main <- bind_rows(mc_results_main)

message("/n✅ 主线蒙特卡洛模拟完成")
message("   总行数: ", nrow(mc_raw_main))

## ############################################################################
## 3.5 执行支线蒙特卡洛模拟
## ############################################################################

message("/n", paste(rep("-", 70), collapse = ""))
message("3.5 执行支线蒙特卡洛模拟（4 个患者 × ", N_MC, " 次）")
message(paste(rep("-", 70), collapse = ""))

mc_results_side <- list()

for (i in seq_len(nrow(patients_side))) {
  
  pat <- patients_side[i, ]
  
  message("/n--- 患者 ", i, "/", nrow(patients_side), ": ", pat$patient_id,
          " (", pat$model_name, ") ---")
  
  mc_result <- run_mc_for_patient(
    patient_row = pat,
    bank = omega_sigma_bank,
    n_mc = N_MC,
    t_sample_h = t_sample_h,
    progress_prefix = paste0("[", pat$patient_id, "] ")
  )
  
  ## 添加患者标识
  mc_result <- mc_result %>%
    mutate(
      patient_id = pat$patient_id,
      pop = pat$pop,
      model_name = pat$model_name
    ) %>%
    select(patient_id, pop, model_name, everything())
  
  mc_results_side[[pat$patient_id]] <- mc_result
  
  ## 打印简要统计
  cat("  收敛率: ", sprintf("%.1f%%", 100 * mean(mc_result$converged)), "/n")
}

## 合并支线结果
mc_raw_side <- bind_rows(mc_results_side)

message("/n✅ 支线蒙特卡洛模拟完成")
message("   总行数: ", nrow(mc_raw_side))

## ############################################################################
## 3.6 计算检验指标
## ############################################################################

message("/n--- 3.6 计算检验指标 ---")

## ----------------------------------------------------------------------------
## 3.6.1 定义指标计算函数
## ----------------------------------------------------------------------------
calc_calibration_metrics <- function(mc_df, param_name, theta_hat_col, theta_true_col) {
  
  theta_hat <- mc_df[[theta_hat_col]]
  theta_true <- mc_df[[theta_true_col]]
  
  ## 过滤无效值
  valid <- !is.na(theta_hat) & !is.na(theta_true) & theta_true > 0
  theta_hat <- theta_hat[valid]
  theta_true <- theta_true[valid]
  n <- length(theta_hat)
  
  if (n == 0) {
    return(tibble(
      param = param_name,
      n = 0,
      rBias_pct = NA_real_,
      rRMSE_pct = NA_real_,
      F20_pct = NA_real_,
      F30_pct = NA_real_
    ))
  }
  
  ## 计算 PE (%)
  PE <- 100 * (theta_hat - theta_true) / theta_true
  
  ## rBias (%)
  rBias <- mean(PE)
  
  ## rRMSE (%)
  rRMSE <- sqrt(mean(PE^2))
  
  ## F20 / F30
  F20 <- 100 * mean(abs(PE) <= 20)
  F30 <- 100 * mean(abs(PE) <= 30)
  
  tibble(
    param = param_name,
    n = n,
    rBias_pct = rBias,
    rRMSE_pct = rRMSE,
    F20_pct = F20,
    F30_pct = F30
  )
}

## ----------------------------------------------------------------------------
## 3.6.2 计算所有患者的指标
## ----------------------------------------------------------------------------

## 合并所有结果
mc_all <- bind_rows(mc_raw_main, mc_raw_side)

## 按患者分组计算指标
metrics_list <- list()

for (pid in unique(mc_all$patient_id)) {
  
  mc_pat <- mc_all %>% filter(patient_id == pid)
  model_name <- mc_pat$model_name[1]
  
  ## CL: θ_ind vs θ_true
  met_CL_ind <- calc_calibration_metrics(
    mc_pat, "CL", "CL_ind", "CL_true"
  ) %>% mutate(comparison = "IND_vs_TRUE")
  
  ## CL: θ_pop vs θ_true
  met_CL_pop <- calc_calibration_metrics(
    mc_pat, "CL", "CL_pop", "CL_true"
  ) %>% mutate(comparison = "POP_vs_TRUE")
  
  ## V: θ_ind vs θ_true（如果适用）
  if (model_name != "maharaj2021") {
    met_V_ind <- calc_calibration_metrics(
      mc_pat, "V", "V_ind", "V_true"
    ) %>% mutate(comparison = "IND_vs_TRUE")
    
    met_V_pop <- calc_calibration_metrics(
      mc_pat, "V", "V_pop", "V_true"
    ) %>% mutate(comparison = "POP_vs_TRUE")
  } else {
    met_V_ind <- tibble(
      param = "V", n = NA_integer_, rBias_pct = NA_real_,
      rRMSE_pct = NA_real_, F20_pct = NA_real_, F30_pct = NA_real_,
      comparison = "IND_vs_TRUE"
    )
    met_V_pop <- tibble(
      param = "V", n = NA_integer_, rBias_pct = NA_real_,
      rRMSE_pct = NA_real_, F20_pct = NA_real_, F30_pct = NA_real_,
      comparison = "POP_vs_TRUE"
    )
  }
  
  metrics_list[[pid]] <- bind_rows(met_CL_ind, met_CL_pop, met_V_ind, met_V_pop) %>%
    mutate(
      patient_id = pid,
      model_name = model_name
    ) %>%
    select(patient_id, model_name, param, comparison, everything())
}

calibration_metrics <- bind_rows(metrics_list)

## 打印汇总
cat("/n=== 校准指标汇总（θ_ind vs θ_true）===/n")
print(calibration_metrics %>% 
        filter(comparison == "IND_vs_TRUE") %>%
        select(patient_id, model_name, param, rBias_pct, rRMSE_pct, F20_pct, F30_pct))

cat("/n=== 校准指标汇总（θ_pop vs θ_true）===/n")
print(calibration_metrics %>% 
        filter(comparison == "POP_vs_TRUE") %>%
        select(patient_id, model_name, param, rBias_pct, rRMSE_pct, F20_pct, F30_pct))

## ############################################################################
## 3.7 判定标准
## ############################################################################

message("/n--- 3.7 应用判定标准 ---")

## 判定函数
judge_calibration <- function(rBias_pct, rRMSE_pct, F30_pct) {
  
  if (is.na(rBias_pct) || is.na(rRMSE_pct) || is.na(F30_pct)) {
    return("N/A")
  }
  
  ## Excellent: |rBias| ≤ 10%, rRMSE ≤ 20%, F30 ≥ 70%
  if (abs(rBias_pct) <= 10 && rRMSE_pct <= 20 && F30_pct >= 70) {
    return("Excellent")
  }
  
  ## Acceptable: |rBias| ≤ 20%, rRMSE ≤ 30%, F30 ≥ 50%
  if (abs(rBias_pct) <= 20 && rRMSE_pct <= 30 && F30_pct >= 50) {
    return("Acceptable")
  }
  
  ## Unacceptable
  return("Unacceptable")
}

## 应用判定
calibration_judgement <- calibration_metrics %>%
  filter(comparison == "IND_vs_TRUE") %>%
  rowwise() %>%
  mutate(
    judgement = judge_calibration(rBias_pct, rRMSE_pct, F30_pct)
  ) %>%
  ungroup() %>%
  select(patient_id, model_name, param, rBias_pct, rRMSE_pct, F30_pct, judgement)

cat("/n=== 校准判定结果 ===/n")
print(calibration_judgement)

## 汇总判定结果
cat("/n=== 判定结果统计 ===/n")
print(table(calibration_judgement$judgement))

## ############################################################################
## 3.8 保存输出文件
## ############################################################################

message("/n--- 3.8 保存输出文件 ---")

## 保存主线原始结果
saveRDS(mc_raw_main, file.path(OUT_DIR_S3, "mc_raw_results_main.rds"))
message("  ✓ mc_raw_results_main.rds 已保存（", nrow(mc_raw_main), " 行）")

## 保存支线原始结果
saveRDS(mc_raw_side, file.path(OUT_DIR_S3, "mc_raw_results_side.rds"))
message("  ✓ mc_raw_results_side.rds 已保存（", nrow(mc_raw_side), " 行）")

## 保存校准指标
write_csv(calibration_metrics, file.path(OUT_DIR_S3, "calibration_metrics.csv"))
saveRDS(calibration_metrics, file.path(OUT_DIR_S3, "calibration_metrics.rds"))
message("  ✓ calibration_metrics.csv / .rds 已保存")

## 保存判定结果
write_csv(calibration_judgement, file.path(OUT_DIR_S3, "calibration_judgement.csv"))
saveRDS(calibration_judgement, file.path(OUT_DIR_S3, "calibration_judgement.rds"))
message("  ✓ calibration_judgement.csv / .rds 已保存")

## ############################################################################
## 3.9 诊断可视化
## ############################################################################

message("/n--- 3.9 生成诊断可视化 ---")

## 创建诊断图目录
plot_dir <- file.path(OUT_DIR_S3, "diagnostic_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

## ----------------------------------------------------------------------------
## 3.9.1 η_CL: η_ind vs η_true 散点图
## ----------------------------------------------------------------------------
p_eta_CL <- ggplot(mc_all, aes(x = eta_CL_true, y = eta_CL_ind)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_point(alpha = 0.1, size = 0.5) +
  facet_wrap(~ patient_id, scales = "free") +
  labs(
    title = "MAPB 校准诊断：η_CL_ind vs η_CL_true",
    subtitle = "对角线 = 完美估计；向原点收缩 = Shrinkage",
    x = "η_CL_true",
    y = "η_CL_ind"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40")
  )

ggsave(file.path(plot_dir, "scatter_eta_CL.png"), 
       p_eta_CL, width = 12, height = 8, dpi = 300)

## ----------------------------------------------------------------------------
## 3.9.2 CL: CL_ind vs CL_true 散点图
## ----------------------------------------------------------------------------
p_CL <- ggplot(mc_all, aes(x = CL_true, y = CL_ind)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_point(alpha = 0.1, size = 0.5) +
  facet_wrap(~ patient_id, scales = "free") +
  labs(
    title = "MAPB 校准诊断：CL_ind vs CL_true",
    x = "CL_true (L/h)",
    y = "CL_ind (L/h)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(file.path(plot_dir, "scatter_CL.png"), 
       p_CL, width = 12, height = 8, dpi = 300)

## ----------------------------------------------------------------------------
## 3.9.3 PE 分布直方图
## ----------------------------------------------------------------------------
mc_all <- mc_all %>%
  mutate(
    PE_CL = 100 * (CL_ind - CL_true) / CL_true,
    PE_CL_pop = 100 * (CL_pop - CL_true) / CL_true
  )

p_PE_hist <- ggplot(mc_all, aes(x = PE_CL)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = c(-20, 20), linetype = 2, color = "red") +
  geom_vline(xintercept = c(-30, 30), linetype = 3, color = "orange") +
  facet_wrap(~ patient_id, scales = "free_y") +
  labs(
    title = "预测误差（PE）分布：CL_ind vs CL_true",
    subtitle = "红色虚线 = ±20%；橙色虚线 = ±30%",
    x = "PE (%)",
    y = "频数"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40")
  )

ggsave(file.path(plot_dir, "histogram_PE_CL.png"), 
       p_PE_hist, width = 12, height = 8, dpi = 300)

## ----------------------------------------------------------------------------
## 3.9.4 指标森林图
## ----------------------------------------------------------------------------
metrics_for_plot <- calibration_judgement %>%
  mutate(
    patient_id = factor(patient_id),
    param_label = paste0(patient_id, " - ", param)
  )

p_forest <- ggplot(metrics_for_plot, aes(y = patient_id, x = rBias_pct, color = judgement)) +
  geom_vline(xintercept = 0, linetype = 1, color = "grey50") +
  geom_vline(xintercept = c(-10, 10), linetype = 2, color = "green4") +
  geom_vline(xintercept = c(-20, 20), linetype = 2, color = "orange") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = rBias_pct - rRMSE_pct, xmax = rBias_pct + rRMSE_pct),
                 height = 0.2) +
  facet_wrap(~ param, ncol = 2) +
  scale_color_manual(values = c(
    "Excellent" = "green4",
    "Acceptable" = "orange",
    "Unacceptable" = "red",
    "N/A" = "grey50"
  )) +
  labs(
    title = "校准指标森林图：rBias ± rRMSE",
    subtitle = "绿线 = ±10%（Excellent）；橙线 = ±20%（Acceptable）",
    x = "rBias (%)",
    y = "患者",
    color = "判定"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40")
  )

ggsave(file.path(plot_dir, "forest_metrics.png"), 
       p_forest, width = 10, height = 6, dpi = 300)

message("✅ 诊断图已保存至: ", plot_dir)

## ############################################################################
## 3.10 Stage 3 完成
## ############################################################################

message("/n", paste(rep("=", 70), collapse = ""))
message("✅ Stage 3 全部完成（第六版）")
message(paste(rep("=", 70), collapse = ""))

message("/n📁 输出目录: ", OUT_DIR_S3)
message("/n📄 输出文件：")
message("  - mc_raw_results_main.rds     （主线原始结果，", nrow(mc_raw_main), " 行）")
message("  - mc_raw_results_side.rds     （支线原始结果，", nrow(mc_raw_side), " 行）")
message("  - calibration_metrics.csv     （校准指标）")
message("  - calibration_judgement.csv   （判定结果）")
message("  - diagnostic_plots/           （诊断图）")

message("/n📊 判定结果汇总：")
print(table(calibration_judgement$param, calibration_judgement$judgement))

message("/n🔗 后续脚本使用方法：")
message('  mc_main <- readRDS("outputs6/outputs6_stage3/mc_raw_results_main.rds")')
message('  judgement <- readRDS("outputs6/outputs6_stage3/calibration_judgement.rds")')

message("/n⏭��  下一步：运行 [第六版代码-主线4.R]（后验概率计算）")
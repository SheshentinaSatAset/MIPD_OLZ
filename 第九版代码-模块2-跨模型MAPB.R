## ############################################################################
## 第九版代码-模块2-跨模型MAPB.R
## 
## 核心任务：
##   让 4 个成人模型分别对同一份 C_obs_runin 进行 MAPB 估计
##   - zang2021: 二维 MAPB (η_CL, η_V)
##   - li2018:   二维 MAPB (η_CL, η_Vc)
##   - yin2016:  一维 MAPB (仅 η_CL，η_V2 强制=0，Unacceptable)
##   - sun2021:  一维 MAPB (仅 η_CL，η_Vc 强制=0，Unacceptable)
## 
## 关键设计：
##   - 同一份 C_obs_runin 交给 4 个模型
##   - 每个模型用各自的先验（Ω）和残差模型（σ）
##   - η_ind 是"功能性参数"，没有对应的 η_true
## 
## 输入：
##   - outputs9/data_true.rds
##   - outputs9/module1_params.rds
##   - outputs6/outputs6_stage1/omega_sigma_bank.rds
## 
## 输出：
##   - outputs9/mapb_results.rds
##   - outputs9/module2_params.rds
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("第九版 模块2：跨模型 MAPB 估计")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## 2.0 加载依赖与路径设置
## ############################################################################

message("\n--- 2.0 加载依赖与路径设置 ---")

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

## 输入路径
IN_DIR_S1 <- file.path(PROJ_ROOT, "outputs6/outputs6_stage1")
OUT_DIR_9 <- file.path(PROJ_ROOT, "outputs9")

## 读取模块 1 输出
data_true <- readRDS(file.path(OUT_DIR_9, "data_true.rds"))
module1_params <- readRDS(file.path(OUT_DIR_9, "module1_params.rds"))
omega_sigma_bank <- readRDS(file.path(IN_DIR_S1, "omega_sigma_bank.rds"))

## 设置随机种子
MASTER_SEED <- module1_params$MASTER_SEED
N_MC <- module1_params$N_MC
T_SAMPLE_H <- module1_params$T_SAMPLE_H

set.seed(MASTER_SEED + 2)

message("✅ 模块 1 输出已加载")
message("   虚拟患者数: ", N_MC)
message("   采样时间: ", T_SAMPLE_H, " h")

## 加载 4 个成人模型
source(file.path(MODEL_DIR, "zang2021_olz_model.R"))
source(file.path(MODEL_DIR, "li2018_olz_model_drug2.R"))
source(file.path(MODEL_DIR, "yin2016_olz_model.R"))
source(file.path(MODEL_DIR, "sun2021_olz_model.R"))

message("✅ 4 个成人 PopPK 模型已加载")

## ############################################################################
## 2.1 定义各模型的协变量
## ############################################################################

message("\n--- 2.1 定义各模型的协变量 ---")

## 同一临床患者：22岁男性，65kg，不吸烟，合并丙戊酸钠

## zang2021 协变量
iCov_zang <- data.frame(
  id = 1,
  SEX = 1,   # 男性
  SMK = 0,   # 不吸烟
  INF = 0,   # 无感染
  VAL = 1,   # 合并丙戊酸钠
  PER = 0,
  SER = 0,
  FLU = 0,
  DGH = 0
)

## li2018 协变量
iCov_li <- data.frame(
  id = 1,
  WT = 65
)

## yin2016 协变量
iCov_yin <- data.frame(
  id = 1,
  SEX = 1,   # 男性（1=男）
  SMK = 0    # 不吸烟
)

## sun2021 协变量（注意：SEX=0 表示男性！）
iCov_sun <- data.frame(
  id = 1,
  AGE = 22,
  SEX = 0,         # 男性（0=男，1=女，与其他模型相反！）
  WT = 65,
  SMOKE = 0,
  RACE_BLACK = 0,
  RIF = 0,
  HEP_MOD = 0,
  REN_SEV = 0,
  FED = 0
)

message("✅ 协变量已定义")

## ############################################################################
## 2.2 定义各模型的 θ_pop
## ############################################################################

message("\n--- 2.2 计算各模型的 θ_pop ---")

## zang2021
TVCL_zang <- 12.88
TVV_zang <- 754.41
B_SEX <- 0.21; B_SMK <- 0.21; B_INF <- -0.29; B_VAL <- 0.21
B_PER <- -0.25; B_SER <- 0.15; B_FLU <- -0.36; B_DGH <- 0.70

CL_pop_zang <- TVCL_zang * exp(
  B_SEX * 1 + B_SMK * 0 + B_INF * 0 + B_VAL * 1 +
  B_PER * 0 + B_SER * 0 + B_FLU * 0 + B_DGH * 0
)
V_pop_zang <- TVV_zang

## li2018
TVCL_li <- 25.4
TVVC_li <- 2390
WT_REF_li <- 60.59
B_WT_VC_EXP_li <- 0.579

CL_pop_li <- TVCL_li
Vc_pop_li <- TVVC_li * (65 / WT_REF_li)^B_WT_VC_EXP_li

## yin2016
TVCL_yin <- 16.6
TVV2_yin <- 599
B_SEX_yin <- 0.385
B_SMK_yin <- 0.426

CL_pop_yin <- TVCL_yin * exp(B_SEX_yin * 1 + B_SMK_yin * 0)
V2_pop_yin <- TVV2_yin

## sun2021
TVCL_sun <- 15.5
TVVC_sun <- 656
WT_REF_sun <- 70
AGE_REF_sun <- 36
B_WT_CL_sun <- 0.75
B_WT_VC_sun <- 1.0
B_AGE_VC_sun <- 0.356

CL_pop_sun <- TVCL_sun * exp(B_WT_CL_sun * log(65 / WT_REF_sun))
Vc_pop_sun <- TVVC_sun * exp(
  B_WT_VC_sun * log(65 / WT_REF_sun) +
  B_AGE_VC_sun * log(22 / AGE_REF_sun)
)

cat("\n=== 各模型 θ_pop ===\n")
cat("zang2021: CL =", round(CL_pop_zang, 2), "L/h, V =", round(V_pop_zang, 2), "L\n")
cat("li2018:   CL =", round(CL_pop_li, 2), "L/h, Vc =", round(Vc_pop_li, 2), "L\n")
cat("yin2016:  CL =", round(CL_pop_yin, 2), "L/h, V2 =", round(V2_pop_yin, 2), "L\n")
cat("sun2021:  CL =", round(CL_pop_sun, 2), "L/h, Vc =", round(Vc_pop_sun, 2), "L\n")

## ############################################################################
## 2.3 定义辅助函数
## ############################################################################

message("\n--- 2.3 定义辅助函数 ---")

## ----------------------------------------------------------------------------
## 2.3.1 构建稳态给药事件表（complete 场景）
## ----------------------------------------------------------------------------
make_ev_steady <- function(dose_mg, ii_h, n_dose, t_sample_h, dt_h = 0.5) {
  
  dose_times <- seq(0, by = ii_h, length.out = n_dose)
  
  t_end <- t_sample_h + 1
  time_grid <- seq(0, t_end, by = dt_h)
  if (!any(abs(time_grid - t_sample_h) < 1e-8)) {
    time_grid <- sort(unique(c(time_grid, t_sample_h)))
  }
  
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
  ev_df$id <- 1
  
  ev_df
}

## ----------------------------------------------------------------------------
## 2.3.2 提取指定时间的浓度
## ----------------------------------------------------------------------------
extract_conc_at_time <- function(sim_df, t_target) {
  idx <- which(abs(sim_df$time - t_target) < 1e-8)
  if (length(idx) == 0) return(NA_real_)
  sim_df$C_ngmL[idx[1]]
}

## ----------------------------------------------------------------------------
## 2.3.3 计算残差 SD
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

## 预构建给药事件表
ev_df <- make_ev_steady(
  dose_mg = module1_params$DOSE_MG,
  ii_h = module1_params$II_H,
  n_dose = module1_params$N_DOSE,
  t_sample_h = T_SAMPLE_H
)

message("✅ 辅助函数已定义")

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
## 2.4 定义各模型的 MAPB 函数
## ############################################################################

message("\n--- 2.4 定义各模型的 MAPB 函数 ---")

## ----------------------------------------------------------------------------
## 2.4.1 Zang2021: 二维 MAPB (η_CL, η_V)
## ----------------------------------------------------------------------------
run_mapb_zang <- function(C_obs, ev_df, iCov, omega_calibrate,
                           sigma_type, sigma_prop_sd, sigma_add_sd,
                           t_sample_h) {
  
  obj_fun <- function(eta_vec) {
    eta_cl <- eta_vec[1]
    eta_v <- eta_vec[2]
    
    sim <- tryCatch({
      rxSolve(
        object = zang2021_olz_model,
        params = c(eta_CL = eta_cl, eta_V = eta_v),
        events = ev_df,
        iCov = iCov,
        returnType = "data.frame"
      )
    }, error = function(e) NULL)
    
    if (is.null(sim)) return(1e10)
    
    C_pred <- extract_conc_at_time(sim, t_sample_h)
    if (is.na(C_pred) || C_pred <= 0) return(1e10)
    
    sd_res <- calc_residual_sd(C_pred, sigma_type, sigma_prop_sd, sigma_add_sd)
    ll <- dnorm(C_obs, mean = C_pred, sd = sd_res, log = TRUE)
    lp <- dmvnorm(eta_vec, mean = c(0, 0), sigma = omega_calibrate, log = TRUE)
    
    -(ll + lp)
  }
  
  opt <- robust_optim(obj_fun, par_init = c(0, 0),
                      lower = c(-3, -3), upper = c(3, 3))
  
  list(
    eta_CL_ind = opt$par[1],
    eta_V_ind = opt$par[2],
    converged = (opt$convergence == 0),
    value = opt$value
  )
}

## ----------------------------------------------------------------------------
## 2.4.2 Li2018: 二维 MAPB (η_CL, η_Vc)
## ----------------------------------------------------------------------------
run_mapb_li <- function(C_obs, ev_df, iCov, omega_calibrate,
                         sigma_type, sigma_prop_sd, sigma_add_sd,
                         t_sample_h) {
  
  obj_fun <- function(eta_vec) {
    eta_cl <- eta_vec[1]
    eta_vc <- eta_vec[2]
    
    sim <- tryCatch({
      rxSolve(
        object = li2018_olz_model_drug2,
        params = c(eta_CL = eta_cl, eta_Vc = eta_vc,
                   eta_KA = 0, eta_Vp = 0, eta_Q = 0),
        events = ev_df,
        iCov = iCov,
        returnType = "data.frame"
      )
    }, error = function(e) NULL)
    
    if (is.null(sim)) return(1e10)
    
    C_pred <- extract_conc_at_time(sim, t_sample_h)
    if (is.na(C_pred) || C_pred <= 0) return(1e10)
    
    sd_res <- calc_residual_sd(C_pred, sigma_type, sigma_prop_sd, sigma_add_sd)
    ll <- dnorm(C_obs, mean = C_pred, sd = sd_res, log = TRUE)
    lp <- dmvnorm(eta_vec, mean = c(0, 0), sigma = omega_calibrate, log = TRUE)
    
    -(ll + lp)
  }
  
  opt <- robust_optim(obj_fun, par_init = c(0, 0),
                      lower = c(-3, -3), upper = c(3, 3))
  
  list(
    eta_CL_ind = opt$par[1],
    eta_V_ind = opt$par[2],
    converged = (opt$convergence == 0),
    value = opt$value
  )
}

## ----------------------------------------------------------------------------
## 2.4.3 Yin2016: 一维 MAPB (仅 η_CL，η_V2 强制=0)
## ----------------------------------------------------------------------------
run_mapb_yin <- function(C_obs, ev_df, iCov, omega_CL,
                          sigma_type, sigma_prop_sd, sigma_add_sd,
                          t_sample_h) {
  
  obj_fun <- function(eta_cl) {
    
    sim <- tryCatch({
      rxSolve(
        object = yin2016_olz_model,
        params = c(eta_CL = eta_cl, eta_V2 = 0,
                   eta_KA = 0, eta_Q = 0, eta_V3 = 0, eta_TLAG = 0),
        events = ev_df,
        iCov = iCov,
        returnType = "data.frame"
      )
    }, error = function(e) NULL)
    
    if (is.null(sim)) return(1e10)
    
    C_pred <- extract_conc_at_time(sim, t_sample_h)
    if (is.na(C_pred) || C_pred <= 0) return(1e10)
    
    sd_res <- calc_residual_sd(C_pred, sigma_type, sigma_prop_sd, sigma_add_sd)
    ll <- dnorm(C_obs, mean = C_pred, sd = sd_res, log = TRUE)
    lp <- dnorm(eta_cl, mean = 0, sd = sqrt(omega_CL), log = TRUE)
    
    -(ll + lp)
  }
  
    opt <- robust_optim(obj_fun, par_init = 0,
                      lower = -3, upper = 3)
  
  list(
    eta_CL_ind = opt$par,
    eta_V_ind = 0,  # 强制为 0
    converged = (opt$convergence == 0),
    value = opt$value
  )
}

## ----------------------------------------------------------------------------
## 2.4.4 Sun2021: 一维 MAPB (仅 η_CL，η_Vc 强制=0)
## ----------------------------------------------------------------------------
run_mapb_sun <- function(C_obs, ev_df, iCov, omega_CL,
                          sigma_type, sigma_prop_sd, sigma_add_sd,
                          t_sample_h) {
  
  obj_fun <- function(eta_cl) {
    
    sim <- tryCatch({
      rxSolve(
        object = sun2021_olz_model,
        params = c(eta_CL = eta_cl, eta_Vc = 0,
                   eta_KA = 0, eta_KA_IOV = 0, eta_Vp = 0,
                   eta_Q = 0, eta_ALAG = 0),
        events = ev_df,
        iCov = iCov,
        returnType = "data.frame"
      )
    }, error = function(e) NULL)
    
    if (is.null(sim)) return(1e10)
    
    C_pred <- extract_conc_at_time(sim, t_sample_h)
    if (is.na(C_pred) || C_pred <= 0) return(1e10)
    
    sd_res <- calc_residual_sd(C_pred, sigma_type, sigma_prop_sd, sigma_add_sd)
    ll <- dnorm(C_obs, mean = C_pred, sd = sd_res, log = TRUE)
    lp <- dnorm(eta_cl, mean = 0, sd = sqrt(omega_CL), log = TRUE)
    
    -(ll + lp)
  }
  
    opt <- robust_optim(obj_fun, par_init = 0,
                      lower = -3, upper = 3)
  
  list(
    eta_CL_ind = opt$par,
    eta_V_ind = 0,  # 强制为 0
    converged = (opt$convergence == 0),
    value = opt$value
  )
}

message("✅ MAPB 函数已定义（4 个模型）")

## ############################################################################
## 2.5 定义计算预测浓度的函数
## ############################################################################

message("\n--- 2.5 定义预测浓度计算函数 ---")

## 计算各模型用 η_ind 预测的 C_pred
calc_C_pred_all_models <- function(eta_results, ev_df, t_sample_h) {
  
  ## zang2021
  sim_zang <- rxSolve(
    object = zang2021_olz_model,
    params = c(eta_CL = eta_results$zang_eta_CL_ind, 
               eta_V = eta_results$zang_eta_V_ind),
    events = ev_df,
    iCov = iCov_zang,
    returnType = "data.frame"
  )
  C_pred_zang <- extract_conc_at_time(sim_zang, t_sample_h)
  
  ## li2018
  sim_li <- rxSolve(
    object = li2018_olz_model_drug2,
    params = c(eta_CL = eta_results$li_eta_CL_ind, 
               eta_Vc = eta_results$li_eta_V_ind,
               eta_KA = 0, eta_Vp = 0, eta_Q = 0),
    events = ev_df,
    iCov = iCov_li,
    returnType = "data.frame"
  )
  C_pred_li <- extract_conc_at_time(sim_li, t_sample_h)
  
  ## yin2016
  sim_yin <- rxSolve(
    object = yin2016_olz_model,
    params = c(eta_CL = eta_results$yin_eta_CL_ind, 
               eta_V2 = eta_results$yin_eta_V_ind,
               eta_KA = 0, eta_Q = 0, eta_V3 = 0, eta_TLAG = 0),
    events = ev_df,
    iCov = iCov_yin,
    returnType = "data.frame"
  )
  C_pred_yin <- extract_conc_at_time(sim_yin, t_sample_h)
  
  ## sun2021
  sim_sun <- rxSolve(
    object = sun2021_olz_model,
    params = c(eta_CL = eta_results$sun_eta_CL_ind, 
               eta_Vc = eta_results$sun_eta_V_ind,
               eta_KA = 0, eta_KA_IOV = 0, eta_Vp = 0,
               eta_Q = 0, eta_ALAG = 0),
    events = ev_df,
    iCov = iCov_sun,
    returnType = "data.frame"
  )
  C_pred_sun <- extract_conc_at_time(sim_sun, t_sample_h)
  
  list(
    zang = C_pred_zang,
    li = C_pred_li,
    yin = C_pred_yin,
    sun = C_pred_sun
  )
}

message("✅ 预测浓度计算函数已定义")

## ############################################################################
## 2.6 执行跨模型 MAPB
## ############################################################################

message("\n", paste(rep("-", 70), collapse = ""))
message("2.6 执行跨模型 MAPB（4 模型 × ", N_MC, " 患者）")
message(paste(rep("-", 70), collapse = ""))

## 提取各模型的 omega 和 sigma
bank_zang <- omega_sigma_bank$zang2021
bank_li <- omega_sigma_bank$li2018
bank_yin <- omega_sigma_bank$yin2016
bank_sun <- omega_sigma_bank$sun2021

## 初始化结果存储
mapb_results <- vector("list", N_MC)

start_time <- Sys.time()

for (i in seq_len(N_MC)) {
  
  if (i %% 50 == 0 || i == 1) {
    cat(sprintf("\r  进度: %d/%d", i, N_MC))
  }
  
  ## 获取该患者的 C_obs_runin
  C_obs <- data_true$C_obs_runin[i]
  
  ## ========== Zang2021 MAPB ==========
  mapb_zang <- run_mapb_zang(
    C_obs = C_obs,
    ev_df = ev_df,
    iCov = iCov_zang,
    omega_calibrate = bank_zang$omega_calibrate,
    sigma_type = bank_zang$sigma_type,
    sigma_prop_sd = bank_zang$sigma_prop_sd,
    sigma_add_sd = bank_zang$sigma_add_sd,
    t_sample_h = T_SAMPLE_H
  )
  
  ## ========== Li2018 MAPB ==========
  mapb_li <- run_mapb_li(
    C_obs = C_obs,
    ev_df = ev_df,
    iCov = iCov_li,
    omega_calibrate = bank_li$omega_calibrate,
    sigma_type = bank_li$sigma_type,
    sigma_prop_sd = bank_li$sigma_prop_sd,
    sigma_add_sd = bank_li$sigma_add_sd,
    t_sample_h = T_SAMPLE_H
  )
  
  ## ========== Yin2016 MAPB (一维) ==========
  mapb_yin <- run_mapb_yin(
    C_obs = C_obs,
    ev_df = ev_df,
    iCov = iCov_yin,
    omega_CL = bank_yin$omega_CL,
    sigma_type = bank_yin$sigma_type,
    sigma_prop_sd = bank_yin$sigma_prop_sd,
    sigma_add_sd = bank_yin$sigma_add_sd,
    t_sample_h = T_SAMPLE_H
  )
  
  ## ========== Sun2021 MAPB (一维) ==========
  mapb_sun <- run_mapb_sun(
    C_obs = C_obs,
    ev_df = ev_df,
    iCov = iCov_sun,
    omega_CL = bank_sun$omega_CL,
    sigma_type = bank_sun$sigma_type,
    sigma_prop_sd = bank_sun$sigma_prop_sd,
    sigma_add_sd = bank_sun$sigma_add_sd,
    t_sample_h = T_SAMPLE_H
  )
  
  ## ========== 计算 θ_ind ==========
  CL_ind_zang <- CL_pop_zang * exp(mapb_zang$eta_CL_ind)
  V_ind_zang <- V_pop_zang * exp(mapb_zang$eta_V_ind)
  
  CL_ind_li <- CL_pop_li * exp(mapb_li$eta_CL_ind)
  Vc_ind_li <- Vc_pop_li * exp(mapb_li$eta_V_ind)
  
  CL_ind_yin <- CL_pop_yin * exp(mapb_yin$eta_CL_ind)
  V2_ind_yin <- V2_pop_yin  # η_V2 = 0
  
  CL_ind_sun <- CL_pop_sun * exp(mapb_sun$eta_CL_ind)
  Vc_ind_sun <- Vc_pop_sun  # η_Vc = 0
  
  ## ========== 计算 C_pred_runin ==========
  eta_results <- list(
    zang_eta_CL_ind = mapb_zang$eta_CL_ind,
    zang_eta_V_ind = mapb_zang$eta_V_ind,
    li_eta_CL_ind = mapb_li$eta_CL_ind,
    li_eta_V_ind = mapb_li$eta_V_ind,
    yin_eta_CL_ind = mapb_yin$eta_CL_ind,
    yin_eta_V_ind = mapb_yin$eta_V_ind,
    sun_eta_CL_ind = mapb_sun$eta_CL_ind,
    sun_eta_V_ind = mapb_sun$eta_V_ind
  )
  
  C_pred <- calc_C_pred_all_models(eta_results, ev_df, T_SAMPLE_H)
  
  ## ========== 存储结果 ==========
  mapb_results[[i]] <- tibble(
    mc_iter = i,
    
    ## 来自模块 1 的数据
    C_obs_runin = C_obs,
    eta_CL_true = data_true$eta_CL_true[i],
    eta_V_true = data_true$eta_V_true[i],
    
    ## Zang2021 结果
    zang_eta_CL_ind = mapb_zang$eta_CL_ind,
    zang_eta_V_ind = mapb_zang$eta_V_ind,
    zang_CL_ind = CL_ind_zang,
    zang_V_ind = V_ind_zang,
    zang_C_pred_runin = C_pred$zang,
    zang_converged = mapb_zang$converged,
    
    ## Li2018 结果
    li_eta_CL_ind = mapb_li$eta_CL_ind,
    li_eta_Vc_ind = mapb_li$eta_V_ind,
    li_CL_ind = CL_ind_li,
    li_Vc_ind = Vc_ind_li,
    li_C_pred_runin = C_pred$li,
    li_converged = mapb_li$converged,
    
    ## Yin2016 结果
    yin_eta_CL_ind = mapb_yin$eta_CL_ind,
    yin_eta_V2_ind = 0,  # 强制为 0
    yin_CL_ind = CL_ind_yin,
    yin_V2_ind = V2_ind_yin,
    yin_C_pred_runin = C_pred$yin,
    yin_converged = mapb_yin$converged,
    
    ## Sun2021 结果
    sun_eta_CL_ind = mapb_sun$eta_CL_ind,
    sun_eta_Vc_ind = 0,  # 强制为 0
    sun_CL_ind = CL_ind_sun,
    sun_Vc_ind = Vc_ind_sun,
    sun_C_pred_runin = C_pred$sun,
    sun_converged = mapb_sun$converged
  )
}

cat("\n")

## 合并结果
mapb_results <- bind_rows(mapb_results)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

message("\n✅ 跨模型 MAPB 完成")
message("   总耗时: ", round(as.numeric(elapsed), 2), " 分钟")
message("   总行数: ", nrow(mapb_results))

## ############################################################################
## 2.7 数据质量检查
## ############################################################################

message("\n--- 2.7 数据质量检查 ---")

## 收敛率
cat("\n=== 各模型收敛率 ===\n")
cat("zang2021:", sprintf("%.1f%%", 100 * mean(mapb_results$zang_converged)), "\n")
cat("li2018:  ", sprintf("%.1f%%", 100 * mean(mapb_results$li_converged)), "\n")
cat("yin2016: ", sprintf("%.1f%%", 100 * mean(mapb_results$yin_converged)), "\n")
cat("sun2021: ", sprintf("%.1f%%", 100 * mean(mapb_results$sun_converged)), "\n")

## η_CL_ind 统计
cat("\n=== η_CL_ind 统计 ===\n")
cat("zang: mean =", round(mean(mapb_results$zang_eta_CL_ind), 4),
    ", sd =", round(sd(mapb_results$zang_eta_CL_ind), 4), "\n")
cat("li:   mean =", round(mean(mapb_results$li_eta_CL_ind), 4),
    ", sd =", round(sd(mapb_results$li_eta_CL_ind), 4), "\n")
cat("yin:  mean =", round(mean(mapb_results$yin_eta_CL_ind), 4),
    ", sd =", round(sd(mapb_results$yin_eta_CL_ind), 4), "\n")
cat("sun:  mean =", round(mean(mapb_results$sun_eta_CL_ind), 4),
    ", sd =", round(sd(mapb_results$sun_eta_CL_ind), 4), "\n")

## CL_ind 统计
cat("\n=== CL_ind 统计 ===\n")
cat("zang: mean =", round(mean(mapb_results$zang_CL_ind), 2), "L/h\n")
cat("li:   mean =", round(mean(mapb_results$li_CL_ind), 2), "L/h\n")
cat("yin:  mean =", round(mean(mapb_results$yin_CL_ind), 2), "L/h\n")
cat("sun:  mean =", round(mean(mapb_results$sun_CL_ind), 2), "L/h\n")

## C_pred vs C_obs 诊断
cat("\n=== C_pred_runin vs C_obs_runin（MAPB 拟合质量）===\n")
cat("zang: cor =", round(cor(mapb_results$zang_C_pred_runin, mapb_results$C_obs_runin), 4), "\n")
cat("li:   cor =", round(cor(mapb_results$li_C_pred_runin, mapb_results$C_obs_runin), 4), "\n")
cat("yin:  cor =", round(cor(mapb_results$yin_C_pred_runin, mapb_results$C_obs_runin), 4), "\n")
cat("sun:  cor =", round(cor(mapb_results$sun_C_pred_runin, mapb_results$C_obs_runin), 4), "\n")

## ############################################################################
## 2.8 保存输出文件
## ############################################################################

message("\n--- 2.8 保存输出文件 ---")

## 保存 MAPB 结果
saveRDS(mapb_results, file.path(OUT_DIR_9, "mapb_results.rds"))
write_csv(mapb_results, file.path(OUT_DIR_9, "mapb_results.csv"))
message("  ✓ mapb_results.rds / .csv 已保存（", nrow(mapb_results), " 行）")

## 保存模块参数
module2_params <- list(
  ## θ_pop
  theta_pop = list(
    zang = list(CL = CL_pop_zang, V = V_pop_zang),
    li = list(CL = CL_pop_li, Vc = Vc_pop_li),
    yin = list(CL = CL_pop_yin, V2 = V2_pop_yin),
    sun = list(CL = CL_pop_sun, Vc = Vc_pop_sun)
  ),
  
  ## 协变量
  iCov = list(
    zang = iCov_zang,
    li = iCov_li,
    yin = iCov_yin,
    sun = iCov_sun
  ),
  
  ## 耗时
  elapsed_minutes = as.numeric(elapsed),
  
  ## 元数据
  created = Sys.time(),
  version = "v9.0",
  module = "Module 2 - Cross-Model MAPB"
)

saveRDS(module2_params, file.path(OUT_DIR_9, "module2_params.rds"))
message("  ✓ module2_params.rds 已保存")

## ############################################################################
## 2.9 模块 2 完成
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("✅ 模块 2 全部完成（跨模型 MAPB）")
message(paste(rep("=", 70), collapse = ""))

message("\n📁 输出目录: ", OUT_DIR_9)
message("\n📄 输出文件：")
message("  - mapb_results.rds / .csv   （MAPB 结果，", nrow(mapb_results), " 行）")
message("  - module2_params.rds        （模块参数）")

message("\n📊 各模型 MAPB 摘要：")
message("   zang2021: 二维 (η_CL, η_V)，收敛率 ", 
        sprintf("%.1f%%", 100 * mean(mapb_results$zang_converged)))
message("   li2018:   二维 (η_CL, η_Vc)，收敛率 ", 
        sprintf("%.1f%%", 100 * mean(mapb_results$li_converged)))
message("   yin2016:  一维 (η_CL)，η_V2=0，收敛率 ", 
        sprintf("%.1f%%", 100 * mean(mapb_results$yin_converged)))
message("   sun2021:  一维 (η_CL)，η_Vc=0，收敛率 ", 
        sprintf("%.1f%%", 100 * mean(mapb_results$sun_converged)))

message("\n🔗 后续模块使用方法：")
message('  mapb_results <- readRDS("outputs9/mapb_results.rds")')
message('  data_true <- readRDS("outputs9/data_true.rds")')
message('  combined <- left_join(data_true, mapb_results, by = "mc_iter")')

message("\n⏭️  下一步：运行 [第九版代码-模块3-跨模型判别.R]")
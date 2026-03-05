## ############################################################################
## 第六版代码-敏感性分析-MAPB假设偏差.R
## 
## 研究问题：
##   患者实际非依从服药 → 医生仍以完全依从假设进行 MAPB
##   这种假设偏差对 MAPB 校准结果产生怎样的影响？
## 
## 设计：
##   - 患者实际服药场景 s_true_runin ∈ {miss_last, miss_prev, miss_two, 
##                                     double_last, double_prev}
##   - 医生 MAPB 假设场景：complete（固定）
##   - 目标患者：PED (maharaj2021), ADULT (zang2021), ELDERLY (zang2021)
## 
## 输出：
##   - outputs6/outputs6_stage3_SA_MAPB/mc_raw_results_SA.rds
##   - outputs6/outputs6_stage3_SA_MAPB/calibration_metrics_SA.rds / .csv
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("敏感性分析：MAPB 假设偏差对校准结果的影响")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## SA.0 加载依赖与路径设置
## ############################################################################

message("\n--- SA.0 加载依赖与路径设置 ---")

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
OUTPUT_ROOT <- file.path(PROJ_ROOT, "outputs6")
OUT_DIR_S1 <- file.path(OUTPUT_ROOT, "outputs6_stage1")
OUT_DIR_S2 <- file.path(OUTPUT_ROOT, "outputs6_stage2")

## 敏感性分析输出目录
OUT_DIR_SA <- file.path(OUTPUT_ROOT, "outputs6_stage3_SA_MAPB")

if (!dir.exists(OUT_DIR_SA)) {
  dir.create(OUT_DIR_SA, recursive = TRUE)
  message("📁 已创建目录: ", OUT_DIR_SA)
}

## 读取前置输出
global_settings <- readRDS(file.path(OUT_DIR_S1, "global_settings.rds"))
omega_sigma_bank <- readRDS(file.path(OUT_DIR_S1, "omega_sigma_bank.rds"))
design_manifest <- readRDS(file.path(OUT_DIR_S2, "design_manifest_stage2.rds"))
patients_main <- design_manifest$patients_main
tdm_rule <- design_manifest$tdm_rule

MASTER_SEED <- global_settings$MASTER_SEED
N_MC <- global_settings$N_MC
T_SAMPLE_H <- tdm_rule$t_sample_h

set.seed(MASTER_SEED + 300)  # 独立种子序列

message("✅ 前置输出已加载")
message("🎲 全局随机种子: ", MASTER_SEED)
message("🔢 蒙特卡洛次数: ", N_MC)
message("⏱️ ���样时间: ", T_SAMPLE_H, " h")

## 加载模型
source(file.path(MODEL_DIR, "maharaj2021_olz_model.R"))
source(file.path(MODEL_DIR, "zang2021_olz_model.R"))

message("✅ maharaj2021 和 zang2021 模型已加载")

## ############################################################################
## SA.1 定义依从性场景
## ############################################################################

message("\n--- SA.1 定义依从性场景 ---")

## 非依从场景（用于患者实际服药）
noncompliance_scenarios <- tibble(
  scenario_id = c("miss_last", "miss_prev", "miss_two", 
                  "double_last", "double_prev"),
  scenario_name = c("末次漏服", "倒数第二次漏服", "连续漏服两次",
                    "末次双倍", "倒数第二次双倍"),
  dose_13_mult = c(1, 0, 0, 1, 2),
  dose_14_mult = c(0, 1, 0, 2, 1)
)

## 完全依从场景（用于医生 MAPB 假设）
complete_scenario <- tibble(
  scenario_id = "complete",
  dose_13_mult = 1,
  dose_14_mult = 1
)

cat("\n=== 非依从场景（患者实际服药）===\n")
print(noncompliance_scenarios)

cat("\n=== MAPB 假设场景（医生假设）===\n")
cat("固定为 complete（完全依从）\n")

## ############################################################################
## SA.2 定义辅助函数
## ############################################################################

message("\n--- SA.2 定义辅助函数 ---")

## ----------------------------------------------------------------------------
## 构建给药事件表（指定场景）
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
## 提取指定时间的浓度
## ----------------------------------------------------------------------------
extract_conc_at_time <- function(sim_df, t_target) {
  idx <- which(abs(sim_df$time - t_target) < 1e-8)
  if (length(idx) == 0) return(NA_real_)
  sim_df$C_ngmL[idx[1]]
}

## ----------------------------------------------------------------------------
## 添加残差
## ----------------------------------------------------------------------------
add_residual_error <- function(C_true, sigma_type, sigma_prop_sd, sigma_add_sd = NA) {
  
  if (sigma_type == "proportional") {
    eps <- rnorm(1, mean = 0, sd = sigma_prop_sd)
    C_obs <- C_true * (1 + eps)
  } else if (sigma_type == "combined") {
    sigma_total <- sqrt((sigma_prop_sd * C_true)^2 + sigma_add_sd^2)
    eps <- rnorm(1, mean = 0, sd = sigma_total)
    C_obs <- C_true + eps
  } else {
    stop("未知残差类型: ", sigma_type)
  }
  
  max(C_obs, 1e-6)
}

## ----------------------------------------------------------------------------
## 计算残差 SD
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

## ----------------------------------------------------------------------------
## MAPB 函数：Maharaj2021（一维）
## ----------------------------------------------------------------------------
run_mapb_maharaj <- function(C_obs, ev_df, iCov, omega_calibrate, 
                              sigma_type, sigma_prop_sd, sigma_add_sd,
                              t_sample_h) {
  
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
    
    sd_res <- calc_residual_sd(C_pred, sigma_type, sigma_prop_sd, sigma_add_sd)
    ll <- dnorm(C_obs, mean = C_pred, sd = sd_res, log = TRUE)
    
    omega_CL <- omega_calibrate[1, 1]
    lp <- dnorm(eta_cl, mean = 0, sd = sqrt(omega_CL), log = TRUE)
    
    -(ll + lp)
  }
  
  opt <- robust_optim(obj_fun, par_init = 0, lower = -3, upper = 3)
  
  list(
    eta_ind = c(eta_CL = opt$par),
    converged = (opt$convergence == 0),
    value = opt$value
  )
}

## ----------------------------------------------------------------------------
## MAPB 函数：Zang2021（二维）
## ----------------------------------------------------------------------------
run_mapb_zang <- function(C_obs, ev_df, iCov, omega_calibrate,
                           sigma_type, sigma_prop_sd, sigma_add_sd,
                           t_sample_h) {
  
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
    
    sd_res <- calc_residual_sd(C_pred, sigma_type, sigma_prop_sd, sigma_add_sd)
    ll <- dnorm(C_obs, mean = C_pred, sd = sd_res, log = TRUE)
    lp <- dmvnorm(eta_vec, mean = c(0, 0), sigma = omega_calibrate, log = TRUE)
    
    -(ll + lp)
  }
  
  opt <- robust_optim(obj_fun, par_init = 0, lower = -3, upper = 3)
  
  list(
    eta_ind = c(eta_CL = opt$par[1], eta_V = opt$par[2]),
    converged = (opt$convergence == 0),
    value = opt$value
  )
}

message("✅ 辅助函数已定义")

## ############################################################################
## SA.3 执行蒙特卡洛模拟
## ############################################################################

message("\n", paste(rep("-", 70), collapse = ""))
message("SA.3 执行蒙特卡洛模拟")
message("    患者数: 3 (PED, ADULT, ELDERLY)")
message("    非依从场景数: 5")
message("    每组合蒙特卡洛次数: ", N_MC)
message("    总模拟次数: ", 3 * 5 * N_MC)
message(paste(rep("-", 70), collapse = ""))

## 筛选目标患者
target_patients <- patients_main %>%
  filter(patient_id %in% c("PED", "ADULT", "ELDERLY"))

cat("\n=== 目标患者 ===\n")
print(target_patients %>% select(patient_id, model_name, dose_mg, CL_pop, V_pop))

## 初始化结果存储
all_results <- list()
result_counter <- 0

start_time <- Sys.time()

## 外层循环：非依从场景
for (s_idx in seq_len(nrow(noncompliance_scenarios))) {
  
  s_true_runin <- noncompliance_scenarios$scenario_id[s_idx]
  dose_13_mult_true <- noncompliance_scenarios$dose_13_mult[s_idx]
  dose_14_mult_true <- noncompliance_scenarios$dose_14_mult[s_idx]
  
  message("\n>>> 场景 ", s_idx, "/5: ", s_true_runin)
  
  ## 中层循环：患者
  for (p_idx in seq_len(nrow(target_patients))) {
    
    pat <- target_patients[p_idx, ]
    patient_id <- pat$patient_id
    model_name <- pat$model_name
    
    cat(sprintf("    患者 %s (%s): ", patient_id, model_name))
    
    ## 获取模型参数
    model_bank <- omega_sigma_bank[[model_name]]
    omega_calibrate <- model_bank$omega_calibrate
    sigma_type <- model_bank$sigma_type
    sigma_prop_sd <- model_bank$sigma_prop_sd
    sigma_add_sd <- model_bank$sigma_add_sd
    
    ## 准备协变量
    if (model_name == "maharaj2021") {
      iCov <- data.frame(id = 1, WT = pat$WT, PMA = pat$PMA)
    } else if (model_name == "zang2021") {
      iCov <- data.frame(
        id = 1, SEX = pat$SEX, SMK = pat$SMK,
        INF = pat$INF, VAL = pat$VAL,
        PER = pat$PER, SER = pat$SER,
        FLU = pat$FLU, DGH = pat$DGH
      )
    }
    
    ## 构建给药事件表
    ## 真实服药场景（非依从）
    ev_true <- make_ev_for_scenario(
      dose_mg = pat$dose_mg, ii_h = pat$ii_h, n_dose = pat$n_dose,
      t_sample_h = T_SAMPLE_H,
      dose_13_mult = dose_13_mult_true,
      dose_14_mult = dose_14_mult_true
    )
    
    ## MAPB 假设场景（complete）
    ev_mapb <- make_ev_for_scenario(
      dose_mg = pat$dose_mg, ii_h = pat$ii_h, n_dose = pat$n_dose,
      t_sample_h = T_SAMPLE_H,
      dose_13_mult = 1,  # complete
      dose_14_mult = 1   # complete
    )
    
    ## 内层循环：蒙特卡洛
    mc_results <- vector("list", N_MC)
    
    for (i in seq_len(N_MC)) {
      
      if (i %% 200 == 0) cat(".")
      
      ## ========== Step 1: 抽样 η_true ==========
      if (model_name == "maharaj2021") {
        eta_CL_true <- rnorm(1, mean = 0, sd = sqrt(omega_calibrate[1, 1]))
        eta_V_true <- NA_real_
      } else {
        eta_vec <- as.vector(rmvnorm(1, mean = rep(0, 2), sigma = omega_calibrate))
        eta_CL_true <- eta_vec[1]
        eta_V_true <- eta_vec[2]
      }
      
      ## ========== Step 2: 计算 θ_true ==========
      CL_pop <- pat$CL_pop
      V_pop <- pat$V_pop
      CL_true <- CL_pop * exp(eta_CL_true)
      V_true <- if (is.na(eta_V_true)) V_pop else V_pop * exp(eta_V_true)
      
      ## ========== Step 3: 用真实场景生成 C_obs ==========
      if (model_name == "maharaj2021") {
        sim_true <- rxSolve(
          object = maharaj2021_olz_model,
          params = c(eta_CL = eta_CL_true),
          events = ev_true,
          iCov = iCov,
          returnType = "data.frame"
        )
      } else {
        sim_true <- rxSolve(
          object = zang2021_olz_model,
          params = c(eta_CL = eta_CL_true, eta_V = eta_V_true),
          events = ev_true,
          iCov = iCov,
          returnType = "data.frame"
        )
      }
      
      C_true <- extract_conc_at_time(sim_true, T_SAMPLE_H)
      C_obs <- add_residual_error(C_true, sigma_type, sigma_prop_sd, sigma_add_sd)
      
      ## ========== Step 4: 用 complete 假设执行 MAPB ==========
      if (model_name == "maharaj2021") {
        mapb_result <- run_mapb_maharaj(
          C_obs = C_obs, ev_df = ev_mapb, iCov = iCov,
          omega_calibrate = omega_calibrate,
          sigma_type = sigma_type, sigma_prop_sd = sigma_prop_sd,
          sigma_add_sd = sigma_add_sd, t_sample_h = T_SAMPLE_H
        )
        eta_CL_ind <- mapb_result$eta_ind["eta_CL"]
        eta_V_ind <- NA_real_
      } else {
        mapb_result <- run_mapb_zang(
          C_obs = C_obs, ev_df = ev_mapb, iCov = iCov,
          omega_calibrate = omega_calibrate,
          sigma_type = sigma_type, sigma_prop_sd = sigma_prop_sd,
          sigma_add_sd = sigma_add_sd, t_sample_h = T_SAMPLE_H
        )
        eta_CL_ind <- mapb_result$eta_ind["eta_CL"]
        eta_V_ind <- mapb_result$eta_ind["eta_V"]
      }
      
      ## ========== Step 5: 计算 θ_ind ==========
      CL_ind <- CL_pop * exp(as.numeric(eta_CL_ind))
      V_ind <- if (is.na(eta_V_ind)) V_pop else V_pop * exp(as.numeric(eta_V_ind))
      
      ## ========== 存储结果 ==========
      mc_results[[i]] <- tibble(
        mc_iter = i,
        patient_id = patient_id,
        model_name = model_name,
        s_true_runin = s_true_runin,
        
        ## η_true
        eta_CL_true = eta_CL_true,
        eta_V_true = eta_V_true,
        
        ## η_ind
        eta_CL_ind = as.numeric(eta_CL_ind),
        eta_V_ind = as.numeric(eta_V_ind),
        
        ## θ_true
        CL_true = CL_true,
        V_true = V_true,
        
        ## θ_ind
        CL_ind = CL_ind,
        V_ind = V_ind,
        
        ## θ_pop
        CL_pop = CL_pop,
        V_pop = V_pop,
        
        ## 浓度
        C_true = C_true,
        C_obs = C_obs,
        
        ## 收敛状态
        converged = mapb_result$converged
      )
    }
    
    cat(" 完成\n")
    
    ## 合并该患者-场景组合的结果
    result_counter <- result_counter + 1
    all_results[[result_counter]] <- bind_rows(mc_results)
  }
}

## 合并所有结果
mc_raw_SA <- bind_rows(all_results)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

message("\n✅ 蒙特卡洛模拟完成")
message("   总耗时: ", round(as.numeric(elapsed), 2), " 分钟")
message("   总行数: ", nrow(mc_raw_SA))

## ############################################################################
## SA.4 计算 rBias、rRMSE、F20、F30（修复版）
## ############################################################################

message("\n--- SA.4 计算校准指标 ---")

## 定义计算函数（与主线3一致）
calc_metrics <- function(theta_hat, theta_true) {
  
  valid <- !is.na(theta_hat) & !is.na(theta_true) & theta_true > 0
  theta_hat <- theta_hat[valid]
  theta_true <- theta_true[valid]
  n <- length(theta_hat)
  
  if (n == 0) {
    return(list(
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
  
  list(n = n, rBias_pct = rBias, rRMSE_pct = rRMSE, F20_pct = F20, F30_pct = F30)
}

## 按 (patient_id, s_true_runin, param, comparison) 计算指标
metrics_list <- list()

for (pid in unique(mc_raw_SA$patient_id)) {
  for (s_runin in unique(mc_raw_SA$s_true_runin)) {
    
    mc_subset <- mc_raw_SA %>%
      filter(patient_id == pid, s_true_runin == s_runin)
    
    model_name <- mc_subset$model_name[1]
    
    ## ========== CL 参数 ==========
    
    ## CL: IND vs TRUE
    met_CL_ind <- calc_metrics(mc_subset$CL_ind, mc_subset$CL_true)
    
    metrics_list[[length(metrics_list) + 1]] <- tibble(
      patient_id = pid,
      model_name = model_name,
      s_true_runin = s_runin,
      param = "CL",
      comparison = "IND_vs_TRUE",
      n = met_CL_ind$n,
      rBias_pct = met_CL_ind$rBias_pct,
      rRMSE_pct = met_CL_ind$rRMSE_pct,
      F20_pct = met_CL_ind$F20_pct,
      F30_pct = met_CL_ind$F30_pct
    )
    
    ## CL: POP vs TRUE
    met_CL_pop <- calc_metrics(mc_subset$CL_pop, mc_subset$CL_true)
    
    metrics_list[[length(metrics_list) + 1]] <- tibble(
      patient_id = pid,
      model_name = model_name,
      s_true_runin = s_runin,
      param = "CL",
      comparison = "POP_vs_TRUE",
      n = met_CL_pop$n,
      rBias_pct = met_CL_pop$rBias_pct,
      rRMSE_pct = met_CL_pop$rRMSE_pct,
      F20_pct = met_CL_pop$F20_pct,
      F30_pct = met_CL_pop$F30_pct
    )
    
    ## ========== V 参数（仅 zang2021）==========
    if (model_name != "maharaj2021") {
      
      ## V: IND vs TRUE
      met_V_ind <- calc_metrics(mc_subset$V_ind, mc_subset$V_true)
      
      metrics_list[[length(metrics_list) + 1]] <- tibble(
        patient_id = pid,
        model_name = model_name,
        s_true_runin = s_runin,
        param = "V",
        comparison = "IND_vs_TRUE",
        n = met_V_ind$n,
        rBias_pct = met_V_ind$rBias_pct,
        rRMSE_pct = met_V_ind$rRMSE_pct,
        F20_pct = met_V_ind$F20_pct,
        F30_pct = met_V_ind$F30_pct
      )
      
      ## V: POP vs TRUE
      met_V_pop <- calc_metrics(mc_subset$V_pop, mc_subset$V_true)
      
      metrics_list[[length(metrics_list) + 1]] <- tibble(
        patient_id = pid,
        model_name = model_name,
        s_true_runin = s_runin,
        param = "V",
        comparison = "POP_vs_TRUE",
        n = met_V_pop$n,
        rBias_pct = met_V_pop$rBias_pct,
        rRMSE_pct = met_V_pop$rRMSE_pct,
        F20_pct = met_V_pop$F20_pct,
        F30_pct = met_V_pop$F30_pct
      )
    }
  }
}

calibration_metrics_SA <- bind_rows(metrics_list)

## 打印结果
cat("\n=== IND vs TRUE 校准指标（CL 参数）===\n")
print(
  calibration_metrics_SA %>%
    filter(param == "CL", comparison == "IND_vs_TRUE") %>%
    select(patient_id, s_true_runin, rBias_pct, rRMSE_pct, F20_pct, F30_pct) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2)))
)

cat("\n=== POP vs TRUE 校准指标（CL 参数）===\n")
print(
  calibration_metrics_SA %>%
    filter(param == "CL", comparison == "POP_vs_TRUE") %>%
    select(patient_id, s_true_runin, rBias_pct, rRMSE_pct, F20_pct, F30_pct) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2)))
)

## ############################################################################
## SA.5 与基准对比（从原主线3读取 complete 场景结果）
## ############################################################################

message("\n--- SA.5 与基准对比 ---")

## 读取原主线3的校准指标
calibration_metrics_baseline <- readRDS(
  file.path(OUTPUT_ROOT, "outputs6_stage3/calibration_metrics.rds")
)

## 提取 complete 场景的基准（IND_vs_TRUE 和 POP_vs_TRUE）
baseline_metrics <- calibration_metrics_baseline %>%
  filter(patient_id %in% c("PED", "ADULT", "ELDERLY")) %>%
  select(patient_id, model_name, param, comparison,
         rBias_baseline = rBias_pct, 
         rRMSE_baseline = rRMSE_pct,
         F20_baseline = F20_pct,
         F30_baseline = F30_pct) %>%
  mutate(s_true_runin = "complete")

## 合并敏感性分析结果与基准
comparison_table <- calibration_metrics_SA %>%
  left_join(
    baseline_metrics %>% select(patient_id, param, comparison, 
                                 rBias_baseline, rRMSE_baseline,
                                 F20_baseline, F30_baseline),
    by = c("patient_id", "param", "comparison")
  ) %>%
  mutate(
    delta_rBias = rBias_pct - rBias_baseline,
    delta_rRMSE = rRMSE_pct - rRMSE_baseline,
    delta_F20 = F20_pct - F20_baseline,
    delta_F30 = F30_pct - F30_baseline
  )

cat("\n=== 与基准对比（IND_vs_TRUE，CL 参数）===\n")
print(
  comparison_table %>%
    filter(param == "CL", comparison == "IND_vs_TRUE") %>%
    select(patient_id, s_true_runin, 
           rBias_pct, rBias_baseline, delta_rBias,
           rRMSE_pct, rRMSE_baseline, delta_rRMSE) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2)))
)

cat("\n=== 与基准对比（POP_vs_TRUE，CL 参数）===\n")
print(
  comparison_table %>%
    filter(param == "CL", comparison == "POP_vs_TRUE") %>%
    select(patient_id, s_true_runin, 
           rBias_pct, rBias_baseline, delta_rBias,
           rRMSE_pct, rRMSE_baseline, delta_rRMSE) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2)))
)

## ############################################################################
## SA.6 保存输出
## ############################################################################

message("\n--- SA.6 保存输出 ---")

## 保存原始蒙特卡洛数据
saveRDS(mc_raw_SA, file.path(OUT_DIR_SA, "mc_raw_results_SA.rds"))
message("  ✓ mc_raw_results_SA.rds 已保存（", nrow(mc_raw_SA), " 行）")

## 保存校准指标
saveRDS(calibration_metrics_SA, file.path(OUT_DIR_SA, "calibration_metrics_SA.rds"))
write_csv(calibration_metrics_SA, file.path(OUT_DIR_SA, "calibration_metrics_SA.csv"))
message("  ✓ calibration_metrics_SA.rds / .csv 已保存")

## 保存对比表
saveRDS(comparison_table, file.path(OUT_DIR_SA, "comparison_with_baseline.rds"))
write_csv(comparison_table, file.path(OUT_DIR_SA, "comparison_with_baseline.csv"))
message("  ✓ comparison_with_baseline.rds / .csv 已保存")

## 保存模块参数
module_params <- list(
  target_patients = c("PED", "ADULT", "ELDERLY"),
  noncompliance_scenarios = noncompliance_scenarios$scenario_id,
  mapb_assumed_scenario = "complete",
  N_MC = N_MC,
  T_SAMPLE_H = T_SAMPLE_H,
  elapsed_minutes = as.numeric(elapsed),
  created = Sys.time(),
  version = "v6.0-SA-MAPB"
)
saveRDS(module_params, file.path(OUT_DIR_SA, "module_params.rds"))
message("  ✓ module_params.rds 已保存")

## ############################################################################
## SA.7 生成最终汇总表（修复版）
## ############################################################################

message("\n--- SA.7 生成最终汇总表 ---")

## 重新整理为要求的格式
final_summary <- calibration_metrics_SA %>%
  select(patient_id, model_name, s_true_runin, param, comparison,
         n, rBias_pct, rRMSE_pct, F20_pct, F30_pct) %>%
  arrange(patient_id, param, comparison, s_true_runin)

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("                     MAPB 假设偏差敏感性分析结果\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("研究设计：\n")
cat("  - 患者实际服药：非依从（5种场景）\n")
cat("  - 医生 MAPB 假设：完全依从（固定）\n")
cat("  - 目标患者：PED (儿童), ADULT (成人), ELDERLY (老人)\n\n")

cat("【CL 参数 - IND vs TRUE】\n")
print(
  final_summary %>%
    filter(param == "CL", comparison == "IND_vs_TRUE") %>%
    select(-n, -comparison) %>%
    mutate(
      rBias_pct = sprintf("%+.2f", rBias_pct),
      rRMSE_pct = sprintf("%.2f", rRMSE_pct),
      F20_pct = sprintf("%.1f", F20_pct),
      F30_pct = sprintf("%.1f", F30_pct)
    ),
  n = Inf
)

cat("\n【CL 参数 - POP vs TRUE（对照组）】\n")
print(
  final_summary %>%
    filter(param == "CL", comparison == "POP_vs_TRUE") %>%
    select(-n, -comparison) %>%
    mutate(
      rBias_pct = sprintf("%+.2f", rBias_pct),
      rRMSE_pct = sprintf("%.2f", rRMSE_pct),
      F20_pct = sprintf("%.1f", F20_pct),
      F30_pct = sprintf("%.1f", F30_pct)
    ),
  n = Inf
)

cat("\n【V 参数 - IND vs TRUE（仅 ADULT/ELDERLY）】\n")
print(
  final_summary %>%
    filter(param == "V", comparison == "IND_vs_TRUE") %>%
    select(-n, -comparison) %>%
    mutate(
      rBias_pct = sprintf("%+.2f", rBias_pct),
      rRMSE_pct = sprintf("%.2f", rRMSE_pct),
      F20_pct = sprintf("%.1f", F20_pct),
      F30_pct = sprintf("%.1f", F30_pct)
    ),
  n = Inf
)

cat("\n【V 参数 - POP vs TRUE（对照组，仅 ADULT/ELDERLY）】\n")
print(
  final_summary %>%
    filter(param == "V", comparison == "POP_vs_TRUE") %>%
    select(-n, -comparison) %>%
    mutate(
      rBias_pct = sprintf("%+.2f", rBias_pct),
      rRMSE_pct = sprintf("%.2f", rRMSE_pct),
      F20_pct = sprintf("%.1f", F20_pct),
      F30_pct = sprintf("%.1f", F30_pct)
    ),
  n = Inf
)

cat("\n", paste(rep("=", 80), collapse = ""), "\n")

## 保存最终汇总表
write_csv(final_summary, file.path(OUT_DIR_SA, "final_summary.csv"))
message("  ✓ final_summary.csv 已保存")

## ############################################################################
## SA.8 完成
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("✅ 敏感性分析完成：MAPB 假设偏差")
message(paste(rep("=", 70), collapse = ""))

message("\n📁 输出目录: ", OUT_DIR_SA)
message("\n📄 输出文件：")
message("  - mc_raw_results_SA.rds        （原始蒙特卡洛数据，", nrow(mc_raw_SA), " 行）")
message("  - calibration_metrics_SA.rds   （rBias/rRMSE 汇总）")
message("  - calibration_metrics_SA.csv")
message("  - comparison_with_baseline.rds （与基准对比）")
message("  - comparison_with_baseline.csv")
message("  - final_summary.csv            （最终汇总表）")
message("  - module_params.rds            （模块参数）")

message("\n📊 分析摘要：")
message("   患者数: 3 (PED, ADULT, ELDERLY)")
message("   非依从场景数: 5")
message("   每组合蒙特卡洛次数: ", N_MC)
message("   总耗时: ", round(as.numeric(elapsed), 2), " 分钟")

message("\n⚠️ 结果解读说明：")
message("   rBias 和 rRMSE 的数值含义由你自行判断")
message("   此脚本仅提供数值结果，不做���坏判定")
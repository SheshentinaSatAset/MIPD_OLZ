## ############################################################################
## 第九版代码-敏感性分析-模块1.R
## 
## 核心任务：
##   参数化的数据生成模块，接受 TRUE_MODEL_NAME 作为输入
##   用于敏感性分析：分别以 li2018, yin2016, sun2021 作为真实模型
## 
## 设计说明：
##   - 此脚本由主控脚本调用，不单独运行
##   - 需要预先设置全局变量：SA_TRUE_MODEL, OUT_DIR_SA
##   - 生成的数据结构与原模块1一致，确保模块2-4无需修改
## 
## 输入（从主控脚本传入）：
##   - SA_TRUE_MODEL: 真实模型名称 ("li2018", "yin2016", "sun2021")
##   - OUT_DIR_SA: 输出目录
## 
## 输出：
##   - {OUT_DIR_SA}/data_true.rds
##   - {OUT_DIR_SA}/scenario_library.rds
##   - {OUT_DIR_SA}/module1_params.rds
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("敏感性分析 模块1：数据生成（真实模型 = ", SA_TRUE_MODEL, "）")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## SA1.0 加载依赖与路径设置
## ############################################################################

message("\n--- SA1.0 加载依赖与路径设置 ---")

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

## 输入路径（第六版输出）
IN_DIR_S1 <- file.path(PROJ_ROOT, "outputs6/outputs6_stage1")

## 创建输出目录
if (!dir.exists(OUT_DIR_SA)) {
  dir.create(OUT_DIR_SA, recursive = TRUE)
  message("📁 已创建目录: ", OUT_DIR_SA)
}

## ============================================================================
## 核心参数设置（与原模块1一致）
## ============================================================================
MASTER_SEED <- 20260222
N_MC <- 5000
T_SAMPLE_H <- 336  # 采样时间 = 第14次给药后24h
DOSE_MG <- 10
II_H <- 24
N_DOSE <- 14

set.seed(MASTER_SEED)

message("\n📊 核心参数：")
message("   真实模型 = ", SA_TRUE_MODEL)
message("   全局随机种子 = ", MASTER_SEED)
message("   虚拟患者数 N_MC = ", N_MC)
message("   采样时间 t = ", T_SAMPLE_H, " h")

## 读取前置输出
omega_sigma_bank <- readRDS(file.path(IN_DIR_S1, "omega_sigma_bank.rds"))

## 加载真实模型
source(file.path(MODEL_DIR, "li2018_olz_model_drug2.R"))
source(file.path(MODEL_DIR, "yin2016_olz_model.R"))
source(file.path(MODEL_DIR, "sun2021_olz_model.R"))

message("\n✅ PopPK 模型已加载")

## ############################################################################
## SA1.1 定义真实模型的参数
## ############################################################################

message("\n--- SA1.1 定义真实模型的参数 ---")

## 获取模型银行中的参数
model_bank <- omega_sigma_bank[[SA_TRUE_MODEL]]

## 提取 Ω 矩阵（用于抽样 η_true）
omega_true <- model_bank$omega_calibrate

## 提取残差参数
sigma_type <- model_bank$sigma_type
sigma_prop_sd <- model_bank$sigma_prop_sd
sigma_add_sd <- model_bank$sigma_add_sd

## ============================================================================
## 根据真实模型定义协变量和计算 θ_pop
## ============================================================================

if (SA_TRUE_MODEL == "li2018") {
  
  ## Li2018 协变量：WT = 65 kg
  PATIENT_COVARIATES <- list(WT = 65)
  
  ## 计算 θ_pop
  TVCL <- 25.4
  TVVC <- 2390
  WT_REF <- 60.59
  B_WT_VC_EXP <- 0.579
  
  CL_POP <- TVCL
  V_POP <- TVVC * (65 / WT_REF)^B_WT_VC_EXP
  
  ## 协变量 data.frame
  iCov_true <- data.frame(id = 1, WT = 65)
  
  ## rxode2 模型对象
  true_model_obj <- li2018_olz_model_drug2
  
  ## η 参数名称
  eta_names <- c("eta_CL", "eta_Vc")
  
  ## 其他 η 参数（设为 0）
  other_etas <- c(eta_KA = 0, eta_Vp = 0, eta_Q = 0)
  
} else if (SA_TRUE_MODEL == "yin2016") {
  
  ## Yin2016 协变量：SEX = 1 (男), SMK = 0 (不吸烟)
  PATIENT_COVARIATES <- list(SEX = 1, SMK = 0)
  
  ## 计算 θ_pop
  TVCL <- 16.6
  TVV2 <- 599
  B_SEX <- 0.385
  B_SMK <- 0.426
  
  CL_POP <- TVCL * exp(B_SEX * 1 + B_SMK * 0)
  V_POP <- TVV2
  
  ## 协变量 data.frame
  iCov_true <- data.frame(id = 1, SEX = 1, SMK = 0)
  
  ## rxode2 模型对象
  true_model_obj <- yin2016_olz_model
  
  ## η 参数名称
  eta_names <- c("eta_CL", "eta_V2")
  
  ## 其他 η 参数（设为 0）
  other_etas <- c(eta_KA = 0, eta_Q = 0, eta_V3 = 0, eta_TLAG = 0)
  
} else if (SA_TRUE_MODEL == "sun2021") {
  
  ## Sun2021 协变量
  ## 注意：SEX=0 表示男性！
  PATIENT_COVARIATES <- list(
    AGE = 22, SEX = 0, WT = 65, SMOKE = 0,
    RACE_BLACK = 0, RIF = 0, HEP_MOD = 0, REN_SEV = 0, FED = 0
  )
  
  ## 计算 θ_pop
  TVCL <- 15.5
  TVVC <- 656
  WT_REF <- 70
  AGE_REF <- 36
  B_WT_CL <- 0.75
  B_WT_VC <- 1.0
  B_AGE_VC <- 0.356
  
  CL_POP <- TVCL * exp(B_WT_CL * log(65 / WT_REF))
  V_POP <- TVVC * exp(
    B_WT_VC * log(65 / WT_REF) +
    B_AGE_VC * log(22 / AGE_REF)
  )
  
  ## 协变量 data.frame
  iCov_true <- data.frame(
    id = 1, AGE = 22, SEX = 0, WT = 65, SMOKE = 0,
    RACE_BLACK = 0, RIF = 0, HEP_MOD = 0, REN_SEV = 0, FED = 0
  )
  
  ## rxode2 模型对象
  true_model_obj <- sun2021_olz_model
  
  ## η 参数名称
  eta_names <- c("eta_CL", "eta_Vc")
  
  ## 其他 η 参数（设为 0）
  other_etas <- c(eta_KA = 0, eta_KA_IOV = 0, eta_Vp = 0, eta_Q = 0, eta_ALAG = 0)
  
} else {
  stop("未知的真实模型: ", SA_TRUE_MODEL)
}

cat("\n=== 真实模型参数 ===\n")
cat("模型: ", SA_TRUE_MODEL, "\n")
cat("θ_pop: CL =", round(CL_POP, 2), "L/h, V =", round(V_POP, 2), "L\n")
cat("Ω 矩阵:\n")
print(round(omega_true, 4))
cat("残差模型:", sigma_type, ", σ_prop =", sigma_prop_sd)
if (!is.na(sigma_add_sd)) cat(", σ_add =", sigma_add_sd)
cat("\n")

## ############################################################################
## SA1.2 定义依从性场景库
## ############################################################################

message("\n--- SA1.2 定义依从性场景库 ---")

scenario_library <- tibble(
  scenario_id = c("complete", "miss_last", "miss_prev", "miss_two", 
                  "double_last", "double_prev"),
  scenario_name = c("完全依从", "末次漏服", "倒数第二次漏服", 
                    "连续漏服两次", "末次双倍", "倒数第二次双倍"),
  dose_13_mult = c(1, 1, 0, 0, 1, 2),
  dose_14_mult = c(1, 0, 1, 0, 2, 1)
)

N_SCENARIOS <- nrow(scenario_library)

cat("\n=== 依从性场景库 ===\n")
print(scenario_library)

## ############################################################################
## SA1.3 定义辅助函��
## ############################################################################

message("\n--- SA1.3 定义辅助函数 ---")

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
add_residual_error <- function(C_pred, sigma_type, sigma_prop_sd, sigma_add_sd) {
  
  if (sigma_type == "proportional") {
    eps <- rnorm(1, mean = 0, sd = sigma_prop_sd)
    C_obs <- C_pred * (1 + eps)
    
  } else if (sigma_type == "combined") {
    sigma_total <- sqrt((sigma_prop_sd * C_pred)^2 + sigma_add_sd^2)
    eps <- rnorm(1, mean = 0, sd = sigma_total)
    C_obs <- C_pred + eps
    
  } else {
    stop("未知残差类型: ", sigma_type)
  }
  
  max(C_obs, 1e-6)
}

message("✅ 辅助函数已定义")

## ############################################################################
## SA1.4 生成虚拟患者数据
## ############################################################################

message("\n--- SA1.4 生成虚拟患者数据 ---")
message("正在生成 ", N_MC, " 个虚拟患者...")

## 预构建先导期事件表（complete 场景）
ev_runin <- make_ev_for_scenario(
  dose_mg = DOSE_MG, ii_h = II_H, n_dose = N_DOSE,
  t_sample_h = T_SAMPLE_H,
  dose_13_mult = 1, dose_14_mult = 1  # complete
)

## 初始化结果存储
data_true <- vector("list", N_MC)

## 蒙特卡洛循环
for (i in seq_len(N_MC)) {
  
  if (i %% 100 == 0 || i == 1) {
    cat(sprintf("\r  进度: %d/%d", i, N_MC))
  }
  
  ## ========== Step 1: 抽样 η_true ==========
  eta_vec <- as.vector(rmvnorm(1, mean = c(0, 0), sigma = omega_true))
  eta_CL_true <- eta_vec[1]
  eta_V_true <- eta_vec[2]
  
  ## ========== Step 2: 计算 θ_true ==========
  CL_true <- CL_POP * exp(eta_CL_true)
  V_true <- V_POP * exp(eta_V_true)
  
  ## ========== Step 3: 构建参数向量 ==========
  ## 根据真实模型构建参数
  params_true <- c(other_etas)  # 先填充其他 η 为 0
  params_true[eta_names[1]] <- eta_CL_true
  params_true[eta_names[2]] <- eta_V_true
  
  ## ========== Step 4: 模拟先导期浓度（complete 场景）==========
  sim_runin <- rxSolve(
    object = true_model_obj,
    params = params_true,
    events = ev_runin,
    iCov = iCov_true,
    returnType = "data.frame"
  )
  
  C_pred_runin <- extract_conc_at_time(sim_runin, T_SAMPLE_H)
  C_obs_runin <- add_residual_error(C_pred_runin, sigma_type, sigma_prop_sd, sigma_add_sd)
  
  ## ========== Step 5: 随机选择测试场景 ==========
  s_true_idx <- sample(1:N_SCENARIOS, 1)
  s_true <- scenario_library$scenario_id[s_true_idx]
  
  ## ========== Step 6: 模拟测试期浓度 ==========
  ev_test <- make_ev_for_scenario(
    dose_mg = DOSE_MG, ii_h = II_H, n_dose = N_DOSE,
    t_sample_h = T_SAMPLE_H,
    dose_13_mult = scenario_library$dose_13_mult[s_true_idx],
    dose_14_mult = scenario_library$dose_14_mult[s_true_idx]
  )
  
  sim_test <- rxSolve(
    object = true_model_obj,
    params = params_true,
    events = ev_test,
    iCov = iCov_true,
    returnType = "data.frame"
  )
  
  C_pred_test <- extract_conc_at_time(sim_test, T_SAMPLE_H)
  C_obs_test <- add_residual_error(C_pred_test, sigma_type, sigma_prop_sd, sigma_add_sd)
  
  ## ========== 存储结果 ==========
  data_true[[i]] <- tibble(
    mc_iter = i,
    
    ## 真实随机效应
    eta_CL_true = eta_CL_true,
    eta_V_true = eta_V_true,
    
    ## 真实参数
    CL_true = CL_true,
    V_true = V_true,
    CL_pop = CL_POP,
    V_pop = V_POP,
    
    ## 先导期浓度（用于 MAPB）
    C_pred_runin = C_pred_runin,
    C_obs_runin = C_obs_runin,
    
    ## 测试期场景和浓度
    s_true = s_true,
    C_pred_test = C_pred_test,
    C_obs_test = C_obs_test
  )
}

cat("\n")

## 合并结果
data_true <- bind_rows(data_true)

message("\n✅ 数据生成完成")
message("   总行数: ", nrow(data_true))

## ############################################################################
## SA1.5 数据质量检查
## ############################################################################

message("\n--- SA1.5 数据质量检查 ---")

## 场景分布
cat("\n=== 场景分布 ===\n")
print(table(data_true$s_true))

## η_true 统计
cat("\n=== η_true 统计 ===\n")
cat("η_CL_true: mean =", round(mean(data_true$eta_CL_true), 4),
    ", sd =", round(sd(data_true$eta_CL_true), 4), "\n")
cat("η_V_true:  mean =", round(mean(data_true$eta_V_true), 4),
    ", sd =", round(sd(data_true$eta_V_true), 4), "\n")

## 预期 SD（从 Ω 矩阵）
expected_sd_CL <- sqrt(omega_true[1, 1])
expected_sd_V <- sqrt(omega_true[2, 2])
cat("\n预期 SD（从 Ω）: η_CL ~", round(expected_sd_CL, 4),
    ", η_V ~", round(expected_sd_V, 4), "\n")

## θ_true 统计
cat("\n=== θ_true 统计 ===\n")
cat("CL_true: mean =", round(mean(data_true$CL_true), 2), "L/h",
    ", range = [", round(min(data_true$CL_true), 2), ",",
    round(max(data_true$CL_true), 2), "]\n")
cat("V_true:  mean =", round(mean(data_true$V_true), 2), "L",
    ", range = [", round(min(data_true$V_true), 2), ",",
    round(max(data_true$V_true), 2), "]\n")

## 浓度统计
cat("\n=== 浓度统计 ===\n")
cat("C_obs_runin: mean =", round(mean(data_true$C_obs_runin), 2), "ng/mL",
    ", range = [", round(min(data_true$C_obs_runin), 2), ",",
    round(max(data_true$C_obs_runin), 2), "]\n")
cat("C_obs_test:  mean =", round(mean(data_true$C_obs_test), 2), "ng/mL",
    ", range = [", round(min(data_true$C_obs_test), 2), ",",
    round(max(data_true$C_obs_test), 2), "]\n")

## 检查 NA
n_na <- sum(is.na(data_true$C_obs_runin) | is.na(data_true$C_obs_test))
if (n_na > 0) {
  warning("发现 ", n_na, " 个 NA 值！")
} else {
  message("✅ 无 NA 值")
}

## ############################################################################
## SA1.6 保存输出文件
## ############################################################################

message("\n--- SA1.6 保存输出文件 ---")

## 保存核心数据
saveRDS(data_true, file.path(OUT_DIR_SA, "data_true.rds"))
write_csv(data_true, file.path(OUT_DIR_SA, "data_true.csv"))
message("  ✓ data_true.rds / .csv 已保存（", nrow(data_true), " 行）")

## 保存场景库
saveRDS(scenario_library, file.path(OUT_DIR_SA, "scenario_library.rds"))
message("  ✓ scenario_library.rds 已保存")

## 保存模块参数（供后续模块使用）
module1_params <- list(
  ## 随机种子
  MASTER_SEED = MASTER_SEED,
  N_MC = N_MC,
  
  ## 真实模型信息
  TRUE_MODEL_NAME = SA_TRUE_MODEL,
  omega_true = omega_true,
  sigma_type = sigma_type,
  sigma_prop_sd = sigma_prop_sd,
  sigma_add_sd = sigma_add_sd,
  
  ## 患者协变量
  PATIENT_COVARIATES = PATIENT_COVARIATES,
  
  ## 给药信息
  DOSE_MG = DOSE_MG,
  II_H = II_H,
  N_DOSE = N_DOSE,
  T_SAMPLE_H = T_SAMPLE_H,
  
  ## θ_pop
  CL_POP = CL_POP,
  V_POP = V_POP,
  
  ## 元数据
  created = Sys.time(),
  version = "v9.0-SA",
  module = "Module 1 - Data Generation (Sensitivity Analysis)"
)

saveRDS(module1_params, file.path(OUT_DIR_SA, "module1_params.rds"))
message("  ✓ module1_params.rds 已保存")

## ############################################################################
## SA1.7 模块 1 完成
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("✅ 敏感性分析模块 1 完成（真实模型 = ", SA_TRUE_MODEL, "）")
message(paste(rep("=", 70), collapse = ""))

message("\n📁 输出目录: ", OUT_DIR_SA)
message("\n📄 输出文件：")
message("  - data_true.rds / .csv    （虚拟患者数据，", nrow(data_true), " 行）")
message("  - scenario_library.rds    （依从性场景库）")
message("  - module1_params.rds      （模块参数）")

message("\n📊 数据生成摘要：")
message("   真实模型: ", SA_TRUE_MODEL)
message("   虚拟患者数: ", N_MC)
message("   θ_pop: CL = ", round(CL_POP, 2), " L/h, V = ", round(V_POP, 2), " L")
message("   η_CL_true SD: ", round(sd(data_true$eta_CL_true), 4), 
        "（预期 ", round(expected_sd_CL, 4), "）")
message("   η_V_true SD: ", round(sd(data_true$eta_V_true), 4),
        "（预期 ", round(expected_sd_V, 4), "）")
## ############################################################################
## 第九版代码-模块1-数据生成.R
## 
## 核心任务：
##   使用 zang2021 模型作为"临床真实模型"，生成唯一的虚拟患者数据
##   这份数据将被交给 4 个模型分别分析，验证跨模型一致性
## 
## 关键设计变化（第九版 vs 之前）：
##   - 之前：4 个模型各自生成 η_true 和 C_obs（4 组不同的虚拟患者）
##   - 现在：仅 zang2021 生成 η_true 和 C_obs（1 组虚拟患者）
##         然后交给 4 个模型分析同一份数据
## 
## 输入：
##   - outputs6/outputs6_stage1/omega_sigma_bank.rds
##   - outputs6/outputs6_stage2/design_manifest_stage2.rds
## 
## 输出：
##   - outputs9/data_true.rds  ← 核心输出！
##   - outputs9/scenario_library.rds
##   - outputs9/module1_params.rds
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("第九版 模块1：数据生成（唯一真实模型：zang2021）")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## 1.0 加载依赖与路径设置
## ############################################################################

message("\n--- 1.0 加载依赖与路径设置 ---")

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
INPUT_ROOT <- file.path(PROJ_ROOT, "outputs6")
IN_DIR_S1 <- file.path(INPUT_ROOT, "outputs6_stage1")
IN_DIR_S2 <- file.path(INPUT_ROOT, "outputs6_stage2")

## 输出路径（第九版）
OUT_DIR_9 <- file.path(PROJ_ROOT, "outputs9")

## 创建输出目录
if (!dir.exists(OUT_DIR_9)) {
  dir.create(OUT_DIR_9, recursive = TRUE)
  message("📁 已创建目录: ", OUT_DIR_9)
}

## ============================================================================
## 核心参数设置
## ============================================================================
MASTER_SEED <- 20260222
N_MC <- 5000
T_SAMPLE_H <- 336  # 采样时间 = 第14次给药后24h

set.seed(MASTER_SEED)

message("\n📊 核心参数：")
message("   全局随机种子 = ", MASTER_SEED)
message("   虚拟患者数 N_MC = ", N_MC)
message("   采样时间 t = ", T_SAMPLE_H, " h")

## 读取前置输出
omega_sigma_bank <- readRDS(file.path(IN_DIR_S1, "omega_sigma_bank.rds"))
design_manifest <- readRDS(file.path(IN_DIR_S2, "design_manifest_stage2.rds"))

## 加载 zang2021 模型
source(file.path(MODEL_DIR, "zang2021_olz_model.R"))

message("\n✅ 前置数据已加载")
message("✅ zang2021 模型已加载")

## ############################################################################
## 1.1 定义"临床真实模型"参数
## ############################################################################

message("\n--- 1.1 定义临床真实模型参数 ---")

## 选择 zang2021 作为"真实模型"
TRUE_MODEL_NAME <- "zang2021"
model_bank <- omega_sigma_bank[[TRUE_MODEL_NAME]]

## 提取 Ω 矩阵（用于抽样 η_true）
omega_true <- model_bank$omega_calibrate  # 2×2 矩阵 (eta_CL, eta_V)

## 提取残差参数
sigma_type <- model_bank$sigma_type  # "proportional"
sigma_prop_sd <- model_bank$sigma_prop_sd  # 0.20

## 定义固定协变量（代表"同一个临床患者"）
## 22岁，男性，体重65kg，不吸烟，无感染，合并丙戊酸钠，其他正常
PATIENT_COVARIATES <- list(
  SEX = 1,   # 男性
  SMK = 0,   # 不吸烟
  INF = 0,   # 无感染
  VAL = 1,   # 合并丙戊酸钠
  PER = 0,   # 无奋乃静
  SER = 0,   # 无舍曲林
  FLU = 0,   # 无氟伏沙明
  DGH = 0    # 无当归龙荟片
)

## 给药信息
DOSE_MG <- 10
II_H <- 24
N_DOSE <- 14

## 计算 θ_pop（使用 zang2021 的协变量效应）
TVCL <- 12.88
TVV <- 754.41
B_SEX <- 0.21
B_SMK <- 0.21
B_INF <- -0.29
B_VAL <- 0.21
B_PER <- -0.25
B_SER <- 0.15
B_FLU <- -0.36
B_DGH <- 0.70

CL_POP <- TVCL * exp(
  B_SEX * PATIENT_COVARIATES$SEX +
  B_SMK * PATIENT_COVARIATES$SMK +
  B_INF * PATIENT_COVARIATES$INF +
  B_VAL * PATIENT_COVARIATES$VAL +
  B_PER * PATIENT_COVARIATES$PER +
  B_SER * PATIENT_COVARIATES$SER +
  B_FLU * PATIENT_COVARIATES$FLU +
  B_DGH * PATIENT_COVARIATES$DGH
)

V_POP <- TVV  # V 没有协变量效应

cat("\n=== 临床真实模型参数 ===\n")
cat("模型: ", TRUE_MODEL_NAME, "\n")
cat("θ_pop: CL =", round(CL_POP, 2), "L/h, V =", round(V_POP, 2), "L\n")
cat("Ω 矩阵:\n")
print(round(omega_true, 4))
cat("残差模型:", sigma_type, ", σ_prop =", sigma_prop_sd, "\n")

## ############################################################################
## 1.2 定义依从性场景库
## ############################################################################

message("\n--- 1.2 定义依从性场景库 ---")

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
## 1.3 定义辅助函数
## ############################################################################

message("\n--- 1.3 定义辅助函数 ---")

## ----------------------------------------------------------------------------
## 1.3.1 构建给药事件表（指定场景）
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
## 1.3.2 提取指定时间的浓度
## ----------------------------------------------------------------------------
extract_conc_at_time <- function(sim_df, t_target) {
  idx <- which(abs(sim_df$time - t_target) < 1e-8)
  if (length(idx) == 0) return(NA_real_)
  sim_df$C_ngmL[idx[1]]
}

## ----------------------------------------------------------------------------
## 1.3.3 添加残差（比例误差）
## ----------------------------------------------------------------------------
add_residual_error <- function(C_pred, sigma_prop_sd) {
  eps <- rnorm(1, mean = 0, sd = sigma_prop_sd)
  C_obs <- C_pred * (1 + eps)
  max(C_obs, 1e-6)  # 确保浓度为正
}

message("✅ 辅助函数已定义")

## ############################################################################
## 1.4 生成虚拟患者数据
## ############################################################################

message("\n--- 1.4 生成虚拟患者数据 ---")
message("正在生成 ", N_MC, " 个虚拟患者...")

## 准备协变量 data.frame
iCov <- data.frame(
  id = 1,
  SEX = PATIENT_COVARIATES$SEX,
  SMK = PATIENT_COVARIATES$SMK,
  INF = PATIENT_COVARIATES$INF,
  VAL = PATIENT_COVARIATES$VAL,
  PER = PATIENT_COVARIATES$PER,
  SER = PATIENT_COVARIATES$SER,
  FLU = PATIENT_COVARIATES$FLU,
  DGH = PATIENT_COVARIATES$DGH
)

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
  
  ## ========== Step 3: 模拟先导期浓度（complete 场景）==========
  sim_runin <- rxSolve(
    object = zang2021_olz_model,
    params = c(eta_CL = eta_CL_true, eta_V = eta_V_true),
    events = ev_runin,
    iCov = iCov,
    returnType = "data.frame"
  )
  
  C_pred_runin <- extract_conc_at_time(sim_runin, T_SAMPLE_H)
  C_obs_runin <- add_residual_error(C_pred_runin, sigma_prop_sd)
  
  ## ========== Step 4: 随机选择测试场景 ==========
  s_true_idx <- sample(1:N_SCENARIOS, 1)
  s_true <- scenario_library$scenario_id[s_true_idx]
  
  ## ========== Step 5: 模拟测试期浓度 ==========
  ev_test <- make_ev_for_scenario(
    dose_mg = DOSE_MG, ii_h = II_H, n_dose = N_DOSE,
    t_sample_h = T_SAMPLE_H,
    dose_13_mult = scenario_library$dose_13_mult[s_true_idx],
    dose_14_mult = scenario_library$dose_14_mult[s_true_idx]
  )
  
  sim_test <- rxSolve(
    object = zang2021_olz_model,
    params = c(eta_CL = eta_CL_true, eta_V = eta_V_true),
    events = ev_test,
    iCov = iCov,
    returnType = "data.frame"
  )
  
  C_pred_test <- extract_conc_at_time(sim_test, T_SAMPLE_H)
  C_obs_test <- add_residual_error(C_pred_test, sigma_prop_sd)
  
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
## 1.5 数据质量检查
## ############################################################################

message("\n--- 1.5 数据质量检查 ---")

## 场景分布
cat("\n=== 场景分布 ===\n")
print(table(data_true$s_true))

## η_true 统计
cat("\n=== η_true 统计 ===\n")
cat("η_CL_true: mean =", round(mean(data_true$eta_CL_true), 4),
    ", sd =", round(sd(data_true$eta_CL_true), 4), "\n")
cat("η_V_true:  mean =", round(mean(data_true$eta_V_true), 4),
    ", sd =", round(sd(data_true$eta_V_true), 4), "\n")

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

## 检查是否有异常值
n_na <- sum(is.na(data_true$C_obs_runin) | is.na(data_true$C_obs_test))
if (n_na > 0) {
  warning("发现 ", n_na, " 个 NA 值！")
} else {
  message("✅ 无 NA 值")
}

## ############################################################################
## 1.6 保存输出文件
## ############################################################################

message("\n--- 1.6 保存输出文件 ---")

## 保存核心数据
saveRDS(data_true, file.path(OUT_DIR_9, "data_true.rds"))
write_csv(data_true, file.path(OUT_DIR_9, "data_true.csv"))
message("  ✓ data_true.rds / .csv 已保存（", nrow(data_true), " 行）")

## 保存场景库
saveRDS(scenario_library, file.path(OUT_DIR_9, "scenario_library.rds"))
message("  ✓ scenario_library.rds 已保存")

## 保存模块参数（供后续模块使用）
module1_params <- list(
  ## 随机种子
  MASTER_SEED = MASTER_SEED,
  N_MC = N_MC,
  
  ## 真实模型信息
  TRUE_MODEL_NAME = TRUE_MODEL_NAME,
  omega_true = omega_true,
  sigma_type = sigma_type,
  sigma_prop_sd = sigma_prop_sd,
  
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
  version = "v9.0",
  module = "Module 1 - Data Generation"
)

saveRDS(module1_params, file.path(OUT_DIR_9, "module1_params.rds"))
message("  ✓ module1_params.rds 已保存")

## ############################################################################
## 1.7 模块 1 完成
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("✅ 模块 1 全部完成（数据生成）")
message(paste(rep("=", 70), collapse = ""))

message("\n📁 输出目录: ", OUT_DIR_9)
message("\n📄 输出文件：")
message("  - data_true.rds / .csv    （虚拟患者数据，", nrow(data_true), " 行）")
message("  - scenario_library.rds    （依从性场景库）")
message("  - module1_params.rds      （模块参数）")

message("\n📊 数据生成摘要：")
message("   真实模型: ", TRUE_MODEL_NAME)
message("   虚拟患者数: ", N_MC)
message("   场景数: ", N_SCENARIOS)
message("   θ_pop: CL = ", round(CL_POP, 2), " L/h, V = ", round(V_POP, 2), " L")

message("\n🔗 后续模块使用方法：")
message('  data_true <- readRDS("outputs9/data_true.rds")')
message('  params <- readRDS("outputs9/module1_params.rds")')

message("\n⏭️  下一步：运行 [第九版代码-模块2-跨模型MAPB.R]")
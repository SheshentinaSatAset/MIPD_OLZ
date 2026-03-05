## ############################################################################
## 第六版代码-主线0&1.R
## 
## Stage 0：前置准备工作
## Stage 1：PopPK 模型的加载与验证
## 
## 核心输出（供后续脚本调用）：
##   - outputs6/outputs6_stage1/model_info.csv
##   - outputs6/outputs6_stage1/omega_sigma_bank.rds  ← 关键！参数银行
##   - outputs6/outputs6_stage1/model_file_paths.rds
## 
## 设计原则：
##   - 后续脚本只需读取 outputs6_stage1/ 中的文件，无需重新 source 模型
##   - 所有模型的 omega、sigma 集中管理，便于维护
##   - omega_sigma_bank 是 Single Source of Truth
## 
## 修复内容（相对于原版）：
##   1.统一 sigma 的存储方式（全部存 SD，不存方差）
##   2.修复 yin2016 的 omega 计算（THETA_F 在 lotri 中无效）
##   3.为每个模型添加 sigma_prop_sd / sigma_add_sd 字段
##   4.修复 model_info 表的 sigma_prop 引用
## 
## ############################################################################

## ############################################################################
## Stage 0：前置准备工作
## ############################################################################

message("/n", paste(rep("=", 70), collapse = ""))
message("Stage 0：前置准备工作（第六版）")
message(paste(rep("=", 70), collapse = ""))

## ----------------------------------------------------------------------------
## 0.1 清空环境 + 设置全局选项
## ----------------------------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = FALSE)
options(digits = 6)

## ----------------------------------------------------------------------------
## 0.2 加载依赖包
## ----------------------------------------------------------------------------
message("/n--- 0.2 加载依赖包 ---")

suppressPackageStartupMessages({
  library(rxode2)
  library(lotri)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(readr)
})

## 检查并安装 mvtnorm（后续 Stage 需要）
if (!requireNamespace("mvtnorm", quietly = TRUE)) {
  message("📦 正在安装 mvtnorm 包...")
  install.packages("mvtnorm")
}
library(mvtnorm)

message("✅ 依赖包加载完成")

## ----------------------------------------------------------------------------
## 0.3 全局种子与常量
## ----------------------------------------------------------------------------
MASTER_SEED <- 20260129
set.seed(MASTER_SEED)

## 蒙特卡洛模拟次数
N_MC <- 5000

message("🎲 全局随机种子: ", MASTER_SEED)
message("🔢 蒙特卡洛次数: ", N_MC)

## ----------------------------------------------------------------------------
## 0.4 路径设置（第六版）
## ----------------------------------------------------------------------------
message("/n--- 0.4 路径设置 ---")

## 项目根目录
PROJ_ROOT <- "D:/Users/YujiaZhang/Desktop/26救赎之道/25.9.29 毕业设计/5正式实施！/多成人模型/体现我卓越的项目管理能力"

## 代码目录
CODE_DIR <- file.path(PROJ_ROOT, "代码这边请/正式代码/第六版代码")

## PopPK 模型目录
MODEL_DIR <- file.path(PROJ_ROOT, "代码这边请/poppk模型")

## 输出总目录（第六版）
OUTPUT_ROOT <- file.path(PROJ_ROOT, "outputs6")

## 各 Stage 输出目录
OUT_DIR_S1 <- file.path(OUTPUT_ROOT, "outputs6_stage1")
OUT_DIR_S2 <- file.path(OUTPUT_ROOT, "outputs6_stage2")
OUT_DIR_S3 <- file.path(OUTPUT_ROOT, "outputs6_stage3")
OUT_DIR_S4 <- file.path(OUTPUT_ROOT, "outputs6_stage4")
OUT_DIR_S5 <- file.path(OUTPUT_ROOT, "outputs6_stage5")
OUT_DIR_S6 <- file.path(OUTPUT_ROOT, "outputs6_stage6")

## 创建输出目录
all_out_dirs <- c(OUT_DIR_S1, OUT_DIR_S2, OUT_DIR_S3, 
                  OUT_DIR_S4, OUT_DIR_S5, OUT_DIR_S6)

for (dir_path in all_out_dirs) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("📁 已创建目录: ", dir_path)
  }
}

message("/n✅ Stage 0 完成：环境准备就绪")
message("   - 项目根目录: ", PROJ_ROOT)
message("   - 模型目录:   ", MODEL_DIR)
message("   - 输出目录:   ", OUTPUT_ROOT)

## ############################################################################
## Stage 1：PopPK 模型的加载与验证
## ############################################################################

message("/n", paste(rep("=", 70), collapse = ""))
message("Stage 1：PopPK 模型的加载与验证（第六版）")
message(paste(rep("=", 70), collapse = ""))

set.seed(MASTER_SEED + 1)

## ----------------------------------------------------------------------------
## 1.1 定义模型文件路径
## ----------------------------------------------------------------------------
message("/n--- 1.1 定义模型文件路径 ---")

model_files <- list(
  ## 主线模型
  zang2021    = file.path(MODEL_DIR, "zang2021_olz_model.R"),
  maharaj2021 = file.path(MODEL_DIR, "maharaj2021_olz_model.R"),
  ## 支线模型
  li2018      = file.path(MODEL_DIR, "li2018_olz_model_drug2.R"),
  yin2016     = file.path(MODEL_DIR, "yin2016_olz_model.R"),
  sun2021     = file.path(MODEL_DIR, "sun2021_olz_model.R")
)

## 检查文件是否存在
for (model_name in names(model_files)) {
  if (!file.exists(model_files[[model_name]])) {
    stop("❌ 模型文件不存在: ", model_files[[model_name]])
  }
}

message("✅ 所有模型文件已确认存在")

## ----------------------------------------------------------------------------
## 1.2 加载所有模型
## ----------------------------------------------------------------------------
message("/n--- 1.2 加载所有模型 ---")

## 主线模型
source(model_files$zang2021)
message("  ✓ zang2021_olz_model 已加载")

source(model_files$maharaj2021)
message("  ✓ maharaj2021_olz_model 已加载")

## 支线模型
source(model_files$li2018)
message("  ✓ li2018_olz_model_drug2 已加载")

source(model_files$yin2016)
message("  ✓ yin2016_olz_model 已加载")

source(model_files$sun2021)
message("  ✓ sun2021_olz_model 已加载")

message("✅ 5 个模型全部加载完成")

## ----------------------------------------------------------------------------
## 1.3 验证模型对象是否正确加载
## ----------------------------------------------------------------------------
message("/n--- 1.3 验证模型对象 ---")

## 检查主线模型对象
stopifnot(exists("zang2021_olz_model"))
stopifnot(exists("zang2021_omega"))
stopifnot(exists("zang2021_sigma"))

stopifnot(exists("maharaj2021_olz_model"))
stopifnot(exists("maharaj2021_omega"))
stopifnot(exists("maharaj2021_sigma"))

## 检查支线模型对象
stopifnot(exists("li2018_olz_model_drug2"))
stopifnot(exists("li2018_omega_drug2"))
stopifnot(exists("li2018_sigma_drug2"))

stopifnot(exists("yin2016_olz_model"))
stopifnot(exists("yin2016_omega"))
stopifnot(exists("yin2016_sigma"))

stopifnot(exists("sun2021_olz_model"))
stopifnot(exists("sun2021_omega"))
stopifnot(exists("sun2021_sigma"))

message("  ✓ 所有模型对象验证通过")

## ----------------------------------------------------------------------------
## 1.4 构建参数银行（omega_sigma_bank）— 核心输出
## ----------------------------------------------------------------------------
message("/n--- 1.4 构建参数银行 ---")

## ========== 预计算各模型的 omega/sigma ==========

## --- Maharaj2021 ---
## omega: eta_CL ~ 0.309（文献直接给出方差）
## sigma: prop_err ~ 0.0772（方差），add_err ~ 1.51（方差）
omega_CL_maharaj <- as.numeric(maharaj2021_omega["eta_CL", "eta_CL"])
sigma_prop_sd_maharaj <- sqrt(as.numeric(maharaj2021_sigma["prop_err", "prop_err"]))
sigma_add_sd_maharaj <- sqrt(as.numeric(maharaj2021_sigma["add_err", "add_err"]))

## --- Zang2021 ---
## omega: eta_CL ~ log((9.76/100)^2 + 1), eta_V ~ log((15.32/100)^2 + 1)
## sigma: prop_err ~ 0.20（这是 SD，不是方差！lotri 直接填 SD）
omega_CL_zang <- as.numeric(zang2021_omega["eta_CL", "eta_CL"])
omega_V_zang <- as.numeric(zang2021_omega["eta_V", "eta_V"])
omega_cov_CL_V_zang <- as.numeric(zang2021_omega["eta_CL", "eta_V"])
## 注意：zang2021_sigma 中 prop_err ~ 0.20 是 SD（lotri 特性）
sigma_prop_sd_zang <- as.numeric(zang2021_sigma["prop_err", "prop_err"])

## --- Li2018 ---
## omega: eta_CL + eta_Vc 有协方差
## sigma: prop_err ~ 0.216（SD），add_err ~ 0.303（SD）
omega_CL_li2018 <- as.numeric(li2018_omega_drug2["eta_CL", "eta_CL"])
omega_Vc_li2018 <- as.numeric(li2018_omega_drug2["eta_Vc", "eta_Vc"])
omega_cov_CL_Vc_li2018 <- as.numeric(li2018_omega_drug2["eta_CL", "eta_Vc"])
sigma_prop_sd_li2018 <- as.numeric(li2018_sigma_drug2["prop_err", "prop_err"])
sigma_add_sd_li2018 <- as.numeric(li2018_sigma_drug2["add_err", "add_err"])

## --- Yin2016 ---
## 注意：模型文件中 THETA_F=1.10 在 lotri 中无法识别，需要手动计算
THETA_F <- 1.10
omega_CL_yin2016 <- log(((36.7 * THETA_F) / 100)^2 + 1)
omega_V2_yin2016 <- log(((42.7 * THETA_F) / 100)^2 + 1)
## 相关系数 = 0.982
omega_cov_CL_V2_yin2016 <- 0.982 * sqrt(omega_CL_yin2016 * omega_V2_yin2016)
## sigma: prop_err ~ 0.311（SD），add_err ~ 1.03（SD）
sigma_prop_sd_yin2016 <- as.numeric(yin2016_sigma["prop_err", "prop_err"])
sigma_add_sd_yin2016 <- as.numeric(yin2016_sigma["add_err", "add_err"])

## --- Sun2021 ---
## omega: 对角矩阵（eta_CL, eta_Vc 独立）
## sigma: prop_err ~ sqrt(0.0462)（已经是 SD）
omega_CL_sun2021 <- as.numeric(sun2021_omega["eta_CL", "eta_CL"])
omega_Vc_sun2021 <- as.numeric(sun2021_omega["eta_Vc", "eta_Vc"])
sigma_prop_sd_sun2021 <- as.numeric(sun2021_sigma["prop_err", "prop_err"])

## ========== 构建参数银行 ==========

omega_sigma_bank <- list(
  
  ## ========== Maharaj2021（主线-儿童）==========
  maharaj2021 = list(
    model_type = "主线-儿童",
    compartment = "一室",
    
    ## omega 矩阵（1×1）
    omega = maharaj2021_omega,
    omega_CL = omega_CL_maharaj,
    
    ## 可校准的 eta
    calibratable_etas = "eta_CL",
    
    ## 用于 MAPB 的 omega 子矩阵（只校准 eta_CL）
    omega_calibrate = matrix(omega_CL_maharaj, nrow = 1, ncol = 1,
                             dimnames = list("eta_CL", "eta_CL")),
    
    ## sigma（残差）— 全部用 SD 表示
    sigma = maharaj2021_sigma,
    sigma_type = "combined",
    sigma_prop_sd = sigma_prop_sd_maharaj,
    sigma_add_sd = sigma_add_sd_maharaj,
    
    ## 协变量要求
    required_covariates = c("WT", "PMA"),
    
    ## 模型对象名称
    model_object_name = "maharaj2021_olz_model",
    
    ## 备注
    note = "Maharaj 2021 - 儿童模型（含成熟函数）"
  ),
  
  ## ========== Zang2021（主线-成人/老人）==========
  zang2021 = list(
    model_type = "主线-成人",
    compartment = "一室",
    
    ## omega 矩阵（2×2）
    omega = zang2021_omega,
    omega_CL = omega_CL_zang,
    omega_V = omega_V_zang,
    omega_cov_CL_V = omega_cov_CL_V_zang,
    
    ## 可校准的 eta
    calibratable_etas = c("eta_CL", "eta_V"),
    
    ## 用于 MAPB 的 omega 子矩阵（校准 eta_CL 和 eta_V）
    omega_calibrate = zang2021_omega[c("eta_CL", "eta_V"), c("eta_CL", "eta_V")],
    
    ## sigma（残差）— 比例误差
    sigma = zang2021_sigma,
    sigma_type = "proportional",
    sigma_prop_sd = sigma_prop_sd_zang,
    sigma_add_sd = NA_real_,
    
    ## 协变量要求
    required_covariates = c("SEX", "SMK", "INF", "VAL", "PER", "SER", "FLU", "DGH"),
    
    ## 模型对象名称
    model_object_name = "zang2021_olz_model",
    
    ## 备注
    note = "Zang 2021 - 中国成人（一室模型）"
  ),
  
  ## ========== Li2018（支线-成人）==========
  li2018 = list(
    model_type = "支线-成人",
    compartment = "两室",
    
    ## omega 矩阵（完整）
    omega = li2018_omega_drug2,
    omega_CL = omega_CL_li2018,
    omega_Vc = omega_Vc_li2018,
    omega_cov_CL_Vc = omega_cov_CL_Vc_li2018,
    
    ## 可校准的 eta
    calibratable_etas = c("eta_CL", "eta_Vc"),
    
    ## 用于 MAPB 的 omega 子矩阵（校准 eta_CL 和 eta_Vc，含协方差）
    omega_calibrate = li2018_omega_drug2[c("eta_CL", "eta_Vc"), c("eta_CL", "eta_Vc")],
    
    ## 所有 eta 名称
    all_etas = rownames(li2018_omega_drug2),
    
    ## sigma（残差）— 组合误差
    sigma = li2018_sigma_drug2,
    sigma_type = "combined",
    sigma_prop_sd = sigma_prop_sd_li2018,
    sigma_add_sd = sigma_add_sd_li2018,
    
    ## 协变量要求
    required_covariates = c("WT"),
    
    ## 模型对象名称
    model_object_name = "li2018_olz_model_drug2",
    
    ## 备注
    note = "Li 2018 - 中国汉族成人（两室模型，制剂#2）"
  ),
  
  ## ========== Yin2016（支线-成人）==========
  yin2016 = list(
    model_type = "支线-成人",
    compartment = "两室",
    
    ## omega 矩阵（手动构建，因为 lotri 中 THETA_F 无效）
    omega = matrix(
      c(omega_CL_yin2016, omega_cov_CL_V2_yin2016,
        omega_cov_CL_V2_yin2016, omega_V2_yin2016),
      nrow = 2, byrow = TRUE,
      dimnames = list(c("eta_CL", "eta_V2"), c("eta_CL", "eta_V2"))
    ),
    omega_CL = omega_CL_yin2016,
    omega_V2 = omega_V2_yin2016,
    omega_cov_CL_V2 = omega_cov_CL_V2_yin2016,
    
    ## 可校准的 eta
    calibratable_etas = c("eta_CL", "eta_V2"),
    
    ## 用于 MAPB 的 omega 子矩阵
    omega_calibrate = matrix(
      c(omega_CL_yin2016, omega_cov_CL_V2_yin2016,
        omega_cov_CL_V2_yin2016, omega_V2_yin2016),
      nrow = 2, byrow = TRUE,
      dimnames = list(c("eta_CL", "eta_V2"), c("eta_CL", "eta_V2"))
    ),
    
    ## 所有 eta 名称
    all_etas = c("eta_CL", "eta_V2", "eta_KA", "eta_Q", "eta_V3", "eta_TLAG"),
    
    ## sigma（残差）— 组合误差
    sigma = yin2016_sigma,
    sigma_type = "combined",
    sigma_prop_sd = sigma_prop_sd_yin2016,
    sigma_add_sd = sigma_add_sd_yin2016,
    
    ## 协变量要求
    required_covariates = c("SEX", "SMK"),
    
    ## 模型对象名称
    model_object_name = "yin2016_olz_model",
    
    ## 备注
    note = "Yin 2016 - 中国精神病患者（两室模型）"
  ),
  
  ## ========== Sun2021（支线-成人）==========
  sun2021 = list(
    model_type = "支线-成人",
    compartment = "两室",
    
    ## omega 矩阵
    omega = sun2021_omega,
    omega_CL = omega_CL_sun2021,
    omega_Vc = omega_Vc_sun2021,
    
    ## 可校准的 eta
    calibratable_etas = c("eta_CL", "eta_Vc"),
    
    ## 用于 MAPB 的 omega 子矩阵（对角，无协方差）
    omega_calibrate = diag(c(omega_CL_sun2021, omega_Vc_sun2021)),
    
    ## 所有 eta 名称
    all_etas = rownames(sun2021_omega),
    
    ## sigma（残差）— 比例误差
    sigma = sun2021_sigma,
    sigma_type = "proportional",
    sigma_prop_sd = sigma_prop_sd_sun2021,
    sigma_add_sd = NA_real_,
    
    ## 协变量要求
    required_covariates = c("AGE", "SEX", "WT", "SMOKE", "RACE_BLACK", 
                            "RIF", "HEP_MOD", "REN_SEV", "FED"),
    
    ## 模型对象名称
    model_object_name = "sun2021_olz_model",
    
    ## 备注
    note = "Sun 2021 - 奥氮平+沙米多芬（两室模型）"
  )
)

## 修复 sun2021 的 omega_calibrate 维度名称
dimnames(omega_sigma_bank$sun2021$omega_calibrate) <- list(
  c("eta_CL", "eta_Vc"), c("eta_CL", "eta_Vc")
)

## 打印参数银行摘要
cat("/n=== 参数银行（omega_sigma_bank）摘要 ===/n")
for (model_name in names(omega_sigma_bank)) {
  bank <- omega_sigma_bank[[model_name]]
  cat("/n【", model_name, "】/n", sep = "")
  cat("  类型:", bank$model_type, "|", bank$compartment, "/n")
  cat("  omega_CL:", round(bank$omega_CL, 4), "/n")
  
  if (!is.null(bank$omega_V) && !is.na(bank$omega_V)) {
    cat("  omega_V:", round(bank$omega_V, 4), "/n")
  }
  if (!is.null(bank$omega_Vc) && !is.na(bank$omega_Vc)) {
    cat("  omega_Vc:", round(bank$omega_Vc, 4), "/n")
  }
  if (!is.null(bank$omega_V2) && !is.na(bank$omega_V2)) {
    cat("  omega_V2:", round(bank$omega_V2, 4), "/n")
  }
  if (!is.null(bank$omega_cov_CL_V) && !is.na(bank$omega_cov_CL_V)) {
    cat("  omega_cov(CL,V):", round(bank$omega_cov_CL_V, 4), "/n")
  }
  if (!is.null(bank$omega_cov_CL_Vc) && !is.na(bank$omega_cov_CL_Vc)) {
    cat("  omega_cov(CL,Vc):", round(bank$omega_cov_CL_Vc, 4), "/n")
  }
  if (!is.null(bank$omega_cov_CL_V2) && !is.na(bank$omega_cov_CL_V2)) {
    cat("  omega_cov(CL,V2):", round(bank$omega_cov_CL_V2, 4), "/n")
  }
  
  cat("  可校准 eta:", paste(bank$calibratable_etas, collapse = ", "), "/n")
  cat("  残差类型:", bank$sigma_type, "/n")
  cat("  sigma_prop_sd:", round(bank$sigma_prop_sd, 4))
  if (!is.na(bank$sigma_add_sd)) {
    cat(" | sigma_add_sd:", round(bank$sigma_add_sd, 4))
  }
  cat("/n")
}

message("/n✅ 参数银行构建完成")

## ----------------------------------------------------------------------------
## 1.5 构建模型信息汇总表
## ----------------------------------------------------------------------------
message("/n--- 1.5 构建模型信息汇总表 ---")

model_info <- tibble(
  model_name = names(omega_sigma_bank),
  model_type = sapply(omega_sigma_bank, function(x) x$model_type),
  compartment = sapply(omega_sigma_bank, function(x) x$compartment),
  sigma_type = sapply(omega_sigma_bank, function(x) x$sigma_type),
  omega_CL = sapply(omega_sigma_bank, function(x) x$omega_CL),
  sigma_prop_sd = sapply(omega_sigma_bank, function(x) x$sigma_prop_sd),
  sigma_add_sd = sapply(omega_sigma_bank, function(x) 
    if (is.na(x$sigma_add_sd)) NA_real_ else x$sigma_add_sd),
  calibratable_etas = sapply(omega_sigma_bank, function(x) 
    paste(x$calibratable_etas, collapse = ", ")),
  note = sapply(omega_sigma_bank, function(x) x$note)
)

cat("/n=== 模型信息汇总 ===/n")
print(model_info)

## ----------------------------------------------------------------------------
## 1.6 简单工程检查：确认模型可运行
## ----------------------------------------------------------------------------
message("/n--- 1.6 工程检查：模型可运行性 ---")

## 辅助函数：创建简单事件表
make_simple_ev <- function(dose_mg, n_dose = 3, ii_h = 24, id = 1) {
  dose_times <- seq(0, by = ii_h, length.out = n_dose)
  t_end <- max(dose_times) + ii_h
  time_grid <- seq(0, t_end, by = 1)
  
  dose_df <- data.frame(time = dose_times, amt = dose_mg, evid = 1, cmt = 1)
  obs_df <- data.frame(time = time_grid, amt = 0, evid = 0, cmt = 0)
  
  ev_df <- rbind(dose_df, obs_df)
  ev_df <- ev_df[order(ev_df$time, -ev_df$evid), ]
  ev_df$id <- id
  ev_df
}

## 检查 Zang 模型
test_ev <- make_simple_ev(dose_mg = 10, id = 1)
test_sim <- rxSolve(
  object = zang2021_olz_model,
  params = c(eta_CL = 0, eta_V = 0),
  events = test_ev,
  iCov = data.frame(id = 1, SEX = 1, SMK = 0, INF = 0, VAL = 0, 
                    PER = 0, SER = 0, FLU = 0, DGH = 0)
)
stopifnot(nrow(test_sim) > 0)
message("  ✓ zang2021_olz_model 运行正常")

## 检查 Maharaj 模型
test_ev_ped <- make_simple_ev(dose_mg = 5, id = 1)
test_sim_ped <- rxSolve(
  object = maharaj2021_olz_model,
  params = c(eta_CL = 0),
  events = test_ev_ped,
  iCov = data.frame(id = 1, WT = 49.3, PMA = 717.8)
)
stopifnot(nrow(test_sim_ped) > 0)
message("  ✓ maharaj2021_olz_model 运行正常")

## 检查 Li2018 模型
test_sim_li <- rxSolve(
  object = li2018_olz_model_drug2,
  params = c(eta_CL = 0, eta_Vc = 0, eta_KA = 0, eta_Vp = 0, eta_Q = 0),
  events = test_ev,
  iCov = data.frame(id = 1, WT = 65)
)
stopifnot(nrow(test_sim_li) > 0)
message("  ✓ li2018_olz_model_drug2 运行正常")

## 检查 Yin2016 模型
test_sim_yin <- rxSolve(
  object = yin2016_olz_model,
  params = c(eta_CL = 0, eta_V2 = 0, eta_KA = 0, eta_Q = 0, eta_V3 = 0, eta_TLAG = 0),
  events = test_ev,
  iCov = data.frame(id = 1, SEX = 1, SMK = 0)
)
stopifnot(nrow(test_sim_yin) > 0)
message("  ✓ yin2016_olz_model 运行正常")

## 检查 Sun2021 模型
test_sim_sun <- rxSolve(
  object = sun2021_olz_model,
  params = c(eta_CL = 0, eta_Vc = 0, eta_KA = 0, eta_KA_IOV = 0, 
             eta_Vp = 0, eta_Q = 0, eta_ALAG = 0),
  events = test_ev,
  iCov = data.frame(id = 1, AGE = 22, SEX = 0, WT = 65, SMOKE = 0, 
                    RACE_BLACK = 0, RIF = 0, HEP_MOD = 0, REN_SEV = 0, FED = 0)
)
stopifnot(nrow(test_sim_sun) > 0)
message("  ✓ sun2021_olz_model 运行正常")

message("/n✅ 所有模型工程检查通过")

## ----------------------------------------------------------------------------
## 1.7 输出文件
## ----------------------------------------------------------------------------
message("/n--- 1.7 保存输出文件 ---")

## 保存模型信息表
write_csv(model_info, file.path(OUT_DIR_S1, "model_info.csv"))
saveRDS(model_info, file.path(OUT_DIR_S1, "model_info.rds"))
message("  ✓ model_info.csv / .rds 已保存")

## 保存参数银行（核心输出！）
saveRDS(omega_sigma_bank, file.path(OUT_DIR_S1, "omega_sigma_bank.rds"))
message("  ✓ omega_sigma_bank.rds 已保存（参数银行）")

## 保存模型文件路径映射
saveRDS(model_files, file.path(OUT_DIR_S1, "model_file_paths.rds"))
message("  ✓ model_file_paths.rds 已保存")

## 保存全局设定
global_settings <- list(
  MASTER_SEED = MASTER_SEED,
  PROJ_ROOT = PROJ_ROOT,
  MODEL_DIR = MODEL_DIR,
  OUTPUT_ROOT = OUTPUT_ROOT,
  N_MC = N_MC,
  version = "v6.0",
  created = Sys.time()
)
saveRDS(global_settings, file.path(OUT_DIR_S1, "global_settings.rds"))
message("  ✓ global_settings.rds 已保存")

## ----------------------------------------------------------------------------
## 1.8 Stage 0 & 1 完成
## ----------------------------------------------------------------------------

message("/n", paste(rep("=", 70), collapse = ""))
message("✅ Stage 0 & Stage 1 全部完成（第六版）")
message(paste(rep("=", 70), collapse = ""))

message("/n📁 输出目录: ", OUT_DIR_S1)
message("/n📄 输出文件：")
message("  - model_info.csv / .rds     （模型信息汇总表）")
message("  - omega_sigma_bank.rds      （参数银行 ← 核心文件）")
message("  - model_file_paths.rds      （模型文件路径映射）")
message("  - global_settings.rds       （全局设定）")

message("/n📦 参数银行包含 5 个��型：")
message("  主线：maharaj2021（儿童）、zang2021（成人/老人）")
message("  支线：li2018、yin2016、sun2021")

message("/n🔗 后续脚本使用方法：")
message('  bank <- readRDS("outputs6/outputs6_stage1/omega_sigma_bank.rds")')
message('  omega_calibrate <- bank$zang2021$omega_calibrate  # 用于 MAPB')
message('  sigma_prop_sd <- bank$zang2021$sigma_prop_sd')

message("/n⏭️  下一步：运行 [第六版代码-主线2.R]")
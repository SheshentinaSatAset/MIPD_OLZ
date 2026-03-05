## ############################################################################
## 第六版代码-主线2.R
## 
## Stage 2：典型患��、给药事件的设置
## 
## 核心变化（第六版 vs 第五版）：
##   - 第五版：每个典型患者设 3-5 个固定的 η_true（baseline/highCL/lowCL/highV/lowV）
##   - 第六版：每个典型患者代表一个亚群，η_true 在 Stage 3 从 Ω 矩阵抽样 N=1000
##   - Stage 2 只定义典型患者的协变量，不再预设 η 的分层
## 
## 核心输出：
##   - outputs6/outputs6_stage2/patients_main.rds / .csv   （主线患者表）
##   - outputs6/outputs6_stage2/patients_side.rds / .csv   （支线患者表）
##   - outputs6/outputs6_stage2/regimens.rds / .csv        （给药方案）
##   - outputs6/outputs6_stage2/tdm_rule.rds               （TDM 采样规则）
##   - outputs6/outputs6_stage2/design_manifest_stage2.rds （完整设计清单）
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("Stage 2：典型患者、给药事件的设置（第六版）")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## 2.0 加载依赖与读取 Stage 1 输出
## ############################################################################

message("\n--- 2.0 加载依赖与读取 Stage 1 输出 ---")

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

## 路径设置（与 Stage 0&1 一致）
PROJ_ROOT <- "D:/Users/YujiaZhang/Desktop/26救赎之道/25.9.29 毕业设计/5正式实施！/多成人模型/体现我卓越的项目管理能力"
CODE_DIR <- file.path(PROJ_ROOT, "代码这边请/正式代码/第六版代码")
MODEL_DIR <- file.path(PROJ_ROOT, "代码这边请/poppk模型")
OUTPUT_ROOT <- file.path(PROJ_ROOT, "outputs6")
OUT_DIR_S1 <- file.path(OUTPUT_ROOT, "outputs6_stage1")
OUT_DIR_S2 <- file.path(OUTPUT_ROOT, "outputs6_stage2")

## 创建输出目录
if (!dir.exists(OUT_DIR_S2)) {
  dir.create(OUT_DIR_S2, recursive = TRUE)
  message("📁 已创建目录: ", OUT_DIR_S2)
}

## 读取 Stage 1 输出
global_settings <- readRDS(file.path(OUT_DIR_S1, "global_settings.rds"))
omega_sigma_bank <- readRDS(file.path(OUT_DIR_S1, "omega_sigma_bank.rds"))
model_files <- readRDS(file.path(OUT_DIR_S1, "model_file_paths.rds"))

MASTER_SEED <- global_settings$MASTER_SEED
N_MC <- global_settings$N_MC

set.seed(MASTER_SEED + 2)

message("✅ Stage 1 输出已加载")
message("🎲 全局随机种子: ", MASTER_SEED)
message("🔢 蒙特卡洛次数: ", N_MC)

## 加载模型（用于后续计算 θ_pop）
source(model_files$zang2021)
source(model_files$maharaj2021)
source(model_files$li2018)
source(model_files$yin2016)
source(model_files$sun2021)

message("✅ 5 个 PopPK 模型已加载")

## ############################################################################
## 2.1 定义辅助函数：计算群体典型值 (θ_pop)
## ############################################################################

message("\n--- 2.1 定义 θ_pop 计算函数 ---")

## ========== Maharaj2021 模型（儿童）==========
calc_theta_pop_maharaj <- function(WT, PMA) {
  ## 模型参数
  TVCL <- 16.8
  TVV <- 663
  B_WT_EXP <- 0.486
  TM50 <- 70
  HILL <- 3.97
  
  ## 成熟函数
  MAT <- (PMA^HILL) / (TM50^HILL + PMA^HILL)
  
  ## 计算 θ_pop（不含 eta）
  CL_pop <- TVCL * (WT / 70)^B_WT_EXP * MAT
  V_pop <- TVV * (WT / 70)
  
  list(CL_pop = CL_pop, V_pop = V_pop)
}

## ========== Zang2021 模型（成人/老人）==========
calc_theta_pop_zang <- function(SEX, SMK, INF = 0, VAL = 0, 
                                 PER = 0, SER = 0, FLU = 0, DGH = 0) {
  ## 模型参数
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
  
  ## 计算 CL_pop
  CL_pop <- TVCL * exp(
    B_SEX * SEX + B_SMK * SMK + B_INF * INF + B_VAL * VAL +
      B_PER * PER + B_SER * SER + B_FLU * FLU + B_DGH * DGH
  )
  
  ## V_pop（无协变量效应）
  V_pop <- TVV
  
  list(CL_pop = CL_pop, V_pop = V_pop)
}

## ========== Li2018 模型（支线-成人）==========
calc_theta_pop_li2018 <- function(WT) {
  TVCL <- 25.4
  TVVC <- 2390
  TVQ <- 8.41
  TVVP <- 168
  WT_REF <- 60.59
  B_WT_VC_EXP <- 0.579
  
  CL_pop <- TVCL
  Vc_pop <- TVVC * (WT / WT_REF)^B_WT_VC_EXP
  Q_pop <- TVQ
  Vp_pop <- TVVP
  
  list(CL_pop = CL_pop, Vc_pop = Vc_pop, Q_pop = Q_pop, Vp_pop = Vp_pop)
}

## ========== Yin2016 模型（支线-成人）==========
calc_theta_pop_yin2016 <- function(SEX, SMK) {
  TVCL <- 16.6
  TVV2 <- 599
  TVQ <- 21.3
  TVV3 <- 174
  B_SEX <- 0.385
  B_SMK <- 0.426
  
  CL_pop <- TVCL * exp(B_SEX * SEX + B_SMK * SMK)
  V2_pop <- TVV2
  Q_pop <- TVQ
  V3_pop <- TVV3
  
  list(CL_pop = CL_pop, V2_pop = V2_pop, Q_pop = Q_pop, V3_pop = V3_pop)
}

## ========== Sun2021 模型（支线-成人）==========
calc_theta_pop_sun2021 <- function(AGE, SEX, WT, SMOKE, RACE_BLACK,
                                    RIF, HEP_MOD, REN_SEV, FED) {
  TVCL <- 15.5
  TVVC <- 656
  TVQ <- 6.15
  TVVP <- 225
  
  WT_REF <- 70
  AGE_REF <- 36
  
  B_WT_CL <- 0.75
  B_WT_VC <- 1.0
  B_AGE_VC <- 0.356
  B_RIF <- log(1.80)
  B_SMK <- log(1.30)
  B_SEX <- log(0.862)
  B_RACE_B <- log(1.10)
  B_HEP_MOD <- log(0.875)
  B_REN_SEV <- log(0.801)
  B_FED <- log(0.943)
  
  CL_pop <- TVCL * exp(
    B_WT_CL * log(WT / WT_REF) +
      B_RIF * RIF +
      B_SMK * SMOKE +
      B_SEX * SEX +
      B_RACE_B * RACE_BLACK +
      B_HEP_MOD * HEP_MOD +
      B_REN_SEV * REN_SEV +
      B_FED * FED
  )
  
  Vc_pop <- TVVC * exp(
    B_WT_VC * log(WT / WT_REF) +
      B_AGE_VC * log(AGE / AGE_REF)
  )
  
  Q_pop <- TVQ
  Vp_pop <- TVVP
  
  list(CL_pop = CL_pop, Vc_pop = Vc_pop, Q_pop = Q_pop, Vp_pop = Vp_pop)
}

message("✅ θ_pop 计算函数已定义")

## ############################################################################
## 2.2 定义主线典型患者（3 个亚群）
## ############################################################################

message("\n--- 2.2 定义主线典型患者 ---")

## ========== 2.2.1 儿童患者 (PED) — Maharaj2021 ==========
## 数据来源：Farhath 2024
## 协变量：WT = 49.3 kg, PMA = 717.8 周

theta_pop_ped <- calc_theta_pop_maharaj(WT = 49.3, PMA = 717.8)

patient_PED <- tibble(
  patient_id = "PED",
  pop = "儿童",
  model_name = "maharaj2021",
  
  ## 协变量
  WT = 49.3,
  PMA = 717.8,
  
  ## Zang 模型协变量（儿童不需要，填 NA）
  SEX = NA_real_,
  SMK = NA_real_,
  INF = NA_real_,
  VAL = NA_real_,
  PER = NA_real_,
  SER = NA_real_,
  FLU = NA_real_,
  DGH = NA_real_,
  
  ## 给药信息
  dose_mg = 5,
  ii_h = 24,
  n_dose = 14,
  
  ## θ_pop
  CL_pop = theta_pop_ped$CL_pop,
  V_pop = theta_pop_ped$V_pop,
  
  ## 数据来源
  data_source = "Farhath 2024"
)

cat("\n【儿童患者 (PED)】\n")
cat("  模型: maharaj2021\n")
cat("  协变量: WT =", 49.3, "kg, PMA =", 717.8, "周\n")
cat("  θ_pop: CL =", round(theta_pop_ped$CL_pop, 2), "L/h, V =", round(theta_pop_ped$V_pop, 2), "L\n")

## ========== 2.2.2 成人患者 (ADULT) — Zang2021 ==========
## 数据来源：Wong 2007
## 协变量：男性(SEX=1), 不吸烟(SMK=0), 无感染(INF=0), 合并丙戊酸钠(VAL=1), 其他=0

theta_pop_adult <- calc_theta_pop_zang(
  SEX = 1, SMK = 0, INF = 0, VAL = 1,
  PER = 0, SER = 0, FLU = 0, DGH = 0
)

patient_ADULT <- tibble(
  patient_id = "ADULT",
  pop = "成人",
  model_name = "zang2021",
  
  ## 协变量
  WT = NA_real_,
  PMA = NA_real_,
  SEX = 1,
  SMK = 0,
  INF = 0,
  VAL = 1,
  PER = 0,
  SER = 0,
  FLU = 0,
  DGH = 0,
  
  ## 给药信息
  dose_mg = 10,
  ii_h = 24,
  n_dose = 14,
  
  ## θ_pop
  CL_pop = theta_pop_adult$CL_pop,
  V_pop = theta_pop_adult$V_pop,
  
  ## 数据来源
  data_source = "Wong 2007"
)

cat("\n【成人患者 (ADULT)】\n")
cat("  模型: zang2021\n")
cat("  协变量: 男性(SEX=1), 不吸烟(SMK=0), 无感染(INF=0), 合并丙戊酸钠(VAL=1)\n")
cat("  θ_pop: CL =", round(theta_pop_adult$CL_pop, 2), "L/h, V =", round(theta_pop_adult$V_pop, 2), "L\n")

## ========== 2.2.3 老人患者 (ELDERLY) — Zang2021 ==========
## 数据来源：Iwahashi 2004
## 协变量：女性(SEX=0), 不吸烟(SMK=0), 无感染(INF=0), 无合并用药(全部=0)

theta_pop_elderly <- calc_theta_pop_zang(
  SEX = 0, SMK = 0, INF = 0, VAL = 0,
  PER = 0, SER = 0, FLU = 0, DGH = 0
)

patient_ELDERLY <- tibble(
  patient_id = "ELDERLY",
  pop = "老人",
  model_name = "zang2021",
  
  ## 协变量
  WT = NA_real_,
  PMA = NA_real_,
  SEX = 0,
  SMK = 0,
  INF = 0,
  VAL = 0,
  PER = 0,
  SER = 0,
  FLU = 0,
  DGH = 0,
  
  ## 给药信息
  dose_mg = 10,
  ii_h = 24,
  n_dose = 14,
  
  ## θ_pop
  CL_pop = theta_pop_elderly$CL_pop,
  V_pop = theta_pop_elderly$V_pop,
  
  ## 数据来源
  data_source = "Iwahashi 2004"
)

cat("\n【老人患者 (ELDERLY)】\n")
cat("  模型: zang2021\n")
cat("  协变量: 女性(SEX=0), 不吸烟(SMK=0), 无感染(INF=0), 无合并用药\n")
cat("  θ_pop: CL =", round(theta_pop_elderly$CL_pop, 2), "L/h, V =", round(theta_pop_elderly$V_pop, 2), "L\n")

## ========== 合并主线患者表 ==========
patients_main <- bind_rows(
  patient_PED,
  patient_ADULT,
  patient_ELDERLY
)

message("\n✅ 主线典型患者定义完成（3 个亚群）")

## ############################################################################
## 2.3 定义支线典型患者（1 个成人 × 4 个模型���角）
## ############################################################################

message("\n--- 2.3 定义支线典型患者 ---")

## 支线的"成人患者"与主线的"成人患者"是同一个人
## 特征：22岁，男性，体重65kg，不吸烟，无感染，非黑人，
##       合并用药（丙戊酸钠），肝功能正常，肾功能正常，空腹
## 只是不同模型需要不同的协变量子集

## ========== 2.3.1 ADULT_zang — Zang2021 ==========
## 与主线成人完全相同

theta_pop_adult_zang <- calc_theta_pop_zang(
  SEX = 1, SMK = 0, INF = 0, VAL = 1,
  PER = 0, SER = 0, FLU = 0, DGH = 0
)

patient_ADULT_zang <- tibble(
  patient_id = "ADULT_zang",
  pop = "成人",
  model_name = "zang2021",
  
  ## Zang 模型协变量
  SEX = 1,
  SMK = 0,
  INF = 0,
  VAL = 1,
  PER = 0,
  SER = 0,
  FLU = 0,
  DGH = 0,
  
  ## 其他协变量（Zang 不需要）
  WT = NA_real_,
  PMA = NA_real_,
  AGE = NA_real_,
  SMOKE = NA_real_,
  RACE_BLACK = NA_real_,
  RIF = NA_real_,
  HEP_MOD = NA_real_,
  REN_SEV = NA_real_,
  FED = NA_real_,
  
  ## 给药信息
  dose_mg = 10,
  ii_h = 24,
  n_dose = 14,
  
  ## θ_pop
  CL_pop = theta_pop_adult_zang$CL_pop,
  V_pop = theta_pop_adult_zang$V_pop,
  Vc_pop = NA_real_,
  V2_pop = NA_real_,
  Q_pop = NA_real_,
  Vp_pop = NA_real_,
  V3_pop = NA_real_,
  
  ## 数据来源
  data_source = "Wong 2007 (同主线成人)"
)

## ========== 2.3.2 ADULT_li — Li2018 ==========
## 协变量：WT = 65 kg

theta_pop_adult_li <- calc_theta_pop_li2018(WT = 65)

patient_ADULT_li <- tibble(
  patient_id = "ADULT_li",
  pop = "成人",
  model_name = "li2018",
  
  ## Li2018 模型协变量
  WT = 65,
  
  ## Zang 模型协变量（Li 不需要，但为统一格式保留）
  SEX = NA_real_,
  SMK = NA_real_,
  INF = NA_real_,
  VAL = NA_real_,
  PER = NA_real_,
  SER = NA_real_,
  FLU = NA_real_,
  DGH = NA_real_,
  
  ## 其他协变量
  PMA = NA_real_,
  AGE = NA_real_,
  SMOKE = NA_real_,
  RACE_BLACK = NA_real_,
  RIF = NA_real_,
  HEP_MOD = NA_real_,
  REN_SEV = NA_real_,
  FED = NA_real_,
  
  ## 给药信息
  dose_mg = 10,
  ii_h = 24,
  n_dose = 14,
  
  ## θ_pop（两室模型）
  CL_pop = theta_pop_adult_li$CL_pop,
  Vc_pop = theta_pop_adult_li$Vc_pop,
  Q_pop = theta_pop_adult_li$Q_pop,
  Vp_pop = theta_pop_adult_li$Vp_pop,
  V_pop = NA_real_,
  V2_pop = NA_real_,
  V3_pop = NA_real_,
  
  ## 数据来源
  data_source = "同主线成人"
)

## ========== 2.3.3 ADULT_yin — Yin2016 ==========
## 协变量：SEX = 1（男）, SMK = 0（不吸烟）

theta_pop_adult_yin <- calc_theta_pop_yin2016(SEX = 1, SMK = 0)

patient_ADULT_yin <- tibble(
  patient_id = "ADULT_yin",
  pop = "成人",
  model_name = "yin2016",
  
  ## Yin2016 模型协变量
  SEX = 1,
  SMK = 0,
  
  ## 其他 Zang 协变量
  INF = NA_real_,
  VAL = NA_real_,
  PER = NA_real_,
  SER = NA_real_,
  FLU = NA_real_,
  DGH = NA_real_,
  
  ## 其他协变量
  WT = NA_real_,
  PMA = NA_real_,
  AGE = NA_real_,
  SMOKE = NA_real_,
  RACE_BLACK = NA_real_,
  RIF = NA_real_,
  HEP_MOD = NA_real_,
  REN_SEV = NA_real_,
  FED = NA_real_,
  
  ## 给药信息
  dose_mg = 10,
  ii_h = 24,
  n_dose = 14,
  
  ## θ_pop（两室模型，Yin 用 V2/V3）
  CL_pop = theta_pop_adult_yin$CL_pop,
  V2_pop = theta_pop_adult_yin$V2_pop,
  Q_pop = theta_pop_adult_yin$Q_pop,
  V3_pop = theta_pop_adult_yin$V3_pop,
  V_pop = NA_real_,
  Vc_pop = NA_real_,
  Vp_pop = NA_real_,
  
  ## 数据来源
  data_source = "同主线成人"
)

## ========== 2.3.4 ADULT_sun — Sun2021 ==========
## 协变量：22岁，男性(SEX=0，注意编码!)，体重65kg，不吸烟，非黑人，
##         无利福平，肝肾正常，空腹
## 注意：Sun2021 模型的 SEX 编码是 0=男性, 1=女性，与 Zang/Yin 相反！

theta_pop_adult_sun <- calc_theta_pop_sun2021(
  AGE = 22, SEX = 0, WT = 65, SMOKE = 0, RACE_BLACK = 0,
  RIF = 0, HEP_MOD = 0, REN_SEV = 0, FED = 0
)

patient_ADULT_sun <- tibble(
  patient_id = "ADULT_sun",
  pop = "成人",
  model_name = "sun2021",
  
  ## Sun2021 模型协变量
  AGE = 22,
  SEX = 0,         # 注意：Sun2021 中 0=男性！
  WT = 65,
  SMOKE = 0,
  RACE_BLACK = 0,
  RIF = 0,
  HEP_MOD = 0,
  REN_SEV = 0,
  FED = 0,
  
  ## Zang/Yin 模型协变量（Sun 不需要）
  SMK = NA_real_,
  INF = NA_real_,
  VAL = NA_real_,
  PER = NA_real_,
  SER = NA_real_,
  FLU = NA_real_,
  DGH = NA_real_,
  
  ## Maharaj 协变量
  PMA = NA_real_,
  
  ## 给药信息
  dose_mg = 10,
  ii_h = 24,
  n_dose = 14,
  
  ## θ_pop（两室模型）
  CL_pop = theta_pop_adult_sun$CL_pop,
  Vc_pop = theta_pop_adult_sun$Vc_pop,
  Q_pop = theta_pop_adult_sun$Q_pop,
  Vp_pop = theta_pop_adult_sun$Vp_pop,
  V_pop = NA_real_,
  V2_pop = NA_real_,
  V3_pop = NA_real_,
  
  ## 数据来源
  data_source = "同主线成人"
)

## ========== 合并支线患者表 ==========
patients_side <- bind_rows(
  patient_ADULT_zang,
  patient_ADULT_li,
  patient_ADULT_yin,
  patient_ADULT_sun
)

cat("\n=== 支线典型患者（同一成人 × 4 个模型视角）===\n")
cat("\n【ADULT_zang】\n")
cat("  模型: zang2021 | θ_pop: CL =", round(theta_pop_adult_zang$CL_pop, 2), "L/h\n")
cat("\n【ADULT_li】\n")
cat("  模型: li2018 | WT = 65 kg | θ_pop: CL =", round(theta_pop_adult_li$CL_pop, 2), "L/h, Vc =", round(theta_pop_adult_li$Vc_pop, 2), "L\n")
cat("\n【ADULT_yin】\n")
cat("  模型: yin2016 | SEX=1(男), SMK=0 | θ_pop: CL =", round(theta_pop_adult_yin$CL_pop, 2), "L/h, V2 =", round(theta_pop_adult_yin$V2_pop, 2), "L\n")
cat("\n【ADULT_sun】\n")
cat("  模型: sun2021 | 22岁, SEX=0(男), 65kg | θ_pop: CL =", round(theta_pop_adult_sun$CL_pop, 2), "L/h, Vc =", round(theta_pop_adult_sun$Vc_pop, 2), "L\n")

cat("\n⚠️  注意：Sun2021 模型 SEX 编码：0=男性, 1=女性（与 Zang/Yin 相反！）\n")

message("\n✅ 支线典型患者定义完成（1 个成人 × 4 个模型视角）")

## ############################################################################
## 2.4 定义给药方案
## ############################################################################

message("\n--- 2.4 定义给药方案 ---")

regimens <- tibble(
  pop = c("儿童", "成人", "老人"),
  dose_mg = c(5, 10, 10),
  ii_h = c(24, 24, 24),
  n_dose = c(14, 14, 14),
  unit = "mg",
  route = "口服",
  note = c(
    "儿童剂量参考 Farhath 2024",
    "成人标准剂量",
    "老人标准剂量"
  )
)

cat("\n=== 给药方案 ===\n")
print(regimens)

## ############################################################################
## 2.5 定义 TDM 采样规则
## ############################################################################

message("\n--- 2.5 定义 TDM 采样规则 ---")

## 采样时间计算：
## 第 14 次给药时间：t = (14-1) × 24 = 312 h
## 采样时间（末次给药后 24h 谷浓度）：t = 312 + 24 = 336 h

tdm_rule <- list(
  type = "single_trough",
  description = "末次给��后的 24h 谷浓度（第 14 次给药后）",
  time_relative_to_last_dose_h = 24,
  t_last_dose_h = (14 - 1) * 24,  # = 312 h
  t_sample_h = (14 - 1) * 24 + 24,  # = 336 h
  n_samples = 1
)

cat("\n=== TDM 采样规则 ===\n")
cat("类型:", tdm_rule$type, "\n")
cat("描述:", tdm_rule$description, "\n")
cat("末次给药时间:", tdm_rule$t_last_dose_h, "h\n")
cat("采样时间:", tdm_rule$t_sample_h, "h\n")

## ############################################################################
## 2.6 构建设计清单（Single Source of Truth）
## ############################################################################

message("\n--- 2.6 构建设计清单 ---")

design_manifest_stage2 <- list(
  
  ## ========== 典型患者 ==========
  patients_main = patients_main,
  patients_side = patients_side,
  
  ## ========== 给药方案 ==========
  regimens = regimens,
  
  ## ========== TDM 采样规则 ==========
  tdm_rule = tdm_rule,
  
  ## ========== 模型-患者映射（主线）==========
  model_patient_mapping_main = tibble(
    patient_id = c("PED", "ADULT", "ELDERLY"),
    model_name = c("maharaj2021", "zang2021", "zang2021"),
    calibratable_etas = c(
      "eta_CL",
      "eta_CL, eta_V",
      "eta_CL, eta_V"
    )
  ),
  
  ## ========== 模型-患者映射（支线）==========
  model_patient_mapping_side = tibble(
    patient_id = c("ADULT_zang", "ADULT_li", "ADULT_yin", "ADULT_sun"),
    model_name = c("zang2021", "li2018", "yin2016", "sun2021"),
    calibratable_etas = c(
      "eta_CL, eta_V",
      "eta_CL, eta_Vc",
      "eta_CL, eta_V2",
      "eta_CL, eta_Vc"
    ),
    has_covariance = c(FALSE, TRUE, TRUE, FALSE)
  ),
  
  ## ========== 全局约定 ==========
  conventions = list(
    concentration_unit = "ng/mL",
    dose_unit = "mg",
    time_unit = "h",
    sex_coding_zang_yin = "1=男, 0=女",
    sex_coding_sun = "0=男, 1=女（注意与 Zang/Yin 相反！）",
    n_mc = N_MC
  ),
  
  ## ========== 备注 ==========
  notes = list(
    design_principle = "第六版核心变化：每个典型患者代表一个亚群，η_true 在 Stage 3 从 Ω 抽样 N=1000",
    main_line = "主线：maharaj2021→儿童；zang2021→成人/老人",
    side_line = "支线：同一成人患者在 4 个模型中的表现（验证一致性）",
    adult_val = "成人患者合并丙戊酸钠（VAL=1），会影响 CL_pop"
  ),
  
  ## ========== 元数据 ==========
  metadata = list(
    created = Sys.time(),
    master_seed = MASTER_SEED,
    version = "v6.0",
    stage = "Stage 2"
  )
)

message("✅ 设计清单构建完成")

## ############################################################################
## 2.7 输出文件
## ############################################################################

message("\n--- 2.7 保存输出文件 ---")

## 保存主线患者表
write_csv(patients_main, file.path(OUT_DIR_S2, "patients_main.csv"))
saveRDS(patients_main, file.path(OUT_DIR_S2, "patients_main.rds"))
message("  ✓ patients_main.csv / .rds 已保存")

## 保存支线患者表
write_csv(patients_side, file.path(OUT_DIR_S2, "patients_side.csv"))
saveRDS(patients_side, file.path(OUT_DIR_S2, "patients_side.rds"))
message("  ✓ patients_side.csv / .rds 已保存")

## 保存给药方案
write_csv(regimens, file.path(OUT_DIR_S2, "regimens.csv"))
saveRDS(regimens, file.path(OUT_DIR_S2, "regimens.rds"))
message("  ✓ regimens.csv / .rds 已保存")

## 保存 TDM 采样规则
saveRDS(tdm_rule, file.path(OUT_DIR_S2, "tdm_rule.rds"))
message("  ✓ tdm_rule.rds 已保存")

## 保存完整设计清单（核心输出）
saveRDS(design_manifest_stage2, file.path(OUT_DIR_S2, "design_manifest_stage2.rds"))
message("  ✓ design_manifest_stage2.rds 已保存（设计清单）")

## ############################################################################
## 2.8 诊断汇总
## ############################################################################

cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("Stage 2 诊断汇总\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

cat("\n=== 主线典型患者（3 个亚群）===\n")
print(patients_main %>% 
        select(patient_id, pop, model_name, dose_mg, CL_pop, V_pop, data_source))

cat("\n=== 支线典型患者（1 成人 × 4 模型）===\n")
print(patients_side %>% 
        select(patient_id, model_name, CL_pop, Vc_pop, V2_pop))

cat("\n=== 给药方案 ===\n")
print(regimens %>% select(pop, dose_mg, ii_h, n_dose))

cat("\n=== TDM 采样规则 ===\n")
cat("采样时间: t =", tdm_rule$t_sample_h, "h（末次给药后 24h 谷浓度）\n")

cat("\n=== 第六版核心设计说明 ===\n")
cat("• 每个典型患者（如 PED）代表一个亚群\n")
cat("• Stage 3 中，每个亚群生成 N=1000 个虚拟个体\n")
cat("• 每个虚拟个体的 η_true 从 Ω 矩阵抽样\n")
cat("• θ_true = θ_pop × exp(η_true)\n")
cat("• C_obs = f(θ_true) + ε\n")
cat("• 然后对每个个体执行 MAPB，得到 η_ind 和 θ_ind\n")

## ############################################################################
## 2.9 Stage 2 完成
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("✅ Stage 2 全部完成（第六版）")
message(paste(rep("=", 70), collapse = ""))

message("\n📁 输出目录: ", OUT_DIR_S2)
message("\n📄 输出文件：")
message("  - patients_main.csv / .rds   （主线患者表，3 个亚群）")
message("  - patients_side.csv / .rds   （支线患者表，1 成人 × 4 模型）")
message("  - regimens.csv / .rds        （给药方案）")
message("  - tdm_rule.rds               （TDM 采样规则）")
message("  - design_manifest_stage2.rds （设计清单 ← 核心文件）")

message("\n🔗 后续��本使用方法：")
message('  manifest <- readRDS("outputs6/outputs6_stage2/design_manifest_stage2.rds")')
message('  patients_main <- manifest$patients_main')
message('  tdm_rule <- manifest$tdm_rule')

message("\n⏭️  下一步：运行 [第六版代码-主线3.R]（MAPB 蒙特卡洛模拟）")
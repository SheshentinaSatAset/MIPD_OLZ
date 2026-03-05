## ############################################################################
## 第八版代码-主线4.R
## 
## Stage 4：后验概率计算（分箱蒙特卡洛法）
## 
## 核心方法变化（第八版 vs 第七版）：
##   - 第七版：直接用正态分布 PDF 计算似然（解析法）
##   - 第八版：分箱蒙特卡洛法，构建条件概率库后查表
## 
## 核心参数：
##   - Δ = 1 ng/mL（分箱步长）
##   - M = 5000（每场景模拟次数）
##   - 并行计算（8 性能核 + 16 能效核）
## 
## η 处理规则（同第七版）：
##   1. Acceptable 参数：全部使用 η_ind
##   2. Unacceptable 参数：强制 η = 0（yin2016-V2, sun2021-Vc）
##   3. MAPB 未收敛的患者：η_ind 回退到 0
## 
## 输入路径：outputs6/（Stage 1-3 的结果）
## 输出路径：outputs8/outputs8_stage4/
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("Stage 4：后验概率计算（第八版 - 分箱蒙特卡洛法）")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## 4.0 加载依赖与路径设置
## ############################################################################

message("\n--- 4.0 加载依赖与路径设置 ---")

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
  library(parallel)
})

## 路径设置
PROJ_ROOT <- "D:/Users/YujiaZhang/Desktop/26救赎之道/25.9.29 毕业设计/5正式实施！/多成人模型/体现我卓越的项目管理能力"
CODE_DIR <- file.path(PROJ_ROOT, "代码这边请/正式代码/第八版代码")
MODEL_DIR <- file.path(PROJ_ROOT, "代码这边请/poppk模型")

## 输入路径（第六版输出）
INPUT_ROOT <- file.path(PROJ_ROOT, "outputs6")
IN_DIR_S1 <- file.path(INPUT_ROOT, "outputs6_stage1")
IN_DIR_S2 <- file.path(INPUT_ROOT, "outputs6_stage2")
IN_DIR_S3 <- file.path(INPUT_ROOT, "outputs6_stage3")

## 输出路径（第八版）
OUTPUT_ROOT <- file.path(PROJ_ROOT, "outputs8")
OUT_DIR_S4 <- file.path(OUTPUT_ROOT, "outputs8_stage4")

## 创建输出目录
if (!dir.exists(OUT_DIR_S4)) {
  dir.create(OUT_DIR_S4, recursive = TRUE)
  message("📁 已创建目录: ", OUT_DIR_S4)
}

## ============================================================================
## 核心参数设置
## ============================================================================
DELTA <- 1          # 分箱步长 (ng/mL)
M_SIM <- 20000       # 每场景模拟次数
C_MIN <- 0          # 浓度下限 (ng/mL)
C_MAX <- 300        # 浓度上限 (ng/mL)
N_CORES <- 8        # 并行核心数（使用性能核）

message("\n📊 核心参数：")
message("   分箱步长 Δ = ", DELTA, " ng/mL")
message("   每场景模拟次数 M = ", M_SIM)
message("   浓度范围 = [", C_MIN, ", ", C_MAX, "] ng/mL")
message("   并行核心数 = ", N_CORES)

## 读取前置 Stage 输出
global_settings <- readRDS(file.path(IN_DIR_S1, "global_settings.rds"))
omega_sigma_bank <- readRDS(file.path(IN_DIR_S1, "omega_sigma_bank.rds"))

## 修复：直接使用当前电脑的路径
model_files <- list(
  zang2021    = file.path(MODEL_DIR, "zang2021_olz_model.R"),
  maharaj2021 = file.path(MODEL_DIR, "maharaj2021_olz_model.R"),
  li2018      = file.path(MODEL_DIR, "li2018_olz_model_drug2.R"),
  yin2016     = file.path(MODEL_DIR, "yin2016_olz_model.R"),
  sun2021     = file.path(MODEL_DIR, "sun2021_olz_model.R")
)

design_manifest <- readRDS(file.path(IN_DIR_S2, "design_manifest_stage2.rds"))
patients_main <- design_manifest$patients_main
patients_side <- design_manifest$patients_side
tdm_rule <- design_manifest$tdm_rule

mc_raw_main <- readRDS(file.path(IN_DIR_S3, "mc_raw_results_main.rds"))
mc_raw_side <- readRDS(file.path(IN_DIR_S3, "mc_raw_results_side.rds"))

MASTER_SEED <- global_settings$MASTER_SEED
N_MC <- global_settings$N_MC

set.seed(MASTER_SEED + 4)

message("\n✅ 前置 Stage 输出已加载")
message("🎲 全局随机种子: ", MASTER_SEED)
message("🔢 虚拟患者数: ", N_MC)

## 加载模型
source(model_files$zang2021)
source(model_files$maharaj2021)
source(model_files$li2018)
source(model_files$yin2016)
source(model_files$sun2021)

message("✅ 5 个 PopPK 模型已加载")

## ############################################################################
## 4.1 合并 Stage 3 结果并应用处理规则
## ############################################################################

message("\n--- 4.1 合并 Stage 3 结果并应用处理规则 ---")

mc_all <- bind_rows(mc_raw_main, mc_raw_side)

message("  原始数据行数: ", nrow(mc_all))
message("  患者组合数: ", length(unique(mc_all$patient_id)))

## ========== 规则1：MAPB 未收敛的患者，η_ind 回退到 0 ==========
## 注意：启用 robust_optim 后，此处理论上不再触发
n_not_converged <- sum(!mc_all$converged)
if (n_not_converged > 0) {
  message("\n⚠️ 发现 ", n_not_converged, " 个未收敛的 MAPB 结果")
  message("   注：已启用 robust_optim 三轮兜底，此情况应极其罕见")
  message("   这些患者的 η_ind 将回退到 0（保底机制）")
  
  mc_all <- mc_all %>%
    mutate(
      eta_CL_ind = ifelse(converged, eta_CL_ind, 0),
      eta_V_ind = ifelse(converged, eta_V_ind, 0),
      CL_ind = ifelse(converged, CL_ind, CL_pop),
      V_ind = ifelse(converged, V_ind, V_pop)
    )
} else {
  message("\n✅ 所有 MAPB 优化均已收敛（robust_optim 三轮兜底生效）")
}

## ========== 规则2：Unacceptable 参数强制 η = 0 ==========
mc_all <- mc_all %>%
  mutate(
    eta_V_ind = case_when(
      model_name == "yin2016" ~ 0,
      model_name == "sun2021" ~ 0,
      TRUE ~ eta_V_ind
    ),
    V_ind = case_when(
      model_name == "yin2016" ~ V_pop,
      model_name == "sun2021" ~ V_pop,
      TRUE ~ V_ind
    )
  )

message("\n✅ Unacceptable 参数已处理:")
message("   - yin2016: η_V2 强制 = 0")
message("   - sun2021: η_Vc 强制 = 0")

## ############################################################################
## 4.2 定义依从性场景库
## ############################################################################

message("\n--- 4.2 定义依从性场景库 ---")

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
## 4.3 定义分箱结构
## ############################################################################

message("\n--- 4.3 定义分箱结构 ---")

## 生成箱边界
bin_edges <- seq(C_MIN, C_MAX, by = DELTA)
n_bins <- length(bin_edges) - 1

message("  箱数量: ", n_bins)
message("  箱边界: [", C_MIN, ", ", C_MAX, "] ng/mL")

## 将浓度值映射到箱索引的函数
get_bin_index <- function(C) {
  ## 边界处理
  C <- pmax(C, C_MIN + 1e-6)
  C <- pmin(C, C_MAX - 1e-6)
  
  ## 计算箱索��（1-based）
  idx <- floor((C - C_MIN) / DELTA) + 1
  idx <- pmax(1, pmin(idx, n_bins))
  idx
}

## ############################################################################
## 4.4 定义辅助函数
## ############################################################################

message("\n--- 4.4 定义辅助函数 ---")

## ----------------------------------------------------------------------------
## 4.4.1 构建指定场景的给药事件表
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
## 4.4.2 提取指定时间的浓度
## ----------------------------------------------------------------------------
extract_conc_at_time <- function(sim_df, t_target) {
  idx <- which(abs(sim_df$time - t_target) < 1e-8)
  if (length(idx) == 0) return(NA_real_)
  sim_df$C_ngmL[idx[1]]
}

## ----------------------------------------------------------------------------
## 4.4.3 生成残差样本（向量化）
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
## 4.4.4 计算预测浓度（给定 η 和场景）
## ----------------------------------------------------------------------------
calc_C_pred <- function(model_name, eta_CL, eta_V, iCov, ev_df, t_sample_h) {
  
  if (model_name == "maharaj2021") {
    sim <- rxSolve(
      object = maharaj2021_olz_model,
      params = c(eta_CL = eta_CL),
      events = ev_df,
      iCov = iCov,
      returnType = "data.frame"
    )
  } else if (model_name == "zang2021") {
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
## 4.4.5 构建条件概率库（核心函数）
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
    
    ## 构建场景 j 的给药事件表
    ev_j <- make_ev_for_scenario(
      dose_mg = dose_mg, ii_h = ii_h, n_dose = n_dose,
      t_sample_h = t_sample_h,
      dose_13_mult = scenario_lib$dose_13_mult[j],
      dose_14_mult = scenario_lib$dose_14_mult[j]
    )
    
    ## 计算预测浓度
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
      ## 如果预测浓度无效，设置均匀分布
      cond_prob_matrix[j, ] <- 1 / n_bins
      next
    }
    
    ## 生成 M 个模拟浓度
    C_sim <- generate_residual_samples(
      C_pred = C_pred,
      sigma_type = sigma_type,
      sigma_prop_sd = sigma_prop_sd,
      sigma_add_sd = sigma_add_sd,
      n_samples = M_sim
    )
    
    ## 分箱统计
    bin_indices <- floor((pmax(pmin(C_sim, C_MIN + n_bins * DELTA - 1e-6), C_MIN + 1e-6) - C_MIN) / DELTA) + 1
    bin_indices <- pmax(1, pmin(bin_indices, n_bins))
    
    bin_counts <- tabulate(bin_indices, nbins = n_bins)
    
    ## 计算条件概率（加拉普拉斯平滑避免零概率）
    cond_prob_matrix[j, ] <- (bin_counts + 1) / (M_sim + n_bins)
  }
  
  list(
    cond_prob = cond_prob_matrix,
    C_pred = C_pred_vec
  )
}

## ----------------------------------------------------------------------------
## 4.4.6 计算后验概率（从条件概率库查表）
## ----------------------------------------------------------------------------
calc_posterior_from_library <- function(C_obs, cond_prob_matrix, n_bins, DELTA, C_MIN) {
  
  ## 确定 C_obs 落入的箱
  bin_idx <- floor((max(min(C_obs, C_MIN + n_bins * DELTA - 1e-6), C_MIN + 1e-6) - C_MIN) / DELTA) + 1
  bin_idx <- max(1, min(bin_idx, n_bins))
  
  ## 从条件概率矩阵中提取该箱的条件概率
  cond_probs <- cond_prob_matrix[, bin_idx]
  
  ## 假设均匀先验
  prior <- rep(1 / nrow(cond_prob_matrix), nrow(cond_prob_matrix))
  
  ## 计算后验概率（贝叶斯公式）
  numerator <- cond_probs * prior
  posterior <- numerator / sum(numerator)
  
  posterior
}

message("✅ 辅助函数已定义")

## ############################################################################
## 4.5 定义单患者处理函数（用于并行）
## ############################################################################

message("\n--- 4.5 定义单患者处理函数 ---")

process_single_patient <- function(patient_id_target, mc_all, patients_all,
                                    omega_sigma_bank, scenario_library,
                                    t_sample_h, n_bins, M_SIM, DELTA, C_MIN,
                                    MASTER_SEED) {
  
  ## 设置该患者的随机种子（确保可重复性）
  set.seed(MASTER_SEED + 4 + which(unique(mc_all$patient_id) == patient_id_target))
  
  ## 获取患者信息
  pat_info <- patients_all %>% filter(patient_id == patient_id_target)
  model_name <- pat_info$model_name[1]
  model_bank <- omega_sigma_bank[[model_name]]
  
  ## 给药信息
  dose_mg <- pat_info$dose_mg[1]
  ii_h <- pat_info$ii_h[1]
  n_dose <- pat_info$n_dose[1]
  
  ## ��备协变量
  if (model_name == "maharaj2021") {
    iCov <- data.frame(id = 1, WT = pat_info$WT[1], PMA = pat_info$PMA[1])
  } else if (model_name == "zang2021") {
    iCov <- data.frame(
      id = 1, SEX = pat_info$SEX[1], SMK = pat_info$SMK[1],
      INF = pat_info$INF[1], VAL = pat_info$VAL[1],
      PER = pat_info$PER[1], SER = pat_info$SER[1],
      FLU = pat_info$FLU[1], DGH = pat_info$DGH[1]
    )
  } else if (model_name == "li2018") {
    iCov <- data.frame(id = 1, WT = pat_info$WT[1])
  } else if (model_name == "yin2016") {
    iCov <- data.frame(id = 1, SEX = pat_info$SEX[1], SMK = pat_info$SMK[1])
  } else if (model_name == "sun2021") {
    iCov <- data.frame(
      id = 1, AGE = pat_info$AGE[1], SEX = pat_info$SEX[1],
      WT = pat_info$WT[1], SMOKE = pat_info$SMOKE[1],
      RACE_BLACK = pat_info$RACE_BLACK[1], RIF = pat_info$RIF[1],
      HEP_MOD = pat_info$HEP_MOD[1], REN_SEV = pat_info$REN_SEV[1],
      FED = pat_info$FED[1]
    )
  }
  
  ## 获取该患者的蒙特卡洛数据
  mc_pat <- mc_all %>% filter(patient_id == patient_id_target)
  n_mc <- nrow(mc_pat)
  n_scenarios <- nrow(scenario_library)
  
  ## 初始化结果存储
  results <- vector("list", n_mc * n_scenarios)
  result_idx <- 0
  
  ## 蒙特卡洛循环
  for (i in seq_len(n_mc)) {
    
    mc_row <- mc_pat[i, ]
    
    eta_CL_true <- mc_row$eta_CL_true
    eta_V_true <- if (is.na(mc_row$eta_V_true)) 0 else mc_row$eta_V_true
    eta_CL_ind <- mc_row$eta_CL_ind
    eta_V_ind <- if (is.na(mc_row$eta_V_ind)) 0 else mc_row$eta_V_ind
    
    ## ========== Step 1: 构建条件概率库（IND 和 POP）==========
    
    ## IND 模式
    cond_lib_ind <- build_cond_prob_library(
      eta_CL = eta_CL_ind, eta_V = eta_V_ind,
      model_name = model_name, model_bank = model_bank,
      iCov = iCov, dose_mg = dose_mg, ii_h = ii_h, n_dose = n_dose,
      t_sample_h = t_sample_h, scenario_lib = scenario_library,
      n_bins = n_bins, M_sim = M_SIM, DELTA = DELTA, C_MIN = C_MIN
    )
    
    ## POP 模式
    cond_lib_pop <- build_cond_prob_library(
      eta_CL = 0, eta_V = 0,
      model_name = model_name, model_bank = model_bank,
      iCov = iCov, dose_mg = dose_mg, ii_h = ii_h, n_dose = n_dose,
      t_sample_h = t_sample_h, scenario_lib = scenario_library,
      n_bins = n_bins, M_sim = M_SIM, DELTA = DELTA, C_MIN = C_MIN
    )
    
    ## ========== Step 2: 对每个真实场景进行测试 ==========
    
    for (s_true_idx in seq_len(n_scenarios)) {
      
      s_true <- scenario_library$scenario_id[s_true_idx]
      
      ## 构建真实场景的给药事件表
      ev_true <- make_ev_for_scenario(
        dose_mg = dose_mg, ii_h = ii_h, n_dose = n_dose,
        t_sample_h = t_sample_h,
        dose_13_mult = scenario_library$dose_13_mult[s_true_idx],
        dose_14_mult = scenario_library$dose_14_mult[s_true_idx]
      )
      
      ## 计算真实预测浓度（使用 θ_true）
      C_true_test <- calc_C_pred(
        model_name = model_name,
        eta_CL = eta_CL_true,
        eta_V = eta_V_true,
        iCov = iCov,
        ev_df = ev_true,
        t_sample_h = t_sample_h
      )
      
      ## 抽取残差，生成观测浓度
      C_obs_test <- generate_residual_samples(
        C_pred = C_true_test,
        sigma_type = model_bank$sigma_type,
        sigma_prop_sd = model_bank$sigma_prop_sd,
        sigma_add_sd = model_bank$sigma_add_sd,
        n_samples = 1
      )
      
      ## ========== Step 3: 从条件概率库计算后验概率 ==========
      
      posterior_ind <- calc_posterior_from_library(
        C_obs = C_obs_test,
        cond_prob_matrix = cond_lib_ind$cond_prob,
        n_bins = n_bins, DELTA = DELTA, C_MIN = C_MIN
      )
      
      posterior_pop <- calc_posterior_from_library(
        C_obs = C_obs_test,
        cond_prob_matrix = cond_lib_pop$cond_prob,
        n_bins = n_bins, DELTA = DELTA, C_MIN = C_MIN
      )
      
      ## ========== Step 4: 存储结果 ==========
      
      for (s_j_idx in seq_len(n_scenarios)) {
        
        result_idx <- result_idx + 1
        
        results[[result_idx]] <- tibble(
          patient_id = patient_id_target,
          pop = pat_info$pop[1],
          model_name = model_name,
          mc_iter = i,
          s_true = s_true,
          s_j = scenario_library$scenario_id[s_j_idx],
          
          C_obs_test = C_obs_test,
          C_true_test = C_true_test,
          C_pred_ind = cond_lib_ind$C_pred[s_j_idx],
          C_pred_pop = cond_lib_pop$C_pred[s_j_idx],
          
          posterior_ind = posterior_ind[s_j_idx],
          posterior_pop = posterior_pop[s_j_idx],
          
          eta_CL_true = eta_CL_true,
          eta_V_true = eta_V_true,
          eta_CL_ind = eta_CL_ind,
          eta_V_ind = eta_V_ind
        )
      }
    }
  }
  
  bind_rows(results)
}

message("✅ 单患者处理函数已定义")

## ############################################################################
## 4.6 执行后验概率计算（并行）
## ############################################################################

message("\n", paste(rep("-", 70), collapse = ""))
message("4.6 执行后验概率计算（并行）")
message(paste(rep("-", 70), collapse = ""))

t_sample_h <- tdm_rule$t_sample_h
message("采样时间: t = ", t_sample_h, " h")

## 合并患者信息表
patients_all <- bind_rows(patients_main, patients_side)

## 获取所有患者 ID
all_patient_ids <- unique(mc_all$patient_id)
message("患者数: ", length(all_patient_ids))

## 创建并行集群
message("\n🚀 启动并行计算...")
message("   核心数: ", N_CORES)

start_time <- Sys.time()

## 使用 parallel 包进行并行计算
cl <- makeCluster(N_CORES)

## 导出必要的对象和函数到集群
clusterExport(cl, c(
  "mc_all", "patients_all", "omega_sigma_bank", "scenario_library",
  "t_sample_h", "n_bins", "M_SIM", "DELTA", "C_MIN", "MASTER_SEED",
  "make_ev_for_scenario", "extract_conc_at_time", "generate_residual_samples",
  "calc_C_pred", "build_cond_prob_library", "calc_posterior_from_library",
  "process_single_patient",
  "maharaj2021_olz_model", "zang2021_olz_model", "li2018_olz_model_drug2",
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
posterior_results <- parLapply(cl, all_patient_ids, function(pid) {
  tryCatch({
    process_single_patient(
      patient_id_target = pid,
      mc_all = mc_all,
      patients_all = patients_all,
      omega_sigma_bank = omega_sigma_bank,
      scenario_library = scenario_library,
      t_sample_h = t_sample_h,
      n_bins = n_bins,
      M_SIM = M_SIM,
      DELTA = DELTA,
      C_MIN = C_MIN,
      MASTER_SEED = MASTER_SEED
    )
  }, error = function(e) {
    message("Error processing patient ", pid, ": ", e$message)
    NULL
  })
})

## 关闭集群
stopCluster(cl)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

message("\n✅ 并行计算完成")
message("   耗时: ", round(as.numeric(elapsed), 2), " 分钟")

## 合并所有结果
posterior_raw <- bind_rows(posterior_results)

message("   总行数: ", nrow(posterior_raw))

## ############################################################################
## 4.7 生成汇总表
## ############################################################################

message("\n--- 4.7 生成汇总表 ---")

posterior_summary <- posterior_raw %>%
  group_by(patient_id, pop, model_name, mc_iter, s_true) %>%
  summarise(
    s_pred_ind = s_j[which.max(posterior_ind)],
    posterior_max_ind = max(posterior_ind),
    posterior_true_ind = posterior_ind[s_j == s_true],
    
    s_pred_pop = s_j[which.max(posterior_pop)],
    posterior_max_pop = max(posterior_pop),
    posterior_true_pop = posterior_pop[s_j == s_true],
    
    C_obs_test = first(C_obs_test),
    C_true_test = first(C_true_test),
    
    .groups = "drop"
  ) %>%
  mutate(
    match_ind = (s_pred_ind == s_true),
    match_pop = (s_pred_pop == s_true)
  )

cat("\n=== 汇总表预览 ===\n")
print(head(posterior_summary, 20))

message("\n✅ 汇总表生成完成")
message("   行数: ", nrow(posterior_summary))

## ############################################################################
## 4.8 保存输出文件
## ############################################################################

message("\n--- 4.8 保存输出文件 ---")

saveRDS(posterior_raw, file.path(OUT_DIR_S4, "posterior_raw.rds"))
message("  ✓ posterior_raw.rds 已保存（", nrow(posterior_raw), " 行）")

saveRDS(posterior_summary, file.path(OUT_DIR_S4, "posterior_summary.rds"))
write_csv(posterior_summary, file.path(OUT_DIR_S4, "posterior_summary.csv"))
message("  ✓ posterior_summary.rds / .csv 已保存（", nrow(posterior_summary), " 行）")

saveRDS(scenario_library, file.path(OUT_DIR_S4, "scenario_library.rds"))
message("  ✓ scenario_library.rds 已保存")

## 保存核心参数
params <- list(
  DELTA = DELTA,
  M_SIM = M_SIM,
  C_MIN = C_MIN,
  C_MAX = C_MAX,
  n_bins = n_bins,
  N_CORES = N_CORES,
  elapsed_minutes = as.numeric(elapsed)
)
saveRDS(params, file.path(OUT_DIR_S4, "params.rds"))
message("  ✓ params.rds 已保存")

## ############################################################################
## 4.9 诊断统计
## ############################################################################

message("\n--- 4.9 诊断统计 ---")

accuracy_by_patient <- posterior_summary %>%
  group_by(patient_id, model_name) %>%
  summarise(
    n_tests = n(),
    n_correct_ind = sum(match_ind),
    n_correct_pop = sum(match_pop),
    accuracy_ind = mean(match_ind) * 100,
    accuracy_pop = mean(match_pop) * 100,
    accuracy_diff = accuracy_ind - accuracy_pop,
    .groups = "drop"
  )

cat("\n=== 各患者准确率 ===\n")
print(accuracy_by_patient)

accuracy_by_scenario <- posterior_summary %>%
  group_by(s_true) %>%
  summarise(
    n_tests = n(),
    accuracy_ind = mean(match_ind) * 100,
    accuracy_pop = mean(match_pop) * 100,
    accuracy_diff = accuracy_ind - accuracy_pop,
    .groups = "drop"
  )

cat("\n=== 各场景准确率 ===\n")
print(accuracy_by_scenario)

overall_accuracy <- posterior_summary %>%
  summarise(
    n_total = n(),
    accuracy_ind = mean(match_ind) * 100,
    accuracy_pop = mean(match_pop) * 100,
    accuracy_diff = accuracy_ind - accuracy_pop
  )

cat("\n=== 总体准确率 ===\n")
print(overall_accuracy)

accuracy_stats <- list(
  by_patient = accuracy_by_patient,
  by_scenario = accuracy_by_scenario,
  overall = overall_accuracy
)
saveRDS(accuracy_stats, file.path(OUT_DIR_S4, "accuracy_stats.rds"))
message("\n  ✓ accuracy_stats.rds 已保存")

## ############################################################################
## 4.10 Stage 4 完成
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("✅ Stage 4 全部完成（第八版 - 分箱蒙特卡洛法）")
message(paste(rep("=", 70), collapse = ""))

message("\n📁 输出目录: ", OUT_DIR_S4)
message("\n📄 输出文件：")
message("  - posterior_raw.rds          （原始后验概率，", nrow(posterior_raw), " 行）")
message("  - posterior_summary.rds/.csv （汇总表，", nrow(posterior_summary), " 行）")
message("  - scenario_library.rds       （场景库定义）")
message("  - params.rds                 （核心参数）")
message("  - accuracy_stats.rds         （准确率统计）")

message("\n📊 核心参数回顾：")
message("   分箱步长 Δ = ", DELTA, " ng/mL")
message("   每场景模拟次数 M = ", M_SIM)
message("   并行核心数 = ", N_CORES)
message("   总耗时 = ", round(as.numeric(elapsed), 2), " 分钟")

message("\n📊 总体准确率：")
message("   方法 A (Ind): ", sprintf("%.2f%%", overall_accuracy$accuracy_ind))
message("   方法 B (Pop): ", sprintf("%.2f%%", overall_accuracy$accuracy_pop))
message("   提升幅度:     ", sprintf("%+.2f%%", overall_accuracy$accuracy_diff))

message("\n🔗 后续脚本使用方法：")
message('  posterior <- readRDS("outputs8/outputs8_stage4/posterior_summary.rds")')
message('  accuracy <- readRDS("outputs8/outputs8_stage4/accuracy_stats.rds")')

message("\n⏭️  下一步：运行 [第八版代码-主线5.R]（依从性判别评估）")
## ############################################################################
## 第九版代码-敏感性分析-主控-修复版.R
## 
## 修复内容：
##   1. 每个真实模型独立运行，一个失败不影响其他
##   2. 支持从指定模块开始（跳过已完成的模块）
##   3. 在汇总前检查所有结果是否存在
## 
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("第九版 敏感性分析：跨真实模型验证（修复版）")
message(paste(rep("=", 70), collapse = ""))

## ############################################################################
## 0. 路径与全局设置
## ############################################################################

## 路径设置
PROJ_ROOT <- "D:/Users/YujiaZhang/Desktop/26救赎之道/25.9.29 毕业设计/5正式实施！/多成人模型/体现我卓越的项目管理能力"
CODE_DIR <- file.path(PROJ_ROOT, "代码这边请/正式代码/第九版代码")

## 敏感性分析输出根目录
OUT_ROOT_SA <- file.path(PROJ_ROOT, "outputs9_SA")

## 创建输出目录
if (!dir.exists(OUT_ROOT_SA)) {
  dir.create(OUT_ROOT_SA, recursive = TRUE)
}

## ============================================================================
## 配置：要运行的模型和起始模块
## ============================================================================

## 定义敏感性分析的真实模型
SENSITIVITY_MODELS <- c("li2018", "yin2016", "sun2021")

## 每个模型的起始模块（1=从头开始，3=从模块3开始，5=跳过）
## 根据已完成的情况调整
START_MODULES <- c(
  li2018 = 1,    # li2018 从头开始
  yin2016 = 1,   # yin2016 从头开始
  sun2021 = 1    # sun2021 从头开始
)

message("\n配置信息：")
for (model in SENSITIVITY_MODELS) {
  message("  ", model, ": 从模块 ", START_MODULES[model], " 开始")
}

## ############################################################################
## 1. 循环执行敏感性分析
## ############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("开始敏感性分析循环")
message(paste(rep("=", 70), collapse = ""))

for (SA_TRUE_MODEL in SENSITIVITY_MODELS) {
  
  message("\n\n", paste(rep("#", 70), collapse = ""))
  message("### 敏感性分析：真实模型 = ", SA_TRUE_MODEL, " ###")
  message(paste(rep("#", 70), collapse = ""))
  
  ## 设置该敏感性分析的输出目录
  OUT_DIR_SA <- file.path(OUT_ROOT_SA, SA_TRUE_MODEL)
  
  if (!dir.exists(OUT_DIR_SA)) {
    dir.create(OUT_DIR_SA, recursive = TRUE)
    message("已创建目录: ", OUT_DIR_SA)
  }
  
  ## 获取该模型的起始模块
  start_module <- START_MODULES[SA_TRUE_MODEL]
  
  ## 使用 tryCatch 包裹，一个模型失败不影响其他
  tryCatch({
    
    ## ========== 模块 1：数据生成 ==========
    if (start_module <= 1) {
      message("\n>>> 运行模块 1：数据生成...")
      source(file.path(CODE_DIR, "第九版代码-敏感性分析-模块1.R"), encoding = "UTF-8")
    } else {
      message("\n>>> 跳过模块 1（已完成）")
    }
    
    ## ========== 模块 2：跨模型 MAPB ==========
    if (start_module <= 2) {
      message("\n>>> 运行模块 2：跨模型 MAPB...")
      source(file.path(CODE_DIR, "第九版代码-敏感性分析-模块2.R"), encoding = "UTF-8")
    } else {
      message("\n>>> 跳过模块 2（已完成）")
    }
    
    ## ========== 模块 3：跨模型判别 ==========
    if (start_module <= 3) {
      message("\n>>> 运行模块 3：跨模型判别...")
      source(file.path(CODE_DIR, "第九版代码-敏感性分析-模块3.R"), encoding = "UTF-8")
    } else {
      message("\n>>> 跳过模块 3（已完成）")
    }
    
        ## ========== 模块 4：一致性评估 ==========
    if (start_module <= 4) {
      message("\n>>> 运行模块 4：一致性评估...")
      source(file.path(CODE_DIR, "第九版代码-敏感性分析-模块4.R"), encoding = "UTF-8")
      
      ## ========== 模块 4E：CL 一致性分析（新增）==========
      message("\n>>> 运行模块 4E：CL 一致性分析...")
      source(file.path(CODE_DIR, "第九版代码-敏感性分析-模块4E-CL一致性.R"), encoding = "UTF-8")
    } else {
      message("\n>>> 跳过模块 4（已完成）")
    }
    
    message("\n✅ 敏感性分析完成：", SA_TRUE_MODEL)
    
  }, error = function(e) {
    message("\n❌ 敏感性分析失败：", SA_TRUE_MODEL)
    message("   错误信息：", e$message)
  })
}

## ############################################################################
## 2. 检查结果并汇总
## ############################################################################

message("\n\n", paste(rep("=", 70), collapse = ""))
message("检查敏感性分析结果")
message(paste(rep("=", 70), collapse = ""))

## 检查每个模型的结果是否存在
results_exist <- sapply(SENSITIVITY_MODELS, function(model) {
  file.exists(file.path(OUT_ROOT_SA, model, "evaluation_results.rds"))
})

message("\n结果检查：")
for (model in SENSITIVITY_MODELS) {
  status <- if (results_exist[model]) "✅ 完成" else "❌ 未完成"
  message("  ", model, ": ", status)
}

## 只有当所有结果都存在时才运行汇总
if (all(results_exist)) {
  message("\n所有敏感性分析已完成，开始汇总...")
  source(file.path(CODE_DIR, "第九版代码-敏感性分析-汇总-独立版.R"), encoding = "UTF-8")
} else {
  message("\n⚠️ 部分敏感性分析未完成，跳过汇总")
  message("   请检查失败的模型并重新运行")
}

message("\n", paste(rep("=", 70), collapse = ""))
message("敏感性分析主控脚本执行完毕")
message(paste(rep("=", 70), collapse = ""))

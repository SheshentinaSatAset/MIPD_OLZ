# ==============================================================================
# 奥氮平PopPK模型（主线）：Maharaj et al. 2021（儿童模型）
# 文件路径：D:/Users/YujiaZhang/Desktop/26救赎之道/25.9.29 毕业设计/5正式实施！/多成人模型/体现我卓越的项目管理能力/代码这边请/poppk模型/maharaj2021_olz_model.R
# 模型类型：一室模型+一级吸收+线性消除+组合残差模型
# 建模工具：rxode2 + nlmixr2（兼容NONMEM FOCE-INTERACTION）
# 命名规范：语义化参数名（B_XXX为系数，TVXXX为典型值，TM50/HILL为成熟模型参数）
# ==============================================================================

# 1. 加载依赖包（工程级检查）
if (!requireNamespace("rxode2", quietly = TRUE)) {
  stop("请安装rxode2包：install.packages('rxode2')")
}
if (!requireNamespace("nlmixr2", quietly = TRUE)) {
  stop("请安装nlmixr2包：install.packages('nlmixr2')")
}
library(rxode2)
library(nlmixr2)

# 2. 定义Maharaj2021儿童PopPK模型（rxode2格式）
maharaj2021_olz_model <- rxode2({
  ## -------------------------------
  ## 群体典型值（文献表3 Final model估计值）
  ## -------------------------------
  KA        <- 0.758    # 吸收速率常数 (h⁻¹)，固定值（来自FDA成人模型）
  TVCL      <- 16.8     # CL/F典型值 (L/h，70kg成熟个体)
  TVV       <- 663      # V/F典型值 (L，70kg成熟个体)
  TM50      <- 70       # 成熟半衰期（PMA，周）：Sigmoid模型中点
  HILL      <- 3.97     # 希尔系数：描述PMA对CL/F的成熟速率
  B_WT_EXP  <- 0.486    # 体重对CL/F的异速生长指数（估计值）
  
  ## -------------------------------
  ## 个体参数（含协变量效应+个体间变异）
  ## 协变量说明：WT=体重(kg)，PMA=出生后年龄(周)
  ## -------------------------------
  # 表观清除率（CL/F）：体重异速生长 + PMA成熟模型 + 个体间变异
  CL <- TVCL * 
    (WT / 70)^B_WT_EXP *  # 体重标度（70kg为参考）
    (PMA^HILL) / (TM50^HILL + PMA^HILL) *  # PMA成熟Sigmoid模型
    exp(eta_CL)           # CL/F个体间变异（指数模型）
  
  # 表观分布容积（V/F）：体重线性标度（指数=1）
  V <- TVV * (WT / 70)^1
  
  ## -------------------------------
  ## 房室微分方程（一级吸收+线性消除）
  ## -------------------------------
  d/dt(depot) <- -KA * depot  # 吸收室
  d/dt(cent)  <-  KA * depot - (CL / V) * cent  # 中央室
  
  ## -------------------------------
  ## 浓度转换（模型内mg/L → 对外ng/mL）
  ## -------------------------------
  C_mgL <- cent / V       # 模型内部浓度单位 (mg/L)
  C_ngmL <- C_mgL * 1000  # 对外输出单位 (ng/mL，1 mg/L = 1000 ng/mL)
})

# 3. 定义个体间变异矩阵（OMEGA，仅CL/F有IIV，文献表3）
maharaj2021_omega <- lotri({
  # IIV-CL/F：ω²=0.309（文献直接给出），严格对应log((CV/100)^2 + 1)推导结果
  eta_CL ~ 0.309
})

# 4. 定义残差模型（组合误差：比例+加性，文献表3）
maharaj2021_sigma <- lotri({
  prop_err ~ 0.0772  # 比例误差方差（文献ε1²=0.0772）
  add_err  ~ 1.51    # 加性误差方差（文献ε2²=1.51）
})
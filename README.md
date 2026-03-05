# MIPD_OLZ

**Model-Informed Precision Dosing for Olanzapine | 奥氮平模型引导精准用药**

## 项目简介 | Overview

本项目聚焦于**奥氮平（Olanzapine）**的模型引导精准用药（MIPD）研究。奥氮平是临床常用的非典型抗精神病药，因其显著的个体间药动学变异性，是开展 MIPD 研究的理想对象。

This project focuses on Model-Informed Precision Dosing (MIPD) for **Olanzapine**, a widely used atypical antipsychotic with substantial pharmacokinetic variability, making it an ideal candidate for individualized dosing strategies.

## 研究方向 | Research Scope

- 奥氮平群体药动学（PPK）建模
- 协变量筛选（吸烟状态、性别、体重、CYP1A2/2D6 基因型等）
- 贝叶斯自适应给药方案设计
- 治疗药物监测（TDM）数据整合
- 临床实践中的剂量个体化工具开发

## AI 协作配置 | Copilot Configuration

本仓库已配置 GitHub Copilot 专属指令（`.github/copilot-instructions.md`），将 Copilot 设置为具备定量药理学专业背景的研究协作者，能够：

- 提供 PK/PD 建模的专业技术支持
- 基于客观评估审视研究方案（而非一味认同）
- 支持中英文学术写作
- 提供数学推导和高质量代码实现

## 技术栈 | Tools & Methods

- **建模软件**: NONMEM, Monolix, nlmixr2
- **数据分析**: R, Python
- **MIPD 方法**: 贝叶斯估算（MAP-Bayesian estimation）
- **统计诊断**: VPC, GOF plots, bootstrap

## 参考框架 | Reference Framework

- FDA Model-Informed Drug Development (MIDD) Guidance
- EMA Modeling and Simulation Qualification Guidelines
- ISOP/PAGE population PK modeling best practices

# PhD Study3: SFC分析

## 数据准备

### 拷贝需要的静息态数据
fsaverage5:`copy_CCNP_fsaverage5.ps1`和`unzip_CCNP.ps1`
fslr10k:数据来自c罗

### 亚型预测
`predict_clusters_abide.py`使用asd研究中的分类器（）预测abide中剩余男性asd（大于等于13岁）的亚型

### 被试汇总
把原来训练得到的分型，和现在预测的分型，和典型对照所有人汇总到一起
`sum_predict_all_abide.R`
/Volumes/Zuolab_XRF/output/abide/ABIDE_cluster_all_subjects.csv

## ROI时间序列提取

### 每个网络选取1个ROI
`select_ROI_fsLR32k.m`
/Users/xuerufan/Downloads/PhD论文材料【师大电脑没有备份】/研究三/ROI/fslr32k/DU15NET_connected_core_ROI_visual.dlabel.nii
/Users/xuerufan/PhD-Study3/output/DU15NET_connected_core_ROI_report.csv

### 空间转换网络ROI
将定义在 **fsLR-32k** 表面空间中的ROIs，转换到 **FreeSurfer 的 fsaverage5** 表面空间，从而与基于 fsaverage5 预处理的功能数据在**顶点层面实现精确对齐**。

`export_ROI_wb.m`  
**fsLR-32k ROI → Workbench metric** 将以顶点索引形式存储的 ROI（fsLR-32k）导出为二值表面 metric 文件（`.func.gii`），作为后续 Workbench 重采样的输入。
output/ROI_metrics_fsLR32k/
lh.DU15network<N>_fsLR32k.func.gii
rh.DU15network<N>_fsLR32k.func.gii

`resample_fsLR32k_to_fsLR10k.sh` 将定义在 **fsLR-32k** 表面空间中的感兴趣区（regions of interest, ROIs），转换到 **fsLR-10k**表面空间。
**fsLR-32k → fsLR-10k（Workbench）**
output/ROI_fslr10k/
lh.DU15Net<N>_fsLR10k.func.gii
rh.DU15Net<N>_fsLR10k.func.gii

`resample_fsLR32k_to_fsaverage164k.sh`  
**fsLR-32k → fsaverage-164k（Workbench）** 将 fsLR-32k ROI metric 重采样至 fsaverage 高分辨率表面（约 164k 顶点/半球）。
output/ROI_metrics_fsaverage164/
lh.DU15network<N>_fsaverage164k.func.gii
rh.DU15network<N>_fsaverage164k.func.gii

`resample_fsaverage164_to_fsaverage5.sh`  
**fsaverage-164k → fsaverage5（FreeSurfer）** 将 fsaverage-164k 的 ROI metric 转换为 FreeSurfer overlay，并映射到 fsaverage5 表面。
output/ROI_fsaverage5/
lh.DU15Net<N>_fsaverage5.mgh
rh.DU15Net<N>_fsaverage5.mgh

`convert_DU15ROI_to_nifti.sh`  
**fsaverage5(.mgh) → fsaverage5（nii.gz）** 格式转换，在windows电脑上（不使用hfreesurfer时）读数据需要。
output/ROI_fsaverage5/
lh.DU15Net<N>_fsaverage5.nii.gz
rh.DU15Net<N>_fsaverage5.nii.gz

### 可视化ROI
`plot_Du15ROI_fsaverage5.py`
在 **fsaverage5 皮层表面空间** 上可视化 **DU15 功能网络（15 个 network）**。将左右半球的 DU15 ROI（`.mgz`）合并为网络标签图，按照DU定义的 RGB 颜色进行渲染，并在膨胀皮层表面上生成多个视角（外侧、内侧、前、后、上、下）的高分辨率图像。
visual/
DU15_fsaverage5_LH_anterior.png
...

`build_DU15Net_fsLR10k_dlabel.sh`
合并多个 ROI 为单一整数 surface metric
output/ROI_fslr10k/
lh.DU15Net_fsLR10k.func.gii
rh.DU15Net_fsLR10k.func.gii
由 surface metric 生成 1个包含所有ROI的CIFTI dense scalar（注意：这一步没有修改每个ROI的颜色和Du15一致，也没有输出ROI得截图，因为已经有fslr32k和fsaverage5两个空间的了，只是分辨率不一样，需要的话再用matlab搞一下吧）
output/ROI_fslr10k/
DU15Net_fsLR10k.dscalar.nii

### 提取时间序列
提取数据中ROI体素的时间序列，然后进行逐个时间序列的去均值和去线性漂移，两个半球的ROI内进行体素平均。

**ABIDE** 
`extract_timeseries_ABIDE.mac.m`
data/abide/timeseries/
29096_DU15_network_ts.mat

**CCNP** 
fsaverage5: `extract_timeseries_CCNP_fsaverage5_mac`和`extract_timeseries_CCNP_fsaverage5_wins.m`
data/timeseries/fsaverage/SITE/SUB/
rest1_DU15_roi_ts.mat
rest2_DU15_roi_ts.mat

fslr10k: `extract_timeseries_CCNP.wins.m`
E:\PhDproject\Study3\data\timeseries\fslr10k
CCNPPEK0001_01_rest01_DU15_network_ts.mat


## 计算SFEI

### 计算汇总SFC
`calculate_SFC_nbwt_ABIDE.m`（不走回头路版）和`calculate_SFC_nbwt_CCNP.m`计算SFC、zSFC和SFEI
`sum_SFC_nbwt_embedding_ABIDE.m`和`sum_SFC_nbwt_embedding_CCNP.m`汇总所有被试的SFEI

### 筛选CCNP被试数据
`select_ccnppek_participant.R`（和`select_ccnpckg_participant_fsaverage5.R`）筛选出meanfd小于0.3的被试数据，并对同一个session内的所有通过qc的run取平均
/Volumes/Zuolab_XRF/output/ccnp/sfc/sfc_nbtw_embedding/pek_step01_fd0.3_sessionAvg.xlsx

### ？？？分subtype × step计算ABIDE的SFEI均值
`mean_SFC_nbtw_by_subtype_step_ABIDE.m`
/Volumes/Zuolab_XRF/output/abide/sfc/stat/meanzSFC
zSFC_mean_TD_step01.csv
zSFC_mean_ASD_L_step03.csv
zSFC_mean_ASD_H_step08.csv

## SFEI的组间比较和相关分析

### SFEI描述性分析
`group_analysis_SFC_nbwt_embedding_ABIDE.R`
整合了多步嵌入计算结果、被试人口学信息、站点信息以及基于机器学习的 ASD 亚型预测结果（仅包含男性被试），并生成可直接用于论文排版的高分辨率图像。

所有分析结果统一输出至以下目录：/Volumes/Zuolab_XRF/output/abide/sfc

其中，表格结果与图像结果分别存放如下。

生成最终纳入分析的被试编号文件/Volumes/Zuolab_XRF/output/abide/sfc/sfc_participant_for_analysis.csv

以及在最后的最后，生成了纳入分析被试的汇总信息，包含编号、站点、分型、年龄。  地址在/Volumes/Zuolab_XRF/output/abide/sfc/sfc_participant_summary.csv

1. 人口学描述统计（按亚型）  
   文件名：  
   sfc_demo.csv  

   内容说明：  
   各亚型被试的样本量、平均年龄及年龄标准差。

2. 功能网络层面的嵌入描述统计  
   文件名：  
   sfc_network.csv  

   内容说明：  
   不同亚型在各功能网络上的 embedding 均值、标准差、中位数及四分位距。

3. 嵌入步骤层面的描述统计  
   文件名：  
   sfc_step.csv  

   内容说明：  
   不同亚型在各嵌入步骤上的 embedding 均值与标准差。

所有图像文件统一保存在以下子目录：/Volumes/Zuolab_XRF/output/abide/sfc/plot
当前版本生成的主要图像为：sfc_heatmap_abide.png 功能网络 × 嵌入步骤的平均 embedding 热图   


### SFEI组间对比分析
`group_contrast_SFC_nbtw_embedding_ABIDE.R`
该脚本在前期 SFC embedding 计算与描述性分析的基础上，进一步对 ABIDE 数据集中不同临床分组相对于 TD 组的 embedding 差异进行系统性统计建模与可视化分析。分析整合了多步嵌入结果、被试人口学信息、扫描站点信息以及基于机器学习预测的 ASD 亚型标签，仅纳入男性被试，以避免性别混杂效应。

## 1. 组间对比统计结果（Network × Step）
/Volumes/Zuolab_XRF/output/abide/sfc/stat/difference/sfc_nbwt_group_constrast.csv
该文件包含在每一个功能网络 × 嵌入步骤组合上进行的线性模型组间对比结果。模型在控制年龄（AGE_AT_SCAN）和扫描站点（site）后，分别估计以下三类对比：
- ASD（合并组） vs TD  
- ASD-L vs TD  
- ASD-H vs TD
- L vs H
输出字段包括对比的估计值（estimate）、标准误、t 值、原始 p 值等统计量，并额外标记是否达到未经多重校正的显著性水平（p < 0.05），用于后续结果展示与探索性解释。

## 2. 功能网络 × 嵌入步骤的组间差异热图
/Volumes/Zuolab_XRF/output/abide/sfc/plot/difference/sfc_nbwt_subtypes_contrasts.png
该图以热图形式展示不同功能网络和嵌入步骤上 embedding 的组间差异估计值（ΔEmbedding），按功能系统组织顺序排列网络维度，并在三个并列面板中分别呈现：
- ASD vs TD  
- ASD-L vs TD  
- ASD-H vs TD
- L vs H
颜色表示回归模型中组别效应的方向与大小，黑色空心圆标记表示在对应 Network × Step 位置达到原始显著性水平（p < 0.05）的对比结果。图像采用透明背景与高分辨率设置，可直接用于博士论文结果章节的核心图像展示。

## 3. 方法学要点说明
- 所有组间对比均基于线性模型（lm），并通过 `emmeans` 计算处理对照（TD）参照的组间对比。  
- ASD 合并组分析与亚型分析在同一数据框架下独立运行，确保统计结果的可比性。  
- 嵌入步骤（Step）与功能网络（Network）均作为分组单元逐一建模，不引入跨步骤或跨网络的平均化假设。  
- 当前结果未进行多重比较校正，显著性标记仅用于结果结构与空间分布的探索性展示，相关解释需在论文中谨慎限定。


### SFEI与认知行为的相关分析
`statistic_SFC_correlations_abide.R`
本脚本用于分析ABIDE数据集中不同分组/分型的SFEI与认知行为量表之间的相关关系。在控制站点效应的前提下，分别对嵌入值与行为指标进行残差相关分析，整体分析流程与此前基于脑形态学百分位数的相关分析保持方法学一致性

本脚本依赖以下数据文件：
1. **SFC embedding 结果（多步）**  `/Volumes/Zuolab_XRF/output/abide/sfc/sfc_embedding/`
   该目录包含 step01–step08 共 8 个 Excel 文件，每个文件对应一步 SFC 嵌入结果。行表示被试，列表示不同功能网络（Net01–Net15）的 embedding 值。

2. **最终纳入分析的被试列表**  `/Volumes/Zuolab_XRF/output/abide/sfc/sfc_participant_for_analysis.csv` 
   本研究中最终用于 SFC 分析的男性被试编号，用于在不同数据源之间进行统一筛选与对齐。

3. **ASD 亚型预测结果（仅男性）**   `/Volumes/Zuolab_XRF/output/abide/abide_cluster_predictions_male.csv`  

4. **认知与临床行为数据**     `/Users/xuerufan/DCM-Project-PhD-Study3-/supplement/abide_A_all_240315.csv`  
   包含 FIQ、ADOS、ADI-R 以及 SRS 等认知与临床行为量表数据，同时包含站点信息（SITE_ID）。

### 分析方法概述
- 自变量  
  不同嵌入步骤（step01–step08）中，各功能网络对应的 SFC embedding 值。
- 因变量  
  认知与临床行为指标，包括 FIQ、ADOS、ADI-R 及 SRS 各分量表。其中 ADOS_2_RRB 作为非正态分布变量单独处理。
- 统计方法  
  对 ADOS_2_RRB 使用 Spearman 相关，对其余连续行为指标使用 Pearson 相关。在相关分析前，分别对自变量与因变量回归掉站点效应（SITE_ID），并基于残差计算相关系数。每一对变量要求有效样本量不少于 30。所有相关结果在每个嵌入步骤及亚型内进行 FDR 多重比较校正。

### 输出结果说明

所有分析结果统一输出至以下目录：`/Volumes/Zuolab_XRF/output/abide/sfc/stat/corr/`

生成的主要结果文件为：

1. **SFC embedding × 认知行为相关结果汇总表（按亚型）**  
   文件名：  
   `/Volumes/Zuolab_XRF/output/abide/sfc/stat/corr/SFEI_nbwt_cognition_LH.csv`
   内容说明：  
   该表格汇总了所有嵌入步骤（step01–step08）、所有功能网络以及所有认知与临床行为指标在 ASD-L 与 ASD-H 亚型中的相关分析结果。表格包含以下字段：  
   - step：嵌入步骤  
   - cluster：亚型（L / H）  
   - name_brain：功能网络  
   - name_cog：认知或行为指标  
   - coef：相关系数  
   - p_value：原始 p 值  
   - P_adj：FDR 校正后的 p 值  
   - df：残差模型自由度
   - 
### 显著相关结果的网络重命名与层级整理
在完成 SFC embedding 与认知行为的相关分析后，本脚本进一步对结果数据进行后处理，以便直接用于论文结果展示与可视化分析。
具体而言，首先将分析结果中以 Net01–Net15 表示的功能网络，映射为对应的功能系统名称缩写（如 VIS-P、DN-A、FPN-B 等），网络定义基于既定的功能网络分区方案。随后，仅保统计显著的结果（P < 0.05）。

在整理后的结果表中，所有显著结果按照以下层级顺序进行排列：  
1）嵌入步骤（step01–step08），  
2）功能网络（按系统名称排序），  
3）ASD 亚型（ASD-L，ASD-H）。

后处理步骤生成的最终显著结果文件为：`/Volumes/Zuolab_XRF/output/abide/sfc/stat/corr/SFEI_nbwt_cognition_LH_significant.csv`

该文件包含以下信息：  
- step：SFC 嵌入步骤  
- cluster：ASD 亚型（L / H）  
- network_abbr：功能网络名称缩写  
- name_cog：认知或临床行为指标  
- coef：相关系数  
- p_value：原始 p 值  
- P_adj：FDR 校正后的 p 值  
- df：残差模型自由度  

该结果表为论文结果部分和补充材料中网络层级与嵌入步骤层级分析的核心输入，可直接用于生成分亚型的相关模式图或热图展示。

### SFC embedding × 认知行为的个体层面可视化（自动生成散点图）

在完成 SFC embedding 与认知行为指标的相关分析后，本脚本基于统计显著结果，自动生成对应的个体层面散点图，用于直观展示功能连接嵌入与行为指标之间的关系模式。该流程严格以相关分析结果为输入，仅对达到统计显著的组合进行作图，避免人工筛选带来的主观偏差，确保结果展示的系统性与可复现性。

对于每一个达到显著性的分析结果，脚本自动遍历以下维度组合：
- SFC 嵌入步骤（step01–step08）
- 功能网络（Net01–Net15）
- 认知或临床行为指标（如 FIQ、SRS、ADOS 等）

针对每一个 *step × network × behavior* 组合，生成一张独立的散点图，具体作图策略如下：
- 每个散点代表一名被试
- 横轴为对应功能网络的 SFC embedding 值（残差）
- 纵轴为对应认知或行为指标（残差）
- 在作图前对所有变量统一控制扫描站点效应（SITE_ID）
- 使用线性回归趋势线（lm）展示组内关系
- 同时展示三组被试：
  - TD（灰色）
  - ASD-L（绿色）
  - ASD-H（红色）
该可视化方式与前述相关分析在数据清洗、协变量控制及统计模型上保持一致，确保图像展示与统计检验具有严格的方法学对应关系。
所有图像文件统一保存至以下目录：/Volumes/Zuolab_XRF/output/abide/sfc/plot/corr


## SFEI的个体偏离分析

### 汇总数据
`sum_SFEI_for_normative.R` 汇总CCNP中meanfd小于0.3，ccnp和abide中的td 筛选年龄不超过18的男性，数据整合起来用于建模
/Volumes/Zuolab_XRF/output/normative/SFEI_normative_demographic.xlsx

### combat去除站点效应
`combat_SFEI_for_combat.R` ABIDE（ASD、TD）和CCNP一起combat
/Volumes/Zuolab_XRF/output/normative/SFEI_normative_data_combat.xlsx

### 建立SFEI的常模
`statistic_FD_correlations_ccnp.R`分析fd和SFEI是否相关，结果是不相关
/Volumes/Zuolab_XRF/output/normative/SFEI_fd_correlation.xlsx

`gamlss_SFEI_ccnp.R` 用CCNP和ABIDE中的TD一起建模GAMLSS

`centile_ASD_SFEI.R`计算 ASD的centile

### 分析ASD整体及亚型的z偏离与网络组织模式
`z_analysis_SFEI_step_ABIDE.R`
`z_analysis_SFEI_network_ABIDE.R`

### centile和认知行为的相关
`statistic_SFC_centile_correlations_abide.R`



# 研究四 rDCM分析

## Step : rDCM 批处理

基于 **TAPAS toolbox**的有效连接估计，所有脚本均采用 **classic rDCM（无 sparsity 约束）** 和 **有sparsity 约束**两种方式，并显式指定采样间隔（TR）。


### 1️⃣ `calculate_rDCM_ABIDE.m`

- 每个被试一个时间序列文件，根据 `siteName` 自动匹配 TR


### 2️⃣ `calculate_rDCM_CCNP.m`

- 每个被试可包含多个 run（如 rest1 / rest2），CKG：2.5 s  ，PEK：2.0 s

### 输出结构示例

rDCM/
├── ABIDE/
│   └── sub-XXXX_rDCM.mat
├── CKG/
│   └── sub-XXXX_ses01/
│       └── sub-XXXX_rest1_rDCM.mat
└── PEK/
    └── sub-XXXX_ses01/
        └── sub-XXXX_rest1_rDCM.mat

## Step : 汇总rDCM结果
`sum_rDCM_ABIDE.m`
`sum_rDCM_CCNP.m`
CCNP_rDCM_summary.xlsx
ABIDE_rDCM_summary.xlsx

`sum_srDCM_ABIDE.m`
ABIDE_rDCM_summary.xlsx


## Step : rDCM组水平分析
`group_analysis_rDCM_ABIDE.R`
/Volumes/Zuolab_XRF/data/abide/stats/ABIDE_rDCM_group_effects.csv
/Volumes/Zuolab_XRF/data/abide/stats/ABIDE_rDCM_group_effects_male_under13.csv

`group_analysis_srDCM_ABIDE.R`
/Volumes/Zuolab_XRF/data/abide/stats/ABIDE_srDCM_density.csv
/Volumes/Zuolab_XRF/data/abide/stats/ABIDE_srDCM_density_male_under13.csv

## Step : 可视化组水平分析结果
`gplot_roup_analysis_rDCM_ABIDE.R`
/Volumes/Zuolab_XRF/data/abide/figures/ABIDE_rDCM_group_effects.png

`group_analysis_srDCM_ABIDE.R`
/Volumes/Zuolab_XRF/data/abide/figures/ABIDE_srDCM_logistic.png

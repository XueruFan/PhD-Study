# PhD Study - SFEI分析

## 数据准备

### 拷贝需要的静息态数据
fsaverage5:`copy_CCNP_fsaverage5.ps1`和`unzip_CCNP.ps1`
fslr10k:数据来自c罗

### 亚型预测
`predict_clusters_abide.py`使用asd研究中的分类器预测abide中剩余男性asd（大于等于13岁）的亚型
/Volumes/Zuolab_XRF/output/abide/abide_cluster_predictions_male.csv

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
`calculate_zSFEI_ABIDE.m`（不走回头路版）和`calculate_zSFEI_CCNP.m`计算SFC、SFEI和zSFEI
`sum_zSFEI_ABIDE.m`和`sum_zSFEI_CCNP.m`汇总所有被试的SFEI

### 筛选CCNP被试数据
`select_ccnppek_participant_zSFEI.R`（和`select_ccnpckg_participant_fsaverage5.R`）筛选出meanfd小于0.3的被试数据，并对同一个session内的所有通过qc的run取平均
/Volumes/Zuolab_XRF/output/ccnp/sfc/zsfei/pek_step01_fd0.3_sessionAvg.xlsx

### ？？？分subtype × step计算ABIDE的SFEI均值
`mean_SFC_nbtw_by_subtype_step_ABIDE.m`
/Volumes/Zuolab_XRF/output/abide/sfc/stat/meanzSFC
zSFC_mean_TD_step01.csv
zSFC_mean_ASD_L_step03.csv
zSFC_mean_ASD_H_step08.csv

## SFEI的组间比较和相关分析

### zSFEI描述性分析
`group_analysis_zSFEI_ABIDE.R`
生成最终纳入分析271人asd的被试编号文件/Volumes/Zuolab_XRF/output/abide/sfc/des/zSFEI_abide_id.csv
人口信息文件/Volumes/Zuolab_XRF/output/abide/sfc/des/zSFEI_abide_demo.csv
年龄人数的分站点统计/Volumes/Zuolab_XRF/output/abide/sfc/des/zSFEI_abide_summary.csv
功能网络层面的zSFEI描述统计 Volumes/Zuolab_XRF/output/abide/sfc/des/zSFEI_abide_network.csv
连接步数层面的zSFEI描述统计  /Volumes/Zuolab_XRF/output/abide/sfc/des/zSFEI_abide_step.csv

/Volumes/Zuolab_XRF/output/abide/sfc/plot/zSFEI_network_abide.png
/Volumes/Zuolab_XRF/output/abide/sfc/plot/zSFEI_network_step_abide.png
/Volumes/Zuolab_XRF/output/abide/sfc/plot/zSFEI_heatmap_abide.png

### zSFEI组间对比分析
`group_contrast_zSFEI_ABIDE.R`
/Volumes/Zuolab_XRF/output/abide/sfc/plot/difference/zSFEI_abide_network_step_difference.png
/Volumes/Zuolab_XRF/output/abide/sfc/stat/difference/zSFEI_abide_network_step.csv

### SFEI与认知行为的相关分析
`statistic_zSFEI_correlations_abide.R`
/Volumes/Zuolab_XRF/output/abide/sfc/stat/corr

## SFEI的个体偏离分析

### 汇总数据
`sum_zSFEI_for_combat.R` 汇总CCNP中meanfd小于0.3，ccnp和abide中的td 筛选年龄不超过18的男性，数据整合起来用于建模
/Volumes/Zuolab_XRF/output/normative/zSFEI_normative_data.xlsx
/Volumes/Zuolab_XRF/output/normative/zSFEI_normative_demo.xlsx

### combat去除站点效应
`combat_zSFEI_for_normative.R` ABIDE（ASD、TD）和CCNP一起combat
/Volumes/Zuolab_XRF/output/normative/zSFEI_normative_data_combat.xlsx

### 建立SFEI的常模
`statistic_FD_zSFEI_correlations_ccnp.R`分析fd和SFEI是否相关，结果是不相关
/Volumes/Zuolab_XRF/output/normative/zSFEI_fd_correlation.xlsx

`gamlss_zSFEI.R` 用CCNP和ABIDE中的TD一起建模GAMLSS
/Volumes/Zuolab_XRF/output/normative/gamlss

### 计算偏离
`centile_ASD_zSFEI.R`计算ASD的centile
/Volumes/Zuolab_XRF/output/normative/centile/ASD_centile_results.xlsx
`centile_plot_ASD_zSFEI.R`
/Volumes/Zuolab_XRF/output/normative/centile/plots

### 分析ASD整体及亚型的z偏离与网络组织模式
`analysis_zSFEI_network_ABIDE.R`
/Volumes/Zuolab_XRF/output/normative/centile/network_stat

### centile和认知行为的相关
`statistic_zSFEI_centile_correlations_abide.R`
/Volumes/Zuolab_XRF/output/abide/sfc/stat/corr/z_score_correlation_LH.csv
/Volumes/Zuolab_XRF/output/abide/sfc/stat/corr/z_score_correlation_LH_significant.csv



# PhD Study - rDCM分析

## 计算与汇总rDCM

### 计算
`calculate_rDCM_ABIDE.m`和`calculate_rDCM_CCNP.m`（calculate_rDCM_CCNP_fsaverage5.m），srDCM和rDCM

### 汇总
`sum_srDCM_ABIDE.m`
/Volumes/Zuolab_XRF/output/abide/dcm/sum/ABIDE_rDCM_summary.xlsx
/Volumes/Zuolab_XRF/output/abide/dcm/sum/ABIDE_srDCM_summary.xlsx

`sum_srDCM_CCNP.m`
/Volumes/Zuolab_XRF/output/ccnp/dcm/sum/CCNP_rDCM_summary.xlsx
/Volumes/Zuolab_XRF/output/ccnp/dcm/sum/CCNP_srDCM_summary.xlsx

### 筛选CCNP被试数据
`select_ccnppek_participant_srDCM.R`筛选出meanfd小于0.3的被试数据，并对同一个session内的所有通过qc的run取平均
/Volumes/Zuolab_XRF/output/ccnp/dcm/sum/pek_srdcm_fd0.3_sessionAvg.xlsx
/Volumes/Zuolab_XRF/output/ccnp/dcm/sum/pek_rdcm_fd0.3_sessionAvg.xlsx

## rDCM发育阶段常模

### 汇总数据
`sum_rDCM_for_combat.R` 汇总CCNP中meanfd小于0.3，ccnp和abide中的td 筛选年龄不超过18的男性，数据整合起来用于建模
/Volumes/Zuolab_XRF/output/norm_rDCM/rDCM_normative_data.xlsx

### combat去除站点效应
`combat_rDCM_for_normative.R` ABIDE（ASD、TD）和CCNP一起combat
/Volumes/Zuolab_XRF/output/norm_rDCM/rDCM_normative_data_combat.xlsx

### 建立rDCM的常模并计算偏离
`statistic_FD_rDCM_correlations_ccnp.R`分析fd和SFEI是否相关
/Volumes/Zuolab_XRF/output/norm_rDCM/rDCM_fd_correlation.xlsx

`gamlss_rDCM.R`
/Volumes/Zuolab_XRF/output/norm_rDCM/gamlss

## 因果调控机制变异轴分析

### 对abide的rDCM结果进行plsda分析
`analysis_rDCM_plsda_ABIDE.R`

### 分析主成分构成的关键边
`analysis_rDCM_plsda_keyNet_ABIDE.R`

### 分析EC和认知行为的相关
`statistic_rDCM_ECcorrelations_ABIDE.R`

### 筛选出来关键边的相关
`pickup_KeyNetwork_corr.R`

## 统计与可视化
把abide和ccnppek中18岁以下的rDCM原始EC、EC的偏离汇总起来，进行描述性统计和绘图
/Volumes/Zuolab_XRF/output/abide/dcm/stat/ABIDE3组人原始EC的显著差异边.xlsx
/Volumes/Zuolab_XRF/output/abide/dcm/plot

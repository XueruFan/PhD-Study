############################################################
# Filter + split by cluster + within-cluster FDR
############################################################

rm(list = ls())

library(readr)
library(readxl)
library(dplyr)
library(openxlsx)

############################################################
# 路径
############################################################

corr_file   <- "/Volumes/Zuolab_XRF/output/rDCM_PLSDA/corr/EC_cognition_correlation.csv"
select_file <- "/Volumes/Zuolab_XRF/output/rDCM_PLSDA/PLSDA_KeyEdge.xlsx"
out_file    <- "/Volumes/Zuolab_XRF/output/rDCM_PLSDA/corr/PLSDA_KeyEdge_corr.xlsx"

############################################################
# 读取数据
############################################################

corr_df   <- read_csv(corr_file, show_col_types = FALSE)
# # 去除 ADI-R 行为
# corr_df <- corr_df %>%
#   filter(!grepl("^ADI_R", name_cog))
select_df <- read_excel(select_file)

############################################################
# 去重（防止 join 放大）
############################################################

select_df <- select_df %>%
  distinct(FromNetwork, ToNetwork)

############################################################
# 筛选
############################################################

select_corr <- corr_df %>%
  inner_join(select_df,
             by = c("FromNetwork", "ToNetwork"))

############################################################
# 在 cluster 内部做 FDR
############################################################

select_corr <- select_corr %>%
  group_by(cluster) %>%
  mutate(
    p_FDR_cluster = p.adjust(p_value, method = "fdr")
  ) %>%
  ungroup()

############################################################
# 写入不同 sheet
############################################################

wb <- createWorkbook()

for (cl in unique(select_corr$cluster)) {
  
  temp_df <- select_corr %>%
    filter(cluster == cl)
  
  addWorksheet(wb, cl)
  writeData(wb, cl, temp_df)
}

saveWorkbook(wb, out_file, overwrite = TRUE)

cat("完成：已分 cluster 保存，并在每个 cluster 内部完成 FDR 校正\n")
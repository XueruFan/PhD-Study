############################################
# Correlation between rDCM and head motion
# XRF
############################################

rm(list = ls())

############################
# Load packages
############################
packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "stringr",
  "openxlsx"
)

sapply(packages, require, character.only = TRUE)

############################
# 1. Read head motion data
############################

fd_path <- "/Volumes/Zuolab_XRF/supplement/ccnp/ccnppek_fdmean0.3.xlsx"

fd_data <- read_excel(fd_path)

colnames(fd_data)[1] <- "ID"

# 取 session-level 平均头动（如果每个session有多个run）
fd_session <- fd_data %>%
  group_by(ID, Session) %>%
  summarise(fd_mean = mean(fd_mean, na.rm = TRUE), .groups = "drop")

sfei_path <- "/Volumes/Zuolab_XRF/output/norm_rDCM/rDCM_normative_data.xlsx"

sfei_data <- read_excel(sfei_path)

############################
# 3. Merge datasets
############################

merged_data <- sfei_data %>%
  left_join(fd_session, by = c("ID", "Session"))

# 检查是否有缺失
cat("Missing fd_mean:", sum(is.na(merged_data$fd_mean)), "\n")

############################
# 4. Correlation analysis
############################

# 找到所有 EC 列
ec_cols <- grep("^EC_", colnames(merged_data), value = TRUE)

# 计算相关
cor <- sapply(ec_cols, function(col) {
  cor(merged_data$fd_mean, merged_data[[col]], use = "complete.obs", method = "pearson")
})

p <- sapply(ec_cols, function(col) {
  p = cor.test(merged_data$fd_mean, merged_data[[col]], method = "pearson")$p.value
})

# 转成数据框
cor_df <- data.frame(
  Edge = ec_cols,
  cor,p
)

cor_df <- cor_df %>%
  mutate(
    p_fdr = p.adjust(p, method = "fdr")
  )

############################
# 6. Save results
############################

write.xlsx(
  cor_df,
  "/Volumes/Zuolab_XRF/output/norm_rDCM/rDCM_fd_correlation.xlsx",
  overwrite = TRUE
)

cat("Analysis finished.\n")

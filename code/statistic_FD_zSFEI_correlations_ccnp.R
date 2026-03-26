############################################
# Correlation between SFEI and head motion
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

############################
# 2. Read SFEI data
############################

sfei_path <- "/Volumes/Zuolab_XRF/output/normative/zSFEI_normative_data.xlsx"

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

# 按 Step × Network 做相关
cor_results <- merged_data %>%
  group_by(Step, Network) %>%
  summarise(
    n = sum(!is.na(fd_mean) & !is.na(zSFEI)),
    cor = cor(zSFEI, fd_mean, use = "complete.obs", method = "pearson"),
    p = cor.test(zSFEI, fd_mean, method = "pearson")$p.value,
    .groups = "drop"
  )

############################
# 5. Multiple comparison correction
############################

cor_results <- cor_results %>%
  mutate(
    p_fdr = p.adjust(p, method = "fdr")
  )

############################
# 6. Save results
############################

write.xlsx(
  cor_results,
  "/Volumes/Zuolab_XRF/output/normative/zSFEI_fd_correlation.xlsx",
  overwrite = TRUE
)

cat("Analysis finished.\n")

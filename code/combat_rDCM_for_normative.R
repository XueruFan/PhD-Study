############################################################
# ComBat harmonization for rDCM Effective Connectivity
# TD + ASD included
# Control covariate: Age
############################################################

rm(list = ls())

library(tidyverse)
library(readxl)
library(openxlsx)
library(neuroCombat)

############################################################
# 1. Load dataset
############################################################

data_path <- "/Volumes/Zuolab_XRF/output/norm_rDCM/rDCM_normative_data.xlsx"

ec_data <- read_excel(data_path)

############################################################
# 2. 预处理
############################################################

ec_data <- ec_data %>%
  mutate(
    subject = as.character(as.numeric(ID)),
    
    # 合并 NYU2
    Site = recode(Site, "NYU2" = "NYU"),
    
    Subtype = factor(Subtype, levels = c("TD","ASD-L","ASD-H"))
  )

############################################################
# 3. 构建 ComBat 输入矩阵
############################################################

meta_data <- ec_data %>%
  select(Cohort, ID, Session, Site, Age, Subtype)

ec_cols <- grep("^EC_", colnames(ec_data), value = TRUE)

dat_matrix <- as.matrix(ec_data[, ec_cols])

# ComBat 需要 features × subjects
dat_matrix <- t(dat_matrix)

############################################################
# 4. 设置 batch 与协变量
############################################################

batch <- as.factor(meta_data$Site)

# 注意：不要加入 Subtype
mod <- model.matrix(~ Age, data = meta_data)

############################################################
# 5. 删除 constant features
############################################################

feature_sd <- apply(dat_matrix, 1, sd, na.rm = TRUE)
non_constant_idx <- which(feature_sd != 0)

cat("Removed", sum(feature_sd == 0), "constant edges.\n")

dat_matrix <- dat_matrix[non_constant_idx, ]
ec_cols    <- ec_cols[non_constant_idx]

############################################################
# 6. 运行 ComBat
############################################################

combat_result <- neuroCombat(
  dat = dat_matrix,
  batch = batch,
  mod = mod,
  parametric = FALSE,
  eb = FALSE,
  mean.only = TRUE
)

############################################################
# 7. 获取 harmonized 数据
############################################################

dat_combat <- t(combat_result$dat.combat)

dat_combat <- as.data.frame(dat_combat)
colnames(dat_combat) <- ec_cols

harmonized_data <- bind_cols(meta_data, dat_combat)

############################################################
# 8. 保存结果
############################################################

output_path <- "/Volumes/Zuolab_XRF/output/norm_rDCM/rDCM_normative_data_combat.xlsx"

write.xlsx(harmonized_data, output_path, overwrite = TRUE)

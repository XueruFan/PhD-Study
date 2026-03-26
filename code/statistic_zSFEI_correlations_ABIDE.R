################################
# SFC embedding × Cognition correlation
# Stratified by subtype (L / H)
################################

rm(list = ls())

packages <- c(
  "tidyverse", "readxl", "openxlsx",
  "dplyr", "stats"
)
sapply(packages, require, character.only = TRUE)

## ==============================
## 路径
## ==============================

sfcDir  <- "/Volumes/Zuolab_XRF/output/abide/sfc/zsfei"
subjFile <- "/Volumes/Zuolab_XRF/output/abide/sfc/des/zSFEI_abide_id.csv"
phenoFile <- "/Volumes/Zuolab_XRF/supplement/abide/abide_A_all_240315.csv"
clusterFile  <- "/Volumes/Zuolab_XRF/output/abide/ABIDE_cluster_all_subjects.csv"

outDir <- "/Volumes/Zuolab_XRF/output/abide/sfc/stat/corr"
# plotDir <- "/Volumes/Zuolab_XRF/output/abide/sfc/plot/corr"

## ==============================
## 读取被试列表
## ==============================
subj_list <- read.csv(subjFile, header = FALSE, stringsAsFactors = FALSE)
colnames(subj_list)[1] <- "participant"
subj_list$participant <- as.character(subj_list$participant)

## ==============================
## 读取亚型
## ==============================
cluster <- read.csv(clusterFile, stringsAsFactors = FALSE)
cluster$subtype <- as.factor(cluster$subtype)
cluster$participant <- as.character(cluster$participant)

## ==============================
## 读取认知行为
## ==============================
pheno <- read.csv(phenoFile, stringsAsFactors = FALSE)
## ==============================
## 统一清洗无效数值（< 0 视为缺失）
## ==============================
pheno <- pheno %>%
  mutate(
    across(
      where(is.numeric),
      ~ ifelse(.x < 0, NA, .x)
    )
  )
colnames(pheno)[colnames(pheno) == "Participant"] <- "participant"
pheno$participant <- as.character(pheno$participant)


## SITE_ID 统一
pheno$SITE_ID <- gsub("ABIDEII-NYU_1|ABIDEII-NYU_2", "NYU", pheno$SITE_ID)
pheno$SITE_ID <- gsub("ABIDEII-KKI_1", "KKI", pheno$SITE_ID)
pheno$SITE_ID <- gsub("ABIDEII-SDSU_1", "SDSU", pheno$SITE_ID)
pheno$SITE_ID <- gsub("ABIDEII-UCLA_1|UCLA_1|UCLA_2", "UCLA", pheno$SITE_ID)
pheno$SITE_ID <- gsub("ABIDEII-GU_1", "GU", pheno$SITE_ID)
pheno$SITE_ID <- gsub("ABIDEII-UCD_1", "UCD", pheno$SITE_ID)
pheno$SITE_ID <- gsub("ABIDEII-EMC_1", "EMC", pheno$SITE_ID)
pheno$SITE_ID <- gsub("TRINITY|ABIDEII-TCD_1", "TCD", pheno$SITE_ID)
pheno$SITE_ID <- gsub("ABIDEII-USM_1", "USM", pheno$SITE_ID)
pheno$SITE_ID <- gsub("ABIDEII-IU_1", "IU", pheno$SITE_ID)
pheno$SITE_ID <- gsub("ABIDEII-U_MIA_1", "UMIA", pheno$SITE_ID)
pheno$SITE_ID <- gsub("ABIDEII-ETH_1", "ETH", pheno$SITE_ID)
pheno$SITE_ID <- gsub("UM_1|UM_2", "UM", pheno$SITE_ID)
pheno$SITE_ID <- gsub("ABIDEII-OHSU_1", "OHSU", pheno$SITE_ID)
pheno$SITE_ID <- gsub("STANFORD", "SU1", pheno$SITE_ID)
pheno$SITE_ID <- gsub("ABIDEII-SU_2", "SU2", pheno$SITE_ID)
pheno$SITE_ID <- gsub("LEUVEN_2", "KUL", pheno$SITE_ID)
pheno$SITE_ID <- gsub("CALTECH", "CALT", pheno$SITE_ID)

## ==============================
## 认知变量
## ==============================

names_cog_p <- c(
  # "FIQ",
  "ADOS_2_SOCAFFECT", "ADOS_2_TOTAL",
  "ADI_R_SOCIAL_TOTAL_A", "ADI_R_RRB_TOTAL_C",
  "SRS_TOTAL_RAW", "SRS_COGNITION_RAW",
  "SRS_AWARENESS_RAW", "SRS_COMMUNICATION_RAW",
  "SRS_MOTIVATION_RAW", "SRS_MANNERISMS_RAW"
)

names_cog_s <- c("ADOS_2_RRB")

## ==============================
## 主循环：step × subtype
## ==============================

all_results <- list()
files <- list.files(sfcDir, pattern = "step[0-9]+\\.xlsx", full.names = TRUE)

for (f in files) {
  # f <- "/Volumes/Zuolab_XRF/output/abide/sfc/sfc_embedding/step01.xlsx"
  step_name <- tools::file_path_sans_ext(basename(f))
  message("Processing ", step_name)
  
  sfc <- read_excel(f)
  colnames(sfc)[1] <- "participant"
  sfc$participant <- as.character(as.integer(sfc$participant))
  
  data_all <- subj_list %>%
    left_join(sfc, by = "participant") %>%
    left_join(pheno, by = "participant") %>%
    left_join(cluster, by = "participant")
  
  data_all$SITE_ID <- as.factor(data_all$SITE_ID)
  
  names_brain <- colnames(sfc)[-1]
  
  for (cl in c("1", "2")) {
    # cl <- "1"
    dat <- subset(data_all, subtype == cl)
    results <- data.frame()
    
    ## -------- Spearman（RRB）
    for (b in names_brain) {
      for (cog in names_cog_s) {
        
        temp <- dat[, c(b, cog, "SITE_ID")]
        temp <- temp[complete.cases(temp), ]
        
        if (nrow(temp) < 40) next
        
        if (length(unique(temp$SITE_ID)) > 1) {
          y_lm <- lm(temp[[cog]] ~ SITE_ID, data = temp)
          x_lm <- lm(temp[[b]] ~ SITE_ID, data = temp)
        } else {
          y_lm <- lm(temp[[cog]] ~ 1, data = temp)
          x_lm <- lm(temp[[b]] ~ 1, data = temp)
        }
        
        cor_test <- cor.test(residuals(y_lm), residuals(x_lm), method = "spearman")
        
        results <- rbind(results, data.frame(
          step = step_name,
          cluster = ifelse(cl == "1", "L", "H"),
          name_brain = b,
          name_cog = cog,
          coef = cor_test$estimate,
          p_value = cor_test$p.value,
          df = df.residual(y_lm)
        ))
      }
    }
    
    ## -------- Pearson（连续）
    for (b in names_brain) {
      for (cog in names_cog_p) {
        
        temp <- dat[, c(b, cog, "SITE_ID")]
        temp <- temp[complete.cases(temp), ]
        
        if (nrow(temp) < 40) next
        
        if (length(unique(temp$SITE_ID)) > 1) {
          y_lm <- lm(temp[[cog]] ~ SITE_ID, data = temp)
          x_lm <- lm(temp[[b]] ~ SITE_ID, data = temp)
        } else {
          y_lm <- lm(temp[[cog]] ~ 1, data = temp)
          x_lm <- lm(temp[[b]] ~ 1, data = temp)
        }
        
        cor_test <- cor.test(residuals(y_lm), residuals(x_lm), method = "pearson")
        
        results <- rbind(results, data.frame(
          step = step_name,
          cluster = ifelse(cl == "1", "L", "H"),
          name_brain = b,
          name_cog = cog,
          coef = cor_test$estimate,
          p_value = cor_test$p.value,
          df = df.residual(y_lm)
        ))
      }
    }
    
    results$P_adj <- p.adjust(results$p_value, method = "fdr")
    all_results[[paste0(step_name, "_", cl)]] <- results
  }
}

final_results <- bind_rows(all_results) %>% arrange(p_value)

write.csv(
  final_results,
  file.path(outDir, "zSFEI_cognition.csv"),
  row.names = FALSE
)


## =========================================================
## 后处理：网络重命名 + 显著性筛选 + 排序
## =========================================================

## -------- 1. 网络索引到缩写的映射表
network_map <- data.frame(
  name_brain = paste0("Net", sprintf("%02d", 1:15)),
  network_abbr = c(
    "VIS-P",   # 1  Visual Peripheral
    "CG-OP",   # 2  Cingulo-Opercular
    "DN-B",    # 3  Default Network-B
    "SMOT-B",  # 4  Somatomotor-B
    "AUD",     # 5  Auditory
    "PM-PPr",  # 6  Premotor–Posterior Parietal rostral
    "dATN-B",  # 7  Dorsal Attention-B
    "SMOT-A",  # 8  Somatomotor-A
    "LANG",    # 9  Language
    "FPN-B",   # 10 Frontoparietal Network-B
    "FPN-A",   # 11 Frontoparietal Network-A
    "dATN-A",  # 12 Dorsal Attention-A
    "VIS-C",   # 13 Visual Central
    "SAL/PMN", # 14 Salience / Parietal Memory Network
    "DN-A"     # 15 Default Network-A
  ),
  stringsAsFactors = FALSE
)

## -------- 2. 合并网络名称
final_results_named <- final_results %>%
  left_join(network_map, by = "name_brain")

## -------- 3. 只保留显著结果
final_results_sig <- final_results_named %>%
  filter(p_value < 0.05)

## -------- 4. 提取 step 数字用于排序
final_results_sig <- final_results_sig %>%
  mutate(
    step_num = as.numeric(gsub("step", "", step)),
    cluster = factor(cluster, levels = c("L", "H"))
  )

## -------- 5. 按 步数 → 网络 → 亚型 排序
final_results_sig <- final_results_sig %>%
  arrange(step_num, network_abbr, cluster)

## -------- 6. 清理辅助列
final_results_sig <- final_results_sig %>%
  select(
    step,
    cluster,
    network_abbr,
    name_cog,
    coef,
    p_value,
    P_adj,
    df
  )

## -------- 7. 保存最终显著结果表
write.csv(
  final_results_sig,
  file.path(outDir, "zSFEI_cognition_significant.csv"),
  row.names = FALSE
)


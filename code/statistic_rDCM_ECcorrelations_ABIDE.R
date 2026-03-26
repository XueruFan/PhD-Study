############################################################
# All EC edges × Cognition Correlation
# Stratified by Subtype (ASD-L / ASD-H)
# GLOBAL FDR correction
############################################################

rm(list = ls())

packages <- c(
  "tidyverse", "readxl", "openxlsx",
  "dplyr", "stats"
)
sapply(packages, require, character.only = TRUE)

set.seed(1205)

############################################################
# 路径
############################################################

ecFile     <- "/Volumes/Zuolab_XRF/output/rDCM_PLSDA/PLSDA_statistics.xlsx"
phenoFile  <- "/Volumes/Zuolab_XRF/supplement/abide/abide_A_all_240315.csv"

outDir <- "/Volumes/Zuolab_XRF/output/rDCM_PLSDA/corr"

###
site <- read_excel("/Volumes/Zuolab_XRF/output/norm_rDCM/rDCM_normative_data.xlsx") %>%
  dplyr::select(ID, Site) %>%
  mutate(Site = recode(Site, "NYU2" = "NYU"))
         
colnames(site) <- c("participant", "site")

############################################################
# 读取 EC 数据
############################################################

scores <- read_excel(ecFile, sheet = "Scores")

dcm_z <- read_excel("/Volumes/Zuolab_XRF/output/norm_rDCM/gamlss/rDCM_deviation_zscores.xlsx")

dcm_z$participant <- as.character(as.numeric(dcm_z$subject))
scores$participant  <- as.character(scores$subject)

# 合并 subtype / site 信息
dcm_all <- scores %>% 
  left_join(
    site, by = "participant"
  ) %>%
  dplyr::select(participant, site) %>% 
  left_join(
    dcm_z, by = "participant"
  )

############################################################
# 提取 EC 边
############################################################

ec_cols <- grep("^EC_", colnames(dcm_all), value = TRUE)

############################################################
# 读取行为数据
############################################################

pheno <- read.csv(phenoFile, stringsAsFactors = FALSE)
colnames(pheno)[colnames(pheno) == "Participant"] <- "participant"
pheno$participant <- as.character(pheno$participant)

pheno <- pheno %>%
  mutate(across(where(is.numeric), ~ ifelse(.x < 0, NA, .x)))

############################################################
# 行为变量
############################################################

names_cog_p <- c(
  "ADOS_2_SOCAFFECT", "ADOS_2_TOTAL",
  "ADI_R_SOCIAL_TOTAL_A", "ADI_R_RRB_TOTAL_C",
  "SRS_TOTAL_RAW", "SRS_COGNITION_RAW",
  "SRS_AWARENESS_RAW", "SRS_COMMUNICATION_RAW",
  "SRS_MOTIVATION_RAW", "SRS_MANNERISMS_RAW"
)

names_cog_s <- c("ADOS_2_RRB")

############################################################
# 合并数据
############################################################

data_all <- dcm_all %>%
  left_join(pheno, by = "participant")

data_all$Site <- as.factor(data_all$site)
data_all$Subtype <- factor(data_all$Subtype)

############################################################
# 主循环
############################################################

all_results <- list()

for (cl in c("ASD-L", "ASD-H")) {
  
  cat("Processing cluster:", cl, "\n")
  
  dat_sub <- data_all %>%
    filter(Subtype == cl)
  
  for (edge in ec_cols) {
    
    ########################################################
    ## Spearman
    ########################################################
    
    for (cog in names_cog_s) {
      
      temp <- dat_sub[, c(edge, cog, "Site")]
      temp <- temp[complete.cases(temp), ]
      if (nrow(temp) < 30) next
      
      if (length(unique(temp$Site)) > 1) {
        y_lm <- lm(temp[[cog]] ~ Site, data = temp)
        x_lm <- lm(temp[[edge]] ~ Site, data = temp)
      } else {
        y_lm <- lm(temp[[cog]] ~ 1, data = temp)
        x_lm <- lm(temp[[edge]] ~ 1, data = temp)
      }
      
      cor_test <- cor.test(
        residuals(x_lm),
        residuals(y_lm),
        method = "spearman"
      )
      
      all_results[[length(all_results)+1]] <- data.frame(
        edge = edge,
        cluster = cl,
        name_cog = cog,
        coef = cor_test$estimate,
        p_value = cor_test$p.value,
        n = nrow(temp),
        df = df.residual(y_lm),
        method = "spearman"
      )
    }
    
    ########################################################
    ## Pearson
    ########################################################
    
    for (cog in names_cog_p) {
      
      temp <- dat_sub[, c(edge, cog, "Site")]
      temp <- temp[complete.cases(temp), ]
      if (nrow(temp) < 30) next
      
      if (length(unique(temp$Site)) > 1) {
        y_lm <- lm(temp[[cog]] ~ Site, data = temp)
        x_lm <- lm(temp[[edge]] ~ Site, data = temp)
      } else {
        y_lm <- lm(temp[[cog]] ~ 1, data = temp)
        x_lm <- lm(temp[[edge]] ~ 1, data = temp)
      }
      
      cor_test <- cor.test(
        residuals(x_lm),
        residuals(y_lm),
        method = "pearson"
      )
      
      all_results[[length(all_results)+1]] <- data.frame(
        edge = edge,
        cluster = cl,
        name_cog = cog,
        coef = cor_test$estimate,
        p_value = cor_test$p.value,
        n = nrow(temp),
        df = df.residual(y_lm),
        method = "pearson"
      )
    }
  }
}

############################################################
# 合并结果 + FDR
############################################################

final_results <- bind_rows(all_results)

final_results <- final_results %>%
  group_by(cluster) %>%
  mutate(P_adj_cluster = p.adjust(p_value, method = "fdr")) %>%
  ungroup()



############################################################
# 建立节点 → 网络映射
############################################################

network_labels_custom <- c(
  5,13,1,
  4,8,
  14,
  6,
  2,
  7,12,
  10,11,
  3,15,
  9
)

network_networks_custom <- c(
  "AUD",
  "VIS-C","VIS-P",
  "SMOT-B","SMOT-A",
  "SAL",
  "PM-PPr",
  "AN",
  "dATN-B","dATN-A",
  "FPN-B","FPN-A",
  "DN-B","DN-A",
  "LANG"
)

network_map <- tibble(
  Label = network_labels_custom,
  Network = network_networks_custom
)

############################################################
# 解析边的 From / To
############################################################

final_results <- final_results %>%
  mutate(
    From = as.integer(stringr::str_extract(edge, "(?<=EC_)\\d+")),
    To   = as.integer(stringr::str_extract(edge, "(?<=to_)\\d+"))
  )

############################################################
# 加入网络名称
############################################################

final_results <- final_results %>%
  left_join(network_map, by = c("From" = "Label")) %>%
  rename(FromNetwork = Network) %>%
  left_join(network_map, by = c("To" = "Label")) %>%
  rename(ToNetwork = Network)

final_results <- final_results %>%
  dplyr::select(
    edge,
    From, FromNetwork,
    To, ToNetwork,
    cluster,
    name_cog,
    coef,
    p_value,
    # P_adj_cluster,
    n,
    df,
    method
  )

final_results <- final_results %>%
  arrange(p_value)

############################################################
# 保存完整结果
############################################################

write.csv(
  final_results,
  file.path(outDir, "EC_cognition_correlation.csv"),
  row.names = FALSE
)

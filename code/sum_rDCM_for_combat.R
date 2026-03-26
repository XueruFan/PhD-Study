rm(list = ls())

library(tidyverse)
library(readxl)
library(openxlsx)

############################################################
# 1. Paths
############################################################

ccnp_path  <- "/Volumes/Zuolab_XRF/output/ccnp/dcm/sum/pek_rdcm_fd0.3_sessionAvg.xlsx"
abide_path <- "/Volumes/Zuolab_XRF/output/abide/dcm/sum/ABIDE_rDCM_summary.xlsx"

abide_demo_path <- "/Volumes/Zuolab_XRF/output/abide/sfc/des/zSFEI_abide_demo.csv"

############################################################
# 2. Read CCNP rDCM
############################################################

ccnp_data <- read.xlsx(ccnp_path)

ccnp_data <- ccnp_data %>%
  rename(ID = Participant) %>%
  mutate(
    ID      = as.character(ID),
    Cohort  = "CCNP",
    Site    = "PEK",
    Sex     = ifelse(Sex == "男", 1, 2),
    Subtype = "TD"
  )

############################################################
# 3. Read ABIDE rDCM
############################################################

abide_data <- read_excel(abide_path)

abide_data <- abide_data %>%
  rename(ID = subject) %>%
  mutate(
    ID     = as.character(as.numeric(ID)),
    Cohort = "ABIDE"
  )

############################################################
# 4. Read ABIDE demo
############################################################

abide_demo <- read.csv(abide_demo_path) %>%
  rename(ID = Subject) %>%
  mutate(
    ID   = as.character(ID),
    Sex  = ifelse(Sex == "Male", 1, 2),
    Site = as.character(site)
  )

############################################################
# 5. Join ABIDE demo
############################################################

abide_data <- abide_data %>%
  left_join(abide_demo, by = "ID")

############################################################
# 6. Filter subjects (保持与你之前一致)
############################################################

ccnp_clean <- ccnp_data %>%
  filter(
    Sex == 1,
    Age <= 18
  )

abide_clean <- abide_data %>%
  filter(
    Sex == 1,
    Age <= 18
  )

############################################################
# 7. Align columns
############################################################

# 找到 EC 列
ec_cols_ccnp  <- grep("^EC_", colnames(ccnp_clean), value = TRUE)
ec_cols_abide <- grep("^EC_", colnames(abide_clean), value = TRUE)

ec_cols <- intersect(ec_cols_ccnp, ec_cols_abide)

ccnp_clean <- ccnp_clean %>%
  select(Cohort, ID, Session, Subtype, Site, Age, all_of(ec_cols))

abide_clean <- abide_clean %>%
  select(Cohort, ID, Subtype, Site, Age, all_of(ec_cols))

############################################################
# 8. Merge
############################################################

normative_data <- bind_rows(ccnp_clean, abide_clean)

############################################################
# 9. Save
############################################################

output_file <- "/Volumes/Zuolab_XRF/output/norm_rDCM/rDCM_normative_data.xlsx"

write.xlsx(normative_data, output_file, overwrite = TRUE)

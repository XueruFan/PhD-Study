rm(list = ls())

library(tidyverse)
library(readxl)
library(openxlsx)

############################################################
# 1. Paths
############################################################

ccnp_path  <- "/Volumes/Zuolab_XRF/output/ccnp/sfc/zsfei"
abide_path <- "/Volumes/Zuolab_XRF/output/abide/sfc/zsfei"

abide_demo_path <- "/Volumes/Zuolab_XRF/output/abide/sfc/des/zSFEI_abide_demo.csv"

############################################################
# 2. Read CCNP (sessionAvg files only)
############################################################

read_ccnp <- function(folder_path) {
  
  files <- list.files(
    folder_path,
    pattern = "^pek_.*_sessionAvg\\.xlsx$",
    full.names = TRUE
  )
  
  map_dfr(files, function(f) {
    
    step_num <- stringr::str_extract(basename(f), "\\d+")
    dat <- read.xlsx(f)
    
    dat$Step   <- as.numeric(step_num)
    dat$Cohort <- "CCNP"
    dat$Site   <- "PEK"   # CCNP 单站点
    
    dat
  })
}

ccnp_embedding <- read_ccnp(ccnp_path)

############################################################
# 3. Read ABIDE embedding
############################################################

read_abide <- function(folder_path) {
  
  files <- list.files(
    folder_path,
    pattern = "^step\\d+\\.xlsx$",
    full.names = TRUE
  )
  
  map_dfr(files, function(f) {
    
    step_num <- stringr::str_extract(basename(f), "\\d+")
    dat <- read.xlsx(f)
    
    dat$Step   <- as.numeric(step_num)
    dat$Cohort <- "ABIDE"
    
    dat
  })
}

abide_embedding <- read_abide(abide_path)

############################################################
# 4. Read ABIDE demo (site already included)
############################################################

abide_demo <- read.csv(abide_demo_path) %>%
  rename(ID = Subject) %>%
  mutate(
    ID   = as.integer(ID),
    Sex  = ifelse(Sex == "Male", 1, 2),
    Site = as.character(site)
  )

############################################################
# 5. Clean CCNP
############################################################

ccnp_clean <- ccnp_embedding %>%
  rename(ID = Participant) %>%
  mutate(
    ID  = as.character(ID),
    Sex = ifelse(Sex == "男", 1, 2),
    Subtype = "TD"
  ) %>%
  filter(
    Sex == 1,
    Age <= 18,
    Step >= 1,
    Step <= 7
  )

############################################################
# 6. Clean ABIDE (直接 join demo)
############################################################

abide_clean <- abide_embedding %>%
  rename(ID = Subject) %>%
  mutate(ID = as.integer(ID)) %>%
  left_join(abide_demo, by = "ID") %>%
  filter(
    Sex == 1,
    Age <= 18,
    Step >= 1,
    Step <= 7
  )

############################################################
# 7. Pivot long
############################################################

ccnp_long <- ccnp_clean %>%
  pivot_longer(
    cols = starts_with("Net"),
    names_to = "Network",
    values_to = "zSFEI"
  ) %>%
  mutate(zSFEI = as.numeric(zSFEI)) %>%
  select(Cohort, ID, Session, Subtype, Site, Age, Step, Network, zSFEI)

abide_long <- abide_clean %>%
  pivot_longer(
    cols = starts_with("Net"),
    names_to = "Network",
    values_to = "zSFEI"
  ) %>%
  mutate(zSFEI = as.numeric(zSFEI),
         ID  = as.character(ID),) %>%
  select(Cohort, ID, Subtype, Site, Age, Step, Network, zSFEI)

abide_long$Session <- NA

############################################################
# 8. Merge
############################################################

normative_data <- bind_rows(ccnp_long, abide_long)

############################################################
# 9. Save
############################################################

output_file <- "/Volumes/Zuolab_XRF/output/normative/zSFEI_normative_data.xlsx"

write.xlsx(normative_data, output_file, overwrite = TRUE)


############################################################
# Demographic summary statistics
# Observation-level & Subject-level
############################################################

library(dplyr)
library(openxlsx)

############################################################
# 1. Observation-level summary
############################################################

# 观测层级：每一行代表一个观测
obs_level <- normative_data %>%
  filter(Step == 1) %>%
  select(ID, Cohort, Site, Age, Subtype) %>%
  distinct()  # 防止 Step × Network 重复

# ---- 按 Cohort ----
obs_cohort <- obs_level %>%
  group_by(Cohort, Subtype) %>%
  summarise(
    N_obs = n(),
    Age_mean = mean(Age, na.rm = TRUE),
    Age_sd   = sd(Age, na.rm = TRUE),
    .groups = "drop"
  )

# ---- 按 Cohort × Site ----
obs_site <- obs_level %>%
  group_by(Cohort, Subtype, Site) %>%
  summarise(
    N_obs = n(),
    Age_mean = mean(Age, na.rm = TRUE),
    Age_sd   = sd(Age, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Cohort, Site)


############################################################
# 3. Save to Excel (multiple sheets)
############################################################

demo_output_file <- "/Volumes/Zuolab_XRF/output/norm_zSFEI/zSFEI_normative_demo.xlsx"

wb <- createWorkbook()

addWorksheet(wb, "Cohort")
writeData(wb, "Cohort", obs_cohort)

addWorksheet(wb, "Site")
writeData(wb, "Site", obs_site)

saveWorkbook(wb, demo_output_file, overwrite = TRUE)

cat("Demographic summary saved to:\n", demo_output_file)

#########保存一份ccnp用于建常模的人口信息

ccnppek <- normative_data %>%
  filter(Step == 1, Cohort == "CCNP") %>%
  select(ID, Age, Session, Subtype) %>%
  distinct()  # 防止 Step × Network 重复

write.csv(ccnppek, "/Volumes/Zuolab_XRF/output/ccnp/CCNPPEK_Participant.csv")

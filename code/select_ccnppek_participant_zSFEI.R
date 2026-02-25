rm(list = ls())

library(readxl)
library(dplyr)
library(stringr)
library(writexl)

# ------------------ paths ------------------
qc_file   <- "/Volumes/Zuolab_XRF/supplement/ccnp/ccnppek_fdmean0.3.xlsx"
demo_file <- "/Volumes/Zuolab_XRF/supplement/ccnp/ccnppek_participant.xlsx"

sfc_dir <- "/Volumes/Zuolab_XRF/output/ccnp/sfc/zsfei"

steps <- sprintf("%02d", 1:6)

# ------------------ PEK QC table (run-level) ------------------
qc <- read_xlsx(qc_file) %>%
  mutate(
    Session = str_pad(Session, 2, pad = "0"),
    Run     = str_pad(Run, 2, pad = "0")
  ) %>%
  distinct(Participant, Session, Run, .keep_all = TRUE)

# ------------------ PEK demographic table (session-level) ------------------
demo <- read_xlsx(demo_file) %>%
  mutate(
    Participant = str_pad(as.numeric(Participant), 4, pad = "0"),
    Session = str_pad(Session, 2, pad = "0")
  ) %>%
  select(Participant, Session, Sex, Age) %>%
  distinct(Participant, Session, .keep_all = TRUE)

# ------------------ loop over steps ------------------
for (s in steps) {
  
  message("Processing PEK step ", s)
  
  sfc <- read_xlsx(file.path(sfc_dir, paste0("step", s, ".xlsx"))) %>%
    filter(Site == "PEK") %>%
    mutate(
      Run = str_pad(Run, 2, pad = "0")
    ) %>%
    rename(
      Participant = Subject
    ) %>%
    mutate(
      across(
        starts_with("Net"),
        ~ as.numeric(.x)
      )
    )
  
  # 1. 严格 run-level QC
  sfc_sel <- sfc %>%
    semi_join(qc, by = c("Participant", "Session", "Run"))
  
  # 2. 对同一个 Participant × Session 内的 run 取平均
  sfc_avg <- sfc_sel %>%
    group_by(Site, Participant, Session) %>%
    summarise(
      across(starts_with("Net"), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  # 3. 合并人口学信息（session-level）
  merged <- sfc_avg %>%
    left_join(demo, by = c("Participant", "Session")) %>%
    select(Site, Participant, Sex, Age, Session, starts_with("Net"))
  
  # 4. 输出
  write_xlsx(
    merged,
    file.path(sfc_dir, paste0("pek_step", s, "_fd0.3_sessionAvg.xlsx"))
  )
}

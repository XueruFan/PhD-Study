# ============================================================
# ABIDE SFC embedding
# Pairwise contrasts (TD, ASD subtypes, and L vs H)
# 4-panel heatmap with raw-p significance
# ============================================================

rm(list = ls())

# ----------------------------
# Libraries
# ----------------------------
library(tidyverse)
library(readxl)
library(emmeans)
library(showtext)

# ----------------------------
# Font (macOS, optional)
# ----------------------------
font_add(
  family = "pingfang",
  regular = "/System/Library/Fonts/PingFang.ttc"
)
showtext_auto()

# ----------------------------
# Paths
# ----------------------------
embedding_dir <- "/Volumes/Zuolab_XRF/output/abide/sfc/zsfei"
demo_path     <- "/Volumes/Zuolab_XRF/supplement/abide/abide_demo.xlsx"
site_dir      <- "/Volumes/Zuolab_XRF/data/abide/sublist"
cluster_path  <- "/Volumes/Zuolab_XRF/output/abide/ABIDE_cluster_all_subjects.csv"

out_dir <- "/Volumes/Zuolab_XRF/output/abide/sfc/stat/difference"
plot_dir <- "/Volumes/Zuolab_XRF/output/abide/sfc/plot/difference"

dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

step_max <- 6
nNet     <- 15

# ----------------------------
# Network labels
# ----------------------------
network_map <- tibble(
  Network = c(
    "Net05", 
    "Net13", "Net01", 
    "Net04", "Net08",
    "Net14", 
    "Net06",
    "Net02", 
    "Net07", "Net12",
    "Net10", "Net11",
    "Net03", "Net15", 
    "Net09"
  ),
  NetworkLabel = c(
    "AUD",
    "VIS-C", "VIS-P",
    "SMOT-B", "SMOT-A",
    "SAL/PMN",
    "PM-PPr",
    "AN",
    "dATN-B", "dATN-A",
    "FPN-B", "FPN-A",
    "DN-B", "DN-A",
    "LANG"
  ),
  Order = 1:15
)

# ----------------------------
# Demographics
# ----------------------------
demo <- read_xlsx(demo_path) %>%
  transmute(
    Subject = as.character(Participant),
    AGE_AT_SCAN,
    SEX = factor(SEX, levels = c(1,2), labels = c("Male","Female"))
  )

# ----------------------------
# Subtype labels
# ----------------------------
cluster <- read_csv(cluster_path, show_col_types = FALSE) %>%
  transmute(
    Subject = as.character(participant),
    subtype = as.integer(subtype)
  ) %>%
  mutate(
    Group = factor(
      subtype,
      levels = c(0, 1, 2),
      labels = c("TD", "L", "H")
    )
  )

# ----------------------------
# Site mapping
# ----------------------------
site_map <- list.files(
  site_dir,
  pattern = "^subjects_.*\\.list$",
  full.names = TRUE
) %>%
  map_dfr(function(f) {
    
    site <- basename(f) %>%
      str_remove("^subjects_") %>%
      str_remove("\\.list$") %>%
      toupper()
    
    subjects <- read_lines(f) %>%
      as.integer() %>%
      as.character()
    
    tibble(Subject = subjects, site = site)
  })

# ----------------------------
# Read embedding
# ----------------------------
read_one_step <- function(step) {
  
  df <- read_excel(file.path(embedding_dir, sprintf("step%02d.xlsx", step)))
  
  df %>%
    pivot_longer(
      cols = starts_with("Net"),
      names_to  = "Network",
      values_to = "zSFEI"
    ) %>%
    mutate(
      Step    = step,
      Subject = as.character(as.integer(Subject)),
      Network = factor(Network, levels = sprintf("Net%02d", 1:nNet))
    )
}

embedding_long <- map_dfr(1:step_max, read_one_step)

# ----------------------------
# Merge all
# ----------------------------
data_all <- embedding_long %>%
  left_join(demo,    by = "Subject") %>%
  left_join(cluster, by = "Subject") %>%
  left_join(site_map,by = "Subject") %>%
  filter(
    SEX == "Male",
    !is.na(Group),
    !is.na(site)
  ) %>%
  mutate(
    Group = relevel(Group, ref = "TD"),
    zSFEI = as.numeric(zSFEI)
  )

# ----------------------------
# Merged ASD group
# ----------------------------
data_all <- data_all %>%
  mutate(
    Group_ASD = factor(
      if_else(Group == "TD", "TD", "ASD"),
      levels = c("TD","ASD")
    )
  )

# ============================================================
# Contrast functions
# ============================================================

fit_asd_vs_td <- function(df) {
  
  m <- lm(zSFEI ~ Group_ASD + AGE_AT_SCAN + site, data = df)
  em <- emmeans(m, ~ Group_ASD)
  
  as.data.frame(
    contrast(em, "trt.vs.ctrl", ref = "TD")
  ) %>%
    mutate(Contrast = "ASD vs TD")
}

fit_subtype_vs_td <- function(df) {
  
  m <- lm(zSFEI ~ Group + AGE_AT_SCAN + site, data = df)
  em <- emmeans(m, ~ Group)
  
  as.data.frame(
    contrast(em, "trt.vs.ctrl", ref = "TD")
  ) %>%
    filter(contrast %in% c("L - TD","H - TD")) %>%
    mutate(
      Contrast = recode(
        contrast,
        "L - TD" = "L vs TD",
        "H - TD" = "H vs TD"
      )
    )
}

## -------- NEW: L vs H (ASD only) --------
fit_L_vs_H <- function(df) {
  
  df_asd <- df %>% filter(Group %in% c("L","H"))
  
  if (n_distinct(df_asd$Group) < 2) return(NULL)
  
  m <- lm(zSFEI ~ Group + AGE_AT_SCAN + site, data = df_asd)
  em <- emmeans(m, ~ Group)
  
  as.data.frame(
    contrast(em, method = list("L vs H" = c(1, -1)))
  ) %>%
    mutate(Contrast = "L vs H")
}

# ============================================================
# Run all contrasts
# ============================================================

stats <- bind_rows(
  
  data_all %>%
    group_by(Network, Step) %>%
    group_modify(~ fit_asd_vs_td(.x)),
  
  data_all %>%
    group_by(Network, Step) %>%
    group_modify(~ fit_subtype_vs_td(.x)),
  
  data_all %>%
    group_by(Network, Step) %>%
    group_modify(~ fit_L_vs_H(.x))
  
) %>%
  ungroup() %>%
  mutate(
    sig_raw = p.value < 0.05,
    Contrast = factor(
      Contrast,
      levels = c("ASD vs TD","L vs TD","H vs TD","L vs H")
    )
  )

write_csv(stats, file.path(out_dir, "zSFEI_abide_network_step.csv"))

# ============================================================
# Plot data
# ============================================================

plot_data <- stats %>%
  left_join(network_map, by = "Network") %>%
  filter(!is.na(NetworkLabel)) %>%
  mutate(
    Step = factor(Step),
    NetworkLabel = factor(NetworkLabel, levels = network_map$NetworkLabel)
  )

# ============================================================
# Plot
# ============================================================

p <- ggplot(
  plot_data,
  aes(x = Step, y = NetworkLabel, fill = estimate)
) +
  geom_tile(color = "lightgray", linewidth = 0.4) +
  geom_point(
    data = subset(plot_data, sig_raw),
    shape = 21, size = 2.6, stroke = 1,
    fill = NA, color = "black"
  ) +
  facet_wrap(~ Contrast, nrow = 1) +
  scale_fill_gradient2(
    low      = "#3d5c6f",
    mid      = "white",
    high     = "#e47159",
    midpoint = 0,
    limits   = c(-0.4, 0.4),
    name = expression(Delta~zSFEI)
  ) +
  labs(x = "SFC步数", y = "功能网络") +
  theme_bw(base_size = 22) +
  theme(
    strip.text   = element_text(size = 32, face = "bold"),
    axis.text.y  = element_text(size = 30),
    axis.text.x  = element_text(size = 30),
    axis.title   = element_text(size = 35),
    legend.title = element_text(size = 35),
    legend.text  = element_text(size = 30),
    panel.grid   = element_blank()
  )

ggsave(
  filename = file.path(plot_dir, "三组人各功能网络在不同SFC连接步数下的zSFEI差异.png"),
  plot     = p,
  width    = 4000,
  height   = 2200,
  dpi      = 300,
  units    = "px"
)

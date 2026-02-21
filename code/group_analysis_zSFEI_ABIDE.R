
rm(list = ls())

# ----------------------------
# 0. Libraries
# ----------------------------
library(tidyverse)
library(readxl)
library(showtext)
library(viridis)
library(ggplot2)
library(dplyr)
# ----------------------------
# 中文字体设置（macOS）
# ----------------------------
font_add(
  family = "pingfang",
  regular = "/System/Library/Fonts/PingFang.ttc"
)
showtext_auto()

theme_set(
  theme_bw(base_family = "pingfang")
)

# ----------------------------
# 1. Paths & parameters
# ----------------------------
embedding_dir <- "/Volumes/Zuolab_XRF/output/abide/sfc/zsfei"
demo_path     <- "/Volumes/Zuolab_XRF/supplement/abide/abide_demo.xlsx"
site_dir      <- "/Volumes/Zuolab_XRF/data/abide/sublist"
cluster_path  <- "/Volumes/Zuolab_XRF/output/abide/ABIDE_cluster_all_subjects.csv"

out_root <- "/Volumes/Zuolab_XRF/output/abide/sfc/des"
plot_dir <- "/Volumes/Zuolab_XRF/output/abide/sfc/plot"

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

step_max <- 6
nNet     <- 15

# ----------------------------
# 2. Read demographic information
# ----------------------------
demo <- read_xlsx(demo_path) %>%
  select(
    Participant,
    AGE_AT_SCAN,
    SEX
  ) %>%
  mutate(
    Participant = as.character(Participant),
    SEX = factor(
      SEX,
      levels = c(1, 2),
      labels = c("Male", "Female")
    )
  )

# ----------------------------
# 3. Read ASD subtype predictions (MALE ONLY)
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
      labels = c("TD", "ASD-L", "ASD-H")
    )
  )

# ----------------------------
# 4. Build subject → site mapping
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
      as.character() %>%
      as.integer() %>%
      as.character()
    
    tibble(
      Subject = subjects,
      site    = site
    )
  })

# ----------------------------
# 5. Read & merge embedding (all steps)
# ----------------------------
read_one_step <- function(step) {
  
  fname <- sprintf("step%02d.xlsx", step)
  fpath <- file.path(embedding_dir, fname)
  
  df <- read_excel(fpath)
  
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
# 6. Merge all information
# ----------------------------
data_all <- embedding_long %>%
  left_join(demo,     by = c("Subject" = "Participant")) %>%
  left_join(cluster,  by = "Subject") %>%
  left_join(site_map, by = "Subject") %>%
  filter(
    SEX == "Male",
    !is.na(Group),
    !is.na(AGE_AT_SCAN),
    !is.na(site)
  ) %>%
  mutate(
    zSFEI = as.numeric(zSFEI)
  )

# ============================================================
# participant for analysis
# ============================================================
participant_for_analysis <- data.frame(unique(data_all$Subject))
write_csv(participant_for_analysis, file.path(out_root, "zSFEI_abide_id.csv"), col_names = F)

participant_summary <- participant_for_analysis %>%
  rename(Subject = unique.data_all.Subject.) %>%  # 如果列名是自动生成的
  left_join(
    demo %>%
      transmute(
        Subject = as.character(Participant),
        Age     = AGE_AT_SCAN,
        Sex     = SEX
      ),
    by = "Subject"
  ) %>%
  left_join(
    cluster %>%
      transmute(
        Subject,
        Subtype = Group
      ),
    by = "Subject"
  ) %>%
  left_join(
    site_map,
    by = "Subject"
  ) %>%
  arrange(Subtype, Subject)

# 查看前几行
head(participant_summary)

# 如需导出
write_csv(
  participant_summary,
  file.path(out_root, "zSFEI_abide_demo.csv")
)


############################################################
# Demographic summary: Overall + Site × Group
############################################################

data_clean <- participant_summary %>%
  mutate(site = recode(site, "NYU2" = "NYU")) %>%
  distinct(Subject, Subtype, Age, site)

overall_summary <- data_clean %>%
  group_by(Subtype) %>%
  summarise(
    N        = n(),
    Age_mean = mean(Age, na.rm = TRUE),
    Age_sd   = sd(Age, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(site = "ALL") %>%
  select(site, everything())

site_summary <- data_clean %>%
  group_by(Subtype, site) %>%
  summarise(
    N        = n(),
    Age_mean = mean(Age, na.rm = TRUE),
    Age_sd   = sd(Age, na.rm = TRUE),
    .groups  = "drop"
  )

final_summary <- bind_rows(
  site_summary,
  overall_summary
) %>%
  arrange(Subtype, site)

final_summary <- final_summary %>%
  mutate(
    Age_mean = round(Age_mean, 2),
    Age_sd   = round(Age_sd, 2)
  )

write_csv(
  final_summary,
  file.path(out_root, "zSFEI_abide_summary.csv")
)


# ============================================================
# Network-wise descriptive stats
# ============================================================
desc_network <- data_all %>%
  group_by(Group, Network) %>%
  summarise(
    N                = n(),
    mean_zSFEI   = mean(zSFEI),
    sd_zSFEI     = sd(zSFEI),
    median_zSFEI = median(zSFEI),
    # iqr_embedding    = IQR(zSFEI),
    .groups = "drop"
  )

write_csv(desc_network, file.path(out_root, "zSFEI_abide_network.csv"))

# ============================================================
# Step-wise descriptive stats
# ============================================================
desc_step <- data_all %>%
  group_by(Group, Step) %>%
  summarise(
    N              = n(),
    mean_zSFEI = mean(zSFEI),
    sd_zSFEI   = sd(zSFEI),
    .groups = "drop"
  )

write_csv(desc_step, file.path(out_root, "zSFEI_abide_step.csv"))

### 按照网络-分组显示zSFEI
library(scico)

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

data_to_plot <- data_all %>%
  left_join(
    network_map %>% select(Network, NetworkLabel),
    by = "Network"
  )

data_to_plot$NetworkLabel <- factor(
  data_to_plot$NetworkLabel,
  levels = network_map$NetworkLabel
)


p2 <- ggplot(data_to_plot,
             aes(x = NetworkLabel, y = zSFEI, fill = Group)) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    linewidth = 0.4,
    position = position_dodge(width = 0.7)
  ) +
  # stat_summary(
  #   fun = median,
  #   geom = "point",
  #   size = 1.6,
  #   colour = "black",
  #   position = position_dodge(width = 0.65)
  # ) +
  scale_fill_manual(
    values = c(
      "TD"    = "white",
      "ASD-L" = "#86b5a1",
      "ASD-H" = "#f9ae78"
    )
  ) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  labs(
    x = "功能网络",
    y = "zSFEI",
    fill = "组别/亚型"
  ) +
  theme_classic(base_size = 22) +
  theme(
    # plot.title   = element_text(size = 40, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 30),
    axis.text.y = element_text(size = 30),
    legend.position = "top",
    plot.title   = element_blank(),
    axis.title   = element_text(size = 35),
    legend.title = element_text(size = 35),
    legend.text  = element_text(size = 35),
  )

p2

ggsave(
  filename = file.path(plot_dir, "zSFEI_network_abide.png"),
  plot     = p2,
  width    = 3000,
  height   = 2000,
  dpi      = 300,
  units    = "px"
)

#######
## network分面图 1-5步，第6步没有zSFEI了

data_to_plot <- data_to_plot %>%
  filter(Step <= 5) %>%
  mutate(
    Step = factor(Step, levels = 1:5)
  )

p3 <- ggplot(data_to_plot,
                  aes(x = Step, y = zSFEI, fill = Group)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    position = position_dodge(width = 0.7),
    linewidth = 0.3
  ) +
  facet_wrap(~ NetworkLabel, ncol = 5) +
  scale_fill_manual(
    values = c(
      "TD"    = "white",
      "ASD-L" = "#86b5a1",
      "ASD-H" = "#f9ae78"
    )
  ) +
  labs(
    x = "连接步数",
    y = "zSFEI",
    fill = "组别"
  ) +
  theme_bw(base_size = 22) +
  theme(
    strip.text = element_text(size = 30, face = "bold", margin = margin(t = 2, b = 2)),
    panel.grid   = element_blank(),
    plot.title   = element_blank(),
    axis.text.y  = element_text(size = 25),
    axis.text.x  = element_text(size = 25),
    axis.title   = element_text(size = 30),
    legend.title = element_text(size = 35),
    legend.text  = element_text(size = 35),
    legend.position = "top"
  )

p3

ggsave(
  filename = file.path(plot_dir, "zSFEI_network_step_abide.png"),
  plot     = p3,
  width    = 3000,
  height   = 2000,
  dpi      = 300,
  units    = "px"
)

# ============================================================
# Network × Step mean zSFEI (for heatmap)
# ============================================================

mean_mat <- data_all %>%
  group_by(Group, Network, Step) %>%
  summarise(
    mean_zSFEI = mean(zSFEI),
    .groups = "drop"
  )

mean_mat_plot <- mean_mat %>%
  left_join(network_map, by = "Network") %>%
  filter(!is.na(NetworkLabel)) %>%
  mutate(
    NetworkLabel = factor(
      NetworkLabel,
      levels = network_map$NetworkLabel
    ),
    Step = factor(Step)
  )

p_heatmap <- ggplot(
  mean_mat_plot,
  aes(
    x = Step,
    y = NetworkLabel,
    fill = mean_zSFEI
  )
) +
  geom_tile(
    color = "lightgray",
    linewidth = 0.4
  ) +
  facet_wrap(
    ~ Group,
    nrow = 1
  ) +
  scale_fill_gradient2(
    low      = "#3d5c6f",
    mid      = "white",
    high     = "#e47159",
    midpoint = 0,
    limits   = c(-0.5, 0.5),
    oob      = scales::squish,
    name     = "zSFEI均值"
  ) +
  labs(
    x     = "SFC步数",
    y     = "功能网络"
  ) +
  theme_bw(base_size = 22) +
  theme(
    # plot.title   = element_text(size = 40, face = "bold", hjust = 0.5),
    plot.title   = element_blank(),
    strip.text   = element_text(size = 32, face = "bold"),
    axis.text.y  = element_text(size = 30),
    axis.text.x  = element_text(size = 30),
    axis.title   = element_text(size = 35),
    legend.title = element_text(size = 35),
    legend.text  = element_text(size = 25),
    panel.grid   = element_blank()
  )


ggsave(
  filename = file.path(plot_dir, "zSFEI_heatmap_abide.png"),
  plot     = p_heatmap,
  width    = 3200,
  height   = 2200,
  dpi      = 300,
  units    = "px"
)

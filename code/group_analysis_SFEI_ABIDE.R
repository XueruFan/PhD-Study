
rm(list = ls())

# ----------------------------
# 0. Libraries
# ----------------------------
library(tidyverse)
library(readxl)
library(showtext)
library(viridis)

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
embedding_dir <- "/Volumes/Zuolab_XRF/output/abide/sfc/sfc_nbtw_embedding"
demo_path     <- "/Volumes/Zuolab_XRF/supplement/abide/abide_demo.xlsx"
site_dir      <- "/Volumes/Zuolab_XRF/data/abide/sublist"
cluster_path  <- "/Volumes/Zuolab_XRF/output/abide/ABIDE_cluster_all_subjects.csv"

out_root <- "/Volumes/Zuolab_XRF/output/abide/sfc/des"
plot_dir <- "/Volumes/Zuolab_XRF/output/abide/sfc/plot"

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

step_max <- 7
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
      values_to = "Embedding"
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
    Embedding = as.numeric(Embedding)
  )

participant_for_analysis <- data.frame(unique(data_all$Subject))
write_csv(participant_for_analysis, file.path(out_root, "sfc_participant_for_analysis.csv"),
          col_names = F)

# ============================================================
# 7. Descriptive statistics (Table 1 by SUBTYPE)
# ============================================================
demo_table <- data_all %>%
  distinct(Subject, Group, AGE_AT_SCAN, site) %>%
  group_by(Group) %>%
  summarise(
    N        = n(),
    Age_mean = mean(AGE_AT_SCAN),
    Age_sd   = sd(AGE_AT_SCAN),
    .groups  = "drop"
  )

write_csv(demo_table, file.path(out_root, "sfc_demo.csv"))

# ============================================================
# 8. Network-wise descriptive stats
# ============================================================
desc_network <- data_all %>%
  group_by(Group, Network) %>%
  summarise(
    N                = n(),
    mean_embedding   = mean(Embedding),
    sd_embedding     = sd(Embedding),
    median_embedding = median(Embedding),
    iqr_embedding    = IQR(Embedding),
    .groups = "drop"
  )

write_csv(desc_network, file.path(out_root, "SFEI_network.csv"))

# ============================================================
# 9. Step-wise descriptive stats
# ============================================================
desc_step <- data_all %>%
  group_by(Group, Step) %>%
  summarise(
    N              = n(),
    mean_embedding = mean(Embedding),
    sd_embedding   = sd(Embedding),
    .groups = "drop"
  )

write_csv(desc_step, file.path(out_root, "SFEI_step.csv"))

# ============================================================
# 9.5 Network × Step mean embedding (for heatmap)
# ============================================================

mean_mat <- data_all %>%
  group_by(Group, Network, Step) %>%
  summarise(
    mean_embedding = mean(Embedding),
    .groups = "drop"
  )

# ============================================================
# 10. Key descriptive visualizations (High-res PNG)
# ============================================================

# ## 10.1 Embedding distribution
# p1 <- ggplot(data_all, aes(x = Embedding, fill = Group)) +
#   geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
#   facet_wrap(~ Step) +
#   labs(
#     title = "不同亚型在各步骤的 Embedding 分布",
#     x = "Embedding 数值",
#     y = "被试数量",
#     fill = "亚型"
#   )
# 
# ggsave(
#   filename = file.path(plot_dir, "sfc_distribution.png"),
#   plot     = p1,
#   width    = 3000,
#   height   = 2000,
#   dpi      = 300,
#   units    = "px"
# )

# ## 10.2 Network-wise boxplot
# p2 <- ggplot(data_all, aes(x = Network, y = Embedding, fill = Group)) +
#   geom_boxplot(outlier.alpha = 0.3) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(
#     title = "不同亚型的网络级 Embedding 分布",
#     x = "功能网络",
#     y = "Embedding 数值",
#     fill = "亚型"
#   )
# 
# ggsave(
#   filename = file.path(plot_dir, "sfc_network.png"),
#   plot     = p2,
#   width    = 3000,
#   height   = 2000,
#   dpi      = 300,
#   units    = "px"
# )


# ============================================================
## 10.3 Network × Step heatmap
# ============================================================
library(scico)

network_map <- tibble(
  Network = c(
    "Net13", "Net01", 
    "Net04", "Net08",
    "Net05",  
    "Net02", 
    "Net07", "Net12",
    "Net06",
    "Net14", 
    "Net10", "Net11",
    "Net03", "Net15", 
    "Net09"
  ),
  NetworkLabel = c(
    "VIS-C", "VIS-P",
    "SMOT-B", "SMOT-A",
    "AUD",
    "AN",
    "dATN-B", "dATN-A",
    "PM-PPr", 
    "SAL/PMN",
    "FPN-B", "FPN-A",
    "DN-B", "DN-A",
    "LANG"
  ),
  Order = 1:15
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
    fill = mean_embedding
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
  scale_fill_gradient(
    low  = "white",
    high = "#ff6666",
    limits = range(mean_mat_plot$mean_embedding),
    name = "SFEI均值"
  ) +
  # scale_fill_scico(
  #   palette = "roma",
  #   direction = -1,
  #   limits = range(mean_mat_plot$mean_embedding),
  #   name = "SFEI均值"
  # ) +
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
    legend.text  = element_text(size = 30),
    panel.grid   = element_blank()
  )


ggsave(
  filename = file.path(plot_dir, "SFEI_heatmap_abide.png"),
  plot     = p_heatmap,
  width    = 3200,
  height   = 2200,
  dpi      = 300,
  units    = "px"
)

# ============================================================
# 11. Build final participant summary table for analysis
# ============================================================

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
  file.path(out_root, "sfc_participant_summary.csv")
)

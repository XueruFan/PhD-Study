############################################################
# Network-level Centile Comparison (No Step)
############################################################

rm(list = ls())

library(tidyverse)
library(readxl)

############################################################
# 路径
############################################################

data_path <- "/Volumes/Zuolab_XRF/output/norm_zSFEI/centile/ASD_centile_results.xlsx"
plot_dir  <- "/Volumes/Zuolab_XRF/output/norm_zSFEI/centile/plots"

dir.create(plot_dir, showWarnings = FALSE)

############################################################
# 读取数据
############################################################

data_all <- read_excel(data_path)

data_all <- data_all %>%
  mutate(
    Network = as.factor(Network),
    Subtype = as.factor(Subtype),
    centile = as.numeric(centile)
  )

############################################################
# 网络标签顺序
############################################################

network_map <- tibble(
  Network = c(
    "Net05","Net13","Net01",
    "Net04","Net08",
    "Net14","Net06","Net02",
    "Net07","Net12",
    "Net10","Net11",
    "Net03","Net15","Net09"
  ),
  NetworkLabel = c(
    "AUD",
    "VIS-C","VIS-P",
    "SMOT-B","SMOT-A",
    "SAL/PMN",
    "PM-PPr",
    "AN",
    "dATN-B","dATN-A",
    "FPN-B","FPN-A",
    "DN-B","DN-A",
    "LANG"
  )
)

data_plot <- data_all %>%
  left_join(network_map, by = "Network")
data_plot$Subtype <- factor(
  data_plot$Subtype,
  levels = c("ASD-L", "ASD-H")
)

data_plot$NetworkLabel <- factor(
  data_plot$NetworkLabel,
  levels = network_map$NetworkLabel
)

############################################################
# 画图
############################################################

p_network <- ggplot(
  data_plot,
  aes(x = NetworkLabel,
      y = centile,
      fill = Subtype)
) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.6) +
  geom_hline(yintercept = c(0.05, 0.95), linetype = 3, linewidth = 0.4) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    linewidth = 0.4,
    position = position_dodge(width = 0.7)
  ) +
  scale_fill_manual(
    values = c(
      "TD"    = "white",
      "ASD-L" = "#86b5a1",
      "ASD-H" = "#f9ae78"
    )
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "功能网络",
    y = "zSFEI的百分位位置",
    fill = "亚型"
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title  = element_text(size = 12),
    # axis.title  = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 12)
  )

p_network

ggsave(
  filename = file.path(plot_dir, "研究三_两个ASD亚型的SFEI百分位位置分布.png"),
  plot     = p_network,
  width    = 3200,
  height   = 2000,
  dpi      = 300,
  units    = "px"
)


############################################################
# Centile Mean Heatmap (Subtype × Network × Step)
############################################################

mean_mat <- data_all %>%
  filter(Subtype %in% c("ASD-L","ASD-H")) %>%
  group_by(Subtype, Network, Step) %>%
  summarise(
    mean_centile = mean(centile, na.rm = TRUE),
    .groups = "drop"
  )

mean_mat_plot <- mean_mat %>%
  left_join(network_map, by = "Network") %>%
  filter(!is.na(NetworkLabel)) %>%
  mutate(
    Subtype = factor(Subtype,
                     levels = c("ASD-L","ASD-H")),
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
    fill = mean_centile
  )
) +
  geom_tile(
    color = "lightgray",
    linewidth = 0.4
  ) +
  facet_wrap(
    ~ Subtype,
    nrow = 1
  ) +
  scale_fill_gradient2(
    low      = "#3d5c6f",
    mid      = "white",
    high     = "#e47159",
    midpoint = 0.5,
    limits   = c(0, 1),
    oob      = scales::squish,
    name     = "zSFEI的Centile均值"
  ) +
  labs(
    x     = "SFC步数",
    y     = "功能网络"
  ) +
  theme_classic(base_size = 22) +
  theme(
    strip.text   = element_text(size = 32, face = "bold"),
    axis.text.y  = element_text(size = 30),
    axis.text.x  = element_text(size = 30),
    axis.title   = element_text(size = 35),
    legend.title = element_text(size = 35),
    legend.text  = element_text(size = 25),
    panel.grid   = element_blank()
  )

ggsave(
  filename = file.path(plot_dir, "Centile_heatmap_subtype.png"),
  plot     = p_heatmap,
  width    = 2500,
  height   = 2200,
  dpi      = 300,
  units    = "px"
)

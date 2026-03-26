rm(list = ls())

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(broom.mixed)

set.seed(1205)

rdcm <- read_xlsx("/Volumes/Zuolab_XRF/output/norm_rDCM/rDCM_normative_data.xlsx") %>%
  filter(Cohort == "ABIDE") %>%
  mutate(subject = as.character(as.numeric(ID)),
    Site = recode(Site, "NYU2" = "NYU"),
    Subtype = factor(Subtype, levels = c("TD","ASD-L","ASD-H"))
  ) %>%
  dplyr::select(-Session, -Cohort, -subject)

data_long <- rdcm %>%
  pivot_longer(
    cols = starts_with("EC_"),
    names_to = "Edge",
    values_to = "EC"
  )

# 8. Model definition
fit_ec <- function(df) {
  m <- lmer(EC ~ Subtype + Age + (1 | Site), data = df, REML = FALSE)  # 这里加了 Age
  em <- emmeans::emmeans(m, pairwise ~ Subtype, adjust = "fdr")
  broom::tidy(em$contrasts)
}

# 9. Fit model for all connections
results <- data_long %>% 
  group_by(Edge) %>% 
  group_modify(~fit_ec(.x)) %>% 
  ungroup()

# 5. 筛选显著结果 adj.p.value < 0.05
significant_results <- results %>%
  filter(adj.p.value < 0.05) %>%
  arrange(adj.p.value) %>%
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 3),
    adj.p.value = round(adj.p.value, 4)
  ) %>%
  select(
    脑边 = Edge,
    比较组 = contrast,
    效应值 = estimate,
    标准误 = std.error,
    统计量 = statistic,
    校正P值 = adj.p.value
  )

# 6. 网络编号 → 英文缩写映射（完全按你的图）
Index <- c(
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

Abbreviation <- c(
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

net_map <- tibble(Index, Abbreviation)

# 7. 给显著表加上 源网络缩写 + 目标网络缩写
significant_final <- significant_results %>%
  mutate(
    src_idx = as.integer(str_extract(脑边, "(?<=EC_)\\d+")),
    tgt_idx = as.integer(str_extract(脑边, "(?<=to_)\\d+"))
  ) %>%
  left_join(net_map, by = c("src_idx" = "Index")) %>% rename(源网络 = Abbreviation) %>%
  left_join(net_map, by = c("tgt_idx" = "Index")) %>% rename(目标网络 = Abbreviation) %>%
  select(
    脑边,
    源网络,
    目标网络,
    比较组,
    效应值,
    标准误,
    统计量,
    校正P值
  )

writexl::write_xlsx(significant_results, "/Volumes/Zuolab_XRF/output/abide/dcm/stat/ABIDE3组人原始EC的显著差异边.xlsx")


############################################################
# ===== 原始EC + 差异热图（带显著性标记） =====
############################################################

library(ggplot2)
library(dplyr)
library(stringr)

##########
# 1. 计算三组 EC 均值
ec_group_mean <- data_long %>%
  group_by(Subtype, Edge) %>%
  summarise(EC_mean = mean(EC, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    From = as.integer(str_extract(Edge,"(?<=EC_)\\d+")),
    To   = as.integer(str_extract(Edge,"(?<=to_)\\d+"))
  ) %>%
  left_join(net_map, by = c("From" = "Index")) %>% rename(FromNetwork = Abbreviation) %>%
  left_join(net_map, by = c("To" = "Index")) %>% rename(ToNetwork = Abbreviation) %>%
  mutate(
    EC_mean = ifelse(From == To, NA, EC_mean)
  )

##########
# 2. 统一网络顺序（关键：和你之前完全一致）
network_order <- net_map$Abbreviation

ec_group_mean$FromNetwork <- factor(
  ec_group_mean$FromNetwork,
  levels = rev(network_order)
)

ec_group_mean$ToNetwork <- factor(
  ec_group_mean$ToNetwork,
  levels = network_order
)

##########
# 3. 统一颜色范围（所有图一致）
ec_range <- range(ec_group_mean$EC_mean, na.rm = TRUE)

##########
# 4. 画三组原始 EC 热图
make_ec_heatmap <- function(subtype_name){
  
  df_plot <- ec_group_mean %>%
    filter(Subtype == subtype_name)
  
  p <- ggplot(df_plot,
              aes(x = ToNetwork, y = FromNetwork, fill = EC_mean)) +
    geom_tile(color="lightgray", linewidth=0.3) +
    scale_fill_gradient2(
      low = "#3d5c6f",
      mid = "white",
      high = "#e47159",
      midpoint = 0,
      limits = ec_range
    ) +
    theme_bw(base_size = 22) +
    labs(fill = "EC") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    paste0("/Volumes/Zuolab_XRF/output/abide/dcm/plot/EC_mean_", subtype_name, ".png"),
    p, width = 3200, height = 2800, dpi = 300, units = "px"
  )
}

make_ec_heatmap("TD")
make_ec_heatmap("ASD-L")
make_ec_heatmap("ASD-H")

############################################################
# ===== 差值热图 + 显著性圈标记（最终版） =====
############################################################

##########
# 1. 转宽格式
ec_wide <- ec_group_mean %>%
  select(Subtype, Edge, EC_mean, FromNetwork, ToNetwork) %>%
  pivot_wider(names_from = Subtype, values_from = EC_mean)

##########
# 2. 计算差值
ec_diff <- ec_wide %>%
  mutate(
    diff_H_vs_TD = `ASD-H` - TD,
    diff_L_vs_TD = `ASD-L` - TD,
    diff_H_vs_L  = `ASD-H` - `ASD-L`
  )

##########
# 3. 去掉自连接（关键）
ec_diff <- ec_diff %>%
  mutate(
    diff_H_vs_TD = ifelse(FromNetwork == ToNetwork, NA, diff_H_vs_TD),
    diff_L_vs_TD = ifelse(FromNetwork == ToNetwork, NA, diff_L_vs_TD),
    diff_H_vs_L  = ifelse(FromNetwork == ToNetwork, NA, diff_H_vs_L)
  )

##########
# 4. 处理显著性（你的 sig_edges）
sig_edges_clean <- sig_edges %>%
  mutate(
    # 统一去括号
    比较组 = str_replace_all(比较组, "[()]", ""),
    
    # 转换为标准方向
    comp_std = case_when(
      比较组 == "TD - ASD-H"  ~ "ASD-H_vs_TD",
      比较组 == "ASD-H - TD"  ~ "ASD-H_vs_TD",
      比较组 == "TD - ASD-L"  ~ "ASD-L_vs_TD",
      比较组 == "ASD-L - TD"  ~ "ASD-L_vs_TD",
      比较组 == "ASD-H - ASD-L" ~ "ASD-H_vs_L",
      比较组 == "ASD-L - ASD-H" ~ "ASD-H_vs_L",
      TRUE ~ NA_character_
    )
  )

##########
# 5. 拆成三个比较
sig_H_vs_TD <- sig_edges_clean %>% filter(comp_std == "ASD-H_vs_TD")
sig_L_vs_TD <- sig_edges_clean %>% filter(comp_std == "ASD-L_vs_TD")
sig_H_vs_L  <- sig_edges_clean %>% filter(comp_std == "ASD-H_vs_L")

##########
# 6. 网络顺序（完全一致）
network_order <- net_map$Abbreviation

##########
# 7. 画差值热图（用圈标记！）
make_diff_heatmap <- function(diff_col, sig_df, name){
  
  df_plot <- ec_diff %>%
    select(Edge, FromNetwork, ToNetwork, !!sym(diff_col)) %>%
    rename(Diff = !!sym(diff_col)) %>%
    mutate(sig = ifelse(Edge %in% sig_df$Edge, 1, 0))
  
  df_plot$FromNetwork <- factor(df_plot$FromNetwork, levels = rev(network_order))
  df_plot$ToNetwork   <- factor(df_plot$ToNetwork, levels = network_order)
  
  diff_range <- max(abs(df_plot$Diff), na.rm = TRUE)
  
  p <- ggplot(df_plot,
              aes(x = ToNetwork, y = FromNetwork, fill = Diff)) +
    
    geom_tile(color="lightgray", linewidth=0.3) +
    
    # 🔴 显著性：用圈！
    geom_point(
      data = df_plot %>% filter(sig == 1 & !is.na(Diff)),
      aes(x = ToNetwork, y = FromNetwork),
      shape = 21,
      size = 5,
      stroke = 1.2,
      color = "black",
      fill = NA
    ) +
    
    scale_fill_gradient2(
      low = "#3d5c6f",
      mid = "white",
      high = "#e47159",
      midpoint = 0,
      limits = c(-diff_range, diff_range)
    ) +
    
    theme_bw(base_size = 22) +
    labs(fill = expression(Delta~EC)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    paste0("/Volumes/Zuolab_XRF/output/abide/dcm/plot/", name, ".png"),
    p,
    width = 3200,
    height = 2800,
    dpi = 300,
    units = "px"
  )
}

##########
# 8. 生成三张图
make_diff_heatmap("diff_H_vs_TD", sig_H_vs_TD, "EC_diff_ASD-H_vs_TD")
make_diff_heatmap("diff_L_vs_TD", sig_L_vs_TD, "EC_diff_ASD-L_vs_TD")
make_diff_heatmap("diff_H_vs_L",  sig_H_vs_L,  "EC_diff_ASD-H_vs_ASD-L")

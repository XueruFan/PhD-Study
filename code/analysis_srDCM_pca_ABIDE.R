rm(list = ls())

library(tidyverse)
library(readxl)
library(openxlsx)
library(neuroCombat)
library(rrcov)
library(effectsize)
library(emmeans)
library(showtext)
library(reshape2)
library(broom)
library(dplyr)
library(stringr)
library(ggplot2)

############################################################
font_add(
  family = "pingfang",
  regular = "/System/Library/Fonts/PingFang.ttc"
)
showtext_auto()
theme_set(theme_bw(base_family = "pingfang"))

############################################################
out_root <- "/Volumes/Zuolab_XRF/output/abide/dcm/des"
plot_dir <- "/Volumes/Zuolab_XRF/output/abide/dcm/plot"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

############################################################

dcm_path  <- "/Volumes/Zuolab_XRF/output/abide/dcm/sum/ABIDE_srDCM_summary.xlsx"
demo_path <- "/Volumes/Zuolab_XRF/output/abide/sfc/des/zSFEI_abide_demo.csv"

dcm_df <- read_excel(dcm_path) %>%
  mutate(subject = as.character(as.numeric(subject)))

demo_df <- read.csv(demo_path) %>%
  mutate(subject = as.character(as.numeric(Subject))) %>%
  dplyr::select(-Subject)

df <- left_join(demo_df, dcm_df, by = "subject") %>%
  mutate(
    site = recode(site, "NYU2" = "NYU"),
    Subtype = factor(Subtype, levels = c("TD","ASD-L","ASD-H"))
  )

############################################################
# 提取有效连接

ec_cols <- grep("^EC_", colnames(df), value = TRUE)
ec_matrix <- as.matrix(df[, ec_cols])

############################################################
# ComBat

combat_input <- t(ec_matrix)
batch <- as.factor(df$site)
mod <- model.matrix(~ Age + Subtype, data = df)

combat_res <- neuroCombat(
  dat = combat_input,
  batch = batch,
  mod = mod,
  parametric = TRUE,
  eb = TRUE,
  mean.only = TRUE
)

ec_matrix_combat <- t(combat_res$dat.combat)

############################################################
# 删除零方差边

edge_variance <- apply(ec_matrix_combat, 2, var)
ec_matrix_combat <- ec_matrix_combat[, edge_variance > 0]

############################################################
# 标准化

ec_matrix_scaled <- scale(ec_matrix_combat)

############################################################
# Robust PCA

rpca_res <- PcaHubert(
  ec_matrix_scaled,
  k = 5,
  scale = FALSE
)

############################################################
# 提取得分

scores_rpca <- as.data.frame(rpca_res@scores)
colnames(scores_rpca)[1:5] <- paste0("PC",1:5)

scores_rpca <- scores_rpca %>%
  mutate(
    subject = df$subject,
    Subtype = df$Subtype,
    Age     = df$Age,
    site    = df$site
  )

############################################################
# 提取载荷

loadings_rpca <- as.data.frame(rpca_res@loadings)
colnames(loadings_rpca)[1:5] <- paste0("PC",1:5)
loadings_rpca$edge <- colnames(ec_matrix_scaled)

############################################################

lm_pc1 <- lm(PC1 ~ Subtype + Age + site, data = scores_rpca)

summary_df <- broom::tidy(lm_pc1)
anova_df   <- broom::tidy(anova(lm_pc1))
eta_df     <- as.data.frame(eta_squared(lm_pc1))

emm_res <- emmeans(lm_pc1, pairwise ~ Subtype)
emm_df  <- as.data.frame(emm_res$emmeans)
contrast_df <- as.data.frame(emm_res$contrasts)

############################################################

wb <- createWorkbook()

addWorksheet(wb,"Scores")
writeData(wb,"Scores",scores_rpca)

addWorksheet(wb,"Loadings")
writeData(wb,"Loadings",loadings_rpca)

addWorksheet(wb,"LM_Summary")
writeData(wb,"LM_Summary",summary_df)

addWorksheet(wb,"ANOVA")
writeData(wb,"ANOVA",anova_df)

addWorksheet(wb,"Effect_Size")
writeData(wb,"Effect_Size",eta_df)

addWorksheet(wb,"EMMEANS")
writeData(wb,"EMMEANS",emm_df)

addWorksheet(wb,"Contrasts")
writeData(wb,"Contrasts",contrast_df)

saveWorkbook(
  wb,
  file.path(out_root,"ABIDE_male_RobustPCA_full_results.xlsx"),
  overwrite = TRUE
)

############################################################
# ================= 图 1 =================
# 机制轴 PC1 分布
############################################################

p_axis <- ggplot(scores_rpca,
                 aes(x = Subtype, y = PC1, fill = Subtype)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  scale_fill_manual(
    values = c(
      "TD"    = "white",
      "ASD-L" = "#86b5a1",
      "ASD-H" = "#f9ae78"
    )
  ) +
  theme_classic(base_size = 28) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text  = element_text(size = 20),
  )
p_axis

ggsave(
  file.path(plot_dir,"RobustPCA_PC1_axis.png"),
  p_axis,
  width = 3000,
  height = 2400,
  dpi = 300,
  units = "px"
)

############################################################
# ================= 图 2 =================
# PC1-PC2 空间
############################################################

p_space <- ggplot(scores_rpca,
                  aes(x = PC1, y = PC2, color = Subtype)) +
  geom_point(size = 4, alpha = 0.85) +
  scale_color_manual(
    values = c(
      "TD"    = "black",
      "ASD-L" = "#86b5a1",
      "ASD-H" = "#f9ae78"
    )
  ) +
  labs(
    x = "PC1",
    y = "PC2"
  ) +
  theme_classic(base_size = 28)
p_space

ggsave(
  file.path(plot_dir,"RobustPCA_space.png"),
  p_space,
  width = 3000,
  height = 2400,
  dpi = 300,
  units = "px"
)

############################################################
# 严格构建 15×15 PC1 载荷矩阵（避免错位）
############################################################

# ---------------------------------------------------------
# 1. 解析 From / To
# ---------------------------------------------------------

load_df <- loadings_rpca %>%
  mutate(
    From = as.integer(str_extract(edge, "(?<=EC_)\\d+")),
    To   = as.integer(str_extract(edge, "(?<=to_)\\d+"))
  ) %>%
  dplyr::select(From, To, PC1)

# ---------------------------------------------------------
# 2. 建立完整 15×15 结构（防止缺边导致错位）
# ---------------------------------------------------------

full_grid <- expand.grid(
  From = 1:15,
  To   = 1:15
)

load_df_full <- full_grid %>%
  left_join(load_df, by = c("From","To"))

# ---------------------------------------------------------
# 3. 加入网络顺序映射
# ---------------------------------------------------------

network_map <- tibble(
  Order = 1:15,
  Network = c(
    "Net01","Net02","Net03","Net04","Net05",
    "Net06","Net07","Net08","Net09","Net10",
    "Net11","Net12","Net13","Net14","Net15"
  )
)

# 这里改成你论文顺序
network_order_custom <- c(
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

network_labels_custom <- c(
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

network_map_custom <- tibble(
  Order = network_order_custom,
  Label = network_labels_custom
)

# ---------------------------------------------------------
# 4. 映射标签
# ---------------------------------------------------------

load_df_full <- load_df_full %>%
  left_join(network_map_custom, by = c("From" = "Order")) %>%
  rename(FromLabel = Label) %>%
  left_join(network_map_custom, by = c("To" = "Order")) %>%
  rename(ToLabel = Label)

# ---------------------------------------------------------
# 5. 设置因子顺序（关键）
# ---------------------------------------------------------

load_df_full$FromLabel <- factor(
  load_df_full$FromLabel,
  levels = network_labels_custom
)

load_df_full$ToLabel <- factor(
  load_df_full$ToLabel,
  levels = network_labels_custom
)

# ---------------------------------------------------------
# 6. 绘图
# ---------------------------------------------------------

p_heat <- ggplot(
  load_df_full,
  aes(x = ToLabel, y = FromLabel, fill = PC1)
) +
  geom_tile(color="lightgray", linewidth=0.4) +
  scale_fill_gradient2(
    low="#3d5c6f",
    mid="white",
    high="#e47159",
    midpoint=0,
    name="PC1载荷"
  ) +
  labs(
    x="输入网络",
    y="输出网络"
  ) +
  theme_bw(base_size = 28) +
  theme(
    axis.text.x  = element_text(angle=45, hjust=1, size=26),
    axis.text.y  = element_text(size=26),
    axis.title   = element_text(size=32),
    legend.title = element_text(size=30),
    legend.text  = element_text(size=26),
    panel.grid   = element_blank()
  )

ggsave(
  file.path(plot_dir,"RobustPCA_PC1_loadings_heatmap_named.png"),
  p_heat,
  width=3400,
  height=3000,
  dpi=300,
  units="px"
)

################################################################################################
# ================= 自动识别关键网络 + 批量验证 =================
################################################################################################

# 网络英文名称映射（与你热图一致）

network_map_custom <- tibble(
  Order = network_order_custom,
  Network = network_labels_custom
)

############################################################
# 1. 计算 PC1 节点贡献
############################################################

node_contribution <- load_df %>%
  group_by(To) %>%
  summarise(
    In_abs = sum(abs(PC1), na.rm = TRUE)
  ) %>%
  rename(Node = To) %>%
  left_join(
    load_df %>%
      group_by(From) %>%
      summarise(
        Out_abs = sum(abs(PC1), na.rm = TRUE)
      ) %>%
      rename(Node = From),
    by = "Node"
  ) %>%
  mutate(
    In_z  = as.numeric(scale(In_abs)),
    Out_z = as.numeric(scale(Out_abs))
  ) %>%
  left_join(network_map_custom, by = c("Node" = "Order"))

############################################################
# 2. 自动筛选关键网络
############################################################

key_nodes <- node_contribution %>%
  filter(abs(In_z) > 2 | abs(Out_z) > 2)

############################################################
# 3. 批量验证函数
############################################################

validate_node <- function(node_id, node_name) {
  
  edge_info <- tibble(edge = colnames(ec_matrix_scaled)) %>%
    mutate(
      From = as.integer(str_extract(edge, "(?<=EC_)\\d+")),
      To   = as.integer(str_extract(edge, "(?<=to_)\\d+"))
    )
  
  node_edges <- edge_info %>%
    filter(To == node_id) %>%
    pull(edge)
  
  node_strength <- ec_matrix_scaled[, node_edges] %>%
    as.data.frame() %>%
    mutate(Node_in_mean = rowMeans(.))
  
  data_tmp <- scores_rpca %>%
    bind_cols(node_strength)
  
  cor_res <- cor.test(data_tmp$PC1, data_tmp$Node_in_mean)
  
  cor_df <- data.frame(
    Network = node_name,
    r       = cor_res$estimate,
    p_value = cor_res$p.value
  )
  
  # 图
  p_cor <- ggplot(data_tmp,
                  aes(x = Node_in_mean,
                      y = PC1,
                      color = Subtype)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    theme_classic(base_size = 24) +
    labs(
      x = paste0(node_name, " Incoming Connectivity"),
      y = "PC1"
    )
  
  ggsave(
    file.path(plot_dir,
              paste0("PC1_vs_", node_name, "_incoming.png")),
    p_cor,
    width = 3000,
    height = 2400,
    dpi = 300,
    units = "px"
  )
  
  return(cor_df)
}

############################################################
# 4. 运行批量验证
############################################################

all_cor_results <- mapply(
  validate_node,
  node_id   = key_nodes$Node,
  node_name = key_nodes$Network,
  SIMPLIFY = FALSE
)

all_cor_results_df <- do.call(rbind, all_cor_results)

############################################################
# 5. 保存汇总文件
############################################################

wb_auto <- createWorkbook()

addWorksheet(wb_auto,"Node_Contribution")
writeData(wb_auto,"Node_Contribution",node_contribution)

addWorksheet(wb_auto,"Key_Networks")
writeData(wb_auto,"Key_Networks",key_nodes)

addWorksheet(wb_auto,"All_Correlations")
writeData(wb_auto,"All_Correlations",all_cor_results_df)

saveWorkbook(
  wb_auto,
  file.path(out_root,"PC1_key_networks_summary.xlsx"),
  overwrite = TRUE
)

################################################################################################
# ================= 自动 PC2 机制轴识别 =================
################################################################################################

############################################################
# 1. 提取 PC2 载荷
############################################################

load_df_pc2 <- loadings_rpca %>%
  mutate(
    From = as.integer(str_extract(edge, "(?<=EC_)\\d+")),
    To   = as.integer(str_extract(edge, "(?<=to_)\\d+"))
  ) %>%
  dplyr::select(From, To, PC2)

############################################################
# 2. 计算节点贡献（PC2）
############################################################

node_contribution_pc2 <- load_df_pc2 %>%
  group_by(To) %>%
  summarise(
    In_abs = sum(abs(PC2), na.rm = TRUE)
  ) %>%
  rename(Node = To) %>%
  left_join(
    load_df_pc2 %>%
      group_by(From) %>%
      summarise(
        Out_abs = sum(abs(PC2), na.rm = TRUE)
      ) %>%
      rename(Node = From),
    by = "Node"
  ) %>%
  mutate(
    In_z  = as.numeric(scale(In_abs)),
    Out_z = as.numeric(scale(Out_abs))
  ) %>%
  left_join(network_map_custom, by = c("Node" = "Order"))

############################################################
# 3. 自动筛选关键网络（PC2）
############################################################

key_nodes_pc2 <- node_contribution_pc2 %>%
  filter(abs(In_z) > 2 | abs(Out_z) > 2)

############################################################
# 4. PC2 批量验证函数
############################################################

validate_node_pc2 <- function(node_id, node_name) {
  
  edge_info <- tibble(edge = colnames(ec_matrix_scaled)) %>%
    mutate(
      From = as.integer(str_extract(edge, "(?<=EC_)\\d+")),
      To   = as.integer(str_extract(edge, "(?<=to_)\\d+"))
    )
  
  node_edges <- edge_info %>%
    filter(To == node_id) %>%
    pull(edge)
  
  node_strength <- ec_matrix_scaled[, node_edges] %>%
    as.data.frame() %>%
    mutate(Node_in_mean = rowMeans(.))
  
  data_tmp <- scores_rpca %>%
    bind_cols(node_strength)
  
  cor_res <- cor.test(data_tmp$PC2, data_tmp$Node_in_mean)
  
  cor_df <- data.frame(
    Network = node_name,
    r       = cor_res$estimate,
    p_value = cor_res$p.value
  )
  
  ############################################################
  # 图
  ############################################################
  
  p_cor <- ggplot(data_tmp,
                  aes(x = Node_in_mean,
                      y = PC2,
                      color = Subtype)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    theme_classic(base_size = 24) +
    labs(
      x = paste0(node_name, " Incoming Connectivity"),
      y = "PC2"
    )
  
  ggsave(
    file.path(plot_dir,
              paste0("PC2_vs_", node_name, "_incoming.png")),
    p_cor,
    width = 3000,
    height = 2400,
    dpi = 300,
    units = "px"
  )
  
  return(cor_df)
}

############################################################
# 5. 运行 PC2 批量验证
############################################################

all_cor_results_pc2 <- mapply(
  validate_node_pc2,
  node_id   = key_nodes_pc2$Node,
  node_name = key_nodes_pc2$Network,
  SIMPLIFY = FALSE
)

all_cor_results_pc2_df <- do.call(rbind, all_cor_results_pc2)

lm_pc2 <- lm(PC2 ~ Subtype + Age + site, data = scores_rpca)

summary_lm_pc2 <- broom::tidy(lm_pc2)
anova_lm_pc2 <- broom::tidy(anova(lm_pc2))

############################################################
# 6. 保存 PC2 汇总文件
############################################################

wb_auto_pc2 <- createWorkbook()

addWorksheet(wb_auto_pc2,"Node_Contribution_PC2")
writeData(wb_auto_pc2,"Node_Contribution_PC2",node_contribution_pc2)

addWorksheet(wb_auto_pc2,"Key_Networks_PC2")
writeData(wb_auto_pc2,"Key_Networks_PC2",key_nodes_pc2)

addWorksheet(wb_auto_pc2,"All_Correlations_PC2")
writeData(wb_auto_pc2,"All_Correlations_PC2",all_cor_results_pc2_df)

addWorksheet(wb_auto_pc2,"LM_Summary_PC2")
writeData(wb_auto_pc2,"LM_Summary_PC2",summary_lm_pc2)

addWorksheet(wb_auto_pc2,"ANOVA_PC2")
writeData(wb_auto_pc2,"ANOVA_PC2",anova_lm_pc2)


saveWorkbook(
  wb_auto_pc2,
  file.path(out_root,"PC2_key_networks_summary.xlsx"),
  overwrite = TRUE
)

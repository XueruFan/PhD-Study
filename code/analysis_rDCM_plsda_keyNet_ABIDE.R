rm(list = ls())
############################################################
# 1️⃣ 加载包
############################################################

library(readxl)
library(openxlsx)
library(dplyr)
library(stringr)

############################################################
# 2️⃣ 读取 PCA 文件
############################################################

file_path <- "/Volumes/Zuolab_XRF/output/rDCM_PLSDA/PLSDA_statistics.xlsx"

loadings_rpca <- read_excel(file_path, sheet = "Loadings")

# 自动获取所在文件夹
out_dir <- dirname(file_path)

############################################################
# 3️⃣ 网络映射（必须与你原分析一致）
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

network_map_custom <- tibble(
  Label = network_labels_custom,
  Network = network_networks_custom
)

############################################################
# 4️⃣ 主函数（自动保存Excel）
############################################################

analyze_network_direction <- function(network_name,
                                      pc_name,
                                      direction = c("out","in")) {
  
  direction <- match.arg(direction)
  
  if (!pc_name %in% colnames(loadings_rpca)) {
    stop("主成分名称不存在，例如 PC1、PC2 ...")
  }
  
  edge_df <- loadings_rpca %>%
    mutate(
      From = as.integer(str_extract(edge, "(?<=EC_)\\d+")),
      To   = as.integer(str_extract(edge, "(?<=to_)\\d+"))
    ) %>%
    dplyr::select(edge, From, To, !!sym(pc_name))
  
  colnames(edge_df)[4] <- "Loading"
  
  edge_df <- edge_df %>%
    left_join(network_map_custom, by = c("From" = "Label")) %>%
    rename(FromNetwork = Network) %>%
    left_join(network_map_custom, by = c("To" = "Label")) %>%
    rename(ToNetwork = Network) %>%
    mutate(
      Loading_z = as.numeric(scale(Loading))
    )
  
  ##########################################################
  # 根据方向筛选
  ##########################################################
  
  if (direction == "out") {
    result_df <- edge_df %>%
      filter(FromNetwork == network_name) %>%
      arrange(desc(abs(Loading_z)))
  } else {
    result_df <- edge_df %>%
      filter(ToNetwork == network_name) %>%
      arrange(desc(abs(Loading_z)))
  }
  
  key_edges <- result_df %>%
    filter(abs(Loading_z) > 2)
  
  ##########################################################
  # 保存 Excel
  ##########################################################
  
  wb <- createWorkbook()
  
  addWorksheet(wb, "All_Edges")
  writeData(wb, "All_Edges", result_df)
  
  addWorksheet(wb, "Key_Edges_Z_gt_2")
  writeData(wb, "Key_Edges_Z_gt_2", key_edges)
  
  out_file <- file.path(
    out_dir,
    paste0("Network_", network_name, "_", pc_name, "_", direction, "_edges.xlsx")
  )
  
  saveWorkbook(wb, out_file, overwrite = TRUE)
  
  cat("=====================================\n")
  cat("已保存文件到:\n")
  cat(out_file, "\n")
  cat("=====================================\n")
  
  return(list(
    All_Edges = result_df,
    Key_Edges = key_edges,
    File_Path = out_file
  ))
}

############################################################
# 5️⃣ 使用示例
############################################################

res <- analyze_network_direction("VIS-C", "Comp1", "out")
res <- analyze_network_direction("dATN-A", "Comp1", "in")
res <- analyze_network_direction("FPN-A", "Comp1", "in")
res <- analyze_network_direction("AUD", "Comp2", "out")
res <- analyze_network_direction("SAL", "Comp2", "in")
res <- analyze_network_direction("FPN-B", "Comp3", "in")
res <- analyze_network_direction("SMOT-A", "Comp3", "out")
res <- analyze_network_direction("DN-A", "Comp3", "out")

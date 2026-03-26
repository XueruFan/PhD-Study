rm(list = ls())

library(tidyverse)
library(readxl)
library(openxlsx)
library(pls)
library(effectsize)
library(emmeans)
library(broom)
library(stringr)
library(ggplot2)

set.seed(1205)

data_path <- "/Volumes/Zuolab_XRF/output/norm_rDCM/gamlss/rDCM_deviation_zscores.xlsx"
out_root  <- "/Volumes/Zuolab_XRF/output/rDCM_PLSDA"
plot_dir  <- "/Volumes/Zuolab_XRF/output/rDCM_PLSDA/plots"

dir.create(out_root, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

data <- read_excel(data_path)

data_asd <- data %>%
  filter(Subtype %in% c("ASD-L","ASD-H")) %>%
  mutate(Subtype = factor(Subtype))

ec_cols <- grep("^EC_", colnames(data_asd), value = TRUE)

X <- as.matrix(data_asd[, ec_cols])
X <- scale(X)

Y <- ifelse(data_asd$Subtype=="ASD-H",1,0)

n_comp <- 3

pls_model <- plsr(
  Y ~ X,
  ncomp = n_comp,
  validation = "none"
)
#######
scores <- as.data.frame(pls_model$scores[,1:n_comp])
colnames(scores) <- paste0("Comp",1:n_comp)

scores <- scores %>%
  mutate(
    subject = data_asd$subject,
    Subtype = data_asd$Subtype
  )

#######
loadings <- as.data.frame(pls_model$loadings[,1:n_comp])
colnames(loadings) <- paste0("Comp",1:n_comp)
loadings$edge <- ec_cols



#######
load_mat <- as.matrix(loadings[,1:n_comp])

scores_classical <- as.matrix(X) %*% load_mat

comp_variances <- apply(scores_classical, 2, var)

total_variance <- sum(apply(X,2,var))

variance_ratio <- comp_variances / total_variance

variance_df <- data.frame(
  Component = paste0("Comp",1:n_comp),
  Variance = comp_variances,
  Variance_Explained = variance_ratio,
  Cumulative = cumsum(variance_ratio)
)


#####
pcs <- paste0("Comp",1:n_comp)

wb <- createWorkbook()

addWorksheet(wb,"Scores")
writeData(wb,"Scores",scores)

addWorksheet(wb,"Loadings")
writeData(wb,"Loadings",loadings)

addWorksheet(wb,"Variance")
writeData(wb,"Variance",variance_df)

for (pc in pcs) {
  
  formula_str <- as.formula(
    paste0(pc," ~ Subtype")
  )
  
  lm_model <- lm(formula_str, data = scores)
  
  summary_df  <- broom::tidy(lm_model)
  anova_df    <- broom::tidy(anova(lm_model))
  eta_df      <- as.data.frame(eta_squared(lm_model))
  
  emm_res     <- emmeans(lm_model, pairwise ~ Subtype)
  
  emm_df      <- as.data.frame(emm_res$emmeans)
  contrast_df <- as.data.frame(emm_res$contrasts)
  
  addWorksheet(wb,paste0(pc,"_LM"))
  writeData(wb,paste0(pc,"_LM"),summary_df)
  
  addWorksheet(wb,paste0(pc,"_ANOVA"))
  writeData(wb,paste0(pc,"_ANOVA"),anova_df)
  
  addWorksheet(wb,paste0(pc,"_EffectSize"))
  writeData(wb,paste0(pc,"_EffectSize"),eta_df)
  
  addWorksheet(wb,paste0(pc,"_EMMEANS"))
  writeData(wb,paste0(pc,"_EMMEANS"),emm_df)
  
  addWorksheet(wb,paste0(pc,"_Contrasts"))
  writeData(wb,paste0(pc,"_Contrasts"),contrast_df)
  
}

saveWorkbook(
  wb,
  file.path(out_root,"PLSDA_statistics.xlsx"),
  overwrite=TRUE
)

########
for (pc in pcs) {
  
  p_axis <- ggplot(scores,
                   aes_string(x="Subtype",y=pc,fill="Subtype")) +
    geom_violin(trim=FALSE,alpha=0.4) +
    geom_boxplot(width=0.15) +
    geom_jitter(width=0.08,alpha=0.6) +
    scale_fill_manual(values=c(
      "ASD-L"="#86b5a1",
      "ASD-H"="#f9ae78"
    )) +
    theme_classic(base_size=28)
  
  ggsave(
    file.path(plot_dir,paste0(pc,"_distribution.png")),
    p_axis,
    width=3000,
    height=2400,
    dpi=300,
    units="px"
  )
  
}

#####
pc_pairs <- combn(pcs,2)

for (i in 1:ncol(pc_pairs)) {
  
  xpc <- pc_pairs[1,i]
  ypc <- pc_pairs[2,i]
  
  p_space <- ggplot(scores,
                    aes_string(x=xpc,y=ypc,color="Subtype")) +
    geom_point(size=4,alpha=0.85) +
    scale_color_manual(values=c(
      "ASD-L"="#86b5a1",
      "ASD-H"="#f9ae78"
    )) +
    theme_classic(base_size=28)
  
  ggsave(
    file.path(plot_dir,
              paste0(xpc,"_",ypc,"_space.png")),
    p_space,
    width=3000,
    height=2400,
    dpi=300,
    units="px"
  )
}

########
network_labels <- c(
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

network_names <- c(
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

network_map <- tibble(
  Label = network_labels,
  Network = network_names
)

make_heatmap <- function(pc_name){
  
  load_df <- loadings %>%
    mutate(
      From = as.integer(str_extract(edge,"(?<=EC_)\\d+")),
      To   = as.integer(str_extract(edge,"(?<=to_)\\d+"))
    ) %>%
    dplyr::select(From,To,!!sym(pc_name))
  
  colnames(load_df)[3] <- "Loading"
  
  load_df_full <- load_df %>%
    left_join(network_map, by = c("From" = "Label")) %>%
    rename(FromNetwork = Network) %>%
    left_join(network_map, by = c("To" = "Label")) %>%
    rename(ToNetwork = Network)
  
  load_df_full$FromNetwork <- factor(
    load_df_full$FromNetwork,
    levels = rev(network_names)
  )
  
  load_df_full$ToNetwork <- factor(
    load_df_full$ToNetwork,
    levels = network_names
  )
  
  p_heat <- ggplot(load_df_full,
                   aes(x=ToNetwork, y=FromNetwork, fill=Loading)) +
    geom_tile(color="lightgray",linewidth=0.4) +
    scale_fill_gradient2(low="#3d5c6f",
                         mid="white",
                         high="#e47159",
                         midpoint=0) +
    theme_bw(base_size=26) +
    theme(axis.text.x=element_text(angle=45,hjust=1))
  
  ggsave(
    file.path(plot_dir,
              paste0(pc_name,"_heatmap.png")),
    p_heat,
    width=3400,
    height=3000,
    dpi=300,
    units="px"
  )
}

make_heatmap("Comp1")
make_heatmap("Comp2")
make_heatmap("Comp3")

##########
compute_node_contribution <- function(comp){
  
  load_df <- loadings %>%
    mutate(
      From = as.integer(str_extract(edge,"(?<=EC_)\\d+")),
      To   = as.integer(str_extract(edge,"(?<=to_)\\d+"))
    ) %>%
    dplyr::select(From,To,!!sym(comp))
  
  colnames(load_df)[3] <- "Loading"
  
  load_df <- load_df %>%
    left_join(network_map,by=c("From"="Label")) %>%
    rename(FromNetwork=Network) %>%
    left_join(network_map,by=c("To"="Label")) %>%
    rename(ToNetwork=Network)
  
  load_df$FromNetwork <- factor(
    load_df$FromNetwork,
    levels = rev(network_names)
  )
  
  load_df$ToNetwork <- factor(
    load_df$ToNetwork,
    levels = network_names
  )
  
  node_df <- load_df %>%
    group_by(ToNetwork) %>%
    summarise(To_abs=sum(abs(Loading),na.rm=TRUE)) %>%
    rename(Node=ToNetwork) %>%
    left_join(
      load_df %>%
        group_by(FromNetwork) %>%
        summarise(From_abs=sum(abs(Loading),na.rm=TRUE)) %>%
        rename(Node=FromNetwork),
      by="Node"
    ) %>%
    mutate(
      To_z=as.numeric(scale(To_abs)),
      From_z=as.numeric(scale(From_abs))
    ) %>%
    left_join(network_map,by=c("Node"="Network"))
  
  return(node_df)
  
}

node_comp1 <- compute_node_contribution("Comp1")
node_comp2 <- compute_node_contribution("Comp2")
node_comp3 <- compute_node_contribution("Comp3")
key_comp1 <- node_comp1 %>% filter(abs(To_z)>2 | abs(From_z)>2)
key_comp2 <- node_comp2 %>% filter(abs(To_z)>2 | abs(From_z)>2)
key_comp3 <- node_comp3 %>% filter(abs(To_z)>2 | abs(From_z)>2)


wb <- createWorkbook()

for (i in c(1,2,3)) {
  
  node_obj <- get(paste0("node_comp", i))
  key_obj  <- get(paste0("key_comp", i))

  addWorksheet(wb, paste0("Node_comp", i))
  writeData(wb, paste0("Node_comp", i), node_obj)
  
  addWorksheet(wb, paste0("Key_comp", i))
  writeData(wb, paste0("Key_comp", i), key_obj)
}

saveWorkbook(
  wb,
  file.path(out_root,
            "PLSDA_KeyNet.xlsx"),
  overwrite=TRUE
)

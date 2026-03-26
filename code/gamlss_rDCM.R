############################################################
# rDCM EC Normative Modeling (GAMLSS)
# TD model → deviation z-score
############################################################

rm(list = ls())

library(tidyverse)
library(readxl)
library(openxlsx)
library(gamlss)
library(gamlss.add)

set.seed(1205)


############################################################
# Developmental trajectory plot
############################################################

plot_trajectory_png <- function(fit, df, out_file){
  
  age_seq <- seq(min(df$Age), max(df$Age), length = 200)
  
  pred <- predictAll(
    fit,
    newdata = data.frame(Age = age_seq),
    type = "response"
  )
  
  png(out_file,
      width = 2400,
      height = 1800,
      res = 300)
  
  plot(df$Age,
       df$EC,
       pch = 16,
       col = rgb(0,0,0,0.3),
       xlab = "年龄",
       ylab = "rDCM")
  
  lines(age_seq, pred$mu, lwd = 3)
  lines(age_seq, pred$mu + pred$sigma, lty = 2)
  lines(age_seq, pred$mu - pred$sigma, lty = 2)
  lines(age_seq, pred$mu + 2*pred$sigma, lty = 3)
  lines(age_seq, pred$mu - 2*pred$sigma, lty = 3)
  
  dev.off()
}

##############
root_dir <- "/Volumes/Zuolab_XRF/output/norm_rDCM/gamlss"

dir.create(root_dir, showWarnings = FALSE)
dir.create(file.path(root_dir,"models"), showWarnings = FALSE)
dir.create(file.path(root_dir,"trajectories"), showWarnings = FALSE)
dir.create(file.path(root_dir,"diagnostics"), showWarnings = FALSE)
dir.create(file.path(root_dir,"deviation"), showWarnings = FALSE)

data <- read_excel("/Volumes/Zuolab_XRF/output/norm_rDCM/rDCM_normative_data_combat.xlsx")

data <- data %>%
  mutate(
    Age = as.numeric(Age),
    Subtype = factor(Subtype,
                     levels=c("TD","ASD-L","ASD-H"))
  ) %>%
  dplyr::select(-Session)

ec_cols <- grep("^EC_", colnames(data), value = TRUE)

td_data <- data %>%
  filter(Subtype == "TD")

deviation_matrix <- data.frame(
  subject = data$ID,
  Age = data$Age,
  Subtype = data$Subtype
)

summary_table <- data.frame()
failure_log <- data.frame()

for (edge in ec_cols) {
  
  cat("Fitting:", edge, "\n")
  
  sub_data <- td_data %>%
    dplyr::select(Age, all_of(edge)) %>%
    rename(EC = all_of(edge)) %>%
    drop_na()
  
  if (nrow(sub_data) < 40) next
  
  fit <- tryCatch({
    
    gamlss(
      EC ~ pb(Age),
      sigma.formula = ~ pb(Age),
      family = NO,
      data = sub_data,
      method = RS(),
      control = gamlss.control(
        n.cyc = 300,
        trace = FALSE
      )
    )
    
  }, error=function(e) NULL)
  
  if (is.null(fit) || !fit$converged) {
    
    failure_log <- rbind(
      failure_log,
      data.frame(Edge=edge)
    )
    
    next
  }
  
  ########################################################
  # 保存模型
  ########################################################
  
  saveRDS(
    fit,
    file.path(root_dir,
              "models",
              paste0(edge,".rds"))
  )
  
  ########################################################
  # Save trajectory
  ########################################################
  plot_trajectory_png(
    fit,
    sub_data,
    file.path(root_dir, "trajectories",
              paste0(edge, "_trajectory.png"))
  )
  
  ########################################################
  # Diagnostic plot
  ########################################################
  
  png(file.path(root_dir, "diagnostics",
                paste0(edge, "_diagnostic.png")),
      width = 2400,
      height = 1800,
      res = 300)
  
  plot(fit)
  
  dev.off()
  
  ########################################################
  # 预测所有人
  ########################################################
  
  pred <- predictAll(
    fit,
    newdata = data.frame(Age=data$Age),
    type="response"
  )
  
  mu <- pred$mu
  sigma <- pred$sigma
  
  observed <- data[[edge]]
  
  zscore <- (observed - mu) / sigma
  
  deviation_matrix[[edge]] <- zscore
  
  ########################################################
  # summary
  ########################################################
  
  summary_table <- rbind(
    summary_table,
    data.frame(
      Edge=edge,
      N=nrow(sub_data),
      DF_total=fit$df.fit,
      AIC=AIC(fit),
      BIC=BIC(fit),
      GlobalDeviance=fit$G.deviance
    )
  )
  
}

write.xlsx(
  deviation_matrix,
  file.path(root_dir,
            "rDCM_deviation_zscores.xlsx"),
  overwrite=TRUE
)

write.xlsx(
  summary_table,
  file.path(root_dir,
            "GAMLSS_summary.xlsx"),
  overwrite=TRUE
)

write.csv(
  failure_log,
  file.path(root_dir,
            "Model_failures.csv"),
  row.names=FALSE
)
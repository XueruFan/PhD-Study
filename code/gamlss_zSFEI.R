############################################################
# GAMLSS Normative Modeling
# Step 1–5 only
# CCNP + ABIDE TD combined
# NO distribution
# Scaled SFEI
# Clean + No environment dependency
############################################################

rm(list = ls())

library(tidyverse)
library(readxl)
library(openxlsx)
library(gamlss)
library(gamlss.add)

set.seed(1234)

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
       df$zSFEI_scaled,
       pch = 16,
       col = rgb(0,0,0,0.3),
       xlab = "年龄",
       ylab = "zSFEI")
  
  lines(age_seq, pred$mu, lwd = 3)
  lines(age_seq, pred$mu + pred$sigma, lty = 2)
  lines(age_seq, pred$mu - pred$sigma, lty = 2)
  lines(age_seq, pred$mu + 2*pred$sigma, lty = 3)
  lines(age_seq, pred$mu - 2*pred$sigma, lty = 3)
  
  dev.off()
}

############################################################
# Directory structure
############################################################

root_dir <- "/Volumes/Zuolab_XRF/output/normative/gamlss"

dir.create(root_dir, showWarnings = FALSE)
dir.create(file.path(root_dir, "models"), showWarnings = FALSE)
dir.create(file.path(root_dir, "trajectories"), showWarnings = FALSE)
dir.create(file.path(root_dir, "diagnostics"), showWarnings = FALSE)

############################################################
# Load data
############################################################

data <- read_excel("/Volumes/Zuolab_XRF/output/normative/zSFEI_normative_data_combat.xlsx")

############################################################
# Preprocess
############################################################

data <- data %>%
  mutate(
    Age = as.numeric(Age),
    Step = as.numeric(as.character(Step)),   # 转 numeric 方便筛选
    Network = as.factor(Network),
    zSFEI_ComBat = as.numeric(zSFEI_ComBat)
  ) %>%
  dplyr::select(-Session) %>%   # 删除 Session
  filter(
    Diagnosis == "TD",
    Step <= 5
  ) %>%
  mutate(
    Step = as.factor(Step),        # 再转回 factor
    zSFEI_scaled = zSFEI_ComBat * 100
  )

############################################################
# Modeling preparation
############################################################

steps <- sort(unique(data$Step))
networks <- sort(unique(data$Network))

summary_table <- data.frame()
failure_log <- data.frame()

############################################################
# Main modeling loop
############################################################

for (s in steps) {
  for (n in networks) {
    
    cat("Fitting:", s, n, "\n")
    
    sub_data <- data %>%
      filter(Step == s,
             Network == n)
    
    if (nrow(sub_data) < 40) next
    
    fit <- tryCatch({
      
      gamlss(
        zSFEI_scaled ~ pb(Age),
        sigma.formula = ~ pb(Age),
        family = NO,             # 固定 NO
        data = sub_data,
        method = RS(),
        control = gamlss.control(
          n.cyc = 300,
          trace = FALSE
        )
      )
      
    }, error = function(e) NULL)
    
    if (is.null(fit) || !fit$converged) {
      
      failure_log <- rbind(
        failure_log,
        data.frame(Step = s, Network = n)
      )
      next
    }
    
    model_name <- paste0("Step", s, "_", n)
    
    ########################################################
    # Save model
    ########################################################
    
    saveRDS(
      fit,
      file.path(root_dir, "models",
                paste0(model_name, ".rds"))
    )
    
    ########################################################
    # Save trajectory
    ########################################################
    
    plot_trajectory_png(
      fit,
      sub_data,
      file.path(root_dir, "trajectories",
                paste0(model_name, "_trajectory.png"))
    )
    
    ########################################################
    # Diagnostic plot
    ########################################################
    
    png(file.path(root_dir, "diagnostics",
                  paste0(model_name, "_diagnostic.png")),
        width = 2400,
        height = 1800,
        res = 300)
    
    plot(fit)
    
    dev.off()
    
    ########################################################
    # Save summary info
    ########################################################
    
    summary_table <- rbind(
      summary_table,
      data.frame(
        Step = s,
        Network = n,
        N = nrow(sub_data),
        DF_total = fit$df.fit,
        AIC = AIC(fit),
        BIC = BIC(fit),
        GlobalDeviance = fit$G.deviance
      )
    )
    
  }
}

############################################################
# Save results
############################################################

write.xlsx(
  summary_table,
  file.path(root_dir, "GAMLSS_summary.xlsx"),
  overwrite = TRUE
)

write.csv(
  failure_log,
  file.path(root_dir, "Model_failures.csv"),
  row.names = FALSE
)

cat("Normative modeling completed successfully.\n")

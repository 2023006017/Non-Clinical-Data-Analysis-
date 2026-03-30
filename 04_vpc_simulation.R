# =============================================================================
# CNS Population PK Project
# Script 04: Visual Predictive Check (VPC) & Simulations
# =============================================================================
# VPC: Gold standard for model evaluation in PopPK
# Method: nlmixr2 vpcPlot or manual via vpc package
# =============================================================================

library(tidyverse)
library(nlmixr2)
library(vpc)
library(ggplot2)
library(patchwork)

theme_pk <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      strip.background = element_rect(fill = "#2C3E50", color = NA),
      strip.text       = element_text(color = "white", face = "bold"),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "#2C3E50"),
      axis.title       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", size = 14)
    )
}

# =============================================================================
# LOAD MODEL AND DATA
# =============================================================================

fit_final <- readRDS("output/fit_cov_model.rds")
pk_obs    <- read_csv("data/pk_obs.csv")

# =============================================================================
# VPC — METHOD 1: nlmixr2 built-in vpcPlot
# =============================================================================

# Run VPC simulation (500 replicates recommended; use 200 for speed)
vpc_result <- vpcPlot(
  fit_final,
  n           = 500,
  bins        = "jenks",   # Jenks natural breaks
  pi          = c(0.05, 0.95),  # 90% PI
  ci          = c(0.025, 0.975),# 95% CI on quantiles
  show        = list(obs_dv = TRUE, obs_ci = TRUE, pi = TRUE, pi_ci = TRUE),
  log_y       = TRUE,
  title       = "Visual Predictive Check — Final Covariate Model",
  subtitle    = "N=500 simulations | 90% Prediction Interval",
  xlab        = "Time (hours post-dose)",
  ylab        = "Plasma Concentration (ng/mL)"
)

ggsave("plots/09_vpc_overall.png", vpc_result, width = 10, height = 7, dpi = 300)

# =============================================================================
# VPC STRATIFIED BY DOSE GROUP
# =============================================================================

vpc_dose <- vpcPlot(
  fit_final,
  n          = 500,
  stratify   = "DOSE",
  bins       = "jenks",
  pi         = c(0.05, 0.95),
  log_y      = TRUE,
  title      = "VPC Stratified by Dose Group",
  xlab       = "Time (hours)",
  ylab       = "Concentration (ng/mL)"
)

ggsave("plots/10_vpc_by_dose.png", vpc_dose, width = 14, height = 5, dpi = 300)

# =============================================================================
# VPC — METHOD 2: Manual (using vpc package — more control)
# For users who want full customization
# =============================================================================

# Simulate from final model
set.seed(123)
n_sim <- 500

sim_data <- nlmixr2(
  fit_final,
  est     = "sim",
  nSim    = n_sim,
  data    = pk_obs |>
              select(ID, TIME, DV, WT, AGE, SEX, eGFR, ALB, SMOK)
)

# Manual VPC with vpc package
manual_vpc <- vpc(
  sim             = sim_data,
  obs             = pk_obs |> rename(time = TIME, dv = DV, id = ID),
  obs_cols        = list(idv = "time", dv = "dv", id = "id"),
  sim_cols        = list(idv = "TIME", dv = "DV", id = "ID"),
  bins            = c(0, 1, 2, 4, 6, 8, 12, 24),
  pi              = c(0.05, 0.95),
  ci              = c(0.025, 0.975),
  lloq            = 1.0,
  log_y           = TRUE,
  plot            = TRUE,
  title           = "VPC with LLOQ Handling — nlmixr2 Final Model",
  xlabel          = "Time (hours post-dose)",
  ylabel          = "Concentration (ng/mL)"
)

ggsave("plots/11_vpc_manual_lloq.png", manual_vpc, width = 10, height = 7, dpi = 300)

# =============================================================================
# COVARIATE EFFECT VISUALIZATION (Forest Plot)
# =============================================================================

# Extract parameter estimates and CIs
param_est <- as_tibble(fit_final$parFixed) |>
  rownames_to_column("Parameter") |>
  filter(str_detect(Parameter, "dCL|dV1")) |>
  rename(Estimate = Estimate, SE = "SE", CI_low = `2.5%`, CI_high = `97.5%`)

# Create forest plot
p_forest <- ggplot(param_est, aes(x = Estimate, y = reorder(Parameter, Estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.3, linewidth = 1.2) +
  geom_point(size = 4, color = "#E74C3C") +
  labs(
    title    = "Covariate Effect Estimates — Forest Plot",
    subtitle = "95% Confidence Intervals | Covariates on PK Parameters",
    x        = "Estimated Coefficient",
    y        = "Covariate Term"
  ) +
  theme_pk()

ggsave("plots/12_forest_plot_covariates.png", p_forest, width = 9, height = 5, dpi = 300)

# =============================================================================
# SIMULATION: DOSE vs EXPOSURE RELATIONSHIP
# =============================================================================

# Simulate typical subject (median demographics)
sim_doses    <- seq(10, 150, by = 5)
time_grid    <- seq(0, 24, by = 0.5)

typical_subs <- expand_grid(
  DOSE = sim_doses,
  TIME = time_grid
) |>
  mutate(
    ID   = as.integer(factor(DOSE)),
    EVID = if_else(TIME == 0, 1L, 0L),
    AMT  = if_else(EVID == 1, DOSE, 0),
    CMT  = 1L,
    MDV  = if_else(EVID == 1, 1L, 0L),
    DV   = NA_real_,
    # Median covariates
    WT   = 70, AGE = 40, SEX = 0, eGFR = 85, ALB = 4.1, SMOK = 0
  )

# Simulate expected exposure
sim_exposure <- nlmixr2(
  fit_final,
  est  = "sim",
  data = typical_subs,
  nSim = 100
)

# Cmax and AUC vs Dose
exposure_summary <- sim_exposure |>
  filter(EVID == 0) |>
  group_by(sim, DOSE) |>
  summarise(
    Cmax = max(DV, na.rm = TRUE),
    AUC  = sum(diff(TIME) * (head(DV, -1) + tail(DV, -1)) / 2),
    .groups = "drop"
  ) |>
  group_by(DOSE) |>
  summarise(
    Cmax_med  = median(Cmax),
    Cmax_lo   = quantile(Cmax, 0.05),
    Cmax_hi   = quantile(Cmax, 0.95),
    AUC_med   = median(AUC),
    AUC_lo    = quantile(AUC, 0.05),
    AUC_hi    = quantile(AUC, 0.95),
    .groups   = "drop"
  )

p_exp1 <- ggplot(exposure_summary, aes(x = DOSE)) +
  geom_ribbon(aes(ymin = Cmax_lo, ymax = Cmax_hi), fill = "#3498DB", alpha = 0.3) +
  geom_line(aes(y = Cmax_med), color = "#3498DB", linewidth = 1.5) +
  geom_vline(xintercept = c(25, 50, 100), linetype = "dashed", color = "#E74C3C") +
  labs(title = "Cmax vs Dose", x = "Dose (mg)", y = "Cmax (ng/mL)",
       subtitle = "Median ± 90% PI | Vertical lines = study doses") +
  theme_pk()

p_exp2 <- ggplot(exposure_summary, aes(x = DOSE)) +
  geom_ribbon(aes(ymin = AUC_lo, ymax = AUC_hi), fill = "#2ECC71", alpha = 0.3) +
  geom_line(aes(y = AUC_med), color = "#2ECC71", linewidth = 1.5) +
  geom_vline(xintercept = c(25, 50, 100), linetype = "dashed", color = "#E74C3C") +
  labs(title = "AUC0-24 vs Dose", x = "Dose (mg)", y = "AUC (ng·hr/mL)",
       subtitle = "Median ± 90% PI | Vertical lines = study doses") +
  theme_pk()

p_exposure <- p_exp1 | p_exp2 +
  plot_annotation(title = "Simulated Exposure-Response Relationship")

ggsave("plots/13_exposure_vs_dose.png", p_exposure, width = 12, height = 5, dpi = 300)

cat("\n=== VPC & Simulation Complete ===\n")
cat("Plots saved:\n")
cat("  09_vpc_overall.png\n  10_vpc_by_dose.png\n")
cat("  11_vpc_manual_lloq.png\n  12_forest_plot_covariates.png\n")
cat("  13_exposure_vs_dose.png\n")

# =============================================================================
# CNS Population PK Project
# Script 03: Population PK Model — Two-Compartment Oral (nlmixr2)
# =============================================================================
# Drug:   Hypothetical CNS Agent
# Model:  2-CMT oral, first-order absorption
# Engine: nlmixr2 (FOCEI estimation)
# Ref:    Fidler et al. (2019) CPT:PSP | NONMEM-compatible parameterization
# =============================================================================

library(tidyverse)
library(nlmixr2)
library(ggplot2)
library(patchwork)
library(xpose)       # Diagnostic plots (xpose4-compatible)
library(xpose.nlmixr2)  # nlmixr2 bridge to xpose

# =============================================================================
# THEME
# =============================================================================

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
# LOAD & FORMAT DATA
# =============================================================================

pk_data <- read_csv("data/pk_data_full.csv") |>
  # nlmixr2 requires: ID, TIME, DV, AMT, EVID, CMT
  mutate(
    DV  = if_else(EVID == 1 | MDV == 1 | BLQ == 1, NA_real_, DV),
    AMT = if_else(EVID == 1, as.numeric(DOSE), 0)
  ) |>
  select(ID, TIME, EVID, AMT, DV, CMT, MDV, WT, AGE, SEX, eGFR, ALB, SMOK) |>
  arrange(ID, TIME, desc(EVID))

cat("Dataset dimensions:", nrow(pk_data), "rows x", ncol(pk_data), "cols\n")

# =============================================================================
# MODEL 1: BASE MODEL (No Covariates)
# Two-compartment, first-order absorption, proportional + additive error
# =============================================================================

base_model <- function() {
  ini({
    # Fixed effects (log-scale for positivity)
    lCL  <- log(8.5)   ; label("Clearance (L/hr)")
    lV1  <- log(45)    ; label("Central volume (L)")
    lQ   <- log(4.2)   ; label("Intercompartmental CL (L/hr)")
    lV2  <- log(80)    ; label("Peripheral volume (L)")
    lKa  <- log(1.2)   ; label("Absorption rate constant (1/hr)")

    # IIV (variance, exponential model)
    eta_CL ~ 0.08   # ~ 28% CV
    eta_V1 ~ 0.10   # ~ 32% CV
    eta_Ka ~ 0.20   # ~ 45% CV
    eta_Q  ~ 0.04
    eta_V2 ~ 0.06

    # Residual error
    prop_err <- 0.15   # Proportional
    add_err  <- 2.0    # Additive (ng/mL)
  })

  model({
    # Individual parameters (exponential IIV)
    CL <- exp(lCL + eta_CL)
    V1 <- exp(lV1 + eta_V1)
    Q  <- exp(lQ  + eta_Q)
    V2 <- exp(lV2 + eta_V2)
    Ka <- exp(lKa + eta_Ka)

    # Two-compartment ODE (linCmt() uses analytical 2CMT solution)
    linCmt()

    # Combined error model
    dv ~ prop(prop_err) + add(add_err)
  })
}

# Fit base model
fit_base <- nlmixr2(
  object  = base_model,
  data    = pk_data,
  est     = "focei",
  control = foceiControl(print = 5, maxOuterIterations = 100)
)

# Model summary
print(summary(fit_base))

# Save base model output
saveRDS(fit_base, "output/fit_base_model.rds")

# =============================================================================
# MODEL 2: COVARIATE MODEL
# Covariates: WT on CL and V1 (allometric), SEX on CL, eGFR on CL
# Selection based on EDA findings and clinical rationale
# =============================================================================

cov_model <- function() {
  ini({
    lCL       <- log(8.5)
    lV1       <- log(45)
    lQ        <- log(4.2)
    lV2       <- log(80)
    lKa       <- log(1.2)
    # Covariate fixed effects
    dCL_WT    <- 0.75    ; label("WT exponent on CL (allometric)")
    dV1_WT    <- 1.00    ; label("WT exponent on V1 (allometric)")
    dCL_SEX   <- -0.20   ; label("Sex effect on CL (Female vs Male)")
    dCL_eGFR  <- 0.30    ; label("eGFR exponent on CL")

    eta_CL ~ 0.06
    eta_V1 ~ 0.09
    eta_Ka ~ 0.18
    eta_Q  ~ 0.03
    eta_V2 ~ 0.05

    prop_err <- 0.12
    add_err  <- 1.5
  })

  model({
    # Standardized covariates (centered at median)
    WT_std   <- WT   / 70
    eGFR_std <- eGFR / 85

    # Covariate-adjusted parameters
    CL <- exp(lCL + dCL_WT * log(WT_std) + dCL_SEX * SEX + dCL_eGFR * log(eGFR_std) + eta_CL)
    V1 <- exp(lV1 + dV1_WT * log(WT_std) + eta_V1)
    Q  <- exp(lQ  + eta_Q)
    V2 <- exp(lV2 + eta_V2)
    Ka <- exp(lKa + eta_Ka)

    linCmt()

    dv ~ prop(prop_err) + add(add_err)
  })
}

fit_cov <- nlmixr2(
  object  = cov_model,
  data    = pk_data,
  est     = "focei",
  control = foceiControl(print = 5, maxOuterIterations = 150)
)

print(summary(fit_cov))
saveRDS(fit_cov, "output/fit_cov_model.rds")

# =============================================================================
# MODEL COMPARISON
# =============================================================================

cat("\n=== MODEL COMPARISON ===\n")
comparison <- tibble(
  Model      = c("Base", "Covariate"),
  OFV        = c(fit_base$objective, fit_cov$objective),
  AIC        = c(AIC(fit_base), AIC(fit_cov)),
  BIC        = c(BIC(fit_base), BIC(fit_cov)),
  N_params   = c(length(coef(fit_base)), length(coef(fit_cov)))
) |>
  mutate(
    ΔOFV = OFV - min(OFV),
    ΔAIC = AIC - min(AIC)
  )

print(comparison)
write_csv(comparison, "output/model_comparison.csv")

# LRT: Covariate model vs Base model
delta_ofv   <- fit_base$objective - fit_cov$objective
delta_df    <- length(coef(fit_cov)) - length(coef(fit_base))
p_val       <- pchisq(delta_ofv, df = delta_df, lower.tail = FALSE)
cat(sprintf("\nLRT: ΔOFV = %.2f, Δdf = %d, p = %.4f\n", delta_ofv, delta_df, p_val))
cat("Conclusion:", if (p_val < 0.05) "Covariate model SIGNIFICANTLY better\n" else "No significant improvement\n")

# =============================================================================
# SELECT FINAL MODEL
# =============================================================================

fit_final <- fit_cov   # Covariate model selected

# =============================================================================
# GOODNESS-OF-FIT PLOTS (Final Model)
# =============================================================================

gof_data <- fit_final |>
  as_tibble() |>
  mutate(CWRES = CWRES, IPRED = IPRED, PRED = PRED)

# DV vs IPRED
p_gof1 <- ggplot(gof_data, aes(x = IPRED, y = DV)) +
  geom_point(alpha = 0.5, color = "#3498DB", size = 1.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 1) +
  geom_smooth(method = "loess", color = "#E74C3C", se = TRUE, linewidth = 1.2) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "DV vs IPRED", x = "Individual Predicted (ng/mL)", y = "Observed (ng/mL)") +
  theme_pk()

# DV vs PRED
p_gof2 <- ggplot(gof_data, aes(x = PRED, y = DV)) +
  geom_point(alpha = 0.5, color = "#2ECC71", size = 1.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 1) +
  geom_smooth(method = "loess", color = "#E74C3C", se = TRUE, linewidth = 1.2) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "DV vs PRED", x = "Population Predicted (ng/mL)", y = "Observed (ng/mL)") +
  theme_pk()

# CWRES vs TIME
p_gof3 <- ggplot(gof_data, aes(x = TIME, y = CWRES)) +
  geom_point(alpha = 0.5, color = "#9B59B6", size = 1.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "red") +
  geom_smooth(method = "loess", color = "#E74C3C", se = TRUE) +
  labs(title = "CWRES vs TIME", x = "Time (hr)", y = "CWRES") +
  theme_pk()

# CWRES vs PRED
p_gof4 <- ggplot(gof_data, aes(x = PRED, y = CWRES)) +
  geom_point(alpha = 0.5, color = "#E67E22", size = 1.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "red") +
  geom_smooth(method = "loess", color = "#E74C3C", se = TRUE) +
  scale_x_log10() +
  labs(title = "CWRES vs PRED", x = "Population Predicted (ng/mL)", y = "CWRES") +
  theme_pk()

gof_panel <- (p_gof1 | p_gof2) / (p_gof3 | p_gof4) +
  plot_annotation(
    title    = "Goodness-of-Fit Diagnostics — Final Covariate Model",
    subtitle = "Two-Compartment Oral PopPK | CNS Agent | FOCEI Estimation"
  )

ggsave("plots/07_goodness_of_fit.png", gof_panel, width = 13, height = 10, dpi = 300)

# =============================================================================
# CWRES DISTRIBUTION (Normality Check)
# =============================================================================

p_cwres <- ggplot(gof_data, aes(x = CWRES)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "#3498DB", alpha = 0.7, color = "white") +
  stat_function(fun = dnorm, color = "#E74C3C", linewidth = 1.5) +
  labs(title = "CWRES Distribution",
       subtitle = "Should approximate N(0,1) — Red line = standard normal",
       x = "Conditional Weighted Residuals", y = "Density") +
  theme_pk()

ggsave("plots/08_cwres_distribution.png", p_cwres, width = 7, height = 5, dpi = 300)

cat("\n=== PopPK Modeling Complete ===\n")
cat("Final model: Two-Compartment Oral with Weight, Sex, eGFR covariates\n")
cat("Results saved to output/\n")

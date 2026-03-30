# =============================================================================
# CNS Population PK Project
# Script 01: Simulate Realistic PK Dataset
# Drug: Hypothetical CNS Agent (e.g., novel anxiolytic / antidepressant)
# Author: [Your Name]
# Date: 2024
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(tidyverse)
library(MASS)   # mvrnorm for correlated random effects

set.seed(42)

# =============================================================================
# STUDY DESIGN PARAMETERS
# =============================================================================

n_subjects  <- 120       # Total subjects
n_dose_grp  <- 3         # Dose groups: 25 mg, 50 mg, 100 mg
doses       <- c(25, 50, 100)  # mg
dose_alloc  <- rep(doses, each = n_subjects / n_dose_grp)

# Sampling times (hours) — sparse PK design typical in PopPK
# Pre-dose + multiple post-dose
sparse_times <- c(0, 0.5, 1, 2, 4, 6, 8, 12, 24)

# =============================================================================
# TRUE POPULATION PK PARAMETERS (Two-Compartment Oral Model)
# Based on realistic CNS small-molecule (MW ~300-400 Da, moderate BBB penetration)
# =============================================================================

# Fixed Effects (theta)
theta <- list(
  CL    = 8.5,    # Clearance (L/hr)      — typical oral CNS drug
  V1    = 45,     # Central volume (L)
  Q     = 4.2,    # Intercompartmental CL (L/hr)
  V2    = 80,     # Peripheral volume (L)
  Ka    = 1.2,    # Absorption rate (1/hr)
  F     = 0.75    # Bioavailability (fraction)
)

# Covariate effects (realistic)
beta_WT_CL  <- 0.75   # Weight effect on CL (allometric)
beta_WT_V1  <- 1.00   # Weight effect on V1 (allometric)
beta_SEX_CL <- -0.20  # Female: 20% lower CL (sex effect)
beta_AGE_CL <- -0.008 # Age effect on CL (per year)
beta_eGFR   <- 0.30   # Renal function effect on CL

# Interindividual variability (IIV) — omega (CV%)
omega_CL  <- 0.28   # ~28% CV
omega_V1  <- 0.32   # ~32% CV
omega_Ka  <- 0.45   # ~45% CV (absorption most variable)
omega_Q   <- 0.20
omega_V2  <- 0.25

# Residual error
sigma_prop <- 0.15  # Proportional error (15%)
sigma_add  <- 0.5   # Additive error (ng/mL)

# =============================================================================
# SIMULATE SUBJECT-LEVEL COVARIATES
# =============================================================================

subjects <- tibble(
  ID     = 1:n_subjects,
  DOSE   = dose_alloc,
  WT     = rnorm(n_subjects, mean = 72, sd = 12) |> pmax(45) |> pmin(130),  # kg
  AGE    = round(runif(n_subjects, 22, 68)),                                  # years
  SEX    = sample(c(0, 1), n_subjects, replace = TRUE, prob = c(0.45, 0.55)), # 0=Male, 1=Female
  eGFR   = rnorm(n_subjects, mean = 85, sd = 18) |> pmax(30) |> pmin(120),   # mL/min/1.73m2
  ALB    = rnorm(n_subjects, mean = 4.1, sd = 0.4) |> pmax(2.5) |> pmin(5.5),# g/dL
  RACE   = sample(c("Asian", "White", "Black", "Other"),
                  n_subjects, replace = TRUE, prob = c(0.20, 0.55, 0.18, 0.07)),
  SMOK   = sample(c(0, 1), n_subjects, replace = TRUE, prob = c(0.75, 0.25))  # 0=No, 1=Yes
)

# =============================================================================
# SIMULATE RANDOM EFFECTS (correlated IIV)
# =============================================================================

Omega <- diag(c(omega_CL^2, omega_V1^2, omega_Ka^2, omega_Q^2, omega_V2^2))
# Small off-diagonal correlation between CL and V1
Omega[1, 2] <- Omega[2, 1] <- 0.15 * omega_CL * omega_V1

eta <- mvrnorm(n_subjects, mu = rep(0, 5), Sigma = Omega)
colnames(eta) <- c("eta_CL", "eta_V1", "eta_Ka", "eta_Q", "eta_V2")

subjects <- bind_cols(subjects, as_tibble(eta))

# =============================================================================
# TWO-COMPARTMENT ORAL PK SOLUTION (Analytical)
# =============================================================================

# Closed-form solution for 2-compartment model with first-order absorption
pk_2cmt_oral <- function(t, DOSE, CL, V1, Q, V2, Ka, F_oral) {
  
  # Micro-rate constants
  k10 <- CL / V1
  k12 <- Q  / V1
  k21 <- Q  / V2
  
  alpha_sum  <- k10 + k12 + k21
  alpha_prod <- k10 * k21
  
  alpha <- (alpha_sum + sqrt(alpha_sum^2 - 4 * alpha_prod)) / 2
  beta  <- (alpha_sum - sqrt(alpha_sum^2 - 4 * alpha_prod)) / 2
  
  # Coefficients
  A <- (Ka * (k21 - alpha)) / (V1 * (alpha - beta) * (Ka - alpha))
  B <- (Ka * (k21 - beta))  / (V1 * (alpha - beta) * (Ka - beta))
  C <- (Ka * (k21 - Ka))    / (V1 * (alpha - Ka)   * (beta - Ka))
  
  amt <- DOSE * F_oral  # Absorbed dose (mg) -> scale to ng/mL via units
  
  # Concentration (ng/mL) — assumes dose in mg, volumes in L, CL in L/hr
  # Multiply by 1000 for mg -> ug/L conversion factor
  Cp <- ifelse(t > 0,
               1000 * amt * (A * exp(-alpha * t) + B * exp(-beta * t) + C * exp(-Ka * t)),
               0)
  
  pmax(Cp, 0)  # floor at 0
}

# =============================================================================
# GENERATE FULL PK DATASET
# =============================================================================

pk_data <- subjects |>
  mutate(
    # Individual PK parameters with covariate effects
    CL_ind = theta$CL *
              (WT / 70)^beta_WT_CL *
              exp(beta_SEX_CL * SEX) *
              exp(beta_AGE_CL * (AGE - 40)) *
              (eGFR / 85)^beta_eGFR *
              exp(eta_CL),
    
    V1_ind = theta$V1 * (WT / 70)^beta_WT_V1 * exp(eta_V1),
    Ka_ind = theta$Ka * exp(eta_Ka),
    Q_ind  = theta$Q  * exp(eta_Q),
    V2_ind = theta$V2 * exp(eta_V2)
  ) |>
  # Expand to all sampling time points
  crossing(TIME = sparse_times) |>
  arrange(ID, TIME) |>
  mutate(
    # Predicted concentration (IPRED)
    IPRED = pmap_dbl(
      list(TIME, DOSE, CL_ind, V1_ind, Q_ind, V2_ind, Ka_ind),
      function(t, d, cl, v1, q, v2, ka) {
        pk_2cmt_oral(t, d, cl, v1, q, v2, ka, F_oral = theta$F)
      }
    ),
    # Add residual error
    EPS1  = rnorm(n(), 0, sigma_prop),
    EPS2  = rnorm(n(), 0, sigma_add),
    DV    = pmax(IPRED * (1 + EPS1) + EPS2, 0.01),  # Observed DV (ng/mL)
    # Log-transform for modeling
    logDV = log(DV),
    # BLQ flag (LLOQ = 1.0 ng/mL)
    LLOQ  = 1.0,
    BLQ   = as.integer(DV < LLOQ),
    MDV   = as.integer(TIME == 0),  # Pre-dose MDV=1
    # EVID: 0=observation, 1=dose
    EVID  = 0L,
    CMT   = 1L,
    AMT   = 0,
    # Standardize covariates
    SEX_label  = ifelse(SEX == 0, "Male", "Female"),
    DOSE_label = paste0(DOSE, " mg")
  ) |>
  # Add dose records (for NONMEM-style format)
  bind_rows(
    subjects |>
      mutate(
        TIME  = 0, EVID = 1L, CMT = 1L,
        AMT   = DOSE, DV = NA, IPRED = NA, MDV = 1L, BLQ = 0L,
        SEX_label  = ifelse(SEX == 0, "Male", "Female"),
        DOSE_label = paste0(DOSE, " mg")
      )
  ) |>
  arrange(ID, TIME, EVID) |>
  select(ID, TIME, EVID, AMT, DV, IPRED, MDV, BLQ, CMT,
         DOSE, WT, AGE, SEX, SEX_label, eGFR, ALB, RACE, SMOK, DOSE_label,
         CL_ind, V1_ind, Ka_ind, Q_ind, V2_ind,
         eta_CL, eta_V1, eta_Ka, LLOQ)

# =============================================================================
# SAVE DATASETS
# =============================================================================

# Full dataset
write_csv(pk_data, "data/pk_data_full.csv")

# Observations only (for analysis)
pk_obs <- pk_data |>
  filter(EVID == 0, BLQ == 0, MDV == 0) |>
  select(-EVID, -AMT)

write_csv(pk_obs, "data/pk_obs.csv")

# NONMEM-style dataset
nm_data <- pk_data |>
  select(ID, TIME, EVID, AMT, DV, MDV, BLQ, CMT, DOSE, WT, AGE, SEX, eGFR, ALB, SMOK) |>
  mutate(across(where(is.numeric), ~ round(.x, 4))) |>
  mutate(DV = ifelse(BLQ == 1 | is.na(DV), ".", as.character(round(DV, 4))))

write_csv(nm_data, "data/nm_dataset.csv")

cat("\n=== Dataset Summary ===\n")
cat("Total subjects:", n_subjects, "\n")
cat("Dose groups:", paste(doses, "mg", collapse = ", "), "\n")
cat("Total observations (non-BLQ):", nrow(pk_obs), "\n")
cat("BLQ rate:", round(mean(pk_data$BLQ[pk_data$EVID == 0], na.rm = TRUE) * 100, 1), "%\n")
cat("Median CL:", round(median(pk_obs$CL_ind, na.rm = TRUE), 2), "L/hr\n")
cat("Median V1:", round(median(pk_obs$V1_ind, na.rm = TRUE), 2), "L\n")
cat("\nFiles saved to data/\n")

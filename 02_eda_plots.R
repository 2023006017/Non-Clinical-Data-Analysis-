# =============================================================================
# CNS Population PK Project
# Script 02: Exploratory Data Analysis (EDA)
# =============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(GGally)
library(ggridges)

# Custom theme for publication-quality plots
theme_pk <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      strip.background = element_rect(fill = "#2C3E50", color = NA),
      strip.text       = element_text(color = "white", face = "bold"),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "#2C3E50"),
      axis.title       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(color = "grey40", size = 10),
      legend.position  = "bottom",
      legend.frame     = element_rect(color = "grey80")
    )
}

# Color palette (colorblind-friendly)
dose_colors <- c("25 mg" = "#3498DB", "50 mg" = "#2ECC71", "100 mg" = "#E74C3C")

# =============================================================================
# LOAD DATA
# =============================================================================

pk_obs <- read_csv("data/pk_obs.csv") |>
  mutate(
    DOSE_label = factor(paste0(DOSE, " mg"), levels = c("25 mg", "50 mg", "100 mg")),
    SEX_label  = factor(SEX_label, levels = c("Male", "Female")),
    logDV      = log(DV)
  )

cat("=== EDA Summary ===\n")
cat("Subjects:", n_distinct(pk_obs$ID), "\n")
cat("Observations:", nrow(pk_obs), "\n")
cat("Time range:", range(pk_obs$TIME), "\n")
cat("Conc range:", round(range(pk_obs$DV), 2), "ng/mL\n")

# =============================================================================
# PLOT 1: Spaghetti Plot — Individual PK Profiles by Dose
# =============================================================================

p1 <- ggplot(pk_obs, aes(x = TIME, y = DV, group = ID, color = DOSE_label)) +
  geom_line(alpha = 0.35, linewidth = 0.5) +
  geom_point(alpha = 0.5, size = 1.2) +
  stat_summary(aes(group = DOSE_label), fun = median,
               geom = "line", linewidth = 1.5, linetype = "solid") +
  scale_color_manual(values = dose_colors, name = "Dose Group") +
  scale_y_log10(labels = scales::label_comma()) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 12, 24)) +
  facet_wrap(~DOSE_label, ncol = 3) +
  labs(
    title    = "Individual Concentration-Time Profiles",
    subtitle = "Each line = 1 subject | Bold line = median | Log10 Y-axis",
    x        = "Time (hours post-dose)",
    y        = "Plasma Concentration (ng/mL)",
    caption  = "Drug: Hypothetical CNS Agent | Two-Compartment Oral Model"
  ) +
  theme_pk() +
  theme(legend.position = "none")

ggsave("plots/01_spaghetti_plot.png", p1, width = 12, height = 5, dpi = 300)

# =============================================================================
# PLOT 2: Observed vs IPRED (Basic Goodness of Fit)
# =============================================================================

p2 <- ggplot(pk_obs, aes(x = IPRED, y = DV, color = DOSE_label)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
  geom_smooth(method = "loess", se = TRUE, color = "#E67E22", fill = "#F39C12", alpha = 0.2) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = dose_colors, name = "Dose Group") +
  labs(
    title    = "Observed vs. Individual Predicted Concentrations",
    subtitle = "Dashed line = line of identity | Log10 scale",
    x        = "Individual Predicted (ng/mL)",
    y        = "Observed DV (ng/mL)"
  ) +
  theme_pk()

ggsave("plots/02_obs_vs_pred.png", p2, width = 7, height = 6, dpi = 300)

# =============================================================================
# PLOT 3: Covariate Distributions
# =============================================================================

cov_data <- pk_obs |>
  distinct(ID, .keep_all = TRUE)

p3a <- ggplot(cov_data, aes(x = WT, fill = SEX_label)) +
  geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("Male" = "#3498DB", "Female" = "#E91E8C"), name = "Sex") +
  labs(title = "Body Weight", x = "Weight (kg)", y = "Count") +
  theme_pk()

p3b <- ggplot(cov_data, aes(x = AGE, fill = DOSE_label)) +
  geom_histogram(bins = 20, alpha = 0.7) +
  scale_fill_manual(values = dose_colors, name = "Dose") +
  labs(title = "Age Distribution", x = "Age (years)", y = "Count") +
  theme_pk()

p3c <- ggplot(cov_data, aes(x = eGFR, fill = SEX_label)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Male" = "#3498DB", "Female" = "#E91E8C"), name = "Sex") +
  labs(title = "eGFR Distribution", x = "eGFR (mL/min/1.73m²)", y = "Density") +
  theme_pk()

p3d <- ggplot(cov_data, aes(x = RACE, fill = RACE)) +
  geom_bar(alpha = 0.85) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Race/Ethnicity", x = "", y = "N Subjects") +
  theme_pk() + theme(legend.position = "none")

p3 <- (p3a + p3b) / (p3c + p3d) +
  plot_annotation(
    title    = "Covariate Distribution Summary",
    subtitle = "N = 120 subjects across 3 dose groups"
  )

ggsave("plots/03_covariate_distributions.png", p3, width = 12, height = 8, dpi = 300)

# =============================================================================
# PLOT 4: Covariate vs Individual PK Parameters
# =============================================================================

param_data <- pk_obs |> distinct(ID, .keep_all = TRUE)

p4a <- ggplot(param_data, aes(x = WT, y = CL_ind, color = SEX_label)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
  scale_color_manual(values = c("Male" = "#3498DB", "Female" = "#E91E8C"), name = "Sex") +
  labs(x = "Weight (kg)", y = "CL (L/hr)", title = "Weight vs. CL") +
  theme_pk()

p4b <- ggplot(param_data, aes(x = eGFR, y = CL_ind, color = DOSE_label)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
  scale_color_manual(values = dose_colors, name = "Dose") +
  labs(x = "eGFR (mL/min/1.73m²)", y = "CL (L/hr)", title = "eGFR vs. CL") +
  theme_pk()

p4c <- ggplot(param_data, aes(x = AGE, y = CL_ind, color = SEX_label)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
  scale_color_manual(values = c("Male" = "#3498DB", "Female" = "#E91E8C"), name = "Sex") +
  labs(x = "Age (years)", y = "CL (L/hr)", title = "Age vs. CL") +
  theme_pk()

p4d <- ggplot(param_data, aes(x = SEX_label, y = CL_ind, fill = SEX_label)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = 21) +
  scale_fill_manual(values = c("Male" = "#3498DB", "Female" = "#E91E8C")) +
  labs(x = "Sex", y = "CL (L/hr)", title = "Sex Effect on CL") +
  theme_pk() + theme(legend.position = "none")

p4 <- (p4a + p4b) / (p4c + p4d) +
  plot_annotation(
    title    = "Covariate Relationships with Individual Clearance",
    subtitle = "Individual PK parameters from true simulation"
  )

ggsave("plots/04_covariate_pk_relationships.png", p4, width = 12, height = 8, dpi = 300)

# =============================================================================
# PLOT 5: Boxplots — Cmax & AUC by Dose Group
# =============================================================================

pk_summary <- pk_obs |>
  group_by(ID, DOSE_label, SEX_label) |>
  summarise(
    Cmax  = max(DV, na.rm = TRUE),
    Ctrough = DV[TIME == 24][1],
    CL_ind  = first(CL_ind),
    .groups = "drop"
  ) |>
  # Approximate AUC by trapezoidal rule per subject
  left_join(
    pk_obs |>
      arrange(ID, TIME) |>
      group_by(ID) |>
      summarise(AUC = sum(diff(TIME) * (head(DV, -1) + tail(DV, -1)) / 2), .groups = "drop"),
    by = "ID"
  )

p5a <- ggplot(pk_summary, aes(x = DOSE_label, y = Cmax, fill = DOSE_label)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = 21, outlier.alpha = 0.5) +
  scale_fill_manual(values = dose_colors) +
  scale_y_log10() +
  labs(x = "Dose Group", y = "Cmax (ng/mL)", title = "Cmax by Dose Group") +
  theme_pk() + theme(legend.position = "none")

p5b <- ggplot(pk_summary, aes(x = DOSE_label, y = AUC, fill = DOSE_label)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = 21, outlier.alpha = 0.5) +
  scale_fill_manual(values = dose_colors) +
  scale_y_log10() +
  labs(x = "Dose Group", y = "AUC0-24 (ng·hr/mL)", title = "AUC by Dose Group") +
  theme_pk() + theme(legend.position = "none")

p5c <- ggplot(pk_summary, aes(x = SEX_label, y = CL_ind, fill = SEX_label)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.2, fill = "white") +
  scale_fill_manual(values = c("Male" = "#3498DB", "Female" = "#E91E8C")) +
  facet_wrap(~DOSE_label) +
  labs(x = "Sex", y = "CL (L/hr)", title = "Clearance by Sex and Dose") +
  theme_pk() + theme(legend.position = "none")

p5 <- (p5a | p5b) / p5c +
  plot_annotation(title = "PK Exposure Metrics by Dose and Demographics")

ggsave("plots/05_pk_exposure_metrics.png", p5, width = 12, height = 9, dpi = 300)

# =============================================================================
# PLOT 6: IIV Eta Distributions
# =============================================================================

eta_data <- pk_obs |>
  distinct(ID, .keep_all = TRUE) |>
  pivot_longer(cols = c(eta_CL, eta_V1, eta_Ka),
               names_to = "Parameter", values_to = "ETA") |>
  mutate(Parameter = recode(Parameter,
    eta_CL = "η_CL (Clearance)",
    eta_V1 = "η_V1 (Central Volume)",
    eta_Ka = "η_Ka (Absorption Rate)"
  ))

p6 <- ggplot(eta_data, aes(x = ETA, fill = Parameter)) +
  geom_histogram(bins = 25, alpha = 0.8, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  facet_wrap(~Parameter, scales = "free_x", ncol = 3) +
  scale_fill_manual(values = c("#3498DB", "#E74C3C", "#2ECC71")) +
  labs(
    title    = "Interindividual Variability (IIV) — Eta Distributions",
    subtitle = "Centered near 0 confirms proper random effects simulation",
    x        = "Eta Value (log-scale)",
    y        = "Count"
  ) +
  theme_pk() + theme(legend.position = "none")

ggsave("plots/06_eta_distributions.png", p6, width = 12, height = 4, dpi = 300)

cat("\n=== EDA Complete ===\n")
cat("Plots saved to plots/\n")
cat("  01_spaghetti_plot.png\n")
cat("  02_obs_vs_pred.png\n")
cat("  03_covariate_distributions.png\n")
cat("  04_covariate_pk_relationships.png\n")
cat("  05_pk_exposure_metrics.png\n")
cat("  06_eta_distributions.png\n")

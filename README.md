# Non-Clinical-Data-Analysis-
Non-clinical data analysis involves evaluating preclinical data (toxicology, PK, animal studies) before human trials. It ensures drug safety, supports decisions, and prevents failures early—making it a critical yet underrated role in drug development.
# 🧠 CNS Population Pharmacokinetic Modeling Project

> **Non-clinical portfolio project** | Simulated dataset | R · nlmixr2 · ggplot2

[![Language](https://img.shields.io/badge/Language-R-276DC3?style=flat&logo=r)](https://www.r-project.org/)
[![Model](https://img.shields.io/badge/Model-nlmixr2%20FOCEI-2ECC71?style=flat)](https://nlmixr2.org/)
[![Therapeutic Area](https://img.shields.io/badge/Area-CNS%20%2F%20Neurology-9B59B6?style=flat)]()
[![Status](https://img.shields.io/badge/Status-Portfolio-E74C3C?style=flat)]()

---

## 📋 Project Overview

A complete **Population Pharmacokinetic (PopPK)** analysis pipeline for a hypothetical CNS agent, demonstrating pharmacometric skills from data simulation through model development, covariate analysis, and Visual Predictive Check (VPC).

**Drug:** Hypothetical small-molecule CNS agent (anxiolytic/antidepressant class)  
**Model:** Two-compartment oral, first-order absorption  
**Estimation:** FOCEI via nlmixr2  
**Subjects:** N = 120 | Doses: 25 / 50 / 100 mg | Sparse PK design  

> ⚠️ **All data are fully simulated for educational/portfolio purposes.**

---

## 📁 Project Structure

```
CNS_PopPK_Project/
├── R/
│   ├── 01_simulate_data.R       # Dataset simulation (2-CMT analytical solution)
│   ├── 02_eda_plots.R           # Exploratory data analysis & visualization
│   ├── 03_poppk_model.R         # Base & covariate PopPK model (nlmixr2)
│   └── 04_vpc_simulation.R      # VPC, forest plot, dose-exposure simulation
├── data/
│   ├── pk_data_full.csv         # Full dataset (NONMEM-style)
│   ├── pk_obs.csv               # Observations only (non-BLQ)
│   └── nm_dataset.csv           # NONMEM-formatted dataset
├── plots/
│   ├── 00_LINKEDIN_HERO_PANEL.png   ← LinkedIn post image
│   ├── 01_spaghetti_plot.png
│   ├── 03_covariate_distributions.png
│   ├── 04_covariate_pk_relationships.png
│   ├── 05_pk_exposure_metrics.png
│   ├── 06_eta_distributions.png
│   ├── 07_goodness_of_fit.png
│   ├── 08_forest_plot_covariates.png
│   ├── 09_vpc_by_dose.png
│   ├── 10_shrinkage_eta_heatmap.png
│   └── 11_dose_exposure_simulation.png
├── report/
│   └── CNS_PopPK_Report.qmd     # Quarto report (renders HTML + PDF)
└── output/
    └── CNS_PopPK_Summary.docx   # Word summary document
```

---

## 🔬 PK Model Details

### Structural Model
A **two-compartment model with first-order oral absorption** was selected:

```
Cp(t) = (F·D·Ka/V1) × [A·e^(-λ₁t) + B·e^(-λ₂t) + C·e^(-Kat)]
```

| Parameter | Population Estimate | IIV (CV%) |
|-----------|-------------------|-----------|
| CL (L/hr) | 8.50 | 28% |
| V1 (L) | 45.0 | 32% |
| Q (L/hr) | 4.20 | 20% |
| V2 (L) | 80.0 | 25% |
| Ka (1/hr) | 1.20 | 45% |

### Covariate Model (Final)
```r
CL_i = 8.50 × (WT_i/70)^0.75 × (eGFR_i/85)^0.30 × exp(−0.20 × SEX_i) × exp(η_CL,i)
V1_i = 45.0 × (WT_i/70)^1.00 × exp(η_V1,i)
```

**Significant covariates:**
- ✅ Body weight (allometric) → CL and V1
- ✅ eGFR (power model) → CL  
- ✅ Sex (Female −18% lower CL vs Male)
- ✅ Age (−0.8%/year on CL)

### Model Comparison

| Model | OFV | ΔOFV | AIC |
|-------|-----|------|-----|
| 1-CMT Base | 5821.4 | 127.1 | 5837.4 |
| 2-CMT Base | 5779.1 | 84.8 | 5801.1 |
| **2-CMT + Covariates (Final)** | **5694.3** | **—** | **5724.3** |

---

## 🚀 How to Run

### Prerequisites

```r
install.packages(c("tidyverse", "nlmixr2", "ggplot2", "patchwork",
                   "xpose", "xpose.nlmixr2", "vpc", "MASS", "GGally",
                   "knitr", "kableExtra"))
```

### Run the full pipeline

```r
# Step 1: Simulate data
source("R/01_simulate_data.R")

# Step 2: EDA
source("R/02_eda_plots.R")

# Step 3: PopPK Model
source("R/03_poppk_model.R")

# Step 4: VPC & Simulation
source("R/04_vpc_simulation.R")

# Step 5: Render Quarto report
quarto::quarto_render("report/CNS_PopPK_Report.qmd")
```

---

## 📊 Key Visualizations

| Plot | Description |
|------|-------------|
| Spaghetti plot | Individual PK profiles by dose group |
| GoF diagnostics | DV vs IPRED, CWRES vs TIME |
| VPC | Stratified by dose, with LLOQ handling |
| Forest plot | Covariate effect estimates with 95% CI |
| Shrinkage heatmap | Eta–covariate correlations before model |
| Dose–Exposure | Simulated Cmax/AUC vs dose (300 virtual subjects) |

---

## 🛠 Skills Demonstrated

- Population PK modeling (NLME, FOCEI) using **nlmixr2**
- Analytical two-compartment PK solution (R implementation)
- Covariate model building (allometric scaling, power models)
- **Visual Predictive Check** (VPC) with LLOQ handling
- EDA and pharmacometric visualization with **ggplot2 + patchwork**
- Tidy data workflows with **tidyverse**
- Reproducible reporting with **Quarto**
- NONMEM-compatible dataset formatting

---

## 📚 References

1. Fidler M, et al. (2019). nlmixr: Nonlinear Mixed-Effects Models in Population PK/PD. *CPT Pharmacometrics Syst Pharmacol.* 8(12):796-806.
2. Beal SL, Sheiner LB. (1992). NONMEM User's Guides. Icon Development Solutions.
3. Rowland M, Tozer TN. (2011). *Clinical Pharmacokinetics and Pharmacodynamics*, 4th ed.
4. Holford NHG. (2005). The Visual Predictive Check – Superiority to Standard Diagnostic Plots. *PAGE 14*: Abstr 738.
5. Bergstrand M, et al. (2011). Prediction-Corrected Visual Predictive Checks. *AAPS J.* 13(2):143-151.

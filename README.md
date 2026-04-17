# Fiscal SVAR: Replication and Sensitivity Analysis

This repository replicates and extends the results from:

**Gregory, McNeil, and Smith (2023)**
*“US Fiscal Policy Shocks: Proxy-SVAR Overidentification via GMM”*

---

## 📌 Overview

This project combines **replication** and **original empirical analysis** of fiscal shock identification using Proxy-SVAR methods.

* First, it replicates the main empirical results and Monte Carlo experiments from the original paper.
* Second, it introduces a **design-based sensitivity analysis** to evaluate how estimation results change across different empirical environments.

The core focus is to compare:

* **Just-identified IV (Proxy-SVAR)**
* **Overidentified GMM (with covariance restrictions)**
* **Local Projections (LP-IV)**

The results show that **GMM improves statistical precision and stabilizes estimates**, especially under weak instruments.

---

## ✨ Original Contribution

Beyond replication, this repository includes several extensions:

### 1. Sensitivity to Instrument Quality

Three scenarios are analyzed:

* **Baseline**: original instruments
* **Trimmed**: instruments set to missing before 1980
* **Noisy**: instruments contaminated with Gaussian noise

This allows a controlled evaluation of **weak identification**.

---

### 2. Efficiency Analysis (IV vs GMM)

The project computes **relative efficiency gains**:

[
Gain = \frac{SE_{IV} - SE_{GMM}}{SE_{IV}}
]

Findings:

* Gains are **larger for government spending shocks**
* Gains **increase when instruments weaken**
* GMM acts as a **stabilizing estimator**

---

### 3. Local Projections Comparison (LP-IV)

An alternative estimation using **local projections** is implemented.

Result:

* Similar point estimates
* **Much wider confidence intervals** → loss of precision

---

### 4. Meta-Regression

A simple meta-analysis relates efficiency gains to:

* Type of shock (spending vs taxes)
* Instrument quality (trimmed / noisy)

This provides a structured summary of the sensitivity results.

---

## 📂 Repository Structure

```
fiscal-svar-replication/
│
├── main.R                      # Original replication
├── main_analysis.R             # Baseline + LP + IV vs GMM comparison
├── stress_test_analysis.R      # Sensitivity (baseline / trim / noise)
├── monteCarlo.R                # Monte Carlo simulations
├── data.csv
├── install_packages.R (optional)
├── README.md
│
├── TablesAndGraphs/
│   └── (Replication outputs: tables and figures)
│
└── Aux_files/
    ├── functions.R
    ├── makeInstData.R
    ├── getStartingValues.R
    ├── estGMM.R
    ├── estGMM_MC.R
    ├── getIRF.R
    ├── getShocks.R
    ├── asymp_ci.R
    ├── makeFigures.R
    └── getTheta.R
```

---

## 📊 Data

The dataset includes:

* **Y**: log real GDP per capita
* **G**: log government spending per capita
* **T**: log tax revenue per capita
* **GSHK_R**: spending shock instrument (Ramey & Zubairy)
* **T90_MMO**: tax instrument (Mertens & Montiel Olea)
* **TFP shock**: FRBNY DSGE-based instrument
* **TBILL3**: short-term interest rate

---

## ⚙️ Requirements

Required R packages:

```
xts, pracma, readxl, vars, MASS, matrixcalc, dynlm, gmm,
ggplot2, scales, gridExtra, AER, expm, gdata, cowplot,
sandwich, lmtest, latex2exp, forcats
```

Optional:

```r
source("install_packages.R")
```

---

## ▶️ How to Run

### 1. Replication

```r
source("main.R")
```

---

### 2. Main Analysis (Recommended)

```r
source("main_analysis.R")
```

This script:

* Replicates baseline results
* Compares IV vs GMM
* Implements LP-IV
* Produces comparison tables and plots

---

### 3. Sensitivity / Stress Test

```r
source("stress_test_analysis.R")
```

This script:

* Runs **baseline / trimmed / noisy** scenarios
* Computes first-stage strength
* Evaluates efficiency gains
* Generates summary tables and plots

---

### 4. Monte Carlo

```r
source("monteCarlo.R")
```

---

## ⚠️ Monte Carlo Note

Default:

```r
mc = 5
```

Paper:

```r
mc = 5000
```

Use higher values for full replication (computationally intensive).

---

## 📄 Project Write-Up

The full project write-up is included in this repository and provides:

* Model derivation
* Data description
* Empirical design
* Results and interpretation

Main finding:

> Overidentification via GMM improves precision without altering point estimates and plays a key stabilizing role under weak identification.

---

## 🧠 Notes

* All paths should be set to **relative paths** for portability
* Some auxiliary functions come from the original authors
* Sensitivity analysis is implemented independently

---

## 📚 Reference

Gregory, A. W., McNeil, J., & Smith, G. W. (2023)
*US Fiscal Policy Shocks: Proxy-SVAR Overidentification via GMM*

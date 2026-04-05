# Fiscal SVAR Replication

This repository replicates the results from:

**Gregory, McNeil, and Smith (2023)**
*“US Fiscal Policy Shocks: Proxy-SVAR Overidentification via GMM”*

---

## 📌 Overview

This project reproduces the main empirical results and Monte Carlo experiments from the paper. The authors study fiscal policy shocks in the United States using a structural VAR (SVAR) identified with external instruments, and propose an extension using GMM that imposes orthogonality (uncorrelated shocks).

The repository includes:

* Replication of the empirical application (impulse responses, multipliers, tables, figures)
* Replication of the Monte Carlo simulations used in the paper

---

## 📂 Repository Structure

```
fiscal-svar-replication/
│
├── main.R
├── monteCarlo.R
├── data.csv
├── install_packages.R (optional)
├── README.md
│
├── TablesAndGraphs/
│   └── (Replications of the original paper’s tables and figures)
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

The file `data.csv` contains the following variables:

1. **Y**: log of real per capita GDP
2. **G**: log of real per capita federal government spending
3. **T**: log of real per capita federal tax revenue
4. **GSHK_R**: instrument for government spending shocks (Ramey and Zubairy, 2018)
5. **T90_MMO**: instrument for tax revenue shocks (Mertens and Montiel Olea, 2018)
6. **TFPSHK_FRBNY_P**: instrument for non-fiscal shocks (NY Fed DSGE model)
7. **TBILL3**: interest rate on three-month Treasury bills

---

## 📁 Additional Folder: TablesAndGraphs

The folder `TablesAndGraphs/` contains the replication outputs corresponding to the original paper’s results. Specifically, it includes:

* Replicated tables presented in the paper
* Replicated figures (e.g., impulse response functions and related plots)

This folder is intended to provide a direct comparison between the reproduced results and those reported by the authors.

---

## ⚙️ Requirements

Before running the code, install the required R packages:

```
xts, pracma, readxl, vars, MASS, matrixcalc, dynlm, gmm,
ggplot2, scales, gridExtra, AER, expm, gdata, cowplot,
sandwich, lmtest, latex2exp, forcats
```

You can optionally run:

```r
source("install_packages.R")
```

---

## ▶️ How to Run

### 1. Set working directory

Make sure your working directory is the root of the repository:

```r
setwd("path/to/fiscal-svar-replication")
```

---

### 2. Run empirical replication

```r
source("main.R")
```

This will:

* Estimate the SVAR
* Compute impulse response functions (IRFs)
* Generate fiscal multipliers
* Reproduce tables and figures from the paper

---

### 3. Run Monte Carlo simulations

```r
source("monteCarlo.R")
```

This will:

* Simulate data according to the paper’s DGP
* Estimate SVAR models under different identification strategies
* Reproduce the distributions shown in Figure 2

---

## ⚠️ Important Note on Monte Carlo Simulations

In this repository, the parameter:

```r
mc = 5
```

was used in `monteCarlo.R` to reduce computation time during testing.

However, the authors use:

```r
mc = 5000
```

Running the simulation with `mc = 5000` fully reproduces the results and figures reported in the paper. The smaller value is only for faster execution.

---

## 🧠 Notes

* All file paths were converted to **relative paths** to ensure reproducibility across different systems.
* The code relies on auxiliary scripts in the `Aux_files/` directory provided by the authors.
* Monte Carlo simulations can be computationally intensive depending on the number of replications.

---

## 📚 Reference

Gregory, A. W., McNeil, J., & Smith, G. W. (2023).
*US Fiscal Policy Shocks: Proxy-SVAR Overidentification via GMM*

---

## 👨‍💻 Author

Replication conducted as part of a course project.

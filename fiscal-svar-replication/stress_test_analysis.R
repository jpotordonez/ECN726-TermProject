# ============================================================
# STRESS TEST: BASE vs TRIM vs NOISE FOR ALL 3 INSTRUMENTS
# ============================================================

library(xts)
library(vars)
library(gmm)
library(AER)
library(lmtest)
library(sandwich)
library(ggplot2)
library(cowplot)
library(dynlm)

lag <- stats::lag
Sys.setenv(TZ = "GMT")

base_dir <- "C:/Users/JuanPabloOrdoñez/Desktop/fiscal-svar-replication"

# Load helper functions
source(file.path(base_dir, "Aux_files/functions.R"))

# -----------------------------
# Parameters
# -----------------------------
p <- 4
irflen <- 21
rng <- "1948/2019"

var.list <- c("G", "T", "Y")
names(var.list) <- c("G", "T", "Y")
N <- length(var.list)

# Baseline instrument names
G.inst0 <- "GSHK_R"
T.inst0 <- "T90_MMO"
Y.inst0 <- "TFPSHK_FRBNY_P_R4"

shks.id <- c("G", "T", "Y")
cov.include <- c("GT", "GY", "TY")

pre.white <- FALSE
demean.inst <- FALSE
use.nwhac <- TRUE
nwhac.lag <- 4
est.method <- "BFGS"

# -----------------------------
# Load data and estimate common VAR once
# -----------------------------
dta <- read.csv(file.path(base_dir, "data.csv"))
df0 <- xts(dta[, 2:ncol(dta)], as.Date(dta$DATE))

df0$t <- seq_len(nrow(df0))
df0$t2 <- df0$t^2
df0$D75 <- ifelse(index(df0) == as.Date("1975-04-01"), 1, 0)
df0$D75L1 <- lag(df0$D75, 1)
df0$D75L2 <- lag(df0$D75, 2)
df0$D75L3 <- lag(df0$D75, 3)
df0$D75L4 <- lag(df0$D75, 4)

df0$LG <- log(df0$G)
df0$LT <- log(df0$T)
df0$LY <- log(df0$Y)

# Residualize the FRBNY shock on its own lags and VAR lags
df0$TFPSHK_FRBNY_P_R4 <- NA
df0["1961/", "TFPSHK_FRBNY_P_R4"] <- residuals(
  dynlm(TFPSHK_FRBNY_P ~ L(TFPSHK_FRBNY_P, 1:4) +
          L(LY, 1:p) + L(LT, 1:p) + L(LG, 1:p),
        data = as.zoo(df0))
)

# Reduced-form VAR, estimated once and reused
endog.var <- log(df0[, var.list])[rng]
exog.var  <- df0[rng, c("t", "t2", "D75", "D75L1", "D75L2", "D75L3", "D75L4")]

BP.var <- VAR(endog.var, p = p, exogen = exog.var)
U <- resid(BP.var)

B <- NULL
D <- NULL
for (v in var.list) {
  eval(parse(text = paste0("B = rbind(B, BP.var$varresult$", v,
                           "$coefficients[1:(N*p)])")))
  eval(parse(text = paste0("D = rbind(D, BP.var$varresult$", v,
                           "$coefficients[c('const', colnames(exog.var))])")))
}
if (p > 1) {
  B <- rbind(B, cbind(eye(N * (p - 1)), zeros(N * (p - 1), N)))
}

var.cov <- t(U) %*% U / nrow(U)

# ============================================================
# Helper functions
# ============================================================

theta_se_matrix <- function(se_vec, var_names) {
  out <- matrix(NA_real_, nrow = length(var_names), ncol = length(var_names),
                dimnames = list(var_names, var_names))
  idx <- 0
  for (i in var_names) {
    for (j in var_names) {
      if (i == j) next
      idx <- idx + 1
      out[i, j] <- se_vec[idx]
    }
  }
  out
}

make_versioned_df <- function(kind = c("base", "trim", "noise"),
                              trim_start = as.Date("1980-01-01"),
                              noise_lambda = 0.35,
                              seed = 123) {
  kind <- match.arg(kind)
  df <- df0
  
  if (kind == "base") {
    return(list(
      df = df,
      G.inst = G.inst0,
      T.inst = T.inst0,
      Y.inst = Y.inst0
    ))
  }
  
  if (kind == "trim") {
    df$GSHK_R_trim <- df$GSHK_R
    df$T90_MMO_trim <- df$T90_MMO
    df$TFPSHK_FRBNY_P_R4_trim <- df$TFPSHK_FRBNY_P_R4
    
    ix <- index(df) < trim_start
    df$GSHK_R_trim[ix] <- NA
    df$T90_MMO_trim[ix] <- NA
    df$TFPSHK_FRBNY_P_R4_trim[ix] <- NA
    
    return(list(
      df = df,
      G.inst = "GSHK_R_trim",
      T.inst = "T90_MMO_trim",
      Y.inst = "TFPSHK_FRBNY_P_R4_trim"
    ))
  }
  
  if (kind == "noise") {
    set.seed(seed)
    
    zG <- as.numeric(df$GSHK_R)
    zT <- as.numeric(df$T90_MMO)
    zY <- as.numeric(df$TFPSHK_FRBNY_P_R4)
    
    df$GSHK_R_noise <- zG + noise_lambda * sd(zG, na.rm = TRUE) * rnorm(length(zG))
    df$T90_MMO_noise <- zT + noise_lambda * sd(zT, na.rm = TRUE) * rnorm(length(zT))
    df$TFPSHK_FRBNY_P_R4_noise <- zY + noise_lambda * sd(zY, na.rm = TRUE) * rnorm(length(zY))
    
    return(list(
      df = df,
      G.inst = "GSHK_R_noise",
      T.inst = "T90_MMO_noise",
      Y.inst = "TFPSHK_FRBNY_P_R4_noise"
    ))
  }
}

run_one_version <- function(kind = c("base", "trim", "noise"),
                            trim_start = as.Date("1980-01-01"),
                            noise_lambda = 0.35,
                            seed = 123) {
  kind <- match.arg(kind)
  
  tmp <- make_versioned_df(
    kind = kind,
    trim_start = trim_start,
    noise_lambda = noise_lambda,
    seed = seed
  )
  
  # Push version-specific objects into the global environment,
  # because the auxiliary scripts rely on globals.
  df <<- tmp$df
  G.inst <<- tmp$G.inst
  T.inst <<- tmp$T.inst
  Y.inst <<- tmp$Y.inst
  
  # Build the instrument dataset
  source(file.path(base_dir, "Aux_files/makeInstData.R"))
  
  # First-stage robust F-statistics
  osls.T <- lm(T ~ ivdat[, shks.list[["T"]]], data = ivdat)
  osls.G <- lm(G ~ ivdat[, shks.list[["G"]]], data = ivdat)
  osls.Y <- lm(Y ~ ivdat[, shks.list[["Y"]]], data = ivdat)
  
  F_T <- unname(waldtest(lm(T ~ 1, data = ivdat), osls.T,
                         vcov = vcovHC(osls.T, type = "HC0"))[2, "F"])
  F_G <- unname(waldtest(lm(G ~ 1, data = ivdat), osls.G,
                         vcov = vcovHC(osls.G, type = "HC0"))[2, "F"])
  F_Y <- unname(waldtest(lm(Y ~ 1, data = ivdat), osls.Y,
                         vcov = vcovHC(osls.Y, type = "HC0"))[2, "F"])
  
  # -------------------------
  # 1) Just-identified IV
  # -------------------------
  use.cov.res <<- FALSE
  
  source(file.path(base_dir, "Aux_files/getStartingValues.R"))
  source(file.path(base_dir, "Aux_files/estGMM.R"))
  source(file.path(base_dir, "Aux_files/getIRF.R"))
  source(file.path(base_dir, "Aux_files/getShocks.R"))
  source(file.path(base_dir, "Aux_files/asymp_ci.R"))
  
  theta_iv <- Thta.gmm
  se_iv <- theta_se_matrix(sqrt(diag(Sig_Thta)), names(var.list))
  eps_iv <- eps
  mult_iv <- MULT
  dmult_iv <- DMULT
  se_mult_iv <- se.mult
  se_dmult_iv <- se.dmult
  
  # -------------------------
  # 2) Overidentified GMM
  # -------------------------
  use.cov.res <<- TRUE
  
  source(file.path(base_dir, "Aux_files/getStartingValues.R"))
  source(file.path(base_dir, "Aux_files/estGMM.R"))
  source(file.path(base_dir, "Aux_files/getIRF.R"))
  source(file.path(base_dir, "Aux_files/getShocks.R"))
  source(file.path(base_dir, "Aux_files/asymp_ci.R"))
  
  theta_gmm <- Thta.gmm
  se_gmm <- theta_se_matrix(sqrt(diag(Sig_Thta)), names(var.list))
  eps_gmm <- eps
  mult_gmm <- MULT
  dmult_gmm <- DMULT
  se_mult_gmm <- se.mult
  se_dmult_gmm <- se.dmult
  
  # -------------------------
  # Summary statistics
  # -------------------------
  shock_cor_iv  <- cor(eps_iv)
  shock_cor_gmm <- cor(eps_gmm)
  
  summary_row <- data.frame(
    scenario = kind,
    n_obs = nrow(ivdat),
    F_G = F_G,
    F_T = F_T,
    F_Y = F_Y,
    J.stat = J.stat,
    J.pval = J.pval,
    
    thetaY_G_iv = theta_iv["Y", "G"],
    seY_G_iv = se_iv["Y", "G"],
    thetaY_G_gmm = theta_gmm["Y", "G"],
    seY_G_gmm = se_gmm["Y", "G"],
    
    thetaY_T_iv = theta_iv["Y", "T"],
    seY_T_iv = se_iv["Y", "T"],
    thetaY_T_gmm = theta_gmm["Y", "T"],
    seY_T_gmm = se_gmm["Y", "T"],
    
    dmultG_h0_iv = dmult_iv[1, "G"],
    dmultG_h0_gmm = dmult_gmm[1, "G"],
    se_dmultG_h0_iv = se_dmult_iv[1, "G"],
    se_dmultG_h0_gmm = se_dmult_gmm[1, "G"],
    
    dmultT_h0_iv = dmult_iv[1, "T"],
    dmultT_h0_gmm = dmult_gmm[1, "T"],
    se_dmultT_h0_iv = se_dmult_iv[1, "T"],
    se_dmultT_h0_gmm = se_dmult_gmm[1, "T"],
    
    rhoGT_iv = shock_cor_iv["G", "T"],
    rhoGY_iv = shock_cor_iv["G", "Y"],
    rhoTY_iv = shock_cor_iv["T", "Y"],
    
    rhoGT_gmm = shock_cor_gmm["G", "T"],
    rhoGY_gmm = shock_cor_gmm["G", "Y"],
    rhoTY_gmm = shock_cor_gmm["T", "Y"],
    
    se_ratio_YG = se_gmm["Y", "G"] / se_iv["Y", "G"],
    se_ratio_YT = se_gmm["Y", "T"] / se_iv["Y", "T"],
    se_ratio_dmultG_h0 = se_dmult_gmm[1, "G"] / se_dmult_iv[1, "G"],
    se_ratio_dmultT_h0 = se_dmult_gmm[1, "T"] / se_dmult_iv[1, "T"]
  )
  
  list(
    summary = summary_row,
    iv = list(theta = theta_iv, se = se_iv, eps = eps_iv,
              MULT = mult_iv, DMULT = dmult_iv,
              se.mult = se_mult_iv, se.dmult = se_dmult_iv),
    gmm = list(theta = theta_gmm, se = se_gmm, eps = eps_gmm,
               MULT = mult_gmm, DMULT = dmult_gmm,
               se.mult = se_mult_gmm, se.dmult = se_dmult_gmm)
  )
}

# ============================================================
# Run the three scenarios
# ============================================================

trim_cutoff <- as.Date("1980-01-01")
noise_lambda <- 0.35

res_base  <- run_one_version("base")
res_trim  <- run_one_version("trim", trim_start = trim_cutoff)
res_noise <- run_one_version("noise", noise_lambda = noise_lambda, seed = 123)

stress_table <- rbind(
  res_base$summary,
  res_trim$summary,
  res_noise$summary
)

print(stress_table)
write.csv(stress_table,
          file = file.path(base_dir, "stress_test_summary.csv"),
          row.names = FALSE)

# ============================================================
# Efficiency comparison
# ============================================================

gain_wide <- rbind(
  data.frame(
    scenario = stress_table$scenario,
    shock = "Spending",
    se_iv = stress_table$seY_G_iv,
    se_gmm = stress_table$seY_G_gmm,
    mult_iv = stress_table$thetaY_G_iv,
    mult_gmm = stress_table$thetaY_G_gmm,
    stringsAsFactors = FALSE
  ),
  data.frame(
    scenario = stress_table$scenario,
    shock = "Taxes",
    se_iv = stress_table$seY_T_iv,
    se_gmm = stress_table$seY_T_gmm,
    mult_iv = stress_table$thetaY_T_iv,
    mult_gmm = stress_table$thetaY_T_gmm,
    stringsAsFactors = FALSE
  )
)

gain_wide$gain_se <- (gain_wide$se_iv - gain_wide$se_gmm) / gain_wide$se_iv
gain_wide$tax_dummy <- ifelse(gain_wide$shock == "Taxes", 1, 0)
gain_wide$weak_instr <- ifelse(gain_wide$scenario == "noise", 1, 0)
gain_wide$trim_instr <- ifelse(gain_wide$scenario == "trim", 1, 0)

gain_wide$scenario <- factor(gain_wide$scenario, levels = c("base", "trim", "noise"))
gain_wide$shock <- factor(gain_wide$shock, levels = c("Spending", "Taxes"))

print(gain_wide)

eff_summary <- aggregate(
  cbind(gain_se) ~ shock + scenario,
  data = gain_wide,
  FUN = mean
)
print(eff_summary)

meta_reg <- lm(gain_se ~ shock + trim_instr + weak_instr, data = gain_wide)
print(summary(meta_reg))

# ============================================================
# Plot: relative efficiency gains
# ============================================================

p_gain <- ggplot(gain_wide,
                 aes(x = interaction(scenario, shock, sep = " / "),
                     y = gain_se,
                     fill = shock)) +
  geom_col(width = 0.7) +
  theme_minimal() +
  labs(
    title = "Relative Efficiency Gain from GMM",
    x = "Scenario / Shock",
    y = "Gain = (SE_IV - SE_GMM) / SE_IV",
    fill = ""
  )

print(p_gain)

# ============================================================
# Optional: descriptive regression table
# ============================================================

meta_long <- gain_wide
meta_long$eff_gain <- meta_long$gain_se

meta_reg_2 <- lm(eff_gain ~ shock + trim_instr + weak_instr, data = meta_long)
print(summary(meta_reg_2))
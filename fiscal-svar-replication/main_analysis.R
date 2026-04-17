############################################################
# LIBRARIES
############################################################

library(xts)
library(pracma)
library(readxl)
library(vars)
library(MASS)
library(matrixcalc)
library(dynlm)
library(gmm)
library(ggplot2)
library(scales)
library(grid)
library(gridExtra)
library(AER)
library(expm)
library(gdata)
library(cowplot)

lag <- stats::lag
Sys.setenv(TZ = 'GMT')

############################################################
# LOAD AUXILIARY FUNCTIONS
############################################################

source("C:/Users/JuanPabloOrdoñez/Desktop/fiscal-svar-replication/Aux_files/functions.R")

############################################################
# MODEL PARAMETERS
############################################################

p       <- 4       # Number of VAR lags
irflen  <- 21      # IRF horizon (quarters)
rng     <- '1948/2019'

var.list <- c("G","T","Y")
names(var.list) <- c("G","T","Y")

N <- length(var.list)

# External instruments
T.inst <- c("T90_MMO")
G.inst <- c("GSHK_R")
Y.inst <- c("TFPSHK_FRBNY_P_R4")

# Identified shocks
shks.id <- c("G","T","Y")

# Covariance restrictions
cov.include <- c("GT","GY","TY")

L <- length(shks.id)

# Estimation options
pre.white   <- FALSE
demean.inst <- FALSE
use.nwhac   <- TRUE
nwhac.lag   <- 4
est.method  <- "BFGS"

############################################################
# LOAD AND PREPARE DATA
############################################################

dta <- read.csv("C:/Users/JuanPabloOrdoñez/Desktop/fiscal-svar-replication/data.csv")
df  <- xts(dta[,2:ncol(dta)], as.Date(dta$DATE))

# Deterministic components
df$t  <- 1:nrow(df)
df$t2 <- df$t^2

# Dummy for 1975Q2 and lags
df$D75   <- ifelse(index(df) == '1975-04-01', 1, 0)
df$D75L1 <- lag(df$D75,1)
df$D75L2 <- lag(df$D75,2)
df$D75L3 <- lag(df$D75,3)
df$D75L4 <- lag(df$D75,4)

# Log variables
df$LG <- log(df$G)
df$LT <- log(df$T)
df$LY <- log(df$Y)

############################################################
# SUMMARY STATISTICS
############################################################

sum_vars <- df[, c("LY","LG","LT",
                   "GSHK_R","T90_MMO","TFPSHK_FRBNY_P","TBILL3")]

sum_df <- data.frame(sum_vars)

summary_stats <- data.frame(
  Variable = colnames(sum_df),
  Mean = sapply(sum_df, mean, na.rm = TRUE),
  SD   = sapply(sum_df, sd, na.rm = TRUE),
  Min  = sapply(sum_df, min, na.rm = TRUE),
  Max  = sapply(sum_df, max, na.rm = TRUE)
)

summary_stats[, -1] <- round(summary_stats[, -1], 3)
print(summary_stats)

############################################################
# PREPROCESSING: TFP SHOCK
############################################################

df$TFPSHK_FRBNY_P_R4 <- NA

df['1961/', 'TFPSHK_FRBNY_P_R4'] <- residuals(
  dynlm(TFPSHK_FRBNY_P ~ L(TFPSHK_FRBNY_P,1:4) +
          L(LY,1:p) + L(LT,1:p) + L(LG,1:p),
        data = as.zoo(df))
)

############################################################
# REDUCED-FORM VAR
############################################################

endog.var <- log(df[,var.list])[rng]
exog.var  <- df[rng, c("t","t2","D75","D75L1","D75L2","D75L3","D75L4")]

BP.var <- VAR(endog.var, p = p, exogen = exog.var)
U <- resid(BP.var)

# Construct companion matrix
B <- NULL
D <- NULL

for (v in var.list) {
  B <- rbind(B, BP.var$varresult[[v]]$coefficients[1:(N*p)])
  D <- rbind(D, BP.var$varresult[[v]]$coefficients[c('const', colnames(exog.var))])
}

if (p > 1) {
  B <- rbind(B, cbind(eye(N*(p-1)), zeros(N*(p-1), N)))
}

var.cov <- t(U) %*% U / nrow(U)

############################################################
# SVAR-IV ESTIMATION
############################################################

source("Aux_files/makeInstData.R")

# First-stage regressions
osls.T <- lm(T ~ ivdat[,shks.list[["T"]]], data = ivdat)
osls.G <- lm(G ~ ivdat[,shks.list[["G"]]], data = ivdat)
osls.Y <- lm(Y ~ ivdat[,shks.list[["Y"]]], data = ivdat)

# Robust F-tests
waldtest(lm(T~1, data=ivdat), osls.T, vcov = vcovHC(osls.T))
waldtest(lm(G~1, data=ivdat), osls.G, vcov = vcovHC(osls.G))
waldtest(lm(Y~1, data=ivdat), osls.Y, vcov = vcovHC(osls.Y))

############################################################
# JUST-IDENTIFIED MODEL (IV)
############################################################

use.cov.res <- FALSE

source("Aux_files/getStartingValues.R")
source("Aux_files/estGMM.R")
source("Aux_files/getIRF.R")
source("Aux_files/getShocks.R")
source("Aux_files/asymp_ci.R")

Thta.iv      <- Thta.gmm
Thta.se.iv   <- sqrt(diag(Sig_Thta))
MULT.iv      <- MULT
MULT.se.iv   <- se.mult

############################################################
# OVERIDENTIFIED MODEL (GMM)
############################################################

use.cov.res <- TRUE

source("Aux_files/estGMM.R")
source("Aux_files/getIRF.R")
source("Aux_files/getShocks.R")
source("Aux_files/asymp_ci.R")

Thta.se.gmm <- sqrt(diag(Sig_Thta))

############################################################
# EFFICIENCY GAINS (META ANALYSIS)
############################################################

Theta_se_iv_mat  <- matrix(Thta.se.iv,  nrow=3, byrow=TRUE)
Theta_se_gmm_mat <- matrix(Thta.se.gmm, nrow=3, byrow=TRUE)

se_iv_g  <- Theta_se_iv_mat[3,1]
se_gmm_g <- Theta_se_gmm_mat[3,1]

se_iv_t  <- Theta_se_iv_mat[3,2]
se_gmm_t <- Theta_se_gmm_mat[3,2]

gain_g <- (se_iv_g - se_gmm_g) / se_iv_g
gain_t <- (se_iv_t - se_gmm_t) / se_iv_t

meta_df <- data.frame(
  shock = c("Spending","Taxes"),
  gain  = c(gain_g, gain_t)
)

meta_df$tax_dummy <- ifelse(meta_df$shock=="Taxes",1,0)

meta_reg <- lm(gain ~ tax_dummy, data = meta_df)
summary(meta_reg)

############################################################
# LOCAL PROJECTIONS (LP-IV)
############################################################

h_target <- 10
df_lp <- as.data.frame(df)

df_lp$LY_h10 <- c(df_lp$LY[(h_target+1):nrow(df_lp)], rep(NA, h_target))

# Create lags
for (i in 1:4) {
  df_lp[[paste0("LY_L", i)]] <- c(rep(NA, i), df_lp$LY[1:(nrow(df_lp)-i)])
  df_lp[[paste0("LG_L", i)]] <- c(rep(NA, i), df_lp$LG[1:(nrow(df_lp)-i)])
  df_lp[[paste0("LT_L", i)]] <- c(rep(NA, i), df_lp$LT[1:(nrow(df_lp)-i)])
}

# LP-IV regression
form_lp <- LY_h10 ~ LG +
  LY_L1 + LY_L2 + LY_L3 + LY_L4 +
  LG_L1 + LG_L2 + LG_L3 + LG_L4 +
  LT_L1 + LT_L2 + LT_L3 + LT_L4 |
  GSHK_R +
  LY_L1 + LY_L2 + LY_L3 + LY_L4 +
  LG_L1 + LG_L2 + LG_L3 + LG_L4 +
  LT_L1 + LT_L2 + LT_L3 + LT_L4

modelo_lp <- ivreg(form_lp, data = df_lp)

lp_coef <- coef(modelo_lp)["LG"]
lp_se   <- sqrt(diag(vcovHAC(modelo_lp)))["LG"]

############################################################
# FINAL COMPARISON TABLE
############################################################

comparison_table <- data.frame(
  Method = c("SVAR-IV",
             "SVAR-GMM",
             "LP-IV"),
  Multiplier = c(MULT.iv[10,"G"], MULT[10,"G"], lp_coef),
  Std_Error  = c(MULT.se.iv[10,"G"], se.mult[10,"G"], lp_se)
)

comparison_table$Rel_Efficiency <- comparison_table$Std_Error /
  comparison_table$Std_Error[2]

print(round(comparison_table, 4))

############################################################
# PLOT: CONFIDENCE INTERVALS
############################################################

comparison_table$lower <- comparison_table$Multiplier -
  1.96 * comparison_table$Std_Error

comparison_table$upper <- comparison_table$Multiplier +
  1.96 * comparison_table$Std_Error

comparison_table$Method <- factor(
  comparison_table$Method,
  levels = rev(c("SVAR-IV","SVAR-GMM","LP-IV"))
)

ggplot(comparison_table,
       aes(x = Multiplier, y = Method, color = Method)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.3, size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Methodological Sensitivity: SVAR vs. Local Projections",
       subtitle = "Cumulative Spending Multiplier (h = 10)",
       x = "Multiplier", y = "") +
  theme(legend.position = "none")
#import libraries
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
Sys.setenv(TZ='GMT')

#load functions
source("Aux_files/functions.R")

#set parameters and options
p = 4 #number of lags
irflen = 21 #maximum impulse response horizon (quarters)
rng = '1948/2019' #estimation range
var.list = c("G","T","Y")
names(var.list) = c("G","T","Y")#name variables (must be a single character)
N = length(var.list) #number of equations in VAR
T.inst = c("T90_MMO") 
G.inst = c("GSHK_R")
Y.inst = c("TFPSHK_FRBNY_P_R4") 

shks.id = c("G","T","Y") #shocks identified by instruments, in same order and names as var.list
cov.include = c("GT","GY",
                "TY") #which covariance restrictions should apply (order important, should match var.list)
L = length(shks.id) #number of shocks we identify
pre.white = FALSE #optionally regress instruments off p lags of the endogenous VAR variables
demean.inst = FALSE #optionally demean instruments (centers non-zero observations)
use.nwhac = TRUE #option to use Newey-West HAC
nwhac.lag = 4
est.method = "BFGS" #choice of estimation method. Can be either BFGS or PORT
  
#load macro data and convert to time series object (xts)  
dta = read.csv("data.csv")
df = xts(dta[,2:ncol(dta)],as.Date(dta$DATE)) #convert to xts dataframe
df$t = c(1:nrow(df)) #add deterministic trends
df$t2 = df$t^2
df$D75 = ifelse(index(df)=='1975-04-01',1,0) #dummy for 1975Q2
df$D75L1 = lag(df$D75,1) #lags of 1975Q2 dummy
df$D75L2 = lag(df$D75,2)
df$D75L3 = lag(df$D75,3)
df$D75L4 = lag(df$D75,4)

#regress FRBNY shock off own lags
df$TFPSHK_FRBNY_P_R4 = NA

df$LG = log(df$G)
df$LT = log(df$T)
df$LY = log(df$Y)

df['1961/','TFPSHK_FRBNY_P_R4'] = residuals(dynlm(TFPSHK_FRBNY_P~L(TFPSHK_FRBNY_P,1:4)+
                                                    L(LY,1:p)+L(LT,1:p)+L(LG,1:p),
                                            data=as.zoo(df)))

#estimate reduced-form VAR 
#first specify endogenous and exogenous variables
endog.var = log(df[,var.list])[rng]
exog.var = df[rng,c("t","t2","D75","D75L1","D75L2","D75L3","D75L4")]

#estimate VAR
BP.var = VAR(endog.var,p=p,exogen=exog.var)
U = resid(BP.var) #save residuals

#put autoregressive coefficients in a matrix B
B = NULL
D = NULL
for (v in var.list){
  eval(parse(text=paste0("B = rbind(B,BP.var$varresult$",v,
                         "$coefficients[1:(N*p)])")))
  eval(parse(text=paste0("D = rbind(D,BP.var$varresult$",v,
                         "$coefficients[c('const',colnames(exog.var))])")))
}
if (p>1){
  B = rbind(B,cbind(eye(N*(p-1)),zeros(N*(p-1),N)))
}

var.cov = t(U)%*%U/nrow(U)

#SVAR IV estimation
source("Aux_files/makeInstData.R")


#first stage regressions (to get first-stage F-statistics)
osls.T = lm(T~ivdat[,shks.list[["T"]]],data=ivdat)
osls.G = lm(G~ivdat[,shks.list[["G"]]],data=ivdat)
osls.Y = lm(Y~ivdat[,shks.list[["Y"]]],data=ivdat)

#robust F-tests 
waldtest(lm(T~1,data=ivdat),osls.T,vcov = vcovHC(osls.T, type = "HC0"))
waldtest(lm(G~1,data=ivdat),osls.G,vcov = vcovHC(osls.G, type = "HC0"))
waldtest(lm(Y~1,data=ivdat),osls.Y,vcov = vcovHC(osls.Y, type = "HC0"))

# First estimate the just-identified model
use.cov.res = FALSE

#use IV/OLS estimates as starting values
source("Aux_files/getStartingValues.R")

# estimate coefficients in Theta
source("Aux_files/estGMM.R")
#impulse response functions
source("Aux_files/getIRF.R")
#get structural shocks
source("Aux_files/getShocks.R")
#calculate confidence intervals
source("Aux_files/asymp_ci.R")

#rename some of this output
Thta.iv = Thta.gmm
Thta.se.iv = sqrt(diag(Sig_Thta)) #standard errors
IRF.iv = IRF
IRF.se.iv = se.irf
MULT.iv = MULT
MULT.se.iv = se.mult
DMULT.iv = DMULT
DMULT.se.iv = se.dmult
eps.iv = eps
eps.sw.iv = eps.sw

# Now estimate the overidentified model
use.cov.res = TRUE
source("Aux_files/estGMM.R")
source("Aux_files/getIRF.R")
source("Aux_files/getShocks.R")
source("Aux_files/asymp_ci.R")

# J-statistic and p-value
J.stat
J.pval

#correlations between just-identified shocks
cor(eps.iv)

#or, using the Stock and Watson partial-identification method (see appendix)
cor(eps.sw.iv)

#Table 1 estimates 
 # Top panel:
Thta.iv #point estimates
Thta.se.iv #standard errors

# Middle panel:
Thta.gmm
sqrt(diag(Sig_Thta))

#Table 2 correlations 
 # Top panel
t(cor(xts(eps.iv,as.Date(rownames(eps.iv)))[index(ivdat)],ivdat[,(N+1):ncol(ivdat)]))
 # bottom panel
t(cor(xts(eps,as.Date(rownames(eps)))[index(ivdat)],ivdat[,(N+1):ncol(ivdat)]))

# produce Figures from paper and appendix:
source("Aux_files/makeFigures.R")

# multipliers (Figure 1):
fig.1

# estimated shocks (Figure A1 from appendix)
fig.A1

# IRFs: (Figures A2A, A2B, and A2C from appendix)
fig.A2A 
fig.A2B 
fig.A2C



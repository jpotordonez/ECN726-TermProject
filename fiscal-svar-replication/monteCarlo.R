#import libraries
library(xts)
library(pracma)
library(vars)
library(MASS)
library(matrixcalc)
library(gmm)
library(ggplot2)
library(AER)
library(expm)
library(forcats)
library(cowplot)
library(latex2exp)
library(gdata)
Sys.setenv(TZ='GMT')

#load functions
source("Aux_files/functions.R")

#specify Monte Carlo experiment
version = "V1" # V1: benchmark Monte Carlo experiment in the main text
               # V2: one instrument violates exogeneity restriction

#set parameters
T = 275 #sample size 
N = 3 #number of endogenous variables
p = 0 #number of lags in VAR
mc = 5000 #number of Monte Carlo replications
nwhac.lag = 4
A.inst = "A.inst" #order of variables and instruments
B.inst = "B.inst"
C.inst = "C.inst"
var.list = c("A","B","C")
names(var.list) = c("A","B","C") #name variables
shks.id = c("A","B","C") #identified shocks, in same order as var.list
#list of all instruments for each identified shocks
shks.list = list(A="A.inst",B="B.inst",C="C.inst") 
cov.include = c("AB","AC",
                "BC") #which covariance restrictions should apply (order important, should match var.list)
L = length(shks.id) #number of identified shocks
cfs.iv = matrix(NA,mc,L*(N-1)) #storage for elements of Theta
cfs.cov = matrix(NA,mc,L*(N-1))
cor.iv = matrix(NA,mc,L*(L+1)/2)
cor.cov = matrix(NA,mc,L*(L+1)/2)
J.mc = matrix(NA,mc,2,dimnames = list(seq(mc),c("iv","cov")))
F.mc = matrix(NA,mc,L,dimnames=list(seq(mc),shks.id))

#--- DGP ---#
Theta = matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),N,N)
Sig.eps = diag(c(1,1,1))
Sig.z = 0.316 #standard deviation of instrument error terms
if (version == "V2"){
  Phi = 0.2*matrix(c(1,0,0.25,0,1,0,0,0,1),N,N)
}else{
  Phi = 0.2*eye(N)
}


#monte carlo draws
for (m in seq(mc)){
  #generate samples
  eps = t(mvrnorm(T,zeros(N,1),Sig.eps)) #structural shocks
  y = Theta%*%eps
  v = t(mvrnorm(T,zeros(N,1),eye(N))) #errors for instruments
  z = Phi%*%eps + Sig.z*v
  rownames(z) = unlist(shks.list)

  rownames(y) = names(var.list)
  u = t(y) 
  var.cov = t(u)%*%u/nrow(u)

  ivdat = xts(cbind(u,t(z[,(p+1):T])),as.Date(seq((p+1):T)))
  
  F.mc[m,"A"] = summary(lm(A~A.inst,data=ivdat))$fstatistic["value"]
  F.mc[m,"B"] = summary(lm(B~B.inst,data=ivdat))$fstatistic["value"]
  F.mc[m,"C"] = summary(lm(C~C.inst,data=ivdat))$fstatistic["value"]
  
  #SVAR-IV estimates
  use.cov.res = FALSE
  source("Aux_files/getStartingValues.R")
  source("Aux_files/estGMM_MC.R")
  cfs.iv[m,] = cfs
  cor.iv[m,] = cor(eps)[upper.tri(eye(L),diag=TRUE)]
  J.mc[m,"iv"] = J.stat
  
  #add covariance restrictions
  use.cov.res = TRUE
  source("Aux_files/getStartingValues.R")
  source("Aux_files/estGMM_MC.R")
  cfs.cov[m,] = cfs
  cor.cov[m,] = cor(eps)[upper.tri(eye(L),diag=TRUE)]
  J.mc[m,"cov"] = J.stat
}

#make figures

fig_thm = theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background=element_blank(),
                panel.border =element_rect(colour='black',fill=NA),
                axis.text=element_text(color='black'),
                axis.ticks=element_line(color='black'),
                plot.title = element_text(face="bold",hjust=0.5),
                text=element_text(family = "Times"),
                legend.title=element_blank()
)


jstat.df = data.frame(method=rep(c("COV","CHI(3)"),each=mc),
                      estimate=c(J.mc[,"cov"],rchisq(mc,3)))

# rejection frequency:
sum(J.mc[,"cov"]>qchisq(0.9,3))/mc

# fig 2.1
ggplot(data=jstat.df, aes(x=estimate,color=method)) + 
  stat_density(geom="line",position="identity",bw = c("nrd"))+
  fig_thm+
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_manual( values = c("grey","red"))+
  labs(y='Density',x=TeX("$\\mathit{J}$"),title = TeX("$\\mathit{J}$-statistic"))+
  theme(legend.position = "none")+
  labs(y='',x='',title = TeX("$J$-statistic and $\\chi^2(3)"))

theta12 = data.frame(method=rep(c("IV","COV"),each=mc),
                    estimate=c(cfs.iv[,1],cfs.cov[,1]))

# fig 2.2
ggplot(data=theta12, aes(x=estimate,color=method)) + 
  stat_density(geom="line",position="identity",bw = c("nrd"))+
  scale_colour_manual( values = c("red","black"))+
  geom_vline(xintercept=Theta[1,2])+fig_thm+
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")+
  labs(y='',x='',title = TeX("$\\Theta_{1,2}$"))


omega12 = data.frame(method=rep(c("IV","COV"),each=mc),
                     estimate=c(cor.iv[,2],cor.cov[,2]))

# fig 2.3
ggplot(data=omega12, aes(x=estimate,color=method)) + 
  stat_density(geom="line",position="identity",bw = c("nrd"))+
  scale_colour_manual( values = c("red","black"))+
  geom_vline(xintercept=Sig.eps[1,2])+fig_thm+
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")+
  labs(y='',x='',title = TeX("$\\rho(\\epsilon_{1,t},\\epsilon_{2,t})$"))



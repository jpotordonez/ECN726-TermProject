#estimate model by GMM and calculate Newey-West HAC

n.mom = size(g.gmm(svals,ivdat),2) #number of moments
n.cfs = size(svals,2) #number of coefficients
#first-stage estimates
gmm.stp1 = optim(svals, gmm.obj, gr=gmm.obj.grad,
                 x=ivdat, 
                 mom.fun=g.gmm,wmat=eye(n.mom),
                 method="BFGS",control = list(maxit = 10000))

#initialize stp2 for algorithm. This is only used in stopping rule, no effect on solution
cfs1=gmm.stp1$par
cfs2=cfs1 + 1000
i=0
while (sqrt(sum((cfs1-cfs2)^2))>0.0001){
  i=i+1
  cfs1=gmm.stp1$par
  grad.at.sol = g.gmm(gmm.stp1$par,ivdat)
  #get weighting matrix
  Sig.nw = 0.5*t(grad.at.sol)%*%grad.at.sol
  if (nwhac.lag>0){
    for (j in seq(1,nwhac.lag)){
      Sig.nw = Sig.nw + (1-j/(nwhac.lag+1))*
        t(grad.at.sol[(j+1):nrow(grad.at.sol),])%*%
        grad.at.sol[1:(nrow(grad.at.sol)-j),]
    }
  }
  Sig.nw = (Sig.nw + t(Sig.nw))/nrow(grad.at.sol)
  
  #second-stage estimates
  gmm.stp2 = optim(cfs1, gmm.obj, gr=gmm.obj.grad,
                   x=ivdat, 
                   mom.fun=g.gmm,
                   wmat=inv(Sig.nw),
                   method="BFGS",control = list(maxit = 10000))
  cfs2=gmm.stp2$par
  gmm.stp1=gmm.stp2
}
cfs = gmm.stp2$par

J.stat = gmm.obj(cfs,ivdat,g.gmm,inv(Sig.nw))*nrow(ivdat)
J.df = n.mom-n.cfs
J.pval = pchisq(J.stat,J.df,lower.tail = FALSE)

source("Aux_files/getTheta.R")
eps = t(inv(Thta.gmm)%*%t(u))



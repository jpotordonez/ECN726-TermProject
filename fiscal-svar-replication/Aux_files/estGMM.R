#estimate model by GMM and calculate Newey-West HAC

n.mom = size(g.gmm(svals,ivdat),2) #number of moments
n.cfs = size(svals,2) #number of coefficients
#first-stage estimates
if (est.method=="BFGS"){
  gmm.stp1 = optim(svals, gmm.obj, gr=gmm.obj.grad,
                   x=ivdat[index(na.omit(ivdat))], 
                   mom.fun=g.gmm,wmat=eye(n.mom),
                   method="BFGS",control = list(maxit = 10000))
}else if (est.method=="PORT"){
  gmm.stp1=nlminb(svals,gmm.obj,x=ivdat[index(na.omit(ivdat))], 
                 gradient=gmm.obj.grad,
                 mom.fun=g.gmm,wmat=eye(n.mom),
                 control = list(maxit = 10000))
}else{
  stop("Error: estimation method must be one of 'BFGS' or 'PORT'")
}

#initialize stp2 for algorithm. This is only used in stopping rule, no effect on solution
cfs1=gmm.stp1$par
cfs2=cfs1 + 1000
i=0
while (sqrt(sum((cfs1-cfs2)^2))>0.0001){
  i=i+1
  cfs1=gmm.stp1$par
  grad.at.sol = g.gmm(gmm.stp1$par,ivdat[index(na.omit(ivdat))])
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
  if (est.method=="BFGS"){
    gmm.stp2 = optim(cfs1, gmm.obj, gr=gmm.obj.grad,
                     x=ivdat[index(na.omit(ivdat))], 
                     mom.fun=g.gmm,
                     wmat=inv(Sig.nw),
                     method="BFGS",control = list(maxit = 10000))
  }else if (est.method=="PORT"){
    gmm.stp2=nlminb(cfs1,gmm.obj,x=ivdat[index(na.omit(ivdat))], 
                    gradient=gmm.obj.grad,
                    mom.fun=g.gmm,wmat=inv(Sig.nw),
                    control = list(maxit = 10000))
  }
  cfs2=gmm.stp2$par
  gmm.stp1=gmm.stp2
}
cfs = gmm.stp2$par

J.stat = gmm.obj(cfs,ivdat,g.gmm,inv(Sig.nw))*nrow(ivdat)
J.df = n.mom-n.cfs
J.pval = pchisq(J.stat,J.df,lower.tail = FALSE)


Thta.gmm = eye(col.max)
if (N>col.max){
  Thta.gmm = cbind(Thta.gmm,zeros(col.max,N-col.max))
}
idx =0
for (i in seq(1,(N*col.max-col.max))){
  if (Thta.gmm[i+idx]==1){
    idx = idx+1
  }
  Thta.gmm[i+idx]=cfs[i]
}
Thta.gmm = t(Thta.gmm)
if (N==col.max){colnames(Thta.gmm)=names(var.list)}else{colnames(Thta.gmm)=shks.id}
rownames(Thta.gmm) = names(var.list)

mom.at.sol = t(g.gmm(cfs,ivdat))
if (use.nwhac){
  S.T = 0.5*mom.at.sol%*%t(mom.at.sol)
  if (nwhac.lag>0){
    for (j in seq(nwhac.lag)){
      S.T = S.T + (1-j/(nwhac.lag+1))*mom.at.sol[,(j+1):ncol(mom.at.sol)]%*%
        t(mom.at.sol[,1:(ncol(mom.at.sol)-j)])
    }
  }
  S.T = (1/(ncol(mom.at.sol)))*(S.T+t(S.T))
}else{
  S.T = mom.at.sol%*%t(mom.at.sol)/ncol(mom.at.sol)
}

#npar = size(cfs,2)
#numgrad = NULL
#for (i in seq(npar)){
#  numgrad = cbind(numgrad,
#                  colMeans((g.gmm((cfs+iota(i,npar)*h),ivdat)-
#                              g.gmm((cfs-iota(i,npar)*h),ivdat))/(2*h))
#  )
#}
sol.grad = G.gmm(cfs,ivdat)

bread = inv(t(sol.grad)%*%inv(Sig.nw)%*%sol.grad)
meat = t(sol.grad)%*%inv(Sig.nw)%*%S.T%*%inv(Sig.nw)%*%sol.grad
Sig_Thta = bread%*%meat%*%bread/ncol(mom.at.sol)


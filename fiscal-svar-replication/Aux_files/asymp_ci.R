#calculate asymptotic confidence intervals for impulse response functions

#first, get covariance matrix for VAR parameters
nobs = BP.var$obs #number of observations
if (p>1){
  X = na.omit(cbind(endog.var[(p+1):(nobs+p)],lag(endog.var,1:(p-1))))
}else{
  X = endog.var[(p+1):(nobs+p)]
}
XL = na.omit(lag(endog.var,1:p))
W = cbind(ones(nobs,1),exog.var[(p+1):(nobs+p)])
Mw = eye(nobs)-W%*%inv(t(W)%*%W)%*%t(W)
MwX = Mw%*%X
MwXL = Mw%*%XL

B.x = t(inv(t(MwXL)%*%MwXL)%*%t(MwXL)%*%MwX) #replicates transition matrix "B"
if (p>1){
  U.x = rbind((t(MwX) - B.x%*%t(MwXL))[1:N,],zeros(N*(p-1),nobs))
}else{
  U.x = (t(MwX) - B.x%*%t(MwXL))[1:N,]
}

Sig.u = U.x%*%t(U.x)/nobs

Sig.b = kron(inv(t(MwXL)%*%MwXL),Sig.u) # covariance matrix for var parameters
 # note that kron(Sig.u,inv(t(MwXL)%*%MwXL)) matches vcov() after adjusting for degrees of freedom
 # vcov(BP.var)*(nobs-N*p-size(W,2))/nobs
 # switch order of kron() because vcov() goes eqn-by-eqn. we vectorize B which goes by regressor
#now, get covariance matrix for structural parameters
if (use.nwhac){
  XU = matrix(repmat(MwXL,N*p,1),nobs,N*p*N*p)*repmat(t(U.x),1,N*p) #estfun(BP.var)
  brd.mat = nobs*kron(inv(t(MwXL)%*%MwXL),eye(N*p)) #multiply by nobs to match Bread() 
  meat.mat = 0.5*t(XU)%*%XU
  if (nwhac.lag>0){
    for (j in seq(1,nwhac.lag)){
      meat.mat = meat.mat + (1-j/(nwhac.lag+1))*t(XU[1:(nobs-j),])%*%XU[(j+1):nobs,]
    }
  }
  meat.mat = (1/(nobs))*(meat.mat+t(meat.mat))
  Sig.b = brd.mat%*%meat.mat%*%brd.mat/nobs
}

#covariance matrix for GMM parameters
if (use.nwhac){
    #can use NeweyWest command for gmm model
    #Sig_Thta = NeweyWest(mod.gmm,lag=nwhac.lag,prewhite = FALSE,adjust=FALSE)
    mom.at.sol = t(g.gmm(cfs,ivdat))
    S.T = 0.5*mom.at.sol%*%t(mom.at.sol)
    if (nwhac.lag>0){
      for (j in seq(nwhac.lag)){
        S.T = S.T + (1-j/(nwhac.lag+1))*mom.at.sol[,1:(ncol(mom.at.sol)-j)]%*%
          t(mom.at.sol[,(j+1):ncol(mom.at.sol)])
      }
    }
    S.T = (1/(ncol(mom.at.sol)))*(S.T+t(S.T))
  }else{
    Sig_Thta = vcov(mod.gmm) 
    S.T = mom.at.sol%*%t(mom.at.sol)/ncol(mom.at.sol)
}
sol.grad = G.gmm(cfs,ivdat)
bread = inv(t(sol.grad)%*%inv(Sig.nw)%*%sol.grad)
meat = t(sol.grad)%*%inv(Sig.nw)%*%S.T%*%inv(Sig.nw)%*%sol.grad
Sig_Thta = bread%*%meat%*%bread/ncol(mom.at.sol)

#choose only Theta parameters (ie. not variances of shocks from over-identifying restrictions)
Sig_Thta = Sig_Thta[1:(N*col.max-col.max),1:(N*col.max-col.max)]
#first letter of label is response variable, second letter causal variable
sig.thta.names = NULL
for (i in names(var.list)){ #rows of Theta
  if (use.cov.res){
    shk.names = names(var.list)
  }else{
    shk.names = shks.id
  }
  for (j in shk.names){ #columns of Theta
    if (i==j){
      next
    } else{
      sig.thta.names = c(sig.thta.names,paste0(i,j))
    }
  }
}
colnames(Sig_Thta) = sig.thta.names
rownames(Sig_Thta) = sig.thta.names

Sig_Thta_Full = diag(0,N*col.max,N*col.max)#covar matrix of vec(Thta.gmm)
#names: response variable, shock variable
colnames(Sig_Thta_Full) = paste0(rep(names(var.list),L),rep(shk.names,each=N))
rownames(Sig_Thta_Full) = paste0(rep(names(var.list),L),rep(shk.names,each=N))

for (i in rownames(Sig_Thta)){
  for (j in colnames(Sig_Thta)){
      Sig_Thta_Full[i,j] = Sig_Thta[i,j]
  }
}

Psi_var = array(NA,c(N*col.max,N*col.max,irflen)) #storage matrix
for (h in seq(0,(irflen-1))){
  Cbar = kron(eye(col.max),J%*%(B%^%h)%*%t(J))
  C = zeros(N*col.max,(N*p)^2)
  if (h >0){
    for (i in seq(0,h-1)){
      C = C + kron(t(Thta.gmm)%*%J%*%(t(B)%^%(h-1-i)),J%*%(B%^%i))
    }
  }
  Psi_var[,,(h+1)] = C%*%Sig.b%*%t(C) + Cbar%*%Sig_Thta_Full%*%t(Cbar)
}
#get standard errors
se.irf = sqrt(t(apply(Psi_var,3,diag)))

#rename columns (response variable, shock variable)
colnames(se.irf) = colnames(Sig_Thta_Full)

#now get standard errors for multipliers

#first, calculate variances of i.bar and x.bar

var.x.bar = NA*eye(3)
rownames(var.x.bar) = c("GY","TY","i3")
colnames(var.x.bar) = c("GY","TY","i3")
var.x.bar["i3","i3"] = summary(lm(df[rng,c("TBILL3")]/100~1))$coefficients[2]^2 
var.x.bar["GY","GY"] = summary(lm(df[rng,var.list["G"]]/df[rng,var.list["Y"]]~1))$coefficients[2]^2 
var.x.bar["TY","TY"] = summary(lm(df[rng,var.list["T"]]/df[rng,var.list["Y"]]~1))$coefficients[2]^2 


#do this one variable at a time
se.mult = array(NA,c(irflen,2),dimnames = list(seq(0,(irflen-1)),c("T","G")))
se.dmult = array(NA,c(irflen,2),dimnames = list(seq(0,(irflen-1)),c("T","G")))

for (shk.var in c("T","G")){
  res.var = "Y" #only need to calculate multipliers for output
  ypos = which(names(var.list)==res.var)
  tpos = which(names(var.list)==shk.var)
  Sig.Thta.mult = Sig_Thta_Full[paste0(names(var.list),rep(shk.var,N)),
                                paste0(names(var.list),rep(shk.var,N))]
  Thta.shk = matrix(Thta.gmm[,shk.var],N,1) #column of Theta associated with shock
  eval(parse(text=paste0("shkY.ratio = ",shk.var,"Y.bar"))) #spending or tax revenue to output ratio
  #first, construct derivative matrices
  
  #storage for C matrices
  Cyt = array(0,c(irflen,(N*p)^2))
  Ctt = array(0,c(irflen,(N*p)^2))
  Cby = array(0,c(irflen,N))
  Cbt = array(0,c(irflen,N))
  
  #iterate through IRF horizons
  for (h in seq(0,(irflen-1))){
    if (h>0){
      for (i in seq(0,h-1)){
        Cyt[h+1,] = Cyt[h+1,] + kron(t(Thta.shk)%*%J%*%(t(B)%^%(h-1-i)),
                                     iota(ypos,N)%*%J%*%(B%^%i))
        Ctt[h+1,] = Ctt[h+1,] + kron(t(Thta.shk)%*%J%*%(t(B)%^%(h-1-i)),
                                     iota(tpos,N)%*%J%*%(B%^%i))
      }
    }
    Cby[h+1,] = iota(ypos,N)%*%J%*%(B%^%h)%*%t(J)
    Cbt[h+1,] = iota(tpos,N)%*%J%*%(B%^%h)%*%t(J)
    
    #now get derivatives
    dmdthta = matrix(
      ((1+i.bar)^(-seq(0,h))%*%Cby[1:(h+1),])*(1/shkY.ratio)/
        (sum(((1+i.bar)^(-seq(0,h)))*IRF[tpos,tpos,1:(h+1)])) -
        (sum((1+i.bar)^(-seq(0,h))*IRF[ypos,tpos,1:(h+1)]))*
        (sum(((1+i.bar)^(-seq(0,h)))*IRF[tpos,tpos,1:(h+1)])^(-2))*(1/shkY.ratio)*
        ((1+i.bar)^(-seq(0,h))%*%Cbt[1:(h+1),]),
      1,N)
    dmdb = matrix(
      ((1+i.bar)^(-seq(0,h))%*%Cyt[1:(h+1),])*(1/shkY.ratio)/
        (sum(((1+i.bar)^(-seq(0,h)))*IRF[tpos,tpos,1:(h+1)])) -
        (sum((1+i.bar)^(-seq(0,h))*IRF[ypos,tpos,1:(h+1)]))*(1/shkY.ratio)*
        (sum(((1+i.bar)^(-seq(0,h)))*IRF[tpos,tpos,1:(h+1)])^(-2))*
        ((1+i.bar)^(-seq(0,h))%*%Ctt[1:(h+1),]),
      1,(N*p)^2)
    
    dmdxybar = -(1/shkY.ratio)*MULT[h+1,shk.var]
    
    if (h>0){
      dmdxi = 
        (sum(-seq(0,h)*(1+i.bar)^(-seq(0,h)-1)*IRF[ypos,tpos,1:(h+1)]))*(1/shkY.ratio)*
        (sum(((1+i.bar)^(-seq(0,h)))*IRF[tpos,tpos,1:(h+1)])^(-1))-
        (sum((1+i.bar)^(-seq(0,h))*IRF[ypos,tpos,1:(h+1)]))*(1/shkY.ratio)*
        (sum(((1+i.bar)^(-seq(0,h)))*IRF[tpos,tpos,1:(h+1)])^(-2))*
        (sum(-seq(0,h)*(1+i.bar)^(-seq(0,h)-1)*IRF[tpos,tpos,1:(h+1)]))
    }else{
      dmdxi = 0
    }
    
    se.mult[h+1,shk.var] = sqrt(
      dmdb%*%Sig.b%*%t(dmdb) + 
        dmdthta%*%Sig.Thta.mult%*%t(dmdthta) + 
        dmdxybar^2*var.x.bar[paste0(shk.var,"Y"),paste0(shk.var,"Y")] +
        dmdxi^2*var.x.bar["i3","i3"]
    )
    #repeat for (simpler) dynamiv multipliers
    dmdb = Cyt[(h+1),]*(1/shkY.ratio)
    dmdthta = Cby[(h+1),]*(1/shkY.ratio)
    dmdxybar = -(1/shkY.ratio)*DMULT[h+1,shk.var]
    se.dmult[(h+1),shk.var] = sqrt(
      dmdb%*%Sig.b%*%(dmdb) + 
      dmdthta%*%Sig.Thta.mult%*%(dmdthta) + 
      dmdxybar^2*var.x.bar[paste0(shk.var,"Y"),paste0(shk.var,"Y")] 
    )
  }
}



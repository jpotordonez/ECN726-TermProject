#calculate impulse response functions and multipliers
# ie. response of N variables to L identified shocks

IRF = array(0,c(size(Thta.gmm),irflen))
if (p>1){
  J = cbind(eye(N), zeros(N,N*(p-1)))
}else{
  J = eye(N)
}

for (j in 1:irflen){
  IRF[,,j] = J%*%(B%^%(j-1))%*%t(J)%*%Thta.gmm
}
if (use.cov.res){
  shk.names = names(var.list)
}else{
  shk.names = shks.id
}
colnames(IRF) = shk.names
rownames(IRF) = names(var.list)

#now, calculate cumulative multipliers. Follow Mountford and Uhlig (2009)
MULT = array(NA,c(irflen,2),dimnames=list(seq(1,irflen),c("T","G")))
i.bar = mean(df[rng,c("TBILL3")])/100
TY.bar = mean(df[rng,var.list["T"]]/df[rng,var.list["Y"]])
GY.bar = mean(df[rng,var.list["G"]]/df[rng,var.list["Y"]])

for (j in 0:(irflen-1)){
  MULT[j+1,"T"]=(1/TY.bar)*sum((1+i.bar)^(-seq(0,j))*IRF["Y","T",1:(j+1)])/
    sum(((1+i.bar)^(-seq(0,j)))*IRF["T","T",1:(j+1)])
  MULT[j+1,"G"]=(1/GY.bar)*sum((1+i.bar)^(-seq(0,j))*IRF["Y","G",1:(j+1)])/
    sum(((1+i.bar)^(-seq(0,j)))*IRF["G","G",1:(j+1)])
}

#now get dynamic multipliers, which are just the IRFs scaled by TY.bar or GY.bar

DMULT = array(NA,c(irflen,2),dimnames=list(seq(1,irflen),c("T","G")))
DMULT[,"T"] = IRF["Y","T",]/TY.bar
DMULT[,"G"] = IRF["Y","G",]/GY.bar

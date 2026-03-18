#make instrument data

shks.list = list() #list of all instruments for each identified shocks
for (shk in shks.id){
  #add instruments to list
  shks.list[[shk]] = c(eval(parse(text=paste0(shk,".inst"))))
}

Z = na.omit(df[rng,unlist(shks.list)][p+1:nrow(U),]) #instruments

if (pre.white){#optionally regress of p lags of the endogenous VAR variables
  Z.tmp = na.omit(merge(Z,lag(df[,var.list],seq(1,p)),exog.var))
  n.insts = length(unlist(shks.list))
  for (i in unlist(shks.list)){
    Z[,i]=residuals(lm(Z.tmp[,i]~Z.tmp[,(n.insts+1):size(Z.tmp,2)]))
  }
}

#data for IV regressions. shrink sample to remove NAs and match GMM sample.
ivdat=na.omit(merge(xts(U,as.Date(rownames(U))),Z))
#rename instruments and residuals
names(ivdat) = c(names(var.list),unlist(shks.list))

#optionally de-mean all series
if (demean.inst){
  for (i in seq(ncol(ivdat))){
    ivdat[,i] = ivdat[,i] - mean(ivdat[,i])
  }
}

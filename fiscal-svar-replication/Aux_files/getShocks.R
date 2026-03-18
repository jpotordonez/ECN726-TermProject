#get columns of theta associated with each identified shock
for (shk in shk.names){
  eval(parse(text=paste0("Thta.",shk,"=matrix(Thta.gmm[,'",shk,"'],N,1)")))
}

if (N==col.max){
  eps = t(inv(Thta.gmm)%*%t(U))
  colnames(eps) = names(var.list) #rename columns 
} 

#compare with Stock and Watson method to uncover shocks
eps.sw = NULL
for (shk in shk.names){
  eval(parse(text=paste0("sig.sw.",shk,shk,"=inv(t(Thta.",shk,")%*%inv(var.cov)",
                         "%*%Thta.",shk,")")))
  eval(parse(text=paste0("eps.sw.",shk,shk,"=sig.sw.",shk,shk,"%*%t(Thta.",shk,
                         ")%*%inv(var.cov)%*%t(U)")))
  eval(parse(text=paste0("eps.sw = cbind(eps.sw,t(eps.sw.",shk,shk,"))")))
}
colnames(eps.sw) = shk.names

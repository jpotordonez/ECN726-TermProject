#starting values for GMM estimation
svals = NULL
for (i in seq(1,N)){#loop through all endogenous variables
  if (use.cov.res){
    col.max = N
  }else{
    col.max = L
  }
  for (j in seq(1,col.max)){#loop through all identified variables
    if (i==j){
      next
    } else{
      #use 0 as starting value if no instrument
      if (is.null(shks.list[[names(var.list[j])]])){
        svals = c(svals,0)
      }else{
        svals = c(svals,ivreg(ivdat[,names(var.list[i])]~
                                ivdat[,names(var.list[j])]-1|
                                ivdat[,shks.list[[names(var.list[j])]]])$coef[1])
      }
    }
  }
}

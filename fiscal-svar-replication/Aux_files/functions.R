g.gmm = function(Thta,x){
  #gmm function
  #relabel data
  for (shk in shks.id){
    eval(parse(text=paste0("z.",shk,"=x[,shks.list[['",shk,"']]]")))
  }
  # IV moment conditions
  f = NULL
  thta.idx = 0
  for (i in seq(N)){#loop through rows
    for (j in seq(col.max)){#loop through (identified) columns
      if (i ==j){#skip diagonal elements
        next
      }else{
        thta.idx = thta.idx + 1
        if (is.na(shks.id[j])){#skip if this shock not identified by an instrument
          next
        }
        ninst = size(eval(parse(text=paste0("z.",shks.id[j]))),2) #number of instruments for shock
        for (k in seq(1,ninst)){
          eval(parse(text=paste0("mom.con=",
                                 "(x[,names(var.list[i])]-Thta[thta.idx]*x[,j])*z.",shks.id[j],"[,k]")))
          f = cbind(f,mom.con)
        }
      }
    }
  }  
  if (use.cov.res){
    Thta.mat = eye(N)
    idx =0
    for (i in seq(1,(N*(N-1)))){
      if (Thta.mat[i+idx]==1){
        idx = idx+1
      }
      Thta.mat[i+idx]=Thta[i]
    }
    Thta.mat = t(Thta.mat)
    eps = t(inv(Thta.mat)%*%t(x[,names(var.list)]))
    colnames(eps) = names(var.list)
    #moment conditions
    for (i in seq(1,length(cov.include))){
      f = cbind(f,eps[,substr(cov.include[i],1,1)]*eps[,substr(cov.include[i],2,2)])
    } 
  }
  names(f) = NULL
  return(f)
}

G.gmm = function(Thta,x){
  #analytical derivatives for gmm estimation
  #if (N>3){
  # warning("error: only working for N=3")
  #}
  #relabel data
  for (shk in shks.id){
    eval(parse(text=paste0("z.",shk,"=x[,shks.list[['",shk,"']]]")))
  }
  npar = size(Thta,2) #number of parameters
  # IV moment conditions
  f = NULL
  thta.idx = 0
  for (i in seq(1,N)){#loop through rows of Theta
    for (j in seq(col.max)){
      if (i ==j){ #skip diagonal elements
        next
      }else if (is.na(shks.id[j])){#skip if parameter not identified by instrument
        thta.idx = thta.idx+1
        next
      }else{
        thta.idx = thta.idx + 1
        ninst = size(eval(parse(text=paste0("z.",shks.id[j]))),2) #num instruments for shock
        for (k in seq(1,ninst)){
          eval(parse(text=paste0("f=rbind(f,",
                                 "cbind(matrix(0,1,(thta.idx-1)),","
            -mean(x[,j]*z.",shks.id[j],"[,k]),matrix(0,1,(npar-thta.idx))))")))
        }
      }
    }
  }  
  #moments for VAR covariance matrix Sig.eps
  if (use.cov.res){
    #build theta matrix
    Thta.mat = eye(N)
    idx =0
    for (i in seq(1,(N*N-N))){
      if (Thta.mat[i+idx]==1){
        idx = idx+1
      }
      Thta.mat[i+idx]=Thta[i]
    }
    Thta.mat = t(Thta.mat)
    Thta.inv = inv(Thta.mat)
    Sig.u = t(x[,names(var.list)])%*%x[,names(var.list)]/nrow(x)
    
    # derivatives of moment conditions
    # moment based on E[eps_ieps_j'] = 0
    # i and j index moments (rows)
    for (i in seq(1,N)){ 
      for (j in seq(1,N)){ 
        if (j<=i){next}
        #make sure restriction is included
        if (!any(paste0(names(var.list)[i],names(var.list)[j])==cov.include)){next}
        thta.grad=NULL
        #k and l index parameters (columns)
        for (k in seq(1,N)){
          for (l in seq(1,N)){
            if (k==l){next}
            thta.grad=c(thta.grad,sum(diag(-t((t(iota(i,N))%*%iota(j,N)+
                                               t(iota(j,N))%*%iota(i,N))%*%
                                                Thta.inv%*%Sig.u)%*%
                                             t(t(Thta.inv[,k]))%*%Thta.inv[l,])))
          }
        }
        f = rbind(f,thta.grad)
      }
    }
  }
  return(f)
}


iota = function(j,k){
  if (j>k){
    return(warning("'j' must be less than or equal to 'k"))
  }
  xx = zeros(1,k)
  xx[j] = 1
  return(xx)
}

gmm.obj = function(theta,x,mom.fun,wmat){
  #GMM objective function, where theta is a vector of parameters, x is the data,
  #mom.fun is a function of moment conditions, and wmat is a weight matrix.
  mom.avg = colMeans(mom.fun(theta,x))
  return(t(mom.avg)%*%wmat%*%mom.avg)
}

gmm.obj.grad = function(theta,x,mom.fun,wmat){
  return(2*colMeans(mom.fun(theta,x))%*%wmat%*%G.gmm(theta,x))
}




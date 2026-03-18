#get Theta matrix
Thta.gmm = eye(L)
if (N>L){
  Thta.gmm = cbind(Thta.gmm,zeros(L,N-L))
}
idx =0
for (i in seq(1,(N*L-L))){
  if (Thta.gmm[i+idx]==1){
    idx = idx+1
  }
  Thta.gmm[i+idx]=cfs[i]
}
Thta.gmm = t(Thta.gmm)
colnames(Thta.gmm) = shks.id
rownames(Thta.gmm) = names(var.list)
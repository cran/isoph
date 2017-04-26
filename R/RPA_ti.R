RPA_ti_inc=function(n, nt, m, z, z.obs, Y, dN){
  Y2=matrix(0,m,nt)
  dN2=matrix(0,m,nt)

  z.obs=c(z.obs,Inf)
  for(h in 1:m){
    idx=which(z.obs[h]<=z & z<z.obs[h+1])
    if(length(idx)==1){
      Y2[h,]=Y[idx,]
      dN2[h,]=dN[idx,]
    }else{
      Y2[h,]=colSums(Y[idx,])
      dN2[h,]=colSums(dN[idx,])
    }
  }
  return(list(Y2=Y2,dN2=dN2));
}

RPA_ti_dec=function(n, nt, m, z, z.obs, Y, dN){
  Y2=matrix(0,m,nt)
  dN2=matrix(0,m,nt)

  z.obs=c(-Inf,z.obs)
  for(h in 1:m){
    idx=which(z.obs[h]<z & z<=z.obs[h+1])
    if(length(idx)==1){
      Y2[h,]=Y[idx,]
      dN2[h,]=dN[idx,]
    }else{
      Y2[h,]=colSums(Y[idx,])
      dN2[h,]=colSums(dN[idx,])
    }
  }
  return(list(Y2=Y2,dN2=dN2));
}

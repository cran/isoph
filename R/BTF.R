BTF_inc=function(m, n, z,z.obs, psi.new){
  psi.full=rep(-Inf,n)
  for(h in 1:m)
    psi.full[z>=z.obs[h]]=psi.new[h]
  return(psi.full);
}

BTF_dec=function(m, n, z,z.obs, psi.new){
  psi.full=rep(-Inf,n);
  for(h in m:1)
    psi.full[z<=z.obs[h]]=psi.new[h]
  return(psi.full);
}


isoph.ti=function(T, STATUS, Z, shape, K, maxiter, eps, maxdec){

  #sorted by z
  n=length(STATUS)
  order.z=order(Z)
  t=T[order.z]
  status=STATUS[order.z]
  z=sort(Z)
  z.obs=unique(z[status==1])
  m=length(z.obs)
  
  t.obs=sort(unique(t))
  nt=length(t.obs)
  
  #anchor
  k=sum(z.obs<K)
  if(k==0) k=1 #for the case when min(z.obs) < K
  zk=z.obs[k]
  
  #counting process
  Y=dN=matrix(0,n,nt)       #row is the subj, col is the time corresponding z_(i);
  for(i in 1:n){
    rank.t=which(t[i]==t.obs)
    Y[i,][1:rank.t]=1
    if(status[i]==1) dN[i,][rank.t]=1
  }  
  
  #RPA (with interval), Y, dN
  Y2=matrix(0,m,nt)
  dN2=matrix(0,m,nt)
  
  rpa.Y=.C('RPA_ti', as.integer(n), as.integer(nt), as.integer(m), as.double(z), as.double(z.obs) , as.double(Y), as.integer(dN), Y2=as.double(Y2), dN2=as.integer(dN2) )
  dN2=matrix(rpa.Y$dN2,m,nt)
  Y2 =matrix(rpa.Y$Y2, m,nt)
  
  #initial value
  zk=z.obs[k]
  z.bar=z-zk
  
  beta.hat=coxph(Surv(t,status)~z.bar)$coefficient      
  if(shape=='increasing'){
    psi= abs(beta.hat)*(z.obs-zk)
  }else if(shape=='decreasing'){
    psi=-abs(beta.hat)*(z.obs-zk)
  }
  
  #picm
  conv="not converged"
  iter=0
  d.e=1
  while(d.e>=eps){  
    iter=iter+1
    
    #picm
    dNsum=colSums(dN2)
    Delta=rowSums(dN2)     
    
    den=colSums(Y2*exp(psi))
    index.zero=which(den>0) #0/0=0
    
    mu=c()
    for(s in 1:m)
      mu[s]=sum( (Y2[s,]/den*dNsum)[index.zero] )
    
    if(shape=='increasing'){
      exp.psi.new=pava(Delta/mu, mu)
    }else if(shape=='decreasing'){
      exp.psi.new=-pava(-Delta/mu, mu)
    }
    psi.new=log(exp.psi.new)    
    
    #distance
    d.e=sum(abs(exp(psi.new)-exp(psi)))
    psi=psi.new
    if(iter>maxiter) break
  }
  if(d.e<eps) conv="converged"
  
  #impose the anchor
  psi.new=psi.new-psi.new[k] #psi.obs is the same as psi.new;
  
  #back to full rank
  psi.obs=round(psi.new, maxdec)
  
  #level sets
  psi.uniq=unique(psi.obs)
  n.lv=length(psi.uniq)
  
  lv.sets=c()
  zmin=formatC( round(min(z),maxdec), format='f', digits=maxdec)
  zmax=formatC( round(max(z),maxdec), format='f', digits=maxdec)
  
  if(n.lv==1){ #only one level sets
    lv.sets[1]=paste('[',zmin,',',zmax,']', sep='')
  }else{
    lv=c()
    for(i in 1:n.lv)
      lv[[i]]=formatC( round(z.obs[which(psi.obs==psi.uniq[i])],maxdec)[1], format='f', digits=maxdec)
    for(i in 1:n.lv){
      if(i==1){
        lv.sets[1]=paste('[',zmin,', ',lv[2],')', sep='')
      }else if(i==n.lv){
        lv.sets[i]=paste('[',lv[i][1],', ',zmax,']', sep='')
      }else{
        lv.sets[i]=paste('[',lv[i],', ',lv[i+1],')', sep='')
      }
    }
  }
  psi.hat=formatC( unique(psi.obs), format='f', digits=maxdec)
  HR.hat=formatC( unique(exp(psi.obs)), format='f', digits=maxdec)
  
  est=data.frame(psi.hat=psi.hat, HR.hat=HR.hat, lv.set=lv.sets)
  names(est)=c("psi.hat","exp(psi.hat)","level set of psi.hat")
  
  #for plot
  psi.obs=c(psi.obs[1],psi.obs,psi.obs[m])
  z.obs=c(z.obs[1],z.obs,z.obs[m])
  
  return(list(est=est, conv=conv, psi=psi.obs, z=z.obs, K=K, shape=shape, iter=iter, dist=d.e, n=n, nevent=sum(STATUS), njump=m))
}

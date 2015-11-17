isoph.td=function(START, STOP, STATUS, Z, shape, K, maxiter, eps, maxdec){

  #data set with observed failure time
  oft.index=(which(START==0)-1)      #observed failure time (oft) for each subj
  oft.index=c(oft.index[-1], length(START)) #-1 for the first subj, with combining oft of the last subj
  
  T.oft=STOP[oft.index]
  STATUS.oft=STATUS[oft.index]
  Z.oft=Z[oft.index]
  n=length(STATUS.oft)
  if( sum(STATUS.oft)!=sum(STATUS) )
    stop("It is not allowed to have an event (status=1) other than the last interval of each subject.")

  #sorted by z
  order.z=order(Z.oft)
  t=T.oft[order.z]
  
    #create ZT
    t.obs=sort(unique(t))
    nt=length(t.obs)    

    ZT=matrix(0,n,nt) #ZT=ZT=matrix(NA,n,nt) does not work for the c code

    for(i in 1:n)
      for(j in oft.index[i]:(max(0,oft.index[i-1])+1))
        ZT[i,which(t.obs<=STOP[j])]=Z[j]
  
  zt=ZT[order.z,]
  status=STATUS.oft[order.z]
  z=sort(Z.oft)
  z.obs=unique(z[status==1])
  m=length(z.obs)
  
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
  
  rpa.Y=.C('RPA_td', as.integer(n), as.integer(nt), as.integer(m), as.double(zt), as.double(z.obs) , as.double(Y), as.integer(dN), Y2=as.double(Y2), dN2=as.integer(dN2) )
  dN2=matrix(rpa.Y$dN2,m,nt)
  Y2 =matrix(rpa.Y$Y2, m,nt)
    
  #initial value
  beta.hat=coxph(Surv(START,STOP,STATUS)~Z)$coefficient  
  if(shape=='increasing'){
     psi= abs(beta.hat[1])*(z.obs-zk)
  }else if(shape=='decreasing'){
     psi=-abs(beta.hat[1])*(z.obs-zk)
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
    
    weight=c()
    for(s in 1:m)
      weight[s]=sum( (Y2[s,]*dNsum/den)[index.zero] )
    
    if(shape=='increasing'){
      exp.psi.new=pava(Delta/weight, weight)
    }else if(shape=='decreasing'){
      exp.psi.new=-pava(-Delta/weight, weight)
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
  z.obs=c(min(Z),z.obs,max(Z))
  
  return(list(est=est, conv=conv, psi=psi.obs, z=z.obs, K=K, shape=shape, iter=iter, dist=d.e, n=n, nevent=sum(STATUS), njump=m))
}
isoph.ti=function(TIME, STATUS, Z, X, shape, K, maxdec, maxiter, eps){
  #the comp speed can be further improved by separating shape, inc or dec

  #sorted by z
  n=length(STATUS)
  data=data.frame(TIME,STATUS,Z)
  if(!is.null(X)){
    data=data.frame(TIME,STATUS,Z,X)
    nparm=1; if(!is.vector(X)) nparm=ncol(X)
    if(nparm>1){; NR.ft=NR2.ft
    }else{;       NR.ft=NR1.ft
    }
  }
  data=data[order(data$Z),]

  if(shape=='increasing'){
    RPA_ti=RPA_ti_inc;
    BTF=BTF_inc;
  }else{
    RPA_ti=RPA_ti_dec;
    BTF=BTF_dec;
  }

  #remove subjects whose cov is less than z^*_(1)
  n.total=n
  if(shape=='increasing'){
    if(sum(data$STATUS[data$Z==data$Z[1]])==0){
      idx=which(status==1)[1]
      idx2=idx:n;      n=length(idx2)
    }
  }else{
    if(sum(data$STATUS[data$Z==data$Z[n]])==0){
      idx=tail(which(data$STATUS==1),1)
      idx2=1:idx;      n=length(idx2)
    }
  }
  if(n.total>n)  data=data[idx2,]

  #re-define variable
  t=data$TIME
  status=data$STATUS
  z=data$Z
  if(!is.null(X)){
    if(nparm==1){ x=as.numeric(data[,4])
    }else{;       x=as.matrix(data[,-c(1,2,3)])
    }
  }

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

  #initial
  zk=z.obs[k]
  Z.BAR=Z-zk

  if(is.null(X)){ #no trt group
    try(beta.hat<-coxph(Surv(TIME,STATUS)~Z.BAR)$coefficient, silent=T)
    if(!is.numeric(beta.hat)) beta.hat=0.01

    if(shape=='increasing'){;       psi= abs(beta.hat)*(z.obs-zk)
    }else if(shape=='decreasing'){; psi=-abs(beta.hat)*(z.obs-zk);  }
  }else{
    try(beta.hat<-coxph(Surv(TIME,STATUS)~Z.BAR+X)$coefficient, silent=T)
    if(!is.numeric(beta.hat)){; beta.hat=0.01; beta=rep(0,nparm);
    }else{ beta=beta.hat[-1];    }
    if(shape=='increasing'){;       psi= abs(beta.hat[1])*(z.obs-zk)
    }else if(shape=='decreasing'){; psi=-abs(beta.hat[1])*(z.obs-zk);  }
  }

  #picm & beta for newton raphson algo
  if(is.null(X)){ #no trt group
    #RPA (with interval), Y, dN
    rpa.Y=RPA_ti(n, nt, m, z, z.obs, Y, dN)
    Y2=rpa.Y$Y2
    dN2=rpa.Y$dN2

    #picm
    dNsum=colSums(dN2)
    Delta=rowSums(dN2)

    dist=0; exp.beta=NA
    picm=picm.ft(psi,m,z.obs,zk,k, dN2,Y2,dNsum,Delta, eps,maxiter, shape)
    psi.new=picm$psi.new

    psi.full=BTF(m, n, z,z.obs, psi.new)
    if(picm$conv==0) stop
  }else{ #trt group
    iter=0;  dist=1; beta.new=rep(NA,nparm)
    while(dist>=eps){
      iter=iter+1
      if(iter>maxiter) break

      #RPA (with interval), Y, dN
      Yest=matrix(NA,n,nt)
      if(nparm>1){ for(j in 1:nt) Yest[,j]=Y[,j]*exp(x%*%beta)
      }else{       for(j in 1:nt) Yest[,j]=Y[,j]*exp(x*beta)
      }

      rpa.Y=RPA_ti(n, nt, m, z, z.obs, Yest, dN)
      Y2=rpa.Y$Y2
      dN2=rpa.Y$dN2

      #picm
      dNsum=colSums(dN2)
      Delta=rowSums(dN2)

      #estimate psi
      picm=picm.ft(psi,m,z.obs,zk,k, dN2,Y2,dNsum,Delta, eps,maxiter, shape)
      psi.new=picm$psi.new

      psi.full=BTF(m, n, z,z.obs, psi.new)
      if(picm$conv==0) stop

      #estimate beta (Y1&w1 or Y2&w2 should be the same);
      beta.new=NR.ft(x,beta,psi.full,n,nt,Y,dN, maxiter,eps)

      #update;
      dist=sqrt(sum((psi.new-psi)^2))+sqrt(sum((beta.new-beta)^2))
      #this can be reduced to: (1-picm$conv)+sqrt(sum((beta.new-beta)^2))

      psi=psi.new
      beta=beta.new
    }
    exp.beta=round(exp(beta.new), maxdec)
    exp.beta=formatC( exp.beta, format='f', digits=maxdec)
  }

  ###
  #picm result
  ###
  conv="converged"
  if(dist>=eps)    conv="not converged"
  if(picm$conv==0) conv="not converged"

  #back to full rank (later)
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
    lv=list()
    if(shape=='increasing'){
      for(i in 1:n.lv)
        lv[[i]]=formatC( round(z.obs[which(psi.obs==psi.uniq[i])],maxdec)[1], format='f', digits=maxdec)
      for(i in 1:n.lv){
        if(i==n.lv){; lv.sets[i]=paste('[',lv[i][1],', ',zmax,']', sep='')
        }else{;       lv.sets[i]=paste('[',lv[i],', ',lv[i+1],')', sep='')
        }
      }
    }else{
      for(i in 1:n.lv)
        lv[i]=formatC( tail(round(z.obs[which(psi.obs==psi.uniq[i])],maxdec),1), format='f', digits=maxdec)
      for(i in 1:n.lv){
        if(i==1){;   lv.sets[i]=paste('[',zmin,', ',lv[i][1],']', sep='')
        }else{;      lv.sets[i]=paste('(',lv[i-1],', ',lv[i],']', sep='')
  };};};};

  psi.hat=formatC( unique(psi.obs), format='f', digits=maxdec)
  HR.hat=formatC( unique(exp(psi.obs)), format='f', digits=maxdec)
  #recover if status[1]=0 associated with z(1)
  if(n.total>n){
    n=n.total
    if(shape=='increasing'){
      psi.hat=c(-Inf,psi.hat)
      HR.hat=c(0,HR.hat)
      lv.sets=c(paste('[',min(Z),', ',zmin,')', sep=''), lv.sets)
    }else{
      psi.hat=c(psi.hat,-Inf)
      HR.hat=c(HR.hat,0)
      lv.sets=c(lv.sets, paste('(',zmax,', ',max(Z),']', sep=''))
    }
  }

  est=data.frame(psi.hat=psi.hat, HR.hat=HR.hat, lv.set=lv.sets)
  names(est)=c("psi.hat","exp(psi.hat)","level set of psi.hat")

  z.range=range(Z)

  #temprary to compute psi.full
  psi.full2=rep(NA,n.total)
  z.full=sort(Z)
  for(i in 1:n){
    psi.full2[which(z.full%in%z[i])]=psi.full[i]
  }
  psi.full2[is.na(psi.full2)]=-Inf
  psi.full=psi.full2

  return(list(est=est, exp.beta=exp.beta, conv=conv,
              psi=psi.obs, z=z.obs, z.range=z.range, K=K, shape=shape, n=n, nevent=sum(STATUS), njump=m,
              psi.full=psi.full, z.full=z.full))
}

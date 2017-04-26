#newton-raphson algo with p>=2
NR2.ft=function(x,beta,psi.full,n,nt,Y,dN, maxiter,eps){
  nparm=length(beta)
  iter=0
  dist=1
  while(dist>=eps){   
    iter=iter+1
    if(iter>maxiter) break    
    
    xb=x%*%beta
    Yest=matrix(0,n,1); S1=matrix(0,nparm,1)
    U=matrix(0,nparm,1);  H=matrix(0,nparm,nparm)
    
    for(i in 1:n){
      for(j in which(dN[i,]>=1)){
        Yest=Y[,j]*exp(psi.full+xb) #Y*estimated exp(psi.hat+xb)
        S0=sum(Yest)
        S1=matrix(t(Yest) %*% x)
        E1=S1/S0
        U=U+(x[i,]-E1)*dN[i,j]
        
        S2=matrix(0,nparm,nparm)
        for(m in 1:n)
          S2=S2+Yest[m]*(x[m,]%*%t(x[m,]))
        
        H=H+(-S2/S0+E1%*%t(E1))*dN[i,j]
      }
    }
    
    beta.new = beta - solve(H)%*%U
    
    #distance
    dist=sqrt(sum((beta.new-beta)^2))
    beta=beta.new    
  }
  return(beta.new)
}
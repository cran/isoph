#newton-raphson algo with p=1
NR1.ft=function(x,beta,psi.full,n,nt,Y,dN, maxiter,eps){
  iter=0
  dist=1

  while(dist>=eps){
    iter=iter+1
    if(iter>maxiter) break

    xb=x*beta
    Yest=matrix(0,n,1)
    U=H=0

    for(i in 1:n){
      for(j in which(dN[i,]>=1)){
        Yest=Y[,j]*exp(psi.full+xb) #Y*estimated exp(psi.hat+xb)
        S0=sum(Yest)
        S1=sum(Yest*x)
        S2=sum(Yest*x^2)
        S2=S1^2
        E1=S1/S0

        U=U+(x[i]-E1)*dN[i,j]
        H=H+(-S2/S0+E1^2)*dN[i,j]
      }
    }

    beta.new = beta - U/H

    #distance
    dist=sqrt(sum((beta.new-beta)^2))
    beta=beta.new
  }
  return(beta.new)
}

plot.isoph=function(x, ...){
  pcex=0.7
  lgcex=0.9
  pch=19

  y=y.obs=x$psi
  z=z.obs=x$z
  z.range=x$z.range
  hr=hr.obs=exp(y.obs)

  #-Inf due to censorship
  if(x$shape=='increasing'){
    if(z.range[1]!=z[1]){
      y=c(-Inf,y)
      hr=c(0,hr)
      z=c(z.range[1],z)
    }
  }else if(x$shape=='decreasing'){
    if(z.range[2]!=max(z)){
      y=c(y,-Inf)
      hr=c(hr,0)
      z=c(z,z.range[2])
    }
  }

  xlab=x$formula[[3]]
  xlim=z.range

  lglab="Potential jump points"
  if(x$shape=='increasing'){
    l.type='s'
    lgloc='topleft'
  }
  if(x$shape=='decreasing'){
    l.type='S'
    lgloc='topright'
  }

  #for psi.hat
  #ylab=expression(hat(psi))
  #ylim=range(y[is.finite(y)])
  #plot(y~z, type=l.type, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
  #points(y.obs~z.obs, pch=pch, cex=pcex)

  #for exp(psi.hat)
  ylab=expression(exp(hat(psi)))
  ylim=range(exp(y))

  plot(hr~z, type=l.type, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
  points(hr.obs~z.obs, pch=pch, cex=pcex)

  #legends
  legend(x=lgloc, pch=pch, legend=lglab, bty='n', cex=lgcex)
}

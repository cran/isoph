isoph=function(formula, trt=NULL, data=NULL, shape='increasing', K=0, maxdec=2, maxiter=10^4, eps=10^-3){

  #3. load data with formula
  #3.1 Surv response and cov
  mf=model.frame(formula=formula, data=data) #data must be data.frame

  #3.2 outcome
  surv.y=model.response(mf)
  if(ncol(surv.y)==2){ #for time-indepnedent
    type='ti'
    TIME=surv.y[,1];     STATUS=surv.y[,2]
  }else if(ncol(surv.y)==3){ #for time-depnedent (or time-independent)
    type='td'
    START=surv.y[,1]; STOP=surv.y[,2];  STATUS=surv.y[,3]
  }

  #3.3 cov
  Z=model.matrix(formula, data = mf, rhs = 2)
  Z=Z[,2] #Z[,1] is intercept;

  #3.5 W (or trt)
  if(!is.null(data))    trt=data$trt

  #4. isoph
  if(type=='ti'){
    est=isoph.ti(TIME=TIME, STATUS=STATUS, Z=Z, X=trt, shape=shape, K=K, maxdec=maxdec, maxiter=maxiter, eps=eps)
  }else if(type=='td'){
    #est=isoph.td(START=START, STOP=STOP, STATUS=STATUS, Z=Z, X=NULL, shape=shape, K=K, maxdec=maxdec, maxiter=maxiter, eps=eps)
    stop("time-depdent cov is not supported for the current version of the isoph function")
  }

  est$call=match.call()
  est$formula=formula

  class(est)="isoph"

  return(est)
}

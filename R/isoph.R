isoph=function(formula, data=NULL, shape='increasing', K=0, maxiter=10^5, eps=10^-3, maxdec=2){
  #1. load/check R packages

  #2. check input parameters
  #2.1. Null
  if( is.null(formula) )  stop("formular is requred")
  if( is.null(shape) )    stop("shape is required")
  if( is.null(K) )        stop("K is required")
  if( is.null(eps) )      stop("eps is required")
  if( is.null(maxiter) )  stop("maxiter is required")
  
  #2.2 Length
  if( length(shape)!=1 )    stop("shape must be scalar")
  if( length(K)!=1 )        stop("K must be scalar")
  if( length(eps)!=1 )      stop("eps must be scalar")
  if( length(maxiter)!=1 )  stop("maxiter must be scalar")
  
  #2.3 NA or Inf
  if( (is.na(K)||is.infinite(K)) )              stop("K must be a finite value")
  if( (is.na(eps)||is.infinite(eps)) )          stop("eps must be a finite value")
  if( (is.na(maxiter)||is.infinite(maxiter)) )  stop("maxiter must be a finite value")
  
  #2.4 class
  if(class(formula)!='formula') stop("formula argument is wrong")
  if(class(shape)!='character') stop("shape must be character")
  if(any(grep("inc",tolower(shape)))==TRUE){
    shape="increasing"
  }else if(any(grep("dec",tolower(shape)))==TRUE){
    shape="decreasing"
  }
  if(shape!='increasing' && shape!='decreasing')
    stop("shape must be either increasing or decreasing")

  if( !(is.numeric(K)) )       stop("K must be numeric")
  if( !(is.numeric(eps)) )     stop("eps must be numeric")
  if( !(is.numeric(maxiter)) ) stop("maxiter must be be numeric")

  #2.5 eps & maxiter
  maxiter=ceiling(maxiter)
  if(eps<0)       stop("eps must be greater than or equal to zero")
  if(maxiter<=0)  stop("maxiter must be greater than zero")
  
  #3. load data with formula 
  #3.1 Surv response and cov
  mf=model.frame(formula=formula, data=data) #data must be data.frame
  
  #outcome
  surv.y=model.response(mf)
  if( !(is.Surv(surv.y)) ) stop("Response must be a survival outcome")
  
  if(ncol(surv.y)==2){ #for time-indepnedent
    type='ti'
    TIME=surv.y[,1];     STATUS=surv.y[,2]
  }else if(ncol(surv.y)==3){ #for time-depnedent (or time-independent)
    type='td'
    START=surv.y[,1]; STOP=surv.y[,2];  STATUS=surv.y[,3]
  }
  
  #cov
  Z=model.matrix(formula, data = mf, rhs = 2)
  if (ncol(Z)>=3) stop("The current version of isoph does not support multivariate explanatory variables")
  Z=Z[,2] #Z[,1] is intercept;  
  
  #3.2. check data entry
  if(type=='ti'){
    #NA or Inf
    if ( any(is.na(TIME)+is.na(STATUS)+is.na(Z)) ) stop("Data included NA")
    if ( any(is.infinite(TIME)+is.infinite(STATUS)+is.infinite(Z)) )  stop("Data included infinite values")
    
    #length
    if( !(length(TIME)==length(STATUS) & length(STATUS)==length(Z)) ) stop("Lengths of data in the formula argument are not matched")
    
    #right censored data
    if( min(TIME)<=0 ) stop("Time must be greater than zero")

  }else if(type=='td'){
    #NA or Inf
    if ( any(is.na(START)+is.na(STOP)+is.na(STATUS)+is.na(Z)) ) stop("Data included NA")
    if ( any(is.infinite(START)+is.infinite(STOP)+is.infinite(STATUS)+is.infinite(Z)) ) stop("Data included infinite values") 
    
    #length
    if( !(length(START)==length(STOP) & length(STOP)==length(STATUS) & length(STATUS)==length(Z)) ) stop("Lengths of data are not matched")
    
    #right censored data
    if( min(STOP)<=0 ) stop("Time must be greater than zero")
  }

  #Censoring
  if ( length(unique(STATUS))>=3 ) stop("status has to be either 0 or 1")
  if ( !all(STATUS %in% c(0,1)) )  stop("status has to be either 0 or 1")
  
  if (sum(STATUS)<=2) stop("At least more than two numbers of event are needed.")

  #4. isoph
  if(type=='ti'){
    est=isoph.ti(TIME=TIME, STATUS=STATUS, Z=Z, shape=shape, K=K, maxiter=maxiter, eps=eps, maxdec=maxdec)
  }else if(type=='td'){
    est=isoph.td(START=START, STOP=STOP, STATUS=STATUS, Z=Z, shape=shape, K=K, maxiter=maxiter, eps=eps, maxdec=maxdec)
  }
    
  est$call=match.call()
  est$formula=formula
  
  class(est)="isoph"
  
  est
}
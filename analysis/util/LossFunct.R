#########################################################################
Integrate <- function(f, time.pnt){
  f.value <- sapply(time.pnt, f)
  result <- diff(time.pnt)%*%(f.value[-length(f.value)] + f.value[-1])/2
  return(result)
}
#########################################################################

Loss.func <- function(Curve, KM, T.pnt){
  suppressWarnings(rm(Extract.surv))
  
  if(substr(Curve,1,3) == "Hot"){
    Mod <- substring(sub("_.*$","",Curve),9)
    Par <- as.numeric(substring(sub("^[^_]*", "", Curve),2))
    
    switch(Mod,
           Exp = {
             Extract.surv <- function(t){1-pexp(t,rate = exp(Par))}
           },
           WebI = {
             Extract.surv <- function(t){1-pweibull(t,shape = 2, scale = 10*exp(Par))}
           },
           WebD = {
             Extract.surv <- function(t){1-pweibull(t,shape = 0.5, scale = 5*exp(Par))}
           },
           Bat = {
             Extract.surv <- function(t){getSurv.Bath(a = exp(Par), t)}
           },
           Lgn = {
             Extract.surv <- function(t){1-plnorm(t, 1.5, exp(Par))}
           }
    )
    
  }else{
    switch(Curve, 
           Bath.T1 = {
             Extract.surv <- function(t){getSurv.Bath(a = 0.01, t)}
           },
           Bath.T2 = {
             Extract.surv <- function(t){getSurv.Bath(a = 0.15, t)}
           },
           Bath.T3 = {
             Extract.surv <- function(t){getSurv.Bath(a = 0.2, t)}
           },
           Bath.T4 = {
             Extract.surv <- function(t){getSurv.Bath(a = 0.9, t)}
           },
           Exp.T1 = {
             Extract.surv <- function(t){1-pexp(t,rate = 0.1)}
           },
           Exp.T2 = {
             Extract.surv <- function(t){1-pexp(t,rate = 0.23)}
           },
           Exp.T3 = {
             Extract.surv <- function(t){1-pexp(t,rate = 0.4)}
           },
           Exp.T4 = {
             Extract.surv <- function(t){1-pexp(t,rate = 0.9)}
           },
           WI.T1 = {
             Extract.surv <- function(t){1-pweibull(t,3,10)}
           },
           WI.T2 = {
             Extract.surv <- function(t){1-pweibull(t,3,6.2)}
           },
           WI.T3 = {
             Extract.surv <- function(t){1-pweibull(t,3,4.3)}
           },
           WI.T4 = {
             Extract.surv <- function(t){1-pweibull(t,3,2)}
           },
           WD.T1 = {
             Extract.surv <- function(t){1-pweibull(t,0.9,7)}
           },
           WD.T2 = {
             Extract.surv <- function(t){1-pweibull(t,0.9,3)}
           },
           WD.T3 = {
             Extract.surv <- function(t){1-pweibull(t,0.9,2.5)}
           },
           WD.T4 = {
             Extract.surv <- function(t){1-pweibull(t,0.9,1)}
           },
           Log.T1 = {
             Extract.surv <- function(t){1-plnorm(t, 2.0, 0.3)}
           },
           Log.T2 = {
             Extract.surv <- function(t){1-plnorm(t, 1.7, 0.2)}
           },
           Log.T3 = {
             Extract.surv <- function(t){1-plnorm(t, 1.3, 0.3)}
           },
           Log.T4 = {
             Extract.surv <- function(t){1-plnorm(t, 0.5, 0.5)}
           }
    )## end of switch
  }
  
  ##
  if ("survfit" %in% class(KM)){
    
    integrand <- function(t){
      ( ipred::getsurv(KM,t) - Extract.surv(t))^2
    }
    
    Int.value <- Integrate(f = integrand, time.pnt = T.pnt)
    return(Int.value/diff(range(T.pnt)))
    
  } else if (class(KM)== "icfit"){
    
    integrand <- function(t){
      (interval::getsurv(t,KM,nonUMLE.method = "right")[[1]]$S - Extract.surv(t))^2
    }
    
    Int.value <- Integrate(f = integrand, time.pnt = T.pnt)
    return(Int.value/diff(range(T.pnt)))
    
  } else if (class(KM)=="sp_curves"||"rfsrc" %in% class(KM)||class(KM) == "ranger"){
    integrand <- function(t){
      (sapply(t, Interpolate, Curve = KM ) - Extract.surv(t))^2
    }
    
    Int.value <- Integrate(f = integrand, time.pnt = T.pnt)
    return(Int.value/diff(range(T.pnt)))
    
  } else if (class(KM) == "list"){
    integrand <- function(t){
      (sapply(t, Interpolate, Curve = KM ) - Extract.surv(t))^2
    }
    
    Int.value <- Integrate(f = integrand, time.pnt = T.pnt)
    return(Int.value/diff(range(T.pnt)))    
  } else if (class(KM) == "ic_ph"||class(KM) == "ic_np"){
    integrand <- function(t){
      (1-icenReg::getFitEsts(fit = KM, q = t) - Extract.surv(t))^2
    }
    
    Int.value <- Integrate(f = integrand, time.pnt = T.pnt)
    return(Int.value/diff(range(T.pnt)))   
  } else if (class(KM) == "numeric"){
    f.value <- (KM - Extract.surv(T.pnt))^2
    Int.value <- diff(T.pnt)%*%(f.value[-length(f.value)] + f.value[-1])/2
    return(Int.value/diff(range(T.pnt)))
  } else{
    stop("wrong class of NPMLE object specified")
  }
}## end of function
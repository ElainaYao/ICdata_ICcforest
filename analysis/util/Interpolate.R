######################################################################################################
Interpolate <- function(Curve,X){
  if(class(Curve)[1]=="sp_curves"){
    if(length(X)!=1)
      stop("X can only be single number")
    
    if (!inherits(Curve, "sp_curves")) 
      stop("obj is not of class sp_curves")
    
    if(length(Curve$S_curves)!=1)
      stop("only one row (one predicted curve) at a time")
    
    intervals <- Curve$Tbull_ints
    surv.rates <- as.vector(unlist(Curve$S_curves))
    surv.rates[is.nan(surv.rates)] = 0
    if(nrow(intervals)!=length(surv.rates))
      stop("dim of matrix and survivals doesn't match")
    
    N <- nrow(intervals)
    Id <- findInterval(X,intervals[,2])
    
    if(Id == N){
      result = surv.rates[N]
    }else if(Id == 0){
      result = approx(c(0,intervals[,2][1]),c(1,surv.rates[1]),xout = X)$y
    }else{
      a <- intervals[,2][Id]
      b <- intervals[,1][Id+1]
      c <- intervals[,2][Id+1]
      s.high <- surv.rates[Id]
      s.low <- surv.rates[Id+1]
      
      if(X>=a && X<b){
        result = s.high
      }else if(X>=b && X<c){
#        if(Id == N-1){
#          result = surv.rates[N]
#        }else{
          result = approx(c(b,c),c(s.high,s.low),xout = X)$y
#        }
      }else{
        stop("something wrong with the interpolating function")
      }
    }
  } else if("rfsrc" %in% class(Curve)) {
#    result = approx(Curve$time.interest, Curve$survival[Curve$i,], xout = X, method="constant",
#                    yleft = 1, yright = 0, f=0)$y
    N = length(Curve$time.interest)
    if (X < Curve$time.interest[1]){
      result = approx(c(0,Curve$time.interest[1]),c(1,Curve$survival[Curve$i,1]),xout = X)$y
    } else if (X > Curve$time.interest[N]){
      result = Curve$survival[Curve$i,N]
    } else {
      result = approx(Curve$time.interest, Curve$survival[Curve$i,], xout = X)$y
    }
    
  } else if ("ranger" %in% class(Curve)){
#    result = approx(Curve$unique.death.times, Curve$survival[Curve$i,], xout = X, method="constant",
#                    yleft = 1, yright = 0, f=0)$y
    N = length(Curv$unique.death.times)
    if (X < Curve$unique.death.times[1]){
      result = approx(c(0,Curve$unique.death.times[1]),c(1,Curve$survival[Curve$i,1]),xout = X)$y
    } else if (X > Curve$unique.death.times[N]){
      result = Curve$survival[Curve$i,N]
    } else {
      result = approx(Curve$unique.death.times, Curve$survival[Curve$i,], xout = X)$y
    }
  } else if ("icfit" %in% class(Curve)){
    pt <- as.numeric(Curve$intmap)
    N <- length(pt)
    s_pt1 <- matrix(0,nrow = N,1)
    
    s_pt0 <- lapply(2:(N-1), function(l){
      ll = floor(l/2)
      resspt = 1-sum(resCurveI1$pf[1:ll])
      return(resspt)
    })
    s_pt1[2:(N-1)] <- as.numeric(s_pt0)
    s_pt1[1] = 1
    
    if (X >= pt[N]){
      result = s_pt1[N]
    } else {
      Id = which(X<pt)
      if (Id[1] == 1){
        result = approx(c(0,pt[1]),c(1,s_pt1[1]),xout = X)$y
      } else {
        result = approx(c(pt[Id[1]-1],pt[Id[1]]),c(s_pt1[Id[1]-1],s_pt1[Id[1]]),xout = X)$y 
      }
    }
  } else {
    stop("wrong class of Interpolcate object specified")
  }
  return(result) 
}

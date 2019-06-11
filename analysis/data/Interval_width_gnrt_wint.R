###########################################################################################
Hothorn_gnrt_width_wint <- function(n=200, Model.type = 2, Dist = "Exp", Censor = 0, tt){
  
  Data <- as.data.frame(matrix(NA,n,13))
  names(Data)<-c("L","R",
                 "X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","T")
  Count = 1
  
  while(Count <= n){
    x1 <- runif(1,0,1)
    x2 <- sample(c(0,1),1)
    x3 <- sample(c(0,1),1)
    x4 <- runif(1,0,1)
    x5 <- runif(1,0,1)
    x6 <- sample(c(0,1),1)
    x7 <- runif(1,0,1)
    x8 <- sample(c(0,1),1)
    x9 <- sample(c(0,1),1)
    x10 <- runif(1,0,1)
    
    if(Model.type == 3){
      Param <- -(cos((x1+x2)*pi)+sqrt(x1+x2))
    }else if(Model.type == 2){
      Param <- -x1-x2
    }else {
      stop("Wrong model type: It's either 2(second setup) or 3 (third setup)")
    }
    
    if(Dist == "Exp"){
      t <- rexp(1,exp(Param))
      Data[Count,"Class"] <- paste("Hothorn-Exp_",Param,sep="")
      if(Censor == 0){
        K = 100
        K.S = 800 
      }else if(Censor == 1){
        K = 4
        K.S = 9 
      }else if(Censor == 2){
        K = 2
        K.S = 5 
      }else{
        stop("wrong censoring type")
      }
    }else if(Dist == "WI"){
      t <- rweibull(1, shape = 2, scale = 10*exp(Param))
      Data[Count,"Class"] <- paste("Hothorn-WebI_",Param,sep="")
      if(Censor == 0){
        K = 100
        K.S = 800 
      }else if(Censor == 1){
        K = 5
        K.S = 11 
      }else if(Censor == 2){
        K = 3
        K.S = 7 
      }else{
        stop("wrong censoring type")
      }
    }else if(Dist == "WD"){
      t <- rweibull(1, shape = 0.5, scale = 5*exp(Param))
      Data[Count,"Class"] <- paste("Hothorn-WebD_",Param,sep="")
      if(Censor == 0){
        K = 200
        K.S = 800 
      }else if(Censor == 1){
        K = 4
        K.S = 9 
      }else if(Censor == 2){
        K = 1
        K.S = 2
      }else{
        stop("wrong censoring type")
      }
    }else if (Dist == "Bat"){
      t <- Bathtub(n=1, a = exp(Param))
      Data[Count,"Class"] <- paste("Hothorn-Bat_",Param,sep="")
      if(Censor == 0){
        K = 200
        K.S = 800 
      }else if(Censor == 1){
        K = 4
        K.S = 9 
      }else if(Censor == 2){
        K = 1
        K.S = 2
      }else{
        stop("wrong censoring type")
      }
    }else if (Dist == "Lgn"){
      t <- rlnorm(1, meanlog = 1.5, sdlog = exp(Param))
      Data[Count,"Class"] <- paste("Hothorn-Lgn_",Param,sep="")
      if(Censor == 0){
        K = 200
        K.S = 800 
      }else if(Censor == 1){
        K = 4
        K.S = 9 
      }else if(Censor == 2){
        K = 1
        K.S = 2
      }else{
        stop("wrong censoring type")
      }
    }else{
      print("Error: Wrong distribution!");
      return(0);
    }
    
    if (tt ==1){
      obs <- round(c(0,cumsum(runif(K.S,0.30,0.35)),Inf),digit = 2)
    } else if (tt ==2){
      obs <- round(c(0,cumsum(runif(K.S,0.45,0.65)),Inf),digit = 2)
    } else if (tt ==3){
      obs <- round(c(0,cumsum(runif(K.S,0.75,0.95)),Inf),digit = 2)
    } else if (tt == 4){
      obs <- round(c(0,cumsum(runif(K.S,1.05,1.25)),Inf),digit = 2)
    } else if (tt == 5){
      obs <- round(c(0,cumsum(runif(K.S,1.35,1.55)),Inf),digit = 2)
    } else if (tt == 6){
      obs <- round(c(0,cumsum(runif(K.S,1.65,1.85)),Inf),digit = 2)
    } else if (tt == 7){
      obs <- round(c(0,cumsum(runif(K.S,1.95,2.15)),Inf),digit = 2)
    }

    ##=============================================
    Data[Count,"X1"] <- x1
    Data[Count,"X2"] <- x2
    Data[Count,"X3"] <- x3
    Data[Count,"X4"] <- x4
    Data[Count,"X5"] <- x5
    Data[Count,"X6"] <- x6
    Data[Count,"X7"] <- x7
    Data[Count,"X8"] <- x8
    Data[Count,"X9"] <- x9
    Data[Count,"X10"] <- x10
    Data[Count,"T"] <- t
    
    Id <- findInterval(t,obs)
    Data[Count,"L"]<- obs[Id]
    Data[Count,"R"]<- obs[Id+1]
    
    Count = Count + 1
  }#end while 
  
  return(Data)
}

###############################
IC.generate_wint <- function(n=200, Dist = "Exp", Censor = 0 , tt){
  Data <- as.data.frame(matrix(NA,n,6))
  names(Data)<-c("L","R",
                 "X1","X2","X3","T")
  Count = 1
  
  while(Count <= n){
    x1 <- sample(1:5,1)
    x2 <- sample(c(1,2),1) 
    x3 <- runif(1,0,2)
    ##============================================
    if(Dist == "Bat"){
      if(x1 < 2.5){
        if(x2 == 1){
          t <- Bathtub(n=1, a = 0.01)
          Data[Count,"Class"] <- "Bath.T1"
        }else{
          t <- Bathtub(n=1, a = 0.15)
          Data[Count,"Class"] <- "Bath.T2"
        }
      }else{
        if(x3 <= 1){
          t <- Bathtub(n=1, a = 0.2)
          Data[Count,"Class"] <- "Bath.T3"
        }else{
          t <- Bathtub(n=1, a = 0.9)
          Data[Count,"Class"] <- "Bath.T4"
        }
      }
      
      if(Censor == 0){
        K = 100
        K.S = 800 
      }else if(Censor == 1){
        K = 3
        K.S = 7 
      }else if(Censor == 2){
        K = 2
        K.S = 5 
      }else{
        stop("wrong censoring type")
      }
    }else if(Dist == "Exp"){
      if(x1 < 2.5){
        if(x2 == 1){
          t <- rexp(1,0.1)
          Data[Count,"Class"] <- "Exp.T1"
        }else{
          t <- rexp(1,0.23)
          Data[Count,"Class"] <- "Exp.T2"
        }
      }else{
        if(x3 <= 1){
          t <- rexp(1,0.4)
          Data[Count,"Class"] <- "Exp.T3"
        }else{
          t <- rexp(1,0.9)
          Data[Count,"Class"] <- "Exp.T4"
        }
      }
      
      if(Censor == 0){
        K = 100
        K.S = 800
      }else if(Censor == 1){
        K = 5
        K.S = 12 
      }else if(Censor == 2){
        K = 2
        K.S = 5 
      }else{
        stop("wrong censoring type")
      }
    }else if(Dist == "WD"){
      if(x1 < 2.5){
        if(x2 == 1){
          t <- rweibull(1, 0.9, 7)
          Data[Count,"Class"] <- "WD.T1"
        }else{
          t <- rweibull(1, 0.9, 3)
          Data[Count,"Class"] <- "WD.T2"
        }
      }else{
        if(x3 <= 1){
          t <- rweibull(1, 0.9, 2.5)
          Data[Count,"Class"] <- "WD.T3"
        }else{
          t <- rweibull(1, 0.9, 1)
          Data[Count,"Class"] <- "WD.T4"
        }
      }
      
      if(Censor == 0){
        K = 100
        K.S = 800 
      }else if(Censor == 1){
        K = 4
        K.S = 9 
      }else if(Censor == 2){
        K = 2
        K.S = 5 
      }else{
        stop("wrong censoring type")
      }
    }else if(Dist == "WI"){
      if(x1 < 2.5){
        if(x2 == 1){
          t <- rweibull(1, 3, 10)
          Data[Count,"Class"] <- "WI.T1"
        }else{
          t <- rweibull(1, 3, 6.2)
          Data[Count,"Class"] <- "WI.T2"
        }
      }else{
        if(x3 <= 1){
          t <- rweibull(1, 3, 4.3)
          Data[Count,"Class"] <- "WI.T3"
        }else{
          t <- rweibull(1, 3, 2)
          Data[Count,"Class"] <- "WI.T4"
        }
      }
      
      if(Censor == 0){
        K = 100
        K.S = 800 
      }else if(Censor == 1){
        K = 6
        K.S = 14 
      }else if(Censor == 2){
        K = 4
        K.S = 9
      }else{
        stop("wrong censoring type")
      }
    }else if(Dist == "Lgn"){
      if(x1 < 2.5){
        if(x2 == 1){
          t <- rlnorm(1, meanlog = 2.0, sdlog = 0.3)
          Data[Count,"Class"] <- "Log.T1"
        }else{
          t <- rlnorm(1, meanlog = 1.7, sdlog = 0.2)
          Data[Count,"Class"] <- "Log.T2"
        }
      }else{
        if(x3 <= 1){
          t <- rlnorm(1, meanlog = 1.3, sdlog = 0.3)
          Data[Count,"Class"] <- "Log.T3"
        }else{
          t <- rlnorm(1, meanlog = 0.5, sdlog = 0.5)
          Data[Count,"Class"] <- "Log.T4"
        }
      }
      
      if(Censor == 0){
        K = 100
        K.S = 800 
      }else if(Censor == 1){
        K = 5
        K.S = 12 
      }else if(Censor == 2){
        K = 4
        K.S = 9
      }else{
        stop("wrong censoring type")
      }
    }
    
    if (tt ==1){
      obs <- round(c(0,cumsum(runif(K.S,0.15,0.35)),Inf),digit = 2)
    } else if (tt ==2){
      obs <- round(c(0,cumsum(runif(K.S,0.45,0.65)),Inf),digit = 2)
    } else if (tt ==3){
      obs <- round(c(0,cumsum(runif(K.S,0.75,0.95)),Inf),digit = 2)
    } else if (tt == 4){
      obs <- round(c(0,cumsum(runif(K.S,1.05,1.25)),Inf),digit = 2)
    } else if (tt == 5){
      obs <- round(c(0,cumsum(runif(K.S,1.35,1.55)),Inf),digit = 2)
    } else if (tt == 6){
      obs <- round(c(0,cumsum(runif(K.S,1.65,1.85)),Inf),digit = 2)
    } else if (tt == 7){
      obs <- round(c(0,cumsum(runif(K.S,0.01,0.05)),Inf),digit = 2)
    } else if (tt == 8){
      obs <- round(c(0,cumsum(runif(K.S,0.55,0.65)),Inf),digit = 2)
    } else if (tt == 9){
      obs <- round(c(0,cumsum(runif(K.S,1.05,1.15)),Inf),digit = 2)
    } else if (tt == 10){
      obs <- round(c(0,cumsum(runif(K.S,5.55,5.65)),Inf),digit = 2)
    } else if (tt == 11){
      obs <- round(c(0,cumsum(runif(K.S,10.55,10.65)),Inf),digit = 2)
    } else if (tt == 12){
      obs <- round(c(0,cumsum(runif(K.S,15.55,15.65)),Inf),digit = 2)
    } else if (tt == 13){
      obs <- round(c(0,cumsum(runif(K.S,20.55,20.65)),Inf),digit = 2)
    }
    ##=============================================
    Data[Count,"X1"] <- x1
    Data[Count,"X2"] <- x2
    Data[Count,"X3"] <- x3
    Data[Count,"T"] <- t
    
    Id <- findInterval(t,obs)
    Data[Count,"L"]<- obs[Id]
    Data[Count,"R"]<- obs[Id+1]
    
    Count = Count + 1
  }
  Data$X4 <- sample(1:5,n,replace = TRUE)
  Data$X5 <- sample(c(1,2),n,replace = TRUE)
  Data$X6 <- runif(n,0,2)
  Data$X7 <- sample(1:5,n,replace = TRUE)
  Data$X8 <- sample(c(1,2),n,replace = TRUE)
  Data$X9 <- runif(n,0,2)
  Data$X10 <- runif(n,0,2)
  return(Data)
}
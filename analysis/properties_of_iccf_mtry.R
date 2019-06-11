library(ipred)
library(survival)
library(icenReg)
library(LTRCtrees)
library(ICcforest)

source("bathtub.R")
source("Interval_width_gnrt_wint.R")
source("LossFunct.R")
source("Interpolate.R")

##############################################################################################
# === Input of the function Pred_width():
##### distribtuion:
##### --- "Bat": Bathtub
##### --- "Exp": Exponential
##### --- "Lgn": Log-normal
##### --- "WI": Weibull-increasing 
##### --- "WD": Weibull-decreasing 
##### model: underlying survival relationship
##### --- 1: tree structure
##### --- 2: linear survival relationship
##### --- 3: nonlinear survival relationship
##### C.rate: right-censoring rate
##### --- 0:  0%
##### --- 1: 20%
##### --- 2: 40%
##### --- 3: 60%
##### tt: interval width, 1-G1; 3-G2; 6-G3
##### M: number of simulations

# === Output of the function Pred_width():
##### mtry: L2 errors of IC cforest with different mtry's and tuned mtry, 
#####       with minsplit, minprob, minbucket set by default

Pred_funct <- function(Nn, distribution, model, C.rate, tt, M){
  C.rate0=0
  RES.L2 <- NULL
  RES.L2$mtry <- data.frame(matrix(0, nrow = M, ncol = 7))
  names(RES.L2$mtry) <- c("1","3","4","6","9","10","T")

  set.seed(101)
  sampleID <- sort(sample(100000000,M))
  for (mm in 1:M){
    mm1 = sampleID[mm]
    set.seed(mm1)
    ## create the simulation dataset 
    if(model == 1){
      DATA <- IC.generate_wint(n=Nn, Dist = distribution, Censor = C.rate0, tt)
    }else{
      DATA <- Hothorn_gnrt_width_wint(n = Nn, Model.type = model, Dist = distribution, Censor = C.rate0, tt)
    }
    
    idx_inf <- (DATA$R == Inf)
    DATA$R[idx_inf] <- 999.
    
    ## make right-censoring 
    if (C.rate == 1){
      DATA$R[sample(Nn,round(Nn*0.20))]<-999.
    } else if (C.rate == 2){
      DATA$R[sample(Nn,round(Nn*0.40))]<-999.
    } else if (C.rate == 3){
      DATA$R[sample(Nn,round(Nn*0.60))]<-999.
    }
    
    ## time points of interest to evaluate the integral
    time.uniq <- unique(sort(c(DATA$T,DATA$L,DATA$R))) 
    time.uniq <- time.uniq[time.uniq <= max(DATA$T)]
    
    Formula = as.formula(Surv(L,R,type="interval2")~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10)
    
    ############################## ------------ mtry ------------- ##############################
    mtrypool <- c(1,3,4,6,9,10)
    for (j in 1:length(mtrypool)){
      IC.cforest <- ICcforest(Formula, data = DATA, mtry = mtrypool[j], 
                              Control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                                minsplit = 20L, minprob = 0.01, minbucket = 7L,
                                                                mincriterion = 0))
      print(sprintf("L2 - IC cforest mtry = %1.0f...",mtrypool[j]))
      Pred.IC.cforest <- predict(IC.cforest, type="prob")
      L2 <- c()
      for(i in 1:nrow(DATA)){
        Km <- Pred.IC.cforest[[i]]
        Cur <- DATA[i,"Class"]
        L2 <- c(L2, Loss.func(Cur, Km, time.uniq))
      }
      RES.L2$mtry[mm,j] <- mean(L2)
      rm(IC.cforest)
      rm(Pred.IC.cforest)
    }    
    
    print(sprintf("L2 - IC cforest mtry is being tuned..."))
    IC.cforest <- ICcforest(Formula, data = DATA, 
                            Control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                              minsplit = 20L, minprob = 0.01, minbucket = 7L,
                                                              mincriterion = 0))
    Pred.IC.cforest <- predict(IC.cforest, type="prob")
    L2 <- c()
    for(i in 1:nrow(DATA)){
      Km <- Pred.IC.cforest[[i]]
      Cur <- DATA[i,"Class"]
      L2 <- c(L2, Loss.func(Cur, Km, time.uniq))
    }
    RES.L2$mtry[mm,7] <- mean(L2)
    rm(IC.cforest)
    rm(Pred.IC.cforest)
    print("=======IC.cforest for different mtry's are done...") 
  }
  return(RES.L2)
}  


############### === Comparison of IC cforest with different mtry's and mtry tuned ==== ################
##### N = 200, no right-censoring, censoring interval width generated from G1
L2.Bat.m1.c0.tt1 <- Pred_funct(Nn = 200, distribution = "Bat", model = 1, C.rate = 0, tt = 1, M = 500)
L2.Bat.m2.c0.tt1 <- Pred_funct(Nn = 200, distribution = "Bat", model = 2, C.rate = 0, tt = 1, M = 500)
L2.Bat.m3.c0.tt1 <- Pred_funct(Nn = 200, distribution = "Bat", model = 3, C.rate = 0, tt = 1, M = 500)
L2.Exp.m1.c0.tt1 <- Pred_funct(Nn = 200, distribution = "Exp", model = 1, C.rate = 0, tt = 1, M = 500)
L2.Exp.m2.c0.tt1 <- Pred_funct(Nn = 200, distribution = "Exp", model = 2, C.rate = 0, tt = 1, M = 500)
L2.Exp.m3.c0.tt1 <- Pred_funct(Nn = 200, distribution = "Exp", model = 3, C.rate = 0, tt = 1, M = 500)
L2.lgn.m1.c0.tt1 <- Pred_funct(Nn = 200, distribution = "Lgn", model = 1, C.rate = 0, tt = 1, M = 500)
L2.Lgn.m2.c0.tt1 <- Pred_funct(Nn = 200, distribution = "Lgn", model = 2, C.rate = 0, tt = 1, M = 500)
L2.Lgn.m3.c0.tt1 <- Pred_funct(Nn = 200, distribution = "Lgn", model = 3, C.rate = 0, tt = 1, M = 500)
L2.WD.m1.c0.tt1 <- Pred_funct(Nn = 200, distribution = "WD", model = 1, C.rate = 0, tt = 1, M = 500)
L2.WD.m2.c0.tt1 <- Pred_funct(Nn = 200, distribution = "WD", model = 2, C.rate = 0, tt = 1, M = 500)
L2.WD.m3.c0.tt1 <- Pred_funct(Nn = 200, distribution = "WD", model = 3, C.rate = 0, tt = 1, M = 500)
L2.WI.m1.c0.tt1 <- Pred_funct(Nn = 200, distribution = "WI", model = 1, C.rate = 0, tt = 1, M = 500)
L2.WI.m2.c0.tt1 <- Pred_funct(Nn = 200, distribution = "WI", model = 2, C.rate = 0, tt = 1, M = 500)
L2.WI.m3.c0.tt1 <- Pred_funct(Nn = 200, distribution = "WI", model = 3, C.rate = 0, tt = 1, M = 500)

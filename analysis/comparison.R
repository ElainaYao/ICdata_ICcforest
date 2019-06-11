library(ipred)
library(survival)
library(icenReg)
library(LTRCtrees)
library(ICcforest)

source("Interval_width_gnrt_wint.R")
source("LossFunct.R")
source("Interpolate.R")

##############################################################################################
.pred_Surv_new <- function(y, w) {
  if (length(y) == 0) return(NA)
  idx = which(w>0)
  y = y[idx]
  w = w[idx]
  yy = as.matrix(y)[,1:2]
  ic_np(yy, weights = w)
}

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
##### iccfD: L2 errors of IC cforest with default settings over M simulations
##### iccfT: L2 errors of IC cforest with tuned mtry and with minsplit, minbucket, minprob set by 
#####        "15%-Default-6% Rule" over M simulations.
##### iccx:  L2 errors of IC Cox over M simulations.
##### ictr:  L2 errors of IC ctree over M simulations.

Pred_funct <- function(Nn, distribution, model, C.rate, tt, M){
  C.rate0=0
  RES.L2 <- NULL
  RES.L2$iccfD <- rep(0,M)
  RES.L2$iccfT <- rep(0,M)
  RES.L2$iccx <- rep(0,M)
  RES.L2$ictr <- rep(0,M)
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
    
    ###############################---------L2 for ICcforest with default settings --------###########################
    ## create a IC cforest object
    IC.cforest <- ICcforest(formula = Formula, data = DATA, mtry = ceiling(sqrt(10)),
                            control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                              minsplit = 20L, minprob = 0.01, minbucket = 7L,
                                                              mincriterion = 0))
    ## make prediction
    Pred.IC.cforest <- predict(IC.cforest, type="prob")
    ## compute the L2 error
    L2 <- c()
    for(i in 1:nrow(DATA)){
      Km <- Pred.IC.cforest[[i]]
      Cur <- DATA[i,"Class"]
      L2 <- c(L2, Loss.func(Cur, Km, time.uniq))
    }
    RES.L2$iccfD[mm] <- mean(L2)
    print("=======IC.cforest is done...")
    
    ###############################---------L2 for ICcforest with mtry tuned--------###########################
    ## create a IC cforest object with tuned mtry
    IC.cforest <- ICcforest(formula = Formula, data = DATA)
    ## make prediction
    Pred.IC.cforest <- predict(IC.cforest, type="prob")
    ## compute the L2 error
    L2 <- c()
    for(i in 1:nrow(DATA)){
      Km <- Pred.IC.cforest[[i]]
      Cur <- DATA[i,"Class"]
      L2 <- c(L2, Loss.func(Cur, Km, time.uniq))
    }
    RES.L2$iccfT[mm] <- mean(L2)
    print("=======IC.cforest is done...")
    
    ###############################---------L2 for IC.tree--------###########################
    ## create a IC ctree object
    IC.tree <- LTRCtrees::ICtree(Formula, DATA)
    km.IC0 <- icenReg::ic_np(DATA[,c("L","R")])
    if(length(IC.tree)==1){
      L2 <- c()
      for(i in 1:nrow(DATA)){ ## computing L2
        Km <- km.IC0
        Cur <- DATA[i,"Class"]
        L2 <- c(L2, Loss.func(Cur, Km, time.uniq)) ## Loss.func is the function computing the integrated L2 loss
      }
      RES.L2$ictr[mm] <- mean(L2)
    }else{
      ## make prediction
      Pred.IC.tree <- predict(IC.tree, type="prob",FUN=.pred_Surv_new)
      L2 <- c()
      for(i in 1:nrow(DATA)){
        Km <- Pred.IC.tree[[i]]
        Cur <- DATA[i,"Class"]
        L2 <- c(L2, Loss.func(Cur, Km, time.uniq))
      }
      RES.L2$ictr[mm] <- mean(L2)
    }
    print("=======IC.tree done...")
    
    #####################------------L2 Score for Cox_IC----------###########################
    ## Fit the Cox model
    Cox_IC <- ic_sp(Formula, DATA)
    L2 <- c() 
    for(i in 1:nrow(DATA)){
      Km <- icenReg::getSCurves(Cox_IC, DATA[i,])
      Cur <- DATA[i,"Class"]
      L2 <- c(L2, Loss.func(Cur, Km, time.uniq))
    }
    RES.L2$iccx[mm] <- mean(L2)
    print("=======IC.Cox done...")
  }
  return(RES.L2)
}  


##### === Comparison of IC Cox, IC ctree, IC cforest with parameters set by default, 
#####                   IC cforest with mtry tuned, parameters set by "15%-Def-6% Rule" ==== ######
##### N = 200, Light right-censoring (20%), censoring interval width generated from G1
L2.Bat.m1.c1.tt1 <- Pred_funct(Nn = 200, distribution = "Bat", model = 1, C.rate = 1, tt = 1, M = 500)
L2.Bat.m2.c1.tt1 <- Pred_funct(Nn = 200, distribution = "Bat", model = 2, C.rate = 1, tt = 1, M = 500)
L2.Bat.m3.c1.tt1 <- Pred_funct(Nn = 200, distribution = "Bat", model = 3, C.rate = 1, tt = 1, M = 500)
L2.Exp.m1.c1.tt1 <- Pred_funct(Nn = 200, distribution = "Exp", model = 1, C.rate = 1, tt = 1, M = 500)
L2.Exp.m2.c1.tt1 <- Pred_funct(Nn = 200, distribution = "Exp", model = 2, C.rate = 1, tt = 1, M = 500)
L2.Exp.m3.c1.tt1 <- Pred_funct(Nn = 200, distribution = "Exp", model = 3, C.rate = 1, tt = 1, M = 500)
L2.lgn.m1.c1.tt1 <- Pred_funct(Nn = 200, distribution = "Lgn", model = 1, C.rate = 1, tt = 1, M = 500)
L2.Lgn.m2.c1.tt1 <- Pred_funct(Nn = 200, distribution = "Lgn", model = 2, C.rate = 1, tt = 1, M = 500)
L2.Lgn.m3.c1.tt1 <- Pred_funct(Nn = 200, distribution = "Lgn", model = 3, C.rate = 1, tt = 1, M = 500)
L2.WD.m1.c1.tt1 <- Pred_funct(Nn = 200, distribution = "WD", model = 1, C.rate = 1, tt = 1, M = 500)
L2.WD.m2.c1.tt1 <- Pred_funct(Nn = 200, distribution = "WD", model = 2, C.rate = 1, tt = 1, M = 500)
L2.WD.m3.c1.tt1 <- Pred_funct(Nn = 200, distribution = "WD", model = 3, C.rate = 1, tt = 1, M = 500)
L2.WI.m1.c1.tt1 <- Pred_funct(Nn = 200, distribution = "WI", model = 1, C.rate = 1, tt = 1, M = 500)
L2.WI.m2.c1.tt1 <- Pred_funct(Nn = 200, distribution = "WI", model = 2, C.rate = 1, tt = 1, M = 500)
L2.WI.m3.c1.tt1 <- Pred_funct(Nn = 200, distribution = "WI", model = 3, C.rate = 1, tt = 1, M = 500)

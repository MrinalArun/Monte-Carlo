library(OptionPricing)

SingleControl <- function(ZArr,Mean,ArrCalls,Strikes){
  VarOut <- vector('numeric')
  Means  <- vector('numeric')
  for(i in c(1:length(Strikes))){
    cov = cov(ZArr,ArrCalls[i,])
    var = var(ZArr)
    c   = -cov/var
    Means = c(Means,mean(ArrCalls[i,]+c*(ZArr-Mean)))
    VarOut  = c(VarOut,var(ArrCalls[i,]+c*(ZArr-Mean)))
  }
  return(c(Means,VarOut))
}

SingleControlCall <- function(YTCalls,Mean,ArrCalls,Strikes){
  VarOut <- vector('numeric')
  Means  <- vector('numeric')
  for(i in c(1:length(Strikes))){
    cov = cov(YTCalls[i,],ArrCalls[i,])
    var = var(YTCalls[i,])
    c   = -cov/var
    Means = c(Means,mean(ArrCalls[i,]+c*(YTCalls[i,]-Mean[i])))
    VarOut  = c(VarOut,var(ArrCalls[i,]+c*(YTCalls[i,]-Mean[i])))
  }
  return(c(Means,VarOut))
  
}

Paths <- 10000
Time  <- 3
S0    <- 100
r     <- 0.04

Lambdas <- c(1,3)
Sigmas   <- c(0.25,0.75)
Strikes  <- c(50,100,150)
CallPricesOld <- matrix(0,nrow = Paths,ncol = length(Strikes))
CallPricesNew <- matrix(0,nrow = Paths,ncol = length(Strikes))

ArrST       <- vector('numeric')
ArrCall1    <- vector('numeric')
ArrCall2    <- vector('numeric')
ArrCall3    <- vector('numeric')
ArrCalls    <- vector('numeric')
ArrYT1      <- vector('numeric')
ArrYT2      <- vector('numeric')
ArrYT1Sq    <- vector('numeric')
ArrYT2Sq    <- vector('numeric')
ArrYT1Exp   <- vector('numeric')
ArrYT2Exp   <- vector('numeric')
ArrYT1Call  <- vector('numeric')
ArrYT2Call  <- vector('numeric')
ArrYT1Call1 <- vector('numeric')
ArrYT2Call1 <- vector('numeric')
ArrYT1Call2 <- vector('numeric')
ArrYT2Call2 <- vector('numeric')
ArrYT1Call3 <- vector('numeric')
ArrYT2Call3 <- vector('numeric')

StatsST      <- vector('numeric')
StatsYT1     <- vector('numeric')
StatsYT2     <- vector('numeric')
StatsYT1Sq   <- vector('numeric')
StatsYT2Sq   <- vector('numeric')
StatsYT1Exp  <- vector('numeric')
StatsYT2Exp  <- vector('numeric')
StatsMulti1  <- vector('numeric')
StatsMulti2  <- vector('numeric')
StatsMulti3  <- vector('numeric')
StatsMulti4  <- vector('numeric')

MeanYT1Call <- vector('numeric')
MeanYT2Call <- vector('numeric')

MeanYT1    <- (r-(Sigmas[1]*Sigmas[1]/2))*Time
MeanYT2    <- (r-(Sigmas[2]*Sigmas[2]/2))*Time
MeanYT1Sq  <- (((r-(Sigmas[1]*Sigmas[1]/2))*Time)^2)+(Sigmas[1]*Sigmas[1]*Time)
MeanYT2Sq  <- (((r-(Sigmas[2]*Sigmas[2]/2))*Time)^2)+(Sigmas[2]*Sigmas[2]*Time)
MeanYT1Exp <- exp(r*Time)
MeanYT2Exp <- exp(r*Time)

for(x in 1:length(Strikes)){
  MeanYT1Call <- c(MeanYT1Call,BS_EC( T = Time, K = Strikes[x], r = r, sigma = Sigmas[1], S0 = S0 )[1][["price"]])
  MeanYT2Call <- c(MeanYT2Call,BS_EC( T = Time, K = Strikes[x], r = r, sigma = Sigmas[2], S0 = S0 )[1][["price"]])
}

for(i in 1:Paths){
  t      <- 0
  T      <- Time
  SNew   <- S0
  YT1    <- 0
  YT2    <- 0
  YT1Exp <- 1
  YT2Exp <- 1
  State  <- 1 # 1 corresponds to low vol. 2 corresponds to high vol
  
  while(TRUE){
    U <- runif(1)
    t <- -log(U)/Lambdas[State]
    T <- T - t
    if(T < 0){
      t <- T + t
    }
    
    Z    <- rnorm(1)
    SNew <- SNew*exp(((r-(Sigmas[State]*Sigmas[State]/2))*t)+(Sigmas[State]*sqrt(t)*Z))
    YT1  <- YT1 + ((r-(Sigmas[1]*Sigmas[1]/2))*t)+(Sigmas[1]*sqrt(t)*Z)
    YT2  <- YT2 + ((r-(Sigmas[2]*Sigmas[2]/2))*t)+(Sigmas[2]*sqrt(t)*Z)
    
    if(State == 1){
      State <- State + 1
    }else{
      State <- State - 1
    }
    
    if(T < 0){
      break
    }
  }
  
  ArrST    <- c(ArrST,SNew)
  ArrCall1 <- c(ArrCall1, exp(-r*Time)*max(0,SNew-Strikes[1]))
  ArrCall2 <- c(ArrCall2, exp(-r*Time)*max(0,SNew-Strikes[2]))
  ArrCall3 <- c(ArrCall3, exp(-r*Time)*max(0,SNew-Strikes[3]))
  ArrCalls <- matrix(c(ArrCall1,ArrCall2,ArrCall3),nrow = 3,ncol=Paths,byrow = TRUE)
  
  ArrYT1Call1 <- c(ArrYT1Call1, exp(-r*Time)*max(0,S0*exp(YT1)-Strikes[1]))
  ArrYT1Call2 <- c(ArrYT1Call2, exp(-r*Time)*max(0,S0*exp(YT1)-Strikes[2]))
  ArrYT1Call3 <- c(ArrYT1Call3, exp(-r*Time)*max(0,S0*exp(YT1)-Strikes[3]))
  ArrYT1Call  <- matrix(c(ArrYT1Call1,ArrYT1Call2,ArrYT1Call3),nrow = 3,ncol=Paths,byrow = TRUE)
  ArrYT2Call1 <- c(ArrYT2Call1, exp(-r*Time)*max(0,S0*exp(YT2)-Strikes[1]))
  ArrYT2Call2 <- c(ArrYT2Call2, exp(-r*Time)*max(0,S0*exp(YT2)-Strikes[2]))
  ArrYT2Call3 <- c(ArrYT2Call3, exp(-r*Time)*max(0,S0*exp(YT2)-Strikes[3]))
  ArrYT2Call  <- matrix(c(ArrYT2Call1,ArrYT2Call2,ArrYT2Call3),nrow = 3,ncol=Paths,byrow = TRUE)

  ArrYT1   <- c(ArrYT1,YT1)
  ArrYT2   <- c(ArrYT2,YT2)
}

ArrYT1Sq <- ArrYT1^2
ArrYT2Sq <- ArrYT2^2
ArrYT1Exp <- exp(ArrYT1)
ArrYT2Exp <- exp(ArrYT2)

StatsST <- c(mean(ArrCall1),mean(ArrCall2),mean(ArrCall3),var(ArrCall1),var(ArrCall2),var(ArrCall3))

StatsYT1      <- SingleControl(ArrYT1,MeanYT1,ArrCalls,Strikes)
StatsYT2      <- SingleControl(ArrYT2,MeanYT2,ArrCalls,Strikes)
StatsYT1Sq    <- SingleControl(ArrYT1Sq,MeanYT1Sq,ArrCalls,Strikes)
StatsYT2Sq    <- SingleControl(ArrYT2Sq,MeanYT2Sq,ArrCalls,Strikes)
StatsYT1Exp   <- SingleControl(ArrYT1Exp,MeanYT1Exp,ArrCalls,Strikes)
StatsYT2Exp   <- SingleControl(ArrYT2Exp,MeanYT2Exp,ArrCalls,Strikes)
StatsYT1Call  <- SingleControlCall(ArrYT1Call,MeanYT1Call,ArrCalls,Strikes)
StatsYT2Call  <- SingleControlCall(ArrYT2Call,MeanYT2Call,ArrCalls,Strikes)

# Multi control var for exp(YT1) and exp(YT2)
ArrMultiVar <- vector()
ArrMultiMean <- vector()
for(i in c(1:length(Strikes))){
  coeffs <- summary(lm(ArrCalls[i,]~ArrYT1Exp + ArrYT2Exp))$coefficients[,1]
  CoeffYT1Exp <- coeffs[2][["ArrYT1Exp"]]
  CoeffYT2Exp <- coeffs[3][["ArrYT2Exp"]]
  ArrNew <- ArrCalls[i,]
  ArrNew <- ArrNew - CoeffYT1Exp*(ArrYT1Exp - MeanYT1Exp)
  ArrNew <- ArrNew - CoeffYT2Exp*(ArrYT2Exp - MeanYT2Exp)
  ArrMultiVar  <- c(ArrMultiVar,var(ArrNew))
  ArrMultiMean <- c(ArrMultiMean,mean(ArrNew))
}
StatsMulti1 <- c(ArrMultiMean,ArrMultiVar)

# Multi control var for YT1^2 and YT2^2
ArrMultiVar <- vector()
ArrMultiMean <- vector()
for(i in c(1:length(Strikes))){
  coeffs <- summary(lm(ArrCalls[i,]~ArrYT1Sq + ArrYT2Sq))$coefficients[,1]
  CoeffYT1Sq <- coeffs[2][["ArrYT1Sq"]]
  CoeffYT2Sq <- coeffs[3][["ArrYT2Sq"]]
  ArrNew <- ArrCalls[i,]
  ArrNew <- ArrNew - CoeffYT1Sq*(ArrYT1Sq - MeanYT1Sq)
  ArrNew <- ArrNew - CoeffYT2Sq*(ArrYT2Sq - MeanYT2Sq)
  ArrMultiVar  <- c(ArrMultiVar,var(ArrNew))
  ArrMultiMean <- c(ArrMultiMean,mean(ArrNew))
}
StatsMulti2 <- c(ArrMultiMean,ArrMultiVar)

# Multi control var for Call option using On S0*exp(YT1) and S0*exp(YT2)
ArrMultiVar <- vector()
ArrMultiMean <- vector()
for(i in c(1:length(Strikes))){
  coeffs <- summary(lm(ArrCalls[i,]~ArrYT1Call[i,] + ArrYT2Call[i,]))$coefficients[,1]
  CoeffYT1Call <- coeffs[2][["ArrYT1Call[i, ]"]]
  CoeffYT2Call <- coeffs[3][["ArrYT2Call[i, ]"]]
  ArrNew <- ArrCalls[i,]
  ArrNew <- ArrNew - CoeffYT1Call*(ArrYT1Call[i,] - MeanYT1Call[i])
  ArrNew <- ArrNew - CoeffYT2Call*(ArrYT2Call[i,] - MeanYT2Call[i])
  ArrMultiVar  <- c(ArrMultiVar,var(ArrNew))
  ArrMultiMean <- c(ArrMultiMean,mean(ArrNew))
}
StatsMulti3 <- c(ArrMultiMean,ArrMultiVar)

# Multi control var for all the control variates
ArrMultiVar <- vector()
ArrMultiMean <- vector()
for(i in c(1:length(Strikes))){
  coeffs <- summary(lm(ArrCalls[i,]~ArrYT1 + ArrYT1Sq + ArrYT1Exp + ArrYT2Exp + ArrYT1Call[i,] + ArrYT2Call[i,]))$coefficients[,1]
  CoeffYT1     <- coeffs[2][["ArrYT1"]]
  CoeffYT1Sq   <- coeffs[3][["ArrYT1Sq"]]
  CoeffYT1Exp  <- coeffs[4][["ArrYT1Exp"]]
  CoeffYT2Exp  <- coeffs[5][["ArrYT2Exp"]]
  CoeffYT1Call <- coeffs[6][["ArrYT1Call[i, ]"]]
  CoeffYT2Call <- coeffs[7][["ArrYT2Call[i, ]"]]
  ArrNew <- ArrCalls[i,]
  ArrNew <- ArrNew - CoeffYT1*(ArrYT1 - MeanYT1)
  ArrNew <- ArrNew - CoeffYT1Sq*(ArrYT1Sq - MeanYT1Sq)
  ArrNew <- ArrNew - CoeffYT1Exp*(ArrYT1Exp - MeanYT1Exp)
  ArrNew <- ArrNew - CoeffYT2Exp*(ArrYT2Exp - MeanYT2Exp)
  ArrNew <- ArrNew - CoeffYT1Call*(ArrYT1Call[i,] - MeanYT1Call[i])
  ArrNew <- ArrNew - CoeffYT2Call*(ArrYT2Call[i,] - MeanYT2Call[i])
  ArrMultiVar  <- c(ArrMultiVar,var(ArrNew))
  ArrMultiMean <- c(ArrMultiMean,mean(ArrNew))
}
StatsMulti4 <- c(ArrMultiMean,ArrMultiVar)

print(StatsST)
print(StatsYT1)
print(StatsYT2)
print(StatsYT1Sq)
print(StatsYT2Sq)
print(StatsYT1Exp)
print(StatsYT2Exp)
print(StatsYT1Call)
print(StatsYT2Call)
print(StatsMulti1)
print(StatsMulti2)
print(StatsMulti3)
print(StatsMulti4)
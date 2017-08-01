# Part a - Pilot Simulation

Paths <- 1000
m     <- 11
r     <- 0.05
sigma <- 0.25
S0    <- 100
Time  <- 1
KArr  <- c(90,100,110,120)
Z     <- 2.58

CovYZ <- vector('numeric')
VarZ  <- vector('numeric')
x     <- c(1:m)
ZAvg  <- S0*mean(exp(r*x*T/m))

for(K in KArr){
  PayArr <- vector('numeric')
  ZArr   <- vector('numeric')
  
  for(i in 1:Paths){
    SArr   <- vector('numeric')
    SOld  <- S0
    
    for(i in 1:m){
      SNew <- SOld*exp(((r-(sigma*sigma/2))/m)+(sigma*sqrt(1/m)*rnorm(1)))
      SArr <- c(SArr,SNew)
      SOld <- SNew
    }
    
    Payoff <- (exp(-r*Time))*max(0,mean(SArr)-K)  
    PayArr <- c(PayArr,Payoff)
    ZArr   <- c(ZArr,mean(SArr))
  }
  
  CovYZ <- c(CovYZ,cov(PayArr,ZArr))
  VarZ  <- c(VarZ,var(ZArr))
}

# Adding Control Variates

Paths <- 10000
Loop  <- c(1:4)

NaivePrice  <- vector('numeric')
NaiveVar    <- vector('numeric')
NaiveLB     <- vector('numeric')
NaiveUB     <- vector('numeric')

NewPrice  <- vector('numeric')
NewVar    <- vector('numeric')
NewLB     <- vector('numeric')
NewUB     <- vector('numeric')

for(y in Loop){
  K      <- KArr[y]
  PayArr <- vector('numeric')
  ZArr   <- vector('numeric')
  
  for(i in 1:Paths){
    SArr   <- vector('numeric')
    SOld  <- S0
    
    for(i in 1:m){
      SNew <- SOld*exp(((r-(sigma*sigma/2))/m)+(sigma*sqrt(1/m)*rnorm(1)))
      SArr <- c(SArr,SNew)
      SOld <- SNew
    }
    
    Payoff <- (exp(-r*Time))*max(0,mean(SArr)-K) 
    PayArr <- c(PayArr,Payoff)
    ZArr   <- c(ZArr,mean(SArr))
  }

  NewArray  <- PayArr + (-1*CovYZ[y]*(ZArr - ZAvg)/VarZ[y])
  
  NaivePrice  <- c(NaivePrice,mean(PayArr))
  NaiveVar    <- c(NaiveVar,var(PayArr))
  NaiveLB     <- c(NaiveLB, mean(PayArr) - Z*sqrt((var(PayArr))/Paths))
  NaiveUB     <- c(NaiveUB, mean(PayArr) + Z*sqrt((var(PayArr))/Paths))
  
  NewPrice  <- c(NewPrice,mean(NewArray))
  NewVar    <- c(NewVar,var(NewArray))
  NewLB     <- c(NewLB, mean(NewArray) - Z*sqrt((var(NewArray))/Paths))
  NewUB     <- c(NewUB, mean(NewArray) + Z*sqrt((var(NewArray))/Paths))
}

print(NaivePrice)
print(NaiveVar)
print(NaiveLB)
print(NaiveUB)

print(NewPrice)
print(NewVar)
print(NewLB)
print(NewUB)
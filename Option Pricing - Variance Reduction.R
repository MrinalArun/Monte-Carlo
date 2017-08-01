library(OptionPricing)

# Pilot Simulation for Control Variate

Paths <- 1000
r     <- 0.05
sigma <- 0.25
S0    <- 100
Time  <- 1
KArr  <- c(110,120)
Z     <- 1.96
L     <- 105

PayArr1 <- vector()
PayArr2 <- vector()
PayArr3 <- vector()
PayArr4 <- vector()

Mean   <- S0*exp(r*Time/2)
Prices <- vector()

for(i in 1:Paths){
  SOld  <- S0
  Rand <- rnorm(1)
  SMid <- SOld*exp(((r-(sigma*sigma/2))*Time/2)+(sigma*sqrt(Time/2)*Rand))

  # Naive calculation of expected Payoff
  if(SMid>=L){
    Payoff2 <- BS_EC( T = Time/2, K = KArr[2], r = r, sigma = sigma, S0 = SMid)[1][["price"]]
  }else{
    Payoff2 <- BS_EC( T = Time/2, K = KArr[1], r = r, sigma = sigma, S0 = SMid)[1][["price"]]
  }  

  PayArr2 <- c(PayArr2,Payoff2)
  Prices  <- c(Prices,SMid) #Using the stock price at T/2 as the control variate
}

Coeff <- -cov(Prices,PayArr2)/var(Prices)

# Actual Simualtion
Paths <- 10000

for(i in 1:Paths){
  SArr   <- vector('numeric')
  SOld  <- S0

  Rand <- rnorm(1)
  SMid <- SOld*exp(((r-(sigma*sigma/2))*Time/2)+(sigma*sqrt(Time/2)*Rand))
  SAnt <- SOld*exp(((r-(sigma*sigma/2))*Time/2)+(sigma*sqrt(Time/2)*-Rand))
  SEnd <- SMid*exp(((r-(sigma*sigma/2))*Time/2)+(sigma*sqrt(Time/2)*rnorm(1)))

  if(SMid>=L){
    Payoff1 <- (exp(-r*Time))*max(0,SEnd-KArr[2])
    Payoff2 <- BS_EC( T = Time/2, K = KArr[2], r = r, sigma = sigma, S0 = SMid)[1][["price"]]
  }else{
    Payoff1 <- (exp(-r*Time))*max(0,SEnd-KArr[1])
    Payoff2 <- BS_EC( T = Time/2, K = KArr[1], r = r, sigma = sigma, S0 = SMid)[1][["price"]]
  }

  if(SAnt>=L){
    Payoff3 <- BS_EC( T = Time/2, K = KArr[2], r = r, sigma = sigma, S0 = SAnt)[1][["price"]]
  }else{
    Payoff3 <- BS_EC( T = Time/2, K = KArr[1], r = r, sigma = sigma, S0 = SAnt)[1][["price"]]
  }

  PayArr1 <- c(PayArr1,Payoff1)
  PayArr2 <- c(PayArr2,Payoff2)
  PayArr3 <- c(PayArr3,mean(c(Payoff2,Payoff3)))
  Payoff4 <- Payoff2 + Coeff*(SMid-Mean)
  PayArr4 <- c(PayArr4,Payoff4)
}

Prices   <- c(mean(PayArr1),mean(PayArr2),mean(PayArr3),mean(PayArr4))
Variance <- c(var(PayArr1),var(PayArr2),var(PayArr3),var(PayArr4))
LowerBnd <- Prices - Z*sqrt(Variance/Paths)
UpperBnd <- Prices + Z*sqrt(Variance/Paths)

print(Prices)
print(Variance)
print(LowerBnd)
print(UpperBnd)
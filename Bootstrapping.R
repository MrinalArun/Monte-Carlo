library(MASS)

Bostondf <- Boston

# Part a - Estimating population mean of medv data
muhat <- mean(Bostondf$medv) # Sample mean is used as the estimate for population mean
muhat

# Part b - Estimating the standard error in the estimator above
sighat <- sd(Bostondf$medv)/sqrt(nrow(Bostondf))
sighat

# Part c- Estimating the standard error using bootstrap
B       = 200
Sum     = 0

for(i in 1:B){
  set.seed(B)
  BootData <- sample(Bostondf$medv,size = nrow(Bostondf),replace = TRUE)
  Sum      <- Sum + (mean(BootData)-muhat)^2
  MeanSum  <- MeanSum + mean(BootData)
}
sigboot <- Sum/B
sigboot

# Part d - providing a 95% confidence interval
lowlmtboot <- muhat - 2*sigboot
uplmtboot  <- muhat + 2*sigboot
c(lowlmtboot,uplmtboot) # CI from bootstrap

t.test(Bostondf$medv) # CI from the given sample

# Part e - Estimate of the median of medv data
medhat <- median(Bostondf$medv)
medhat

# Part f - Estimate the standard error of medhat
Sum = 0
for(i in 1:B){
  set.seed(B)
  BootData <- sample(Bostondf$medv,size = nrow(Bostondf),replace = TRUE)
  Sum      <- Sum + (median(BootData)-medhat)^2
}
medsigma <- Sum/B
medsigma

# Part g - Estimate of the quantile of medv data
quanthat <- quantile(Bostondf$medv,0.1)
quanthat

# Part h - Estimate the standard error of quanthat
Sum = 0
for(i in 1:B){
  set.seed(B)
  BootData <- sample(Bostondf$medv,size = nrow(Bostondf),replace = TRUE)
  Sum      <- Sum + (quantile(BootData,0.1)-quanthat)^2
}
quantsigma <- Sum/B
quantsigma


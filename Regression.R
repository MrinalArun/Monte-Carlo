library(boot)

# Part a - Generate simulated data

set.seed(1)
y=rnorm(100)
x=rnorm(100)
y=x-2*x^2+rnorm(100)

# Part b - Scatter plot of x vs y
plot(x,y)

# Part c & d - Computing LOOCV for different models

seeds = c(13,50)

for(seed in seeds){
  set.seed(seed)
  x1=rnorm(100)
  y=x1-2*x1^2+rnorm(100)
  x2 = x1^2
  x3 = x1^3
  x4 = x1^4
  
  table  = cbind(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4)
  df     = data.frame(table)
  errors = numeric()

  glm.fit1 = glm(y~x1,data=df)  # Linear fit
  errors   = c(errors, cv.glm(df,glm.fit1)$delta[1])   # LOOCV Errors
  
  glm.fit2 = glm(y~x1+x2,data=df)  # Quadratic fit
  errors   = c(errors, cv.glm(df,glm.fit2)$delta[1])  
  
  glm.fit3 = glm(y~x1+x2+x3,data=df)  # Cubic fit
  errors   = c(errors, cv.glm(df,glm.fit3)$delta[1])  
  
  glm.fit4 = glm(y~x1+x2+x3+x4,data=df)  # Quartic fit
  errors   = c(errors, cv.glm(df,glm.fit4)$delta[1])  
  
  print(seed)
  print(errors)  
}

summary(glm.fit3)
summary(glm.fit4)
library(glmnet)
# Part a - Reading the data into a data frame
setwd("G:/Studies/Columbia/Sem 2 - Spring 2017/IEOR E4525 - Machine Learning for OR & FE/Assignments/Assignment 3")
Df <-read.table("ceo.csv", header = TRUE,sep = ",")

# Part b - Dropping the last column - X
Df$X <- NULL

# Part d - Standardizing the independent columns
NewDf <- Df
NewDf$salary  <- scale(Df$salary)[,1]
NewDf$tenure  <- scale(Df$tenure)[,1]
NewDf$age     <- scale(Df$age)[,1]
NewDf$sales   <- scale(Df$sales)[,1]
NewDf$profits <- scale(Df$profits)[,1]
NewDf$assets  <- scale(Df$assets)[,1]

# Part c - Splitting the data into training and test
set.seed(1)

X <- model.matrix(totcomp~.,NewDf)[,-1]
Y <- NewDf$totcomp

rows      <- as.numeric(nrow(X)) # No of rows
shuffled <- sample(1:rows) 
trainrows <- shuffled[1:round(0.75*rows)] # 75% of rows for training data
testrows  <- shuffled[(round(0.75*rows)+1):rows]
Xtrain    <- X[trainrows,]
Xtest     <- X[testrows,] # Rest of the 25% for test data
Ytrain    <- Y[trainrows]
Ytest     <- Y[testrows] 

LY         <- log(Y)
LYtrain    <- LY[trainrows]
LYtest     <- LY[testrows]

# Part e - Best LASSO Model for predicting TotComp
grid    = 10^seq(10,-2,length = 121)
Model1  = glmnet(Xtrain,Ytrain,alpha=1,lambda = grid)
CVOut1  = cv.glmnet(Xtrain,Ytrain,alpha=1,lambda = grid)
Lambda1 = CVOut1$lambda.min
plot(CVOut1)
print(Lambda1)

# Part f - Best LASSO Model for predicting log(TotComp)

Model2  = glmnet(Xtrain,LYtrain,alpha=1,lambda = grid)
CVOut2  = cv.glmnet(Xtrain,LYtrain,alpha=1,lambda = grid)
Lambda2 = CVOut2$lambda.min
print(Lambda2)
plot(CVOut2)
print(Lambda2)

# Part g - Predict the test data using both the models
Pred1 = predict(Model1,s=Lambda1,newx=Xtest)
MSE1  = mean((Pred1-Ytest)^2)
print(MSE1)

Pred2 = predict(Model2,s=Lambda2,newx=Xtest)
MSE2  = mean((Pred2-LYtest)^2)
print(MSE2)
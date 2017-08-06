library(e1071)
library(caret)

# Part a - Estimating the generalization error
spamdata <- read.csv("G:/Studies/Columbia/Sem 2 - Spring 2017/IEOR E4525 - Machine Learning for OR & FE/Assignments/Assignment 2/spam.csv")
names(spamdata)

# Set columns timeofday, special, isuid and id to categorical variables
spamdata$time.of.day <- as.factor(spamdata$time.of.day)
spamdata$special     <- as.factor(spamdata$special)
spamdata$isuid       <- as.factor(spamdata$isuid)
spamdata$id          <- as.factor(spamdata$id)

# Add a column which contains a binary variable which checks if spampct is na or not
spamdata$spampctna   <- as.factor(is.na(spamdata$spampct))

rows      <- as.numeric(nrow(spamdata)) # No of rows
trainrows <- round(0.8*rows) # 80% of rows for training data
train <- spamdata[1:trainrows,]
test  <- spamdata[(trainrows+1):rows,] # Rest of the 20% for test data

model  <- naiveBayes(spam ~ ., data = train) # Fit the naive bayes
predict <- predict(model, test) # Predict on the test data
table  <- table(predict, test$spam) # Build the confusion matrix
print('The confusion matrix of the model is')
confusionMatrix(table)
print('The prediction error in the test data is')
error  <- sum(table[row(table)!=col(table)])/sum(table)
error

# Part c - Reducing the number of cases when the model predicts spam, when it actually isn't

tablenew <- table(predict(model, test, type="raw")[,"yes"]>=0.995,test$spam == "yes")
print('The old confusion matrix is')
table     # Old table
print('The new confusion matrix is')
tablenew  # New table

# Part b - Estimating the test error by repeating multiple times

K = 10
errors   <- numeric()

for(x in 1:K){ # Repeat the experiment in part a 10 times
  set.seed(x)
  array    <- c(1:rows)
  shuffled <- sample(array) # shuffle the array c(1:2171)

  trainnew <- spamdata[shuffled[1:trainrows],] # Pick the first 80% of the data in the shuffled array for training
  testnew  <- spamdata[shuffled[(trainrows+1):rows],] # Pick the remaining 20% for testing
 
  model  <- naiveBayes(spam ~ ., data = trainnew) # Fit the naive bayes
  predict <- predict(model, testnew) # Predict on the test data
  table  <- table(predict, testnew$spam) # Build the confusion matrix
  error  <- sum(table[row(table)!=col(table)])/sum(table)
  errors <- c(errors, error) # Store the errors in an array
}
print('The summary of the error % across the 10 runs is')
summary(errors)
print('The variance of the errors is')
var(errors)
print('The standard deviation of the errors is')
sd(errors)
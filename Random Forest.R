library(caret)
library(gbm) 
library(MASS)
library(tree)

MedianCrime <- median(Boston$crim) # Median crime rate
Y <- Boston$crim >= MedianCrime # An array indicating if the crime is greater than median or not 

NewDf <- Boston
NewDf$crim    <- NULL  #Dropping the crim column
NewDf$Y       <- Y     #Adding the indicator vector to the Data Frame

# Split the data into training and test
Test  <- NewDf[401:nrow(NewDf),]
Train <- NewDf[1:400,]

# # Part a - Logistic Regression
# LogFit  <- glm(Y~.,family = binomial,data = Train)
# LogPred <- predict(LogFit,newdata = Test,type= "response")
# LogPred <- LogPred >= 0.5
# 
# #Get the confusion matrix for training and test data
# confusionMatrix(table(predict(LogFit,newdata = Train,type= "response")>=0.5,Train$Y))
# confusionMatrix(table(predict(LogFit,newdata = Test,type= "response")>=0.5,Test$Y))
# confusionMatrix(table(LogPred,Test$Y)) 
# # Accuracy of 89% for training and 90.57% for test data
# 
# # Part b - Classificaton trees
# TreeFit <- tree(as.factor(Y)~.,Train)
# confusionMatrix(table(predict(TreeFit,Train,type = "class"),Train$Y))
# confusionMatrix(table(predict(TreeFit,Test,type = "class"),Test$Y))
# # Accuracy of 96.25% for training and 70.75% for test data
# 
# #Prune the tree to improve the performance
# TreeCV <- cv.tree(TreeFit,FUN = prune.misclass)
# par(mfrow=c(1,2))
# plot(TreeCV$size,TreeCV$dev ,type="b")
# plot(TreeCV$k,TreeCV$dev ,type="b")
# 
# #Size 6 produces the best results
# PruneTree = prune.misclass(TreeFit,best=6)
# # par(mfrow=c(1,1))
# # plot(PruneTree)
# # text(PruneTree,pretty=0)
# 
# confusionMatrix(table(predict(PruneTree,Train,type = "class"),Train$Y))
# confusionMatrix(table(predict(PruneTree,Test,type = "class"),Test$Y))
# # Accuracy of 96% for training and 70% for test data

# Part c - Bagging

#Split the training data into training and validation
Sample <- sample(1:400)
TrainNew <- Train[Sample[1:320],]
Valid    <- Train[Sample[321:400],]
IntDepth  <- c(1:10)
LeastError <- 1
BestVal  <- -1

# Doing a cross validation to determine the Max number of nodes
for(i in IntDepth)
{
  BoostFit = gbm(as.factor(Y)~.,data = TrainNew, distribution = "multinomial", n.trees = 5000, interaction.depth = i)
  # confusionMatrix(table((predict(BoostFit,Train,n.trees = 5000)<=0)[1:400],Train$Y))
  # confusionMatrix(table((predict(BoostFit,Test,n.trees = 5000)<=0)[1:106],Test$Y))
  # # Accuracy of 99.5% for training and 87.74% for test data
  Error  = sum((predict(BoostFit,Valid,n.trees = 5000)<=0)[1:80]!=Valid$Y)/80
  if(Error < LeastError){
    LeastError <- Error
    BestVal  <- i
  }
  print(Error)
}


# # Part d - Random Forest
# RFFit = randomForest(as.factor(Y)~., data = Train, mtry = 4, importance = TRUE)
# confusionMatrix(table(predict(RFFit,Train,type = "class"),Train$Y))
# confusionMatrix(table(predict(RFFit,Test,type = "class"),Test$Y))
# importance(BagFit)
# # Accuracy of 100% for training and 91.51% for test data
# 
# Part e - Boosting
BoostFit = gbm(as.factor(Y)~.,data = Train, distribution = "multinomial", n.trees = 5000, interaction.depth = IntDepth)
par(mfrow=c(1,1))
summary(BoostFit)
confusionMatrix(table((predict(BoostFit,Train,n.trees = 5000)<=0)[1:400],Train$Y))
confusionMatrix(table((predict(BoostFit,Test,n.trees = 5000)<=0)[1:106],Test$Y))
# Accuracy of 99.5% for training and 87.74% for test data

# setwd("G:/Studies/Columbia/Sem 2 - Spring 2017/IEOR E4525 - Machine Learning for OR & FE/Assignments/Assignment 7")
# Train = read.csv("Train.csv", header = FALSE)
# Test  = read.csv("Test.csv", header = FALSE)
# 
# # Part a - Average rating as baseline estimator
# Ratings = Train$V3
# Average = mean(Ratings)
# 
# Pred = rep(Average,nrow(Test))
# RMSE = sqrt(mean((Pred-Test$V3)^2))
# print(RMSE)

# # Part b - Adding Biases
# MovieBias = rep(0,1682)
# UserBias  = rep(0,943)
# 
# for(i in 1:length(MovieBias)){
#   if(sum(Train$V2 == i) == 0){
#     MovieBias[i] = 0
#   }else{
#     MovieBias[i] = mean(Train$V3[Train$V2==i]) - Average
#   }
# }
# 
# for(i in 1:length(UserBias)){
#   if(sum(Train$V1 == i) == 0){
#     UserBias[i] = 0
#   }else{
#     UserBias[i] = mean(Train$V3[Train$V1==i]) - Average
#   }
# }
# 
# Pred = rep(Average,nrow(Test))
# for(i in 1:length(Pred)){
#   Pred[i] = Pred[i] + MovieBias[Test$V2[i]] + UserBias[Test$V1[i]]
# }
# 
# RMSE = sqrt(mean((Pred-Test$V3)^2))
# print(RMSE)

# # Part c - Regularization
# # LambdaArr = c(3:8)
# LambdaArr = c(4)
# RMSEArr   = vector()
# for(Lambda in LambdaArr){
#   MovieBias = rep(0,1682)
#   UserBias  = rep(0,943)
#   for(j in 1:20){
#     for(i in 1:length(MovieBias)){
#       Count = sum(Train$V2==i)
#       MovieBias[i] = (sum(Train$V3[Train$V2==i]) - sum(UserBias[Train$V1[Train$V2==i]]) - Count*Average)/(Count + Lambda)
#     }
# 
#     for(i in 1:length(UserBias)){
#       Count = sum(Train$V1==i)
#       UserBias[i] = (sum(Train$V3[Train$V1==i]) - sum(MovieBias[Train$V2[Train$V1==i]]) - Count*Average)/(Count + Lambda)
#     }
#   }
# 
#   Pred = rep(Average,nrow(Test))
#   for(i in 1:length(Pred)){
#     Pred[i] = Pred[i] + MovieBias[Test$V2[i]] + UserBias[Test$V1[i]]
#   }
# 
#   RMSE = sqrt(mean((Pred-Test$V3)^2))
#   print(c(Lambda,RMSE))
#   RMSEArr = c(RMSEArr,RMSE)
# }

# # Part d - Construct the residual matrix
# 
# Res = matrix(0,nrow = 1682,943)
# Test$Pred = Pred
# TrainPred = rep(Average,nrow(Train))
# for(i in 1:length(TrainPred)){
#   TrainPred[i] = TrainPred[i] + MovieBias[Train$V2[i]] + UserBias[Train$V1[i]]
# }
# 
# Train$Pred = TrainPred
# for(i in 1:nrow(Train)){
#   Res[Train$V2[i],Train$V1[i]] = Train$V3[i]-Train$Pred[i]
# }


# Part e - Use neighbourhood method
DMat = matrix(0,nrow = 1682, ncol = 1682)
for(i in 1:1682){
  if(sum(Res[i,]^2)==0){
    DMat[i,] = 0
  }else{
    for(j in 1:1682){
      if(sum(Res[j,]^2)==0){
        DMat[i,j] = 0
      }else{
        DMat[i,j] = sum(Res[i,]*Res[j,])/(sqrt(sum(Res[i,]^2)*sum(Res[j,]^2)))
      }
    }
  }
}


# Get the neighbours for each movie
RMSEVect    = vector()
Neighbours = c(1:50)
for(NeighCount in Neighbours){
  NeighMat = matrix(0,nrow = 1682, ncol = NeighCount)
  for(Row in 1:1682){
    kthMax     = sort(abs(DMat[Row,]),partial = 1682 - NeighCount)[1682-NeighCount]
    Indices    = c(1:1682)[abs(DMat[Row,]) > kthMax]
    Indices    = c(Indices,c(1:1682)[abs(DMat[Row,]) == kthMax][1:(NeighCount + 1 - length(Indices))])
    Indices    = Indices[Indices != Row]
    if(length(Indices)!=NeighCount){
      Indices = Indices[1:NeighCount]
    }
    NeighMat[Row,] = Indices
  }


  TestPredNew = vector()
  for(i in 1:nrow(Test)){
    User  = Test$V1[i]
    Movie = Test$V2[i]
    Pred = Average + MovieBias[Movie] + UserBias[User]
    Neig = NeighMat[Movie,]
    if(sum(abs(DMat[Movie,Neig])!=0)){
      Pred = Pred + sum(DMat[Movie,Neig]*Res[Neig,User])/sum(abs(DMat[Movie,Neig]))
    }
    TestPredNew = c(TestPredNew,Pred)
  }
  RMSE = sqrt(mean((TestPredNew-Test$V3)^2))
  RMSEVect = c(RMSEVect,RMSE)
}

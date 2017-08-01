# Part a

require(DAAG)
require(ggplot2)
require(MASS)
# This code implements reduced rank LDA (Fisher Discriminant Analysis)
# It can reproduce the subplots of Figure 4.8 in HTF by specifing coordinates a,b 
# For example, a=1,b=3 reproduces the top-left sub-figure of Figure 4.8

a=1 # First Fisher coordinate to plot
b=3 # second Fisher coordinate to plot

################################################################
# First download the training data from the HTF website
url<-"http://www-stat.stanford.edu/~tibs/ElemStatLearn/datasets/vowel.train"
vtrain<-read.table(url,header=TRUE,sep=',')
vtrain<-as.data.frame(vtrain)
# columns are row.names,y,x.1,x.2,x.3,x.4,x.5,x.6,x.7,x.8,x.9,x.10
# y is the class, and x.1 to x.10 are predictors

# Now download the test data
url<-"http://www-stat.stanford.edu/~tibs/ElemStatLearn/datasets/vowel.test"
vtest<-read.table(url,header=TRUE,sep=',')
vtest<-as.data.frame(vtest)
#################################################################

# Now find the Fisher discriminant directions using the "lda" function
ldatrain<-lda(y~x.1+x.2+x.3+x.4+x.5+x.6+x.7+x.8+x.9+x.10, data=vtrain)
scaling <- ldatrain$scaling

factors1 <- unname(scaling[,a])
factors2 <- unname(scaling[,b])

# Create a Data Frame with only the x's
NewDf <- vtrain
NewDf$row.names <- NULL
NewDf$y <- NULL
Y <- vtrain$y

# Get the coordinates along the 2 dimensions
CoordA <- vector()
CoordB <- vector()

for(i in 1:nrow(NewDf)){
  CoordA <- c(CoordA, sum(NewDf[i,]*factors1))
  CoordB <- c(CoordB, sum(NewDf[i,]*factors2))
}

#Scaled Coords
ScaledCoordA <- scale(CoordA)
ScaledCoordB <- scale(CoordB)

#Set colors for the scatter plot
colors <- c("black","blue2","chocolate4","blueviolet","cyan1","red","darkgray","bisque","darkgreen","chocolate1","chartreuse")
plot(ScaledCoordA[Y==1],ScaledCoordB[Y==1],col = colors[Y[1]], pch = 20, xaxt = 'n', yaxt = 'n', xlab = "Coordinate 1", ylab = "Coordinate 3", xlim = c(-2.5,2.5), ylim = c(-3.5,3.5))

# Plot the points
for(i in 2:11){
  points(ScaledCoordA[Y==i],ScaledCoordB[Y==i],col = colors[Y[i]], pch = 20)
}

# Get the means along the 2 dimensions
MeansA <- vector()
MeansB <- vector()

for(i in 1:11){
  MeansA <- c(MeansA, (sum(unname(ldatrain$means)[i,]*factors1)-mean(CoordA))/sd(CoordA))
  MeansB <- c(MeansB, (sum(unname(ldatrain$means)[i,]*factors2)-mean(CoordB))/sd(CoordB))
}

for(i in 1:11){
  points(MeansA[i],MeansB[i],col = colors[Y[i]], pch = "O")
}

# Part b
ArrTestError  <- vector()
ArrTrainError <- vector()

for(i in 1:10){
  ldanew     <- ldatrain
  
  # Choose just i directions
  ldanew$svd <- ldanew$svd[1:i]
  ldanew$scaling <- ldatrain$scaling[,1:i]
  
  trainpred <- predict(ldanew,vtrain)
  testpred  <- predict(ldanew,vtest)
  
  TestError  <- sum(testpred$class!=vtest$y)/462
  TrainError <- sum(trainpred$class!=vtrain$y)/528
  
  ArrTestError  <- c(ArrTestError,TestError)
  ArrTrainError <- c(ArrTrainError,TrainError)
}

plot(c(1:10),ArrTestError,ylim = c(0.25,0.75), xlab = "Dimension", ylab = "Misclassification Rate")
lines(c(1:10),ArrTestError)
lines(c(1:10),ArrTrainError, col = "blue")

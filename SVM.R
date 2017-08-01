require(mlbench)
require(e1071)

data(BreastCancer)
summary(BreastCancer)
na_idx<-rowSums(is.na(BreastCancer))>0
BreastCancer<-BreastCancer[!na_idx,]
L<-length(BreastCancer[,1])
tcon<-tune.control(sampling="cross",cross=5)

tvals<-(10:19)*0.05
N<-50
M<-length(tvals)

er<-matrix(0,nrow=N,ncol=M)

for (m in 1:M)
{
  for (n in 1:N)
  {
    train_idx<-sample(L,floor(L*tvals[m]))
    train_data<-BreastCancer[train_idx,2:11]
    test_data<-BreastCancer[-train_idx,2:11]
      
    bsvm <- best.svm(Class~.,data=train_data,kernel="radial", 
                     gamma = 2^(-2:2), cost = 2^(0:4),
                     tunecontrol=tcon)
    svm.pred<-predict(bsvm,test_data[,-10])
    er[n,m]<-sum(svm.pred!=test_data[,10])/length(test_data[,10])
  }
}  

er_sd<-apply(er,2,sd)
er_m<-apply(er,2,mean)

plot(tvals,er_m,ylim=c(0.017,0.035),xlab="t",ylab="error rate")
points(tvals,er_m-er_sd*1.96/sqrt(50),col=2)
points(tvals,er_m+er_sd*1.96/sqrt(50),col=3)

boxplot(er,names=tvals)
title(xlab="t",ylab="error rate")

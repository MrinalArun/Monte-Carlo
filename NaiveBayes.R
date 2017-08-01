# Q1 -- Naive Bayes

# check if we can predict missing spampct using observed variables

data<-read.csv(file.choose(),header=T) #choose spam2.csv
data_spam = data[which(data[,"spam"]==1),] # spam data
data_nspam = data[which(data[,"spam"]==0),] #not spam data
# subset of data_nspam with no entries missing
data_spam_nms = data_spam[which(is.finite(data_spam[,"spampct"])),] 
# subset of data_spam with no entries missing
data_nspam_nms = data_nspam[which(is.finite(data_nspam[,"spampct"])),]

# regression of spampct on other features for spam data
par(mfrow=c(1,2))
fit1 <- lm(spampct~., data=data_spam_nms)
summary(fit1) # very small R-square
# plotting the fit of spampct against precited spampct
plot(fitted(fit1),data_spam_nms[,"spampct"], main = 'Fitted vs True spampct for spam data', cex.main = 0.7) # not a good fit
#title('Fitted vs True spampct for spam data', cex = 0.5)

# regression of spampct on other features for not-spam data
fit2 <- lm(spampct~., data=data_nspam_nms)
summary(fit2) # very small R-square
# plotting the fit of spampct against precited spampct
plot(fitted(fit2),data_nspam_nms[,"spampct"], main = 'Fitted vs True spampct for spam data', cex.main = 0.7) # not a good fit

# Therefore, we just ignore the missing spampct value....see rest of the code in matlab
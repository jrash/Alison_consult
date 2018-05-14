#install and load any packages I used on your system
list.of.packages <- c("randomForest","tree","psych","ggplot2","gplots","stats","devtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

#clear workspace
rm(list=ls())

processed.dir <- 'C:/Users/Vestige/Google Drive/AlisonConsult/clinvitro/analyses/regressionTree/'
setwd(processed.dir)

d = read.csv("AllVarProcessReg.csv")[,-1]
head(d)
colnames(d)
ac50d = d[,c(-1,-(23:34))]
hist(d$AC50)


colnames(ac50d)
summary(ac50d)
head(ac50d)

##### RANDOMLY DIVIDE DATA INTO TRAINING (75%) & TEST (25%) SET #####
set.seed(123)
ind <- sample(2, nrow(ac50d),replace=T, prob=c(0.75,0.25))
train <- ac50d[ind==1,]
test <- ac50d[ind==2,]
attach(ac50d)

dim(train)
dim(ac50d)

#### Tree Based Methods
par(mfrow=c(1,1))
set.seed(1)
tree.d = tree(AC50~., data=ac50d)
summary(tree.d)
plot(tree.d)
text(tree.d, pretty=0)
cv.d = cv.tree(tree.d, K=10)
cv.d

plot(cv.d$size, cv.d$dev,type="b", cex=1.5, cex.lab=1.5, xlab="Tree Size", ylab="CV SSE")

prune.d=prune.tree(tree.d, best=9)
plot(prune.d)
text(prune.d,pretty =0)

#Test MSE
cv.d$dev[1]/10

#Average Deviation between predicted and actual
sqrt(cv.d$dev[1]/10)


#Test set R^2
tree.d = tree(AC50~., data=train)
yhat=predict(tree.d, newdata = test)
plot(yhat, test$AC50)
abline(0,1, col="blue", lty="dashed", lwd=2)
cor(yhat, test$AC50)^2



#remove variables with lots of missing data
#random forests cant handle missing data
ac50d2 = ac50d
summary(ac50d)
ac50d2$Grade = NULL
ac50d2$StageFinal = NULL
#remove missing data
ac50d2 = na.omit(ac50d2)

# #split into test and train set
# set.seed(1)
# ind <- sample(2, nrow(ac50d2),replace=T, prob=c(0.75,0.25))
# train <- ac50d2[ind==1,]
# test <- ac50d2[ind==2,]
# 
# 
# #Random Forests
# 
# dim(ac50d)
# 
# set.seed(1)
# p = c(1,2,3,9,13,18)
# 
# MSEs = data.frame(matrix(nrow=length(p),ncol=1000))
# 
# for(j in 1:length(p)){
#   rf.d=randomForest(x=train[,-19],y=train[,19],xtest= test[,-19], 
#                        ytest=test[,19], mtry=p[j],ntree=1000)
#   MSEs[j,] = rf.d$mse
# }
# 
# ?randomForest
# 
# #min Test MSE
# which(MSEs==min(MSEs), arr.ind=T)
# MSEs[1, 38]
# sqrt(MSEs[1,38])
# 
# cols = c("red","darkgreen", "black", "cyan", "purple", "brown")
# 
# plot(1, type="n", xlab="Number of Trees", ylab="Test MSE", xlim=c(0, 500), ylim=c(400,700))
# for(i in 1:length(p)){
#   lines(1:500,MSEs[i,1:500],type="l", lwd=2, col=cols[i])
# }
# legend("topright",legend=p,col=cols,lty=1,lwd=4, title="P",pt.cex=1.5)
# 
# rf.d=randomForest(AC50~.,data=ac50d2,mtry=1,ntree=1000)
# importance(rf.d)
# varImpPlot(rf.d, main="Variable Importance")
# 
# #test MSE for Random Forest
# rf.d=randomForest(AC50~.,data=train,mtry=3,ntree=960)
# yhat.rf = predict(rf.d,newdata=test)
# mean((yhat.rf-test$AC50)^2)
# cor(yhat.rf,test$AC50)

#Random Forests Take 2
set.seed(1)
p = c(1,2,3,9,13,18)
ntrees = seq(1,1000,10)

MSEs = data.frame(matrix(nrow=length(p),ncol=1000))

for(j in 1:length(p)){
  rf.d=randomForest(AC50~., data=ac50d2,mtry=p[j],ntree=1000)
  #extract OOB MSE
  MSEs[j,] = rf.d$mse
}

?randomForest

#min Test MSE
which(MSEs==min(MSEs), arr.ind=T)
MSEs[1, 38]

cols = c("red","darkgreen", "black", "cyan", "purple", "brown")

plot(1, type="n", xlab="Number of Trees", ylab="OOB MSE", xlim=c(0, 500), ylim=c(400,700))
for(i in 1:length(p)){
  lines(1:500,MSEs[i,1:500],type="l", lwd=2, col=cols[i])
}
legend("topright",legend=p,col=cols,lty=1,lwd=4, title="P",pt.cex=1.5)

rf.d=randomForest(AC50~.,data=ac50d2,mtry=1,ntree=300)
importance(rf.d)
varImpPlot(rf.d, main="Variable Importance")
sqrt(rf.d$mse[300])


#OOB R^2
names(rf.d)
plot(rf.d)
plot(rf.d$predicted, ac50d2$AC50)
abline(0,1, col="blue", lty="dashed", lwd=2)

cor(rf.d$predicted, ac50d2$AC50)
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
noClind = d[,-c(1,6,12:14,17,19)]
colnames(noClind)

summary(noClind)
noClind = na.omit(noClind)

summary(d$ResponseTaxane)
respon = d$ResponseTaxane

##### RANDOMLY DIVIDE DATA INTO TRAINING (75%) & TEST (25%) SET #####
set.seed(123)
ind <- sample(2, nrow(noClind),replace=T, prob=c(0.75,0.25))
train <- noClind[ind==1,]
test <- noClind[ind==2,]
attach(noClind)

dim(train)
dim(noClind)


####### Classification Tree
par(mfrow=c(1,1))
set.seed(1)
tree.d = tree(ResponseTaxane~., data=train)
summary(tree.d)
plot(tree.d)
text(tree.d, pretty=0)
cv.d = cv.tree(tree.d, K=93)
cv.d

?cv.tree

png("ClassTreePrune.png")
prune.d=prune.tree(tree.d, best=4)
plot(prune.d)
text(prune.d,pretty =0)
dev.off()

png("classTreeTune.png")
plot(cv.d$size, cv.d$dev,type="b", cex=1.5, cex.lab=1.5, xlab="Tree Size", ylab="CV SSE")
dev.off()


yhat=predict(tree.d, newdata = test, type="class")
?predict.tree
plot(yhat, ac50d$AC50)
abline(0,1, col="blue", lty="dashed", lwd=2)

#CV error rate
mean(yhat!=test$ResponseTaxane)



####### RANDOM FOREST FOR CLASSIFICATION #######
#### Random Forest for Classification


dim(noClind)


set.seed(4)
p = seq(1,26,3)
ntrees = seq(1,1000,10)

missclass_errors = data.frame(matrix(nrow=length(p),ncol=1000))

j=10
for(j in 1:length(p)){
  rf.d=randomForest(ResponseTaxane~., data=noClind,mtry=p[j],ntree=1000)
  missclass_errors[j,] = rf.d$err.rate[,1]
}

cols = c("red","darkgreen", "black", "cyan", "purple", "brown","orange","pink","darkgrey")


png("classRFtune.png")
plot(1, type="n", xlab="Number of Trees", ylab="OOB Misclass Error", xlim=c(0, 1000), ylim=c(.4,1))
for(i in 1:length(p)){
  lines(1:1000,missclass_errors[i,1:1000],type="l", lwd=2, col=cols[i])
}
legend("topright",legend=p,col=cols,lty=1,lwd=4, title="P",pt.cex=1.5)
dev.off()

#min misclassification error
which(missclass_errors==min(missclass_errors), arr.ind=T)
missclass_errors[8,1000]

set.seed(4)
rf.d=randomForest(ResponseTaxane~.,data=noClind,mtry=22,ntree=10000)
importance(rf.d)



png("classRFVarImport.png")
varImpPlot(rf.d, main="Variable Importance")
dev.off()

rf.d$err.rate[10000,1]


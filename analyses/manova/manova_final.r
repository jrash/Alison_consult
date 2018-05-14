list.of.packages <- c("arules","car", "psych")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

###############################################
###### clinvitro scripts: regression ##########
############################# by JJ and AAMR ##
rm(list=ls()) #################################


# Load general libraries
outdir <- 'C:/Users/Vestige/Google Drive/AlisonConsult/clinvitro/analyses/manova/'
setwd(outdir)
processed.dir <- 'C:/Users/Vestige/Google Drive/clinvitro/analyses/processed/clinical/'
curves.dir <- ("C:/Users/Vestige/Google Drive/clinvitro/analyses/curvefits.100/")


d = read.csv("C:/Users/Vestige/Google Drive/AlisonConsult/clinvitro/analyses/manova/AllVarProcessMan.csv")[,-1]
setwd("C:/Users/Vestige/Google Drive/AlisonConsult/clinvitro/analyses/manova/")

colnames(d)
results = d[,25:34]

#discretize continuous variables
d$Age <-discretize(d$Age,method="frequency")
d$TotalWeeks <-discretize(d$TotalWeeks,method="frequency")
d$NumCycles <-discretize(d$NumCycles,method="frequency")


d$AC50 <- NULL
d$EMAX <- NULL
d$mean_Score <- NULL

colnames(d)
head(d)

combo = d
colnames(combo)[22:31]

attach(combo)

#all dose response variables
resp1 = as.matrix(combo[,22:31])

library(caret)
findCorrelation(cor(resp1), cutoff = .7)

matmat <- data.matrix(cor(resp1) > .7)
matmat <- apply(matmat, 2, as.numeric)
gplots::heatmap.2(matmat, dendrogram = "none", Rowv = F, Colv = F, trace = "none")
as.numeric(cor(resp1) < .7)
cor.plot(resp1)

#less correlated dose response variables
resp2 = resp1[,c(1,4,6,10)]
cor(resp2)

#even less correlated dose response variables
resp3 = resp1[,c(3,6,10)]

pairs(resp1)
head(resp1)
cor(resp1)

pairs(resp2)
cor(resp2)

pairs(resp3)
cor(resp3)

summary(manova(resp1~ Smoke), test="Pillai")
summary(manova(resp2~ Smoke), test="Pillai")
summary(manova(resp3~ Smoke), test="Pillai")


save.image("manovaStart.R")
load("manovaStart.R")

vars = colnames(combo[,c(-1,-(22:31))])

#MANOVA for all dose response variables

colnames(combo[2:(length(vars)+1)])
resp1y <- vector("list", length(vars))
resp1p_ls <- vector("list", length(vars))
for (i in 2:(length(vars)+1)){
  print(colnames(combo)[i])
  lmfit <- manova(resp1 ~ combo[,i])
  resp1y[[i-1]] <- lmfit
  #pvalue for regressing each variable in df on resp1
  resp1p_ls[[i-1]] <- summary(manova(resp1 ~ combo[,i]))[[4]][1,6]
}
resp1ps = unlist(resp1p_ls)
p.adjust(resp1ps, method = "BH")
vars[p.adjust(resp1ps, method = "BH") < .2]

sigresp1 = cbind(vars[resp1ps<.05],resp1ps[resp1ps<.05])
sigresp1
resp1res = data.frame(variables = vars, pvalues = resp1ps)

#plot MANOVA results

pdf("Response1results.pdf",20,7)
plot(resp1res,cex.axis=.6)
abline(h=.05,lty="dashed", col="red")
dev.off()

#test the pairwise differences in means for significant effects
lmfit <- lm(resp1 ~ Smoke)
summary(Anova(lmfit))
anova(lmfit)
linearHypothesis(lmfit, "SmokeFormer Smoker = SmokeNever Smoked",verbose=T)
linearHypothesis(lmfit, "SmokeFormer Smoker",verbose=T)
linearHypothesis(lmfit, "SmokeNever Smoked",verbose=T)


#Boxplots of dose response vars

par(mfrow=c(3,4))
for(i in 1:10){
  boxplot(resp1[,i] ~ Smoke,main=colnames(resp1)[i], cex=.75)
}

# MANOVA for less correlated response variables
resp2y <- vector("list", dim(d)[2]-1)
resp2p_ls <- vector("list", dim(d)[2]-1)
i <- 3
for (i in 2:(length(vars)+1)){
  print(colnames(combo)[i])
  lmfit <- manova(resp2 ~ combo[,i])
  resp2y[[i-1]] <- lmfit
  #pvalue for regressing each variable in df on resp2
  resp2p_ls[[i-1]] <- summary(manova(resp2 ~ combo[,i]))[[4]][1,6]
}

# dont now why so many of the BH corrected p-values are the
# same, but I know the procedure Alex used in our paper was valid
resp2ps = unlist(resp2p_ls)
p.adjust(resp2ps, method = "BH")
vars[p.adjust(resp2ps, method = "BH") < .2]

sigresp2 = cbind(vars[resp2ps<.05],resp2ps[resp2ps<.05])
sigresp2
resp2res = data.frame(variables = vars, pvalues = resp2ps)

#plot MANOVA results

pdf("Response2results.pdf",20,7)
plot(resp2res,cex.axis=.6)
abline(h=.05,lty="dashed", col="red")
dev.off()

#test the pairwise differences in means for significant effects
lmfit <- lm(resp2 ~ Smoke)
summary(Anova(lmfit))
linearHypothesis(lmfit, "SmokeFormer Smoker = SmokeNever Smoked",verbose=T)
linearHypothesis(lmfit, "SmokeFormer Smoker",verbose=T)
linearHypothesis(lmfit, "SmokeNever Smoked",verbose=T)

lmfit <- lm(resp2 ~ Race)
summary(Anova(lmfit))
linearHypothesis(lmfit, "RaceWhite",verbose=T)

linearHypothesis(lmfit, "SmokeFormer Smoker",verbose=T)
linearHypothesis(lmfit, "SmokeNever Smoked",verbose=T)

#Boxplots of dose response vars

par(mfrow=c(1,4))
for(i in 1:4){
  boxplot(resp2[,i] ~ Smoke,main=colnames(resp2)[i],cex=.75)
} 

fs <- na.omit(resp2[Smoke == "Former Smoker", ])
ns <- na.omit(resp2[Smoke == "Never Smoked", ])
cs <- na.omit(resp2[Smoke == "Current Smoker", ])
apply(fs, 2, mean) - apply(ns, 2, mean)

fs_mean <- apply(fs, 2, mean)
ns_mean <- apply(ns, 2, mean)
cs_mean <- apply(cs, 2, mean)

plot(cs_mean, type = "l", col = "red")
lines(fs_mean, type = "l", col = "blue")
lines(ns_mean, type = "l", col = "black")


par(mfrow=c(1,4))
for(i in 1:4){
  boxplot(resp2[,i] ~ Race,main=colnames(resp2)[i],cex=.75)
}

white <- na.omit(resp2[Race == "White", ])
black <- na.omit(resp2[Race == "Black", ])
apply(white, 2, mean) - apply(black, 2, mean)

white_mean <- apply(white, 2, mean)
black_mean <- apply(black, 2, mean)

par(mfrow = c(1,1))
plot(white_mean, type = "l", col = "red")
lines(black_mean, type = "l", col = "blue")

mean(na.omit(resp2[Race == "White", ]))

#2way ANOVA with smoking and race

mod.resp <- lm(resp2 ~ Smoke*Race)
manova.resp = Anova(mod.resp)
manova.resp



#focus on X600
summary(aov(resp2[,2]~Smoke))
TukeyHSD(aov(resp2[,2]~Smoke))
summary(aov(resp2[,2]~Race))
TukeyHSD(aov(resp2[,2]~Race))


#Redoing John Jack's experiment


curves.dir <- ("C:/Users/Vestige/Google Drive/AlisonConsult/clinvitro/analyses/curvefits.100/")
fits <- read.table(file=paste(curves.dir,'ALL_LINES_Model_Output_File_ConfInt-Boot_PredInt-TRUE_Alpha-0.9_ConcTrans-FALSE_06May2015.csv',sep=''),sep=',',header=T,stringsAsFactors=FALSE)

hist(fits$AC50)
hist(fits$LEC)

preds = d[,c("ID",vars)]
head(preds)

# preds$AC50_Score <- NA
# preds$EMAX_Score <- NA
setdiff(preds$ID,fits$Cell_Line)
setdiff(fits$Cell_Line,preds$ID)
for (i in preds$ID){
  if (dim(fits[fits$Cell_Line==i,])[1]>0){
    preds$AC50_Score[preds$ID==i] <- fits$AC50[fits$Cell_Line==i]
    preds$EMAX_Score[preds$ID==i] <- fits$EMAX[fits$Cell_Line==i]
  }        
}
preds$mean_Score <- (preds$AC50_Score+preds$EMAX_Score)/2

# nocurves <- preds[-which(preds$AC50_Score%in%NA),]
nocurves <- preds

dim(nocurves)
vars = colnames(nocurves)[2:21]

AC50modelsy <- vector("list", 20)
AC50modelspvalsy <- vector("list", 20)
for (i in 2:(dim(nocurves)[2]-3)){
  lmfit <- lm(nocurves$AC50_Score~nocurves[,i])
  AC50modelsy[[i-1]] <- lmfit
  #pvalue for regressing each variable in df on AC50
  AC50modelspvalsy[[i-1]] <- summary(aov(lmfit))[[1]][1,5]
}

AC50ps = unlist(AC50modelspvalsy)
sigAC50 = cbind(vars[AC50ps<.05],AC50ps[AC50ps<.05])
sigAC50
AC50res = data.frame(variables = vars, pvalues = AC50ps)

#plot ANOVA results

pdf("AC50results.pdf",20,7)
plot(AC50res,cex.axis=.6)
abline(h=.05,lty="dashed", col="red")
dev.off()

EMAXmodelsy <- vector("list", 20)
EMAXmodelspvalsy <- vector("list", 20)
for (i in 2:(dim(nocurves)[2]-3)){
  lmfit <- lm(nocurves$EMAX_Score~nocurves[,i])
  EMAXmodelsy[[i-1]] <- lmfit
  EMAXmodelspvalsy[[i-1]] <- summary(aov(lmfit))[[1]][1,5]
}
EMAXps = unlist(EMAXmodelspvalsy)
sigEMAX = cbind(vars[EMAXps<.05],EMAXps[EMAXps<.05])
sigEMAX
EMAXres = data.frame(variables = vars, pvalues = EMAXps)

#plot ANOVA results

pdf("EMAXresults.pdf",20,7)
plot(EMAXres,cex.axis=.6)
abline(h=.05,lty="dashed", col="red")
dev.off()


par(mfrow=c(1,1))
boxplot(nocurves$EMAX_Score~nocurves$Smoke, ylab="EMAX score",cex=.75)

TukeyHSD(aov(nocurves$EMAX_Score~nocurves$Smoke))

Meanmodelsy <- vector("list", 20)
Meanmodelspvalsy <- vector("list", 20)
for (i in 2:(dim(nocurves)[2]-3)){
  lmfit <- lm(nocurves$mean_Score~nocurves[,i])
  Meanmodelsy[[i-1]] <- lmfit
  Meanmodelspvalsy[[i-1]] <- summary(aov(lmfit))[[1]][1,5]
}
Meanps = unlist(Meanmodelspvalsy)
sigMean = cbind(vars[Meanps<.05],Meanps[Meanps<.05])
sigMean
Meanres = data.frame(variables = vars, pvalues = Meanps)

#plot ANOVA results

pdf("Meanresults.pdf",20,7)
plot(Meanres,cex.axis=.6, ylim = 0:1)
abline(h=.05,lty="dashed", col="red")
dev.off()


#####ALEX, load data with the below command 
#####and you will have access to pvalues

load("C:/Users/Vestige/Google Drive/AlisonConsult/clinvitro/analyses/manova/MANOVA_finished.RDATA")

pvalues = data.frame(variable=vars,FullResponse = resp1ps, UncorrelatedResponse = resp2ps, AC50= AC50ps, EMAX = EMAXps ,MeanScore=Meanps)
write.csv('C:/Users/Vestige/Google Drive/AlisonConsult/clinvitro/analyses/manova/AnovaPvalues.csv', row.names = F)

# save.image(paste(outdir,"MANOVA_finished.RDATA", sep=""))



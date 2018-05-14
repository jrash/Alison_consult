list.of.packages <- c("compare","reshape","psych")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# rm(list=ls()) #################################
setwd("C:/Users/vestige/Google Drive/AlisonConsult/clinvitro/regressionTree/")
processed.dir <- 'C:/Users/Vestige/Google Drive/AlisonConsult/clinvitro/processed/clinical/'


###### WORKING WITH d5 ######

## With_genos.tsv
d <- read.table(file = paste(processed.dir,'With_genos.tsv',sep=''),sep='\t',header=TRUE,stringsAsFactors=FALSE,comment.char = "$")
## Remove unnecessary cols
d <- d[,-c(1,3,5,17,18)]
d <- d[,-c(3:13)]

cols <- colnames(d)
colnames(d) <- c('ID','Race',
                 #"CYP1B1.3","X2B6.6.516","X2B6.6.785","X2C8.3.139","X2C8.3.399","X3A4.1B","X3A5.3C","ABCB1.3435","ABCB1.2677","ABCB1.1236","GSTA1.B",
                 'Age',
                 'Menopause',
                 'Smoke',
                 'Grade',
                 'ER',
                 'ERorPR',
                 'HER2',
                 'IHC',
                 'Regimen',
                 'Adjuvant',
                 'NumCycles',
                 'StartDate',
                 'EndDate',
                 'Toxtype',
                 'ToxGrade',
                 'ToxCall',
                 'Neutropenia',
                 'Myalgia',
                 'Neuropathy',
                 'DoseInterval',
                 'TotalWeeks',
                 'FollowupStatus',
                 'StagePre',
                 'StageFinal',
                 'PercChangeRegimens',
                 'RespRegimens',
                 'PercChangeTaxane',
                 'ResponseTaxane',
                 'PercChangeNonTaxane',
                 'ResponseNonTaxane')

vars = c("Race","Age","Menopause","Smoke","Grade","ER","ERorPR","HER2","Regimen","NumCycles","Neutropenia","Myalgia","Neuropathy","DoseInterval","TotalWeeks","FollowupStatus","StagePre","StageFinal","ResponseTaxane","ResponseNonTaxane")

d = d[,c("ID",vars)]

## QC STEP ##
##### Change ER NA to 0



d$ER[which(d$ER %in% NA)] <- 0


############ Load clincal data#####################
d$ID <- paste('x',d$ID,sep='')

## Just Black and White
table(d$Race)
# Black White 
# 27    80 


table(d$Menopause)
table(d$Smoke)

#ordinal or categorical variable?
table(d$Grade)

table(d$ERorPR)
d$ERorPR[d$ERorPR=='Y'] <- 1
d$ERorPR[d$ERorPR=='N'] <- 0
table(d$ERorPR)

d$ER = as.character(d$ER)
table(d$ER)


table(d$HER2)
d$HER2[d$HER2=='Positive'] <- 2
d$HER2[d$HER2=='Borderline'] <- 1
d$HER2[d$HER2=='Negative'] <- 0
table(d$HER2)

# Not sure if I need this variable

table(d$Regimen)
# d$Regimen[!d$Regimen=='1st'] <- 1
# d$Regimen[d$Regimen=='1st'] <- 0
# table(d$Regimen)

# Not sure if I need this variable

table(d$Neutropenia)
d$Neutropenia[d$Neutropenia=='Unknown'] <- NA
d$Neutropenia[d$Neutropenia=='N'] <- 0
d$Neutropenia[d$Neutropenia=='Y'] <- 1

table(d$Myalgia)
d$Myalgia[d$Myalgia=='Unknown'] <- NA
d$Myalgia[d$Myalgia=='N'] <- 0
d$Myalgia[d$Myalgia=='Y'] <- 1

table(d$Neuropathy)
d$Neuropathy[d$Neuropathy=='Unknown'] <- NA
d$Neuropathy[d$Neuropathy=='N'] <- 0
d$Neuropathy[d$Neuropathy=='Y'] <- 1


#categorical or ordinal?
table(d$DoseInterval)
d$DoseInterval[d$DoseInterval=='Weekly'] <- 0
d$DoseInterval[d$DoseInterval=='Q2wks'] <- 1
d$DoseInterval[d$DoseInterval=='Q2.5wks'] <- 2
d$DoseInterval[d$DoseInterval=='Q3wks'] <- 3
table(d$DoseInterval)

table(d$FollowupStatus)

table(d$StagePre)
d$StagePre[d$StagePre=='IIA'] <- 0
d$StagePre[d$StagePre=='IIB'] <- 1
d$StagePre[d$StagePre=='IIIA'] <- 2
d$StagePre[d$StagePre=='IIIB'] <- 3
d$StagePre[d$StagePre=='IIIC'] <- 4
d$StagePre[d$StagePre=='IV'] <- 5
table(d$StagePre)


table(d$StageFinal)
d$StageFinal[d$StageFinal=='I'] <- 0
d$StageFinal[d$StageFinal=='IIA'] <- 1
d$StageFinal[d$StageFinal=='IIB'] <- 2
d$StageFinal[d$StageFinal=='IIIA'] <- 3
d$StageFinal[d$StageFinal=='IIIB'] <- 4
d$StageFinal[d$StageFinal=='IIIC'] <- 5
d$StageFinal[d$StageFinal=='IV'] <- 6
table(d$StageFinal)


d$ResponseTaxane = sub("UE|SD", "SDorUE", d$ResponseTaxane)
d$ResponseTaxane = sub("PD|PR", "PDorPR", d$ResponseTaxane)
table(d$ResponseTaxane)
# d$ResponseTaxane[d$ResponseTaxane=='SDorUE'] <- 0
# d$ResponseTaxane[d$ResponseTaxane=='PDorPR'] <- 1
# d$ResponseTaxane[d$ResponseTaxane=='CR'] <- 2

d$ResponseNonTaxane = sub("UE|SD", "SDorUE", d$ResponseNonTaxane)
d$ResponseNonTaxane = sub("PD|PR", "PDorPR", d$ResponseNonTaxane)
table(d$ResponseNonTaxane)
# d$ResponseNonTaxane[d$ResponseNonTaxane=='SDorUE'] <- 0
# d$ResponseNonTaxane[d$ResponseNonTaxane=='PDorPR'] <- 1
# d$ResponseNonTaxane[d$ResponseNonTaxane=='CR'] <- 2

############ Load clincal data#####################

curves.dir <- ("C:/Users/vestige/Google Drive/AlisonConsult/clinvitro/curvefits.100/")
fits <- read.table(file=paste(curves.dir,'ALL_LINES_Model_Output_File_ConfInt-Boot_PredInt-TRUE_Alpha-0.9_ConcTrans-FALSE_06May2015.csv',sep=''),sep=',',header=T,stringsAsFactors=FALSE)

hist(fits$AC50)
boxplot(fits$AC50)
hist(fits$LEC)

d$AC50 <- NA
d$EMAX <- NA
setdiff(d$ID,fits$Cell_Line)
setdiff(fits$Cell_Line,d$ID)
for (i in d$ID){
  if (dim(fits[fits$Cell_Line==i,])[1]>0){
    d$AC50[d$ID==i] <- fits$AC50[fits$Cell_Line==i]
    d$EMAX[d$ID==i] <- fits$EMAX[fits$Cell_Line==i]
  }        
}
d$mean_Score <- (d$AC50+d$EMAX)/2

nocurves <- d[-which(d$AC50%in%NA),]

dim(nocurves)
head(nocurves)
bd <- nocurves

qc <- read.csv(file="C:/Users/Vestige/Google Drive/clinvitro/qc/conc.resp/mean.qc.response.csv")
response <- qc
response$linesrepid <- paste('x',response$linesrepid,sep='')

bd = merge(bd,response,by.x="ID", by.y="linesrepid")
dim(bd)
write.csv(bd, file="C:/Users/Vestige/Google Drive/AlisonConsult/clinvitro/regressionTree/AllVarProcessReg.csv")
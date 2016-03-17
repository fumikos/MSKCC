library(Biobase)
library(CePa)
library(dplyr)
library(stats)
library(glmnet)
library(parallel)
library(doParallel)
parallel:::detectCores()
library(stringr)
library(hgu133plus2.db)
library(GenomicFeatures)
library(AnnotationDbi)
library(ggplot2)
library(caret)
library(rafalib)
library(genefilter)
library(BiocInstaller)
library(genefilter)
library(gplots)
library(RColorBrewer)
library(devtools)
library(Heatplus) 
library(vegan)
library(marray)
library(Rcpp)
library(lme4)
library(nnet)
library(quantreg)
install.packages("")
biocLite("")

#CCLE_Drug_all = read.csv("G:/My Documents/CCLE/CCLE_analysis/CCLE_NP24.2009_Drug_data_2015.02.24.csv")
#CCLE_Exp_all =read.gct("G:/My Documents/CCLE/CCLE_analysis/CCLE_Expression_Entrez_2012-09-29.gct")
#CCLE_CELL_all <- read.csv("G:/My Documents/CCLE/CCLE_analysis/Cell_tissue.csv")

idx0 <- match(CCLE_CELL_all$CCLE.name,colnames(CCLE_Exp_all))
length(which(is.na(idx0)))#17 colnames of CCLE_Exp_all have letter "X" before cell name
(correct_name<-CCLE_CELL_all$CCLE.name[which(is.na(idx0))]%>%as.vector)
CCLE_Exp_minus1 <- CCLE_Exp_all[,-991] # remove "NCIH292_LUNG.1" from CCLE_Exp_all
colnames(CCLE_Exp_minus1)[which(is.na(idx0))]
colnames(CCLE_Exp_minus1)[which(is.na(idx0))]<-correct_name
colnames(CCLE_Exp_minus1)[which(is.na(idx0))]

glioma_HSP90_index = which(grepl("CENTRAL_NERVOUS_SYSTEM", CCLE_Drug_all$CCLE.Cell.Line.Name) 
                           & CCLE_Drug_all$Target=="HSP90")# 29 glioma cells have 17AAG data
glioma_HSP90 = CCLE_Drug_all[glioma_HSP90_index,c(1,13)]

glioma_Exp_index = which(grepl("CENTRAL_NERVOUS_SYSTEM", colnames(CCLE_Exp_minus1)))
glioma_Exp = CCLE_Exp_minus1[,glioma_Exp_index]
dim(glioma_Exp)#69 glioma cells have microarray data
colnames(glioma_Exp)

(idx = match(glioma_HSP90$CCLE.Cell.Line.Name,colnames(glioma_Exp)))
length(which(is.na(idx)))#All 29 cells have microarray data
rearranged = glioma_Exp[,idx]
dim(rearranged)
colnames(rearranged)

(idx2 = match(colnames(rearranged), glioma_HSP90$CCLE.Cell.Line.Name))
length(idx2)
Y <- matrix(glioma_HSP90[idx2,2], nrow = length(idx2), ncol = 1) #29 responses
idx4<-sample(1:29,14,replace=F)
Y_train<-Y[idx4]
Y_test<-Y[-idx4]
mypar(1,3)
plot(sort(Y))
plot(sort(Y_train))
plot(sort(Y_test))

Y_sh<-sample(Y_train, replace=F)
plot(Y_train)
plot(Y_sh)

X <- matrix(rearranged, nrow = 18988, ncol = length(idx2)) 
dim(X)#18988 predictors
X<-t(X)%>%scale
dim(X)
X_train<-X[idx4,]
X_test<-X[-idx4,]
dim(X_train)
dim(X_test)

R <- cor(X,Y, method="pearson")
mypar(1,1)
hist(R,breaks=60)
idx_by_R = which(R>0.525 | R<(-0.525))
X = X[,idx_by_R]
X_train<-X_train[,idx_by_R]
X_test<-X_test[,idx_by_R]

gene_list= row.names(CCLE_Exp_all)[idx_by_R]
length(gene_list)# 103 predictors

a<-0.05
mol<-0.206
grid=exp(seq(3,-8,length=250))

fit <- glmnet(X_train, Y_train, alpha = a, lambda=grid)
Y_pred<-predict(fit, X_test, type = "response", s = mol)

length(Y_pred)
cor(Y_test, Y_pred)
cor.test(Y_test,Y_pred,method="kendall")

Y2 <-as.character(glioma_HSP90[idx2,1])
label_idx<-match(Y2, CCLE_CELL_all$CCLE.name)
(Y_names<-CCLE_CELL_all$Cell.line.primary.name[label_idx])
Y_names[-idx4]

mypar(1,1)
plot(Y_test,Y_pred,xlim=c(1.5,5.5),ylim=c(1.5,5.5),xlab="Observed",ylab="Predicted")
abline(0,1)
text(Y_test,Y_pred,labels=Y_names[-idx4],cex=0.7,pos=3)
idx4<-c(14,5, 3, 29,  7, 15, 13, 23, 24, 17, 16, 22, 12, 28)

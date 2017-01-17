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

CCLE_Drug_all = read.csv("G:/My Documents/CCLE/CCLE_analysis/CCLE_NP24.2009_Drug_data_2015.02.24.csv")
CCLE_Exp_all =read.gct("G:/My Documents/CCLE/CCLE_analysis/CCLE_Expression_Entrez_2012-09-29.gct")
CCLE_CELL_all <- read.csv("G:/My Documents/CCLE/CCLE_analysis/Cell_tissue.csv")

idx0 <- match(CCLE_CELL_all$CCLE.name,colnames(CCLE_Exp_all))
length(which(is.na(idx0)))#17 colnames of CCLE_Exp_all have letter "X" before cell name
(correct_name<-CCLE_CELL_all$CCLE.name[which(is.na(idx0))])
correct_name<-as.vector(correct_name)
correct_name
CCLE_Exp_minus1 <- CCLE_Exp_all[,-991] # remove "NCIH292_LUNG.1" from CCLE_Exp_all
colnames(CCLE_Exp_minus1)[which(is.na(idx0))]
colnames(CCLE_Exp_minus1)[which(is.na(idx0))]<-correct_name
colnames(CCLE_Exp_minus1)[which(is.na(idx0))]

(HSP90_index = which(CCLE_Drug_all$Target=="HSP90"))#503 cells have 17AAG data
HSP90 = CCLE_Drug_all[HSP90_index,c(1,13)]
(idx = match(HSP90$CCLE.Cell.Line.Name, colnames(CCLE_Exp_minus1)))
length(which(is.na(idx)))#of the 503 cells 13 do not have microarray data
HSP90[is.na(idx),]
CCLE_Exp_all[1,which(grepl("SF", colnames(CCLE_Exp_all)))]
idx <- idx[!is.na(idx)]
length(idx)#of the 503, 490 have microarray data
idx
rearranged = CCLE_Exp_minus1[,idx]
dim(rearranged)
rearranged<-rearranged[,-(which(grepl("CENTRAL_NERVOUS_SYSTEM", colnames(rearranged))))]
dim(rearranged)

glioma_HSP90_index = which(grepl("CENTRAL_NERVOUS_SYSTEM", CCLE_Drug_all$CCLE.Cell.Line.Name) 
                           & CCLE_Drug_all$Target=="HSP90")# 29 glioma cells have 17AAG data
glioma_HSP90 = CCLE_Drug_all[glioma_HSP90_index,c(1,13)]

glioma_Exp_index = which(grepl("CENTRAL_NERVOUS_SYSTEM", colnames(CCLE_Exp_minus1)))
glioma_Exp = CCLE_Exp_minus1[,glioma_Exp_index]
dim(glioma_Exp)#69 glioma cells have microarray data
colnames(glioma_Exp)

(idx_2 = match(glioma_HSP90$CCLE.Cell.Line.Name,colnames(glioma_Exp)))
length(which(is.na(idx_2)))#All 29 cells have microarray data
rearranged_g = glioma_Exp[,idx_2]
dim(rearranged_g)
colnames(rearranged_g)

(idx2 = match(colnames(rearranged), HSP90$CCLE.Cell.Line.Name))
length(idx2)
Y <- matrix(HSP90[idx2,2], nrow = length(idx2), ncol = 1) 

(idx2_2 = match(colnames(rearranged_g), glioma_HSP90$CCLE.Cell.Line.Name))
length(idx2_2)
Y_g <- matrix(glioma_HSP90[idx2_2,2], nrow = length(idx2_2), ncol = 1) #29 responses

dim(Y)
plot(sort(Y))
dim(Y_g)
plot(sort(Y_g))

X <- matrix(rearranged, nrow = 18988, ncol = length(idx2)) 
dim(X)#18988 predictors
X <- t(X)
X <- scale(X)
dim(X)

X_g <- matrix(rearranged_g, nrow = 18988, ncol = length(idx2_2)) 
dim(X_g)#18988 predictors
X_g <- t(X_g)
X_g <- scale(X_g)
dim(X_g)

R <- cor(X,Y, method="pearson")
idx_by_R = which(R>0.1 | R<(-0.1))
length(idx_by_R)
X = X[,idx_by_R]
X_g = X_g[,idx_by_R]
gene_list= row.names(CCLE_Exp_all)[idx_by_R]
length(gene_list)# 2841 predictors

a<-0.05
mol<-0.521
grid=exp(seq(3,-8,length=250))

fit <- glmnet(X, Y, alpha = a, lambda=grid)
Y_pred<-predict(fit, X_g, type = "response", s = mol)
cor(Y_g, Y_pred)
cor.test(Y_g,Y_pred,method="kendall")

Y2 <-as.character(glioma_HSP90[idx2,1])
label_idx<-match(Y2, CCLE_CELL_all$CCLE.name)

mypar(1,1)
plot(Y_g,Y_pred,xlim=c(1.5,5.5),ylim=c(1.5,5.5),xlab="Observed",ylab="Predicted")
abline(0,1)
text(Y_g,Y_pred,labels=CCLE_CELL_all$Cell.line.primary.name[label_idx],cex=0.7,pos=3)


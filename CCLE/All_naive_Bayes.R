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
library(e1071)
library(pROC)
library(OptimalCutpoints)
library(limma)
library(GEOquery)
source("http://www.bioconductor.org/biocLite.R")
install.packages("limma")
biocLite("GEOquery")

version
#CCLE_Drug_all = read.csv("G:/My Documents/CCLE/CCLE_analysis/CCLE_NP24.2009_Drug_data_2015.02.24.csv")
#CCLE_Exp_all =read.gct("G:/My Documents/CCLE/CCLE_analysis/CCLE_Expression_Entrez_2012-09-29.gct")
#CCLE_CELL_all <- read.csv("G:/My Documents/CCLE/CCLE_analysis/Cell_tissue.csv")

idx0 <- match(CCLE_CELL_all$CCLE.name,colnames(CCLE_Exp_all))
length(which(is.na(idx0)))#17 colnames of CCLE_Exp_all have letter "X" before cell name
correct_name<-as.vector(CCLE_CELL_all$CCLE.name[which(is.na(idx0))])
length(correct_name)
CCLE_Exp_minus1 <- CCLE_Exp_all[,-991] # remove "NCIH292_LUNG.1" from CCLE_Exp_all
colnames(CCLE_Exp_minus1)[which(is.na(idx0))]<-correct_name

HSP90_index = which(CCLE_Drug_all$Target=="HSP90")#503 cells have 17AAG data
HSP90 = CCLE_Drug_all[HSP90_index,c(1,13)]

idx = match(HSP90$CCLE.Cell.Line.Name,colnames(CCLE_Exp_minus1))
length(which(is.na(idx)))#of the 503 cells 13 do not have microarray data
idx <- idx[!is.na(idx)]
length(idx)#of the 503, 490 have microarray data
rearranged = CCLE_Exp_minus1[,idx]
rearranged[1:5,1:5]
HSP90[1:5,]

 idx2 = match(colnames(rearranged), HSP90$CCLE.Cell.Line.Name)
length(which(is.na(idx2)))
HSP90<-HSP90[idx2,]

summary(HSP90)
res<-HSP90[order(HSP90[,2])[1:50],]
sen<-HSP90[order(HSP90[,2])[441:490],]
res_idx<-match(res$CCLE.Cell.Line.Name,colnames(CCLE_Exp_minus1))
sen_idx<-match(sen$CCLE.Cell.Line.Name,colnames(CCLE_Exp_minus1))
mypar(1,2)
plot(sort(HSP90[,2]))
plot(c(res[,2],sen[,2]))

match(colnames(CCLE_Exp_minus1), CCLE_CELL_all$CCLE.name)
cell_name =CCLE_CELL_all$Cell.line.primary.name
colnames(CCLE_Exp_minus1)[1:5]
res_exp<-CCLE_Exp_minus1[,res_idx]
sen_exp<-CCLE_Exp_minus1[,sen_idx]
res_exp[1:5,1:5]
sen_exp[1:5,1:5]
colnames(res_exp)<-cell_name[res_idx]
colnames(sen_exp)<-cell_name[sen_idx]

# Wilcoxon Sum rank test
dim(res_exp)
dim(sen_exp)
pvals<-vector("numeric",length(res_exp[,1]))
for(i in 1:length(res_exp[,1])) {pvals[i]<-wilcox.test(res_exp[i,],sen_exp[i,])$p.value}
summary(pvals)
hist(pvals,breaks=50)
fdrs = p.adjust(pvals, method="fdr")
summary(fdrs)
plot(fdrs)
hist(fdrs)
length(which(fdrs<0.05))

# naiveBayes
X<-cbind(res_exp,sen_exp)
dim(X)
Y<-as.factor(c(rep(1, 50),rep(2, 50)))#1:res, 2:sen
X<-X[which(fdrs<0.05),]
row.names(X) <- str_replace_all(row.names(X),"_at","")
k<-as.character(row.names(X))
gene_name <- select(hgu133plus2.db, keys=k, columns=c("SYMBOL"), keytype="ENTREZID")
gene_table <- select(hgu133plus2.db, keys=k, columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
which(is.na(gene_name$SYMBOL))
row.names(X)<-gene_name[,2]
X<-X[-which(is.na(gene_name$SYMBOL)),] #remove genes without SYMBOL
head(X)
X<-t(X)

?naiveBayes
dim(X)
length(Y)
CV <- matrix(data=NA, nrow=length(Y), ncol=5)
for(j in 1:5){
  idx <- createFolds(Y, k=10)
  for(i in 1:10){
    model <- naiveBayes(X[-idx[[i]],], Y[-idx[[i]]])
    pred <- predict(model, X[idx[[i]],])
    CV[idx[[i]],j]<-pred
  }}

cor=numeric(5)
for (i in 1:5) {cor[i]<-mean(CV[,i]==Y)}
mean(cor)

# ROC analysis
X2 <- rearranged[which(fdrs<0.05),]
X2[1:5,1:5]

 <- matrix(data=NA, nrow=length(Y), ncol=5)
for(j in 1:5){
  idx <- createFolds(Y, k=10)
  for(i in 1:10){
    model <- naiveBayes(X[-idx[[i]],], Y[-idx[[i]]])
    pred <- predict(model, X[idx[[i]],])
    CV[idx[[i]],j]<-pred
  }}


?pROC
Y_pred<-round(rowMeans(CV))
roc<-roc(Y, Y_pred, plot=TRUE, ci=TRUE)
?coords
coords(roc,"l", ret=c("threshold", "specificity", "sensitivity", "ppv", "npv"))
coords(roc,"b", ret=c("threshold", "specificity", "sensitivity", "ppv", "npv"))

# GSEA
library(GSEABase)
library(GSE5859Subset)
data(GSE5859Subset)
library(sva)
library(limma)
X = sampleInfo$group
mod<-model.matrix(~X)
svafit <- sva(geneExpression,mod)

?limma
library(GEOquery)
g <- getGEO("GSE34313")
e <- g[[1]]




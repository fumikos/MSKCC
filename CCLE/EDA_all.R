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
library(caret)
library(rafalib)
library(genefilter)
library(BiocInstaller)
library(genefilter)
library(gplots)
library(RColorBrewer)
library(devtools)


#CCLE_Drug_all = read.csv("CCLE_NP24.2009_Drug_data_2015.02.24.csv")
#CCLE_Exp_all =read.gct("CCLE_Expression_Entrez_2012-09-29.gct")
#CCLE_CELL_all <- read.csv("Cell_tissue.csv")

# annotate expression data with tissue type
(idx0 <- match(CCLE_CELL_all$CCLE.name,colnames(CCLE_Exp_all)))
CCLE_CELL_all$CCLE.name[which(is.na(idx0))]
colnames(CCLE_Exp_all)[which(is.na(idx0))] #17 colnames of CCLE_Exp_all have letter "X" before cell name
match(colnames(CCLE_Exp_all),CCLE_CELL_all$CCLE.name) #which cell line is missing in CCLE_CELL_all?
colnames(CCLE_Exp_all)[991]#"NCIH292_LUNG.1" is missing

which(grepl("NCIH292", CCLE_CELL_all$CCLE.name))
CCLE_CELL_all[which(grepl("NCIH292", CCLE_CELL_all$CCLE.name)),]
which(grepl("NCIH292", colnames(CCLE_Exp_all)))# this object has two columns for NCIH292 cell
CCLE_Exp_all[1,which(grepl("NCIH292", colnames(CCLE_Exp_all)))]
cor(CCLE_Exp_all[,963],CCLE_Exp_all[,991])# the expression data for NCIH292 are similar
CCLE_Exp_minus1 <- CCLE_Exp_all[,-991] # remove "NCIH292_LUNG.1" from CCLE_Exp_all
(idx0 <- match(CCLE_CELL_all$CCLE.name, colnames(CCLE_Exp_minus1)))
match(colnames(CCLE_Exp_minus1), CCLE_CELL_all$CCLE.name)
CCLE_CELL_all[1:10,]
CCLE_Exp_minus1[1:3,1:10]
table(CCLE_CELL_all[,3])
(tissue <- CCLE_CELL_all[,3])
CCLE_Exp_minus1[1:3,1:10]
length(tissue)
table(tissue[idx])

#PCA on all tissue
s <- svd(CCLE_Exp_minus1-rowMeans(CCLE_Exp_minus1))
PC1 <- s$d[1]*s$v[,1]
PC2 <- s$d[2]*s$v[,2]
mypar(1,1)
plot(PC1,PC2,pch=21,bg=as.numeric(tissue))
legend("bottomright",levels(tissue),col=seq(along=levels(tissue)),pch=15,cex=0.5)

#PCA on selected tissue
tissues_x8 <- tissue%in%c("breast","central_nervous_system","lung","pancreas",
                      "haematopoietic_and_lymphoid_tissue","stomach","oesophagus","large_intestine")
s2 <- svd(CCLE_Exp_minus1[,tissues_x8]-rowMeans(CCLE_Exp_minus1[,tissues_x8]))
PC1 <- s2$d[1]*s2$v[,1]
PC2 <- s2$d[2]*s2$v[,2]
mypar(1,1)
plot(PC1,PC2,pch=21,bg=as.numeric(factor(tissue[tissues_x8])))
legend("bottomright",levels(factor(tissue[tissues_x8])),col=seq(along=levels(factor(tissue[tissues_x8]))),pch=15,cex=1)

#PCA on gliomas, breast, and blood Ca
gbb<-tissue%in%c("breast","central_nervous_system","haematopoietic_and_lymphoid_tissue")
s3 <- svd(CCLE_Exp_minus1[,gbb]-rowMeans(CCLE_Exp_minus1[,gbb]))
PC1 <- s3$d[1]*s3$v[,1]
PC2 <- s3$d[2]*s3$v[,2]
mypar(1,1)
plot(PC1,PC2,pch=21,bg=as.numeric(factor(tissue[gbb])))
legend("bottomright",levels(factor(tissue[gbb])),col=seq(along=levels(factor(tissue[gbb]))),pch=15,cex=1)
text(PC1,PC2,labels=CCLE_CELL_all$Cell.line.primary.name[gbb],cex=0.5,pos=3)

#PCA on gliomas only
gliomas<-tissue%in%c("central_nervous_system")
s4 <- svd(CCLE_Exp_minus1[,gliomas]-rowMeans(CCLE_Exp_minus1[,gliomas]))
PC1 <- s4$d[1]*s4$v[,1]
PC2 <- s4$d[2]*s4$v[,2]
mypar(1,1)
plot(PC1,PC2,pch=21,bg=as.numeric(factor(tissue[gliomas])))
legend("bottomright",levels(factor(tissue[gliomas])),col=seq(along=levels(factor(tissue[gliomas]))),pch=15,cex=1)
text(PC1,PC2,labels=CCLE_CELL_all$Cell.line.primary.name[gliomas],cex=0.5,pos=3)

# PCA on all cells that have 17AAG data
CCLE_Drug_all[1:10,]
HSP90_index = which(CCLE_Drug_all$Target=="HSP90")
length(HSP90_index) #503 cells were tested with 17AAG
CCLE_Drug_all[HSP90_index,13]

(idx = match(CCLE_Drug_all[HSP90_index,]$CCLE.Cell.Line.Name, CCLE_CELL_all$CCLE.name))
length(idx) 
idx <- idx[!is.na(idx)]
length(idx) # of the 503 cells tested with 17AAG, 490 has microarray data
idx # rownums for CCLE_CELL_all & colnums for CCLE_Exp_minus1

s5 <- svd(CCLE_Exp_minus1[,idx]-rowMeans(CCLE_Exp_minus1[,idx]))
PC1 <- s5$d[1]*s5$v[,1]
PC2 <- s5$d[2]*s5$v[,2]
mypar(1,1)
plot(PC1,PC2,pch=21,bg=as.numeric(factor(tissue[idx])))
legend("bottomleft",levels(factor(tissue[idx])),col=seq(along=levels(factor(tissue[idx]))),pch=15,cex=0.5)

# Boxplot for cells with 17AAG data
X <- matrix(CCLE_Exp_minus1[,idx], nrow = 18988, ncol = length(idx)) #18988 predictors
boxplot(X, range=0)

# Boxplot for glioma cells
gliomas<-tissue%in%c("central_nervous_system")
X <- matrix(CCLE_Exp_minus1[,gliomas], nrow = 18988, ncol = length(which(gliomas==T))) 
boxplot(X, range=0)

# Take a look at NQO1
HSP90_index = which(CCLE_Drug_all$Target=="HSP90")
HSP90 = CCLE_Drug_all[HSP90_index,c(1,13)]

match(colnames(CCLE_Exp_minus1[,idx]),CCLE_CELL_all$CCLE.name[idx])
(idx2 = match(CCLE_CELL_all$CCLE.name[idx], HSP90$CCLE.Cell.Line.Name))# rownums for HSP90
length(idx2)
Y <- matrix(HSP90[idx2,2], nrow = length(idx2), ncol = 1) # Activity area
#Y <- sample(Y, replace=F)
dim(Y) #490 responses
plot(Y)
class(Y)

(e_1728 <- rearranged["1728_at",])
qqnorm(e_1728)
qqline(e_1728)
plot(e_1728,Y)
cor(e_1728,Y)




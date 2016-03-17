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
#17 colnames of CCLE_Exp_all have letter "X" before cell name
correct_name<-CCLE_CELL_all$CCLE.name[which(is.na(idx0))]%>%as.vector
CCLE_Exp_minus1 <- CCLE_Exp_all[,-991] # remove "NCIH292_LUNG.1" from CCLE_Exp_all
colnames(CCLE_Exp_minus1)[which(is.na(idx0))]<-correct_name

glioma_HSP90_index = which(grepl("CENTRAL_NERVOUS_SYSTEM",CCLE_Drug_all$CCLE.Cell.Line.Name) 
                                & CCLE_Drug_all$Target=="HSP90")# 29 glioma cells have 17AAG data
glioma_HSP90 = CCLE_Drug_all[glioma_HSP90_index,c(1,13)]

glioma_Exp_index = which(grepl("CENTRAL_NERVOUS_SYSTEM", colnames(CCLE_Exp_minus1)))
glioma_Exp = CCLE_Exp_minus1[,glioma_Exp_index]
dim(glioma_Exp)#69 glioma cells have microarray data
colnames(glioma_Exp)
(glioma_name =CCLE_CELL_all$Cell.line.primary.name[glioma_Exp_index])
(idx = match(glioma_HSP90$CCLE.Cell.Line.Name,colnames(glioma_Exp)))
length(which(is.na(idx)))#All 29 cells have microarray data
rearranged = glioma_Exp[,idx]
dim(rearranged)
colnames(rearranged)

(idx2 = match(colnames(rearranged), glioma_HSP90$CCLE.Cell.Line.Name))
length(idx2)
Y <- matrix(glioma_HSP90[idx2,2], nrow = length(idx2), ncol = 1) #29 responses

dim(Y)
par(mfrow=c(1,2))
plot(sort(Y))

dim(glioma_HSP90)
res<-glioma_HSP90[order(glioma_HSP90[,2])[1:10],]
sen<-glioma_HSP90[order(glioma_HSP90[,2])[20:29],]
res_idx<-match(res$CCLE.Cell.Line.Name,colnames(glioma_Exp))
sen_idx<-match(sen$CCLE.Cell.Line.Name,colnames(glioma_Exp))
colnames(glioma_Exp)<-glioma_name
res_exp<-glioma_Exp[,res_idx]
sen_exp<-glioma_Exp[,sen_idx]
colnames(res_exp)
colnames(sen_exp)

# PCA for resistant and sensitive glioma cells
pheno<-vector("numeric",69)
pheno[res_idx]<-rep(1,length(res_idx))
pheno[sen_idx]<-rep(2,length(sen_idx))
s <- svd(glioma_Exp-rowMeans(glioma_Exp))
PC1 <- s$d[1]*s$v[,1]
PC2 <- s$d[2]*s$v[,2]
mypar(1,1)
plot(PC1,PC2,pch=21,bg=pheno)
legend("bottomright",c("others","resistant","sensitive"),col=c(0,1,2),pch=15,cex=0.8)
text(PC1,PC2,labels=glioma_name,cex=0.5,pos=3)

#Spearman correlation
X <- matrix(rearranged, nrow = 18988, ncol = length(idx2)) 
X <- t(X)
X <- scale(X)
R <- cor(X,Y, method="spearman")
hist(R, breaks=60)
idx_by_R = which(R>0.1 | R<(-0.1))
length(idx_by_R)
res_exp<-res_exp[idx_by_R,]
sen_exp<-sen_exp[idx_by_R,]

# Wilcoxon Sum rank test
dim(res_exp)
dim(sen_exp)
pvals<-vector("numeric",length(idx_by_R))
for(i in 1:length(idx_by_R)) {pvals[i]<-wilcox.test(res_exp[i,],sen_exp[i,])$p.value}
summary(pvals)
hist(pvals,breaks=50)
fdrs = p.adjust(pvals,method="fdr")
summary(fdrs)
plot(fdrs)
which(fdrs<0.45)# no hit

which(pvals<=0.05)

sort(pvals)[1:40]
which(pvals<=0.0006)
gene_list= row.names(res_exp)
length(gene_list)
gene_list <- str_replace_all(gene_list,"_at","")
length(gene_list)
gene_table <- select(hgu133plus2.db, keys=gene_list, columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
gene_table[which(fdrs<0.45),]

exp<-glioma_Exp[,c(res_idx,sen_idx)]
dim(exp)
exp[1:5,]
rownames(exp)<-gene_table[,2]
exp<-t(scale(t(exp)))
head(rowMeans(exp))
head(rowSds(exp))
dim(exp)

mypar(1,1)
cols<-colorRampPalette(c("green","black","red"))
brks<-unique(c(seq(-5,-1,length=100),seq(-1,1,length=50),seq(1,5,length=100)))
brks2<-unique(c(seq(-3.5,3.5,length=100)))

heatmap.2(exp[which(pvals<=0.0001),],dendrogram="none",
          col=cols,breaks=brks2,trace="none")

fdrs = p.adjust(pvals,method="fdr")
plot(fdrs)
which(fdrs<0.25)# no hit

# try t-test
pvals<-vector("numeric",18988)
for(i in 1:18988) {pvals[i]<-t.test(res_exp[i,],sen_exp[i,])$p.value}
summary(pvals)
hist(pvals,breaks=50)
pvals_0.05<-pvals[which(pvals<=0.05)]


fdrs = p.adjust(pvals_0.05,method="fdr")
summary(fdrs)
hist(fdrs)
which(fdrs<0.25)

?p.adjust

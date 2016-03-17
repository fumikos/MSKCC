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
#(Y <- sample(Y, replace=F))
dim(Y)
par(mfrow=c(1,2))
plot(sort(Y))
plot(sort(glioma_HSP90[,2]))

X <- matrix(rearranged, nrow = 18988, ncol = length(idx2)) 
dim(X)#18988 predictors
hist(rowMeans(X),breaks=60)
idx_by_rm<-which(rowMeans(X)>4.5)
length(idx_by_rm)
X <- X[idx_by_rm,]%>%t%>%scale
#X_s<-t(X)%>%scale
dim(X)
X[1:10,1:10]

R <- cor(X,Y, method="pearson")
hist(R,breaks=60)
idx_by_R = which(R>0.52 | R<(-0.52))
X = X[,idx_by_R]
dim(X)
#sort(R)[50]

gene_list= row.names(CCLE_Exp_all)[idx_by_rm]
gene_list<-gene_list[idx_by_R]
length(gene_list)# 12292 predictors

# Regression analysis
# Optimize alpha below
grid=exp(seq(3,-8,length=250))
fold_id=sample(rep(seq(8),length=length(Y)))
table(fold_id)

alphaslist<-seq(0.05, 0.95, length=10)
alphaslist
registerDoParallel(4)
elasticnet<-lapply(alphaslist, function(a){cv.glmnet(X, Y, alpha=a, nfolds=8, foldid=fold_id, lambda=grid, parallel = TRUE)})
mypar(3,4)
for (i in 1:10) {plot(elasticnet[[i]])}
for (i in 1:10) {print(min(elasticnet[[i]]$cvm))}
# -> optimal alpha = 0.05
a<-0.05
for (i in 1:10) {print(elasticnet[[i]]$lambda.min)}
for (i in 1:10) {print(elasticnet[[i]]$lambda.1se)}

# Optimize lambda below
registerDoParallel(4)
lambdas <- replicate (50, {
  elasticnet<-cv.glmnet(X, Y, alpha = a, nfolds = 8, lambda=grid, parallel = TRUE)
  elasticnet$lambda.min
})
lambdas
sd(lambdas)
(mol <- mean(lambdas))
mol<-0.521

# boostrap 200 times
grid=exp(seq(3,-8,length=250))
boot_betas <- matrix(data=NA, nrow=ncol(X)+1, ncol=200)
X_resample <- matrix(data=NA, nrow=length(Y), ncol=ncol(X))
Y_resample <- matrix(data=NA, nrow=length(Y), ncol=1)
for(i in 1:200){idx <- sample(length(Y), replace = T)
for(j in 1:length(Y)){X_resample[j,] <- X[idx[j],]; Y_resample[j] <- Y[idx[j]]}
fit <- glmnet(X_resample, Y_resample, alpha = a, lambda=grid)
boot_betas[,i] <- as.vector(coef(fit, s = mol, exact = FALSE))
}
head(boot_betas)
dim(boot_betas)
boot_betas <- boot_betas[2:(ncol(X)+1),]
head(boot_betas)
dim(boot_betas)

top_feats <- matrix(data=NA, nrow=ncol(X), ncol=1)
for (i in 1:ncol(X)) { top_feats[i] <- sum(boot_betas[i,]!=0)>160}
dim(top_feats)
head(top_feats)
sum(top_feats)

length(gene_list)
(k <- gene_list[which(top_feats==T)])
(k <- str_replace_all(k,"_at",""))
length(k)
(k <- as.character(k))
gene_table <- select(hgu133plus2.db, keys=k, columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")

(betas_mean <- apply(boot_betas,1,function(x) mean(x[x != 0])))
(betas_sd <- apply(boot_betas,1,function(x) sd(x[x != 0])))
(betas_weight <- apply(boot_betas,1,function(x) sum(x != 0)))
table(betas_weight)

betas_weight[which(top_feats==T)]
betas_mean[which(top_feats==T)]

dim(gene_table)
gene_table[,4] <- betas_weight[which(top_feats==T)]
colnames(gene_table)[4] <- "hits_BS"
gene_table[,5] <- betas_mean[which(top_feats==T)]
colnames(gene_table)[5] <- "beta_mean"

gene_table[,4]<-gene_table[,4]/200
(result<-gene_table[order(gene_table[,5]),])
write.csv(result, file="glioma_vs_HSP90_result.csv")

# boostrap 200 times again
boot_betas2 <- matrix(data=NA, nrow=ncol(X)+1, ncol=200)
X_resample2 <- matrix(data=NA, nrow=length(Y), ncol=ncol(X))
Y_resample2 <- matrix(data=NA, nrow=length(Y), ncol=1)
for(i in 1:200){idx <- sample(length(Y), replace = T)
for(j in 1:length(Y)){X_resample2[j,] <- X[idx[j],]; Y_resample2[j] <- Y[idx[j]]}
fit <- glmnet(X_resample2, Y_resample2, alpha = a, lambda=grid)
boot_betas2[,i] <- as.vector(coef(fit, s = mol, exact = FALSE))
}
head(boot_betas2)
dim(boot_betas2)
boot_betas2 <- boot_betas2[2:(ncol(X)+1),]
head(boot_betas2)
dim(boot_betas2)

top_feats2 <- matrix(data=NA, nrow=ncol(X), ncol=1)
for (i in 1:ncol(X)) { top_feats2[i] <- sum(boot_betas2[i,]!=0)>160}
dim(top_feats2)
head(top_feats2)
sum(top_feats2)

length(gene_list)
(k2 <- gene_list[which(top_feats2==T)])
(k2 <- str_replace_all(k2,"_at",""))
length(k2)
(k2 <- as.character(k2))
gene_table2 <- select(hgu133plus2.db, keys=k2, columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")

(betas_mean2 <- apply(boot_betas2,1,function(x) mean(x[x != 0])))
(betas_sd2 <- apply(boot_betas2,1,function(x) sd(x[x != 0])))
(betas_weight2 <- apply(boot_betas2,1,function(x) sum(x != 0)))
table(betas_weight2)

betas_weight2[which(top_feats2==T)]
betas_mean2[which(top_feats2==T)]

dim(gene_table2)
gene_table2[,4] <- betas_weight2[which(top_feats2==T)]
colnames(gene_table2)[4] <- "hits_BS"
gene_table2[,5] <- betas_mean2[which(top_feats2==T)]
colnames(gene_table2)[5] <- "beta_mean"

gene_table2[,4]<-gene_table2[,4]/200
result2<-gene_table2[order(gene_table2[,5]),]
write.csv(result2, file="glioma_vs_HSP90_result2.csv")

# boostrap 200 times yet again
boot_betas3 <- matrix(data=NA, nrow=ncol(X)+1, ncol=200)
X_resample3 <- matrix(data=NA, nrow=length(Y), ncol=ncol(X))
Y_resample3 <- matrix(data=NA, nrow=length(Y), ncol=1)
for(i in 1:200){idx <- sample(length(Y), replace = T)
for(j in 1:length(Y)){X_resample3[j,] <- X[idx[j],]; Y_resample3[j] <- Y[idx[j]]}
fit <- glmnet(X_resample3, Y_resample3, alpha = a, lambda=grid)
boot_betas3[,i] <- as.vector(coef(fit, s = mol, exact = FALSE))
}
head(boot_betas3)
dim(boot_betas3)
boot_betas3 <- boot_betas3[2:(ncol(X)+1),]
head(boot_betas3)
dim(boot_betas3)

top_feats3 <- matrix(data=NA, nrow=ncol(X), ncol=1)
for (i in 1:ncol(X)) { top_feats3[i] <- sum(boot_betas3[i,]!=0)>160}
dim(top_feats3)
head(top_feats3)
sum(top_feats3)

length(gene_list)
(k3 <- gene_list[which(top_feats3==T)])
(k3 <- str_replace_all(k3,"_at",""))
length(k3)
(k3 <- as.character(k3))
gene_table3 <- select(hgu133plus2.db, keys=k3, columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")

(betas_mean3 <- apply(boot_betas3,1,function(x) mean(x[x != 0])))
(betas_sd3 <- apply(boot_betas3,1,function(x) sd(x[x != 0])))
(betas_weight3 <- apply(boot_betas3,1,function(x) sum(x != 0)))
table(betas_weight3)

betas_weight3[which(top_feats3==T)]
betas_mean3[which(top_feats3==T)]

dim(gene_table3)
gene_table3[,4] <- betas_weight3[which(top_feats3==T)]
colnames(gene_table3)[4] <- "hits_BS"
gene_table3[,5] <- betas_mean3[which(top_feats3==T)]
colnames(gene_table3)[5] <- "beta_mean"

gene_table3[,4]<-gene_table3[,4]/200
result3<-gene_table3[order(gene_table3[,5]),]
write.csv(result3, file="glioma_vs_HSP90_result3.csv")

# heatmap
result_at<-read.csv("glioma_vs_HSP90_result_at.csv")
result_at
dim(result_at)

match(CCLE_CELL_all$CCLE.name,colnames(CCLE_Exp_minus1))# 100% match
(glioma_Exp_index = which(grepl("CENTRAL_NERVOUS_SYSTEM", colnames(CCLE_Exp_minus1))))
glioma_Exp = CCLE_Exp_minus1[,glioma_Exp_index]
colnames(glioma_Exp)
(glioma_name =CCLE_CELL_all$Cell.line.primary.name[glioma_Exp_index])
(idx = match(glioma_HSP90$CCLE.Cell.Line.Name,colnames(glioma_Exp)))
rearranged = glioma_Exp[,idx]
colnames(rearranged)
(glioma_name<-glioma_name[idx])# matches "rearranged"
colnames(rearranged)<-glioma_name
colnames(rearranged)

order_gene<-match(result_at[,1],rownames(rearranged))
order_gene

mypar(1,1)
plot(Y)
plot(Y[order(Y)])

exp<-rearranged[order_gene,order(Y)]
dim(exp)
exp<-t(scale(t(exp)))
rownames(exp)<-result_at[,2]
exp[1:5,1:5]
tail(exp)

rowMeans(exp)
rowSds(exp)

mypar(1,1)
cols<-colorRampPalette(c("green","black","red"))
brks<-unique(c(seq(-5,-1,length=100),seq(-1,1,length=50),seq(1,5,length=100)))
brks2<-unique(c(seq(-3.5,3.5,length=100)))

heatmap.2(exp,Rowv=FALSE,Colv=FALSE,dendrogram="none",
          col=cols,breaks=brks2,trace="none")
?heatmap.2

result_at
barplot(result_at[,5])
?plot

head(HSP90[order(HSP90[,2],decreasing = TRUE),])
tail(exp)
sd(exp[1,])
mean(exp[1,])
hist(log10(exp[1,]))

mypar(1,1)
boxplot(t(exp))

# compare with r
R <- cor(X,Y, method="pearson")
dim(R)
R[which(top_feats==T)]
summary(R)

k4<-c("MSI1","ATP10A","TM2D2","LHFPL1","CDHR3","C12orf57","INMT","CSNK1E","CYCS","NQO1")
(published <- select(hgu133plus2.db, keys=k4, columns=c("ENTREZID","GENENAME","SYMBOL"), keytype="SYMBOL"))

cor(rearranged["4440_at",],Y)
cor(rearranged["57194_at",],Y)
cor(rearranged["83877_at",],Y)
cor(rearranged["340596_at",],Y)
cor(rearranged["222256_at",],Y)#CDHR3 has R=-0.09591676
cor(rearranged["113246_at",],Y)
cor(rearranged["11185_at",],Y)
cor(rearranged["1454_at",],Y)
cor(rearranged["54205_at",],Y)
cor(rearranged["1728_at",],Y)

# cross validation (Barretina et al.)
CV_B <- matrix(data=NA, nrow=length(Y), ncol=8)
for(j in 1:8){
  idx <- createFolds(Y, k=8)
  for(i in 1:8){
    fit <- glmnet(X[-idx[[i]], ], Y[-idx[[i]]], alpha = a, lambda=grid)
    pred<-predict(fit, X[idx[[i]], ], type = "response", s = mol)
    CV_B[idx[[i]],j]<-pred
  }}

head(CV_B)
dim(CV_B)
Y_pred <- rowMeans(CV_B)
head(Y_pred)
length(Y_pred)
cor(Y, Y_pred)
cor.test(Y,Y_pred,method="kendall")

Y2 <-as.character(glioma_HSP90[idx2,1])
label_idx<-match(Y2, CCLE_CELL_all$CCLE.name)

mypar(1,1)
plot(Y,Y_pred,xlim=c(1.5,5.5),ylim=c(1.5,5.5),xlab="Observed",ylab="Predicted")
abline(0,1)
text(Y,Y_pred,labels=CCLE_CELL_all$Cell.line.primary.name[label_idx],cex=0.7,pos=3)
?plot



library(rafalib)
library(stringr)
library(hgu133plus2.db)
library(e1071)
library(caret)
library(pROC)

source("http://www.bioconductor.org/biocLite.R")
install.packages("")
biocLite("")

#Extract gene expression and drug response data from CCLE data

library(rafalib)
library(stringr)
library(hgu133plus2.db)

#CCLE_Drug_all = read.csv("G:/My Documents/CCLE/CCLE_data/CCLE_NP24.2009_Drug_data_2015.02.24.csv")
#CCLE_Exp_all = read.gct("G:/My Documents/CCLE/CCLE_data/CCLE_Expression_Entrez_2012-09-29.gct")
#CCLE_CELL_all <- read.csv("G:/My Documents/CCLE/CCLE_data/Cell_tissue.csv")

#First, check if the cell ID matches in CCLE_CELL_all and CCLE_Exp_all
idx0 <- match(CCLE_CELL_all$CCLE.name,colnames(CCLE_Exp_all))
length(which(is.na(idx0)))
CCLE_CELL_all$CCLE.name[which(is.na(idx0))]
colnames(CCLE_Exp_all)[which(is.na(idx0))] #17 colnames of CCLE_Exp_all have letter "X" before cell name
length(CCLE_CELL_all$CCLE.name)# 1036
length(colnames(CCLE_Exp_all))# 1037 CCLE_Exp_all also has an extra column
CCLE_Exp_minus1 <- CCLE_Exp_all[,-991] # remove "NCIH292_LUNG.1" from CCLE_Exp_all
correct_name<-as.vector(CCLE_CELL_all$CCLE.name[which(is.na(idx0))])
length(correct_name)
colnames(CCLE_Exp_minus1)[which(is.na(idx0))]<-correct_name# colname for CCLE_Exp_minus1 now matches CCLE_CELL_all$CCLE.name

#Next, extract Activity Area data for HSP90 inhibitor
head(CCLE_Drug_all)
HSP90_index = which(CCLE_Drug_all$Target=="HSP90")#503 cells have 17AAG data
HSP90 = CCLE_Drug_all[HSP90_index,c(1,13)]
head(HSP90)# cell name and ActArea

#Next, match cell name in HSP90 and CCLE_Exp_minus1
idx = match(HSP90$CCLE.Cell.Line.Name,colnames(CCLE_Exp_minus1))
length(which(is.na(idx)))#of the 503 cells 13 cell names do not exist in CCLE_Exp_minus1
HSP90$CCLE.Cell.Line.Name[which(is.na(idx))]
sum(colnames(CCLE_Exp_minus1)=="BGC823_STOMACH")#confirm above
idx <- idx[!is.na(idx)]
length(idx)#of the 503, 490 have microarray data
rearranged = CCLE_Exp_minus1[,idx]#extract expression data of cells that have HSP90i data and reorder to match HSP90 data order
dim(rearranged)
rearranged[1:5,1:5]
dim(HSP90)
HSP90[1:5,]

#Remove cell lines with no microarray data from HSP90 data
idx2 = match(colnames(rearranged), HSP90$CCLE.Cell.Line.Name)
length(which(is.na(idx2)))
HSP90<-HSP90[idx2,]
dim(HSP90)

#Assign resistant and sensitive cells
summary(HSP90)
res<-HSP90[order(HSP90[,2])[1:50],]
sen<-HSP90[order(HSP90[,2])[441:490],]
res_idx<-match(res$CCLE.Cell.Line.Name,colnames(CCLE_Exp_minus1))
sen_idx<-match(sen$CCLE.Cell.Line.Name,colnames(CCLE_Exp_minus1))
mypar(1,2)
plot(sort(HSP90[,2]))
plot(c(res[,2],sen[,2]))

#Simplify cell name
class(colnames(CCLE_Exp_minus1))
class(CCLE_CELL_all$CCLE.name)
identical(colnames(CCLE_Exp_minus1), as.character(CCLE_CELL_all$CCLE.name))
colnames(CCLE_Exp_minus1)<-CCLE_CELL_all$Cell.line.primary.name
colnames(CCLE_Exp_minus1)[1:5]

#Extract microarray data for resistant and sensitive cells and combine into a matrix
res_exp<-CCLE_Exp_minus1[,res_idx]
sen_exp<-CCLE_Exp_minus1[,sen_idx]
res_exp[1:5,1:5]
sen_exp[1:5,1:5]
X<-cbind(res_exp,sen_exp)
X[1:5,1:5]
dim(X)

# Wilcoxon Sum rank test

?wilcox.test

dim(res_exp)
dim(sen_exp)
pvals_sen<-vector("numeric",length(sen_exp[,1]))
for(i in 1:length(sen_exp[,1])) {
  pvals_sen[i]<-wilcox.test(sen_exp[i,],res_exp[i,],alternative ="greater")$p.value
}
mypar(1,1)
hist(pvals_sen,breaks=50)
plot(pvals_sen)
fdrs_sen = p.adjust(pvals_sen, method="fdr")
hist(fdrs_sen)
plot(fdrs_sen)
cutoffs<-seq(0,0.25,0.01)
gene_len_sen<-vector("numeric",length(cutoffs))
for(i in 1:length(cutoffs)){
  gene_len_sen[i]<-length(which(fdrs_sen<cutoffs[i]))
}
plot(cutoffs,gene_len_sen)
abline(100,0)

length(which(fdrs_sen<0.05))# 100 genes
length(order(fdrs_sen)[1:100])# index of top 100 genes
fdrs_sen[order(fdrs_sen)[1:101]]# 100th and 101st gene not a tie

pvals_res<-vector("numeric",length(res_exp[,1]))
for(i in 1:length(res_exp[,1])) {
  pvals_res[i]<-wilcox.test(res_exp[i,],sen_exp[i,],alternative ="greater")$p.value
}
hist(pvals_res,breaks=50)
plot(pvals_res)
fdrs_res = p.adjust(pvals_res, method="fdr")
hist(fdrs_res)
plot(fdrs_res)
gene_len_res<-vector("numeric",length(cutoffs))
for(i in 1:length(cutoffs)){
  gene_len_res[i]<-length(which(fdrs_res<cutoffs[i]))
}
plot(cutoffs,gene_len_res)
abline(100,0)

length(order(fdrs_res)[1:100])# index of top 100 genes
fdrs_res[order(fdrs_res)[1:101]]# 100th and 101st gene not a tie

# Save gene symbol

k <- str_replace_all(row.names(CCLE_Exp_minus1),"_at","")# probe ID to Entrez ID
gene_name <- select(hgu133plus2.db, keys=k, columns=c("SYMBOL"), keytype="ENTREZID")# Entrez ID to symbol
gene_name<-cbind(gene_name, fdrs_sen, fdrs_res)
head(gene_name)

sen_sig<-gene_name[order(fdrs_sen)[1:100],]
sen_sig
write.table(sen_sig$SYMBOL, file="sen_sig_100_wilcox.txt", quote=F, row.names=F, col.names=F)

res_sig<-gene_name[order(fdrs_res)[1:100],]
res_sig
write.table(res_sig$SYMBOL, file="res_sig_100_wilcox.txt", quote=F, row.names=F, col.names=F)

# naive Bayes classification

library(e1071)
library(stringr)
library(hgu133plus2.db)
library(caret)

?naiveBayes
dim(X)
Y<-as.factor(c(rep(1, 50),rep(2, 50)))#1:res, 2:sen
X1<-X[c(order(fdrs_sen)[1:100], order(fdrs_res)[1:100]),]
dim(X1)# 200 genes
X1[1:10,1:5]
X1<-t(X1)

# Cross Validation

dim(X1)
length(Y)
CV1 <- matrix(data=NA, nrow=length(Y), ncol=5)
for(j in 1:5){
  idx <- createFolds(Y, k=10)
  for(i in 1:10){
    model <- naiveBayes(X1[-idx[[i]],], Y[-idx[[i]]])
    pred <- predict(model, X1[idx[[i]],], type = "raw")
    CV1[idx[[i]],j]<-pred[,2]
  }}

# ROC analysis

library(pROC)

?pROC
Y_pred<-rowMeans(CV1)
roc1<-roc(Y, Y_pred, plot=TRUE, ci=TRUE)
auc(roc1)
# Area under the curve: 0.9044
ci(roc1)
# 95% CI: 0.8432-0.9656 (DeLong)

?coords
coords(roc1,"l", ret=c("threshold", "specificity", "sensitivity", "ppv", "npv"))
coords(roc1,"b", ret=c("threshold", "specificity", "sensitivity", "ppv", "npv"))

# Predict all cell lines (including sen and res)
# Removing NQO1 (1728_at) did not significantly change the results

X2 <- rearranged[c(order(fdrs_sen)[1:100], order(fdrs_res)[1:100]),]
dim(X2)
X2<-t(X2)
dim(X2)
X1[1:5,1:5]
X2[1:5,1:5]
#select(hgu133plus2.db, keys="NQO1", columns=c("ENTREZID"), keytype="SYMBOL")
#X1<-X1[,-which(colnames(X1)=="1728_at")]
#X2<-X2[,-which(colnames(X2)=="1728_at")]
#dim(X1)
#dim(X2)
identical(colnames(X2),colnames(X1))
model2 <- naiveBayes(X1, Y)
pred2 <- predict(model2, X2, type = "raw")
sen_pred<-pred2[,2]
res_pred<-pred2[,1]
mypar(1,2)
plot(res_pred, col=(as.numeric(res_pred==1)+1))
plot(sen_pred, col=(as.numeric(sen_pred==1)+1))
length(HSP90[,2])
mypar(1,1)
plot(HSP90[,2], col=(as.numeric(sen_pred==1)+1))
plot(HSP90[,2], col=(as.numeric(res_pred==1)+1))
plot(HSP90[,2], col=(as.numeric(sen_pred!=1&res_pred!=1)+1))
plot(density(HSP90[,2]))
lines(density(HSP90[sen_pred==1,2]), col=2)
lines(density(HSP90[res_pred==1,2]), col=3)
lines(density(HSP90[sen_pred!=1&res_pred!=1,2]), col=4)
t.test(HSP90[sen_pred==1,2],HSP90[res_pred==1,2])$p.value
# [1] 2.53978e-27
t.test(HSP90[sen_pred==1,2],HSP90[sen_pred!=1&res_pred!=1,2])$p.value
# [1] 9.341802e-06
t.test(HSP90[res_pred==1,2],HSP90[sen_pred!=1&res_pred!=1,2])$p.value
# [1] 1.735112e-14

mypar(1,2)
plot(sort(HSP90[,2]),col=(as.numeric(sen_pred==1)[order(HSP90[,2])]))
plot(sort(HSP90[,2]),col=(as.numeric(res_pred==1)[order(HSP90[,2])]))

# Predict all cell lines (excluding sen and res)
# Removing NQO1 (1728_at) did not significantly change the results

(res_idx2<-match(res$CCLE.Cell.Line.Name,colnames(rearranged)))
(match(res$CCLE.Cell.Line.Name,HSP90$CCLE.Cell.Line.Name))# confirm match
(sen_idx2<-match(sen$CCLE.Cell.Line.Name,colnames(rearranged)))
(match(sen$CCLE.Cell.Line.Name,HSP90$CCLE.Cell.Line.Name))
X3 <- rearranged[c(order(fdrs_sen)[1:100], order(fdrs_res)[1:100]),-c(res_idx2,sen_idx2)]
HSP90_2<-HSP90[-c(res_idx2,sen_idx2),]
mypar(1,1)
plot(sort(HSP90[,2]))
points(sort(HSP90_2[,2]),col=1)
X3[1:5,1:5]
dim(X3)
dim(HSP90_2)
X3<-t(X3)
dim(X3)
#X3<-X3[,-which(colnames(X3)=="1728_at")]
#dim(X3)
identical(colnames(X3),colnames(X1))
pred3 <- predict(model2, X3, type = "raw")
sen_pred2<-pred3[,2]
res_pred2<-pred3[,1]
length(sen_pred2)
mypar(1,2)
plot(res_pred2, col=(as.numeric(res_pred2==1)+1))
plot(sen_pred2, col=(as.numeric(sen_pred2==1)+1))
length(HSP90_2[,2])
mypar(1,1)
plot(HSP90_2[,2], col=(as.numeric(sen_pred2==1)+1))
plot(HSP90_2[,2], col=(as.numeric(res_pred2==1)+1))
plot(HSP90_2[,2], col=(as.numeric(sen_pred2!=1&res_pred2!=1)+1))
plot(density(HSP90_2[,2]))
lines(density(HSP90_2[sen_pred2==1,2]), col=2)
lines(density(HSP90_2[res_pred2==1,2]), col=3)
lines(density(HSP90_2[sen_pred2!=1&res_pred2!=1,2]), col=4)
t.test(HSP90_2[sen_pred2==1,2],HSP90_2[res_pred2==1,2])$p.value
# [1] 3.727346e-08
t.test(HSP90_2[sen_pred2==1,2],HSP90_2[sen_pred2!=1&res_pred2!=1,2])$p.value
# [1] 0.0685431
t.test(HSP90_2[res_pred2==1,2],HSP90_2[sen_pred2!=1&res_pred2!=1,2])$p.value
# [1] 2.85561e-05

mypar(1,2)
plot(sort(HSP90_2[,2]),col=(as.numeric(sen_pred2==1)[order(HSP90_2[,2])]))
plot(sort(HSP90_2[,2]),col=(as.numeric(res_pred2==1)[order(HSP90_2[,2])]))

# Predict glioma cell lines (including res and sen)

glioma_HSP90_idx = which(grepl("CENTRAL_NERVOUS_SYSTEM",CCLE_Drug_all$CCLE.Cell.Line.Name) 
                         & CCLE_Drug_all$Target=="HSP90")
length(glioma_HSP90_idx)# 29 glioma cells have 17AAG data
glioma_HSP90 = CCLE_Drug_all[glioma_HSP90_idx,c(1,2,13)]
glioma_HSP90[1:5,]
dim(glioma_HSP90)# 29 glioma cells with 17AAD data
dim(rearranged)
glioma_Exp_idx = which(grepl("CENTRAL_NERVOUS_SYSTEM", colnames(rearranged)))
glioma_Exp = rearranged[,glioma_Exp_idx]
dim(glioma_Exp)
match(colnames(glioma_Exp), glioma_HSP90[,1])# the order of cells in glioma_Exp and glioma_HSP90 match
match(glioma_HSP90[,2],colnames(X))# 3 glioma cells were included in the traiuning set

X4 <- glioma_Exp[c(order(fdrs_sen)[1:100], order(fdrs_res)[1:100]),]
dim(X4)
X4<-t(X4)
dim(X4)
X4[1:5,1:5]

identical(colnames(X4),colnames(X1))
pred4 <- predict(model2, X4, type = "raw")
sen_pred3<-pred4[,2]
res_pred3<-pred4[,1]
mypar(1,2)
plot(res_pred3, col=(as.numeric(res_pred3==1)+1))
plot(sen_pred3, col=(as.numeric(sen_pred3==1)+1))
length(glioma_HSP90[,3])
mypar(1,1)
plot(glioma_HSP90[,3], col=(as.numeric(sen_pred3==1)+1))
plot(glioma_HSP90[,3], col=(as.numeric(res_pred3==1)+1))
plot(glioma_HSP90[,3], col=(as.numeric(sen_pred3!=1&res_pred3!=1)+1))
plot(density(glioma_HSP90[,3]))
lines(density(glioma_HSP90[sen_pred3==1,3]), col=2)
lines(density(glioma_HSP90[res_pred3==1,3]), col=3)
lines(density(glioma_HSP90[sen_pred3!=1&res_pred3!=1,3]), col=4)
t.test(glioma_HSP90[sen_pred3==1,3],glioma_HSP90[res_pred3==1,3])$p.value
# Error
t.test(glioma_HSP90[sen_pred3==1,3],glioma_HSP90[sen_pred3!=1&res_pred3!=1,3])$p.value
# Error
t.test(glioma_HSP90[res_pred3==1,3],glioma_HSP90[sen_pred3!=1&res_pred3!=1,3])$p.value
# 0.02652505

mypar(1,2)
plot(sort(glioma_HSP90[,3]),col=(as.numeric(sen_pred3==1)[order(glioma_HSP90[,3])]))
plot(sort(glioma_HSP90[,3]),col=(as.numeric(res_pred3==1)[order(glioma_HSP90[,3])]))



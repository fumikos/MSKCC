### Sensitivity prediction using naive Bayes classifier

## Load gene expression and drug response data from CCLE data
library(CePa)
CCLE_Drug_all <- read.csv("G:/My Documents/CCLE/CCLE_data/CCLE_NP24.2009_Drug_data_2015.02.24.csv")
CCLE_Exp_all <- read.gct("G:/My Documents/CCLE/CCLE_data/CCLE_Expression_Entrez_2012-09-29.gct")
CCLE_CELL_all <- read.csv("G:/My Documents/CCLE/CCLE_data/Cell_tissue.csv")
head(CCLE_Drug_all)
CCLE_Exp_all[1:5,1:5]
head(CCLE_CELL_all)

## Run code in "CCLE_17AAG_data.R"
dim(CCLE_Exp_minus1)# 18988  1036
CCLE_Exp_minus1[1:5,1:5]# expression data for all cells
dim(rearranged)# 18988   492
rearranged[1:5,1:5]# expression data for 17AAG study
dim(HSP90)# 492   2
head(HSP90)# 17AAG response data (activity area)

## Assign resistant and sensitive cells
dim(HSP90)
set.seed(0)
res <- HSP90[sample(order(HSP90[,2])[1:60],50),]# ramdomly select 50 cells from 60 most resistant
set.seed(8)
sen <- HSP90[sample(order(HSP90[,2])[433:492],50),]# ramdomly select 50 cells from 60 most sensitive
res_idx <- match(res$Primary.Cell.Line.Name,colnames(CCLE_Exp_minus1))
sen_idx <- match(sen$Primary.Cell.Line.Name,colnames(CCLE_Exp_minus1))
library(rafalib)
mypar(1,2)
plot(sort(HSP90[,2]))
plot(c(res[,2],sen[,2]))

## Extract microarray data for resistant and sensitive cells and combine into a matrix
res_exp <- CCLE_Exp_minus1[,res_idx]
sen_exp <- CCLE_Exp_minus1[,sen_idx]
res_exp[1:5,1:5]
sen_exp[1:5,1:5]
X <- cbind(res_exp,sen_exp)
X[1:5,1:5]
dim(X)

## Wilcoxon Sum rank test
dim(res_exp)
dim(sen_exp)
pvals <- vector("numeric",nrow(sen_exp))
for(i in 1:nrow(sen_exp)) {
  pvals[i] <- wilcox.test(sen_exp[i,],res_exp[i,],alternative ="two.sided")$p.value
}
mypar(1,1)
hist(pvals,breaks=50)# anti-conservative
plot(pvals)
fdrs <- p.adjust(pvals, method="fdr")
hist(fdrs)
plot(fdrs)
sum(fdrs<0.25)# 1222 genes

## See how gene number changes with different cutoffs
cutoffs <- seq(0,0.5,0.05)
gene_len <- vector("numeric",length(cutoffs))
for(i in 1:length(cutoffs)){
  gene_len[i] <- length(which(fdrs<cutoffs[i]))
}
plot(cutoffs,gene_len)

## Save top 30 gene symbol
library(hgu133plus2.db)
gene_name <- select(hgu133plus2.db, keys=rownames(X), columns=c("SYMBOL"), keytype="ENTREZID")# Entrez ID to symbol
top30 <- gene_name[order(fdrs)[1:30],]
top30
# write.table(top30$SYMBOL, file="top30_wilcox.txt", quote=F, row.names=F, col.names=F)

## For the rest of the significant genes, reduce the number of correlated features by clustering them
library(apcluster)
X0 <- X[order(fdrs)[31:length(which(fdrs<0.25))],]
dim(X0)
apres <- apcluster(negDistMat(r=2), X0, details=TRUE)
apres # 22 cluster representatives
apres2 <- apcluster(negDistMat(r=2), X0, q=0.1, details=TRUE)
apres2 # 9 cluster representatives
list <- c(top30$ENTREZID, names(apres@exemplars))

## Save the "cluster representatives"
exemplars1<-select(hgu133plus2.db, keys=names(apres@exemplars), columns=c("SYMBOL"), keytype="ENTREZID")
# write.table(exemplars1$SYMBOL, file="exemplars_q0.5_wilcox.txt", quote=F, row.names=F, col.names=F)
exemplars2<-select(hgu133plus2.db, keys=names(apres2@exemplars), columns=c("SYMBOL"), keytype="ENTREZID")
# write.table(exemplars2$SYMBOL, file="exemplars_q0.1_wilcox.txt", quote=F, row.names=F, col.names=F)

## Naive Bayes classification
library(e1071)
library(stringr)
library(hgu133plus2.db)
library(caret)
dim(X)
Y <- as.factor(c(rep(1, ncol(res_exp)),rep(2, ncol(sen_exp))))# 1:res, 2:sen
X1 <- X[list,]
dim(X1)
X1[1:10,1:5]
X1 <- t(X1)
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

## ROC analysis
library(pROC)
Y_pred <- rowMeans(CV1)
roc1 <- roc(Y, Y_pred, plot=TRUE, ci=TRUE)
auc(roc1)# Area under the curve: 0.9168
ci(roc1)# 95% CI: 0.8591-0.9745 (DeLong)
coords(roc1,"b", ret=c("threshold", "specificity", "sensitivity", "ppv", "npv"))
coords(roc1,"l", ret=c("threshold", "specificity", "sensitivity", "ppv", "npv"))
# threshold for sen > 0.9994029
# threshold for res < 1.609626e-05

## Predict all cell lines (including sen and res)
X2 <- rearranged[list,]
dim(X2)
X2<-t(X2)
dim(X2)
X1[1:5,1:5]
X2[1:5,1:5]
identical(colnames(X2),colnames(X1))
model2 <- naiveBayes(X1, Y)
pred2 <- predict(model2, X2, type = "raw")
head(pred2)# 1:res, 2:sen

## Use probability thresholds suggested by the ROC analysis to define sen, res and mid
sen1 <- pred2[,2] > 0.9994029
res1 <- pred2[,2] < 1.609626e-05
mid1 <- pred2[,2] <= 0.9994029 & pred2[,2] >= 1.609626e-05

## Visualize the prediction result
mypar(1,2)
plot(pred2[,2], col=(as.numeric(sen1)+1))
plot(pred2[,2], col=(as.numeric(res1)+1))
# most cells have probabilities very close to 0 or 1 (not ideal..?)

## ROC analysis
mypar(1,1)
roc_sen1 <- roc(as.numeric(sen1), HSP90[,2], plot=TRUE, ci=TRUE)
auc(roc_sen1)# Area under the curve: 0.7399
ci(roc_sen1)# 95% CI: 0.6963-0.7836 (DeLong)
coords(roc_sen1,"b", ret=c("threshold", "specificity", "sensitivity", "ppv", "npv"))
# threshold specificity sensitivity         ppv         npv 
# 3.2585500   0.5986842   0.7819149   0.5464684   0.8161435 
roc_res1 <- roc(as.numeric(res1), HSP90[,2], plot=TRUE, ci=TRUE)
auc(roc_res1)# Area under the curve: 0.8
ci(roc_res1)# 95% CI: 0.7507-0.8493 (DeLong)
coords(roc_res1,"b", ret=c("threshold", "specificity", "sensitivity", "ppv", "npv"))
# threshold specificity sensitivity         ppv         npv 
# 2.8583500   0.7896104   0.7009346   0.4807692   0.9047619 
mypar(1,2)
plot(HSP90[,2], col=(as.numeric(sen1)+1), ylab="17AAG activity area")
abline(3.2585500,0)
plot(HSP90[,2], col=(as.numeric(res1)+1), ylab="17AAG activity area")
abline(2.8583500,0)

## t-test
mypar(1,1)
plot(density(HSP90[,2]))
lines(density(HSP90[sen1,2]), col=2)
lines(density(HSP90[res1,2]), col=3)
lines(density(HSP90[mid1,2]), col=4)
t.test(HSP90[sen1,2],HSP90[res1,2])$p.value# sen vs res
# [1] 9.241815e-28
t.test(HSP90[sen1,2],HSP90[mid1,2])$p.value# sen vs mid
# [1] 3.297282e-10
t.test(HSP90[res1,2],HSP90[mid1,2])$p.value# res vs mid
# [1] 3.348326e-12

## Predict all cell lines (excluding sen and res)
(res_idx2 <- match(res$Primary.Cell.Line.Name, colnames(rearranged)))
(match(res$Primary.Cell.Line.Name, HSP90$Primary.Cell.Line.Name))# confirm match
(sen_idx2 <- match(sen$Primary.Cell.Line.Name, colnames(rearranged)))
(match(sen$Primary.Cell.Line.Name, HSP90$Primary.Cell.Line.Name))
X3 <- rearranged[list,-c(res_idx2,sen_idx2)]
HSP90_2 <- HSP90[-c(res_idx2,sen_idx2),]
mypar(1,1)
plot(sort(HSP90[,2]))
points(sort(HSP90_2[,2]),col=1)
X3[1:5,1:5]
dim(X3)
dim(HSP90_2)
X3 <- t(X3)
dim(X3)
identical(colnames(X3),colnames(X1))
pred3 <- predict(model2, X3, type = "raw")
head(pred3)# 1:res, 2:sen

## Use thresholds suggested by the ROC analysis to define sen, res and mid
sen2 <- pred3[,2] > 0.9994029
res2 <- pred3[,2] < 1.609626e-05
mid2 <- pred3[,2] <= 0.9994029 & pred3[,2] >= 1.609626e-05

## Visualize the prediction result
mypar(1,2)
plot(pred3[,2], col=(as.numeric(sen2)+1))
plot(pred3[,2], col=(as.numeric(res2)+1))
# most cells have probabilities very close to 0 or 1 (not ideal..?)
mypar(1,3)
plot(HSP90_2[,2], col=(as.numeric(sen2)+1), ylab="17AAG activity area")
plot(HSP90_2[,2], col=(as.numeric(res2)+1), ylab="17AAG activity area")
plot(HSP90_2[,2], col=(as.numeric(mid2)+1), ylab="17AAG activity area")

## ROC analysis
mypar(1,1)
roc_sen2 <- roc(as.numeric(sen2), HSP90_2[,2], plot=TRUE, ci=TRUE)
auc(roc_sen2)# Area under the curve: 0.646
ci(roc_sen2)# 95% CI: 0.5907-0.7013 (DeLong)
coords(roc_sen2,"b", ret=c("threshold", "specificity", "sensitivity", "ppv", "npv"))
# threshold specificity sensitivity         ppv         npv 
# 3.2585500   0.5450820   0.7297297   0.4931507   0.7687861 
roc_res2 <- roc(as.numeric(res2), HSP90_2[,2], plot=TRUE, ci=TRUE)
auc(roc_res2)# Area under the curve: 0.6996
ci(roc_res2)# 95% CI: 0.6282-0.7711 (DeLong)
coords(roc_res2,"b", ret=c("threshold", "specificity", "sensitivity", "ppv", "npv"))
# threshold specificity sensitivity         ppv         npv 
# 3.1345000   0.6843750   0.6666667   0.3221477   0.9012346
mypar(1,2)
plot(HSP90_2[,2], col=(as.numeric(sen2)+1), ylab="17AAG activity area")
abline(3.2585500,0)
plot(HSP90_2[,2], col=(as.numeric(res2)+1), ylab="17AAG activity area")
abline(3.1345000,0)

## t-test
mypar(1,1)
plot(density(HSP90_2[,2]))
lines(density(HSP90_2[sen2,2]), col=2)
lines(density(HSP90_2[res2,2]), col=3)
lines(density(HSP90_2[mid2,2]), col=4)
t.test(HSP90_2[sen2,2],HSP90_2[res2,2])$p.value# sen vs res
# [1] 1.506088e-08
t.test(HSP90_2[sen2,2],HSP90_2[mid2,2])$p.value# sen vs mid
# [1] 0.001372112
t.test(HSP90_2[res2,2],HSP90_2[mid2,2])$p.value# res vs mid
# [1] 0.0002651682
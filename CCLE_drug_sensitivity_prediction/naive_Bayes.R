### Reconstruction of a procedure predicting sensitivity using naive Bayes classifier
### described in the CCLE drug sensitivity paper (Barretina et al, Nature, 2012)
### Data downloaded from CCLE website (https://portals.broadinstitute.org/ccle/home)

## Run code in "CCLE_data.R" to load and preprocess the data objects
dim(rearranged)# expression data for the drug of interest
rearranged[1:5,1:5]
dim(drug_data)# drug response data (activity area)
head(drug_data)

## Assign resistant and sensitive cells
dim(drug_data)
res <- drug_data[order(drug_data[,2])[1:50],]
sen <- drug_data[order(drug_data[,2])[(nrow(drug_data)-49):nrow(drug_data)],]
library(rafalib)
mypar(1,2)
plot(sort(drug_data[,2]))
plot(sort(c(res[,2],sen[,2])))
res_idx <- match(res$Primary.Cell.Line.Name, colnames(CCLE_Exp_minus1))
sen_idx <- match(sen$Primary.Cell.Line.Name, colnames(CCLE_Exp_minus1))

## Extract microarray data for resistant and sensitive cells and combine into a matrix
res_exp <- CCLE_Exp_minus1[,res_idx]
sen_exp <- CCLE_Exp_minus1[,sen_idx]
res_exp[1:5,1:5]
sen_exp[1:5,1:5]
X <- cbind(res_exp,sen_exp)
X[1:5,1:5]
dim(X)

## Naive Bayes classification with 10-fold cross validation
dim(sen_exp)
dim(res_exp)
Y <- as.factor(c(rep(1, ncol(res_exp)),rep(2, ncol(sen_exp))))# 1:res, 2:sen
length(Y)
dim(X)
X <- t(X)
dim(X)
CV <- matrix(data=NA, nrow=length(Y), ncol=5)
library(hgu133plus2.db)
gene_name <- select(hgu133plus2.db, keys=colnames(X), columns=c("SYMBOL"), keytype="ENTREZID")
head(gene_name)
library(caret)
library(apcluster)
library(e1071)
for(j in 1:5){
  idx <- createFolds(Y, k=10)
  for(i in 1:10){
    pvals <- vector("numeric", nrow(sen_exp))
    for(k in 1:nrow(sen_exp)) {
      pvals[k] <- wilcox.test(sen_exp[k,-idx[[i]][idx[[i]]<50]],
                              res_exp[k,-(idx[[i]][idx[[i]]>50]-50)],
                              alternative ="two.sided")$p.value
    }
    fdrs <- p.adjust(pvals, method="fdr")
    top30 <- gene_name[order(fdrs)[1:30],]# select top 30 significant genes
    m <- t(X[-idx[[i]],order(fdrs)[31:length(which(fdrs<0.25))]])# cluster the rest
    apres <- apcluster(negDistMat(r=2), m, details=TRUE)
    features<- c(top30$ENTREZID, names(apres@exemplars))
    model <- naiveBayes(X[-idx[[i]],features], Y[-idx[[i]]])
    pred <- predict(model, X[idx[[i]],features], type = "raw")
    CV[idx[[i]],j]<-pred[,2]
  }
}

## ROC analysis
library(pROC)
Y_pred <- rowMeans(CV)
roc1 <- roc(Y, Y_pred, plot=TRUE, ci=TRUE)
auc(roc1)# Area under the curve
ci(roc1)# 95% CI
coords(roc1,"b", ret=c("threshold", "specificity", "sensitivity", "ppv", "npv"))
roc1_tab <- coords(roc1,"l", ret=c("threshold", "specificity", "sensitivity", "ppv", "npv"))
sen_local <- min(which(roc1_tab["specificity",]>=0.98))
res_local <- max(which(roc1_tab["sensitivity",]>=0.98))
sen_thresh <- roc1_tab["threshold",sen_local]
res_thresh <- roc1_tab["threshold",res_local]
plot(Y_pred)
boxplot(t(CV))

## Predict category of all cell lines (excluding training data) 
## Re-assign resistant and sensitive cells for training set
dim(drug_data)
res2 <- drug_data[sample(order(drug_data[,2])[1:60],50),]# ramdomly select 50 cells from 60 most resistant
# do the same with 60 most resistant below
sen2 <- drug_data[sample(order(drug_data[,2])[(nrow(drug_data)-59):nrow(drug_data)],50),]
res_idx2 <- match(res2$Primary.Cell.Line.Name,colnames(CCLE_Exp_minus1))
sen_idx2 <- match(sen2$Primary.Cell.Line.Name,colnames(CCLE_Exp_minus1))
library(rafalib)
mypar(1,2)
plot(sort(drug_data[,2]))
plot(c(res2[,2],sen2[,2]))

## Extract microarray data for resistant and sensitive cells and combine into a matrix
res_exp2 <- CCLE_Exp_minus1[,res_idx2]
sen_exp2 <- CCLE_Exp_minus1[,sen_idx2]
X2 <- cbind(res_exp2,sen_exp2)
X2[1:5,1:5]
dim(X2)
X2 <- t(X2)

## Select features using naive Bayes on training set with 10-fold cross validation
idx <- createFolds(Y2, k=10)
features <- matrix(data = 0, nrow = nrow(sen_exp2), ncol = 10)
row.names(features) <- row.names(sen_exp2)
for(i in 1:10){
  pvals <- vector("numeric", nrow(sen_exp2))
  for(k in 1:nrow(sen_exp2)) {
    pvals[k] <- wilcox.test(sen_exp2[k,-idx[[i]][idx[[i]]<50]],
                            res_exp2[k,-(idx[[i]][idx[[i]]>50]-50)],
                            alternative ="two.sided")$p.value
  }
  fdrs <- p.adjust(pvals, method="fdr")
  top30 <- gene_name[order(fdrs)[1:30],]# select top 30 significant genes
  m <- t(X[-idx[[i]],order(fdrs)[31:length(which(fdrs<0.25))]])# cluster the rest
  apres <- apcluster(negDistMat(r=2), m, details=TRUE)
  features[c(top30$ENTREZID, names(apres@exemplars)),i] <- 1
}
features2 <- rownames(features)[rowMeans(features)>0.8]

## Naive Bayes classification on all cell lines (test set)
(res_idx3 <- match(res2$Primary.Cell.Line.Name, colnames(rearranged)))
(match(res2$Primary.Cell.Line.Name, drug_data$Primary.Cell.Line.Name))# confirm match
(sen_idx3 <- match(sen2$Primary.Cell.Line.Name, colnames(rearranged)))
(match(sen2$Primary.Cell.Line.Name, drug_data$Primary.Cell.Line.Name))
X3 <- rearranged[features2,-c(res_idx3,sen_idx3)]
drug_data_2 <- drug_data[-c(res_idx3,sen_idx3),]
mypar(1,1)
plot(sort(drug_data[,2]))
points(sort(drug_data_2[,2]),col=1)
X3[1:5,1:5]
dim(X3)
dim(drug_data_2)
X3 <- t(X3)
dim(X3)
identical(colnames(X2[,features2]),colnames(X3))
Y2 <- as.factor(c(rep(1, ncol(res_exp2)),rep(2, ncol(sen_exp2))))# 1:res, 2:sen
model2 <- naiveBayes(X2[,features2], Y2)
pred2 <- predict(model2, X3, type = "raw")
head(pred2)# 1:res, 2:sen
# Use probability thresholds suggested by the ROC analysis to define sen, res and mid
sen_pred2 <- pred2[,2] > sen_thresh
res_pred2 <- pred2[,2] < res_thresh
mid_pred2 <- (pred2[,2] <= sen_thresh) & (pred2[,2] >= res_thresh)
# Visualize the prediction result
mypar(1,2)
plot(pred2[,2], col=(as.numeric(sen_pred2)+1))
plot(pred2[,2], col=(as.numeric(res_pred2)+1))
mypar(1,3)
plot(drug_data_2[,2], col=(as.numeric(sen_pred2)+1), ylab="activity area")
plot(drug_data_2[,2], col=(as.numeric(res_pred2)+1), ylab="activity area")
plot(drug_data_2[,2], col=(as.numeric(mid_pred2)+1), ylab="activity area")
# ROC analysis
mypar(1,1)
roc_sen1 <- roc(as.numeric(sen_pred2), drug_data_2[,2], plot=TRUE, ci=TRUE)
auc(roc_sen1)# Area under the curve
ci(roc_sen1)# 95% CI
sen_pred2_thresh <- coords(roc_sen1,"b", ret = "threshold")
roc_res1 <- roc(as.numeric(res_pred2), drug_data_2[,2], plot=TRUE, ci=TRUE)
auc(roc_res1)# Area under the curve
ci(roc_res1)# 95% CI
res_pred2_thresh <- coords(roc_res1,"b", ret = "threshold")
mypar(1,2)
plot(drug_data_2[,2], col=(as.numeric(sen_pred2)+1), ylab="activity area")
abline(sen_pred2_thresh,0)
plot(drug_data_2[,2], col=(as.numeric(res_pred2)+1), ylab="activity area")
abline(res_pred2_thresh,0)
# t-test
mypar(1,1)
plot(density(drug_data_2[,2]))
lines(density(drug_data_2[sen_pred2,2]), col=2)
lines(density(drug_data_2[res_pred2,2]), col=3)
lines(density(drug_data_2[mid_pred2,2]), col=4)
t.test(drug_data_2[sen_pred2,2],drug_data_2[res_pred2,2])$p.value# sen vs res
t.test(drug_data_2[sen_pred2,2],drug_data_2[mid_pred2,2])$p.value# sen vs mid
t.test(drug_data_2[res_pred2,2],drug_data_2[mid_pred2,2])$p.value# res vs mid
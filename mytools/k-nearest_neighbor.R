### Predicting drug sensitivity using k-nearest neighbors (kNN)
### Method is adapted from EdX course PH525x series - Biomedical Data Science
### http://genomicsclass.github.io/book/pages/machine_learning.html
### Data downloaded from CCLE website (https://portals.broadinstitute.org/ccle/home)

## Run code in "CCLE_drug_sensitivity_prediction\CCLE_data.R" to load and preprocess the data objects
dim(rearranged)# expression data for the drug of interest
rearranged[1:5,1:5]
dim(drug_data)# drug response data (activity area)
head(drug_data)

## Assign resistant and sensitive cells
# Select top 60 resistant or sensitive cells
dim(drug_data)
res <- drug_data[order(drug_data[,2])[1:60],]
sen <- drug_data[order(drug_data[,2])[(nrow(drug_data)-59):nrow(drug_data)],]
res_idx <- match(res$Primary.Cell.Line.Name,colnames(CCLE_Exp_minus1))
sen_idx <- match(sen$Primary.Cell.Line.Name,colnames(CCLE_Exp_minus1))
library(rafalib)
mypar(1,2)
plot(sort(drug_data[,2]))
plot(c(res[,2],sen[,2]))

## Extract microarray data for resistant and sensitive cells and combine into a matrix
res_exp <- CCLE_Exp_minus1[,res_idx]
sen_exp <- CCLE_Exp_minus1[,sen_idx]
X <- cbind(res_exp,sen_exp)
X[1:5,1:5]
dim(X)

## Devide X into training set and test set
idx_train <- c(sample(1:60,48),sample(61:120,48))
X_train <- X[,idx_train]
X_test <- X[,-idx_train]
dim(X_train)
dim(X_test)
Y_train <- as.factor(c(rep("res", 48),rep("sen", 48)))
Y_test <- as.factor(c(rep("res", 12),rep("sen", 12)))

## Find best value for tunable parameter k using training set with 10-fold CV
train_error <- matrix(data = NA, nrow = 48, ncol = 10)
test_error <- matrix(data = NA, nrow = 48, ncol = 10)
library(caret)
idx <- createFolds(Y_train, k=10)
library(limma)
library(class)
for(k in 1:48){
  for(i in 1:10){
    mod <- model.matrix(~ Y_train[-idx[[i]]])
    fit <- lmFit(X_train[,-idx[[i]]], design = mod)
    fit <- eBayes(fit)
    tt <- topTable(fit, coef=2, number=1000, p.value = 0.05)
    tt <- tt[abs(tt$logFC) > 0.2,]
    ##predict on train
    yhat <- knn(t(X_train[row.names(tt),-idx[[i]]]),t(X_train[row.names(tt),-idx[[i]]]),Y_train[-idx[[i]]],k=k)
    train_error[k,i] <- 1-mean(as.numeric(yhat)==as.numeric(Y_train[-idx[[i]]]))
    ##make plot
    yhat <- knn(t(X_train[row.names(tt),-idx[[i]]]),t(X_train[row.names(tt),idx[[i]]]),Y_train[-idx[[i]]],k=k)
    test_error[k,i] <- 1-mean(as.numeric(yhat)==as.numeric(Y_train[idx[[i]]]))
  }
}
plot(rowMeans(train_error))
plot(rowMeans(test_error))
min(rowMeans(test_error))

## Identify differentially expressed genes in the training set and use signature to
## classify samples in the test set using the value for k that minimizes the error rate
mod <- model.matrix(~ Y_train)
fit <- lmFit(X_train, design = mod)
fit <- eBayes(fit)
tt <- topTable(fit, coef=2, number=1000, p.value = 0.05)
tt <- tt[abs(tt$logFC) > 0.2,]
k = which.min(rowMeans(test_error))
pred <- knn(t(X_train[row.names(tt),]),t(X_test[row.names(tt),]),Y_train,k=k)
1-mean(as.numeric(pred)==as.numeric(Y_test))# Error rate
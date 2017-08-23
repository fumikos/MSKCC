### Clustering samples and predicting drug sensitivity using basic machine learning methods
### Method is adapted from EdX course PH525x series - Biomedical Data Science
### http://genomicsclass.github.io/book/pages/clustering_and_heatmaps.html
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

## Differential gene expression analysis and classification with 10-fold cross validation
## using the training set to compare 2 different clustering methods
res_train <- X_train[, 1:48]
sen_train <- X_train[, 49:96]
dim(res_train)
dim(sen_train)
error_rates <- matrix(data = NA, nrow = 10, ncol = 2)
colnames(error_rates) <- c("hclust", "kmeans")
?matrix
library(hgu133plus2.db)
gene_name <- select(hgu133plus2.db, keys=rownames(X), columns=c("SYMBOL"), keytype="ENTREZID")
head(gene_name)
library(caret)
library(limma)
idx <- createFolds(Y_train, k=10)
for(i in 1:10){
  mod <- model.matrix(~ Y_train[-idx[[i]]])
  fit <- lmFit(X_train[,-idx[[i]]], design = mod)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef=2, number=1000, p.value = 0.05)
  tt <- tt[abs(tt$logFC) > 0.2,]
  d <- dist(t(X_train[row.names(tt),idx[[i]]]))
  hc <- hclust(d)
  hclusters <- cutree(hc, k=2)
  tab <- table(true=Y_train[idx[[i]]], cluster=hclusters)
  error_rates[i, 1] <- min((tab[1,1]+tab[2,2])/sum(tab),(tab[1,2]+tab[2,1])/sum(tab))
  km <- kmeans(t(X_train[row.names(tt),idx[[i]]]), centers=2)
  tab <- table(true=Y_train[idx[[i]]],cluster=km$cluster)
  error_rates[i, 2] <- min((tab[1,1]+tab[2,2])/sum(tab),(tab[1,2]+tab[2,1])/sum(tab))
  }
colMeans(error_rates)# Select clustering method that gives lower mean error rate

## Identify differentially expressed genes in the training set and use signature to
## classify samples in the test set
mod <- model.matrix(~ Y_train)
fit <- lmFit(X_train, design = mod)
fit <- eBayes(fit)
tt <- topTable(fit, coef=2, number=1000, p.value = 0.05)
tt <- tt[abs(tt$logFC) > 0.2,]
dim(tt)
# Using hclust
d <- dist(t(X_test[row.names(tt),]))
hc <- hclust(d)
hclusters <- cutree(hc, k=2)
tab <- table(true=Y_test, cluster=hclusters)
min((tab[1,1]+tab[2,2])/sum(tab),(tab[1,2]+tab[2,1])/sum(tab))# error rate
# Using kmeans
km <- kmeans(t(X_test), centers=2)
tab <- table(true=Y_test,cluster=km$cluster)
min((tab[1,1]+tab[2,2])/sum(tab),(tab[1,2]+tab[2,1])/sum(tab))# error rate
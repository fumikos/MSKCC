### Reconstruction of an elastic net regression combined with bootstrapping procedure
### described in the CCLE drug sensitivity paper (Barretina et al, Nature, 2012)
### Data downloaded from CCLE website (https://portals.broadinstitute.org/ccle/home)

## Run code in "CCLE_data.R" to load and preprocess the data objects
dim(rearranged)# expression data for the drug of interest
rearranged[1:5,1:5]
dim(drug_data)# drug response data (activity area)
head(drug_data)

## Input to the algorithm
Y <- drug_data[,2]
plot(sort(Y))
X <- rearranged
X[1:5,1:5]
X <- scale(t(X))
X[1:5,1:5]
# Use only features in X that correlates with Y with |R| > 0.1
R <- cor(X, Y, method="pearson")
hist(R, breaks=60)
idx_by_R <- which(R > 0.1 | R < -0.1)
length(idx_by_R)
X <- X[,idx_by_R]

## Parameter optimization with 10-fold cross validation
grid <- exp(seq(5,-6,length=250))
fold_id <- sample(rep(seq(10),length=length(Y)))
table(fold_id)
alphaslist <- seq(0.2, 1, length=10)
alphaslist
library(parallel)
detectCores()# 4 cores
library(doParallel)
registerDoParallel(4)
library(glmnet)
elasticnet <- lapply(alphaslist, function(a){
  cv.glmnet(X, Y, alpha = a, foldid = fold_id, lambda = grid, parallel = TRUE)
})
library(rafalib)
mypar(3,4)
for (i in 1:10) {plot(elasticnet[[i]])}
alphas <- 0# initiate alphas
for(i in 1:10) {alphas[i] <- min(elasticnet[[i]]$cvm)}
a <- alphaslist[which.min(alphas)]
l <- elasticnet[[which.min(alphas)]]$lambda.min

## Identification of predictive features with bootstrapping 200 times
boot_betas <- matrix(data=NA, nrow=ncol(X), ncol=200)
X_resample <- matrix(data=NA, nrow=length(Y), ncol=ncol(X))
Y_resample <- matrix(data=NA, nrow=length(Y), ncol=1)
for(i in 1:200){
  idx <- sample(length(Y), replace = T)
  for(j in 1:length(Y)){
    X_resample[j,] <- X[idx[j],]
    Y_resample[j] <- Y[idx[j]]
  }
  fit <- glmnet(X_resample, Y_resample, alpha = a, lambda = grid)
  boot_betas[,i] <- as.vector(coef(fit, s = l, exact = FALSE))[-1]# remove first row "(Intercept)"
}
head(boot_betas)
dim(boot_betas)
sum(rowSums(boot_betas != 0)/200 > 0.8)
top_feats <- rowSums(boot_betas != 0)/200 > 0.8
k <- colnames(X)[top_feats]
library(hgu133plus2.db)
gene_table <- select(hgu133plus2.db, keys=k, columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
betas_weight <- rowSums(boot_betas != 0)/200
betas_mean <- rowMeans(boot_betas, na.rm = TRUE)
gene_table$hits_BS <- betas_weight[top_feats]
gene_table$beta_mean <- betas_mean[top_feats]
gene_table <- gene_table[order(gene_table[,"beta_mean"]),]
gene_table

## Heatmap
order_gene <- match(gene_table$ENTREZID, rownames(rearranged))
order_gene
exp <- rearranged[order_gene,order(Y)]
exp <- t(scale(t(exp)))
gene_table$SYMBOL[is.na(gene_table$SYMBOL)] <- "<NA>"
rownames(exp) <- gene_table$SYMBOL
exp[1:5,1:5]
dim(exp)
mypar(1,1)
cols<-colorRampPalette(c("green","black","red"))
brks<-unique(c(seq(-4,4,length=100)))
top50s_idx <- c(1:50, (ncol(exp)-50):ncol(exp))
library(gplots)
heatmap.2(exp[,top50s_idx], Rowv = FALSE, Colv = FALSE, dendrogram = "none",
          col = cols, breaks = brks, trace = "none")# plot all cells
top_50s <- Y[order(Y)]
top_50s <- top_50s[top50s_idx]
plot(top_50s, ylab = "activity area")# curve indicating drug respose
barplot(gene_table$beta_mean)# barplot indicating weight of top predictive features

## Evaluating prediction performance by cross validation
library(caret)
CV_B <- matrix(data = NA, nrow = length(Y), ncol = 10)
for(j in 1:10){
  idx <- createFolds(Y, k = 10)
  for(i in 1:10){
    fit <- glmnet(X[-idx[[i]], ], Y[-idx[[i]]], alpha = a, lambda = grid)
    pred <- predict(fit, X[idx[[i]], ], type = "response", s = l)
    CV_B[idx[[i]],j] <- pred
  }
}
head(CV_B)
dim(CV_B)
Y_pred <- rowMeans(CV_B)
head(Y_pred)
length(Y_pred)
cor.test(Y, Y_pred, method = "pearson")
cor.test(Y, Y_pred, method = "kendall")
mypar(1,1)
plot(Y,Y_pred)
abline(0,1)

### Feature identification with elastic net with 10-fold cross validation for 100 times 
### described in the GCP paper (Garnett et al., Nature, 2012)
CV_G <- matrix(data = NA, nrow = ncol(X), ncol = 100)
for(j in 1:100){
  idx <- createFolds(Y, k = 10)
  for(i in 1:10){
    fit <- glmnet(X[-idx[[i]], ], Y[-idx[[i]]], alpha = a, lambda = grid)
    CV_G[,j] <- as.vector(coef(fit, s = l))[-1]# remove first row "(Intercept)"
  }
}
dim(CV_G)
head(CV_G)
feats <- rowSums(CV_G != 0)/100 != 0
sum(feats)
k2 <- colnames(X)[feats]
gene_table2 <- select(hgu133plus2.db, keys=k2, columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
head(gene_table2)
betas_weight2 <- rowSums(CV_G != 0)/100
betas_mean2 <- rowMeans(CV_G, na.rm = TRUE)
gene_table2$hits_BS <- betas_weight2[feats]
gene_table2$beta_mean <- betas_mean2[feats]
head(gene_table2)
M <- mean(gene_table2$beta_mean)
SD <- sd(gene_table2$beta_mean)
M2 <- mean(gene_table2$hits_BS)
SD2 <- sd(gene_table2$hits_BS)
sig_feats <- (gene_table2$beta_mean > (M + 2*SD) | gene_table2$beta_mean < (M - 2*SD)) & gene_table2$hits_BS > (M2 + 2*SD2)
sum(sig_feats)
gene_table2 <- gene_table2[sig_feats,]
gene_table2 <- gene_table2[order(gene_table2[,"beta_mean"]),]
gene_table2
# Heatmap
order_gene2 <- match(gene_table2$ENTREZID, rownames(rearranged))
order_gene2
exp2 <- rearranged[order_gene2,order(Y)]
exp2 <- t(scale(t(exp2)))
gene_table2$SYMBOL[is.na(gene_table2$SYMBOL)] <- "<NA>"
rownames(exp2) <- gene_table2$SYMBOL
exp2[1:5,1:5]
dim(exp2)
mypar(1,1)
cols<-colorRampPalette(c("green","black","red"))
brks<-unique(c(seq(-4,4,length=100)))
top50s_idx <- c(1:50, (ncol(exp)-50):ncol(exp))
library(gplots)
heatmap.2(exp2[,top50s_idx], Rowv = FALSE, Colv = FALSE, dendrogram = "none",
          col = cols, breaks = brks, trace = "none")# plot all cells
top_50s <- Y[order(Y)]
top_50s <- top_50s[top50s_idx]
plot(top_50s, ylab = "activity area")# curve indicating drug respose
barplot(gene_table2$beta_mean)# barplot indicating weight of top predictive features

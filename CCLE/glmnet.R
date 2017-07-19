### Sensitivity prediction using elastic net regression 

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

## Sensitivity prediction using elastic net regression analysis
Y <- HSP90[,2]
plot(sort(Y))
X <- rearranged
X[1:5,1:5]
X <- t(X)
X <- scale(X)
X[1:5,1:5]
# Use only features in X that correlates with Y with |R| > 0.1
R <- cor(X, Y, method="pearson")
hist(R, breaks=60)
idx_by_R <- which(R > 0.1 | R < -0.1)
length(idx_by_R)# 2857 features
X <- X[,idx_by_R]
# Parameter optimization 1: alpha
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
for (i in 1:10) {print(min(elasticnet[[i]]$cvm))}# optimal alpha = 0.2
a <- 0.2
for (i in 1:10) {print(elasticnet[[i]]$lambda.min)}
for (i in 1:10) {print(elasticnet[[i]]$lambda.1se)}
# Parameter optimization 2: lambda
registerDoParallel(4)
lambdas <- replicate(50, {
  elasticnet <- cv.glmnet(X, Y, alpha = a, nfolds = 10, lambda = grid, parallel = TRUE)
  elasticnet$lambda.min
  })
lambdas
sd(lambdas)
mol <- mean(lambdas)
mol# lambda = 0.07229556
# Identification of predictive features with bootstrapping
boot_betas <- matrix(data=NA, nrow=ncol(X), ncol=200)
X_resample <- matrix(data=NA, nrow=length(Y), ncol=ncol(X))
Y_resample <- matrix(data=NA, nrow=length(Y), ncol=1)
dim(boot_betas)# 2858  200
dim(X_resample)# 492 2857
dim(Y_resample)# 492   1
for(i in 1:200){
  idx <- sample(length(Y), replace = T)
  for(j in 1:length(Y)){
    X_resample[j,] <- X[idx[j],]
    Y_resample[j] <- Y[idx[j]]
    }
  fit <- glmnet(X_resample, Y_resample, alpha = a, lambda = grid)
  boot_betas[,i] <- as.vector(coef(fit, s = mol, exact = FALSE))[-1]# remove first row "(Intercept)"
}
head(boot_betas)
dim(boot_betas)
sum(rowSums(boot_betas != 0) > 160)# 20 features significant in > 80% of the bootstrap datasets
top_feats <- rowSums(boot_betas != 0) > 160
k <- colnames(X)[top_feats]
library(hgu133plus2.db)
gene_table <- select(hgu133plus2.db, keys=k, columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
head(gene_table)
betas_weight <- rowSums(boot_betas != 0)/200
gene_table$hits_BS <- betas_weight[top_feats]
betas_mean <- rowMeans(boot_betas, na.rm = TRUE)
gene_table$beta_mean <- betas_mean[top_feats]
head(gene_table)
gene_table <- gene_table[order(gene_table[,"beta_mean"]),]
gene_table
# write.csv(gene_table, file="All_vs_HSP90_BS.csv", row.names = FALSE)

## Heatmap
library(gplots)
order_gene <- match(gene_table$ENTREZID, rownames(rearranged))
order_gene
exp <- rearranged[order_gene,order(Y)]
exp <- t(scale(t(exp)))
rowMeans(exp)
rowSds(exp)
rownames(exp) <- gene_table$SYMBOL
exp[1:5,1:5]
dim(exp)
mypar(1,1)
cols<-colorRampPalette(c("green","black","red"))
brks<-unique(c(seq(-5,-1,length=100),seq(-1,1,length=50),seq(1,5,length=100)))
brks2<-unique(c(seq(-3.5,3.5,length=100)))
idx3<-51:442
length(idx3)# 392
heatmap.2(exp[,-idx3],Rowv = FALSE, Colv = FALSE, dendrogram = "none",
          col = cols, breaks = brks2, trace = "none")# plot top 50s only
heatmap.2(exp, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
          col = cols, breaks = brks2, trace = "none")# plot all cells
top_50s <- Y[order(Y)]
top_50s <- top_50s[-idx3]
plot(top_50s, ylab = "17AAG activity area")# curve indicating drug respose
barplot(gene_table$beta_mean)# barplot indicating weight of top predictive features

## Compare result with Pearson's correlation coefficient r
R <- cor(t(rearranged), Y, method = "pearson")
mypar(1,2)
plot(R)
R2 <- cor(X, Y, method = "pearson")
plot(R2[top_feats])

## Evaluating prediction performance by cross validation (Barretina et al.)
library(caret)
CV_B <- matrix(data = NA, nrow = length(Y), ncol = 10)
for(j in 1:10){
  idx <- createFolds(Y, k = 10)
  for(i in 1:10){
    fit <- glmnet(X[-idx[[i]], ], Y[-idx[[i]]], alpha = a, lambda = grid)
    pred <- predict(fit, X[idx[[i]], ], type = "response", s = mol)
    CV_B[idx[[i]],j] <- pred
    }
  }
head(CV_B)
dim(CV_B)
Y_pred <- rowMeans(CV_B)
head(Y_pred)
length(Y_pred)
cor(Y, Y_pred)
cor.test(Y, Y_pred, method = "pearson")
cor.test(Y, Y_pred, method = "kendall")
mypar(1,1)
plot(Y,Y_pred)
abline(0,1)

## Elastic net with cross validation (Garnett et al.)
CV_G <- matrix(data = NA, nrow = ncol(X), ncol = 100)
for(j in 1:100){
  idx <- createFolds(Y, k = 10)
  for(i in 1:10){
    fit <- glmnet(X[-idx[[i]], ], Y[-idx[[i]]], alpha = a, lambda = grid)
    CV_G[,j] <- as.vector(coef(fit, s = mol))[-1]# remove first row "(Intercept)"
    }
  }
dim(CV_G)
head(CV_G)
feats <- rowSums(CV_G != 0)/100 != 0
sum(feats)# 1449
k2 <- colnames(X)[feats]
gene_table2 <- select(hgu133plus2.db, keys=k2, columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
head(gene_table2)
betas_weight2 <- rowSums(CV_G != 0)/100
gene_table2$hits_BS <- betas_weight2[feats]
betas_mean2 <- rowMeans(CV_G, na.rm = TRUE)
gene_table2$beta_mean <- betas_mean2[feats]
head(gene_table)
M <- mean(gene_table2$beta_mean)#-0.0003600625
SD <- sd(gene_table2$beta_mean)# 0.01164365
M2 <- mean(gene_table2$hits_BS)# 0.273568
SD2 <- sd(gene_table2$hits_BS)# 0.3146805
sig_feats <- (gene_table2$beta_mean > M + 2*SD | gene_table2$beta_mean < M - 2*SD) & gene_table2$hits_BS > M2 + 2*SD2
sum(sig_feats)
gene_table2 <- gene_table2[sig_feats,]
gene_table2 <- gene_table2[order(gene_table2[,"beta_mean"]),]
gene_table2
# write.csv(gene_table2, file="All_vs_HSP90_CV.csv", row.names = FALSE)

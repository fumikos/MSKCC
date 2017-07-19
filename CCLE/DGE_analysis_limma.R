### Test for differentially expressed genes using limma

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
res <- HSP90[order(HSP90[,2])[1:50],]# top 50 resistant
sen <- HSP90[order(HSP90[,2])[443:492],]# top 50 sensitive
library(rafalib)
mypar(1,2)
plot(sort(HSP90[,2]))
plot(sort(c(res[,2],sen[,2])))
res_idx <- match(res$Primary.Cell.Line.Name, colnames(CCLE_Exp_minus1))
sen_idx <- match(sen$Primary.Cell.Line.Name, colnames(CCLE_Exp_minus1))
X <- CCLE_Exp_minus1[,c(res_idx,sen_idx)]
dim(X)
X[1:5,1:5]

## Exploratory data analysis (EDA) on X
## Have a look at distributions
mypar(1,1)
boxplot(X)
mypar(2,1)
hist(rowMeans(X), breaks=50, xlim = c(2,15))
library(genefilter)
hist(rowSds(X), breaks=50, xlim = c(0,3))

## Heatmap
mypar(1,1)
dim(X[rowMeans(X)>5 & rowSds(X)>1,])# 3175  100
heatmap(X[rowMeans(X)>5 & rowSds(X)>1,])

## Correlation matrix
library(RColorBrewer)
n <- ncol(X)
cors <- cor(X - rowMeans(X))
cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
mypar(1,1)
image(cors,xaxt="n",yaxt="n",col=cols,xlab="",ylab="",zlim=c(-1,1))

## Principal component analysis (PCA)
s0 <- svd(X - rowMeans(X))
PC1 <- s0$d[1]*s0$v[,1]
PC2 <- s0$d[2]*s0$v[,2]
mypar(1,1)
Y <- as.factor(c(rep("res", 50),rep("sen", 50)))
plot(PC1,PC2,pch=21,bg=Y)
legend("bottomright",levels(Y), col=seq(along=levels(Y)), pch=15, cex=0.8)
plot(s0$d^2/sum(s0$d^2), ylab="% variance explained",xlab="Principal component")

## Test for differentially expressed genes using limma
## http://genomicsclass.github.io/book/pages/bioc1_geneset_1.html
library(limma)
library(sva)
# Fit all genes without adjusting with surrogate variables
dim(X)
mod <- model.matrix(~ Y)
fit <- lmFit(X, design = mod)
fit <- eBayes(fit)# adjust variance estimates using a hierarchical model
dim(fit$p.value)
head(fit$p.value)
pvals <- fit$p.value[,2]
hist(pvals)# anti-conservative
qvals<-p.adjust(pvals, method = "BH")
hist(qvals)# conservative
fc <- fit$coefficients[,2]
plot(fc, -log10(pvals), xlab = "Difference in means", 
     col = as.factor(qvals < 0.25 & abs(fc) > 0.15))# too permissive
plot(fc, -log10(pvals), xlab = "Difference in means", 
     col = as.factor(qvals < 0.1 & abs(fc) > 0.15))# better
library(hgu133plus2.db)
symbols <- select(hgu133plus2.db, keys=row.names(fit), columns=c("SYMBOL"), keytype="ENTREZID")
symbols$SYMBOL[!(qvals < 0.1 & abs(fc) > 0.15)] <- ""
text(fc, -log10(pvals), labels=symbols$SYMBOL, cex = 0.5, pos = 3)
sum(qvals < 0.1 & abs(fc) > 0.15)# 530 significant genes
tt <- topTable(fit, coef=2, number=1000, p.value = 0.1)
dim(tt)
tt <- tt[abs(tt$logFC) > 0.15,]
dim(tt)# 530   6
tt <- tt[order(abs(tt$logFC),decreasing=TRUE),]
tt[1:10,]
gene_table <- select(hgu133plus2.db, keys=row.names(tt), columns=c("SYMBOL"), keytype="ENTREZID")
tt <- cbind(gene_table, tt)
tt[1:10,]
# write.table(tt, file="limma_tt.csv", quote=F, row.names=F)
# Fit all genes with adjustment using surrogate variables (SVs)
s <- sva(X, mod)
svamod <- model.matrix(~ Y + s$sv)# 18 significant SVs
svafit <- lmFit(X, design = svamod)
svafit <- eBayes(svafit)
svapvals <- svafit$p.value[,2]
hist(svapvals)# anti-conservative
svaqvals<-p.adjust(svapvals, method = "BH")
hist(svaqvals)# conservative
svafc <- svafit$coefficients[,2]
plot(svafc, -log10(svapvals), xlab = "Difference in means", 
     col = as.factor(svaqvals < 0.25 & abs(svafc) > 0.15))
svasymbols <- select(hgu133plus2.db, keys=row.names(svafit), columns=c("SYMBOL"), keytype="ENTREZID")
svasymbols$SYMBOL[!(svaqvals < 0.25 & abs(svafc) > 0.15)] <- ""
text(svafc,-log10(svapvals),labels=svasymbols$SYMBOL,cex=0.5,pos=3)
sum(svaqvals < 0.25 & abs(svafc) > 0.15)# 36 genes
svatt <- topTable(svafit, coef=2, number=100, p.value = 0.25)
dim(svatt)
svatt <- svatt[abs(svatt$logFC) > 0.15,]
dim(svatt)# 36  6
svatt <- svatt[order(abs(svatt$logFC),decreasing=TRUE),]
svatt[1:10,]
svagene_table <- select(hgu133plus2.db, keys=row.names(svatt), columns=c("SYMBOL"), keytype="ENTREZID")
svatt <- cbind(svagene_table, svatt)
svatt[1:10,]
# write.table(svatt, file="limma_svatt.csv", quote=F, row.names=F)
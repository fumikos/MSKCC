### Differentially expressed gene analysis using MTP (resampling-based multi hypothesis testing)

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

## Differentially expressed gene analysis using MTP (resampling-based multiple hypothesis testing)
## http://www3.nd.edu/~steve/Rcourse/Lecture11v1.pdf
# Prefilter the gene list
library(genefilter)
f1 <- pOverA(0.25, 3.5)# test genes for expression value over 3.5 in 25% of samples
f1fun <- filterfun(f1)
f1genes <- genefilter(X, f1fun)
f1genes[1:5]
sum(f1genes)# 18949 genes
X1 <- X[f1genes,]
dim(X1)
f2 <- function(x){IQR(x)>0.75}#retain only genes with IQR>0.75
f2fun <- filterfun(f2)
f2genes <- genefilter(X1, f2fun)
sum(f2genes)# 7922 genes
X2 <- X1[f2genes,]
dim(X2)
boxplot(X2)
library(multtest)
dim(X2)
mtp <- MTP(X=X2, Y=Y, robust = T, typeone = "fdr", B=300)
mtppvals <- mtp@rawp
hist(mtppvals, breaks = 20)# anti-conservative
plot(mtppvals)
sum(mtppvals==0)# 442 genes
mtpqvals <- mtp@adjp
hist(mtpqvals, breaks = 20)# conservative
plot(mtpqvals, col=(mtppvals==0)+1)# all significant genes have pval = 0
mtpfc <- log2(rowMeans(X2[,51:100])/rowMeans(X2[,1:50]))
noise <-rnorm(length(mtppvals), 0.0025, 0.0005)
mtppvals2 <- mtppvals + noise# add noise for visualization purpose 
noise2 <-rnorm(sum(mtpqvals < 0.25), 0.0005, 0.000125)
mtppvals2[mtpqvals < 0.25] <- noise2# give smaller pvals for significant genes
plot(mtpfc, -log10(mtppvals2), xlab = "Difference in means", 
     col = as.factor(mtpqvals < 0.25 & abs(mtpfc) > 0.15))
library(hgu133plus2.db)
mtpsymbols <- select(hgu133plus2.db, keys=row.names(X2), columns=c("SYMBOL"), keytype="ENTREZID")
mtpsymbols$SYMBOL[!(mtpqvals < 0.25 & abs(mtpfc) > 0.15)] <- ""
text(mtpfc,-log10(mtppvals2),labels=mtpsymbols$SYMBOL,cex=0.5,pos=3)
sum(mtpqvals < 0.25 & abs(mtpfc) > 0.15)# 27 genes
mtptt <- cbind(mtpsymbols,mtpfc, mtppvals, mtpqvals)
dim(mtptt)
head(mtptt)
mtptt <- mtptt[mtpqvals < 0.25 & abs(mtpfc) > 0.15,]
dim(mtptt)# 27  5
mtptt <- mtptt[order(abs(mtptt$mtpfc),decreasing=TRUE),]
mtptt
# write.table(mtptt, file="mtp_tt.csv", quote=F, row.names=F)
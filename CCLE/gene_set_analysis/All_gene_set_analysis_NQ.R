library(CePa)
library(stringr)
library(hgu133plus2.db)
library(rafalib)
library(genefilter)
library(RColorBrewer)
library(limma)
library(org.Hs.eg.db)
library(GO.db)
library(multtest)
library(sva)

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

library(CePa)
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
HSP90 = CCLE_Drug_all[HSP90_index,c(1,2,13)]
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

#Change rownames(probeID) in rearranged to gene symbol
row.names(rearranged) <- str_replace_all(row.names(rearranged),"_at","")
k<-as.character(row.names(rearranged))
gene_name <- select(hgu133plus2.db, keys=k, columns=c("SYMBOL"), keytype="ENTREZID")
row.names(rearranged)[1:10]
gene_name[1:10,]
row.names(rearranged)<-gene_name[,2]
rearranged[1:5,1:5]
dim(rearranged)

# Subset 'rearranged' based on NQO1 expression
mypar(1,1)
hist(rearranged["NQO1",], breaks = 50)
plot(rearranged["NQO1",order(HSP90[,3])])
abline(10,0)
NQ_high<-rearranged[,rearranged["NQO1",]>10]
dim(NQ_high)# 358 cells have NQO1 expression higher than 10
row.names(NQ_high)<-row.names(CCLE_Exp_minus1)# change back row names to probe ID
HSP90_NQ_high<-HSP90[rearranged["NQO1",]>10,]
dim(HSP90_NQ_high)
mypar(1,2)
plot(sort(HSP90_NQ_high[,3]))
plot(sort(HSP90[,3]))

#Assign resistant and sensitive cells
res<-HSP90_NQ_high[order(HSP90_NQ_high[,3])[1:35],]
sen<-HSP90_NQ_high[order(HSP90_NQ_high[,3])[324:358],]
res_idx<-match(res$CCLE.Cell.Line.Name,colnames(CCLE_Exp_minus1))
sen_idx<-match(sen$CCLE.Cell.Line.Name,colnames(CCLE_Exp_minus1))
mypar(1,2)
plot(sort(HSP90_NQ_high[,3]))
plot(c(res[,3],sen[,3]))
mypar(2,1)
table(CCLE_CELL_all$Site.Primary[res_idx])
table(CCLE_CELL_all$Site.Primary[sen_idx])
plot(CCLE_CELL_all$Site.Primary[res_idx])
plot(CCLE_CELL_all$Site.Primary[sen_idx])

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

#save X as cvs file for other use
#write.csv(X, "expression_100each_SYNBL.csv")

# Have a look at distribution of X

dim(X)
mypar(1,1)
boxplot(X)
mypar(2,1)
summary(colMeans(X))
hist(rowMeans(X), breaks=50,xlim = c(2,15))
library(genefilter)
summary(rowSds(X))
hist(rowSds(X), breaks=50,xlim = c(0,3))
sum(rowMeans(X)<4.1 & rowSds(X)<0.3)#2098 genes--not much impact
sum(rowMeans(X)<4.1)#3141 genes
sum(rowSds(X)<0.3)#3788 genes

# PCA

s0<-svd(X-rowMeans(X))
PC1 <- s0$d[1]*s0$v[,1]
PC2 <- s0$d[2]*s0$v[,2]
mypar(1,1)
Y<-as.factor(c(rep(0, 35),rep(1, 35)))
plot(PC1,PC2,pch=21,bg=Y)
legend("bottomright",levels(Y),col=seq(along=levels(Y)),pch=15,cex=0.5)
plot(s0$d^2/sum(s0$d^2),ylab="% variance explained",xlab="Principal component")
# Less than 20% of variation is explained by the first two PCs

# Correlation matrix

library(RColorBrewer)#

n <- ncol(X)
cors=cor(X-rowMeans(X))
cols=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
mypar()
image(cors,xaxt="n",yaxt="n",col=cols,xlab="",ylab="",zlim=c(-1,1))

# Prefilter the gene list(from http://www3.nd.edu/~steve/Rcourse/Lecture11v1.pdf)

library(genefilter)#

f1<-pOverA(0.25, 3.5)#test genes for expression value over 3.5 in 25% of samples
f1fun<-filterfun(f1)
f1genes<-genefilter(X, f1fun)
f1genes[1:5]
sum(f1genes)#18936 genes
X1<-X[f1genes,]
dim(X1)

f2<-function(x){IQR(x)>0.75}#retain only genes with IQR>0.75
f2fun<-filterfun(f2)
f2genes<-genefilter(X1, f2fun)
sum(f2genes)
X2<-X1[f2genes,]
dim(X2)#7399 genes

boxplot(X2)
mypar(2,1)
summary(rowMeans(X2))
hist(rowMeans(X2), breaks=50,xlim = c(2,15))
summary(rowSds(X2))
hist(rowSds(X2), breaks=50,xlim = c(0,3))

s2<-svd(X2-rowMeans(X2))
PC1 <- s2$d[1]*s2$v[,1]
PC2 <- s2$d[2]*s2$v[,2]
mypar(1,1)
plot(PC1,PC2,pch=21,bg=Y)
legend("bottomright",levels(Y),col=seq(along=levels(Y)),pch=15,cex=0.5)

# Test for differentially expressed genes using limma
# http://genomicsclass.github.io/book/pages/bioc1_geneset_1.html

library(limma)
library(sva)

# Fit all genes without adjusting with surrogate variables

dim(X)
mod<- model.matrix(~Y)
fit <- lmFit(X, design=mod)
fit <- eBayes(fit)
tt <- topTable(fit, coef=2, number=200)
tt_up<-tt[tt$logFC>0,]
tt_dn<-tt[tt$logFC<0,]
dim(tt_up)
dim(tt_dn)
head(tt_up)
head(tt_dn)

#plot(X["NQO1",],col=(X["NQO1",]<10)+1)

# Fit all genes adjusting with surrogate variables
#s <- sva(X,mod)# Takes a while to run
#svamod<-model.matrix(~Y+s$sv)# 18 significant SVs -- this data probably does not require correction
#svafit <- lmFit(X, design=svamod)
#svafit <- eBayes(svafit)
#svatt <- topTable(svafit, coef=2, number=50)
#svatt$ID

#Fit filtered genes without adjusting with surrogate variables
#dim(X2)
#fit2 <- lmFit(X2, design=mod)
#fit2 <- eBayes(fit2)
#tt2 <- topTable(fit2, coef=2, number=200)
#tt2_up<-tt2[tt2$logFC>0,]
#tt2_dn<-tt2[tt2$logFC<0,]
#dim(tt2_up)
#dim(tt2_dn)

#k2_up <- str_replace_all(row.names(tt2_up),"_at","")
#sen_sig2 <- select(hgu133plus2.db, keys=k2_up, columns=c("SYMBOL"), keytype="ENTREZID")
#k2_dn <- str_replace_all(row.names(tt2_dn),"_at","")
#res_sig2 <- select(hgu133plus2.db, keys=k2_dn, columns=c("SYMBOL"), keytype="ENTREZID")

#write.table(sen_sig2$SYMBOL, file="sen_sig_filtered_limma_NQ.txt", quote=F, row.names=F, col.names=F)
#write.table(res_sig2$SYMBOL, file="res_sig_filtered_limma_NQ.txt", quote=F, row.names=F, col.names=F)

# Fit filtered genes adjusting with surrogate variables
#s2 <- sva(X2,mod)# Takes a while to run
#svamod2<-model.matrix(~Y+s2$sv)# 16 significant SVs
#svafit2 <- lmFit(X, design=svamod)
#svafit2 <- eBayes(svafit2)
#svatt2 <- topTable(svafit2, coef=2, number=50)
#svatt2$ID

# Save gene symbol (from non-filtered result)

k <- str_replace_all(row.names(CCLE_Exp_minus1),"_at","")# probe ID to Entrez ID
gene_name <- select(hgu133plus2.db, keys=k, columns=c("SYMBOL"), keytype="ENTREZID")# Entrez ID to symbol
head(gene_name)
dim(gene_name)

sen_sig<-gene_name[fit$coefficients[,2]>0,]
dim(sen_sig)
sen_sig<-sen_sig[order(fit$p.value[fit$coefficients[,2]>0, 2])[1:100],]
head(sen_sig)
write.table(sen_sig$SYMBOL, file="sen_limma_NQ.txt", quote=F, row.names=F, col.names=F)

res_sig<-gene_name[fit$coefficients[,2]<0,]
dim(res_sig)
res_sig<-res_sig[order(fit$p.value[fit$coefficients[,2]<0, 2])[1:100],]
head(res_sig)
write.table(res_sig$SYMBOL, file="res_limma_NQ.txt", quote=F, row.names=F, col.names=F)

# Test for differentially expressed genes using MTP (permutated null distribution)
# http://www3.nd.edu/~steve/Rcourse/Lecture11v1.pdf

library(multtest)#

mtp<-MTP(X=X2, Y=Y, robust = T, typeone = "fdr", B=300)#takes a long time to run (~30 min)
print(mtp)
summary(mtp)
hist(mtp@rawp, breaks = 50)
hist(mtp@adjp, breaks = 50)

cutoffs<-seq(0,0.5,0.01)
gene_len<-vector("numeric",length(cutoffs))
for(i in 1:length(cutoffs)){
  gene_len[i]<-sum(mtp@adjp<cutoffs[i])
}
plot(cutoffs,gene_len)

k3 <- str_replace_all(names(mtp@adjp)[mtp@adjp<0.20],"_at","")
multi_sig <- select(hgu133plus2.db, keys=k3, columns=c("SYMBOL"), keytype="ENTREZID")

write.table(multi_sig$SYMBOL, file="multtest_NQ.txt", quote=F, row.names=F, col.names=F)

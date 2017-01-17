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
head(CCLE_Drug_all)
CCLE_Exp_all[1:5,1:5]
head(CCLE_CELL_all)

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

#Next, extract Activity Area data for HSP90 inhibitor in gliomas
glioma_HSP90_idx = which(grepl("CENTRAL_NERVOUS_SYSTEM",CCLE_Drug_all$CCLE.Cell.Line.Name) 
                                & CCLE_Drug_all$Target=="HSP90")
length(glioma_HSP90_idx)# 29 glioma cells have 17AAG data
glioma_HSP90 = CCLE_Drug_all[glioma_HSP90_idx,c(1,2,13)]
glioma_HSP90[1:5,]

#Next, match cell name in glioma_HSP90 and CCLE_Exp_minus1
glioma_Exp_idx = which(grepl("CENTRAL_NERVOUS_SYSTEM", colnames(CCLE_Exp_minus1)))
length(glioma_Exp_idx)#69 glioma cells have microarray data
sum(grepl("central_nervous_system",CCLE_CELL_all$Site.Primary))#confirm this in CCLE_CELL_all
glioma_Exp = CCLE_Exp_minus1[,glioma_Exp_idx]
colnames(glioma_Exp)
(glioma_names =CCLE_CELL_all$Cell.line.primary.name[glioma_Exp_idx])
colnames(glioma_Exp)<-glioma_names# Simplify the cell names
(idx = match(glioma_HSP90$Primary.Cell.Line.Name, colnames(glioma_Exp)))
length(which(is.na(idx)))#All 29 cells which have 17AAG data also have microarray data
rearranged = glioma_Exp[,idx]#extract expression data and reorder to match glioma_HSP90 data order
dim(rearranged)
colnames(rearranged)

#Assign resistant and sensitive cells
dim(glioma_HSP90)
HSP90_P1<-glioma_HSP90[order(glioma_HSP90[,3])[seq(1,29,2)],]# split data into 2 partitions
HSP90_P2<-glioma_HSP90[order(glioma_HSP90[,3])[seq(2,28,2)],]
dim(HSP90_P1)
dim(HSP90_P2)
res<-HSP90_P2[order(HSP90_P2[,3])[1:5],]
sen<-HSP90_P2[order(HSP90_P2[,3])[10:14],]
res_idx<-match(res$Primary.Cell.Line.Name,colnames(glioma_Exp))
sen_idx<-match(sen$Primary.Cell.Line.Name,colnames(glioma_Exp))
res_exp<-glioma_Exp[,res_idx]
sen_exp<-glioma_Exp[,sen_idx]
colnames(res_exp)
colnames(sen_exp)
mypar(1,2)
plot(sort(glioma_HSP90[,3]))
plot(c(res[,3],sen[,3]))

# PCA for resistant and sensitive glioma cells
#pheno<-vector("numeric",69)
#pheno[res_idx]<-rep(1,length(res_idx))
#pheno[sen_idx]<-rep(2,length(sen_idx))
#s <- svd(glioma_Exp-rowMeans(glioma_Exp))
#PC1 <- s$d[1]*s$v[,1]
#PC2 <- s$d[2]*s$v[,2]
#mypar(1,1)
#plot(PC1,PC2,pch=21,bg=pheno)
#legend("topleft",c("others","resistant","sensitive"),col=c(0,1,2),pch=15,cex=0.8)
#text(PC1,PC2,labels=glioma_names,cex=0.5,pos=3)

# Prefilter the gene list(from http://www3.nd.edu/~steve/Rcourse/Lecture11v1.pdf)

library(genefilter)#

X<-cbind(res_exp,sen_exp)
X[1:5,1:5]
dim(X)

f1<-pOverA(0.25, 3.5)#test genes for expression value over 3.5 in 25% of samples
f1fun<-filterfun(f1)
f1genes<-genefilter(X, f1fun)
f1genes[1:5]
sum(f1genes)#18866 genes
X1<-X[f1genes,]
dim(X1)

f2<-function(x){IQR(x)>0.75}#retain only genes with IQR>0.75
f2fun<-filterfun(f2)
f2genes<-genefilter(X1, f2fun)
sum(f2genes)
X2<-X1[f2genes,]
dim(X2)#4678 genes

boxplot(X2)
mypar(2,1)
summary(rowMeans(X2))
hist(rowMeans(X2), breaks=50,xlim = c(2,15))
summary(rowSds(X2))
hist(rowSds(X2), breaks=50,xlim = c(0,3))

# Correlation matrix

#library(RColorBrewer)#

#n <- ncol(X)
#cors=cor(X-rowMeans(X))
#cols=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
#mypar()
#image(cors,xaxt="n",yaxt="n",col=cols,xlab="",ylab="",zlim=c(-1,1))

#n2 <- ncol(X2)
#cors2=cor(X2-rowMeans(X2))
#image(cors2,xaxt="n",yaxt="n",col=cols,xlab="",ylab="",zlim=c(-1,1))

# Test for differentially expressed genes using limma
# http://genomicsclass.github.io/book/pages/bioc1_geneset_1.html

library(limma)
library(sva)

# Fit all genes
Y<-as.factor(c(rep(0, 5),rep(1, 5)))#0:res, 1:sen
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
tail(tt_up)
tail(tt_dn)

# Fit filtered genes
dim(X2)
fit2 <- lmFit(X2, design=mod)
fit2 <- eBayes(fit2)
tt2 <- topTable(fit2, coef=2, number=200)
tt2_up<-tt2[tt2$logFC>0,]
tt2_dn<-tt2[tt2$logFC<0,]
dim(tt2_up)
dim(tt2_dn)
head(tt2_up)
head(tt2_dn)# does not help much with the adjusted p-values

# Save gene symbol

k <- str_replace_all(row.names(CCLE_Exp_minus1),"_at","")# probe ID to Entrez ID
gene_name <- select(hgu133plus2.db, keys=k, columns=c("SYMBOL"), keytype="ENTREZID")# Entrez ID to symbol
head(gene_name)

sen_sig<-gene_name[fit$coefficients[,2]>0,]
dim(sen_sig)
sen_sig<-sen_sig[order(fit$p.value[fit$coefficients[,2]>0, 2])[1:100],]
head(sen_sig)
write.table(sen_sig$SYMBOL, file="glioma_sen_limma.txt", quote=F, row.names=F, col.names=F)

res_sig<-gene_name[fit$coefficients[,2]<0,]
dim(res_sig)
res_sig<-res_sig[order(fit$p.value[fit$coefficients[,2]<0, 2])[1:100],]
head(res_sig)
write.table(res_sig$SYMBOL, file="glioma_res_limma.txt", quote=F, row.names=F, col.names=F)

# naive Bayes classification

library(e1071)
library(stringr)
library(hgu133plus2.db)
library(caret)

?naiveBayes
dim(X)
Y<-as.factor(c(rep(1, 5),rep(2, 5)))#1:res, 2:sen
sig<-order(fit$p.value[,2])[1:50]
#sig<-c(row.names(sen_sig)[1:25],row.names(res_sig)[1:25])
class(sig)
sig<-as.numeric(sig)
X3<-X[sig,]# just use top 50 genes from res and sen signatures
dim(X3)# 100 genes
X3[1:10,1:5]
X3<-t(X3)

# Cross Validation

dim(X3)
length(Y)
CV1 <- matrix(data=NA, nrow=length(Y), ncol=5)
for(j in 1:5){
  idx <- createFolds(Y, k=5)
  for(i in 1:length(idx)){
    model <- naiveBayes(X3[-idx[[i]],], Y[-idx[[i]]])
    pred <- predict(model, X3[idx[[i]],], type = "raw")
    CV1[idx[[i]],j]<-pred[,2]
  }}
CV1

# ROC analysis

library(pROC)

Y_pred<-rowMeans(CV1)
roc1<-roc(Y, Y_pred, plot=TRUE, ci=TRUE)
auc(roc1)
# Area under the curve: 1
ci(roc1)
# 95% CI: 1-1 (DeLong)

# Predict all cell lines (including sen and res)

X4 <- rearranged[sig,]
dim(X4)
X4<-t(X4)
dim(X4)
X3[1:5,1:5]
X4[1:5,1:5]

identical(colnames(X4),colnames(X3))
model2 <- naiveBayes(X3, Y)
pred2 <- predict(model2, X4, type = "raw")
sen_pred<-pred2[,2]
res_pred<-pred2[,1]
mypar(1,2)
plot(res_pred, col=(as.numeric(res_pred==1)+1))
plot(sen_pred, col=(as.numeric(sen_pred==1)+1))
length(glioma_HSP90[,3])
mypar(1,1)
plot(glioma_HSP90[,3], col=(as.numeric(sen_pred==1)+1))
plot(glioma_HSP90[,3], col=(as.numeric(res_pred==1)+1))
plot(glioma_HSP90[,3], col=(as.numeric(sen_pred!=1&res_pred!=1)+1))
plot(density(glioma_HSP90[,3]))
lines(density(glioma_HSP90[sen_pred==1,3]), col=2)
lines(density(glioma_HSP90[res_pred==1,3]), col=3)
lines(density(glioma_HSP90[sen_pred!=1&res_pred!=1,3]), col=4)
t.test(glioma_HSP90[sen_pred==1,3],glioma_HSP90[res_pred==1,3])$p.value
# [1] 8.646668e-06
t.test(glioma_HSP90[sen_pred==1,3],glioma_HSP90[sen_pred!=1&res_pred!=1,3])$p.value
# [1] 0.3928271
t.test(glioma_HSP90[res_pred==1,3],glioma_HSP90[sen_pred!=1&res_pred!=1,3])$p.value
# [1] 0.1611635

mypar(1,2)
plot(sort(glioma_HSP90[,3]),col=(as.numeric(sen_pred==1)[order(glioma_HSP90[,3])]))
plot(sort(glioma_HSP90[,3]),col=(as.numeric(res_pred==1)[order(glioma_HSP90[,3])]))

# Predict all cell lines (excluding sen and res)

(res_idx2<-match(res$Primary.Cell.Line.Name, colnames(rearranged)))
(match(res$Primary.Cell.Line.Name, glioma_HSP90$Primary.Cell.Line.Name))# confirm match
(sen_idx2<-match(sen$Primary.Cell.Line.Name, colnames(rearranged)))
(match(sen$Primary.Cell.Line.Name, glioma_HSP90$Primary.Cell.Line.Name))
X5 <- rearranged[sig, -c(res_idx2,sen_idx2)]
HSP90_2<-glioma_HSP90[-c(res_idx2,sen_idx2),]
mypar(1,1)
plot(sort(glioma_HSP90[,3]))
points(sort(HSP90_2[,3]),col=1)
X5[1:5,1:5]
dim(X5)
dim(HSP90_2)
X5<-t(X5)
dim(X5)
identical(colnames(X5),colnames(X3))
pred3 <- predict(model2, X5, type = "raw")
sen_pred2<-pred3[,2]
res_pred2<-pred3[,1]
length(sen_pred2)
mypar(1,2)
plot(res_pred2, col=(as.numeric(res_pred2==1)+1))
plot(sen_pred2, col=(as.numeric(sen_pred2==1)+1))
length(HSP90_2[,3])
mypar(1,1)
plot(HSP90_2[,3], col=(as.numeric(sen_pred2==1)+1))
plot(HSP90_2[,3], col=(as.numeric(res_pred2==1)+1))
plot(HSP90_2[,3], col=(as.numeric(sen_pred2!=1&res_pred2!=1)+1))
plot(density(HSP90_2[,3]))
lines(density(HSP90_2[sen_pred2==1,3]), col=2)
lines(density(HSP90_2[res_pred2==1,3]), col=3)
lines(density(HSP90_2[sen_pred2!=1&res_pred2!=1,3]), col=4)
t.test(HSP90_2[sen_pred2==1,3],HSP90_2[res_pred2==1,3])$p.value
# [1] 0.003839564
t.test(HSP90_2[sen_pred2==1,3],HSP90_2[sen_pred2!=1&res_pred2!=1,3])$p.value
# [1] 0.5785514
t.test(HSP90_2[res_pred2==1,3],HSP90_2[sen_pred2!=1&res_pred2!=1,3])$p.value
# [1] 0.1877027

mypar(1,2)
plot(sort(HSP90_2[,3]),col=(as.numeric(sen_pred2==1)[order(HSP90_2[,3])]))
plot(sort(HSP90_2[,3]),col=(as.numeric(res_pred2==1)[order(HSP90_2[,3])]))

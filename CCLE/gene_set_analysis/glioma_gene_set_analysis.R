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
res<-glioma_HSP90[order(glioma_HSP90[,3])[1:5],]
sen<-glioma_HSP90[order(glioma_HSP90[,3])[25:29],]
res_idx<-match(res$Primary.Cell.Line.Name,colnames(glioma_Exp))
sen_idx<-match(sen$Primary.Cell.Line.Name,colnames(glioma_Exp))
res_exp<-glioma_Exp[,res_idx]
sen_exp<-glioma_Exp[,sen_idx]
colnames(res_exp)
colnames(sen_exp)
mypar(1,2)
plot(sort(glioma_HSP90[,3]))
plot(c(res[,3],sen[,3]))

# Prefilter the gene list(from http://www3.nd.edu/~steve/Rcourse/Lecture11v1.pdf)

library(genefilter)

X<-cbind(res_exp,sen_exp)
X[1:5,1:5]
dim(X)

f1<-pOverA(0.25, 3.5)#test genes for expression value over 3.5 in 25% of samples
f1fun<-filterfun(f1)
f1genes<-genefilter(X, f1fun)
f1genes[1:5]
sum(f1genes)#18889 genes
X1<-X[f1genes,]
dim(X1)

f2<-function(x){IQR(x)>0.75}#retain only genes with IQR>0.75
f2fun<-filterfun(f2)
f2genes<-genefilter(X1, f2fun)
sum(f2genes)
X2<-X1[f2genes,]
dim(X2)#4422 genes

boxplot(X2)
mypar(2,1)
summary(rowMeans(X2))
hist(rowMeans(X2), breaks=50,xlim = c(2,15))
summary(rowSds(X2))
hist(rowSds(X2), breaks=50,xlim = c(0,3))

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

# Save gene symbol (from non-filtered result)

k <- str_replace_all(row.names(CCLE_Exp_minus1),"_at","")# probe ID to Entrez ID
gene_name <- select(hgu133plus2.db, keys=k, columns=c("SYMBOL"), keytype="ENTREZID")# Entrez ID to symbol
head(gene_name)

sen_sig<-gene_name[fit$coefficients[,2]>0,]
dim(sen_sig)
sen_sig<-sen_sig[order(fit$p.value[fit$coefficients[,2]>0, 2])[1:100],]
head(sen_sig)
write.table(sen_sig$SYMBOL, file="glioma_sen_sig.txt", quote=F, row.names=F, col.names=F)

res_sig<-gene_name[fit$coefficients[,2]<0,]
dim(res_sig)
res_sig<-res_sig[order(fit$p.value[fit$coefficients[,2]<0, 2])[1:100],]
head(res_sig)
write.table(res_sig$SYMBOL, file="glioma_res_sig.txt", quote=F, row.names=F, col.names=F)

# Test for differentially expressed genes using MTP (permutated null distribution)
# http://www3.nd.edu/~steve/Rcourse/Lecture11v1.pdf

library(multtest)#


dim(X)
X[1:5,1:5]

length(fit$p.value[,2])
X3<-X[order(fit$p.value[,2])[1:200],]
dim(X3)

mtp<-MTP(X=X3, Y=Y, typeone = "fdr", B=500)#takes a long time to run (~30 min)
plot(mtp@adjp)
hist(mtp@rawp, breaks = 50)
sum(mtp@reject)
names(mtp@adjp)[mtp@reject]

k3 <- str_replace_all(names(mtp@rawp)[order(mtp@rawp)[1:20]],"_at","")
multi_sig <- select(hgu133plus2.db, keys=k3, columns=c("SYMBOL"), keytype="ENTREZID")

write.table(multi_sig$SYMBOL, file="glioma_multitest.txt", quote=F, row.names=F, col.names=F)

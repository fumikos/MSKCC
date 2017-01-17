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
install.packages("")
source("http://www.bioconductor.org/biocLite.R")
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
dim(CCLE_CELL_all)

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
HSP90 = CCLE_Drug_all[HSP90_index,c(1,13)]
head(HSP90)# cell name and ActArea

#Next, match cell name in HSP90 and CCLE_Exp_minus1
idx = match(HSP90$CCLE.Cell.Line.Name,colnames(CCLE_Exp_minus1))
length(which(is.na(idx)))#of the 503 cells 13 cell names do not exist in CCLE_Exp_minus1
HSP90$CCLE.Cell.Line.Name[which(is.na(idx))]
sum(colnames(CCLE_Exp_minus1)=="BGC823_STOMACH")#confirm above
idx <- idx[!is.na(idx)]
length(idx)#of the 503, 490 have microarray data
rearranged = CCLE_Exp_minus1[,idx]#extract expression data of cells that have HSP90i data and reorder to match HSP90 data order
rearranged[1:5,1:5]
HSP90[1:5,]

#Remove cell lines with no microarray data from HSP90 data
idx2 = match(colnames(rearranged), HSP90$CCLE.Cell.Line.Name)
length(which(is.na(idx2)))
HSP90<-HSP90[idx2,]
dim(HSP90)

#Assign resistant and sensitive cells
summary(HSP90)
res<-HSP90[order(HSP90[,2])[1:50],]
sen<-HSP90[order(HSP90[,2])[441:490],]
res_idx<-res_idx<-match(res$CCLE.Cell.Line.Name,colnames(CCLE_Exp_minus1))
sen_idx<-match(sen$CCLE.Cell.Line.Name,colnames(CCLE_Exp_minus1))
mypar(1,2)
plot(sort(HSP90[,2]))
plot(c(res[,2],sen[,2]))
match(res$CCLE.Cell.Line.Name,colnames(CCLE_Exp_minus1))
match(res$CCLE.Cell.Line.Name,CCLE_CELL_all$CCLE.name)
mypar(2,1)
table(CCLE_CELL_all$Site.Primary[res_idx])
table(CCLE_CELL_all$Site.Primary[sen_idx])
plot(CCLE_CELL_all$Site.Primary[res_idx])
plot(CCLE_CELL_all$Site.Primary[sen_idx])

     
#Simplify cell name
match(colnames(CCLE_Exp_minus1), CCLE_CELL_all$CCLE.name)
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
sum(rowMeans(X)<4.1 & rowSds(X)<0.3)#1693 genes--not much impact
sum(rowMeans(X)<4.1)#3025 genes
sum(rowSds(X)<0.3)#3248 genes

# PCA

s0<-svd(X-rowMeans(X))
PC1 <- s0$d[1]*s0$v[,1]
PC2 <- s0$d[2]*s0$v[,2]
mypar(1,1)
Y<-as.factor(c(rep(0, 50),rep(1, 50)))
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
sum(f1genes)#18949 genes
X1<-X[f1genes,]
dim(X1)

f2<-function(x){IQR(x)>0.75}#retain only genes with IQR>0.75
f2fun<-filterfun(f2)
f2genes<-genefilter(X1, f2fun)
sum(f2genes)
X2<-X1[f2genes,]
dim(X2)#7922 genes

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
tail(tt_up)
tail(tt_dn)

#plot(X["NQO1",],col=(X["NQO1",]<10)+1)

# Fit all genes adjusting with surrogate variables
#s <- sva(X,mod)# Takes a while to run
#svamod<-model.matrix(~Y+s$sv)# 18 significant SVs -- this data probably does not require correction
#svafit <- lmFit(X, design=svamod)
#svafit <- eBayes(svafit)
#svatt <- topTable(svafit, coef=2, number=50)
#svatt$ID

# Fit filtered genes without adjusting with surrogate variables
dim(X2)
fit2 <- lmFit(X2, design=mod)
fit2 <- eBayes(fit2)
tt2 <- topTable(fit2, coef=2, number=200)
tt2_up<-tt2[tt2$logFC>0,]
tt2_dn<-tt2[tt2$logFC<0,]
dim(tt2_up)
dim(tt2_dn)
head(tt2_up)
head(tt2_dn)
tail(tt2_up)
tail(tt2_dn)

# Fit filtered genes adjusting with surrogate variables
#s2 <- sva(X2,mod)# Takes a while to run
#svamod2<-model.matrix(~Y+s2$sv)# 16 significant SVs
#svafit2 <- lmFit(X, design=svamod)
#svafit2 <- eBayes(svafit2)
#svatt2 <- topTable(svafit2, coef=2, number=50)
#svatt2$ID

# Save gene symbol

k <- str_replace_all(row.names(CCLE_Exp_minus1),"_at","")# probe ID to Entrez ID
gene_name <- select(hgu133plus2.db, keys=k, columns=c("SYMBOL"), keytype="ENTREZID")# Entrez ID to symbol
head(gene_name)
dim(gene_name)

sen_sig<-gene_name[fit$coefficients[,2]>0,]
dim(sen_sig)
sen_sig<-sen_sig[order(fit$p.value[fit$coefficients[,2]>0, 2])[1:100],]
head(sen_sig)
write.table(sen_sig$SYMBOL, file="sen_limma.txt", quote=F, row.names=F, col.names=F)

res_sig<-gene_name[fit$coefficients[,2]<0,]
dim(res_sig)
res_sig<-res_sig[order(fit$p.value[fit$coefficients[,2]<0, 2])[1:100],]
head(res_sig)
write.table(res_sig$SYMBOL, file="res_limma.txt", quote=F, row.names=F, col.names=F)

gene_name2<-gene_name[f1genes,]
gene_name2<-gene_name2[f2genes,]
dim(gene_name2)
dim(fit2$coefficients)
sen_sig2<-gene_name2[fit2$coefficients[,2]>0,]
dim(sen_sig2)
sen_sig2<-sen_sig2[order(fit2$p.value[fit2$coefficients[,2]>0, 2])[1:100],]
head(sen_sig2)
write.table(sen_sig2$SYMBOL, file="sen_limma_filtered.txt", quote=F, row.names=F, col.names=F)

res_sig2<-gene_name2[fit2$coefficients[,2]<0,]
dim(res_sig2)
res_sig2<-res_sig2[order(fit2$p.value[fit2$coefficients[,2]<0, 2])[1:100],]
head(res_sig2)
write.table(res_sig2$SYMBOL, file="res_limma_filtered.txt", quote=F, row.names=F, col.names=F)

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

write.table(multi_sig$SYMBOL, file="multtest.txt", quote=F, row.names=F, col.names=F)

# mroast
# http://genomicsclass.github.io/book/pages/bioc1_roast.html

library(org.Hs.eg.db)# Gene set database

org.Hs.egGO2ALLEGS
go2eg <- org.Hs.egGO2ALLEGS
head(go2eg)

govector <- unlist(go2eg)# unlist the GO list into vector
head(govector)
golengths <- sapply(go2eg, length)

X3<-cbind(res_exp,sen_exp)# reconstruct X to recover probe ID
row.names(X3) <- str_replace_all(row.names(X3),"_at","")# convert to Entrez ID

head(row.names(X3))# Entrez ID
idxvector <- match(govector, row.names(X3))# get matches for each Entrez ID to the index in X3
head(idxvector)
table(is.na(idxvector))
idx_list <- split(idxvector, rep(names(go2eg), golengths))# rebuild list
go2eg[[1]]
row.names(X3)[idx_list[[1]]]# confirm that genes in the first gene set matches

idxclean <- lapply(idx, function(x) x[!is.na(x)])# remove no values
idxlengths <- sapply(idxclean, length)
idxsub <- idxclean[idxlengths > 10]# remove gene sets which has less than 10 genes
length(idxsub)# 7354 gene sets are left

?mroast
r <- mroast(X3, idxsub, design)
head(r)
r <- r[order(r$PValue),]

r_up<-r[r$Direction == "Up",]
r_up$PValue[1:15]
sum(r_up$PValue==0.001)# 240 gene sets

r_dn<-r[r$Direction == "Down",]
r_dn$PValue[1:50]
sum(r_dn$PValue<=0.005)# 56 gene sets

library(GO.db)#

rtab_up <- select(GO.db, keys=rownames(r_up)[1:15],
                columns=c("GOID","TERM","DEFINITION"), keytype="GOID")
rtab_up[,1:2]

rtab_dn <- select(GO.db, keys=rownames(r_dn)[1:15],
                columns=c("GOID","TERM","DEFINITION"), keytype="GOID")
rtab_dn[,2]

#make a res file for GSEA (no need?)
colnames(CCLE_Exp_minus1[,res_idx])
colnames(CCLE_Exp_minus1[,sen_idx])

keytypes(hgu133plus2.db)
k<-as.character(row.names(X))
refseqID <- select(hgu133plus2.db, keys=k, columns=c("REFSEQ"), keytype="ENTREZID")
dim(refseqID)

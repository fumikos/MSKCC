library(Biobase)
library(CePa)
library(dplyr)
library(stats)
library(glmnet)
library(parallel)
library(doParallel)
parallel:::detectCores()
library(stringr)
library(hgu133plus2.db)
library(GenomicFeatures)
library(AnnotationDbi)
library(ggplot2)
library(caret)
library(rafalib)
library(genefilter)
library(BiocInstaller)
library(genefilter)
library(gplots)
library(RColorBrewer)
library(devtools)
library(Heatplus) 
library(vegan)
library(marray)
library(Rcpp)
library(lme4)
library(nnet)
library(quantreg)
library(e1071)
library(pROC)
library(OptimalCutpoints)
library(limma)
library(GEOquery)
library(org.Hs.eg.db)
library(GO.db)
library(rJava)
install.packages("rJava")
source("http://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")

#CCLE_Drug_all = read.csv("G:/My Documents/CCLE/CCLE_analysis/CCLE_NP24.2009_Drug_data_2015.02.24.csv")
#CCLE_Exp_all =read.gct("G:/My Documents/CCLE/CCLE_analysis/CCLE_Expression_Entrez_2012-09-29.gct")
#CCLE_CELL_all <- read.csv("G:/My Documents/CCLE/CCLE_analysis/Cell_tissue.csv")

idx0 <- match(CCLE_CELL_all$CCLE.name,colnames(CCLE_Exp_all))
length(which(is.na(idx0)))#17 colnames of CCLE_Exp_all have letter "X" before cell name
correct_name<-as.vector(CCLE_CELL_all$CCLE.name[which(is.na(idx0))])
length(correct_name)
CCLE_Exp_minus1 <- CCLE_Exp_all[,-991] # remove "NCIH292_LUNG.1" from CCLE_Exp_all
colnames(CCLE_Exp_minus1)[which(is.na(idx0))]<-correct_name

HSP90_index = which(CCLE_Drug_all$Target=="HSP90")#503 cells have 17AAG data
HSP90 = CCLE_Drug_all[HSP90_index,c(1,13)]

idx = match(HSP90$CCLE.Cell.Line.Name,colnames(CCLE_Exp_minus1))
length(which(is.na(idx)))#of the 503 cells 13 do not have microarray data
idx <- idx[!is.na(idx)]
length(idx)#of the 503, 490 have microarray data
rearranged = CCLE_Exp_minus1[,idx]
rearranged[1:5,1:5]
HSP90[1:5,]

idx2 = match(colnames(rearranged), HSP90$CCLE.Cell.Line.Name)
length(which(is.na(idx2)))
HSP90<-HSP90[idx2,]

summary(HSP90)
res<-HSP90[order(HSP90[,2])[1:100],]
sen<-HSP90[order(HSP90[,2])[391:490],]
res_idx<-match(res$CCLE.Cell.Line.Name,colnames(CCLE_Exp_minus1))
sen_idx<-match(sen$CCLE.Cell.Line.Name,colnames(CCLE_Exp_minus1))
mypar(1,2)
plot(sort(HSP90[,2]))
plot(c(res[,2],sen[,2]))

match(colnames(CCLE_Exp_minus1), CCLE_CELL_all$CCLE.name)
cell_name =CCLE_CELL_all$Cell.line.primary.name
colnames(CCLE_Exp_minus1)[1:5]
res_exp<-CCLE_Exp_minus1[,res_idx]
sen_exp<-CCLE_Exp_minus1[,sen_idx]
res_exp[1:5,1:5]
sen_exp[1:5,1:5]



colnames(res_exp)<-cell_name[res_idx]
colnames(sen_exp)<-cell_name[sen_idx]
res_exp[1:5,1:5]
sen_exp[1:5,1:5]


X<-cbind(res_exp,sen_exp)
X[1:5,1:5]

row.names(X) <- str_replace_all(row.names(X),"_at","")
k<-as.character(row.names(X))
gene_name <- select(hgu133plus2.db, keys=k, columns=c("SYMBOL"), keytype="ENTREZID")
row.names(X)[1:10]
gene_name[1:10,]
row.names(X)<-gene_name[,2]
X[1:5,1:5]
X[1:5,96:100]

write.csv(X, "expression_100each_SYNBL.csv")

# lmFit
dim(X)
Y<-as.factor(c(rep(0, 100),rep(1, 100)))#0:res, 1:sen
row.names(X) <- str_replace_all(row.names(X),"_at","")
X[1:5,1:5]
design <- model.matrix(~ Y)
fit <- lmFit(X, design=design)
fit <- eBayes(fit)
?topTable
tt <- topTable(fit, coef=2)
tt

# mroast
org.Hs.egGO2ALLEGS
go2eg <- as.list(org.Hs.egGO2ALLEGS)
head(go2eg)

govector <- unlist(go2eg)
head(govector)
golengths <- sapply(go2eg, length)
head(row.names(X))
idxvector <- match(govector, row.names(X))
head(idxvector)
table(is.na(idxvector))
idx <- split(idxvector, rep(names(go2eg), golengths))
go2eg[[1]]
row.names(X)[idx[[1]]]

idxclean <- lapply(idx, function(x) x[!is.na(x)])
idxlengths <- sapply(idxclean, length)
idxsub <- idxclean[idxlengths > 10]
length(idxsub)

?mroast
r <- mroast(X, idxsub, design)
head(r)
r <- r[order(r$PValue),]

r_up<-r[r$Direction == "Up",]
r_up$PValue[1:50]
sum(r_up$PValue==0.001)

r_dn<-r[r$Direction == "Down",]
r_dn$PValue[1:50]
sum(r_dn$PValue<=0.005)

rtab_up <- select(GO.db, keys=rownames(r_up)[1:232],
                columns=c("GOID","TERM","DEFINITION"), keytype="GOID")
rtab_up[,1:2]

rtab_dn <- select(GO.db, keys=rownames(r_dn)[1:60],
                columns=c("GOID","TERM","DEFINITION"), keytype="GOID")
rtab_dn[,2]

#make a res file for GSEA (no need?)
colnames(CCLE_Exp_minus1[,res_idx])
colnames(CCLE_Exp_minus1[,sen_idx])

keytypes(hgu133plus2.db)
k<-as.character(row.names(X))
refseqID <- select(hgu133plus2.db, keys=k, columns=c("REFSEQ"), keytype="ENTREZID")
dim(refseqID)

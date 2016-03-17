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
install.packages("")
biocLite("")

CCLE_Exp_all =read.gct("G:/My Documents/CCLE/CCLE_analysis/CCLE_Expression_Entrez_2012-09-29.gct")
CCLE_CELL_all <- read.csv("G:/My Documents/CCLE/CCLE_analysis/Cell_tissue.csv")

ls("package:hgu133plus2.db")

class(hgu133plus2.db)
keytypes(hgu133plus2.db)
columns(hgu133plus2.db)
length(keys(hgu133plus2.db))

idx0 <- match(CCLE_CELL_all$CCLE.name,colnames(CCLE_Exp_all))
length(which(is.na(idx0)))#17 colnames of CCLE_Exp_all have letter "X" before cell name
(correct_name<-CCLE_CELL_all$CCLE.name[which(is.na(idx0))])
correct_name<-as.vector(correct_name)
correct_name
CCLE_Exp_minus1 <- CCLE_Exp_all[,-991] # remove "NCIH292_LUNG.1" from CCLE_Exp_all
colnames(CCLE_Exp_minus1)[which(is.na(idx0))]
colnames(CCLE_Exp_minus1)[which(is.na(idx0))]<-correct_name
colnames(CCLE_Exp_minus1)[which(is.na(idx0))]

CCLE_Exp_minus1[1:10,1:4]
head(CCLE_CELL_all)

idx1 <- match(CCLE_CELL_all$CCLE.name,colnames(CCLE_Exp_minus1))
colnames(CCLE_Exp_minus1)<-CCLE_CELL_all$Cell.line.primary.name[idx1]

cell_idx=c("NCI-H1975", "MDA-MB-468", "Daudi", "MOLM-13","AsPC-1","MCF7","MDA-MB-415")
cell_exp=CCLE_Exp_minus1[,cell_idx]
head(cell_exp)
row.names(cell_exp) <- str_replace_all(row.names(cell_exp),"_at","")
k<-as.character(row.names(cell_exp))
gene_name <- select(hgu133plus2.db, keys=k, columns=c("SYMBOL"), keytype="ENTREZID")
cell_exp[1:10,]

dim(gene_name)
head(gene_name)
tail(gene_name)
head(k)
tail(k)
idx2 <- match(uniprot_name$ ENTREZID,k)
cell_exp<-cbind(gene_name,cell_exp)


write.csv(cell_exp, file="microarray_log2.csv")

file<-read.csv("microarray_log2.csv")
file[1:10,]

map <- select(hgu133plus2.db, keys=k, columns=c("SYMBOL","UNIPROT"), keytype="ENTREZID")
head(map)
dim(map)
write.csv(map, file="map_to_uniprot.csv")

map_file<-read.csv("map_to_uniprot.csv")
head(map_file)

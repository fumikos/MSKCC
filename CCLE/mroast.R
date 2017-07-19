### Roatation gene set analysis (MROAST)

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
Y <- as.factor(c(rep("res", 50),rep("sen", 50)))

## mroast
## http://genomicsclass.github.io/book/pages/bioc1_roast.html
library(org.Hs.eg.db)# Gene set database
library(AnnotationDbi)
org.Hs.egGO2ALLEGS# a later version than org.Hs.egGO2EG
go2eg <- AnnotationDbi::as.list(org.Hs.egGO2ALLEGS)
head(go2eg)
govector <- unlist(go2eg)# unlist the GO list into vector
head(govector)
golengths <- sapply(go2eg, length)
head(row.names(X))# Entrez ID
idxvector <- match(govector, row.names(X))# get matches for each Entrez ID to the index in X
head(idxvector)
table(is.na(idxvector))
idx_list <- split(idxvector, rep(names(go2eg), golengths))# rebuild list
go2eg[[1]]
row.names(X)[idx_list[[1]]]# confirm that genes in the first gene set matches
idxclean <- lapply(idx_list, function(x) x[!is.na(x)])# remove no values
idxlengths <- sapply(idxclean, length)
idxsub <- idxclean[idxlengths > 10]# remove gene sets which has less than 10 genes
length(idxsub)# 7354 gene sets are left
library(limma)
mod <- model.matrix(~Y)
r <- mroast(X, idxsub, design=mod)
head(r)
r <- r[order(r$PValue),]
library(GO.db)
columns(GO.db)
keytypes(GO.db)
GOTERM[[rownames(r)[1]]]
rtab_up <- select(GO.db, keys=rownames(r)[r$Direction == "Up"][1:15],
                columns=c("GOID","TERM","DEFINITION"), keytype="GOID")
rtab_up[,1:2]
rtab_dn <- select(GO.db, keys=rownames(r)[r$Direction == "Down"][1:15],
                columns=c("GOID","TERM","DEFINITION"), keytype="GOID")
rtab_dn[,1:2]



### Gene set testing using the limma package
### Method is based on course materials from EdX PH525x series - Biomedical Data Science 
### Data downloaded from CCLE website (https://portals.broadinstitute.org/ccle/home)

## Run code in "CCLE_drug_sensitivity_prediction\CCLE_data.R" to load and preprocess the data
dim(rearranged)# expression data for the drug of interest
rearranged[1:5,1:5]
dim(drug_data)# drug response data (activity area)
head(drug_data)

## Run code in "CCLE_drug_sensitivity_prediction\CCLE_EDA.R" to do EDA (optional)

## Assign resistant and sensitive cells (top 50 for each)
dim(drug_data)
res <- drug_data[order(drug_data[,2])[1:50],]
sen <- drug_data[order(drug_data[,2])[(nrow(drug_data)-49):nrow(drug_data)],]
res_idx <- match(res$Primary.Cell.Line.Name,colnames(CCLE_Exp_minus1))
sen_idx <- match(sen$Primary.Cell.Line.Name,colnames(CCLE_Exp_minus1))
library(rafalib)
mypar(1,2)
plot(sort(drug_data[,2]))
plot(c(res[,2],sen[,2]))
X <- CCLE_Exp_minus1[,c(res_idx,sen_idx)]
dim(X)
X[1:5,1:5]
Y <- as.factor(c(rep("res", 50),rep("sen", 50)))

## Load gene set database and remove genes that are not present in X
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
head(idx_list)
go2eg[[1]]
row.names(X)[idx_list[[1]]]# confirm that genes in the first gene set matches
idxclean <- lapply(idx_list, function(x) x[!is.na(x)])# remove no values
head(idxclean)
idxlengths <- sapply(idxclean, length)
idxsub <- idxclean[idxlengths > 10]# remove gene sets which has less than 10 genes
length(idxsub)# Number of gene sets left

## Run a self-contained rotation-based gene set testing method (MROAST)
library(limma)
mod <- model.matrix(~Y)
r <- mroast(X, idxsub, design=mod)
head(r)
sum(r$FDR[!is.na(r$FDR)]<0.25)
r <- r[which(r$FDR<0.25),]
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

## Run a competitive gene set testing method (CAMERA)
c <- camera(X, idxsub, design=mod)
head(c)
sum(c$FDR[!is.na(c$FDR)]<0.25)
c <- c[which(c$FDR<0.25),]
c <- c[order(c$PValue),]
ctab_up <- select(GO.db, keys=rownames(c)[c$Direction == "Up"][1:15],
                  columns=c("GOID","TERM","DEFINITION"), keytype="GOID")
ctab_up[,1:2]
ctab_dn <- select(GO.db, keys=rownames(c)[c$Direction == "Down"][1:15],
                  columns=c("GOID","TERM","DEFINITION"), keytype="GOID")
ctab_dn[,1:2]

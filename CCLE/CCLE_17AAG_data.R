### Prepare gene expression and drug response data for downstream analysis

## Load gene expression and drug response data from CCLE data
library(CePa)
#CCLE_Drug_all <- read.csv("G:/My Documents/CCLE/CCLE_data/CCLE_NP24.2009_Drug_data_2015.02.24.csv")
#CCLE_Exp_all <- read.gct("G:/My Documents/CCLE/CCLE_data/CCLE_Expression_Entrez_2012-09-29.gct")
#CCLE_CELL_all <- read.csv("G:/My Documents/CCLE/CCLE_data/Cell_tissue.csv")
head(CCLE_Drug_all)
CCLE_Exp_all[1:5,1:5]
head(CCLE_CELL_all)

## Check if the cell ID matches in CCLE_Exp_all and CCLE_CELL_all
idx0 <- match(CCLE_CELL_all$CCLE.name, colnames(CCLE_Exp_all))
length(which(is.na(idx0)))# 17 do not match
CCLE_CELL_all$CCLE.name[which(is.na(idx0))]
colnames(CCLE_Exp_all)[which(is.na(idx0))]# 17 cells have letter "X" before its name
length(CCLE_CELL_all$CCLE.name)# 1036
length(colnames(CCLE_Exp_all))# 1037: CCLE_Exp_all also has an extra column
which(is.na(match(CCLE_CELL_all$CCLE.name, colnames(CCLE_Exp_all))))# Which one is the extra?
which(is.na(match(colnames(CCLE_Exp_all), CCLE_CELL_all$CCLE.name)))# 991 is the extra one
CCLE_Exp_minus1 <- CCLE_Exp_all[,-991] # remove row 991 (NCIH292_LUNG.1) from CCLE_Exp_all
correct_name <- as.vector(CCLE_CELL_all$CCLE.name[which(is.na(idx0))])
colnames(CCLE_Exp_minus1)[which(is.na(idx0))] <- correct_name
identical(colnames(CCLE_Exp_minus1), as.character(CCLE_CELL_all$CCLE.name))# the two now match

## Change probeID to EntrezID and simplify the cell names
library(stringr)
rownames(CCLE_Exp_minus1) <- str_replace_all(row.names(CCLE_Exp_minus1),"_at","")
colnames(CCLE_Exp_minus1) <- CCLE_CELL_all$Cell.line.primary.name
CCLE_Exp_minus1[1:5,1:5]

## Next, extract Activity Area data for HSP90 inhibitor
head(CCLE_Drug_all)
HSP90_index <- which(CCLE_Drug_all$Target=="HSP90")
length(HSP90_index)# 503 cells have 17AAG data
HSP90 <- CCLE_Drug_all[HSP90_index,c(2,13)]
head(HSP90)

## Match cell name in HSP90 and CCLE_Exp_minus1
idx1 <- match(HSP90$Primary.Cell.Line.Name, colnames(CCLE_Exp_minus1))
length(which(is.na(idx1)))# of the 503 cells 11 cell names do not exist in CCLE_Exp_minus1
idx1 <- idx1[!is.na(idx1)]
length(idx1)# of the 503, 492 have microarray data
rearranged <- CCLE_Exp_minus1[,idx1]# extract/reorder expression data to match HSP90 data
dim(rearranged)# 18988   492
rearranged[1:5,1:5]

## Remove cell lines with no microarray data from HSP90 data
idx2 <- match(colnames(rearranged), HSP90$Primary.Cell.Line.Name)
length(idx2)
length(which(is.na(idx2)))
HSP90 <- HSP90[idx2,]
dim(HSP90)# 492   2
head(HSP90)

## Confirm cell name match between expression and drug response data
identical(colnames(rearranged), as.character(HSP90$Primary.Cell.Line.Name))# TRUE
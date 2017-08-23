### Preprocess gene expression and drug response data for downstream analysis
### Data downloaded from CCLE website (https://portals.broadinstitute.org/ccle/home)
### Enter the name of your drug of interst on line 46 before running this code

## Load gene expression and drug response data from CCLE data
library(CePa)
CCLE_Drug_all <- read.csv("CCLE_NP24.2009_Drug_data_2015.02.24.csv")
CCLE_Exp_all <- read.gct("CCLE_Expression_Entrez_2012-09-29.gct")
CCLE_CELL_all <- read.csv("Cell_tissue.csv")
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
colnames(CCLE_Exp_all)[991]# 991 is "NCIH292_LUNG.1"
which(grepl("NCIH292", CCLE_CELL_all$CCLE.name))# 963
which(grepl("NCIH292", colnames(CCLE_Exp_all)))# 963 991
head(CCLE_Exp_all[,c(963,991)])
cor(CCLE_Exp_all[,963],CCLE_Exp_all[,991])# 0.9804931
CCLE_Exp_minus1 <- CCLE_Exp_all[,-991] # remove row 991 (NCIH292_LUNG.1) from CCLE_Exp_all
correct_name <- as.vector(CCLE_CELL_all$CCLE.name[which(is.na(idx0))])
colnames(CCLE_Exp_minus1)[which(is.na(idx0))] <- correct_name
identical(colnames(CCLE_Exp_minus1), as.character(CCLE_CELL_all$CCLE.name))# the two now match

## Change probeID to EntrezID and simplify the cell names
library(stringr)
rownames(CCLE_Exp_minus1) <- str_replace_all(row.names(CCLE_Exp_minus1),"_at","")
colnames(CCLE_Exp_minus1) <- CCLE_CELL_all$Cell.line.primary.name
CCLE_Exp_minus1[1:5,1:5]

## Next, extract Activity Area data for the drug of interest
head(CCLE_Drug_all)
levels(CCLE_Drug_all$Compound)
# [1] "17-AAG"       "AEW541"       "AZD0530"      "AZD6244"      "Erlotinib"    "Irinotecan"  
# [7] "L-685458"     "Lapatinib"    "LBW242"       "Nilotinib"    "Nutlin-3"     "Paclitaxel"  
# [13] "Panobinostat" "PD-0325901"   "PD-0332991"   "PF2341066"    "PHA-665752"   "PLX4720"     
# [19] "RAF265"       "Sorafenib"    "TAE684"       "TKI258"       "Topotecan"    "ZD-6474" 
drug_index <- which(CCLE_Drug_all$Compound=="DRUG")
length(drug_index)
drug_data <- CCLE_Drug_all[drug_index,c(2,13)]# save cell line name and activity area
head(drug_data)

## Match cell name in HSP90 and CCLE_Exp_minus1
idx1 <- match(drug_data$Primary.Cell.Line.Name, colnames(CCLE_Exp_minus1))
length(which(is.na(idx1)))# these cell names do not exist in CCLE_Exp_minus1
idx1 <- idx1[!is.na(idx1)]
length(idx1)# number of cells which have both microarray and drug response data
rearranged <- CCLE_Exp_minus1[,idx1]# extract/reorder expression data to match drug_data
dim(rearranged)
rearranged[1:5,1:5]

## Remove cell lines with no microarray data from drug_data
idx2 <- match(colnames(rearranged), drug_data$Primary.Cell.Line.Name)
length(idx2)
length(which(is.na(idx2)))
drug_data <- drug_data[idx2,]
dim(drug_data)
head(drug_data)

## Confirm cell name match between expression and drug response data
identical(colnames(rearranged), as.character(drug_data$Primary.Cell.Line.Name))# TRUE

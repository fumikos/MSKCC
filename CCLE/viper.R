### Network-based inference of protein activity using ARACNE-AP/VIPER algorithm

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

## ARACNE-AP
## Expression profile for lung cancer cells
X <- CCLE_Exp_minus1[,which(CCLE_CELL_all$Site.Primary=="lung")]
dim(X)
X <- data.frame(X)
X[1:5,1:5]
gene <- rownames(X)
X <- cbind(gene, X)
X[1:5,1:5]
library(data.table)
# fwrite(X, "expression_lung.txt", quote = FALSE, sep = "\t", row.names = FALSE)

## Extract gene list of regulators
library(aracne.networks)
data(regulongbm)
class(regulongbm)
regulongbm# Object of class regulon with 6056 regulators, 19858 targets and 563850 interactions
regulators <- names(regulongbm)
length(regulators)# 6056
head(regulators)
# write.table(regulators,"regulators.txt", quote=F, row.names=F, col.names=F)

## On the command line, run "run_aracne". This generates an aracne network file ("network.txt")
## The first line (column label) had to be removed in order for the aracne2regulon function to work

## Create regulon objects
library(viper)
class(X)# "data.frame"
X[1:5,1:5]
X <- CCLE_Exp_minus1[,which(CCLE_CELL_all$Site.Primary=="lung")]
class(X)# "matrix"
X[1:5,1:5]
regs<- aracne2regulon("network.txt", X, format = "3col")
regs# Object of class regulon with 5408 regulators, 18298 targets and 110112 interactions

## Assign resistant and sensitive lung cancer cells
idx3 <- match(HSP90$Primary.Cell.Line.Name, CCLE_CELL_all$Cell.line.primary.name)
CCLE_CELL_HSP90 <- CCLE_CELL_all[idx3,]
head(CCLE_CELL_HSP90)
sum(CCLE_CELL_HSP90$Site.Primary=="lung")#91 lung cells have 17-AAG data
HSP90_lung <- HSP90[CCLE_CELL_HSP90$Site.Primary=="lung",]
dim(HSP90_lung)
head(HSP90_lung)
res_lung <- HSP90_lung[order(HSP90_lung[,2])[1:15],]
sen_lung <- HSP90_lung[order(HSP90_lung[,2])[77:91],]
mid_lung <- HSP90_lung[sample(order(HSP90_lung[,2])[30:62],15),]
res_lung_idx <- match(res_lung$Primary.Cell.Line.Name,colnames(X))
sen_lung_idx <- match(sen_lung$Primary.Cell.Line.Name,colnames(X))
mid_lung_idx <- match(mid_lung$Primary.Cell.Line.Name,colnames(X))
library(rafalib)
mypar(1,3)
plot(sort(HSP90[,2]))
plot(sort(HSP90_lung[,2]))
plot(sort(c(res_lung[,2],mid_lung[,2],sen_lung[,2])))

## Create an ExpressionSet object
eset <- ExpressionSet(X, phenoData = annotatedDataFrameFrom(X, byrow = FALSE))
description <- rep("N", ncol(X))
class(description)
description[sen_lung_idx] <- "sen"
description[res_lung_idx] <- "res"
description[mid_lung_idx] <- "mid"
pData(eset)$description <- description
table(pData(eset)$description)

description2 <- rep("N", ncol(X))
pData(eset)$description2 <- description2

## Generating the gene expression signatures
sig1 <- rowTtest(eset, "description", "sen", "res")
mypar(1,1)
hist(sig1$p.value)# uniform
sig1 <- (qnorm(sig1$p.value/2, lower.tail = FALSE)*sign(sig1$statistic))[, 1]
hist(sig1)

## Generating the NULL model by sample permutations
null1 <- ttestNull(eset, "description", "sen", "res", per = 1000, repos = TRUE, verbose = FALSE)
dim(null1)
null1[1:5,1:5]

## msVIPER
mrs <- msviper(sig1, regs, null1)
summary(mrs)
plot(mrs, cex=.7)
topten <- as.character(summary(mrs)$Regulon)
library(hgu133plus2.db)
gene_name <- select(hgu133plus2.db, keys=topten, columns=c("SYMBOL"), keytype="ENTREZID")
gene_name$SYMBOL

## Bootstrap msVIPER
sig_bs <- bootstrapTtest(eset, "description", "sen", "res")
dim(sig_bs)
sig_bs[1:5,1:5]
mrs_bs <- msviper(sig_bs, regs, null1)
mrs_bs <- bootstrapmsviper(mrs_bs, "mode")
summary(mrs_bs)
plot(mrs_bs, cex = .7)
topten2 <- as.character(summary(mrs_bs)$Regulon)
gene_name2 <- select(hgu133plus2.db, keys=topten2, columns=c("SYMBOL"), keytype="ENTREZID")
gene_name2$SYMBOL

## VIPER
vpres <- viper(eset, regs)
dim(vpres)#1535 186
exprs(vpres)[1:5,1:5]
head(pData(vpres))
sig2 <- rowTtest(vpres, "description", "sen", "res")
length(sig2)
hist(sig2$p.value)# uniform
df <- data.frame(Gene = rownames(sig2$p.value), t = round(sig2$statistic, 2),
               "p-value" = signif(sig2$p.value, 3))[order(sig2$p.value)[1:30], ]
df
topten3 <- as.character(df$Gene)
gene_name3 <- select(hgu133plus2.db, keys=topten3, columns=c("SYMBOL"), keytype="ENTREZID")
gene_name3$SYMBOL

## Running VIPER with a null model
## An additional column had to be added to pData(eset) in order for the viperSignature function to work
description2 <- rep("N", ncol(X))
pData(eset)$description2 <- description2
vpsig <- viperSignature(eset, "description", c("mid","N"))
vpres <- viper(vpsig, regs)

pos <- pData(vpres)[["description"]] %in% c("sen", "res")
d1 <- exprs(vpres)[, pos]
colnames(d1) <- pData(vpres)[["description"]][pos]
dd <- dist(t(d1), method = "euclidean")
heatmap(as.matrix(dd), Rowv = as.dendrogram(hclust(dd, method = "average")), symm = T)
dd <- viperSimilarity(d1)
heatmap(as.matrix(as.dist(dd)), Rowv = as.dendrogram(hclust(as.dist(dd), method = "average")), symm = T)

### Network-based inference of protein activity using ARACNE-AP/VIPER algorithm
### Method is based on VIPER's vignette
### https://www.bioconductor.org/packages/devel/bioc/vignettes/viper/inst/doc/viper.pdf
### Data downloaded from CCLE website (https://portals.broadinstitute.org/ccle/home)

## Run code in "CCLE_drug_sensitivity_prediction\CCLE_data.R" to load and preprocess the data
dim(rearranged)# expression data for the drug of interest
rearranged[1:5,1:5]
dim(drug_data)# drug response data (activity area)
head(drug_data)
sort(table(CCLE_CELL_all$Site.Primary), decreasing = T)[1]# Tissue with largest number of data
# lung 
# 186

## Extract expression data for lung cancer cells to build an context-specific interactome
X <- CCLE_Exp_minus1[,which(CCLE_CELL_all$Site.Primary=="lung")]
dim(X)
X <- data.frame(X)
X[1:5,1:5]
gene <- rownames(X)
X <- cbind(gene, X)
X[1:5,1:5]
library(data.table)
fwrite(X, "expression_lung.txt", quote = FALSE, sep = "\t", row.names = FALSE)

## Extract a gene list of regulators (transcription factors and signaling proteins)
## from a pre-assembled ARACNe interactome
library(aracne.networks)
data(regulongbm)
regulongbm# Object of class regulon with 6056 regulators, 19858 targets and 563850 interactions
regulators <- names(regulongbm)
length(regulators)# 6056
head(regulators)
write.table(regulators,"regulators.txt", quote=F, row.names=F, col.names=F)

## On the command line, run "run_aracne.sh". This generates an aracne network file ("network.txt")
## The first line (column label) may need to be removed in order for the aracne2regulon function 
## below to work.

## Create regulon objects
library(viper)
X <- CCLE_Exp_minus1[,which(CCLE_CELL_all$Site.Primary=="lung")]# Convert X back to a matrix
class(X)# "matrix"
regs<- aracne2regulon("network.txt", X, format = "3col")
regs# Object of class regulon with 5408 regulators, 18298 targets and 110112 interactions

## Assign resistant and sensitive lung cancer cells
idx3 <- match(drug_data$Primary.Cell.Line.Name, CCLE_CELL_all$Cell.line.primary.name)
CCLE_CELL_drug_data <- CCLE_CELL_all[idx3,]
head(CCLE_CELL_drug_data)
sum(CCLE_CELL_drug_data$Site.Primary=="lung")
drug_data_lung <- drug_data[CCLE_CELL_drug_data$Site.Primary=="lung",]
dim(drug_data_lung)
head(drug_data_lung)
res <- 1:15
sen <- (length(drug_data_lung$ActArea)-14):length(drug_data_lung$ActArea)
mid <- (round(length(drug_data_lung$ActArea)/2)-15):(round(length(drug_data_lung$ActArea)/2)+15)
res_lung <- drug_data_lung[order(drug_data_lung$ActArea)[res],]
sen_lung <- drug_data_lung[order(drug_data_lung$ActArea)[sen],]
mid_lung <- drug_data_lung[sample(order(drug_data_lung$ActArea)[mid],15),]
res_lung_idx <- match(res_lung$Primary.Cell.Line.Name,colnames(X))
sen_lung_idx <- match(sen_lung$Primary.Cell.Line.Name,colnames(X))
mid_lung_idx <- match(mid_lung$Primary.Cell.Line.Name,colnames(X))
library(rafalib)
mypar(1,3)
plot(sort(drug_data[,2]))
plot(sort(drug_data_lung[,2]))
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
hist(sig1$p.value)
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

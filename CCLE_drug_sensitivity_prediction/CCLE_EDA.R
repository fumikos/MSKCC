### Exploratory data analysis (EDA) for CCLE data
### Data downloaded from CCLE website (https://portals.broadinstitute.org/ccle/home)
### EDA on 1) all cells included in a drug sensitivity study (~ 500 cells), or
### 2) top 50 sensitive or resistant (100 cells) is shown below

## Run code in "CCLE_data.R" to load and preprocess the data objects
dim(rearranged)# expression data for the drug of interest
rearranged[1:5,1:5]
dim(CCLE_CELL_all)
head(CCLE_CELL_all)
dim(drug_data)# drug response data (activity area)
head(drug_data)

## Have a look at distributions
X1 <- rearranged
CELL_X1 <- CCLE_CELL_all[match(colnames(X1), CCLE_CELL_all$Cell.line.primary.name),]
library(rafalib)
mypar(1,1)
boxplot(X1)
mypar(2,1)
hist(rowMeans(X1), breaks=50, xlim = c(2,15))
library(genefilter)
hist(rowSds(X1), breaks=50, xlim = c(0,3))

## Heatmap
dim(X1[rowMeans(X1)>5 & rowSds(X1)>1.5,])
levels(CELL_X1$Site.Primary)# 24 levels
cc <- rainbow(24)
sidecols1 <- as.character(CELL_X1$Site.Primary)
for(i in 1:24){
  sidecols1[sidecols1==levels(CELL_X1$Site.Primary)[i]] <- cc[i]
}
mypar(1, 1)
heatmap(X1[rowMeans(X1)>5 & rowSds(X1)>1.5,], labRow = "", labCol = CCLE_CELL_all$Site.Primary, 
        ColSideColors = sidecols1)# Cells generally cluster according to tissue of origin

## Correlation matrix
library(RColorBrewer)
n1 <- ncol(X1)
cors1 <- cor(X1 - rowMeans(X1))
cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
mypar(1,1)
image(cors1, xaxt = "n", yaxt = "n", col = cols, xlab = "", ylab = "", zlim = c(-1,1))

## Principal component analysis (PCA)
# Create labels for top 7 tissues for plotting purpose (below)
Tissues1 <- rep("others", nrow(CELL_X1))
sort(table(CELL_X1$Site.Primary), decreasing = TRUE)
top_tissues1 <- names(sort(table(CELL_X1$Site.Primary), decreasing = TRUE))[1:7]
top_tissues1
for (i in 1:7) {
  Tissues1[CELL_X1$Site.Primary==top_tissues1[i]] <- top_tissues1[i]
}
s1 <- svd(X1 - rowMeans(X1))
PC1_1 <- s1$d[1]*s1$v[,1]
PC2_1 <- s1$d[2]*s1$v[,2]
mypar(1,1)
plot(PC1_1, PC2_1, pch = 21, bg = as.numeric(as.factor(Tissues1)))
legend("topright", levels(factor(Tissues1)), col = seq(along = levels(factor(Tissues1))),
       pch = 15, cex = 0.7, ncol = 1)
# The first PC is mainly associated with tissue type
plot(s1$d^2/sum(s1$d^2), ylab="% variance explained",xlab="Principal component")

## EDA on top 50 drug sensitive or resistant subset of CCLE data
## Assign resistant and sensitive cells
dim(drug_data)
res <- drug_data[order(drug_data[,2])[1:50],]
sen <- drug_data[order(drug_data[,2])[(nrow(drug_data)-49):nrow(drug_data)],]
res_idx <- match(res$Primary.Cell.Line.Name,colnames(CCLE_Exp_minus1))
sen_idx <- match(sen$Primary.Cell.Line.Name,colnames(CCLE_Exp_minus1))
library(rafalib)
mypar(1,2)
plot(sort(drug_data[,2]))
plot(c(res[,2],sen[,2]))
X2 <- CCLE_Exp_minus1[,c(res_idx,sen_idx)]
CELL_X2 <- CCLE_CELL_all[match(colnames(X2), CCLE_CELL_all$Cell.line.primary.name),]

## Look at distributions
mypar(1,1)
boxplot(X2)
mypar(2,1)
hist(rowMeans(X2), breaks=50, xlim = c(2,15))
library(genefilter)
hist(rowSds(X2), breaks=50, xlim = c(0,3))

## Heatmap
dim(X2[rowMeans(X2)>5 & rowSds(X2)>1.5,])
levels(CELL_X2$Site.Primary)# 24 levels
cc <- rainbow(24)
sidecols2 <- as.character(CELL_X2$Site.Primary)
for(i in 1:24){
  sidecols2[sidecols2==levels(CELL_X2$Site.Primary)[i]] <- cc[i]
}
mypar(1,1)
heatmap(X2[rowMeans(X2)>5 & rowSds(X2)>1.5,], labRow = "", labCol = CELL_X2$Site.Primary, 
        ColSideColors = sidecols2)# Cells are clustered according to tissue
heatmap(X2[rowMeans(X2)>5 & rowSds(X2)>1.5,], labRow = "", labCol = CELL_X2$Site.Primary, 
        ColSideColors = c(rep(cc[8],50),rep(cc[20],50)))# res: green, sen: purple

## Correlation matrix
library(RColorBrewer)
n2 <- ncol(X2)
cors2 <- cor(X2 - rowMeans(X2))
cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
mypar(1,1)
image(cors2, xaxt = "n", yaxt = "n", col = cols, xlab = "", ylab = "", zlim = c(-1,1))

## Principal component analysis (PCA)
# Create labels for top 7 tissues for plotting purpose (below)
Tissues2 <- rep("others", nrow(CELL_X2))
sort(table(CELL_X2$Site.Primary), decreasing = TRUE)
top_tissues2 <- names(sort(table(CELL_X2$Site.Primary), decreasing = TRUE))[1:7]
top_tissues2
for (i in 1:7) {
  Tissues2[CELL_X2$Site.Primary==top_tissues2[i]] <- top_tissues2[i]
}
s2 <- svd(X2 - rowMeans(X2))
PC1_2 <- s2$d[1]*s2$v[,1]
PC2_2 <- s2$d[2]*s2$v[,2]
mypar(1,1)
Y <- as.factor(c(rep("res", 50),rep("sen", 50)))
plot(PC1_2, PC2_2, pch = as.numeric(Y) + 20, bg = as.numeric(factor(Tissues2)))
legend("topright", levels(factor(Tissues2)), col = seq(along = levels(factor(Tissues2))),
       pch = 15, cex = 0.7, ncol = 1)
legend("bottomright", levels(Y), pch = seq(along = levels(Y)) + 20, cex = 0.8, ncol = 1)
# The first PC is mainly associated with tissue type
plot(s2$d^2/sum(s2$d^2), ylab="% variance explained",xlab="Principal component")

### Differentially expressed gene analysis using MTP (resampling-based multi hypothesis testing)
### Followed http://www3.nd.edu/~steve/Rcourse/Lecture11v1.pdf
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

## Differentially expressed gene analysis using MTP 
# Prefilter the gene list
library(genefilter)
f1 <- pOverA(0.25, 3.5)# test genes for expression value over 3.5 in 25% of samples
f1fun <- filterfun(f1)
f1genes <- genefilter(X, f1fun)
X1 <- X[f1genes,]
dim(X1)
f2 <- function(x){IQR(x)>0.75}#retain only genes with IQR>0.75
f2fun <- filterfun(f2)
f2genes <- genefilter(X1, f2fun)
X2 <- X1[f2genes,]
dim(X2)
# Run MTP
library(multtest)
Y <- as.factor(c(rep("res", 50),rep("sen", 50)))
mtp <- MTP(X=X2, Y=Y, robust = T, typeone = "fdr", B=300)
mtppvals <- mtp@rawp
hist(mtppvals, breaks = 20)# anti-conservative
plot(mtppvals)# lots of genes are on the y = 0 axis
sum(mtppvals==0)
mtpqvals <- mtp@adjp
hist(mtpqvals, breaks = 20)# conservative
plot(mtpqvals, col=(mtppvals==0)+1)# all significant genes (adjp < 0.25) have pval = 0
mtpfc <- log2(rowMeans(X2[,51:100])/rowMeans(X2[,1:50]))# calculate fold changes
noise <-rnorm(length(mtppvals), 0.0025, 0.0005)# create noise
mtppvals2 <- mtppvals + noise# add noise to rawp for visualization purpose 
noise2 <-rnorm(sum(mtpqvals < 0.25), 0.0005, 0.000125)
mtppvals2[mtpqvals < 0.25] <- noise2# give smaller values for significant genes
plot(mtpfc, -log10(mtppvals2), xlab = "Difference in means", 
     col = as.factor(mtpqvals < 0.25 & abs(mtpfc) > 0.15))
library(hgu133plus2.db)
mtpsymbols <- select(hgu133plus2.db, keys=row.names(X2), columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
mtpsymbols$SYMBOL[!(mtpqvals < 0.25 & abs(mtpfc) > 0.15)] <- ""# mask non-significant genes
text(mtpfc,-log10(mtppvals2),labels=mtpsymbols$SYMBOL,cex=0.5,pos=3)
mtptt <- cbind(mtpsymbols,mtpfc, mtppvals, mtpqvals)
mtptt <- mtptt[mtpqvals < 0.25 & abs(mtpfc) > 0.15,]
dim(mtptt)
mtptt <- mtptt[order(abs(mtptt$mtpfc),decreasing=TRUE),]
mtptt
cc <- rainbow(24)
CELL_X <- CCLE_CELL_all[c(res_idx,sen_idx),]
heatmap(X[mtptt$ENTREZID,], labRow = "", labCol = CELL_X$Site.Primary, 
        ColSideColors = c(rep(cc[8],50),rep(cc[20],50)))# res: green, sen: purple

### Test for differentially expressed genes using the 'limma' and 'sva' packages
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

## Test for differentially expressed genes using limma
## First, fit all genes using a hierarchical model without adjusting 
## http://genomicsclass.github.io/book/pages/hierarchical_models.html
library(limma)
dim(X)
Y <- as.factor(c(rep("res", 50),rep("sen", 50)))
mod <- model.matrix(~ Y)
fit <- lmFit(X, design = mod)
fit <- eBayes(fit)# adjust variance estimates using a hierarchical model
dim(fit$p.value)
head(fit$p.value)
pvals <- fit$p.value[,2]
hist(pvals)# anti-conservative
qvals<-p.adjust(pvals, method = "BH")
hist(qvals)# conservative
fc <- fit$coefficients[,2]
mypar(1,1)
plot(fc, -log10(pvals), xlab = "Difference in means", col = as.factor(qvals < 0.05 & abs(fc) > 0.2))
library(hgu133plus2.db)
symbols <- select(hgu133plus2.db, keys=row.names(fit), columns=c("SYMBOL"), keytype="ENTREZID")
symbols$SYMBOL[!(qvals < 0.05 & abs(fc) > 0.2)] <- ""# mask labels for non significant genes
text(fc, -log10(pvals), labels=symbols$SYMBOL, cex = 0.6, pos = 3)
tt <- topTable(fit, coef=2, number=1000, p.value = 0.05)
tt <- tt[abs(tt$logFC) > 0.2,]
dim(tt)
tt <- tt[order(abs(tt$logFC),decreasing=TRUE),]
gene_table <- select(hgu133plus2.db, keys=row.names(tt), columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
tt <- cbind(gene_table, tt)
head(tt)
cc <- rainbow(24)
CELL_X <- CCLE_CELL_all[c(res_idx,sen_idx),]
heatmap(X[gene_table$ENTREZID[1:50],], labRow = "", labCol = CELL_X$Site.Primary, 
        ColSideColors = c(rep(cc[8],50),rep(cc[20],50)))# res: green, sen: purple

## Next, fit all genes adjusting for tissue type with linear models
## http://genomicsclass.github.io/book/pages/adjusting_with_linear_models.html
mod2 <- model.matrix(~ Y + as.character(CELL_X$Site.Primary))
dim(mod2)
fit2 <- lmFit(X, design = mod2)
fit2 <- eBayes(fit2)# adjust variance estimates using a hierarchical model
dim(fit2$p.value)
head(fit2$p.value)
pvals2 <- fit2$p.value[,2]
hist(pvals2)# anti-conservative
qvals2<-p.adjust(pvals2, method = "BH")
hist(qvals2)# conservative
fc2 <- fit2$coefficients[,2]
mypar(1,1)
plot(fc2, -log10(pvals2), xlab = "Difference in means", col = as.factor(qvals2 < 0.05 & abs(fc2) > 0.2))
library(hgu133plus2.db)
symbols2 <- select(hgu133plus2.db, keys=row.names(fit2), columns=c("SYMBOL"), keytype="ENTREZID")
symbols2$SYMBOL[!(qvals2 < 0.05 & abs(fc2) > 0.2)] <- ""# mask labels for non significant genes
text(fc2, -log10(pvals2), labels=symbols2$SYMBOL, cex = 0.6, pos = 3)
tt2 <- topTable(fit2, coef=2, number=1000, p.value = 0.05)
tt2 <- tt2[abs(tt2$logFC) > 0.2,]
dim(tt2)
tt2 <- tt2[order(abs(tt2$logFC),decreasing=TRUE),]
gene_table2 <- select(hgu133plus2.db, keys=row.names(tt2), columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
tt2 <- cbind(gene_table2, tt2)
head(tt2)

## Finally, fit all genes with adjustment using surrogate variables (SVs)
## http://genomicsclass.github.io/book/pages/adjusting_with_factor_analysis.html
library(sva)
s <- sva(X, mod)
svamod <- model.matrix(~ Y + s$sv)
svafit <- lmFit(X, design = svamod)
svafit <- eBayes(svafit)
svapvals <- svafit$p.value[,2]
hist(svapvals)# anti-conservative
svaqvals<-p.adjust(svapvals, method = "BH")
hist(svaqvals)# conservative
svafc <- svafit$coefficients[,2]
plot(svafc, -log10(svapvals), xlab = "Difference in means", col = as.factor(svaqvals < 0.05 & abs(svafc) > 0.2))
svasymbols <- select(hgu133plus2.db, keys=row.names(svafit), columns=c("SYMBOL"), keytype="ENTREZID")
svasymbols$SYMBOL[!(svaqvals < 0.05 & abs(svafc) > 0.2)] <- ""
text(svafc,-log10(svapvals),labels=svasymbols$SYMBOL,cex=0.6,pos=3)
svatt <- topTable(svafit, coef=2, number=100, p.value = 0.05)
svatt <- svatt[abs(svatt$logFC) > 0.2,]
dim(svatt)
svatt <- svatt[order(abs(svatt$logFC),decreasing=TRUE),]
svagene_table <- select(hgu133plus2.db, keys=row.names(svatt), columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
svatt <- cbind(svagene_table, svatt)
head(svatt)
library(GEOquery)
library(limma)
library(org.Hs.eg.db)
library(GO.db)
source("https://bioconductor.org/biocLite.R") 
biocLite("GO.db")

g <- getGEO("GSE34313")
e <- g[[1]]
e$condition <- e$characteristics_ch1.2
levels(e$condition) <- c("dex24","dex4","control")
table(e$condition)
boxplot(exprs(e), range=0)
names(fData(e))
lvls <- c("control", "dex4")
es <- e[,e$condition %in% lvls]
es$condition <- factor(es$condition, levels=lvls)

design <- model.matrix(~ es$condition)

?lmFit
fit <- lmFit(es, design=design)
?eBayes
fit <- eBayes(fit)
?topTable
tt <- topTable(fit, coef=2, genelist=fData(es)$GENE_SYMBOL)
tt

idx <- grep("GO:0006955", fData(es)$GO_ID)
length(idx)

r1 <- roast(es, idx, design)
?roast
r1

org.Hs.egGO2ALLEGS # a later version than org.Hs.egGO2EG
go2eg <- as.list(org.Hs.egGO2ALLEGS)
head(go2eg)

govector <- unlist(go2eg)
head(govector)
golengths <- sapply(go2eg, length)
head(fData(es)$GENE)
idxvector <- match(govector, fData(es)$GENE)
head(idxvector)
table(is.na(idxvector))
idx <- split(idxvector, rep(names(go2eg), golengths))
go2eg[[1]]
fData(es)$GENE[idx[[1]]]

idxclean <- lapply(idx, function(x) x[!is.na(x)])
idxlengths <- sapply(idxclean, length)
idxsub <- idxclean[idxlengths > 10]
length(idxsub)

?mroast
r2 <- mroast(es, idxsub, design)
head(r2)
r2 <- r2[order(r2$PValue.Mixed),]
columns(GO.db)
keytypes(GO.db)
GOTERM[[rownames(r2)[1]]]
r2tab <- select(GO.db, keys=rownames(r2)[1:10],
                columns=c("GOID","TERM","DEFINITION"), 
                keytype="GOID")
r2tab[,1:2]

r2 <- r2[order(r2$PValue),]
r2tab <- select(GO.db, keys=rownames(r2)[r2$Direction == "Up"][1:10],
                columns=c("GOID","TERM","DEFINITION"), 
                keytype="GOID")
r2tab[,1:2]

r2tab <- select(GO.db, keys=rownames(r2)[r2$Direction == "Down"][1:5],
                columns=c("GOID","TERM","DEFINITION"), 
                keytype="GOID")
r2tab[,1:2]

# from gene set testing
library(rafalib)
library(GSEABase)
library(devtools)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)
library(sva)
library(limma)
X = sampleInfo$group
mod<-model.matrix(~X)
svafit <- sva(geneExpression,mod)

svaX<-model.matrix(~X+svafit$sv)
lmfit <- lmFit(geneExpression,svaX)
tt<-lmfit$coef[,2]*sqrt(lmfit$df.residual)/(2*lmfit$sigma)
pval<-2*(1-pt(abs(tt),lmfit$df.residual[1]))
qval <- p.adjust(pval,"BH")

gsets <- getGmt("c1.all.v4.0.entrez.gmt") 
length(gsets)
head(names(gsets))
gsets[["chryq11"]]
head(geneIds(gsets[["chryq11"]]))

mapGMT2Affy <- function(object,gsets){
  ann<-annotation(object)
  dbname<-paste(ann,"db",sep=".")
  require(dbname,character.only=TRUE)
  gns<-featureNames(object)
  ##This call may generate warnings
  map<-select(get(dbname), keys=gns,columns=c("ENTREZID", "PROBEID"))
  map<-split(map[,1],map[,2])
  indexes<-sapply(gsets,function(ids){
    gns2<-unlist(map[geneIds(ids)])
    match(gns2,gns)
  })
  names(indexes)<-names(gsets)
  return(indexes)
}
##create an Expression Set
rownames(sampleInfo)<- colnames(geneExpression)
e=ExpressionSet(assay=geneExpression,
                phenoData=AnnotatedDataFrame(sampleInfo),
                annotation="hgfocus")
##can safely ignore the warning
gsids <- mapGMT2Affy(e,gsets) 

tab <- table(ingenset=1:nrow(e) %in% gsids[["chryq11"]],signif=qval<0.05)
chisq.test(tab)$p.val

library(rafalib)
mypar(1,1)
qs<-seq(0,1,len=length(tt)+1)-1/(2*length(tt));qs<-qs[-1]
qqplot(qt(qs,lmfit$df.resid),tt,ylim=c(-10,10),xlab="Theoretical quantiles",ylab="Observed") ##note we leave some of the obvious ones out
abline(0,1)






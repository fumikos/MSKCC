### Exploratory data analysis for 17AAG study in CCLE data

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
colnames(CCLE_Exp_all)[991]#"NCIH292_LUNG.1" is the extra
which(grepl("NCIH292", CCLE_CELL_all$CCLE.name))# 963
which(grepl("NCIH292", colnames(CCLE_Exp_all)))# 963 991
head(CCLE_Exp_all[,c(963,991)])
cor(CCLE_Exp_all[,963],CCLE_Exp_all[,991])# 0.9804931
CCLE_Exp_minus1 <- CCLE_Exp_all[,-991] # remove row 991 (NCIH292_LUNG.1) from CCLE_Exp_all
correct_name <- as.vector(CCLE_CELL_all$CCLE.name[which(is.na(idx0))])
colnames(CCLE_Exp_minus1)[which(is.na(idx0))] <- correct_name
identical(colnames(CCLE_Exp_minus1), as.character(CCLE_CELL_all$CCLE.name))# the two now match

## Simplify the cell names
colnames(CCLE_Exp_minus1) <- CCLE_CELL_all$Cell.line.primary.name
CCLE_Exp_minus1[1:5,1:5]

## Creat a vector containing tissue data
tissue <- CCLE_CELL_all[,3]
length(tissue)
table(tissue)
length(table(tissue))# 24 tissue types
sort(table(tissue), decreasing = TRUE)

## Have a look at distributions
library(rafalib)
mypar(1,1)
boxplot(CCLE_Exp_minus1)
mypar(2,1)
hist(rowMeans(CCLE_Exp_minus1), breaks=50, xlim = c(2,15))
library(genefilter)
hist(rowSds(CCLE_Exp_minus1), breaks=50, xlim = c(0,3))

## Principal component analysis (PCA)
# PCA on all cells
s <- svd(CCLE_Exp_minus1-rowMeans(CCLE_Exp_minus1))
PC1 <- s$d[1]*s$v[,1]
PC2 <- s$d[2]*s$v[,2]
mypar(1,1)
plot(PC1,PC2,pch=21,bg=as.numeric(tissue))
legend("topright",levels(tissue),col=seq(along=levels(tissue)),pch=15,cex=0.6, ncol = 1)
plot(s$d^2/sum(s$d^2),ylab="% variance explained",
     xlab="Principal component")
# PCA on selected tissue (8 most abundant)
tissues_x8 <- tissue%in%c("lung","haematopoietic_and_lymphoid_tissue","central_nervous_system","skin",
                      "large_intestine","breast","ovary","pancreas")
s2 <- svd(CCLE_Exp_minus1[,tissues_x8]-rowMeans(CCLE_Exp_minus1[,tissues_x8]))
PC1 <- s2$d[1]*s2$v[,1]
PC2 <- s2$d[2]*s2$v[,2]
mypar(1,1)
plot(PC1,PC2,pch=21,bg=as.numeric(factor(tissue[tissues_x8])))
legend("bottomright",levels(factor(tissue[tissues_x8])),col=seq(along=levels(factor(tissue[tissues_x8]))),pch=15,cex=0.6)
plot(s2$d^2/sum(s2$d^2),ylab="% variance explained",
     xlab="Principal component")
# PCA on gliomas, breast, and blood Ca
gbb<-tissue%in%c("breast","central_nervous_system","haematopoietic_and_lymphoid_tissue")
s3 <- svd(CCLE_Exp_minus1[,gbb]-rowMeans(CCLE_Exp_minus1[,gbb]))
PC1 <- s3$d[1]*s3$v[,1]
PC2 <- s3$d[2]*s3$v[,2]
mypar(1,1)
plot(PC1,PC2,pch=21,bg=as.numeric(factor(tissue[gbb])))
legend("bottomright",levels(factor(tissue[gbb])),col=seq(along=levels(factor(tissue[gbb]))),pch=15,cex=0.6)
text(PC1,PC2,labels=CCLE_CELL_all$Cell.line.primary.name[gbb],cex=0.5,pos=3)
# PCA on gliomas only
gliomas<-tissue%in%c("central_nervous_system")
s4 <- svd(CCLE_Exp_minus1[,gliomas]-rowMeans(CCLE_Exp_minus1[,gliomas]))
PC1 <- s4$d[1]*s4$v[,1]
PC2 <- s4$d[2]*s4$v[,2]
mypar(1,1)
plot(PC1,PC2,pch=21,bg=as.numeric(factor(tissue[gliomas])))
legend("bottomright",levels(factor(tissue[gliomas])),col=seq(along=levels(factor(tissue[gliomas]))),pch=15,cex=0.6)
text(PC1,PC2,labels=CCLE_CELL_all$Cell.line.primary.name[gliomas],cex=0.5,pos=3)
# PCA on all cells that have 17AAG data
CCLE_Drug_all[1:10,]
HSP90_index = which(CCLE_Drug_all$Target=="HSP90")
length(HSP90_index) #503 cells were tested with 17AAG
idx = match(CCLE_Drug_all[HSP90_index,]$Primary.Cell.Line.Name, colnames(CCLE_Exp_minus1))
sum(!is.na(idx)) 
idx <- idx[!is.na(idx)]
length(idx) # of the 503 cells tested with 17AAG, 492 has microarray data
s5 <- svd(CCLE_Exp_minus1[,idx]-rowMeans(CCLE_Exp_minus1[,idx]))
PC1 <- s5$d[1]*s5$v[,1]
PC2 <- s5$d[2]*s5$v[,2]
mypar(1,1)
plot(PC1,PC2,pch=21,bg=as.numeric(factor(tissue[idx])))
legend("bottomright",levels(factor(tissue[idx])),col=seq(along=levels(factor(tissue[idx]))),pch=15,cex=0.6)

## More boxplots
# Boxplot for cells with 17AAG data
boxplot(CCLE_Exp_minus1[,idx], range=0)
# Boxplot for glioma cells
boxplot(CCLE_Exp_minus1[,gliomas], range=0)

## Take a look at NQO1 (1728_at) expression
HSP90_index = which(CCLE_Drug_all$Target=="HSP90")
HSP90 = CCLE_Drug_all[HSP90_index,c(2,13)]
idx2 = match(colnames(CCLE_Exp_minus1)[idx], HSP90$Primary.Cell.Line.Name)# rownums for HSP90
length(idx2)# 492 responses
e_1728 <- CCLE_Exp_minus1["1728_at",idx]
head(e_1728)
head(HSP90[idx2,])
qqnorm(e_1728)
qqline(e_1728)
plot(e_1728,HSP90[idx2,2])
cor(e_1728,HSP90[idx2,2])# 0.2642656




source("http://bioconductor.org/biocLite.R")
biocLite("mixtools")
biocLite("bcellViper")
biocLite("viper")
#biocLite("aracne.networks")# not working (4/13/17)->download 1.0.0 tar ball

library(viper)
data(bcellViper, package="bcellViper")#loads 2 objects, "dset" and "regulon"
class(dset)
class(regulon)
dset
dim(dset)
exprs(dset)[1:5,1:5]
head(pData(dset))
regulon
dim(regulon)
regulon[1]
length(regulon)

#Generating the regulon object
adjfile <- system.file("aracne", "bcellaracne.adj", package = "bcellViper")
adjfile#Full path for ARACNe network file
?aracne2regulon
regul <- aracne2regulon(adjfile, dset, verbose = T)
regul
regul[1]
length(regul)

#Generating the gene expression signatures
?rowTtest
signature <- rowTtest(dset, "description", c("CB", "CC"), "N")
hist(signature$p.value)
signature <- (qnorm(signature$p.value/2, lower.tail = FALSE)*sign(signature$statistic))[, 1]
hist(signature)

#Generating the NULL model by sample permutations
?ttestNull
nullmodel <- ttestNull(dset, "description", c("CB", "CC"), "N", per = 1000, repos = TRUE, verbose = FALSE)
dim(nullmodel)
nullmodel[1:5,1:5]

#msVIPER
regulon
?msviper
mrs <- msviper(signature, regulon, nullmodel)
summary(mrs)
plot(mrs, cex=.7)

#Leading-edge analysis
mrs <- ledge(mrs)
summary(mrs)

#Bootstrap msVIPER
signature <- bootstrapTtest(dset, "description", c("CB", "CC"), "N", verbose = FALSE)
dim(signature)
signature[1:5,1:5]
mrs <- msviper(signature, regulon, nullmodel, verbose = FALSE)
plot(mrs, cex = .7)
?bootstrapmsviper
mrs <- bootstrapmsviper(mrs, "mode")
plot(mrs, cex = .7)

#Shadow analysis
mrshadow <- shadow(mrs, regulators = 25, verbose = FALSE)
summary(mrshadow)
plot(mrshadow, cex = .7)

#Synergy analysis
mrs <- msviperCombinatorial(mrs, regulators = 25, verbose = FALSE)
mrs <- msviperSynergy(mrs, verbose = FALSE)
summary(mrs)
plot(mrs, 25, cex = .7)

#VIPER
vpres <- viper(dset, regulon, verbose = FALSE)
dim(vpres)
exprs(vpres)[1:5,1:5]
dim(dset)
exprs(dset)[1:5,1:5]
tmp <- rowTtest(vpres, "description", c("CB", "CC"), "N")
length(tmp)
hist(tmp$p.value)
data.frame(Gene = rownames(tmp$p.value), t = round(tmp$statistic, 2), "p-value" = signif(tmp$p.value, 3))[order(tmp$p.value)[1:10], ]

#Running VIPER with a null model
?viperSignature
vpsig <- viperSignature(dset, "description", "N", verbose = FALSE)
#> Warning message:
#> Not enough reference samples to compute null model, all samples will be used
vpsig
vpsig$signature
class(vpsig$signature)
dim(vpsig$signature)
exprs(vpsig$signature)[1:5,1:5]
vpsig$nullmodel
class(vpsig$nullmodel)
dim(vpsig$nullmodel)

vpres <- viper(vpsig, regulon, verbose = FALSE)

pos <- pData(vpres)[["description"]] %in% c("M", "CB", "CC")
d1 <- exprs(vpres)[, pos]
colnames(d1) <- pData(vpres)[["description"]][pos]
dd <- dist(t(d1), method = "euclidean")
heatmap(as.matrix(dd), Rowv = as.dendrogram(hclust(dd, method = "average")), symm = T)
dd <- viperSimilarity(d1)
heatmap(as.matrix(as.dist(dd)), Rowv = as.dendrogram(hclust(as.dist(dd), method = "average")), symm = T)

#practice
library(aracne.networks)
data("regulongbm")
write.regulon(regulongbm,n=10)

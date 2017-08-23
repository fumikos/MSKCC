library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
library(limma)

# Setting up data
pheno = pData(bladderEset)
pheno
class(pheno)
# [1] "data.frame"
edata = exprs(bladderEset)
dim(edata)
class(edata)
# [1] "matrix"
mod = model.matrix(~as.factor(cancer), data=pheno)
mod
mod0 = model.matrix(~1,data=pheno)
mod0

# Estimate batch and other artifacts 
n.sv = num.sv(edata,mod,method="leek")
n.sv
svobj = sva(edata,mod,mod0,n.sv=n.sv)

# calqulate F-test pvals without adjusting for SVs
pValues = f.pvalue(edata,mod,mod0)
hist(pValues)
qValues = p.adjust(pValues,method="BH")
sum(qValues<0.05)/length(qValues)
#[1] 0.6818202

# calqulate F-test pvals adjusting for SVs
modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(edata,modSv,mod0Sv)
hist(pValuesSv)
qValuesSv = p.adjust(pValuesSv,method="BH")
sum(qValuesSv<0.05)/length(qValuesSv)
# [1] 0.6647669

# Adjusting for SVs using limma package
fit = lmFit(edata,modSv)
contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),"C2"=c(0,-1,1,rep(0,svobj$n.sv)),"C3"=c(-1,0,1,rep(0,svobj$n.sv)))
fitContrasts = contrasts.fit(fit,contrast.matrix)
eb = eBayes(fitContrasts)
topTableF(eb, adjust="BH")

# Applying the ComBat function
batch = pheno$batch
class(batch)
class(edata)
edata[1:5,1:5]
modcombat = model.matrix(~1, data=pheno)
class(modcombat)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=FALSE, prior.plots=TRUE)
pValuesComBat = f.pvalue(combat_edata,mod,mod0)
library(rafalib)
mypar(1,1)
hist(pValuesComBat)
qValuesComBat = p.adjust(pValuesComBat,method="BH")
sum(qValuesComBat<0.05)/length(qValuesComBat)
# [1] 0.4963425

# Removing known batch effects with a linear model
modBatch = model.matrix(~as.factor(cancer) + as.factor(batch),data=pheno)
modBatch
dim(modBatch)
mod0Batch = model.matrix(~as.factor(batch),data=pheno)
mod0Batch
pValuesBatch = f.pvalue(edata,modBatch,mod0Batch)
hist(pValuesBatch)
qValuesBatch = p.adjust(pValuesBatch,method="BH")
sum(qValuesBatch<0.05)/length(qValuesBatch)
# [1] 0.6635552

# Variance fltering to speed computations when the number of features is large (m > 100,000)
n.sv = num.sv(edata,mod,vfilter=2000,method="leek")
svobj = sva(edata,mod,mod0,n.sv=n.sv,vfilter=2000)

# Applying the fsva function to remove batch effects for prediction
set.seed(12354)
trainIndicator = sample(1:57,size=30,replace=FALSE)
testIndicator = (1:57)[-trainIndicator]
trainData = edata[,trainIndicator]
testData = edata[,testIndicator]
trainPheno = pheno[trainIndicator,]
testPheno = pheno[testIndicator,]

mydata = list(x=trainData,y=trainPheno$cancer)
mytrain = pamr.train(mydata)

table(pamr.predict(mytrain,testData,threshold=2),testPheno$cancer)

trainMod = model.matrix(~cancer,data=trainPheno)
trainMod0 = model.matrix(~1,data=trainPheno)
trainSv = sva(trainData,trainMod,trainMod0)
?sva
fsvaobj = fsva(trainData,trainMod,trainSv,testData)
?fsva
mydataSv = list(x=fsvaobj$db,y=trainPheno$cancer)
mytrainSv = pamr.train(mydataSv)

table(pamr.predict(mytrainSv,fsvaobj$new,threshold=1),testPheno$cancer)

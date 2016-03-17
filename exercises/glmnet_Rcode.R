# http://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
# 100 samples, 20 predictors

library(devtools)
library(glmnet)
load("QuickStartExample.RData")
head(x)
head(y)
fit = glmnet(x, y)
plot(fit, label=T)
class(x)
class(y)
print(fit)
coef(fit,s=0.1)
nx = matrix(rnorm(10*20),10,20)
predict(fit,newx=nx,s=c(0.1,0.05))
cvfit = cv.glmnet(x, y)
plot(cvfit)
cvfit$lambda.min
coef(cvfit, s = "lambda.min")
w = c(rep(1,50),rep(2,50))
w
fit = glmnet(x, y, alpha = 0.2, weights = c(rep(1,50),rep(2,50)), nlambda = 20)
print(fit)
plot(fit, xvar = "lambda", label = TRUE)
plot(fit, xvar = "dev", label = TRUE)
any(fit$lambda == 0.5)
coef.exact = coef(fit, s = 0.5, exact = TRUE)
coef.apprx = coef(fit, s = 0.5, exact = FALSE)
cbind2(coef.exact, coef.apprx)
predict(fit, newx = x[1:5,], type = "response", s = 0.05)
cvfit = cv.glmnet(x, y, type.measure = "mse", nfolds = 20)

library(parallel)
parallel:::detectCores()

foldid=sample(1:10,size=length(y),replace=TRUE)
foldid
table(foldid)

foldid=sample(rep(seq(10),length=length(y)))

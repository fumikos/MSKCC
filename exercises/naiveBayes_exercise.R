library(e1071)
library(pROC)
install.packages("mlbench")
data(HouseVotes84, package = "mlbench")
head(HouseVotes84)

## Categorical data only:
model <- naiveBayes(Class ~ ., data = HouseVotes84)
head(model)
dim(HouseVotes84)
predict(model, HouseVotes84[1:10,])
predict(model, HouseVotes84[1:10,], type = "raw")

pred <- predict(model, HouseVotes84)
pred2 <- predict(model, HouseVotes84[,-1])
table(pred, HouseVotes84$Class)
mean(pred==HouseVotes84[,1])
mean(pred==pred2)

## using laplace smoothing:
model <- naiveBayes(Class ~ ., data = HouseVotes84, laplace = 3)
pred <- predict(model, HouseVotes84[,-1])
table(pred, HouseVotes84$Class)
mean(pred==HouseVotes84$Class)

## Example of using a contingency table:
data(Titanic)
m <- naiveBayes(Survived ~ ., data = Titanic)
m
predict(m, as.data.frame(Titanic))

## Example with metric predictors:
data(iris)
m <- naiveBayes(Species ~ ., data = iris)
predict(m, as.data.frame(iris))

# ROC analysis
data(aSAH)
# Basic example
roc(aSAH$outcome, aSAH$s100b, levels=c("Good", "Poor"))
roc(aSAH$outcome, aSAH$s100b, levels=c("Poor", "Good"))

# Plot and CI (see plot.roc and ci for more options):
roc(aSAH$outcome, aSAH$s100b, percent=TRUE, plot=TRUE, ci=TRUE)

# Smoothed ROC curve
roc(aSAH$outcome, aSAH$s100b, smooth=TRUE, plot=TRUE, ci=TRUE)









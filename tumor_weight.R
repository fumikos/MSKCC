library(genefilter)

data=read.csv("2-1-2015.csv")
data
data=data[data$L.2.W.2!=0, c(1,5)]
names<-data$ID
names
data
data=data[,2]
data
names(data)<-names
data
data=sort(data)
data
S<-data[1:6]
M<-data[7:12]
L<-data[13:18]
S
M
L
groups<-matrix(data = NA, nrow = 6, ncol = 3)
rownames(groups)<-c("G1", "G2","G3","G4","G5","G6")
groups
seeds<-seq(0, 9, 1)
results<-matrix(data = NA, nrow = 10, ncol = 2)
for(seed in seeds){
  set.seed(seed*10)
  groups[,1]<-sample(S, 6, replace = F)
  set.seed(seed*10+1)
  groups[,2]<-sample(M, 6, replace = F)
  set.seed(seed*10+2)
  groups[,3]<-sample(L, 6, replace = F)
  sd(rowMeans(groups))
  sd(rowSds(groups))
  results[seed+1,]<-c(sd(rowMeans(groups)),sd(rowSds(groups)))
}
results[,1]+results[,2]
colMeans(results)
results_norm<-results
results_norm[,1]<-results[,1]/14.64185
results_norm[,2]<-results[,2]/18.56203
colMeans(results_norm)
results_norm[,1]+results_norm[,2]


S_ID<-names(data)[1:6]
M_ID<-names(data)[7:12]
L_ID<-names(data)[13:18]
set.seed(60)
groups[,1]<-sample(S_ID, 6, replace = F)
set.seed(61)
groups[,2]<-sample(M_ID, 6, replace = F)
set.seed(62)
groups[,3]<-sample(L_ID, 6, replace = F)
sd(rowMeans(groups))
sd(rowSds(groups))
groups
groups<-t(groups)
groups
write.csv(groups, "groups.csv")

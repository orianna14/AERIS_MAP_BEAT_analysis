# This script performs community type analysis using two different methods
library(textmineR)
library(clValid)
library(cluster) 
library(clusterSim)
library(fpc)
library(ade4)

#Method 1: hclust based on JSD and clValid calculation of Silhouette

data<-read.table("L6_batchrm.txt",header=TRUE,row.names=1,check.names=F,sep="\t")
jsd<-CalcJSDivergence(as.matrix(data))
plot(hclust(dist(jsd), method="ward"))
# use a local version to allow for JSD metrics
clValid<-(as.matrix(data),2:10,clMethods="hierarchical",validation="internal",metric="JSD",method="ward",neighbSize=10)

#Method 2: PAM clustering based on JSD and calculation of CH index and Silhouette

data=read.table("L6_batchrm.txt",header=T,row.names=1,dec=".",sep="\t",check.names=F)
data=data[-1,]
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}
data.dist=dist.JSD(data)
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}
nclusters=NULL
for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist,
                          centrotypes = "medoids")
  }
}
plot(nclusters, type="h", xlab="k clusters", ylab="CH index",main="Optimal number of clusters")

data.cluster<-pam.clustering(data.dist,2)
obs.silhouette=mean(silhouette(data.cluster, data.dist)[,2])
cat(obs.silhouette) 
silhouette(pam(data.dist,2))
obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F,nf=k-1) 
dev.new()
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=T,sub="Between-class analysis")

obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
dev.new()
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F)
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(3,2,4))

clust<-pam(data.dist,2)
write.table(clust$clustering,"clusters.txt",sep="\t",append=FALSE,quote=FALSE)
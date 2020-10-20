# This script performs Combat analysis for microbiome data
library(sva)
library(funrar)

#read in genus-level microbiome profile
data<-read.table("L6.txt",header=TRUE,row.names=1,check.names=F,sep="\t")
data=log(data+0.5*min(data),base=2) # data transformation

#read in metadata file with batch effect information in the 'batch' column
meta<-read.table("all_metadata.txt",header=TRUE,row.names=1,check.name=F,sep="\t")
batch<-meta$batch

combat<-ComBat(as.matrix(data),batch,mod=NULL,par.prior=TRUE,prior.plots=FALSE)
## exponential and total sum scaling tranformation convert back to relative abundance
combat_rel<-make_relative(2^^(combat-max(combat)))

write.table(combat_rel,file="L6_batchrm.txt",append=FALSE,quote=FALSE,sep="\t")

# optional: perform PCA to assess removal of batch effects
pca<-prcomp(t(all_combat))
write.table(pca$x,file="pca.txt",append=FALSE,quote=FALSE,sep="\t")
pca<-prcomp(t(data))
write.table(pca$x,file="pca2.txt",append=FALSE,quote=FALSE,sep="\t")

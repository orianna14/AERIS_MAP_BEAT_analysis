args=commandArgs(trailingOnly=TRUE)
data<-read.table(args[1],row.names=1,sep="\t",header=T)
cor<-ccf(data$Abundance,data$Neutro,type="correlation")
capture.output(cor,file=args[2],append=FALSE)


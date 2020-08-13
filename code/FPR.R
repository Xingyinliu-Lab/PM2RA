#setwd("/home/liuZ/projects/1-PM/COMBO")
###generating synthetic dataset by randomly sampling samples into case or control groups
data=read.table("table.from_biom.txt",head=T,row.names=1)
fil=apply(data,1,function(x){length(x[which(x>0)])/100})
comp=apply(data,2,function(x){x/sum(x)})
df=comp[which(fil>=0.5),]
df2=data[which(fil>=0.5),]
out=NULL
for(i in seq(1,100)){
samp=sample(seq(1,100),50,replace=FALSE)
dif=setdiff(seq(1,100),samp)
lis=c(samp,dif)
out=rbind(out,lis)
case=df2[,samp];contr=df2[,-samp]
write.table(case,paste(i,".sampling.case.txt",sep=""),col.names=T,quote=F,sep="\t")
write.table(contr,paste(i,".sampling.control.txt",sep=""),col.names=T,quote=F,sep="\t")
case=df[,samp];contr=df[,-samp]
tc=t(case)
tc=data.frame("cases",tc)
tt=t(contr)
tt=data.frame("control",tt)
colnames(tt)=c("lable",seq(1,123))
colnames(tc)=c("lable",seq(1,123))
md=rbind(tc,tt)
colnames(md)=c("label",row.names(case))
write.table(md,paste(i,".table.4.PM2RA.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
}
write.table(out,"sampling.record.txt",row.names=F,col.names=F,quote=F,sep="\t") ##50VS50 samples
###generating synthetic dataset by randomly sampling samples into case or control groups


####Detecting association alteration

###PM2RA
sh pm.sh

####co-occurrence-based methods
args=command(TRUE)
input<-args[1]
library(SpiecEasi)
library(Matrix)
library(reshape2)
data=read.table(paste(input,".sampling.case.txt",sep=""),head=T,row.names=1)
mtdata<-matrix(as.numeric(unlist(data)),nrow = length(rownames(data)),ncol = length(colnames(data)))
rownames(mtdata)=rownames(data)
colnames(mtdata)=colnames(data)
case=t(mtdata)
se.mb.case <- spiec.easi(case, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
se.gl.case <- spiec.easi(case, method='glasso', lambda.min.ratio=1e-2,nlambda=20, pulsar.params=list(rep.num=50))
sparcc.case <- sparcc(case)
sparcc.graph <- abs(sparcc.case$Cor) >= 0.3
diag(sparcc.graph) <- 0
sparcc.case <- Matrix(sparcc.graph, sparse=TRUE)
mb.case<-getRefit(se.mb.case)
gl.case<-getRefit(se.gl.case)
mb.case<-matrix(as.numeric(unlist(mb.case)),nrow = 123,ncol = 123)
rownames(mb.case)=colnames(case)
colnames(mb.case)=colnames(case)
mb.rcase=melt(mb.case)
gl.case<-matrix(as.numeric(unlist(gl.case)),nrow = 123,ncol = 123)
rownames(gl.case)=colnames(case)
colnames(gl.case)=colnames(case)
gl.rcase=melt(gl.case)
sparcc.case<-matrix(as.numeric(unlist(sparcc.case)),nrow = 123,ncol = 123)
rownames(sparcc.case)=colnames(case)
colnames(sparcc.case)=colnames(case)
sparcc.rcase=melt(sparcc.case)

name=apply(mb.rcase,1,function(x){paste(max(x[1],x[2]),min(x[1],x[2]),sep=":")})
ind=which(duplicated(name))
mb.rcase=mb.rcase[-ind,];gl.rcase=gl.rcase[-ind,];sparcc.rcase=sparcc.rcase[-ind,]


data=read.table(paste(input,".sampling.control.txt",sep=""),head=T,row.names=1)
mtdata<-matrix(as.numeric(unlist(data)),nrow = length(rownames(data)),ncol = length(colnames(data)))
rownames(mtdata)=rownames(data)
colnames(mtdata)=colnames(data)
ncase=t(mtdata)
se.mb.ncase <- spiec.easi(ncase, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
se.gl.ncase <- spiec.easi(ncase, method='glasso', lambda.min.ratio=1e-2,nlambda=20, pulsar.params=list(rep.num=50))
sparcc.ncase <- sparcc(ncase)
sparcc.graph <- abs(sparcc.ncase$Cor) >= 0.3
diag(sparcc.graph) <- 0
sparcc.ncase <- Matrix(sparcc.graph, sparse=TRUE)
mb.ncase<-getRefit(se.mb.ncase)
gl.ncase<-getRefit(se.gl.ncase)
rownames(mb.ncase)=colnames(ncase)
colnames(mb.ncase)=colnames(ncase)
mb.ncase<-matrix(as.numeric(unlist(mb.ncase)),nrow = 123,ncol = 123)
rownames(mb.ncase)=colnames(ncase)
colnames(mb.ncase)=colnames(ncase)
mb.rncase=melt(mb.ncase)

gl.ncase<-matrix(as.numeric(unlist(gl.ncase)),nrow = 123,ncol = 123)
rownames(gl.ncase)=colnames(ncase)
colnames(gl.ncase)=colnames(ncase)
gl.rncase=melt(gl.ncase)

sparcc.ncase<-matrix(as.numeric(unlist(sparcc.ncase)),nrow = 123,ncol = 123)
rownames(sparcc.ncase)=colnames(ncase)
colnames(sparcc.ncase)=colnames(ncase)
sparcc.rncase=melt(sparcc.ncase)

name=apply(mb.rncase,1,function(x){paste(max(x[1],x[2]),min(x[1],x[2]),sep=":")})
ind=which(duplicated(name))
mb.rncase=mb.rncase[-ind,];gl.rncase=gl.rncase[-ind,];sparcc.rncase=sparcc.rncase[-ind,]
mb=cbind(mb.rcase,mb.rncase)
mb.diff=mb[which(mb[,3] != mb[,6]),]
gls=cbind(gl.rcase,gl.rncase)
gl.diff=gls[which(gls[,3] != gls[,6]),]
sparcc=cbind(sparcc.rcase,sparcc.rncase)
sparcc.diff=sparcc[which(sparcc[,3] != sparcc[,6]),]
out1=paste(i,".case.mb.res.txt",sep="")
out2=paste(i,".case.gl.res.txt",sep="")
out3=paste(i,".case.sparcc.res.txt",sep="")
out11=paste(i,".control.mb.res.txt",sep="")
out22=paste(i,".control.gl.res.txt",sep="")
out33=paste(i,".control.sparcc.res.txt",sep="")
out111=paste(i,".mb.diff.txt",sep="")
out222=paste(i,".gl.diff.txt",sep="")

out333=paste(i,".sparcc.diff.txt",sep="")
write.table(mb.rcase,out1,quote=F,sep="\t")
write.table(gl.rcase,out2,quote=F,sep="\t")
write.table(sparcc.rcase,out3,quote=F,sep="\t")
write.table(mb.rncase,out11,quote=F,sep="\t")
write.table(gl.rncase,out22,quote=F,sep="\t")
write.table(sparcc.rncase,out33,quote=F,sep="\t")
write.table(mb.diff,out111,quote=F,sep="\t")
write.table(gl.diff,out222,quote=F,sep="\t")
write.table(sparcc.diff,out333,quote=F,sep="\t")
####Detecting association alteration

########result processing######
#1. PM2RA
out=NULL
for(i in seq(1,100)){
	data=read.table(paste(i,".iteration.2Dscan.res.txt",sep=""),head=T,sep="\t")
	sig=data[which(data$FDR < 0.05),]
	g=dim(sig)[1]
	out=rbind(out,c(i,g))
}
write.table(out,"100.sampling.pm2ca.res.stat.txt",quote=F,sep="\t")

#2. co-occurrence-based methods	
out=NULL
for(i in seq(1,100)){
	in1=paste(i,"gl.diff.txt",sep=".")
	in2=paste(i,"mb.diff.txt",sep=".")
	in3=paste(i,"sparcc.diff.txt",sep=".")
	glf=read.table(in1,sep="\t",head=T)
	mb=read.table(in2,sep="\t",head=T)
	sparcc=read.table(in3,sep="\t",head=T)
	sg=dim(glf)[1]
	sm=dim(mb)[1]
	ss=dim(sparcc)[1]
	out=rbind(out,c(i,sg,sm,ss))
}
write.table(out,"100.sampling.net.res.stat.txt",quote=F,sep="\t")
########result processing######











	


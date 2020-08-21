data=read.table("/home/liuZ/projects/1-PM/COMBO/table.from_biom.txt",head=T,row.names=1)
fil=apply(data,1,function(x){length(x[which(x>0)])/100})
df=data[which(fil>=0.8),]
comp=apply(df,2,function(x){x/sum(x)})
td=t(comp)
write.table(td,"biom.filter.txt",quote=F,sep="\t")

###calcuate pairwise correlation and difference in abundance
for(i in seq(1,36)){
	for(j in seq(i+1,37)){
			pv=cor.test(td[,i],td[,j],method="spearman")$p.value
			cr=cor.test(td[,i],td[,j],method="spearman")$estimate
			pdiff=wilcox.test(td[,i],td[,j])$p.value
			fd=log2(sum(td[,i])/sum(td[,j]))
			out=rbind(out,c(i,j,pv,cr,pdiff,fd))
	}
}
sig=out
write.table(out,"simu.cor.diff.pvalue.txt",row.names=F,col.names=F,sep="\t",quote=F)
###calcuate pairwise correlation and difference in abundance

###generating synthetic dataset by exchange each pair of the OTU
ct=td
df=out
for(i in seq(1,666)){
	c1=df[i,1]
	c2=df[i,2]
	temp=td
	temp[,c1]=td[,c2]
	temp[,c2]=td[,c1]
	colnames(temp)=colnames(ct)
	df=rbind(temp,ct)
	Disease=c(rep("cases",100),rep("control",100))
	md=cbind(Disease,df)
	write.table(md,paste(i,".table.4.PM2RA.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
	write.table(t(temp),paste(i,".sampling.case.txt",sep=""),col.names=T,quote=F,sep="\t")
	write.table(t(ct),paste(i,".sampling.control.txt",sep=""),col.names=T,quote=F,sep="\t")
}
###generating synthetic dataset by exchange each pair of the OTU

####Detect association alteration
1. run PM2RA
2. run other methods:

library(SpiecEasi)
library(Matrix)
library(reshape2)
for(i in seq(1,666)){
data=read.table(paste(i,".sampling.case.txt",sep=""),head=T,row.names=1)
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
mb.case<-matrix(as.numeric(unlist(mb.case)),nrow = 37,ncol = 37)
rownames(mb.case)=colnames(case)
colnames(mb.case)=colnames(case)
mb.rcase=melt(mb.case)
gl.case<-matrix(as.numeric(unlist(gl.case)),nrow = 37,ncol = 37)
rownames(gl.case)=colnames(case)
colnames(gl.case)=colnames(case)
gl.rcase=melt(gl.case)
sparcc.case<-matrix(as.numeric(unlist(sparcc.case)),nrow = 37,ncol = 37)
rownames(sparcc.case)=colnames(case)
colnames(sparcc.case)=colnames(case)
sparcc.rcase=melt(sparcc.case)

name=apply(mb.rcase,1,function(x){paste(max(x[1],x[2]),min(x[1],x[2]),sep=":")})
ind=which(duplicated(name))
mb.rcase=mb.rcase[-ind,];gl.rcase=gl.rcase[-ind,];sparcc.rcase=sparcc.rcase[-ind,]

data=read.table(paste(i,".sampling.control.txt",sep=""),head=T,row.names=1)
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
mb.ncase<-matrix(as.numeric(unlist(mb.ncase)),nrow = 37,ncol = 37)
rownames(mb.ncase)=colnames(ncase)
colnames(mb.ncase)=colnames(ncase)
mb.rncase=melt(mb.ncase)

gl.ncase<-matrix(as.numeric(unlist(gl.ncase)),nrow = 37,ncol = 37)
rownames(gl.ncase)=colnames(ncase)
colnames(gl.ncase)=colnames(ncase)
gl.rncase=melt(gl.ncase)

sparcc.ncase<-matrix(as.numeric(unlist(sparcc.ncase)),nrow = 37,ncol = 37)
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

sparcc =cbind(sparcc.rcase,sparcc.rncase)
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
}
####Detect association alteration


###result processing######
#1. PM2RA
data=read.table("simu.cor.diff.pvalue.txt")
cmp=read.table("biom.filter.txt",head=T,row.names=1)
bact=colnames(cmp)
out=NULL
for(i in seq(1,666)){
	c1=data[i,1]
	c2=data[i,2]
	b1=bact[c1];b2=bact[c2]
	input=paste(i,".iteration.2Dscan.res.txt",sep="")
	res=read.table(input,head=T,sep="\t")
	sig=res[which(res$FDR<0.05),]
	sig2=sig
	index=which((sig2$OTU1 %in% c(b1,b2)) & (sig2$OTU2 %in% c(b1,b2)))
	if(length(index>0)){
	sig2=sig2[-index,]
	}
	sig2s=apply(sig2,1,function(x){a=max(as.character(x[1]),as.character(x[2]));b=min(as.character(x[1]),as.character(x[2]));out=c(a,b,x[3],x[4])})
	sig2s=t(sig2s)
	sig2s=sig2s[!duplicated(sig2s),]
	s=dim(sig2s)[1]
	FNR=(70-s)/70
	out=rbind(out,c(i,b1,b2,FNR))
}
write.table(out,"PM2RA.FNR.res.txt",row.names=F,quote=F,sep="\t")

ab=read.table("simu.cor.diff.pvalue.name.txt",head=TRUE)
res=read.table("PM2RA.FNR.all.res.txt")
res=res[!duplicated(res),]
res$name=apply(res,1,function(x){a=max(as.character(x[2]),as.character(x[3]));b=min(as.character(x[2]),as.character(x[3]));paste(a,b,sep=":")})

ind1=which((ab$dfdr<0.05)&(abs(ab$V6)>1)) ###Difference in abundance
ind2=which((ab$rfdr>0.05)|(ab$V4 <0.2)) ###not correlated
index2=intersect(ind1,ind2)
fil=ab[index2,]
res1=res[which(res$name %in% fil$name),]
mean(res1$V4)

fil=ab[ind1,] #Difference in abundance
res2=res[which(res$name %in% fil$name),]
mean(res2$V4)

fil=ab[ind2,] ###not correlated
res3=res[which(res$name %in% fil$name),]
mean(res3$V4) 


FNR=c(res1$V4,res2$V4,res3$V4)
label=c(rep("A",length(res1$V4)),rep("B",length(res2$V4)),rep("C",length(res3$V4)))
df=data.frame(label,FNR)
write.table(df,"PM2RA.FNR.stat.txt",quote=F,sep="\t",row.names=F)

#2. Co-coccurrence-based methods
data=read.table("simu.cor.diff.pvalue.name.txt",head=T)
cmp=read.table("biom.filter.txt",head=T,row.names=1)
bact=colnames(cmp)
sparcc=NULL
for(i in seq(1,666)){
	c1=data[i,1]
	c2=data[i,2]
	b1=bact[c1];b2=bact[c2]
	input=paste(i,".sparcc.diff.txt",sep="")
	res=read.table(input,head=T,sep="\t")
	s=dim(res)[1]
	FNR=(70-s)/70
	sparcc=rbind(sparcc,c(i,b1,b2,FNR))
}
write.table(sparcc,"Sparcc.FNR.res.txt",row.names=F,col.names=F,quote=F,sep="\t")

gl=NULL
for(i in seq(1,666)){
	c1=data[i,1]
	c2=data[i,2]
	b1=bact[c1];b2=bact[c2]
	input=paste(i,".gl.diff.txt",sep="")
	res=read.table(input,head=T,sep="\t")
	s=dim(res)[1]
	FNR=(70-s)/70
	gl=rbind(gl,c(i,b1,b2,FNR))
}
write.table(gl,"gl.FNR.res.txt",row.names=F,col.names=F,quote=F,sep="\t")

mb=NULL
for(i in seq(1,666)){
	c1=data[i,1]
	c2=data[i,2]
	b1=bact[c1];b2=bact[c2]
	input=paste(i,".mb.diff.txt",sep="")
	res=read.table(input,head=T,sep="\t")
	s=dim(res)[1]
	FNR=(70-s)/70
	mb=rbind(mb,c(i,b1,b2,FNR))
}
write.table(mb,"mb.FNR.res.txt",row.names=F,col.names=F,quote=F,sep="\t")

###result processing######


####clr transformation
library(psych)
data=read.table("table.from_biom.txt",head=T,row.names=1)
fil=apply(data,1,function(x){length(x[which(x>0)])/100})
df=data[which(fil>=0.8),]
df=df+1
gm<-apply(df,2,function(x){log2(x/geometric.mean(x))})





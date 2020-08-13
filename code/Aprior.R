library(arules)
library(arulesViz)

data=read.table("crc_baxter.CRC.H.id.match.ready.txt",head=T,row.names=1)
a=apply(data,1,function(x){length(x[which(x>0)])/length(x)})
data=data[which(a>0.1),]
##case
case=data[,seq(1,120)]
mtdata<-matrix(as.numeric(unlist(case)),nrow = length(rownames(case)),ncol = length(colnames(case)))
rownames(mtdata)=rownames(case)
colnames(mtdata)=colnames(case)
case=mtdata
case[which(case>0)] <-1
case[which(case ==0)] <-0
case=t(case)
a=case[1,]
b=names(a[which(a==1)])
MyList = list()
for( i in seq(1,dim(case)[1])){
	a=case[i,]
	b=names(a[which(a==1)])
	MyList[[i]]=b
	}
	
MyTrans<-as(MyList,"transactions")
MyRules<-apriori(data=MyTrans,parameter=list(support=0.6,confidence=0.99,target="rules"))
MyRules.sorted<-sort(x=MyRules,by="lift",decreasing=TRUE)
a<-is.significant(MyRules.sorted, MyTrans, method = "fisher", alpha = 0.01, adjust = "fdr")
rules=MyRules.sorted[which(a=="TRUE"),]
write(rules, file = "CRC.H.s06c099.csv", sep = ",")
subrules<-rules[seq(1,10)]
write(subrules, file = "s06c099.top10.csv", sep = ",")

pdf("CRC.H.rules.top10.pdf")
for(i in seq(1,10)){
plot(subrules[i], method = "graph", control = list(verbose = TRUE))
}
dev.off()

##case
case=data[,seq(121,292)]
mtdata<-matrix(as.numeric(unlist(case)),nrow = length(rownames(case)),ncol = length(colnames(case)))
rownames(mtdata)=rownames(case)
colnames(mtdata)=colnames(case)
case=mtdata
case[which(case>0)] <-1
case[which(case ==0)] <-0
case=t(case)
a=case[1,]
b=names(a[which(a==1)])
MyList = list()
for( i in seq(1,dim(case)[1])){
	a=case[i,]
	b=names(a[which(a==1)])
	MyList[[i]]=b
	}
	
MyTrans<-as(MyList,"transactions")
MyRules<-apriori(data=MyTrans,parameter=list(support=0.6,confidence=0.99,target="rules"))
MyRules.sorted<-sort(x=MyRules,by="lift",decreasing=TRUE)
a<-is.significant(MyRules.sorted, MyTrans, method = "fisher", alpha = 0.01, adjust = "fdr")
rules=MyRules.sorted[which(a=="TRUE"),]
write(rules, file = "CRC.case.s06c099.csv", sep = ",")
subrules<-rules[seq(1,10)]
write(subrules, file = "CRC.case.s06c099.top10.csv", sep = ",")
pdf("CRC.case.rules.top10.pdf")
for(i in seq(1,10)){
plot(subrules[i], method = "graph")
}
dev.off()

##### 
library(arules)
library(arulesViz)

data=read.table("ob_goodrich.HOB.id.match.ready.txt",head=T,row.names=1)
a=apply(data,1,function(x){length(x[which(x>0)])/length(x)})
data=data[which(a>0.1),]
##case
case=data[,seq(452,644)]
mtdata<-matrix(as.numeric(unlist(case)),nrow = length(rownames(case)),ncol = length(colnames(case)))
rownames(mtdata)=rownames(case)
colnames(mtdata)=colnames(case)
case=mtdata
case[which(case>0)] <-1
case[which(case ==0)] <-0
case=t(case)
a=case[1,]
b=names(a[which(a==1)])
MyList = list()
for( i in seq(1,dim(case)[1])){
	a=case[i,]
	b=names(a[which(a==1)])
	MyList[[i]]=b
	}
	
MyTrans<-as(MyList,"transactions")
MyRules<-apriori(data=MyTrans,parameter=list(support=0.6,confidence=0.99,target="rules"))
MyRules.sorted<-sort(x=MyRules,by="lift",decreasing=TRUE)
a<-is.significant(MyRules.sorted, MyTrans, method = "fisher", alpha = 0.01, adjust = "fdr")
rules=MyRules.sorted[which(a=="TRUE"),]
write(rules, file = "HOB.case.s06c099.csv", sep = ",")
subrules<-rules[seq(1,10)]
write(subrules, file = "HOB.case.s06c099.top10.csv", sep = ",")
pdf("HOB.case.rules.top10.pdf")
for(i in seq(1,10)){
plot(subrules[i], method = "graph")
}
dev.off()

###control
case=data[,seq(1,451)]
mtdata<-matrix(as.numeric(unlist(case)),nrow = length(rownames(case)),ncol = length(colnames(case)))
rownames(mtdata)=rownames(case)
colnames(mtdata)=colnames(case)
case=mtdata
case[which(case>0)] <-1
case[which(case ==0)] <-0
case=t(case)
a=case[1,]
b=names(a[which(a==1)])
MyList = list()
for( i in seq(1,dim(case)[1])){
	a=case[i,]
	b=names(a[which(a==1)])
	MyList[[i]]=b
	}
	
MyTrans<-as(MyList,"transactions")
MyRules<-apriori(data=MyTrans,parameter=list(support=0.6,confidence=0.99,target="rules"))
MyRules.sorted<-sort(x=MyRules,by="lift",decreasing=TRUE)
a<-is.significant(MyRules.sorted, MyTrans, method = "fisher", alpha = 0.01, adjust = "fdr")
rules=MyRules.sorted[which(a=="TRUE"),]
write(rules, file = "HOB.control.s06c099.csv", sep = ",")
subrules<-rules[seq(1,10)]
write(subrules, file = "HOB.control.s06c099.top10.csv", sep = ",")
pdf("HOB.control.rules.top10.pdf")
for(i in seq(1,10)){
plot(subrules[i], method = "graph")
}
dev.off()





############################github
library(arules)
library(arulesViz)
apq<-function(df,name){
case=df
mtdata<-matrix(as.numeric(unlist(case)),nrow = length(rownames(case)),ncol = length(colnames(case)))
rownames(mtdata)=rownames(case)
colnames(mtdata)=colnames(case)
case=mtdata
case[which(case>0.001)] <-1
case[which(case <=0.001)] <-0
#case=t(case)
a=case[1,]
b=names(a[which(a==1)])
MyList = list()
for( i in seq(1,dim(case)[1])){
	a=case[i,]
	b=names(a[which(a==1)])
	MyList[[i]]=b
	}
	
MyTrans<-as(MyList,"transactions")
MyRules<-apriori(data=MyTrans,parameter=list(support=0.6,confidence=0.99,target="rules"))
MyRules.sorted<-sort(x=MyRules,by="lift",decreasing=TRUE)
a<-is.significant(MyRules.sorted, MyTrans, method = "fisher", alpha = 0.05, adjust = "fdr")
rules=MyRules.sorted[which(a=="TRUE"),]
write(rules, file = paste(name,".s06c099.csv",sep=""), sep = ",")
subrules<-rules[seq(1,10)]
write(subrules, file = paste(name,".s06c099.top10.csv",sep=""), sep = ",")
pdf(paste(name,".s06c099.top10.pdf",sep=""))
for(i in seq(1,10)){
plot(subrules[i], method = "graph")
}
dev.off()

}


#####feng
orifile<-'feng.fil.txt'
cd1<-46
cd2<-63
a=1
b=cd1
cc=b+1
d=cd1+cd2
tdata<-read.table(file=orifile,header=T)
rownames(tdata)=seq(1,dim(tdata)[1])
data=t(tdata)
index=apply(data,1,function(x){length(x[which(x>0)])/length(x)})
data=data[which(index>0.1),]
case=t(data[,seq(a,b)])
apq(case,"feng.case")
contrs=t(data[,seq(cc,d)])
apq(contrs,"feng.control")

###Zeller
orifile<-'zeller.fil.txt'
tdata<-read.table(file=orifile,header=T)
rownames(tdata)=seq(1,dim(tdata)[1])
cd1<-88
cd2<-64
a=1
b=cd1
cc=b+1
d=cd1+cd2
data=t(tdata)
index=apply(data,1,function(x){length(x[which(x>0)])/length(x)})
data=data[which(index>0.1),]
case=t(data[,seq(a,b)])
apq(case,"Zeller.case")
contrs=t(data[,seq(cc,d)])
apq(contrs,"Zeller.control")
####

###YU
orifile<-'yu.fil.txt'
tdata<-read.table(file=orifile,header=T)
rownames(tdata)=seq(1,dim(tdata)[1])
cd1<-73
cd2<-92
a=1
b=cd1
cc=b+1
d=cd1+cd2
data=t(tdata)
index=apply(data,1,function(x){length(x[which(x>0)])/length(x)})
data=data[which(index>0.1),]
case=t(data[,seq(a,b)])
apq(case,"yu.case")
contrs=t(data[,seq(cc,d)])
apq(contrs,"yu.control")
###YU

################HOW
library(arules)
library(arulesViz)
apq<-function(df,name){
case=df
mtdata<-matrix(as.numeric(unlist(case)),nrow = length(rownames(case)),ncol = length(colnames(case)))
rownames(mtdata)=rownames(case)
colnames(mtdata)=colnames(case)
case=mtdata
case[which(case>0)] <-1
case[which(case ==0)] <-0
case=t(case)
a=case[1,]
b=names(a[which(a==1)])
MyList = list()
for( i in seq(1,dim(case)[1])){
	a=case[i,]
	b=names(a[which(a==1)])
	MyList[[i]]=b
	}
	
MyTrans<-as(MyList,"transactions")
MyRules<-apriori(data=MyTrans,parameter=list(support=0.6,confidence=0.99,target="rules"))
MyRules.sorted<-sort(x=MyRules,by="lift",decreasing=TRUE)
a<-is.significant(MyRules.sorted, MyTrans, method = "fisher", alpha = 0.01, adjust = "fdr")
rules=MyRules.sorted[which(a=="TRUE"),]
write(rules, file = paste(name,".s06c099.csv",sep=""), sep = ",")
subrules<-rules[seq(1,10)]
write(subrules, file = paste(name,".s06c099.top10.csv",sep=""), sep = ",")
pdf(paste(name,".s06c099.top10.pdf",sep=""))
for(i in seq(1,10)){
plot(subrules[i], method = "graph")
}
dev.off()

}


data=read.table("ob_goodrich.HOW.id.match.ready.txt",head=T,row.names=1)
a=apply(data,1,function(x){length(x[which(x>0)])/length(x)})
data=data[which(a>0.1),]
case=data[,seq(1,451)]
apq(case,"HOW.control")
case=data[,seq(452,786)]
apq(case,"HOW.case")

################pcb


library(arules)
library(arulesViz)
apq<-function(df,name){
case=df
mtdata<-matrix(as.numeric(unlist(case)),nrow = length(rownames(case)),ncol = length(colnames(case)))
rownames(mtdata)=rownames(case)
colnames(mtdata)=colnames(case)
case=mtdata
case[which(case>0)] <-1
case[which(case <=0)] <-0
#case=t(case)
a=case[1,]
b=names(a[which(a==1)])
MyList = list()
for( i in seq(1,dim(case)[1])){
	a=case[i,]
	b=names(a[which(a==1)])
	MyList[[i]]=b
	}
	
MyTrans<-as(MyList,"transactions")
MyRules<-apriori(data=MyTrans,parameter=list(support=0.6,confidence=0.99,target="rules"))
MyRules.sorted<-sort(x=MyRules,by="lift",decreasing=TRUE)
a<-is.significant(MyRules.sorted, MyTrans, method = "fisher", alpha = 0.01, adjust = "fdr")
rules=MyRules.sorted[which(a=="TRUE"),]
write(rules, file = paste(name,".s06c099.csv",sep=""), sep = ",")
subrules<-rules[seq(1,10)]
write(subrules, file = paste(name,".s06c099.top10.csv",sep=""), sep = ",")
pdf(paste(name,".s06c099.top10.pdf",sep=""))
for(i in seq(1,10)){
plot(subrules[i], method = "graph")
}
dev.off()

}


orifile<-'t2dmeta.stage1.format.fil.csv'
tdata<-read.csv(file=orifile,header=T)
tdata<-tdata[,-1]
cd1<-74
cd2<-71
a=1
b=cd1
cc=b+1
d=cd1+cd2
rownames(tdata)=seq(1,dim(tdata)[1])
data=t(tdata)
index=apply(data,1,function(x){length(x[which(x>0)])/length(x)})
data=data[which(index>0.1),]
case=t(data[,seq(a,b)])
apq(case,"t2d.control")
case=t(data[,seq(cc,d)])
apq(case,"t2d.case")


orifile<-'t2dmeta.stage2.format.fil.csv'
tdata<-read.csv(file=orifile,header=T)
tdata<-tdata[,-1]
cd1<-99
cd2<-99
a=1
b=cd1
cc=b+1
d=cd1+cd2
rownames(tdata)=seq(1,dim(tdata)[1])
data=t(tdata)
index=apply(data,1,function(x){length(x[which(x>0)])/length(x)})
data=data[which(index>0.1),]
case=t(data[,seq(a,b)])
apq(case,"t2d2.control")
case=data[,seq(cc,d)]
apq(case,"t2d2.case")









args=commandArgs(T)
orifile<<-args[1]
label1<<-args[2]
label2<<-args[3]
pct<<-as.numeric(args[4])
pjname<<-args[5]
type<<-args[6]
Test<<-args[7]

tdata<<-read.csv(file=orifile,header=T)
d1<-tdata[which(tdata[,1]==label1),]
d2<-tdata[which(tdata[,1]==label2),]
tdata<<-rbind(d1,d2)[,-1]
a<<-1
b<<-dim(d1)[1]
cc<<-b+1
d<<-dim(tdata)[1]
g<-apply(tdata,2,function(x){length(which(x>0))})
ind<-which(g/d>=pct)

if(length(ind)<2){
	info<-paste("No enough(<2) microbe passed the prevelence filter of ", pct)
	write.table(info,"Error.txt",row.names=F,col.names=F,quote=F)
	quit()
	}
	tdata<<-tdata[,ind]
	n<<-length(colnames(tdata))
	if(type=="MultiSearch"){
		listfile<-args[8]
		microbe<<-as.character(read.table(file=listfile)$V1)
		tdata<<-tdata[,colnames(tdata)%in% microbe]
		M<<-length(microbe)
	}
library(gtools,quietly=TRUE)
library(MASS,quietly = TRUE)
library(sfsmisc,quietly = TRUE)
library(pracma,quietly = TRUE)
library(DMwR,quietly = TRUE)
library(qcc,quietly = TRUE)
library(reshape2,quietly = TRUE)
library(gtools,quietly = TRUE)
library(stringr,quietly = TRUE)
library(parallel,quietly = TRUE)

###remove outliers
remove_outliers <- function(X){
	gAB<-quantile(X,probs=c(0.25,0.5,0.75))
	upper<-gAB[2]+1.5*(gAB[3]-gAB[1])	
	outl<-c(1:length(X))
	cirteia<-(X>upper)
	s<--1
	s<-sum(cirteia,na.rm=T)
	Y<-X
	if(s==0){
		Y<-X } else {
		Y<-X[-outl[cirteia]]
			}
		return(Y)
	}

###calculate non-overlapped area of two probability distributions
diffXY <- function(X,Y)
{
	densX<-density(X,from=0,bw=0.1)
	densY<-density(Y,from=0,bw=0.1)
	densXx<-densX$x
	densXy<-densX$y
	densYx<-densY$x
	densYy<-densY$y
	aera<-0
	if((max(densXx)<min(densYx))|(min(densXx)>max(densYx))){
		aera<-0}else{
		Xxi <- linspace(min(densXx), max(densXx), 100)
		Xppcub <- ppfit(densXx, densXy, Xxi, method = "cubic")
		Yxi <- linspace(min(densYx), max(densYx), 100)
		Yppcub <- ppfit(densYx, densYy, Yxi, method = "cubic")
		minx<-max(min(densXx),min(densYx))	
		maxx<-min(max(densXx),max(densYx))
		xs <- linspace(minx, maxx, 100)
		yX<-ppval(Xppcub, xs)
		yY<-ppval(Yppcub, xs)
		yscub <- apply(cbind(yX, yY), 1, min)
		aera<-integrate.xy(xs,yscub)
		aera<-max(0,aera)
	}
	diff=1-aera
	return(diff)
}

###functions-figures

###set ggplot theme
theme_zg <- function(..., bg='white'){
    require(grid)
    theme_classic(...) +
        theme(rect=element_rect(fill=bg),
              plot.margin=unit(rep(0.5,4), 'lines'),
              panel.background=element_rect(fill='transparent', color='black'),
              panel.border=element_rect(fill='transparent', color='transparent'),
              panel.grid=element_blank(),
              axis.title = element_text(color='black', vjust=0.1),
              axis.ticks.length = unit(0.4,"lines"),
              axis.ticks = element_line(color='black'),
              #axis.ticks.margin = unit(0.8,"lines"),
              legend.title=element_blank(),
              legend.key=element_rect(fill='transparent', color='transparent'))
}

###plot density distribution of Hotelling’s T2 statistics for a defined microbial community
Dens <- function(s,title)
{
	library(ggplot2,quietly = TRUE);library(RColorBrewer,quietly = TRUE)
	pm<-round(unlist(s[1]),2)
	pvalue<-round(unlist(s[2]),4)	
	tit= 	paste(title,paste("PM=",pm,sep=""),paste("pvalue=",pvalue,sep=""),sep=";")
    pj<-c(s[[3]],s[[4]])
    label<-c(rep(label1,length(s[[3]])),rep(label2,length(s[[4]])))
    df<-data.frame(pj,label)
	p=ggplot(df,aes(x=pj))+geom_density(aes(fill=label),alpha=0.9,adjust=0.5)+ scale_fill_manual(values=brewer.pal(3, "Set1")[c(1,2)])+theme_zg() + labs(title=tit) + theme(plot.title=element_text(hjust=0.5))
	print(p)
 }

### plot CA network for significant community alterations with FDR <0.05
Graph<-function(fa){
	library(igraph,quietly = TRUE)
	library(RColorBrewer,quietly = TRUE)
	ndscale = function(x,a,b,c,d){c + (x-a)/(b-a)*(d-c)}
	pmvalue=read.table(fa,head=T,sep="\t")
	pmvalue=pmvalue[which(pmvalue$FDR<0.05),]
	if(dim(pmvalue)[1]>0){
		g <- graph_from_data_frame(pmvalue[,c(1,2)], directed=FALSE)
		d<-degree(g, v=V(g), mode = c("all", "out", "in", "total"),loops = TRUE, normalized = TRUE)
		if(max(d)>min(d)){d<-ndscale(d, min(d), max(d), 1, 20)}
		else{d<-d/min(d)*10}
		weight<-pmvalue$PMscore
		weight<-ndscale(weight,min(weight),max(weight),0.1,3)
		V(g)$color<-brewer.pal(3,"Set1")[1]
		V(g)$size<-d
		V(g)$frame.color<-NA
		E(g)$width<-weight
		E(g)$color<-brewer.pal(3,"Set1")[3]
		dsize<-ndscale(d, min(d), max(d), 0.05, 1)
		V(g)$label.cex<-dsize
		V(g)$label.dist=1
		plot(g,layout = layout.sphere,main="CA network for 2D PM")
	}else{
		plot.new()
		text(0.5,1,"No significantly changed assoction ( FDR<0.05 )")
	}
}

### plot interactions within a defined module
Module<-function(df,n){
	library(igraph,quietly = TRUE)
	library(RColorBrewer,quietly = TRUE)
	ndscale = function(x,a,b,c,d){c + (x-a)/(b-a)*(d-c)}
	pmvalue=df
	pmvalue=pmvalue[which(pmvalue$FDR<0.05),]
	node=unique(c(as.character(pmvalue[,1]),as.character(pmvalue[,2])))
	diff=setdiff(n,node)
	if(dim(pmvalue)[1]>0){
		g <- graph_from_data_frame(pmvalue[,c(1,2)], directed=FALSE)
		g<-g+vertices(diff)
		weight<-pmvalue$PMscore
		V(g)$color<-brewer.pal(6,"Paired")[6]
		V(g)$frame.color<-NA
		E(g)$width<-weight*3
		E(g)$color<-brewer.pal(3,"Set1")[3]
		V(g)$label.cex<-0.8
		V(g)$label.dist=1
		plot(g,mark.groups=(list(1:length(V(g)))),mark.col =brewer.pal(7,"Blues")[1], alpha = 0.3,mark.border = NA,main="CA network for Sub-community")
	}else{
		g<-make_empty_graph() + vertices(diff)
		V(g)$color<-brewer.pal(6,"Paired")[6]
		V(g)$frame.color<-NA
		V(g)$label.cex<-0.8
		V(g)$label.dist=1
		plot(g,mark.groups=(list(1:length(V(g)))),mark.col =brewer.pal(7,"Blues")[1], alpha = 0.3,mark.border = NA,main="CA network for Sub-community")
	}
}

### plot heatmap for pairwise PM2CA  score in the 2DScaning
heat<-function(fb){
	library(pheatmap,quietly = TRUE)
	library(reshape2,quietly = TRUE)
	df=read.table(fb,head=T,sep="\t")
	df<-df[,c(1,2,3)]
	revdf=data.frame(df$OTU2,df$OTU1,df$PMscore)
	colnames(revdf)=c("OTU1","OTU2","PMscore")
	df=rbind(df,revdf)
	df=df[!duplicated(df),]
	mat=dcast(df, OTU1~ OTU2,value.var="PMscore")
	mat[is.na(mat)]<-0
	rownames(mat)<-mat[,1]
	mat=mat[,-1]
	pheatmap(mat, gaps_col=0.2, gaps_row=0.2, fontsize_row = 5,fontsize_col=5,main="Heatmap for pairwise PM score",fontsize=10)
}

###functions-figures

###functions-main function
### calculate 2D PM score
MSD2 <- function(W){
	A<-tdata[a:b,W]
	B<-tdata[cc:d,W]
	A<-A[which(apply(A,1,sum)>0),]
	B<-B[which(apply(B,1,sum)>0),]
	if((dim(A)[1]>2)&(dim(B)[1]>2)){
		SA<-stats.T2.single(A)
		SAcov<-SA$cov
		SAcenter<-SA$center
		Astatistics<-SA$statistics
		SB<-stats.T2.single(B)
		SBcov<-SB$cov
		SBcenter<-SB$center
		Bstatistics<-SB$statistics		
		AB<-sweep(A,2,SBcenter, "-" )
		AB<-as.matrix(AB)
		ABstatistics<-diag(AB%*%solve(SBcov)%*%(t(AB)))
		BA<-sweep(B,2,SAcenter, "-" )
		BA<-as.matrix(BA)
		BAstatistics<-diag(BA%*%solve(SAcov)%*%(t(BA)))
		ABstatistics<-remove_outliers(ABstatistics)
		BAstatistics<-remove_outliers(BAstatistics)
		Astatistics<-remove_outliers(Astatistics)
		Bstatistics<-remove_outliers(Bstatistics)
		if((length(Astatistics)>1)&(length(Bstatistics)>1)&(length(BAstatistics)>1)&(length(ABstatistics)>1)){
			diffs<-max(diffXY(ABstatistics,Bstatistics),diffXY(BAstatistics,Astatistics))
			if(Test =="wc"){
				pval<-min(wilcox.test(ABstatistics,Bstatistics)$p.value,wilcox.test(BAstatistics,Astatistics)$p.value)
			}else{
				pval<-min(ks.test(ABstatistics,Bstatistics)$p.value,ks.test(BAstatistics,Astatistics)$p.value)
			}
		}
		else{diffs<-NA;pval<-NA}
	}else{diffs<-NA;pval<-NA}
	cm<-paste(str_c(colnames(tdata)[W],collapse='\t'),diffs,pval,sep="\t")
	return (cm)
}

### calculate the 1D PM score 
MSD1 <- function(Z){
	A<-tdata[a:b,Z]
	B<-tdata[cc:d,Z]
	SA<-stats.T2.single(A)
	SAcov<-SA$cov
	SAcenter<-SA$center
	Astatistics<-SA$statistics
	SB<-stats.T2.single(B)
	SBcov<-SB$cov
	SBcenter<-SB$center
	Bstatistics<-SB$statistics	
    	AB<-data.frame(A-SBcenter)		
	AB<-as.matrix(AB)
	ABstatistics<-diag(AB%*%solve(SBcov)%*%(t(AB)))
	BA<-data.frame(B-SAcenter)
	BA<-as.matrix(BA)
	BAstatistics<-diag(BA%*%solve(SAcov)%*%(t(BA)))
	ABstatistics<-remove_outliers(ABstatistics)
	BAstatistics<-remove_outliers(BAstatistics)
	Astatistics<-remove_outliers(Astatistics)
	Bstatistics<-remove_outliers(Bstatistics)
	if((length(Astatistics)>1)&(length(Bstatistics)>1)&(length(BAstatistics)>1)&(length(ABstatistics)>1)){
		diffs<-max(diffXY(ABstatistics,Bstatistics),diffXY(BAstatistics,Astatistics))
		if(Test =="wc"){
			pval<-min(wilcox.test(ABstatistics,Bstatistics)$p.value,wilcox.test(BAstatistics,Astatistics)$p.value)
		}else{
			pval<-min(ks.test(ABstatistics,Bstatistics)$p.value,ks.test(BAstatistics,Astatistics)$p.value)
		}	
	}
	else{diffs<-NA;pval<-NA}
	cm<-paste(colnames(tdata)[Z],colnames(tdata)[Z],diffs,pval,sep="\t")
	return (cm)
}

### calculate the PM score for a defined module
MSDs <- function(S){
	A<-tdata[a:b,S]
	B<-tdata[cc:d,S]
	A<-A[which(apply(A,1,sum)>0),]
	B<-B[which(apply(B,1,sum)>0),]
	SA<-stats.T2.single(A)
	SAcov<-SA$cov
	SAcenter<-SA$center
	Astatistics<-SA$statistics
	SB<-stats.T2.single(B)
	SBcov<-SB$cov
	SBcenter<-SB$center
	Bstatistics<-SB$statistics	
    AB<-data.frame(A-SBcenter)		
	AB<-as.matrix(AB)
	ABstatistics<-diag(AB%*%solve(SBcov)%*%(t(AB)))
	BA<-data.frame(B-SAcenter)
	BA<-as.matrix(BA)
	BAstatistics<-diag(BA%*%solve(SAcov)%*%(t(BA)))
	ABstatistics<-remove_outliers(ABstatistics)
	BAstatistics<-remove_outliers(BAstatistics)
	Astatistics<-remove_outliers(Astatistics)
	Bstatistics<-remove_outliers(Bstatistics)
	if((length(Astatistics)>1)&(length(Bstatistics)>1)&(length(BAstatistics)>1)&(length(ABstatistics)>1)){
		Bdiff<-diffXY(ABstatistics,Bstatistics)
		Adiff<-diffXY(BAstatistics,Astatistics)
		diffs<-max(Adiff,Bdiff)
		if(Test =="wc"){
			pval<-min(wilcox.test(ABstatistics,Bstatistics)$p.value,wilcox.test(BAstatistics,Astatistics)$p.value)
		}else{
			pval<-min(ks.test(ABstatistics,Bstatistics)$p.value,ks.test(BAstatistics,Astatistics)$p.value)
		}
	}
	else{Adiff<-NA;diffs<-NA;pval<-NA}	
	pm<-paste(str_c(colnames(tdata)[S],collapse='\t'),diffs,pval,sep="\t")
	if(diffs==Adiff){
		Aprojection <-Astatistics
		Bprojection<-BAstatistics		
	}else{
		Aprojection <-ABstatistics
		Bprojection<-Bstatistics
	}
	res<-list(pm=pm,Apj=Aprojection,Bpj=Bprojection)
	return (res)
}


###functions-main function
###main
if(type=="2DScan"){
	outtable<-paste(pjname,".2Dscan.res.txt",sep="")
	outpdf<-paste(pjname,".2Dscan.res.pdf",sep="")	
	cb1<-combinations(n,1)
	cl <- makeCluster(detectCores(),type="FORK")
	results1 <- parApply(cl,cb1,1,MSD1)
	res1<-unlist(results1)
	stopCluster(cl)
	cb2<-combinations(n,2)
	cl <- makeCluster(detectCores(),type="FORK")
	results <- parApply(cl,cb2,1,MSD2)
	res2<-unlist(results)
	stopCluster(cl)
	res.df<-c(res1,res2)
	write.table(res.df,outtable,quote=F,sep="\t",row.names=F,col.names=F)
	temp<-read.table(outtable,sep="\t")
	temp$fdr<-p.adjust(temp$V4,method="fdr",n=length(temp$V4))
	colnames(temp)=c("OTU1","OTU2","PMscore","Raw pvalue","FDR")
	write.table(temp, outtable,quote=F,sep="\t",row.names=F)
	pdf(outpdf)
	heat(outtable)
	Graph(outtable)
	dev.off()
}else{	
	outtable<-paste(pjname,".MultiSearch.pairwise.res.txt",sep="")
	outpdf<-paste(pjname,".MultiSearch.res.pdf",sep="")		
	Ms<-MSDs(seq(1,M))
	cmb<-combinations(M,2)
	resmD <- apply(cmb,1,MSD2)
	write.table(resmD,outtable,quote=F,row.names=F,col.names=F,sep="\t")
	temp<-read.table(outtable,sep="\t")
	temp$fdr<-p.adjust(temp$V4,method="fdr",n=length(temp$V4))
	colnames(temp)=c("OTU1","OTU2","PMscore","Raw.pvalue","FDR")
	write.table(temp, outtable,quote=F,sep="\t",row.names=F)
	pdf(outpdf)
	temp<-unlist(strsplit(Ms$pm,split="\t"))
	Ms.pm<-as.numeric(temp[length(temp)-1])
	Ms.pvalue<-as.numeric(temp[length(temp)])
	Apj<-Ms$Apj;Bpj<-Ms$Bpj
	Dens(list(Ms.pm,Ms.pvalue,Apj,Bpj),"PM distribution")
	temp=read.table(outtable,head=T)
	Module(temp,microbe)
	dev.off()
}
	
###main
	
	

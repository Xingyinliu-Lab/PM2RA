import argparse
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--Input', type=str, default = None)
parser.add_argument('--Case', type=int, default = None)
parser.add_argument('--Control', type=int, default=None)
parser.add_argument('--TopN', type=int, default=10)
parser.add_argument('--Thread', type=int, default=1)
#parser.add_argument('--Output', type=str)

args = parser.parse_args()

import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
numpy2ri.activate()

base = importr('base')
utils=importr('utils')
### datafile place
datafile = args.Input

pddata=pd.read_csv(filepath_or_buffer=datafile,header=0)
bioName=list(pddata.columns.values)
bioName=np.delete(bioName,0,0)
rscript_readdata = '''
tdata <- read.csv('%s',header=T)
n=length(colnames(tdata))-1
#n=20
'''%datafile
robjects.r(rscript_readdata)
#print(robjects.r['summary'](robjects.r['tdata']))
# bio kind number
bioN=int(robjects.r['n'][0])
# top N
N=args.TopN
### gloable var
a=1
b=args.Case
# c=b+1
c=b+1
d=args.Case+args.Control

robjects.globalenv['a']=a
robjects.globalenv['b']=b
robjects.globalenv['c']=c
robjects.globalenv['d']=d
robjects.globalenv['N']=N


rscript_lib = '''
library(parallel)
library(sfsmisc)
library(pracma)
library(gtools)
library(qcc)
'''
robjects.r(rscript_lib)

rscript_fun = '''
func<-function(x){
	cb=combinations(n,x)
	# res <- apply(cb,1,MSD)
	cl <- makeCluster(8,type="FORK")
	results <- parApply(cl,cb,1,MSD)
	res.df<-unlist(results,recursive =FALSE)
	stopCluster(cl)
	return(res.df)
}

func_bycb<-function(cb){
	cl <- makeCluster(8,type="FORK")
	results <- parApply(cl,cb,1,MSDforIteration)
	res.df<-unlist(results,recursive =FALSE)
	stopCluster(cl)
	return(res.df)
}

removeoutlier3 <- function(X,upper)
{
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

diffXY <- function(X,Y)
{
	densX<-density(X,from=0,bw=0.1)
	densY<-density(Y,from=0,bw=0.1)
	densXx<-densX$x
	densXy<-densX$y
	densYx<-densY$x
	densYy<-densY$y
	aera<-0
	if((max(densXx)<min(densYx))|(min(densXx)>max(densYx)))
	{
		aera<-0
	}
	else{
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
	}
	diff=1-aera
	return(diff)
}
changetovec<-function(x){
	res.df<-data.frame(id=c(NA))
	for(i in 1:length(x))
		res.df<-data.frame(res.df,x[i])
	res.df<-subset(res.df,select=-id)
	colnames(res.df)<-c()
	td<-t(res.df)
	return(td)
}
MSD <- function(cmb){
	A<-tdata[a:b,cmb+1]
	B<-tdata[c:d,cmb+1]
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
	AB<-sweep(A,2,SBcenter, "-" )
	AB<-as.matrix(AB)
	ABstatistics<-diag(AB%*%solve(SBcov)%*%(t(AB)))
	BA<-sweep(B,2,SAcenter, "-" )
	BA<-as.matrix(BA)
	BAstatistics<-diag(BA%*%solve(SAcov)%*%(t(BA)))
		
	gAB=unname(quantile(ABstatistics,probs=c(0.25,0.5,0.75)))
	ABstatistics<-removeoutlier3(ABstatistics,(gAB[2]+1.5*(gAB[3]-gAB[1])))
	gBA=unname(quantile(BAstatistics,probs=c(0.25,0.5,0.75)))
	BAstatistics<-removeoutlier3(BAstatistics,(gBA[2]+1.5*(gBA[3]-gBA[1])))
	gA=unname(quantile(Astatistics,probs=c(0.25,0.5,0.75)))
	Astatistics<-removeoutlier3(Astatistics,(gA[2]+1.5*(gA[3]-gA[1])))
	gB=unname(quantile(Bstatistics,probs=c(0.25,0.5,0.75)))
	Bstatistics<-removeoutlier3(Bstatistics,(gB[2]+1.5*(gB[3]-gB[1])))

	pval<-min(ks.test(ABstatistics,Bstatistics)$p.value,ks.test(BAstatistics,Astatistics)$p.value)
	if(pval>0.05){
		diffs<-0
	}
		else{
		diffs<-max(diffXY(ABstatistics,Bstatistics),diffXY(BAstatistics,Astatistics))
		}
	#cm<-paste(str_c(colnames(tdata)[cmb+1],collapse=';'),diffs,pval,colapse="\t")
	v<-rep(0,time=n)
	v[cmb]=1
	cm<-c(v,diffs,pval)
	return (cm)
}

MSDforIteration <- function(rawcmb){
    sq<-seq(1,n)
    cmb<-sq[rawcmb>0]
	A<-tdata[a:b,cmb+1]
	B<-tdata[c:d,cmb+1]
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
	AB<-sweep(A,2,SBcenter, "-" )
	AB<-as.matrix(AB)
	ABstatistics<-diag(AB%*%solve(SBcov)%*%(t(AB)))
	BA<-sweep(B,2,SAcenter, "-" )
	BA<-as.matrix(BA)
	BAstatistics<-diag(BA%*%solve(SAcov)%*%(t(BA)))
		
	gAB=unname(quantile(ABstatistics,probs=c(0.25,0.5,0.75)))
	ABstatistics<-removeoutlier3(ABstatistics,(gAB[2]+1.5*(gAB[3]-gAB[1])))
	gBA=unname(quantile(BAstatistics,probs=c(0.25,0.5,0.75)))
	BAstatistics<-removeoutlier3(BAstatistics,(gBA[2]+1.5*(gBA[3]-gBA[1])))
	gA=unname(quantile(Astatistics,probs=c(0.25,0.5,0.75)))
	Astatistics<-removeoutlier3(Astatistics,(gA[2]+1.5*(gA[3]-gA[1])))
	gB=unname(quantile(Bstatistics,probs=c(0.25,0.5,0.75)))
	Bstatistics<-removeoutlier3(Bstatistics,(gB[2]+1.5*(gB[3]-gB[1])))

	pval<-min(ks.test(ABstatistics,Bstatistics)$p.value,ks.test(BAstatistics,Astatistics)$p.value)
	if(pval>0.05){
		diffs<-0
	}
		else{
		diffs<-max(diffXY(ABstatistics,Bstatistics),diffXY(BAstatistics,Astatistics))
		}
	#cm<-paste(str_c(colnames(tdata)[cmb+1],collapse=';'),diffs,pval,colapse="\t")
	v<-rep(0,time=n)
	v[cmb]=1
	cm<-c(v,diffs,pval)
	return (cm)
}
'''
robjects.r(rscript_fun)
#### cal all pm

# rscript_calALL = '''
# x=seq(2,n)
# starttime<-Sys.time()
# res<-sapply(x,func,simplify = "array")
# res<-changetovec(res)
# endtime<-Sys.time()
# ctime<-endtime-starttime
# '''
# robjects.r(rscript_calALL)
# print(robjects.r['summary']('res'))
# print(robjects.r['ctime'])
# npres=np.array(robjects.r['res'])
# print(npres.shape)
# print(npres)

# ### cal pm on a certain dimension
# dimx=str(2)
#
# rscript_calDimX = '''
# x=%s
# starttimeX<-Sys.time()
# resX<-sapply(x,func,simplify = "array")
# endtimeX<-Sys.time()
# ctimeX<-endtimeX-starttimeX
# '''%dimx
# robjects.r(rscript_calDimX)
#
#
#
# print(robjects.r['summary']('resX'))
# print(robjects.r['ctimeX'])
# npres=np.array(robjects.r['resX'])
# npres=np.reshape(npres,newshape=(npres.shape[0],npres.shape[1]))
# npres=np.transpose(npres)
# print(npres)


#### find module
## first step traversing 2nd dimension

dimx=str(2)

rscript_calDimX = '''
x=%s
starttimeX<-Sys.time()
resX<-sapply(x,func,simplify = "array")
endtimeX<-Sys.time()
ctimeX<-endtimeX-starttimeX
'''%dimx
robjects.r(rscript_calDimX)

#print(robjects.r['summary']('resX'))
#print(robjects.r['ctimeX'])
npres=np.array(robjects.r['resX'])
npres=np.reshape(npres,newshape=(npres.shape[0],npres.shape[1]))
npres=np.transpose(npres)
tmpres=-1*npres
npres=npres[tmpres[:,npres.shape[1]-2].argsort()]
del(tmpres)
#print(npres)

def getTopN(dataset,N,bioN):
    j=1
    res=dataset[0,:]
    for i in range(N-1):
        notout = True
        while(notout):
            if(j+1>dataset.shape[0]):
                break
            if(i==0):
                tmp = np.concatenate(([res],[dataset[j,:]]),axis=0)
                tmp = np.delete(tmp, [bioN, bioN + 1], axis=1)
                if(sum(sum(tmp)>1)==0):
                    notout=False
                    break
                j=j+1
            if (i > 0):
                tmp = np.concatenate((res, [dataset[j, :]]), axis=0)
                tmp = np.delete(tmp, [bioN, bioN + 1], axis=1)
                if (sum(sum(tmp) > 1) == 0):
                    notout = False
                    break
                j = j + 1
        if(i==0):
            res = np.concatenate(([res],[dataset[j,:]]),axis=0)
            j = j + 1
        if (i > 0):
            res = np.concatenate((res, [dataset[j, :]]), axis=0)
            j = j + 1
    return res

topres=getTopN(npres,N,bioN)
print('Start Iter: 1  ----------')
print(topres)



### step2 first iteration
## prepare pm combination

def prepareCombination( top,bioN):
    #N total number of micorbio
    top=np.delete(top,[bioN,bioN+1],axis=1)
    x=np.arange(1,bioN+1)
    comb=None
    for i in range(top.shape[0]-1):
        for j in range(i+1,top.shape[0]):
            combination = np.zeros(bioN)
            tmp=top[j, :]+top[i, :]
            tmpselected = x[tmp > 0]
            combination[tmpselected-1] = 1
            if comb is None:
                comb = [combination]
            else:
                comb= np.concatenate(([combination],comb))
    for i in range(top.shape[0]):
        selected = x[top[i, :] > 0]
        remained = x[-top[i, :] > -1]
        for j in range(len(remained)):
            combination = np.zeros(bioN)
            tmpselected=np.append(selected,remained[j])
            combination[tmpselected-1] = 1
            comb= np.concatenate(([combination],comb))
    return comb

## implentment iteration
thresholdN=10
iterationT=2
while(iterationT<thresholdN):
    iterationT=iterationT+1
    comb = prepareCombination(topres, bioN)
    rcomb=numpy2ri.py2ri(comb)
    rcomb=robjects.Matrix(rcomb)
    robjects.globalenv['rcomb']=rcomb
    rscript_calC = '''
    rcomb <- data.frame(rcomb)
    starttimeC<-Sys.time()
    resC<-func_bycb(rcomb)
    endtimeC<-Sys.time()
    ctimeC<-endtimeC-starttimeC
    '''
    robjects.r(rscript_calC)
    #print(robjects.r['head']('rcomb'))
    npresC=np.array(robjects.r['resC'])
    npresC=np.reshape(npresC,newshape=(npresC.shape[0],npresC.shape[1]))
    npresC=np.transpose(npresC)
    npres=np.concatenate((npres,npresC),axis=0)
    tmpres=-1*npres
    npres=npres[tmpres[:,npres.shape[1]-2].argsort()]
    del(tmpres)
    del(npresC)
    resTopres=getTopN(npres,N,bioN)
    if((topres==resTopres).all()):
        break
    else:
        topres=resTopres


topindex=np.delete(topres,[bioN,bioN+1],axis=1)
for i in range(topres.shape[0]):
    print('------------------------------------------')
    print('Module %i: \t Diff:%f \t Pvalue:%f' %(i,topres[i,bioN],topres[i, bioN+1]))
    print('BioName:')
    print(bioName[topindex[i,:]>0])

    

#print(comb)

## cal combination



numpy2ri.deactivate()

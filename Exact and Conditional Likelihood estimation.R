### Analysis of monhtly log returns of IBM and KO stocks
### Sample period: 2001-2011 for 132 observations
### Purpose: to demonstrate the difference between exact and conditional 
### likelihood estimation in small sample.
###

da=read.table("m-ibmko-0111.txt",header=T)
head(da)

lrtn=log(da[,2:3]+1)*100
dim(da)

tdx=c(1:132)/12+2001
colnames(lrtn) <- c("ibm","ko")
MTSplot(lrtn,tdx)
ccm(lrtn)

mq(zt,20)

yt=diffM(lrtn)
mm=ccm(yt)


m1=VMA(lrtn,q=1,include.mean=F)


m2=VMAe(lrtn,q=1,include.mean=F)

m1=VMA(yt,q=1,include.mean=F)

m2=VMAe(yt,q=1,include.mean=F)



t1=m1$Theta; t2=m2$Theta
eigen(t1)

eigen(t2)

#### Compute the theoretical covariance and CCM matrices 
####
phi=matrix(c(.816,-1.116,-.623,1.074,-.643,.615,.592,-.133),2,4)
phi

theta=matrix(c(0,-.801,-1.248,0),2,2)
sig=matrix(c(4,2,2,5),2,2)
VARMAcov(Phi=phi,Theta=theta,Sigma=sig,lag=2)


### Simulation
m1=VARMAsim(400,arlags=c(1,2),malags=c(1),phi=phi,theta=theta,sigma=sig)

zt=m1$series
m2=Eccm(zt,maxp=5,maxq=6)

### U.S. Hog data
da=read.table("ushog.txt",header=T)
head(da)

m1=Eccm(da,maxp=5,maxq=6)

VARorder(da,maxp=9)  ## VAR order
### Analysis of the growth rates, in percentages, of PCE and DSPI
da1=read.table("m-pce.txt",header=T)
da2=read.table("m-dspi.txt",header=T)
head(da1)
head(da2)
x=cbind(da1$pce,da2$dspi)
x=log(x)
zt=diffM(x)*100
colnames(zt) <- c("pceg","dspig")
tdx=da1[,1]+da1[,2]/12
dim(da1)


MTSplot(zt,tdx[2:639]) ### 1 less data point due to differencing

VARorder(zt)

m1=VAR(zt,3)

MTSdiag(m1a)  ## Model checking

### VARMA modeling
Eccm(zt,maxp=6,maxq=6)

m2=VARMA(zt,p=3,q=1)
m2a=refVARMA(m2,thres=0.8)  ## refine the model

m2b=refVARMA(m2a,thres=1) ## Further refinement

MTSdiag(m2b) ### Model checking
names(m2b)

phi=m2b$Phi; theta=m2b$Theta; sig=m2b$Sigma
VARMAirf(Phi=phi,Theta=theta,Sigma=sig,orth=F)
da=read.table("m-hsmort7112.txt",header=T)
head(da)

zt=da[,3:4]
colnames(zt) <- c("hs","mort")
dzt=diffM(zt)
dim(da)
tdx=da[,1]+da[,2]/12
dzt[,1]=dzt[,1]/1000   ### Make the scale of the two series similar
MTSplot(dzt,tdx[2:492])


VARorder(dzt)

m1=VAR(dzt,4)
MTSdiag(m1a)

Eccm(dzt,maxp=6,maxq=6)

### VARMA(2,1) model
m2=VARMA(dzt,p=2,q=1)

m2a=refVARMA(m2,thres=0.8)

m2b=refVARMA(m2a,thres=1)
m2c=refVARMA(m2b,thres=1)

MTSdiag(m2c)

##### VARMA(1,2) model
m3=VARMA(dzt,p=1,q=2)

m3a=refVARMA(m3,thres=0.6)

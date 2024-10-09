#######################
# CAP mediation
#######################

library("mvtnorm")

library("mediation")

rm(list=ls())

######################################
# parameter setting
p<-10
q<-0
n<-100
nT.mean<-100
nT<-rep(nT.mean,n)

# M
set.seed(100)
gamma.mat0<-matrix(runif(p),nrow=p,ncol=p)
gamma.mat<-qr.Q(qr(gamma.mat0))
for(j in 1:p)
{
  if(gamma.mat[which.max(abs(gamma.mat[,j])),j]<0)
  {
    gamma.mat[,j]<-(-gamma.mat[,j])
  }
}
Gamma<-gamma.mat

# eigenvalue \neq 0 index
eigenv.nz.idx<-c(2,4)

# W
if(q>0)
{
  W.name<-paste0("W",1:q)
  
  phi1<-cbind(c(0.5,-0.5),c(-0.5,0.5))
  colnames(phi1)<-paste0("Dim",eigenv.nz.idx)
  rownames(phi1)<-W.name
  
  phi2<-c(-0.5,0.5)
}else
{
  phi1<-NULL
  phi2<-NULL
  W.name<-NULL
}

# alpha
alpha0<-c(1,-1)
alpha<-rbind(cbind(c(1),c(1)),phi1)
colnames(alpha)<-paste0("Dim",eigenv.nz.idx)
rownames(alpha)<-c("X",W.name)

# beta
beta<-c(1,1)

gamma0<-1
gamma<-c(1,phi2)

eigenv.base<-seq(3,-1,length.out=5)
eigenv.base.vec<-c(eigenv.base,seq(-1.5,-3,length.out=p-length(eigenv.base)))

# eigenvalue=0 corresponding component SD
eigenv.nr.sd<-0.1
tau<-0.1
sigma<-0.1
######################################

######################################
# generate X and W: (X,W) is a n by q data matrix
set.seed(100)
X0<-rbinom(n,size=1,prob=0.5)
if(q>0)
{
  W<-cbind(rnorm(n,mean=0,sd=0.5),rbinom(n,size=1,prob=0.5))
}else
{
  W<-NULL
}
X<-cbind(X0,W)
colnames(X)<-c("X",W.name)

# delta: a n×p matrix, the true projection matrix used to generate the data.
delta<-matrix(NA,n,p)

# Sigma:a p×p×n array, the covariance matrix of the n subjects.
Sigma<-array(NA,c(p,p,n))
for(j in 1:p)
{
  if(sum(j==eigenv.nz.idx)!=0)
  {
    itmp<-which(eigenv.nz.idx==j)
    eta<-rnorm(n,mean=0,sd=tau)
    # for some eigenvalue, it satisfy the following log-linear model 
    delta[,j]<-exp(alpha0[itmp]+X%*%alpha[,itmp]+eta)
  }else
  {
    delta[,j]<-exp(eigenv.base.vec[j])
  }
}
# generate M: a list of length n by n, where each list element is a T by p matrix, the data matrix of T observations from p features.
M<-vector("list",length=n)
for(i in 1:n)
{
  Sigma[,,i]<-Gamma%*%diag(delta[i,])%*%t(Gamma)
  M[[i]]<-rmvnorm(n=nT[i],mean=rep(0,p),sigma=Sigma[,,i])
}

# generate Y:a n×1 outcome vector
set.seed(100)
epsilon<-rnorm(n,mean=0,sd=sigma)
Y<-gamma0+X%*%gamma+log(delta[,eigenv.nz.idx])%*%beta+ epsilon
######################################

######################################
# method parameters
max.itr<-1000
tol<-1e-4
score.return<-TRUE
trace<-FALSE

stop.crt<-c("DfD")
nD<-3
DfD.thred<-2

# verbose<-TRUE
verbose<-FALSE

boot<-TRUE
sims<-500
boot.ci.type<-c("se")
conf.level<-0.95
######################################

######################################
# run function
source("CAPMediation.R")

# sample covariance of M
M.cov<-array(NA,c(p,p,n))
for(i in 1:n)
{
  M.cov[,,i]<-cov(M[[i]])
}
H<-apply(M.cov,c(1,2),mean,na.rm=TRUE)


# call CAPMediation function 
re<-CAPMediation(X,M,Y,H=H,stop.crt=stop.crt,DfD.thred=DfD.thred,Y.remove=FALSE)

# get bootstrap result
re.boot<-vector("list",length=ncol(re$theta))
names(re.boot)<-paste0("C",1:ncol(re$theta))
for(jj in 1:ncol(re$theta))
{
  re.boot[[jj]]<-CAPMediation_boot(X,M,Y,theta=re$theta[,jj],H=H,boot=boot,sims=sims,boot.ci.type=boot.ci.type,conf.level=conf.level,verbose=verbose)
}
######################################

# save output to R data file 
save.image("example.RData")


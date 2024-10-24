library(mvtnorm)

## read data, you need to make sure the data file is in your current working directory 
earthquakes.dat <- read.delim("earthquakes.txt")
earthquakes.dat$Quakes=as.numeric(earthquakes.dat$Quakes)

y.dat=earthquakes.dat$Quakes[1:100] ## this is the training data
y.new=earthquakes.dat$Quakes[101:103] ## this is the test data

# Plot Data

par(mfrow = c(1,1))
plot(earthquakes.dat$Quakes~earthquakes.dat$Year, xlab = "Year", ylab = "Earthquakes", main="Earthquakes by Year")
cor(earthquakes.dat$Quakes,earthquakes.dat$Year)

# ACF and PCF plot

par(mfrow=c(1,2))

acf(y.dat,main="ACF",xlab='Lag')

pacf(y.dat,main="PACF",xlab='Lag')

# Check best Hyperparameters

## prior sensitivity analysis
## plot posterior distribution of phi_1 and phi_2 on a grid 

library(colorRamps)
library(leaflet)
library(fields)
library(mvtnorm)

## generate grid
coordinates_1=seq(-3,3,length.out = 100)
coordinates_2=seq(-3,3,length.out = 100)
coordinates=expand.grid(coordinates_1,coordinates_2)
coordinates=as.matrix(coordinates)

## set up
N=100
p=2  ## order of AR process
n.all=length(y.dat) ## T, total number of data

Y=matrix(y.dat[(p+1):n.all],ncol=1)
Fmtx=matrix(c(y.dat[2:(n.all-1)],y.dat[1:(n.all-2)]),nrow=p,byrow=TRUE)
n=length(Y)

## function to compute parameters for the posterior distribution of phi_1 and phi_2
## the posterior distribution of phi_1 and phi_2 is a multivariate t distribution

cal_parameters=function(m0=matrix(c(0,0),nrow=2),C0=diag(p),n0,d0){
  e=Y-t(Fmtx)%*%m0
  Q=t(Fmtx)%*%C0%*%Fmtx+diag(n)
  Q.inv=chol2inv(chol(Q))  ## similar as solve, but more robust
  A=C0%*%Fmtx%*%Q.inv
  m=m0+A%*%e
  C=C0-A%*%Q%*%t(A)
  n.star=n+n0
  d.star=t(Y-t(Fmtx)%*%m0)%*%Q.inv%*%(Y-t(Fmtx)%*%m0)+d0
  
  params=list()
  params[[1]]=n.star
  params[[2]]=d.star
  params[[3]]=m
  params[[4]]=C
  
  return(params)
}

## evaluate density at the grid points
get_density=function(param){
  location=param[[3]]
  scale=as.numeric(param[[2]]/param[[1]])*param[[4]]
  density=rep(0,N^2)
  
  for (i in 1:N^2) {
    xi=coordinates[i,]
    density[i]=dmvt(xi,delta=location,sigma=scale,df=param[[1]])
  }
  
  density_expand=matrix(density,nrow=N)
  return(density_expand)
}

## calculate density for three sets of hyperparameters
params1=cal_parameters(n0=2,d0=2)
params2=cal_parameters(n0=6,d0=1)
params3=cal_parameters(m0=matrix(c(-0.5,-0.5),nrow=2),n0=6,d0=1)

col.list=matlab.like2(N)
Z=list(get_density(params1),get_density(params2),get_density(params3))

op <- par(mfrow = c(1,3),
          oma = c(5,4,0,0) + 0.1,
          mar = c(4,4,0,0) + 0.2)
image(coordinates_1,coordinates_2,Z[[1]],col=col.list,
      zlim=range(unlist(Z)),xlab=expression(phi[1]),ylab=expression(phi[2]))
image(coordinates_1,coordinates_2,Z[[2]],col=col.list,
      zlim=range(unlist(Z)),xlab=expression(phi[1]),ylab=expression(phi[2]))
image(coordinates_1,coordinates_2,Z[[3]],col=col.list,
      zlim=range(unlist(Z)),xlab=expression(phi[1]),ylab=expression(phi[2]))

# Select Model Order

n.all=length(y.dat)
p.star=5
Y=matrix(y.dat[(p.star+1):n.all],ncol=1)
sample.all=matrix(y.dat,ncol=1)
n=length(Y)
p=seq(1,p.star,by=1)

design.mtx=function(p_cur){
  Fmtx=matrix(0,ncol=n,nrow=p_cur)
  for (i in 1:p_cur) {
    start.y=p.star+1-i
    end.y=start.y+n-1
    Fmtx[i,]=sample.all[start.y:end.y,1]
  }
  return(Fmtx)
}

criteria.ar=function(p_cur){
  Fmtx=design.mtx(p_cur)
  beta.hat=chol2inv(chol(Fmtx%*%t(Fmtx)))%*%Fmtx%*%Y
  R=t(Y-t(Fmtx)%*%beta.hat)%*%(Y-t(Fmtx)%*%beta.hat)
  sp.square=R/(n-p_cur)
  aic=2*p_cur+n*log(sp.square)
  bic=log(n)*p_cur+n*log(sp.square)
  result=c(aic,bic)
  return(result)
}

criteria=sapply(p,criteria.ar)

par(mfrow = c(1,1))
plot(p,criteria[1,],type='p',pch='a',col='red',xlab='AR order p',ylab='Criterion',main='',
     ylim=c(min(criteria)-10,max(criteria)+10))
points(p,criteria[2,],pch='b',col='blue')

## set up
N=100
p=3 ## order of AR process
n.all=length(y.dat) ## T, total number of data

Y=matrix(y.dat[(p+1):n.all],ncol=1)
n=length(Y)

if (p == 1) {
  Fmtx=matrix(c(y.dat[1:(n.all-1)]),nrow=p,byrow=TRUE)
  } else if (p == 2){
  Fmtx=matrix(c(y.dat[2:(n.all-1)],y.dat[1:(n.all-2)]),nrow=p,byrow=TRUE)
  }else if (p == 3){
    Fmtx=matrix(c(y.dat[3:(n.all-1)],y.dat[2:(n.all-2)],y.dat[1:(n.all-3)]),nrow=p,byrow=TRUE)}

## posterior inference

## set the prior
m0=matrix(rep(0,p),ncol=1)
C0=diag(p)
n0=0.02
d0=0.02

## calculate parameters that will be reused in the loop
e=Y-t(Fmtx)%*%m0
Q=t(Fmtx)%*%C0%*%Fmtx+diag(n)
Q.inv=chol2inv(chol(Q))
A=C0%*%Fmtx%*%Q.inv
m=m0+A%*%e
C=C0-A%*%Q%*%t(A)
n.star=n+n0
d.star=t(Y-t(Fmtx)%*%m0)%*%Q.inv%*%(Y-t(Fmtx)%*%m0)+d0

n.sample=5000

## store posterior samples
nu.sample=rep(0,n.sample)
phi.sample=matrix(0,nrow=n.sample,ncol=p)

for (i in 1:n.sample) {
  set.seed(i)
  nu.new=1/rgamma(1,shape=n.star/2,rate=d.star/2)
  nu.sample[i]=nu.new
  phi.new=rmvnorm(1,mean=m,sigma=nu.new*C)
  phi.sample[i,]=phi.new
}

par(mfrow=c(1,4))
hist(phi.sample[,1],freq=FALSE,xlab=expression(phi[1]),main="",ylim=c(0,6.4))
lines(density(phi.sample[,1]),type='l',col='red')
hist(phi.sample[,2],freq=FALSE,xlab=expression(phi[2]),main="",ylim=c(0,6.4))
lines(density(phi.sample[,2]),type='l',col='red')
hist(phi.sample[,3],freq=FALSE,xlab=expression(phi[3]),main="",ylim=c(0,6.4))
lines(density(phi.sample[,3]),type='l',col='red')
hist(nu.sample,freq=FALSE,xlab=expression(nu),main="")
lines(density(nu.sample),type='l',col='red')

## get in sample prediction
post.pred.y=function(s){
  
  beta.cur=matrix(phi.sample[s,],ncol=1)
  nu.cur=nu.sample[s]
  mu.y=t(Fmtx)%*%beta.cur
  sapply(1:length(mu.y), function(k){rnorm(1,mu.y[k],sqrt(nu.cur))})
  
  
}

y.post.pred.sample=sapply(1:5000, post.pred.y)

s=1

beta.cur=matrix(phi.sample[s,],ncol=1)
nu.cur=nu.sample[s]
mu.y=t(Fmtx)%*%beta.cur
sapply(1:length(mu.y), function(k){rnorm(1,mu.y[k],sqrt(nu.cur))})

## show the result
summary.vec95=function(vec){
  c(unname(quantile(vec,0.025)),mean(vec),unname(quantile(vec,0.975)))
}

summary.y=apply(y.post.pred.sample,MARGIN=1,summary.vec95)

par(mfrow=c(1,1))

plot(Y,type='b',xlab='Time',ylab='',ylim=c(0,35),pch=16)
lines(summary.y[2,],type='b',col='grey',lty=2,pch=4)
lines(summary.y[1,],type='l',col='purple',lty=3)
lines(summary.y[3,],type='l',col='purple',lty=3)
legend("topright",legend=c('Truth','Mean','95% C.I.'),lty=1:3,col=c('black','grey','purple'),
       horiz = T,pch=c(16,4,NA))

# Get Model Predictions

## the prediction function

y_pred_h_step=function(h.step,s){
  phi.cur=matrix(phi.sample[s,],ncol=1)
  nu.cur=nu.sample[s]
  y.cur=c(y.post.pred.sample[n],y.post.pred.sample[n-1],y.post.pred.sample[n-2])
  y.pred=rep(0,h.step)
  for (i in 1:h.step) {
    mu.y=sum(y.cur*phi.cur)
    y.new=rnorm(1,mu.y,sqrt(nu.cur))
    y.pred[i]=y.new
    y.cur=c(y.new,y.cur)
    y.cur=y.cur[-length(y.cur)]
  }
  return(y.pred)
}

set.seed(1)
y.post.pred.ahead=sapply(1:5000, function(s){y_pred_h_step(h.step=3,s=s)})

par(oma=c(0,0,2,0), mfrow=c(1,1))
plot(y.new~c(2018,2019,2020), xlab = "Year", ylab = "Earthquakes", type='b')
lines(rowMeans(y.post.pred.ahead)~c(2018,2019,2020),type='b',col='red',lty=3)
legend("topright",legend=c('Truth','Single Component','Mixture'),lty=c(1,3,3),col=c('black','red','blue'),
       horiz = T,pch=c(1,1,1))
title(main = 'Earthquake Prediction vs. Actual Values', outer = TRUE)

# AR Model Mixture

# install.packages("MCMCpack")
library(MCMCpack)
library(mvtnorm)

##
p=3 ## order of AR process
K=2 ## number of mixing component
# Y=matrix(y[3:200],ncol=1) ## y_{p+1:T}

# Fmtx=matrix(c(y[2:199],y[1:198]),nrow=2,byrow=TRUE) ## design matrix F
# Fmtx=matrix(c(y.dat[2:(n.all-1)],y.dat[1:(n.all-2)]),nrow=p,byrow=TRUE)

# n=length(Y) ## T-p

## prior hyperparameters
m0=matrix(rep(0,p),ncol=1)
C0=10*diag(p)
C0.inv=0.1*diag(p)
n0=2
d0=2
a=rep(1,K)

# Sample Functions

sample_omega=function(L.cur){
  n.vec=sapply(1:K, function(k){sum(L.cur==k)})
  rdirichlet(1,a+n.vec)
}

sample_L_one=function(beta.cur,omega.cur,nu.cur,y.cur,Fmtx.cur){
  prob_k=function(k){
    beta.use=beta.cur[((k-1)*p+1):(k*p)]
    omega.cur[k]*dnorm(y.cur,mean=sum(beta.use*Fmtx.cur),sd=sqrt(nu.cur[k]))
  }
  prob.vec=sapply(1:K, prob_k)
  L.sample=sample(1:K,1,prob=prob.vec/sum(prob.vec))
  return(L.sample)
}

sample_L=function(y,x,beta.cur,omega.cur,nu.cur){
  L.new=sapply(1:n, function(j){sample_L_one(beta.cur,omega.cur,nu.cur,y.cur=y[j,],Fmtx.cur=x[,j])})
  return(L.new)
}

sample_nu=function(k,L.cur){
  idx.select=(L.cur==k)
  n.k=sum(idx.select)
  if(n.k==0){
    d.k.star=d0
    n.k.star=n0
  }else{
    y.tilde.k=Y[idx.select,]
    Fmtx.tilde.k=Fmtx[,idx.select]
    e.k=y.tilde.k-t(Fmtx.tilde.k)%*%m0
    Q.k=t(Fmtx.tilde.k)%*%C0%*%Fmtx.tilde.k+diag(n.k)
    Q.k.inv=chol2inv(chol(Q.k))
    d.k.star=d0+t(e.k)%*%Q.k.inv%*%e.k
    n.k.star=n0+n.k
  }
  
  1/rgamma(1,shape=n.k.star/2,rate=d.k.star/2)
}

sample_beta=function(k,L.cur,nu.cur){
  nu.use=nu.cur[k]
  idx.select=(L.cur==k)
  n.k=sum(idx.select)
  if(n.k==0){
    m.k=m0
    C.k=C0
  }else{
    y.tilde.k=Y[idx.select,]
    Fmtx.tilde.k=Fmtx[,idx.select]
    e.k=y.tilde.k-t(Fmtx.tilde.k)%*%m0
    Q.k=t(Fmtx.tilde.k)%*%C0%*%Fmtx.tilde.k+diag(n.k)
    Q.k.inv=chol2inv(chol(Q.k))
    A.k=C0%*%Fmtx.tilde.k%*%Q.k.inv
    m.k=m0+A.k%*%e.k
    C.k=C0-A.k%*%Q.k%*%t(A.k)
  }
  
  rmvnorm(1,m.k,nu.use*C.k)
}

# The Gibbs Sampler

## number of iterations
nsim=20000

## store parameters

beta.mtx=matrix(0,nrow=p*K,ncol=nsim)
L.mtx=matrix(0,nrow=n,ncol=nsim)
omega.mtx=matrix(0,nrow=K,ncol=nsim)
nu.mtx=matrix(0,nrow=K,ncol=nsim)

## initial value

beta.cur=rep(0,p*K)
L.cur=rep(1,n)
omega.cur=rep(1/K,K)
nu.cur=rep(1,K)

# Run Gibbs Sampler

## Gibbs Sampler
for (i in 1:nsim) {
  set.seed(i)
  
  ## sample omega
  omega.cur=sample_omega(L.cur)
  omega.mtx[,i]=omega.cur
  
  ## sample L
  L.cur=sample_L(Y,Fmtx,beta.cur,omega.cur,nu.cur)
  L.mtx[,i]=L.cur
  
  ## sample nu
  nu.cur=sapply(1:K,function(k){sample_nu(k,L.cur)})
  nu.mtx[,i]=nu.cur
  
  ## sample beta
  beta.cur=as.vector(sapply(1:K, function(k){sample_beta(k,L.cur,nu.cur)}))
  beta.mtx[,i]=beta.cur
  
  ## show the numer of iterations 
  if(i%%1000==0){
    print(paste("Number of iterations:",i))
  }
  
}

RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

# Check the Posterior Result

sample.select.idx=seq(10001,20000,by=1)

post.pred.y.mix=function(idx){
  
  k.vec.use=L.mtx[,idx]
  beta.use=beta.mtx[,idx]
  nu.use=nu.mtx[,idx]
  
  
  get.mean=function(s){
    k.use=k.vec.use[s]
    sum(Fmtx[,s]*beta.use[((k.use-1)*p+1):(k.use*p)])
  }
  get.sd=function(s){
    k.use=k.vec.use[s]
    sqrt(nu.use[k.use])
  }
  mu.y=sapply(1:n, get.mean)
  sd.y=sapply(1:n, get.sd)
  sapply(1:length(mu.y), function(k){rnorm(1,mu.y[k],sd.y[k])})
  
}


y.post.pred.sample=sapply(sample.select.idx, post.pred.y.mix)

summary.vec95=function(vec){
  c(unname(quantile(vec,0.025)),mean(vec),unname(quantile(vec,0.975)))
}

summary.y_2=apply(y.post.pred.sample,MARGIN=1,summary.vec95)

par(mfrow=c(1,1))
plot(Y,type='b',xlab='Time',ylab='',ylim=c(0,35),pch=16)
lines(summary.y[2,],type='b',col='grey',lty=2,pch=4)
lines(summary.y[1,],type='l',col='purple',lty=3)
lines(summary.y[3,],type='l',col='purple',lty=3)
legend("topright",legend=c('Truth','Mean','95% C.I.'),lty=1:3,col=c('black','grey','purple'),
       horiz = T,pch=c(16,4,NA))

# Fit Parameter Distributions

par(mfrow=c(1,4))
hist(phi.sample[,1],freq=FALSE,xlab=expression(phi[1]),main="",ylim=c(0,6.4))
lines(density(phi.sample[,1]),type='l',col='red')
hist(phi.sample[,2],freq=FALSE,xlab=expression(phi[2]),main="",ylim=c(0,6.4))
lines(density(phi.sample[,2]),type='l',col='red')
hist(phi.sample[,3],freq=FALSE,xlab=expression(phi[3]),main="",ylim=c(0,6.4))
lines(density(phi.sample[,3]),type='l',col='red')
hist(nu.sample,freq=FALSE,xlab=expression(nu),main="")
lines(density(nu.sample),type='l',col='red')

# Get Model Predictions

## the prediction function

y_pred_h_step_mixture=function(h.step,s){
  phi.cur=matrix(t(beta.mtx[1:3,s]),ncol=1)
  nu.cur=nu.mtx[1,][s]
  y.cur=c(y.post.pred.sample[n],y.post.pred.sample[n-1],y.post.pred.sample[n-2])
  y.pred=rep(0,h.step)
  for (i in 1:h.step) {
    mu.y=sum(y.cur*phi.cur)
    y.new=rnorm(1,mu.y,sqrt(nu.cur))
    y.pred[i]=y.new
    y.cur=c(y.new,y.cur)
    y.cur=y.cur[-length(y.cur)]
  }
  return(y.pred)
}

set.seed(1)
y.post.pred.ahead_mixture=sapply(1:5000, function(s){y_pred_h_step_mixture(h.step=3,s=s)})

mix2.steps = c(1)
i = 1
if (i %in% mix2.steps) {
  m = 2
} else {
  m = 1
}

y_pred_h_step_mixture_2=function(h.step, s, mix2.steps = c(2)){
  i = 1
  if (i %in% mix2.steps) {
    m = 2
  } else {
    m = 1
  }
  phi.cur=matrix(t(beta.mtx[(1+(m-1)*3):(3*m),s]),ncol=1)
  nu.cur=nu.mtx[1*m,][s]
  y.cur=c(y.post.pred.sample[n],y.post.pred.sample[n-1],y.post.pred.sample[n-2])
  y.pred=rep(0,h.step)
  for (i in 1:h.step) {
    if (i %in% mix2.steps) {
      m = 2
    } else {
      m = 1
    }
    phi.cur=matrix(t(beta.mtx[(1+(m-1)*3):(3*m),s]),ncol=1)
    nu.cur=nu.mtx[1*m,][s]
    mu.y=sum(y.cur*phi.cur)
    y.new=rnorm(1,mu.y,sqrt(nu.cur))
    y.pred[i]=y.new
    y.cur=c(y.new,y.cur)
    y.cur=y.cur[-length(y.cur)]
  }
  return(y.pred)
}

set.seed(1)

par(oma=c(0,0,2,0), mfrow=c(1,1))
plot(y.new~c(2018,2019,2020), xlab = "Year", ylab = "Earthquakes", type='b', ylim = c(0,20), xaxt="n")
lines(rowMeans(y.post.pred.ahead)~c(2018,2019,2020),type='b',col='blue',lty=3)
lines(rowMeans(y.post.pred.ahead_mixture)~c(2018,2019,2020),type='b',col='red',lty=3)
lines(rowMeans(y.post.pred.ahead_mixture_2)~c(2018,2019,2020),type='b',col='orange',lty=3)
legend("topright",legend=c('Truth','Single Component','Mixture Component 1 Only','Mixture Component 2 Year 2'),lty=c(1,3,3,3),col=c('black','blue','red','orange'),
       horiz = T,pch=c(1,1,1,1))
title(main = 'Earthquake Prediction vs. Actual Values', outer = TRUE)
x=c("2018", "2019", "2020")
axis(1, at=floor(seq(2018,2020,length=3)), labels=x)


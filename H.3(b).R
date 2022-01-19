#install.packages("MCMCpack")
library(MCMCpack)

#------------ H3.a  Bootstrap particle filter:START         ------------#
#sigma2: sigma^2
#beta2 : beta^2
BootstrapF.H3=function(Y,N,T,phi,sigma2,beta2){
  sigma=sqrt(sigma2)
  beta=sqrt(beta2)
  Xt=matrix(NA,nrow=N,ncol=T)
  Wt=matrix(NA,nrow=N,ncol=T)
  At=matrix(NA,nrow=N,ncol=T)
  log.p=matrix(NA,nrow=N,ncol=T)
  # Initialization t=0
  X0=rnorm(N,mean=0,sd=sigma)
  w0=rep(1/N,N)
  #--- t=1:start
  #resample
  At[,1]=sample(1:N,size=N,prob=w0,replace = TRUE)
  #propogate
  Xt[,1]=rnorm(N,mean=phi*X0[At[,1]],sd=sigma) # X_1{i}|X_0^{a1,i} ~N(AX0{a1,i},Q)
  #Weight & normalize
  sd.yt=beta*sqrt(exp(Xt[,1]))
  log.p[,1]=dnorm(Y[1],mean=0,sd=sd.yt,log=TRUE) #logP(Yt|Xt)
  log.Pmax=max(log.p[,1])                        #max(logP(Yt|Xt))
  Wt[,1]=exp(log.p[,1]-log.Pmax)/sum(exp(log.p[,1]-log.Pmax))
  #--- t=1:end  
  #---------- for t=2 to T :start ------
  if( T>1){
    for (t in 2:T){
      #resample
      At[,t]=sample(1:N,size=N,prob=Wt[,t-1],replace = TRUE)
      #propogate
      Xt[,t]=rnorm(N,mean=phi*Xt[At[,t],t-1],sd=sigma) # X_1{i}|X_0^{a1,i} ~N(AX0{a1,i},Q)
      #Weight & normalize
      sd.yt=beta*sqrt(exp(Xt[,t]))
      log.p[,t]=dnorm(Y[t],mean=0,sd=sd.yt,log=TRUE) #logP(Yt|Xt)
      log.Pmax=max(log.p[,t])                        #max(logP(Yt|Xt))
      Wt[,t]=exp(log.p[,t]-log.Pmax)/sum(exp(log.p[,t]-log.Pmax))
    }  
  }
  #---------- for t=2 to T :end ------
  H3.p=exp(log.p)
  log.z.est=sum(log(apply(H3.p,MARGIN=2,FUN=mean)))  #log(z)
  result=list(At=At,Wt=Wt,Xt=Xt,log.p=log.p,log.z.est=log.z.est)
  return(result)  
}


#///////////////////////////////////////////////////////////////////////#
#------------ H3.b  PMH algorithm:START                     ------------#
#   iter  : number of iterations 
# sigma20 : initial value of sigma^2
#  beta20 : initial value of beta^2
#    var1 : variance of gussain random walk of proposal for sigma^2
#    var2 : variance of gussain random walk of proposal for beta^2
#       a : shape prameter of Inverse Gamma for prior of sigma^2/beta^2
#       b : scale prameter of Inverse Gamma for prior of sigma^2/beta^2
PMH.alg=function(Y,N,T,phi,sigma20,beta20,iter, var1,var2,a,b){  
  sd1=sqrt(var1)
  sd2=sqrt(var2)
  sigma2.t=rep(NA,iter)  
  beta2.t=rep(NA,iter)
  log.z.est.t=rep(NA,iter) #log(z)
  # Initialize m=1
  sigma2.t[1]=sigma20  #set sigma2[1]
  beta2.t[1]=beta20   #set beta2[1]
  log.z.est.t[1]=BootstrapF.H3(Y,N,T,phi,sigma2.t[1], beta2.t[1])$log.z.est #bootstrap particle filter.
  
  # create progress bar
  pb <- txtProgressBar(min = 2, max = iter, style = 3) 
  # for m=2 to M(iter)
  for (m in 2:iter){
    if(m==2){
      ptm=proc.time()[3]
    }
    setTxtProgressBar(pb, m)
    sigma2.t[m]=rnorm(1,mean=sigma2.t[m-1],sd=sd1)  #sample sigma2' from Gaussian random walk
    beta2.t[m]=rnorm(1,mean= beta2.t[m-1],sd=sd2)  #sample beta2'  from Gaussian random walk
    if (sigma2.t[m]<=0 | beta2.t[m]<=0){
      sigma2.t[m]= sigma2.t[m-1]
      beta2.t[m]=  beta2.t[m-1]  
      log.z.est.t[m]=log.z.est.t[m-1]
    }else{
      log.z.est.t[m]=BootstrapF.H3(Y,N,T,phi,sigma2.t[m], beta2.t[m])$log.z.est #bootstrap particle filter.
      log.accept.p=log.z.est.t[m]+log(dinvgamma(sigma2.t[m],a,b))+log(dinvgamma(beta2.t[m],a,b))
      log.accept.p=log.accept.p-log.z.est.t[m-1]-log(dinvgamma(sigma2.t[m-1],a,b))-log(dinvgamma(beta2.t[m-1],a,b))# accept.alpha
      log.u=log(runif(1,0,1))
      if (log.u>log.accept.p){               # if u>accept.alpha, reject new proposal
        sigma2.t[m]=   sigma2.t[m-1]
        beta2.t[m]=    beta2.t[m-1] 
        log.z.est.t[m]=log.z.est.t[m-1]
      }
    }
    if(m==12){
      ptm=round((proc.time()[3]-ptm)*iter/60/60/10,2)
      print(list("Running time will be about ",ptm, "hours"))
    }
  }
  result=list(log.z.est.t=log.z.est.t,beta2.t=beta2.t,sigma2.t=sigma2.t)
  return(result)
}

#--------------**********************------------------------------------------------------------
N=200
T=length(Y)
phi=0.985
sigma20=0.1^2
beta20=1^2
iter=5000
var1=0.01^2
var2=0.3^2
a=0.01
b=0.01
set.seed(12389)    
H3b=PMH.alg(Y,N,T,phi,sigma20,beta20,iter, var1,var2,a,b)
save.image("H3b_test10.RData")  
mean(H3b$sigma2.t[-(1:1000)])
mean(H3b$beta2.t[-(1:1000)])

pdf(paste("../Figures/H.3(b)TracePlot.pdf"))
par(mfrow=c(2,1))
plot(1000:iter, H3b$sigma2.t[1000:5000],type="l",  ylab = expression(sigma^2), xlab = "")
abline(h=0.1^2, col="red")
plot(1000:iter, H3b$beta2.t[1000:5000],type="l", ylab = expression(beta^2), xlab = "")
abline(h=1^2, col="red")
dev.off()

pdf(paste("../Figures/H.3(b)Histogram.pdf"))
par(mfrow=c(2,1))
hist(H3b$sigma2.t[1000:5000], main = expression(sigma^2), xlab = "", ylab = "")
abline(v=0.1^2, col="red")
hist(H3b$beta2.t[1000:5000], main = expression(beta^2), xlab = "", ylab = "")
abline(v=1^2, col="red")
dev.off()

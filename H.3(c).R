#///////////////////////////////////////////////////////////////////////#
#------------ H3.c PGA: Bootstrap particle filter:START     ------------#
# x.ref : refrence trajectory x_{0:T}
# length(x.ref)=T+1
PGA.BF=function(Y,N,T,phi,sigma,beta,x.ref){
  Xt=matrix(NA,nrow=N,ncol=T)
  Wt=matrix(NA,nrow=N,ncol=T)
  At=matrix(NA,nrow=N,ncol=T)
  log.p=matrix(NA,nrow=N,ncol=T)
  #Xt2=matrix(NA,nrow=N,ncol=T+1)
  # Initialization t=0
  X0=rnorm(N,mean=0,sd=sigma)
  X0[N]=x.ref[1]
  w0=rep(1/N,N)
  #--- t=1:start
  #resample
  At[,1]=sample(1:N,size=N,prob=w0,replace = TRUE)
  #propogate
  Xt[,1]=rnorm(N,mean=phi*X0[At[,1]],sd=sigma) # X_1{i}|X_0^{a1,i} ~N(AX0{a1,i},Q)
  At[N,1]=N
  Xt[N,1]=x.ref[2]
  #Xt2[,1:2]=cbind(X0[At[,1]],Xt[,1])
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
      At[N,t]=N
      Xt[N,t]=x.ref[t+1]
      #Xt2[,1:(t+1)]=cbind(Xt2[At[,t],1:t],Xt[,t])
      #Weight & normalize
      sd.yt=beta*sqrt(exp(Xt[,t]))
      log.p[,t]=dnorm(Y[t],mean=0,sd=sd.yt,log=TRUE) #logP(Yt|Xt)
      log.Pmax=max(log.p[,t])                        #max(logP(Yt|Xt))
      Wt[,t]=exp(log.p[,t]-log.Pmax)/sum(exp(log.p[,t]-log.Pmax))
    }  
  }
  
  #sample new refrence trajectory
  index.b=sample(1:N,size=1,prob=Wt[,T])
  x.ref.new=rep(NA,T)
  t.ref.new=rep(NA,T)
  t.ref.new[T]=At[index.b,T] 
  x.ref.new[T]=Xt[index.b,T]
  for (t in (T-1):1){
    t.ref.new[t]=At[t.ref.new[t+1],t] 
    x.ref.new[t]=Xt[t.ref.new[t+1],t]
  }
  x.ref.new=c(X0[t.ref.new[1]],x.ref.new)
  #sum(abs(x.ref.new-Xt2[index.b,]))
  #---------- for t=2 to T :end ------
  #H3.p=exp(log.p)
  #z.est=sum(log(apply(H3.p,MARGIN=2,FUN=mean))) 
  #result=list(At=At,Wt=Wt,Xt=Xt,log.p=log.p,x.ref.new=x.ref.new)
  #return(result)  
  return(x.ref.new)  
}



#///////////////////////////////////////////////////////////////////////#
#------------ H3.c  PGS algorithm:START                     ------------#
#   iter  : number of iterations 
# sigma20 : initial value of sigma^2
#  beta20 : initial value of beta^2
#       a : shape prameter of Inverse Gamma for prior of sigma^2/beta^2
#       b : scale prameter of Inverse Gamma for prior of sigma^2/beta^2
#   x.ref : refrence trajectory x_{0:T}
PGS.alg=function(Y,N,T,phi,sigma20,beta20,x.ref,iter,a,b){  
  sigma2.t=rep(NA,iter)  
  beta2.t=rep(NA,iter)
  
  # Initialize m=1
  sigma2.t[1]=sigma20  #set sigma2[1]
  beta2.t[1]=beta20   #set beta2[1]
  
  # create progress bar
  pb <- txtProgressBar(min = 2, max = iter, style = 3) 
  # for m=2 to M(iter)
  for (m in 2:iter){
    setTxtProgressBar(pb, m)
    a1=a+T/2
    b1=b+sum((x.ref[-1]-phi*x.ref[-(T+1)])^2)/2
    b2=b+sum(Y^2*exp(-x.ref[-1]))/2
    sigma2.t[m]=rinvgamma(1,a1,b1)                #sample sigma2' from Gaussian random walk
    beta2.t[m]=rinvgamma(1,a1,b2)                #sample beta2'  from Gaussian random walk
    x.ref=PGA.BF(Y,N,T,phi,sqrt(sigma2.t[m]),sqrt(beta2.t[m]),x.ref)
    #PGA.BF(Y,N,T,phi,sigma,beta,x.ref)
  }
  result=list(beta2.t=beta2.t,sigma2.t=sigma2.t)
  return(result)
}

#--------------**********************------------------------------------------------------------
sigma20=0.1^2
beta20=1^2
phi=0.985
iter=5000
N=200
T=length(Y)
T
a=0.01
b=0.01
set.seed(123867)
x.ref=rnorm((T+1),0,0.16)
H3c=PGS.alg(Y,N,T,phi,sigma20,beta20,x.ref,iter,a,b)
save.image("H3c_test1new.RData")

pdf(paste("../Figures/H.3(c)TracePlot.pdf"))
par(mfrow=c(2,1))
plot(1000:iter, H3c$sigma2.t[1000:iter], type = "l", ylab = expression(sigma^2), xlab = "")
abline(h=0.1^2, col="red")
plot(1000:iter, H3c$beta2.t[1000:iter], type = "l", ylab = expression(beta^2), xlab = "")
abline(h=1^2, col="red")
dev.off()

pdf(paste("../Figures/H.3(c)Histogram.pdf"))
par(mfrow=c(2,1))
hist(H3c$sigma2.t[1000:iter], main = expression(sigma^2), xlab = "")
abline(v = 0.1^2, col = "red")
hist(H3c$beta2.t[1000:iter], main = expression(beta^2), xlab = "")
abline(v = 1, col = "red")
dev.off()


plot(H3c$beta2.t[500:2000])
plot(H3c$sigma2.t[500:2000])
mean(H3c$beta2.t[500:2000])
mean(H3c$sigma2.t[500:2000])

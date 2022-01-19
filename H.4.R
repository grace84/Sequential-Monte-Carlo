# H.4 SMC sampler
install.packages("plot3D")
library("plot3D")

x1 <- seq(-1,1,length.out=500)
x2 <- seq(-1,1,length.out=500)
z <- sin(x1*pi)^2*cos(2*pi*x2)^8*exp(-5*(x1^2+x2^2))
scatter3D(x1, x2, z, colvar = NULL, col = "blue",
          pch = 19, cex = 0.5)




## r_{k}(x) propotion to pi(x)^{\tau_k}*p(x)
## take Guassian random walk as proposal, set variance as 0.02^2
## take p(x)=p(x1)*p(x2) ,p(x1)=unif(0,1),p(x2)=unif(0,1)
rm(list=ls())
## pi.tilde(x1,x2) propotional to pi(x1,x2)
pi.tilde=function(x1,x2){
  y=(sin(x1*pi)^2)*(cos(x2*2*pi)^8)*exp(-5*(x1^2+x2^2))  
  #y=0
  #if((x1*(1-x1))>=0 & (x2*(1-x2))>=0){
  return(y)
}

# sig1: standard deviation of proposal of x1
# sig2: standard deviation of proposal of x2
SMC.sampler=function(N,K,sig1,sig2){
  X1t=matrix(NA,nrow=N,ncol=K)
  X2t=matrix(NA,nrow=N,ncol=K)
  Wt=matrix(NA,nrow=N,ncol=K)
  At=matrix(NA,nrow=N,ncol=K)
  log.p=matrix(NA,nrow=N,ncol=K)
  log.z=rep(NA,K)
  Ess1=rep(NA,K)
  Ess2=rep(NA,K)
  # Initialization t=0
  X01=runif(N,0,1)
  X02=runif(N,0,1)
  w0=rep(1/N,N)
  #--- k=1:start
  #resample : update weights & resample
  k=1
  tk=k/K
  log.p[,k]=log(w0)+log(pi.tilde(X01,X02))*tk # log.p=log(w0)+log(incre)
  log.z[k]=log(sum(exp(log.p[,k])))
  #log.z[k]=sum(log.p[,k])+log.z[k-1]
  log.Pmax=max(log.p[,k])
  Wt[,k]=exp(log.p[,k]-log.Pmax)/sum(exp(log.p[,k]-log.Pmax))
  Ess=1/sum(Wt[,k]^2)
  Ess1[k]=Ess
  Ess2[k]=Ess
  if(Ess/N<0.7){                #if ture, resample
    At[,k]=sample(1:N,size=N,prob=Wt[,k],replace = TRUE)
    Wt[,k]=rep(1/N,N)
    Ess2[k]=1/sum(Wt[,k]^2)
  }else{
    At[,k]=1:N
  }
  #propogate
  for(n in 1:N){
    X1t[n,k]=rnorm(1,X01[At[n,k]],sd=sig1)
    X2t[n,k]=rnorm(1,X02[At[n,k]],sd=sig2) 
    index.rsp=((X2t[n,k]*(1-X2t[n,k]))>0)*1+((X1t[n,k]*(1-X1t[n,k]))>0)*1
    if(index.rsp!=2){
      X1t[n,k]=X01[At[n,k]]  
      X2t[n,k]=X02[At[n,k]]
    }else{
      x1=X1t[n,k]
      x2=X2t[n,k]
      x1p=X01[At[n,k]]
      x2p=X02[At[n,k]]
      log.a=tk*log(pi.tilde(x1,x2))-tk*log(pi.tilde(x1p,x2p))
      log.u=log(runif(1,0,1))
      if(log.u>log.a) #reject the proposal
        X1t[n,k]=X01[At[n,k]]  
      X2t[n,k]=X02[At[n,k]] 
    }
  }
  #end k=1-----------------------
  
  # for k=2 to K
  # create progress bar
  pb <- txtProgressBar(min = 2, max = K, style = 3) 
  for (k in 2:K){
    if(k==2){
      ptm=proc.time()[3]
    }
    setTxtProgressBar(pb, k)
    tk=k/K
    log.p[,k]=log(Wt[,k-1])+log(pi.tilde(X1t[,k-1],X2t[,k-1]))/K
    
    
    # log.p=log(w0)+log(incre)
    log.z[k]=log(sum(exp(log.p[,k])))+log.z[k-1]
    log.Pmax=max(log.p[,k])
    Wt[,k]=exp(log.p[,k]-log.Pmax)/sum(exp(log.p[,k]-log.Pmax))
    Ess=1/sum(Wt[,k]^2)
    Ess1[k]=Ess
    Ess2[k]=Ess
    if(Ess/N<0.7){                #if ture, resample
      At[,k]=sample(1:N,size=N,prob=Wt[,k],replace = TRUE)
      Wt[,k]=rep(1/N,N)
      Ess2[k]=1/sum(Wt[,k]^2)
    }else{
      At[,k]=1:N
    }
    #start: propogate ---MH-------------
    for(n in 1:N){
      X1t[n,k]=rnorm(1,X1t[At[n,k],k-1],sd=sig1)
      X2t[n,k]=rnorm(1,X2t[At[n,k],k-1],sd=sig2) 
      # index.rsp=((X2t[n,k]*(1-X2t[n,k]))>0)*1+((X1t[n,k]*(1-X1t[n,k]))>0)*1
      if ((X1t[n,k] > -1 & X1t[n,k] < 1) & (X2t[n,k] > -1 & X2t[n,k] < 1)){
        index.rsp = 2
      }
      #tip : if x1 not in (0,1) or x2 not in (0,1) ,reject the proposal directly
      if(index.rsp!=2){ 
        X1t[n,k]=X1t[At[n,k],k-1]    
        X2t[n,k]=X2t[At[n,k],k-1]
      }else{
        x1=X1t[n,k]
        x2=X2t[n,k]
        x1p=X1t[At[n,k],k-1]
        x2p=X2t[At[n,k],k-1]
        log.a=tk*log(pi.tilde(x1,x2))-tk*log(pi.tilde(x1p,x2p))
        log.u=log(runif(1,0,1))
        if(log.u>log.a) # if u > accept.alpha , reject the proposal
          X1t[n,k]=X1t[At[n,k],k-1] 
        X2t[n,k]=X2t[At[n,k],k-1]
      }# end if-else
    } #end: propogate ---MH-------------
    if(k==12){
      ptm=round((proc.time()[3]-ptm)*K/60/60/10,2)
      print(list("Running time will be about ",ptm, "hours"))
    }
  }#---------------end k from 2 to K
  
  result=list(X01=X01,X02=X02,At=At,Wt=Wt,X1t=X1t,X2t=X2t,log.z=log.z,Ess1=Ess1,Ess2=Ess2)
  return(result)  
}
#-----------------------H4: SMC.sampler end -----------------
N=100
K=8000
sig1=sig2=0.02
set.seed(123456)
H4=SMC.sampler(N,K,sig1,sig2)
tail(exp(H4$log.z))
tail((H4$log.z))

At=H4$At
X1t=H4$X1t
X2t=H4$X2t

Track4.T=matrix(NA,nr=N,nc=K)
Track4.X1=matrix(NA,nr=N,nc=K)
Track4.X2=matrix(NA,nr=N,nc=K)
Track4.T[,K]=At[,K] 
Track4.X1[,K]=X1t[,K]
Track4.X2[,K]=X2t[,K]
for (t in (K-1):1){
  Track4.T[,t]=At[Track4.T[,t+1],t] 
  Track4.X1[,t]=X1t[Track4.T[,t+1],t]
  Track4.X2[,t]=X2t[Track4.T[,t+1],t]
}

X1.all=cbind(H4$X01,X1t)
X2.all=cbind(H4$X02,X2t)
#plot location of particles at time 0,1000,4000, and K
index.piece=c(0,1000,4000,K)+1
range(X1.all[,index.piece])

# gname = c("H4a.eps",sep="")  
# postscript(gname,width=12,height=12,horizontal = FALSE, onefile = FALSE, paper = "special")
# par(mfrow=c(2,2),oma=c(0.2,3.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1.5,las=1,mgp=c(1,0.5,0),adj=0.5)
pdf(paste("./Figures/H.4(a)LocationOfParticles_X1.pdf"), width = 8, height = 6)
par(mfrow = c(2,1))
plot(X1.all[1,index.piece],xlab="",ylab="",
     ylim=c(0,1), type="p",pch=20,col="black",xaxt="n", cex.axis=1)
for (n in 2:N){
  lines(X1.all[n,index.piece],
        type="p",pch=20,col="black")
}
mtext("X1", las=0,side = 2, cex = 1, line = 2.5, outer = FALSE)
mtext(expression(k), side = 1, cex = 1, line = 2, outer = FALSE)
axis(1, at =1:4, labels =(index.piece-1), las=1, cex.axis=1)

plot(density(X1.all[,1]),col=1,ylim=c(0,5),xlim=c(-1,1),
     xlab="",ylab="",main="", cex.axis=1)
mtext("est. density of X1", las=0,side = 2, cex = 1, line = 2, outer = FALSE)
lines(density(X1.all[,1001]),col=2)
lines(density(X1.all[,4001]),col=3)
lines(density(X1.all[,K+1]),col=4)
legend("topright",c("k=0","k=1000","k=4000","k=8000"),
       cex=1,lty=c(1,1,1,1), bty="n",lwd=c(1.2,1.2,1.2,1.2),col=c(1,2,3,4))
dev.off()


pdf(paste("./Figures/H.4(a)LocationOfParticles_X2.pdf"), width = 8, height = 6)
par(mfrow = c(2,1))
plot(X2.all[1,index.piece], xlab="", ylab="",
     ylim=c(0,1),type="p",pch=20,col="black",xaxt="n", cex.axis=1)
for (n in 2:N){
  lines(X2.all[n,index.piece],
        type="p",pch=20,col="black")
}
mtext("X2", las=0,side = 2, cex = 1, line = 2.5, outer = FALSE)
mtext(expression(k), side = 1, cex = 1, line = 2, outer = FALSE)
axis(1, at =1:4, labels =(index.piece-1), las=1, cex.axis=1)

plot(density(X2.all[,1]),col=1,ylim=c(0,5),xlim=c(-1,1),
     xlab="",ylab="",main="", cex.axis=1)
mtext("est. density of X2", las=0,side = 2, cex = 1, line = 2.5, outer = FALSE)
lines(density(X2.all[,1001]),col=2)
lines(density(X2.all[,4001]),col=3)
lines(density(X2.all[,K+1]),col=4)
dev.off()

#--------------------------------------------
par(mfrow=c(2,4))
hist(X1.all[,1],xlab="",ylab="",main="",xaxt="n",xlim=c(0,1),ylim=c(0,50))
mtext("X1", las=0,side = 2, cex = 1.2, line = 2.5, outer = FALSE)
mtext("k=0", side = 1, cex = 1.5, line = 2, outer = FALSE)

hist(X1.all[,1000],xlab="",ylab="",main="",xaxt="n",xlim=c(0,1),ylim=c(0,50))
mtext("k=1000", side = 1, cex = 1.5, line = 2, outer = FALSE)

hist(X1.all[,4000],xlab="",ylab="",main="",xaxt="n",xlim=c(0,1),ylim=c(0,50))
mtext("k=4000", side = 1, cex = 1.5, line = 2, outer = FALSE)

hist(X1.all[,8000],xlab="",ylab="",main="",xaxt="n",xlim=c(0,1),ylim=c(0,50))
mtext("k=8000", side = 1, cex = 1.5, line = 2, outer = FALSE)


hist(X2.all[,1],xlab="",ylab="",main="",xaxt="n",xlim=c(0,1),ylim=c(0,50))
mtext("X2", las=0,side = 2, cex = 1.2, line = 2.5, outer = FALSE)
mtext("k=0", side = 1, cex = 1.5, line = 2, outer = FALSE)

hist(X2.all[,1001],xlab="",ylab="",main="",xaxt="n",xlim=c(0,1),ylim=c(0,50))
mtext("k=1000", side = 1, cex = 1.5, line = 2, outer = FALSE)

hist(X2.all[,4001],xlab="",ylab="",main="",xaxt="n",xlim=c(0,1),ylim=c(0,50))
mtext("k=4000", side = 1, cex = 1.5, line = 2, outer = FALSE)

hist(X2.all[,8001],xlab="",ylab="",main="",xaxt="n",xlim=c(0,1),ylim=c(0,50))
mtext("k=8000", side = 1, cex = 1.5, line = 2, outer = FALSE)



# gname = c("H4bc.eps",sep="")  
# postscript(gname,width=12,height=12,horizontal = FALSE, onefile = FALSE, paper = "special")
# par(mfrow=c(2,1),oma=c(0.2,3.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1.5,las=1,mgp=c(1,0.5,0),adj=0.5)
pdf(paste("./Figures/H.4(b)ESSvsK.pdf"), width = 12, height = 5)
plot(x=1:K,y=H4$Ess1,xlab="",ylab="",
     ylim=c(60,100),type="l",pch=20, col="black", cex.axis = 1)
mtext("ESS", side = 2,las=0, cex = 1, line = 3, outer = FALSE)
mtext(expression(k), side = 1, cex = 1, line = 2, outer = FALSE)
abline(h=70,col="red",lty=4)
dev.off()

plot(x=1:K,y=exp(H4$log.z),xlab="",ylab="",
     type="l",pch=20,col="black")
mtext("Est. Z", side = 2,las=0, cex = 1.5, line = 3, outer = FALSE)
mtext(expression(k), side = 1, cex = 1.5, line = 2, outer = FALSE)
abline(h=70,col="red")


tail(round(exp(H4$log.z),5)) # Estimated normalizing constant
range(H4$Ess1)
plot(H4$Ess1,type="l")
plot(exp(H4$log.z),type="l")
logZ.est=rep(NA,10)
for(j in 1:10){
 H4=SMC.sampler(N,K,sig1,sig2)
 logZ.est[j]=H4$log.z[K]
}
boxplot(logZ.est)
#plot(exp(H4$log.z),type="l")
#tail((H4$log.z))
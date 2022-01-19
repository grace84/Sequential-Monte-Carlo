
set.seed(666)

### (a) Simulate the model to produce T = 2 000 synthetic measurements y1:T

T     <- 2000
Y     <- matrix(0, nrow=1, ncol=T)
X     <- matrix(0, nrow=1, ncol=T + 1)

x[1] <- rnorm(1, 0, 1)  # X0 ~ N(0, 1)

#process noise (standard deviation)
Vt <- sqrt(0.5)

#measurement noise (standard deviation)
Et <- sqrt(0.1)

t = 1

while (t < T){
  
  t = t + 1
  
  X[t] <- 0.9*X[t-1] + rnorm(1, 0, Vt) 
  Y[t] <- 1.3*X[t] + rnorm(1, 0, Et)

}

pdf(paste("../Figures/H.2(a)_SimulatedY&X.pdf"), width = 12, height = 6)
plot(1:T, Y, type = "l", ylab = expression(y[1:T]), lwd =1)
lines(0:T, X, type = "l", col = "red")
legend("topleft", lty=c(1, 2, 2),
       col=c("black", "red"),
       legend=c(expression(y[1:T]), expression(x[0:T])),
       bty="n", y.intersp=1, cex=1)
dev.off()

pdf(paste("../Figures/H.2(a)_SimulatedX.pdf"), width = 12, height = 6)
plot(0:T, X, type = "l",  ylab = expression(x[0:T]))
dev.off()

### (b) The Kalman Filter for the model

A <- matrix(0.9, 1, 1)
C <- matrix(1.3, 1, 1)
Q <- matrix(0.5, 1, 1)
R <- matrix(0.1, 1, 1)

x.hat <- matrix(0, 1, T+1) # Start from x_{0|0}
P.hat <- matrix(0, 1, T+1) # Start from P_{0|0}
 
t = 1
x.hat[1] <- 0
P.hat[1] <- 0

while (t < T){
  
  t = t + 1
  
  P.hat.t.tm1 <- A %*% P.hat[t-1] %*% t(A) + Q
  K <- P.hat.t.tm1 %*% t(C) %*% solve(C %*% P.hat.t.tm1 %*% t(C) + R)
  
  P.hat[t] <- P.hat.t.tm1 - K %*% C %*% P.hat.t.tm1
  x.hat[t] <- A %*% x.hat[t-1] + K %*% (Y[t] - C %*% A %*% x.hat[t-1])
}

pdf(paste("../Figures/H.2(b)KF_StateMeanEstimation.pdf"), width = 12, height = 5)
par(mfrow=c(1,1))
t <- 1999
plot(0:t, X[1:(t+1)], type = "l", main = "Kalman filter solution", 
     ylab ="X", xlab = "Time", col=c(gray(level=.5)))
lines(0:t, x.hat[1:(t+1)], col = "blue", lwd = 1, lty = 2)
legend("topleft", lty=c(1, 2, 2),
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "KF"),
       bty="n", y.intersp=1.2, cex=1)
dev.off()

### (c) Implement the bootstrap particle filter (with multinomial resampling) and compare its mean and variance
### estimates to the Kalman filter results.
set.seed(666)
N = 5000 # number of particles

x.sample <- matrix(0, N, T)
norm.weights <- matrix(0, N, T)
ancestor.ind <- matrix(0, N, T)

# t=0
x.sample[,1] <- rnorm(N, 0, 1)
norm.weights[,1] <- 1/N

# t=1 to T

t = 1

while(t < T){
  
  t = t+1
  
  # Resample: the ancestor index a_t^i represents the ancestor of partivle x_t^i at time t-1
  ancestor.ind[,t] <- sample(1:N, N, replace = T, prob = norm.weights[, t-1]) 
  
  # Propogate
  x.sample[,t] <- rnorm(N, mean = 0.9*x.sample[ancestor.ind[,t], t-1], sd = sqrt(0.5))
  
  # Weight
  weight <- dnorm(Y[t], mean = 1.3*x.sample[,t], sd = sqrt(0.1))
  norm.weights[,t] <-  weight/sum(weight)

}

x.means.BPF <- apply((norm.weights*x.sample),MARGIN=2,FUN=sum)
x.variance.BPF <- apply(((x.sample-matrix(x.means.BPF,nr=N,nc=T,byrow=TRUE))^2*norm.weights),MARGIN=2,FUN=sum)

# x.means.BPF <- apply(x.sample, 2, mean)
# x.variance.BPF <- apply(x.sample, 2, var)

## Report the average absolute difference of mean and variance
mean(abs(x.means.BPF-x.hat[1:2000]))
mean(abs(x.variance.BPF - P.hat[1:2000]))

# mean(abs(x.hat[1:(T-1)] - x.means.BPF[2:T]))
# mean(abs(P.hat[1:(T-1)] - x.variance.BPF[2:T]))

par(mfrow=c(1,2))
plot(2:T, abs(x.hat[1:(T-1)] - x.means.BPF[2:T]), type = "l", xlab ="1:T",
     ylim = c(0,3), ylab = "Absolute difference of mean", main = paste("N =",N))
plot(2:T, abs(P.hat[1:(T-1)] - x.variance.BPF[2:T]), type = "l", xlab ="1:T",
     ylim = c(0,3), ylab = "Absolute difference of variance", main = paste("N =",N))

pdf(paste("../Figures/H.2(c)BootsrapPF_StateMeanEstimation.pdf"), width = 12, height = 5)
par(mfrow=c(1,1))
t <- 1998
plot(0:t, X[1:(t+1)], type = "l", main = "Bootstrap Particle Filter", 
     ylab ="X", xlab = "Time", col=c(gray(level=.5)))
# lines(0:t, x.hat[1:(t+1)], col = "blue", lwd = 1, lty = 2) # KF
lines(0:t, x.means.BPF[2:(t+2)], col = "blue", lwd = 1, lty = 2) # PF
legend("topleft", lty=c(1, 2, 2),
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "Bootstrap PF"),
       bty="n", y.intersp=1.2, cex=1)
dev.off()

pdf(paste("../Figures/H.2(c)BootsrapPF_StateVarianceEstimation.pdf"), width = 12, height = 5)
par(mfrow=c(1,1))
t <- 1999
plot(1:t, P.hat[2:(t+1)], type = "l", main = "Kalman filter vs. Bootstrap Particle Filter", 
     ylab ="Variance",  xlab = "Time", col= "red", ylim=c(0.0, 0.1))
lines(1:t, x.variance.BPF[3:(t+2)], col = "blue", lty = 2) # PF
legend("topleft", lty=c(1, 2),
       col=c("red", "blue"),
       legend=c("KF", "Bootstrap PF"),
       bty="n", y.intersp=1.2, cex=1)
dev.off()

### (d) Fully adapted particle filter
set.seed(666)
N = 5000 # number of particles

x.sample <- matrix(0, N, T)
norm.weights <- matrix(0, N, T)
ancestor.ind <- matrix(0, N, T)

# t=0
x.sample[,1] <- rnorm(N, 0, 1)
norm.weights[,1] <- 1/N

# t=1 to T

t = 1

while(t < T){
  
  t = t+1
  
  # Resample: the ancestor index a_t^i represents the ancestor of particle x_t^i at time t-1
  ancestor.ind[,t] <- sample(1:N, N, replace = T, prob = norm.weights[, t-1]) 
  
  # Propogate
  x.sample[,t] <- rnorm(N,
                        mean = 0.687831*Y[t] - 0.095238*x.sample[ancestor.ind[,t], t-1],
                        sd = sqrt(0.0529))

  # Weight
  weight <- dnorm(Y[t], mean = 1.17*x.sample[,t-1], sd = sqrt(0.945))
  norm.weights[,t] <-  weight/sum(weight)
  
}

x.means.FAPF <- apply((norm.weights*x.sample),MARGIN=2,FUN=sum)
x.variance.FAPF <- apply(((x.sample-matrix(x.means.FAPF,nr=N,nc=T,byrow=TRUE))^2*norm.weights),MARGIN=2,FUN=sum)


## Report the average absolute difference of mean and variance

# FAPF vs. KF
mean(abs(x.hat[1:2000] - x.means.FAPF))
mean(abs(P.hat[1:2000] - x.variance.FAPF))

# FAPF vs. BPF
mean(abs(x.means.BPF - x.means.FAPF))
mean(abs(x.variance.BPF - x.variance.FAPF))

## Estimation of mean of state trajectory
pdf(paste("../Figures/H.2(d)FullyAdaptedPF_StateMeanEstimation.pdf"), width = 12, height = 5)
plot(0:t, X[1:(t+1)], type = "l", main = "Fully Adapted Particle Filter", 
     ylab ="X", xlab = "Time", col=c(gray(level=.5)))
lines(0:t, x.means.FAPF[2:(t+2)], col = "blue", lwd = 1, lty = 2) 
legend("topleft", lty=c(1, 2, 2),
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "FAPF"),
       bty="n", y.intersp=1.2, cex=1)
dev.off()


## Estimation of variance of state trajectory
pdf(paste("../Figures/H.2(d)FullyAdaptedPF_StateVarianceEstimation.pdf"), width = 12, height = 5)
t <- 1999
plot(1:t, P.hat[2:(t+1)], type = "l", main = "KF, Bootstrap PF, Fully adapted PF", 
     ylab ="Variance",  xlab = "Time", col= "blue",  lty = 1, ylim=c(0.0, 0.2))
lines(1:t, x.variance.BPF[3:(t+2)], col = "red", lwd = 1, lty = 2) # Bootsrap PF
lines(1:t, x.variance.FAPF[3:(t+2)], col = "orange", lwd = 1, lty = 2) # Fully adapted PF
legend("topleft", lty=c(1, 2, 2),
       col=c("blue", "orange", "red"),
       legend=c("KF", "FAPF", "BPF"),
       bty="n", y.intersp=1.2, cex=1)
dev.off()


### (e) Trace the particle genealogy in fully adapted particle filter

N = 100 # number of particles

x.sample <- matrix(0, N, T)
x.resample <- matrix(0, N, T)
norm.weights <- matrix(0, N, T)
ancestor.ind <- matrix(0, N, T)

# t=0
x.sample[,1] <- rnorm(N, 0, 1)
norm.weights[,1] <- 1/N

# t=1 to T

t = 1

while(t < T){
  
  t = t+1
  
  # Resample: the ancestor index a_t^i represents the ancestor of particle x_t^i at time t-1
  ancestor.ind[,t] <- sample(1:N, N, replace = T, prob = norm.weights[, t-1]) 
  
  # Store the resampled x at time t-1
  x.resample[,t] <- x.sample[ancestor.ind[,t], t-1] 
  
  # Propogate
  x.sample[,t] <- rnorm(N, 
                        mean = 0.687831*Y[t] + 0.095238*x.resample[,t] , 
                        sd = sqrt(0.156085))
  
  # Weight
  weight <- dnorm(Y[t], mean = 1.17*x.sample[,t-1], sd = sqrt(0.945))
  norm.weights[,t] <-  weight/sum(weight)
  
}


### Visualize joint filtering density
pdf(paste("./Figures/H.2(e)FullyAdaptedPF_ParticleGenealogy.pdf"), width = 12, height = 5)

end.time = 100 # End time point for visualizing path degeneracy
plot(x = 1:end.time, y = x.sample[1,1:end.time], pch = 16, col="grey", 
     cex = 0.5, xlab = "Time", ylab = "State", ylim = c(-5,4))

# Add dots
for(i in 2:N){
  points(x = 1:end.time, y = x.sample[i,1:end.time], pch = 16, col="grey", cex = 0.5)
}
# 
# Connect each particle with its ancestor
for (j in 2:end.time){
    segments(x0 = j-1, y0 = x.sample[ancestor.ind[,j], j-1], x1 = j, y1 = x.sample[,j], col = "grey")
  }

# Trace resampled particles
segments(x0 = end.time, y0 = x.sample[, end.time], 
         x1 = end.time-1, y1 = x.sample[ancestor.ind[,end.time], end.time-1], col = "black") 
prev.unique <- unique(ancestor.ind[,end.time])

for (j in (end.time-1):2){
  next.unique <- unique(ancestor.ind[prev.unique,j])
  segments(x0 = j, y0 = x.sample[prev.unique, j], 
           x1 = j-1, y1 = x.sample[next.unique,j-1])
  prev.unique <- next.unique
}
dev.off()

### (f) Systematic resampling

devtools::install_github("jarad/smcUtils")
library("smcUtils")

source("check.weights.R")
source("inverse.cdf.weights.R")

N = 100 # number of particles

x.sample <- matrix(0, N, T)
x.resample <- matrix(0, N, T)
norm.weights <- matrix(0, N, T)
ancestor.ind <- matrix(0, N, T)

# t=0
x.sample[,1] <- rnorm(N, 0, 1)
norm.weights[,1] <- 1/N

# t=1 to T

t = 1

while(t < T){
  
  t = t+1
  
  # Resample: the ancestor index a_t^i represents the ancestor of particle x_t^i at time t-1
  ancestor.ind[,t] <- systematic.resample(weights = norm.weights[, t-1], num.samples = N) 
  
  # Store the resampled x at time t-1
  x.resample[,t] <- x.sample[ancestor.ind[,t], t-1] 
  
  # Propogate
  x.sample[,t] <- rnorm(N, 
                        mean = 0.687831*Y[t] + 0.095238*x.resample[,t] , 
                        sd = sqrt(0.156085))
  
  # Weight
  weight <- dnorm(Y[t], mean = 1.17*x.sample[,t-1], sd = sqrt(0.945))
  norm.weights[,t] <-  weight/sum(weight)
  
}

### Visualize joint filtering density
pdf(paste("./Figures/H.2(f)FullyAdaptedPF_ParticleGenealogy.pdf"), width = 12, height = 5)

end.time = 100 # End time point for visualizing path degeneracy
plot(x = 1:end.time, y = x.sample[1,1:end.time], pch = 16, col="grey", 
     cex = 0.5, xlab = "Time", ylab = "State", ylim = c(-5,4))

# Add dots
for(i in 2:N){
  points(x = 1:end.time, y = x.sample[i,1:end.time], pch = 16, col="grey", cex = 0.5)
}
# 
# Connect each particle with its ancestor
for (j in 2:end.time){
  segments(x0 = j-1, y0 = x.sample[ancestor.ind[,j], j-1], x1 = j, y1 = x.sample[,j], col = "grey")
}

# Trace resampled particles
segments(x0 = end.time, y0 = x.sample[, end.time], 
         x1 = end.time-1, y1 = x.sample[ancestor.ind[,end.time], end.time-1], col = "black") 
prev.unique <- unique(ancestor.ind[,end.time])

for (j in (end.time-1):2){
  next.unique <- unique(ancestor.ind[prev.unique,j])
  segments(x0 = j, y0 = x.sample[prev.unique, j], 
           x1 = j-1, y1 = x.sample[next.unique,j-1])
  prev.unique <- next.unique
}

dev.off()


### (g) ESS-adaptive resampling

N = 100 # number of particles

x.sample <- matrix(0, N, T)
x.resample <- matrix(0, N, T)
norm.weights <- matrix(0, N, T)
ancestor.ind <- matrix(0, N, T)
N.eff <- c()

# t=0
x.sample[,1] <- rnorm(N, 0, 1)
norm.weights[,1] <- 1/N

t = 1

while(t < T){
  
  t = t+1
  
  # Compute effective sample size
  N.eff <- c(1/sum(norm.weights[,t-1]^2),N.eff)
  
  # ESS-updated resample: 
  if(N.eff < N/2){
    ancestor.ind[,t] <- sample(1:N, N, replace = T, prob = norm.weights[, t-1])
    norm.weights[, t-1] <- 1/N
  }else{
    ancestor.ind[,t] <- 1:N
  }
  
  # Store the resampled x at time t-1
  x.resample[,t] <- x.sample[ancestor.ind[,t], t-1] 
  
  # Propogate
  x.sample[,t] <- rnorm(N, 
                        mean = 0.687831*Y[t] + 0.095238*x.resample[,t] , 
                        sd = sqrt(0.156085))
  
  # Weight
  weight <- dnorm(Y[t], mean = 1.17*x.sample[,t-1], sd = sqrt(0.945))
  norm.weights[,t] <-  weight/sum(weight)
  
}

pdf(paste("./Figures/H.2(g)FullyAdaptedPF_ParticleGenealogy.pdf"), width = 12, height = 5)

end.time = 100 # End time point for visualizing path degeneracy
plot(x = 1:end.time, y = x.sample[1,1:end.time], pch = 16, col="grey", 
     cex = 0.5, xlab = "Time", ylab = "State", ylim = c(-5,4))

# Add dots
for(i in 2:N){
  points(x = 1:end.time, y = x.sample[i,1:end.time], pch = 16, col="grey", cex = 0.5)
}
# 
# Connect each particle with its ancestor
for (j in 2:end.time){
  segments(x0 = j-1, y0 = x.sample[ancestor.ind[,j], j-1], x1 = j, y1 = x.sample[,j], col = "grey")
}

# Trace resampled particles
segments(x0 = end.time, y0 = x.sample[, end.time], 
         x1 = end.time-1, y1 = x.sample[ancestor.ind[,end.time], end.time-1], col = "black") 
prev.unique <- unique(ancestor.ind[,end.time])

for (j in (end.time-1):2){
  next.unique <- unique(ancestor.ind[prev.unique,j])
  segments(x0 = j, y0 = x.sample[prev.unique, j], 
           x1 = j-1, y1 = x.sample[next.unique,j-1])
  prev.unique <- next.unique
}
dev.off()


### Plot N.eff/N over time
pdf(paste("./Figures/H.2(g)FullyAdaptedPF_N.effOverN.pdf"), width = 12, height = 5)

plot(1:1999, N.eff/N, type = "l", xlab = "Time")

dev.off()

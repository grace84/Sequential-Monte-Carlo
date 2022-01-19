# H.3. Parameter estimation in the stochastic volatility model

# (a) Compute 10 estimates for each grid point and plot all of them in one figure

Y <- read.csv("seOMXlogreturns2012to2014.csv", header = F)
Y <- Y$V1

grid <- seq(0.1, 1, 0.1)
log.likelihood <- matrix(0, 10, length(grid))    
N = 100 # number of particles
T <- length(Y)

for (i in 1:length(grid)){
  
  phi <- grid[i]
  
  for (j in 1: 10){
    
    sigma <- 0.16
    beta <- 0.7
    
    x.sample <- matrix(0, N, T)
    weight <- matrix(0, N, T)
    norm.weights <- matrix(0, N, T)
    ancestor.ind <- matrix(0, N, T)
    
    # t=0
    x.sample[,1] <- rnorm(N, 0, sigma)
    norm.weights[,1] <- 1/N
    
    # t=1 to T
    
    t = 1
    
    while(t < T){
      
      t = t+1
      
      # Resample: the ancestor index a_t^i represents the ancestor of particle x_t^i at time t-1
      ancestor.ind[,t] <- sample(1:N, N, replace = T, prob = norm.weights[, t-1]) 
      
      # Propogate
      x.sample[,t] <- rnorm(N, 
                            mean = phi*x.sample[ancestor.ind[,t], t-1], 
                            sd = sigma)
      
      # Weight
      weight[,t] <- dnorm(Y[t], mean = 0, sd = sqrt(beta^2*exp(x.sample[,t])))
      norm.weights[,t] <-  weight[,t]/sum(weight[,t])
      
    }
    
    log.likelihood[j, i] <- sum(log(colSums(weight[, 2:T])) - log(N))
  }
}

pdf(paste("./Figures/H.3(a)Boxplot_Loglikelihood.pdf"), width = 12, height = 5)

colnames(log.likelihood) <- grid
boxplot(log.likelihood, xlab = expression(phi),  ylab = "Log-likelihood")

dev.off()

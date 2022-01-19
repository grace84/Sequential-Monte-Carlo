lambda <- seq(0, 2, 0.01)
r <-matrix(0, nrow = 1, ncol = length(lambda))

for (l in 1:length(lambda)){
  ## define the integrated function
  integrand <- function(x) {exp(-(x^2*(1-lambda[l]/2)))}
  
  ## integrate the function from 0 to infinity
  r[1,l] <- integrate(integrand, lower = 0, upper = Inf)$value
  
}
par(mfrow=c(1,1))
plot(lambda[1:(length(lambda)-1)], r[1,1:(length(lambda)-1)], xlab = "Lambda", ylab = "Evaluation of integral",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

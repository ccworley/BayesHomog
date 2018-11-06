#rm(list=ls())
require(rjags)
require(coda)
library(xtable)

set.seed(4444)
Beta0 <- 2
Beta1 <- 8
Beta2 <- -0.8
N <- 100
Sigma <- 5

x <- xbm  #runif(n=N, min=0, max=10) 
y <- ylumsrt  #Beta0 + Beta1 * x + Beta2 * x^2 + rnorm(n=N, mean=0, sd=Sigma)


Data <- list(x=x,  y=y)


lmfit <- lm(y~x+I(x^2))
summary(lmfit)
plot(x, y)
lines(x[order(x)], predict(lmfit)[order(x)], col="blue") # add fit line

jags.script <- "
model{
# likelihood
for( i in 1:length(x[])) {
y[i] ~ dnorm(mu[i], tau)
mu[i] <- beta0 + beta1 * x[i] + beta2 * pow(x[i], 2)
}

# priors
beta0 ~ dnorm(0.0, 1.0E-6)
beta1  ~ dnorm(0.0, 1.0E-6)
beta2  ~ dunif(-10000, -0.01) # exclude zero to permit turning point
tau  ~ dgamma(0.1,0.1)

# transformations
sigma <- pow(tau, 0.5)
delta <- (0 - beta1) / (2 * beta2)
}
" 

jags.fit <- jags.model(textConnection(jags.script), 
                       data=Data, n.chains=4, n.adapt=1000)

update(jags.fit, n.iter=1000) # burnin

jags.samples <- coda.samples(model=jags.fit,
                             variable.names=c('beta0', 'beta1', 
                                              'beta2', 'sigma', 'delta'),
                             n.iter=2000)
plot(jags.samples)
summary(jags.samples) 

beta0.posterior.mean <- summary(jags.samples)$statistics["beta0", "Mean"]
beta1.posterior.mean <- summary(jags.samples)$statistics["beta1", "Mean"]
beta2.posterior.mean <- summary(jags.samples)$statistics["beta2", "Mean"]

plot(x, y)
ypred <- beta0.posterior.mean + beta1.posterior.mean*x + 
  beta2.posterior.mean*x^2
lines(x[order(x)], ypred[order(x)], col="red")
lines(x[order(x)], predict(lmfit)[order(x)], col="blue") # add fit line
abline(v=summary(jags.samples)$statistics['delta', 'Mean']) # turning point
legend(x=min(x[order(x)]), y=max(y)*0.95, legend=c("posterior mean", "lm fit"), 
       col=c("red", "blue"), lty=1)

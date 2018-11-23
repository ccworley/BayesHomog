# R. Smiljanic
# (start, end) = (Apr 2018, Oct 2018)
# Gaia-ESO iDR6
# Homogenisation code for the atmospheric parameters
# Bayesian MCMC inference of the CovMatrix of the Nodes (i.e., sigmas and correlations of the atmospheric parameters) and the bias functions
# The CovMatrix and bias functions can then be applied on the "other stars" to estimate the best value of their atmospheric parameters
#################################
# Load the libraries needed to the MCMC Bayesian analysis and do some plots

library(rjags);library(runjags);library(coda);library(jagstools);library(lattice);library(ggplot2);library(ggthemes)#;library(plot3D)

# nodes.param and bench.param come from previous steps

# Oct 2018 - Try using only the Fe1 scatter as Node error

############################################
# Function to add a column with E_FEH in the bench.param matrix

adds.feh.error <- function(data.of.bench=bench.param,error.columns=c('SIG1_MH0')) {# From Fe1 scatter
  if (length(error.columns) == 0) {
    stop(paste('No choice of columns to adopt as E_FEH was made'))
  } else if (length(error.columns) == 1) {
    num.col <- which(colnames(data.of.bench) == error.columns)
    errors.metallicity <- round(data.of.bench[,num.col],digits=2)   
  } else if (length(error.columns) > 1) {
    squared.errors <- apply(data.of.bench[,(colnames(data.of.bench) %in% error.columns)],2,function(x) x^2)
    errors.metallicity <- apply(squared.errors,1,sum)
    errors.metallicity <- sqrt(errors.metallicity)
  }
  new.col.names <- c(colnames(data.of.bench),'E_FEH')
  data.of.bench <- cbind(data.of.bench,round(errors.metallicity,digits=2))         # This is the matrix that holds the parameters of the benchmarks
  colnames(data.of.bench) <- new.col.names
  
  return(data.of.bench)
}

############################################
# Function to create the matrix with SNR ratio of the individual benchmark spectra

create.snr.vector <- function(snr.data=metadata.of.bench.spectra$SNR) {
  snr.spec <- as.numeric(as.vector(snr.data))
  filter.high.snr <- (snr.spec > 900)
  snr.spec[filter.high.snr] <- trunc(runif(snr.spec[filter.high.snr],900,1000))
  snr.spec.vec <- matrix(rep(snr.spec,num.nodes),ncol=num.nodes,nrow=length(snr.spec))
  return(snr.spec.vec)
}

############################################
# Function to separate the measurement of the benchmark spectra from the general node results

node.measurements.of.reference.stars <- function(data.of.nodes=node.measured.param.bench) {
  
  observed.node.teff.spectrum <- matrix(NA,ncol=ncol(data.of.nodes),nrow=nrow(data.of.nodes))
  for (ik in seq(1,ncol(observed.node.teff.spectrum))) {
    observed.node.teff.spectrum[,ik] <- trunc(as.numeric(as.vector(data.of.nodes[,ik,'TEFF']))) # To be sure they are numbers
  }
  
  observed.node.logg.spectrum <- matrix(NA,ncol=ncol(data.of.nodes),nrow=nrow(data.of.nodes))
  for (ik in seq(1,ncol(observed.node.logg.spectrum))) {
    observed.node.logg.spectrum[,ik] <- round(as.numeric(as.vector(data.of.nodes[,ik,'LOGG'])),digits=2) # To be sure they are numbers
  }
  
  observed.node.feh.spectrum <- matrix(NA,ncol=ncol(data.of.nodes),nrow=nrow(data.of.nodes))
  for (ik in seq(1,ncol(observed.node.feh.spectrum))) {
    observed.node.feh.spectrum[,ik] <- round(as.numeric(as.vector(data.of.nodes[,ik,'FEH'])),digits=2) # To be sure they are numbers
  }
  
  all.data <- list(observed.node.teff.spectrum,observed.node.logg.spectrum,observed.node.feh.spectrum)
  return(all.data)
}


############################################
# Function to create some data structures for FEH, because some benchmarks do not have FEH measurmeent (despite having TEFF and LOGG)

correct.data.for.feh <- function(bench.data=bench.param,orig.star.code=star.code,orig.metadata=metadata.of.bench.spectra,observed.feh=observed.node.feh.spectrum,observed.snr=snr.spec.vec) {
  bench.without.feh <- which(is.na(bench.data$FEH))  
  bench.with.feh <- which(!is.na(bench.data$FEH))  
  positions.with.feh <- which(!(orig.star.code %in% bench.without.feh))  
  
  corr.bench.param <- bench.data[bench.with.feh,]
  corr.metadata.of.bench.spectra <- orig.metadata[positions.with.feh,]
  corr.observed.node.feh.spectrum <- observed.feh[positions.with.feh,]
  corr.observed.snr <- observed.snr[positions.with.feh,]
  
  corr.star.code <- vector('numeric',length=nrow(corr.observed.node.feh.spectrum))
  for (ik in seq(1,nrow(corr.observed.node.feh.spectrum))) {
    corr.star.code[ik] <- which(corr.bench.param$ID1 == corr.metadata.of.bench.spectra$CNAME[ik])
    #corr.star.code[ik] <- which(corr.bench.param$GES_FLD == corr.metadata.of.bench.spectra$GES_FLD[ik])
  }
  
  corr.given.feh.bench <- bench.data$FEH[bench.with.feh]
  corr.given.sigma.feh.bench <- bench.data$E_FEH[bench.with.feh]
  corr.given.sigma.feh.bench[is.na(corr.given.sigma.feh.bench)] <- 0.2 # If there is no official error, assume a large value
  
  all.data <- list(corr.bench.param,corr.metadata.of.bench.spectra,corr.observed.node.feh.spectrum,corr.star.code,corr.given.feh.bench,corr.given.sigma.feh.bench,positions.with.feh,corr.observed.snr)
  return(all.data)
}

############################################
# Function to create the "manipulated" version of TEFF and LOGG and FEH data

create.manipulated.teff.logg.feh <- function(observed.teff=observed.node.teff.spectrum,observed.logg=observed.node.logg.spectrum,observed.feh=observed.node.feh.spectrum,bench.data=bench.param,the.stars=star.code) {
  # What we are trying to do boils down to inferring "sigma" given that a node measured "x[i]" and knowing that the real values were "y[i]" for each star "i", or in R jags notation: x[i] ~ dnorm ( mu = y[i] , tau = 1/sigma^2 ); but it happens that sometimes x[i] is missing because the node failed to analyse that star. What I would like to do is, in these cases of missing measurement, simply assume a broad non informative prior for the missing x[i]: x[i] ~ dunif(3000,8000); in other words: ok, the node did not provide a measurement, but I know the value should be between 3000 and 8000 K. **BUT**, I can not write these two lines of code together: x[i] ~ dnorm ( mu = y[i] , tau = 1/sigma^2 ) and x[i] ~ dunif(3000,8000) because JAGS correctly thinks that I am trying to define the same "node" x[i] twice and throws an error. On the other hand, JAGS knows how to deal with the situation where it is "y[i]" that has the missing values. If we write x[i] ~ dnorm ( mu = y[i] , tau = 1/sigma^2 ) and y[i] ~ dunif(3000,9000) it automatically knows that it should treat **ONLY** the missing values of y[i] with the non-informative priors (really! I tested it. If y = c(1,2,NA,3,4), with the y[i] ~ dunif(3000,9000) it will substitute only the third element of y for the prior, the rest remains unchanged.
  # So, what is the numerical solution for this? I am proposing the following trick where the values to be randomized are inverted. Let's say that the measured values are x = c(1.1,1.9,NA,4.2,4.8,6.3) and the true values are y = c(1,2,3,4,5,6) and I want to infer the sigma as before in x[i] ~ dnorm ( mu = y[i] , tau = 1/sigma^2 ). I would like to treat the missing value with a prior x[i] ~ dunif (0,10), but as I said, this can not be done directly. So, what I do is, I manipulate x and y by taking out the NA from "x" and inputing there the corresponding true value; and removing the corresponding true value from y and inputing there the NA: new.x = c(1.1,1.9,3,4.2,4.8,6.3) and new.y = c(1,2,NA,4,5,6) telling JAGS that when it sees the NA in new.y it should use new.y ~ dunif(0,10). So, instead of randomizing the measured value, I am randomizing the true value. In essence, I think the two formulations are equivalent in what concerns the determination of the "sigma" that we are interested in. Before, I was trying to say that when the node measured y[3], it got a random value x[3] between 0 and 10. Now, what I am saying that it measured x[3] out of a random value of y[3] between 0 and 10.
  # This is the new matrix without missing values, where the NAs and NaNs where substituted by the reference benchmak values
  
  manipulated.node.teff.spectrum <- matrix(NA,ncol=ncol(observed.teff),nrow=nrow(observed.teff))
  manipulated.node.logg.spectrum <- matrix(NA,ncol=ncol(observed.logg),nrow=nrow(observed.logg))
  manipulated.node.feh.spectrum <- matrix(NA,ncol=ncol(observed.feh),nrow=nrow(observed.feh))
  
  # And I also need to save the position (i,j) in the matrix of each element that was changed from NA to a value, those that were not changed, and the star.code of the changed values. (Perhaps there is a better and less verbose way to do this, but this was the solution I could come up with)
  i.of.missing.teff <- vector()
  j.of.missing.teff <- vector()
  i.not.missing.teff <- vector()
  j.not.missing.teff <- vector()
  code.bench.missing <- vector()     # The star.code of the benchmark for which the measurement was missing
  code.bench.not.missing <- vector() # The star.code of the benchmark for which measurements were not missing
  keep.track.of.missing <- matrix(NA,nrow=nrow(observed.teff),ncol=ncol(observed.teff))
  
  ik <- 1
  jk <- 1
  for (i in seq(1,nrow(observed.teff))) {
    for (j in seq(1,ncol(observed.teff))) {
      if (is.na(observed.teff[i,j])) {
        manipulated.node.teff.spectrum[i,j] <- bench.data$TEFF[the.stars[i]]
        manipulated.node.logg.spectrum[i,j] <- bench.data$LOGG[the.stars[i]]
        manipulated.node.feh.spectrum[i,j] <- bench.data$FEH[the.stars[i]]
        i.of.missing.teff[ik] <- i
        j.of.missing.teff[ik] <- j
        code.bench.missing[ik] <- the.stars[i]
        keep.track.of.missing[i,j] <- 2
        ik <- ik+1
      } else {
        if (is.na(observed.logg[i,j])) {
          manipulated.node.teff.spectrum[i,j] <- bench.data$TEFF[the.stars[i]]
          manipulated.node.logg.spectrum[i,j] <- bench.data$LOGG[the.stars[i]]
          manipulated.node.feh.spectrum[i,j] <- bench.data$FEH[the.stars[i]]
          i.of.missing.teff[ik] <- i
          j.of.missing.teff[ik] <- j
          code.bench.missing[ik] <- the.stars[i]
          keep.track.of.missing[i,j] <- 2
          ik <- ik+1
        } else {      
          if (is.na(observed.feh[i,j])) {
            manipulated.node.teff.spectrum[i,j] <- bench.data$TEFF[the.stars[i]]
            manipulated.node.logg.spectrum[i,j] <- bench.data$LOGG[the.stars[i]]
            manipulated.node.feh.spectrum[i,j] <- bench.data$FEH[the.stars[i]]
            i.of.missing.teff[ik] <- i
            j.of.missing.teff[ik] <- j
            code.bench.missing[ik] <- the.stars[i]
            keep.track.of.missing[i,j] <- 2
            ik <- ik+1
          } else {      
            manipulated.node.teff.spectrum[i,j] <- observed.teff[i,j]
            manipulated.node.logg.spectrum[i,j] <- observed.logg[i,j]
            manipulated.node.feh.spectrum[i,j] <- observed.feh[i,j]
            i.not.missing.teff[jk] <- i
            j.not.missing.teff[jk] <- j
            code.bench.not.missing[jk] <- the.stars[i]
            keep.track.of.missing[i,j] <- 1
            jk <- jk+1
          }
        }
      }
    }
  }
  
  num.missing.values <- ik-1     # save the total number of missing entries
  num.not.missing.values <- jk-1 # save the total number of non-missing entries
  
  all.data <- list(manipulated.node.teff.spectrum,manipulated.node.logg.spectrum,manipulated.node.feh.spectrum,i.of.missing.teff,j.of.missing.teff,i.not.missing.teff,
                   j.not.missing.teff,code.bench.missing,code.bench.not.missing,keep.track.of.missing,num.missing.values,num.not.missing.values)
  return(all.data)
  
}


############################################
# Function to create the "manipulated" version of TEFF and LOGG data

create.manipulated.teff.logg <- function(observed.teff=observed.node.teff.spectrum,observed.logg=observed.node.logg.spectrum,bench.data=bench.param,the.stars=star.code) {
# What we are trying to do boils down to inferring "sigma" given that a node measured "x[i]" and knowing that the real values were "y[i]" for each star "i", or in R jags notation: x[i] ~ dnorm ( mu = y[i] , tau = 1/sigma^2 ); but it happens that sometimes x[i] is missing because the node failed to analyse that star. What I would like to do is, in these cases of missing measurement, simply assume a broad non informative prior for the missing x[i]: x[i] ~ dunif(3000,8000); in other words: ok, the node did not provide a measurement, but I know the value should be between 3000 and 8000 K. **BUT**, I can not write these two lines of code together: x[i] ~ dnorm ( mu = y[i] , tau = 1/sigma^2 ) and x[i] ~ dunif(3000,8000) because JAGS correctly thinks that I am trying to define the same "node" x[i] twice and throws an error. On the other hand, JAGS knows how to deal with the situation where it is "y[i]" that has the missing values. If we write x[i] ~ dnorm ( mu = y[i] , tau = 1/sigma^2 ) and y[i] ~ dunif(3000,9000) it automatically knows that it should treat **ONLY** the missing values of y[i] with the non-informative priors (really! I tested it. If y = c(1,2,NA,3,4), with the y[i] ~ dunif(3000,9000) it will substitute only the third element of y for the prior, the rest remains unchanged.
# So, what is the numerical solution for this? I am proposing the following trick where the values to be randomized are inverted. Let's say that the measured values are x = c(1.1,1.9,NA,4.2,4.8,6.3) and the true values are y = c(1,2,3,4,5,6) and I want to infer the sigma as before in x[i] ~ dnorm ( mu = y[i] , tau = 1/sigma^2 ). I would like to treat the missing value with a prior x[i] ~ dunif (0,10), but as I said, this can not be done directly. So, what I do is, I manipulate x and y by taking out the NA from "x" and inputing there the corresponding true value; and removing the corresponding true value from y and inputing there the NA: new.x = c(1.1,1.9,3,4.2,4.8,6.3) and new.y = c(1,2,NA,4,5,6) telling JAGS that when it sees the NA in new.y it should use new.y ~ dunif(0,10). So, instead of randomizing the measured value, I am randomizing the true value. In essence, I think the two formulations are equivalent in what concerns the determination of the "sigma" that we are interested in. Before, I was trying to say that when the node measured y[3], it got a random value x[3] between 0 and 10. Now, what I am saying that it measured x[3] out of a random value of y[3] between 0 and 10.
# This is the new matrix without missing values, where the NAs and NaNs where substituted by the reference benchmak values

manipulated.node.teff.spectrum <- matrix(NA,ncol=ncol(observed.teff),nrow=nrow(observed.teff))
manipulated.node.logg.spectrum <- matrix(NA,ncol=ncol(observed.logg),nrow=nrow(observed.logg))

# And I also need to save the position (i,j) in the matrix of each element that was changed from NA to a value, those that were not changed, and the star.code of the changed values. (Perhaps there is a better and less verbose way to do this, but this was the solution I could come up with)
i.of.missing.teff <- vector()
j.of.missing.teff <- vector()
i.not.missing.teff <- vector()
j.not.missing.teff <- vector()
code.bench.missing <- vector()     # The star.code of the benchmark for which the measurement was missing
code.bench.not.missing <- vector() # The star.code of the benchmark for which measurements were not missing
keep.track.of.missing <- matrix(NA,nrow=nrow(observed.teff),ncol=ncol(observed.teff))

ik <- 1
jk <- 1
for (i in seq(1,nrow(observed.teff))) {
  for (j in seq(1,ncol(observed.teff))) {
    if (is.na(observed.teff[i,j])) {
      manipulated.node.teff.spectrum[i,j] <- bench.data$TEFF[the.stars[i]]
      manipulated.node.logg.spectrum[i,j] <- bench.data$LOGG[the.stars[i]]
      i.of.missing.teff[ik] <- i
      j.of.missing.teff[ik] <- j
      code.bench.missing[ik] <- the.stars[i]
      keep.track.of.missing[i,j] <- 2
      ik <- ik+1
    } else {
      if (is.na(observed.logg[i,j])) {
        manipulated.node.teff.spectrum[i,j] <- bench.data$TEFF[the.stars[i]]
        manipulated.node.logg.spectrum[i,j] <- bench.data$LOGG[the.stars[i]]
        i.of.missing.teff[ik] <- i
        j.of.missing.teff[ik] <- j
        code.bench.missing[ik] <- the.stars[i]
        keep.track.of.missing[i,j] <- 2
        ik <- ik+1
      } else {      
        manipulated.node.teff.spectrum[i,j] <- observed.teff[i,j]
        manipulated.node.logg.spectrum[i,j] <- observed.logg[i,j]
        i.not.missing.teff[jk] <- i
        j.not.missing.teff[jk] <- j
        code.bench.not.missing[jk] <- the.stars[i]
        keep.track.of.missing[i,j] <- 1
        jk <- jk+1
      }
    }
  }
}

num.missing.values <- ik-1     # save the total number of missing entries
num.not.missing.values <- jk-1 # save the total number of non-missing entries

all.data <- list(manipulated.node.teff.spectrum,manipulated.node.logg.spectrum,i.of.missing.teff,j.of.missing.teff,i.not.missing.teff,
                 j.not.missing.teff,code.bench.missing,code.bench.not.missing,keep.track.of.missing,num.missing.values,num.not.missing.values)
return(all.data)

}

############################################
# Function to create the manipulated version of FEH data

create.manipulated.feh <- function(observed.feh=observed.node.feh.spectrum,bench.data=corr.bench.param,the.stars=corr.star.code,positions=positions.with.feh) {

# For FEH is different because some of the benchmarks do not have [Fe/H]
i.of.missing.feh <- vector()
j.of.missing.feh <- vector()
i.not.missing.feh <- vector()
j.not.missing.feh <- vector()
corr.code.bench.missing <- vector()     # The corr.star.code of the benchmark for which the measurement was missing
corr.code.bench.not.missing <- vector() # The corr.star.code of the benchmark for which measurements were not missing
manipulated.node.feh.spectrum <- matrix(NA,ncol=ncol(observed.feh),nrow=nrow(observed.feh))
corr.manipulated.node.feh.spectrum <- manipulated.node.feh.spectrum[positions,]
corr.keep.track.of.missing <- matrix(NA,nrow=nrow(corr.observed.node.feh.spectrum),ncol=ncol(corr.observed.node.feh.spectrum))
ik <- 1
jk <- 1
for (i in seq(1,nrow(corr.observed.node.feh.spectrum))) {
  for (j in seq(1,ncol(corr.observed.node.feh.spectrum))) {
    if (is.na(corr.observed.node.feh.spectrum[i,j])) {
      corr.manipulated.node.feh.spectrum[i,j] <- bench.data$FEH[the.stars[i]]
      i.of.missing.feh[ik] <- i
      j.of.missing.feh[ik] <- j
      corr.code.bench.missing[ik] <- the.stars[i]
      corr.keep.track.of.missing[i,j] <- 2
      ik <- ik+1
    } else {
      corr.manipulated.node.feh.spectrum[i,j] <- corr.observed.node.feh.spectrum[i,j]
      i.not.missing.feh[jk] <- i
      j.not.missing.feh[jk] <- j
      corr.code.bench.not.missing[jk] <- the.stars[i]
      corr.keep.track.of.missing[i,j] <- 1
      jk <- jk+1
    }
  }
}
corr.num.missing.values <- ik-1    # save the total number of missing entries
corr.num.not.missing.values <- jk-1 # save the total number of non-missing entries

all.data <- list(corr.manipulated.node.feh.spectrum,i.of.missing.feh,j.of.missing.feh,i.not.missing.feh,j.not.missing.feh,
                 corr.code.bench.missing,corr.code.bench.not.missing,corr.keep.track.of.missing,corr.num.missing.values,corr.num.not.missing.values)
return(all.data)
}

############################################
# Function to create the priors for the covariance matrix

create.prior.covariance <- function(n.nodes=num.nodes,typical.sigma) {
  # 5) The Inverse Covariance Matrix
  #  We will model the problem using a multi variate Gaussian where we know the vector of measured values of the node, know the vector with true.teff of the benchmark and want to infer the "Inverse Covariance Matrix". Actually we want the Covariance Matrix, but in the JAGS notation we write the problem with the Inverse Covariance Matrix:
  #  observed.node.teff[j,1:K] ~ dmnorm( true.teff.benchmark[j,1:K] , teff.InvCovMat[1:K,1:K] )
  
  # The Inverse Covariance Matrix needs a prior.
  # From what I read in the internet, the recommended prior on the inverse covariance matrix is a generic Wishart distribution, which is a multidimensional generalization of a gamma distribution. W( V, n ) takes as input the degrees of freedom "n" and V = Sigma_0 * n^-1 where Sigma_0 is a first guess of the covariance matrix. The non-informative choice for the degrees of freedom is to set it equal to the number of rows (in this case the number of nodes).
  
  #prior.deg.freedom = num.nodes
  prior.matrix = (1/n.nodes)*diag(typical.sigma^2,ncol=n.nodes,nrow=n.nodes)  # Usually assuming that the typical sigma of Teff = 100 K, sigma of logg = 0.2 dex, sigma of FEH = 0.2 dex
  
  # Some links:
  #http://doingbayesiandataanalysis.blogspot.com/2017/06/bayesian-estimation-of-correlations-and.html
  #http://www.themattsimpson.com/2012/08/20/prior-distributions-for-covariance-matrices-the-scaled-inverse-wishart-prior/
  #https://stats.stackexchange.com/questions/66338/how-to-specify-the-wishart-distribution-scale-matrix
  #https://en.wikipedia.org/wiki/Wishart_distribution
  return(prior.matrix)
}

############################################
# Function to create the vector of setups of each spectrum of the benchmarks

create.setup.vector <- function(column.setups=metadata.of.bench.spectra$SETUP) {
#6) We are going to fit one bias function per setup (i.e., U520, U580 - this can eventually be extended for Giraffe, but maybe not easilly. 
#   In WG11 all nodes analyze all spectra of all setups. That is not the case in WG10. This would break the symmetry of the data structures per setup). 
#   We create a vector with a code for each setup:

setups.used <- levels(as.factor((column.setups)))    # Finds out all the setups. Below it converts each setup "name" to a nummeric code
vector.of.setups <- as.character(as.vector(column.setups))
ik <- 0
for (each.setup in setups.used) {
  ik <- ik+1
  vector.of.setups[(vector.of.setups == each.setup)] <- as.character(ik)
}
vector.of.setups <- as.numeric(as.vector(vector.of.setups)) 
num.setups <- length(setups.used) # Total number of setups

all.data <- list(vector.of.setups,num.setups)
return(all.data)
}

############################################
# And now let's write the model itself

# TEFF MODEL FOR BENCHMARKS

model.teff.matrix.bias <- "model{
# We set up the model to use multivariate normal distributions. For a given spectrum of the benchmarks, the vector of Teff measured by the nodes teff.nodes.vector = (Teff1, Teff2,...,Teffn) is a random drawing from a multivariate normal distribution with a mean = vector of true.teff of that benchmark (true.teff, true.teff,...,true.teff) - where obviously the true.teff is the same irrespective of node - and there is a covariance matrix (JAGS expects the inverse covariance matrix actually): teff.nodes.vector ~ dmnorm ( mu = true.teff.vector, InvCovMatrix ). So we need to first set up all vectors accordingly.

# 1)
# This first part is to define the vector with TRUE values of the benchmark Teff
# We do not use the given Teff value directly, but use the given value as a prior for the true value

# Prior on the true Teff of the benchmarks
for (i in 1:N) { # Running over N which is the num.bench

    the.true.teff.bench[i] ~ dnorm(given.teff.benchmarks[i], given.tau.teff.bench[i])  # The given benchmark Teff value is a prior on the true value of the star's Teff
    given.tau.teff.bench[i] <- pow(given.error.teff.benchmarks[i],-2)                  # dnorm uses the precision, tau, which is the inverse of the variance: tau = 1 / (sigma^2). So, from here we also believe on the error of the Teff value

}

# Now we will create the vector of true Teffs: vector.true.teff = (true.teff, ..., true.teff)
# Here we have to remember that we use directly true.teff only when the node provided a measurement. If the value was missing, we actually use a broad non-informative prior to simulate that the value of the node measurement was unknown. This means that the actual vector.true.teffs is not one per benchmark, but one per spectrum! Because for each spectrum there is a different number of node values missing.

# true.teff.benchmark.per.spec is a matrix with nrow = number of spectrum and ncol = number of nodes, equivalent to the input node.teff matrix

for (i in 1:n.not.missing) { # When the node gave a value (n.not.missing), we include the true.teff of the benchmark estimated above
    true.teff.benchmark.per.spec[i.not.missing[i],j.not.missing[i]] <- the.true.teff.bench[code.bench.not.missing[i]]
}

#Prior on the missing values

for (l in 1:n.missing) { # When the node was missing a value (n.missing), we include a value drawn from the broad prior

    true.teff.benchmark.per.spec[i.of.missing.values[l],j.of.missing.values[l]] ~ dunif(3000,8000) #dnorm(offset.teff,scale.teff) I(3000,8000)
}

# Prior on the teff.InvCovMat

teff.InvCovMat ~ dwish( prior.matrix[1:K,1:K] , prior.deg.freedom )

# Other Priors needed for coefficients
# Priors on the coefficients alpha for the fit of the bias function. We fit one bias function per node. The bias.teff is for now a quadratic function of the true.teff
#for (i in 1:2) { # For missing and non-missing values
    for (l in 1:K) {       # On number of nodes
       for (m in 1:P) {    # On number of setups
          for (j in 1:ncoeff) { # 3 coefficients needed for fitting a quadratic function to 1 parameter
              alpha[j,l,m] ~ dnorm(0,pow(0.1,-2)) # Remember here that dnorm(mu,tau) where mu = mean and tau = 1/sigma^2, and sigma is the standard deviation
           }
        }
    }
#}

# I will eventually need also a vector that really has only the the.true.teff.bench[i], without taking into account the missing values, for the fitting of the bias function (see below)

for (i in 1:M) {  #M=num.spectra, K=num.nodes
   for (j in 1:K) {
        only.the.true.teff.benchmark[i,j] <- (given.teff.benchmarks[star.code[i]]-offset.teff)/scale.teff
       #only.the.true.teff.benchmark[i,j] <- the.true.teff.bench[star.code[i]]

   }
}

# 2)
# Now, we finally write that what the nodes give is a noisy measurement of the true.teff.benchmark

#Likelihood - per spectrum (vectorized over K = num.nodes)

for (j in 1:M) {    #M=num.spectra, K=num.nodes
#Fitting the true value
  observed.node.teff.per.spec[j,1:K] ~ dmnorm( (true.teff.benchmark.per.spec[j,1:K] + (norm.bias.vector.teff.node[j,1:K]*scale.teff)) , teff.InvCovMat[1:K,1:K] )

  #Based on TEFF only -> ncoeff = 3 coefficients    
  norm.bias.vector.teff.node[j,1:K] <- alpha[1,1:K,setup.code[j]] + alpha[2,1:K,setup.code[j]]*only.the.true.teff.benchmark[j,1:K] + alpha[3,1:K,setup.code[j]]*pow(only.the.true.teff.benchmark[j,1:K],2) 

}

# Notice what was done with the bias function above. Let's assume, for example, that I know that everytime I measure the Teff of the Sun, my value comes out with a bias of +50 K. I could write that my measurement minus the bias was drawn from the normal distribution: (my.measurement - 50 K) ~ dnorm(mu = 5777 K, sigma = my.random.error). What I am assuming is that this is equivalent to say that I am instead drawing from a normal distribution with mu = 5827 K (5777 + 50): my.measurement ~ dnorm(mu = 5827 K + my.bias , sigma = my.random.error)

# 3)
# This is extra. Since the quantities I know how to understand are the sigmas and correlations out of the Covariance Matrix:
# Convert teff.invCovMat to node.sd and correlations:

  teff.CovMat <- inverse( teff.InvCovMat )
  for ( varIdx in 1:K ) {
      node.sd.teff[varIdx] <- sqrt(teff.CovMat[varIdx,varIdx])
   }
  for ( varIdx1 in 1:K ) {
       for ( varIdx2 in 1:K ) {
           teff.Rho[varIdx1,varIdx2] <- ( teff.CovMat[varIdx1,varIdx2] / (node.sd.teff[varIdx1]*node.sd.teff[varIdx2]) )
    }
  }

}"

########################################
# TEFF MODEL with LOGG FOR BENCHMARKS

model.teff.logg.matrix.bias <- "model{
# We set up the model to use multivariate normal distributions. For a given spectrum of the benchmarks, the vector of Teff measured by the nodes teff.nodes.vector = (Teff1, Teff2,...,Teffn) is a random drawing from a multivariate normal distribution with a mean = vector of true.teff of that benchmark (true.teff, true.teff,...,true.teff) - where obviously the true.teff is the same irrespective of node - and there is a covariance matrix (JAGS expects the inverse covariance matrix actually): teff.nodes.vector ~ dmnorm ( mu = true.teff.vector, InvCovMatrix ). So we need to first set up all vectors accordingly.

# 1)
# This first part is to define the vector with TRUE values of the benchmark Teff
# We do not use the given Teff value directly, but use the given value as a prior for the true value

# Prior on the true Teff of the benchmarks
for (i in 1:N) { # Running over N which is the num.bench

the.true.teff.bench[i] ~ dnorm(given.teff.benchmarks[i], given.tau.teff.bench[i])  # The given benchmark Teff value is a prior on the true value of the star's Teff
given.tau.teff.bench[i] <- pow(given.error.teff.benchmarks[i],-2)                  # dnorm uses the precision, tau, which is the inverse of the variance: tau = 1 / (sigma^2). So, from here we also believe on the error of the Teff value

the.true.logg.bench[i] ~ dnorm(given.logg.benchmarks[i], given.tau.logg.bench[i])  # The given benchmark Teff value is a prior on the true value of the star's Teff
given.tau.logg.bench[i] <- pow(given.error.logg.benchmarks[i],-2)                  # dnorm uses the precision, tau, which is the inverse of the variance: tau = 1 / (sigma^2). So, from here we also believe on the error of the Teff value

#the.true.feh.bench[i] ~ dnorm(given.feh.benchmarks[i], given.tau.feh.bench[i])  # The given benchmark Teff value is a prior on the true value of the star's Teff
#given.tau.feh.bench[i] <- pow(given.error.feh.benchmarks[i],-2)                  # dnorm uses the precision, tau, which is the inverse of the variance: tau = 1 / (sigma^2). So, from here we also believe on the error of the Teff value

}

# Now we will create the vector of true Teffs: vector.true.teff = (true.teff, ..., true.teff)
# Here we have to remember that we use directly true.teff only when the node provided a measurement. If the value was missing, we actually use a broad non-informative prior to simulate that the value of the node measurement was unknown. This means that the actual vector.true.teffs is not one per benchmark, but one per spectrum! Because for each spectrum there is a different number of node values missing.

# true.teff.benchmark.per.spec is a matrix with nrow = number of spectrum and ncol = number of nodes, equivalent to the input node.teff matrix

for (i in 1:n.not.missing) { # When the node gave a value (n.not.missing), we include the true.teff of the benchmark estimated above

true.teff.benchmark.per.spec[i.not.missing[i],j.not.missing[i]] <- the.true.teff.bench[code.bench.not.missing[i]]
true.logg.benchmark.per.spec[i.not.missing[i],j.not.missing[i]] <- the.true.logg.bench[code.bench.not.missing[i]]
#true.feh.benchmark.per.spec[i.not.missing[i],j.not.missing[i]] <- the.true.feh.bench[code.bench.not.missing[i]]

}

#Prior on the missing values

for (l in 1:n.missing) { # When the node was missing a value (n.missing), we include a value drawn from the broad prior

true.teff.benchmark.per.spec[i.of.missing.values[l],j.of.missing.values[l]] ~ dunif(3000,8000) #dnorm(offset.teff,scale.teff) I(3000,8000)
true.logg.benchmark.per.spec[i.of.missing.values[l],j.of.missing.values[l]] ~ dunif(0,6)   #dnorm(5300,pow(800,-2))   #dunif(-1,1)  #dunif(3000,8000)   #dnorm(0,1) #when scaled by meanTbm and sigTbm  #
#true.feh.benchmark.per.spec[i.of.missing.values[l],j.of.missing.values[l]] ~ dunif(-3.5,0.5)   #dnorm(5300,pow(800,-2))   #dunif(-1,1)  #dunif(3000,8000)   #dnorm(0,1) #when scaled by meanTbm and sigTbm  #

}

# Prior on the teff.InvCovMat

teff.InvCovMat ~ dwish( prior.matrix[1:K,1:K] , prior.deg.freedom )

# Other Priors needed for coefficients
# Priors on the coefficients alpha for the fit of the bias function. We fit one bias function per node. The bias.teff is for now a quadratic function of the true.teff
#for (i in 1:2) { # For missing and non-missing values
for (l in 1:K) {       # On number of nodes
for (m in 1:P) {    # On number of setups
for (j in 1:ncoeff) { # 3 coefficients needed for fitting a quadratic function to 1 parameter
alpha[j,l,m] ~ dnorm(0,pow(0.1,-2)) # Remember here that dnorm(mu,tau) where mu = mean and tau = 1/sigma^2, and sigma is the standard deviation
}
}
}
#}

# I will eventually need also a vector that really has only the the.true.teff.bench[i], without taking into account the missing values, for the fitting of the bias function (see below)

for (i in 1:M) {  #M=num.spectra, K=num.nodes
for (j in 1:K) {
only.the.true.teff.benchmark[i,j] <- (given.teff.benchmarks[star.code[i]]-offset.teff)/scale.teff
#only.the.true.teff.benchmark[i,j] <- the.true.teff.bench[star.code[i]]

only.the.true.logg.benchmark[i,j] <- (given.logg.benchmarks[star.code[i]]-offset.logg)/scale.logg
#only.the.true.feh.benchmark[i,j] <- (given.feh.benchmarks[star.code[i]]-offset.feh)/scale.feh
}
}

# 2)
# Now, we finally write that what the nodes give is a noisy measurement of the true.teff.benchmark

#Likelihood - per spectrum (vectorized over K = num.nodes)

for (j in 1:M) {    #M=num.spectra, K=num.nodes
#Fitting the true value
observed.node.teff.per.spec[j,1:K] ~ dmnorm( (true.teff.benchmark.per.spec[j,1:K] + (norm.bias.vector.teff.node[j,1:K]*scale.teff)) , teff.InvCovMat[1:K,1:K] )

#Based on TEFF only -> ncoeff = 3 coefficients    
#norm.bias.vector.teff.node[j,1:K] <- alpha[1,1:K,setup.code[j]] + alpha[2,1:K,setup.code[j]]*only.the.true.teff.benchmark[j,1:K] + alpha[3,1:K,setup.code[j]]*pow(only.the.true.teff.benchmark[j,1:K],2) 

#Based on TEFF & LOGG only -> ncoeff = 6 coefficients
norm.bias.vector.teff.node[j,1:K] <- alpha[1,1:K,setup.code[j]] + alpha[2,1:K,setup.code[j]]*only.the.true.teff.benchmark[j,1:K] + alpha[3,1:K,setup.code[j]]*pow(only.the.true.teff.benchmark[j,1:K],2) + alpha[4,1:K,setup.code[j]]*only.the.true.logg.benchmark[j,1:K] + alpha[5,1:K,setup.code[j]]*pow(only.the.true.logg.benchmark[j,1:K],2) + alpha[6,1:K,setup.code[j]]*only.the.true.teff.benchmark[j,1:K]*only.the.true.logg.benchmark[j,1:K]

}

# Notice what was done with the bias function above. Let's assume, for example, that I know that everytime I measure the Teff of the Sun, my value comes out with a bias of +50 K. I could write that my measurement minus the bias was drawn from the normal distribution: (my.measurement - 50 K) ~ dnorm(mu = 5777 K, sigma = my.random.error). What I am assuming is that this is equivalent to say that I am instead drawing from a normal distribution with mu = 5827 K (5777 + 50): my.measurement ~ dnorm(mu = 5827 K + my.bias , sigma = my.random.error)

# 3)
# This is extra. Since the quantities I know how to understand are the sigmas and correlations out of the Covariance Matrix:
# Convert teff.invCovMat to node.sd and correlations:

teff.CovMat <- inverse( teff.InvCovMat )
for ( varIdx in 1:K ) {
node.sd.teff[varIdx] <- sqrt(teff.CovMat[varIdx,varIdx])
}
for ( varIdx1 in 1:K ) {
for ( varIdx2 in 1:K ) {
teff.Rho[varIdx1,varIdx2] <- ( teff.CovMat[varIdx1,varIdx2] / (node.sd.teff[varIdx1]*node.sd.teff[varIdx2]) )
}
}

}"


########################################
# TEFF with FEH MODEL FOR BENCHMARKS

model.teff.feh.matrix.bias <- "model{
# We set up the model to use multivariate normal distributions. For a given spectrum of the benchmarks, the vector of Teff measured by the nodes teff.nodes.vector = (Teff1, Teff2,...,Teffn) is a random drawing from a multivariate normal distribution with a mean = vector of true.teff of that benchmark (true.teff, true.teff,...,true.teff) - where obviously the true.teff is the same irrespective of node - and there is a covariance matrix (JAGS expects the inverse covariance matrix actually): teff.nodes.vector ~ dmnorm ( mu = true.teff.vector, InvCovMatrix ). So we need to first set up all vectors accordingly.

# 1)
# This first part is to define the vector with TRUE values of the benchmark Teff
# We do not use the given Teff value directly, but use the given value as a prior for the true value

# Prior on the true Teff of the benchmarks
for (i in 1:N) { # Running over N which is the num.bench

  the.true.teff.bench[i] ~ dnorm(given.teff.benchmarks[i], given.tau.teff.bench[i])  # The given benchmark Teff value is a prior on the true value of the star's Teff
  given.tau.teff.bench[i] <- pow(given.error.teff.benchmarks[i],-2)                  # dnorm uses the precision, tau, which is the inverse of the variance: tau = 1 / (sigma^2). So, from here we also believe on the error of the Teff value
  
  #the.true.logg.bench[i] ~ dnorm(given.logg.benchmarks[i], given.tau.logg.bench[i])  # The given benchmark Teff value is a prior on the true value of the star's Teff
  #given.tau.logg.bench[i] <- pow(given.error.logg.benchmarks[i],-2)                  # dnorm uses the precision, tau, which is the inverse of the variance: tau = 1 / (sigma^2). So, from here we also believe on the error of the Teff value
  
  the.true.feh.bench[i] ~ dnorm(given.feh.benchmarks[i], given.tau.feh.bench[i])  # The given benchmark Teff value is a prior on the true value of the star's Teff
  given.tau.feh.bench[i] <- pow(given.error.feh.benchmarks[i],-2)                  # dnorm uses the precision, tau, which is the inverse of the variance: tau = 1 / (sigma^2). So, from here we also believe on the error of the Teff value

}

# Now we will create the vector of true Teffs: vector.true.teff = (true.teff, ..., true.teff)
# Here we have to remember that we use directly true.teff only when the node provided a measurement. If the value was missing, we actually use a broad non-informative prior to simulate that the value of the node measurement was unknown. This means that the actual vector.true.teffs is not one per benchmark, but one per spectrum! Because for each spectrum there is a different number of node values missing.

# true.teff.benchmark.per.spec is a matrix with nrow = number of spectrum and ncol = number of nodes, equivalent to the input node.teff matrix

for (i in 1:n.not.missing) { # When the node gave a value (n.not.missing), we include the true.teff of the benchmark estimated above

  true.teff.benchmark.per.spec[i.not.missing[i],j.not.missing[i]] <- the.true.teff.bench[code.bench.not.missing[i]]
  #true.logg.benchmark.per.spec[i.not.missing[i],j.not.missing[i]] <- the.true.logg.bench[code.bench.not.missing[i]]
  true.feh.benchmark.per.spec[i.not.missing[i],j.not.missing[i]] <- the.true.feh.bench[code.bench.not.missing[i]]

}

#Prior on the missing values

for (l in 1:n.missing) { # When the node was missing a value (n.missing), we include a value drawn from the broad prior

  true.teff.benchmark.per.spec[i.of.missing.values[l],j.of.missing.values[l]] ~ dunif(3000,8000) #dnorm(offset.teff,scale.teff) I(3000,8000)
  #true.logg.benchmark.per.spec[i.of.missing.values[l],j.of.missing.values[l]] ~ dunif(0,6)   #dnorm(5300,pow(800,-2))   #dunif(-1,1)  #dunif(3000,8000)   #dnorm(0,1) #when scaled by meanTbm and sigTbm  #
  true.feh.benchmark.per.spec[i.of.missing.values[l],j.of.missing.values[l]] ~ dunif(-3.5,0.5)   #dnorm(5300,pow(800,-2))   #dunif(-1,1)  #dunif(3000,8000)   #dnorm(0,1) #when scaled by meanTbm and sigTbm  #

}

# Prior on the teff.InvCovMat

teff.InvCovMat ~ dwish( prior.matrix[1:K,1:K] , prior.deg.freedom )

# Other Priors needed for coefficients
# Priors on the coefficients alpha for the fit of the bias function. We fit one bias function per node. The bias.teff is for now a quadratic function of the true.teff
#for (i in 1:2) { # For missing and non-missing values
for (l in 1:K) {       # On number of nodes
  for (m in 1:P) {    # On number of setups
    for (j in 1:ncoeff) { # 3 coefficients needed for fitting a quadratic function to 1 parameter
      alpha[j,l,m] ~ dnorm(0,pow(0.1,-2)) # Remember here that dnorm(mu,tau) where mu = mean and tau = 1/sigma^2, and sigma is the standard deviation
    }
  }
}
#}

# I will eventually need also a vector that really has only the the.true.teff.bench[i], without taking into account the missing values, for the fitting of the bias function (see below)

for (i in 1:M) {
  for (j in 1:K) {
    only.the.true.teff.benchmark[i,j] <- (given.teff.benchmarks[star.code[i]]-offset.teff)/scale.teff
    #only.the.true.teff.benchmark[i,j] <- the.true.teff.bench[star.code[i]]
    
    #only.the.true.logg.benchmark[i,j] <- (given.logg.benchmarks[star.code[i]]-offset.logg)/scale.logg
    only.the.true.feh.benchmark[i,j] <- (given.feh.benchmarks[star.code[i]]-offset.feh)/scale.feh
  }
}

# 2)
# Now, we finally write that what the nodes give is a noisy measurement of the true.teff.benchmark

#Likelihood - per spectrum (vectorized over K = num.nodes)

for (j in 1:M) {
#Fitting the true value
observed.node.teff.per.spec[j,1:K] ~ dmnorm( (true.teff.benchmark.per.spec[j,1:K] + (norm.bias.vector.teff.node[j,1:K]*scale.teff)) , teff.InvCovMat[1:K,1:K] )

#Based on TEFF only -> ncoeff = 3 coefficients    
#norm.bias.vector.teff.node[j,1:K] <- alpha[1,1:K,setup.code[j]] + alpha[2,1:K,setup.code[j]]*only.the.true.teff.benchmark[j,1:K] + alpha[3,1:K,setup.code[j]]*pow(only.the.true.teff.benchmark[j,1:K],2) 

#Based on TEFF & LOGG only -> ncoeff = 6 coefficients
#norm.bias.vector.teff.node[j,1:K] <- alpha[1,1:K,setup.code[j]] + alpha[2,1:K,setup.code[j]]*only.the.true.teff.benchmark[j,1:K] + alpha[3,1:K,setup.code[j]]*pow(only.the.true.teff.benchmark[j,1:K],2) + alpha[4,1:K,setup.code[j]]*only.the.true.logg.benchmark[j,1:K] + alpha[5,1:K,setup.code[j]]*pow(only.the.true.logg.benchmark[j,1:K],2) + alpha[6,1:K,setup.code[j]]*only.the.true.teff.benchmark[j,1:K]*only.the.true.logg.benchmark[j,1:K]

#Based on TEFF & FEH only -> ncoeff = 6 coefficients
norm.bias.vector.teff.node[j,1:K] <- alpha[1,1:K,setup.code[j]] + alpha[2,1:K,setup.code[j]]*only.the.true.teff.benchmark[j,1:K] + alpha[3,1:K,setup.code[j]]*pow(only.the.true.teff.benchmark[j,1:K],2) + alpha[4,1:K,setup.code[j]]*only.the.true.feh.benchmark[j,1:K] + alpha[5,1:K,setup.code[j]]*pow(only.the.true.feh.benchmark[j,1:K],2) + alpha[6,1:K,setup.code[j]]*only.the.true.teff.benchmark[j,1:K]*only.the.true.feh.benchmark[j,1:K]

}

# Notice what was done with the bias function above. Let's assume, for example, that I know that everytime I measure the Teff of the Sun, my value comes out with a bias of +50 K. I could write that my measurement minus the bias was drawn from the normal distribution: (my.measurement - 50 K) ~ dnorm(mu = 5777 K, sigma = my.random.error). What I am assuming is that this is equivalent to say that I am instead drawing from a normal distribution with mu = 5827 K (5777 + 50): my.measurement ~ dnorm(mu = 5827 K + my.bias , sigma = my.random.error)

# 3)
# This is extra. Since the quantities I know how to understand are the sigmas and correlations out of the Covariance Matrix:
# Convert teff.invCovMat to node.sd and correlations:

teff.CovMat <- inverse( teff.InvCovMat )
for ( varIdx in 1:K ) {
node.sd.teff[varIdx] <- sqrt(teff.CovMat[varIdx,varIdx])
}
for ( varIdx1 in 1:K ) {
for ( varIdx2 in 1:K ) {
teff.Rho[varIdx1,varIdx2] <- ( teff.CovMat[varIdx1,varIdx2] / (node.sd.teff[varIdx1]*node.sd.teff[varIdx2]) )
}
}

}"


############################################
# LOGG MODEL FOR BENCHMARKS

model.logg.matrix.bias <- "model{

# Prior on the true logg of the benchmarks
for (i in 1:N) { # Running over N which is the num.bench

the.true.logg.bench[i] ~ dnorm(given.logg.benchmarks[i], given.tau.logg.bench[i]) 
given.tau.logg.bench[i] <- pow(given.error.logg.benchmarks[i],-2)                 

}

for (i in 1:n.not.missing) { 

true.logg.benchmark.per.spec[i.not.missing[i],j.not.missing[i]] <- the.true.logg.bench[code.bench.not.missing[i]]

}

#Prior on the missing values

for (l in 1:n.missing) { 

true.logg.benchmark.per.spec[i.of.missing.values[l],j.of.missing.values[l]] ~ dunif(0.0,5.0) 

}

# Prior on the logg.InvCovMat

logg.InvCovMat ~ dwish( prior.matrix[1:K,1:K] , prior.deg.freedom )

# Other Priors needed for coefficients

for (l in 1:K) {       # On number of nodes
    for (m in 1:P) {    # On number of setups
        for (j in 1:ncoeff) { # 3 coefficients needed for fitting a quadratic function
            alpha[j,l,m] ~ dnorm(0,pow(0.1,-2)) # Remember here that dnorm(mu,tau) where mu = mean and tau = 1/sigma^2, and sigma is the standard deviation
        }
    }
}


# I will eventually need also a vector that really has only the the.true.logg.bench[i], without taking into account the missing values, for the fitting of the bias function (see below)

for (i in 1:M) {
    for (j in 1:K) {
        only.the.true.logg.benchmark[i,j] <- (given.logg.benchmarks[star.code[i]]-offset.logg)/scale.logg
    }
}

# 2)

#Likelihood - per spectrum (vectorized over K = num.nodes)

for (j in 1:M) {
#Fitting the true value
observed.node.logg.per.spec[j,1:K] ~ dmnorm( (true.logg.benchmark.per.spec[j,1:K] + (norm.bias.vector.logg.node[j,1:K]*scale.logg)) , logg.InvCovMat[1:K,1:K] )

#Bias function as a quadratic function of logg  -> ncoeff = 3 coefficients    
norm.bias.vector.logg.node[j,1:K] <- alpha[1,1:K,setup.code[j]] + alpha[2,1:K,setup.code[j]]*only.the.true.logg.benchmark[j,1:K] + alpha[3,1:K,setup.code[j]]*pow(only.the.true.logg.benchmark[j,1:K],2) 

  #Based on TEFF & LOGG only -> ncoeff = 6 coefficients
  #norm.bias.vector.logg.node[j,1:K] <- alpha[1,1:K,setup.code[j]] + alpha[2,1:K,setup.code[j]]*only.the.true.logg.benchmark[j,1:K] + alpha[3,1:K,setup.code[j]]*pow(only.the.true.logg.benchmark[j,1:K],2) + alpha[4,1:K,setup.code[j]]*only.the.true.teff.benchmark[j,1:K] + alpha[5,1:K,setup.code[j]]*pow(only.the.true.teff.benchmark[j,1:K],2) + alpha[6,1:K,setup.code[j]]*only.the.true.logg.benchmark[j,1:K]*only.the.true.teff.benchmark[j,1:K]

}


# 3)
# This is extra. Since the quantities I know how to understand are the sigmas and correlations out of the Covariance Matrix:
# Convert logg.invCovMat to node.sd and correlations:

logg.CovMat <- inverse( logg.InvCovMat )
for ( varIdx in 1:K ) {
node.sd.logg[varIdx] <- sqrt(logg.CovMat[varIdx,varIdx])
}
for ( varIdx1 in 1:K ) {
for ( varIdx2 in 1:K ) {
logg.Rho[varIdx1,varIdx2] <- ( logg.CovMat[varIdx1,varIdx2] / (node.sd.logg[varIdx1]*node.sd.logg[varIdx2]) )
}
}

}"

########################################
# LOGG with FEH MODEL FOR BENCHMARKS

model.logg.feh.matrix.bias <- "model{
# We set up the model to use multivariate normal distributions. For a given spectrum of the benchmarks, the vector of logg measured by the nodes logg.nodes.vector = (logg1, logg2,...,loggn) is a random drawing from a multivariate normal distribution with a mean = vector of true.logg of that benchmark (true.logg, true.teff,...,true.teff) - where obviously the true.teff is the same irrespective of node - and there is a covariance matrix (JAGS expects the inverse covariance matrix actually): teff.nodes.vector ~ dmnorm ( mu = true.teff.vector, InvCovMatrix ). So we need to first set up all vectors accordingly.

# 1)
# This first part is to define the vector with TRUE values of the benchmark logg
# We do not use the given logg value directly, but use the given value as a prior for the true value

# Prior on the true logg of the benchmarks
for (i in 1:N) { # Running over N which is the num.bench

  the.true.logg.bench[i] ~ dnorm(given.logg.benchmarks[i], given.tau.logg.bench[i])  # The given benchmark logg value is a prior on the true value of the star's Teff
  given.tau.logg.bench[i] <- pow(given.error.logg.benchmarks[i],-2)                  # dnorm uses the precision, tau, which is the inverse of the variance: tau = 1 / (sigma^2). So, from here we also believe on the error of the Teff value
  
  the.true.feh.bench[i] ~ dnorm(given.feh.benchmarks[i], given.tau.feh.bench[i])  # The given benchmark feh value is a prior on the true value of the star's Teff
  given.tau.feh.bench[i] <- pow(given.error.feh.benchmarks[i],-2)                  # dnorm uses the precision, tau, which is the inverse of the variance: tau = 1 / (sigma^2). So, from here we also believe on the error of the Teff value

}

# Now we will create the vector of true loggs: vector.true.teff = (true.logg, ..., true.logg)
# Here we have to remember that we use directly true.logg only when the node provided a measurement. If the value was missing, we actually use a broad non-informative prior to simulate that the value of the node measurement was unknown. This means that the actual vector.true.teffs is not one per benchmark, but one per spectrum! Because for each spectrum there is a different number of node values missing.

# true.logg.benchmark.per.spec is a matrix with nrow = number of spectrum and ncol = number of nodes, equivalent to the input node.teff matrix

for (i in 1:n.not.missing) { # When the node gave a value (n.not.missing), we include the true.teff of the benchmark estimated above

  true.logg.benchmark.per.spec[i.not.missing[i],j.not.missing[i]] <- the.true.logg.bench[code.bench.not.missing[i]]
  true.feh.benchmark.per.spec[i.not.missing[i],j.not.missing[i]] <- the.true.feh.bench[code.bench.not.missing[i]]

}

#Prior on the missing values

for (l in 1:n.missing) { # When the node was missing a value (n.missing), we include a value drawn from the broad prior

  true.logg.benchmark.per.spec[i.of.missing.values[l],j.of.missing.values[l]] ~ dunif(0,6)   #dnorm(5300,pow(800,-2))   #dunif(-1,1)  #dunif(3000,8000)   #dnorm(0,1) #when scaled by meanTbm and sigTbm  #
  true.feh.benchmark.per.spec[i.of.missing.values[l],j.of.missing.values[l]] ~ dunif(-3.5,0.5)   #dnorm(5300,pow(800,-2))   #dunif(-1,1)  #dunif(3000,8000)   #dnorm(0,1) #when scaled by meanTbm and sigTbm  #

}

# Prior on the teff.InvCovMat

logg.InvCovMat ~ dwish( prior.matrix[1:K,1:K] , prior.deg.freedom )

# Other Priors needed for coefficients
# Priors on the coefficients alpha for the fit of the bias function. We fit one bias function per node. The bias.teff is for now a quadratic function of the true.teff
#for (i in 1:2) { # For missing and non-missing values
for (l in 1:K) {       # On number of nodes
  for (m in 1:P) {    # On number of setups
    for (j in 1:ncoeff) { # 3 coefficients needed for fitting a quadratic function to 1 parameter
      alpha[j,l,m] ~ dnorm(0,pow(0.1,-2)) # Remember here that dnorm(mu,tau) where mu = mean and tau = 1/sigma^2, and sigma is the standard deviation
    }
  }
}
#}

# I will eventually need also a vector that really has only the the.true.teff.bench[i], without taking into account the missing values, for the fitting of the bias function (see below)

for (i in 1:M) {
  for (j in 1:K) {
    only.the.true.logg.benchmark[i,j] <- (given.logg.benchmarks[star.code[i]]-offset.logg)/scale.logg
    only.the.true.feh.benchmark[i,j] <- (given.feh.benchmarks[star.code[i]]-offset.feh)/scale.feh
  }
}

# 2)
# Now, we finally write that what the nodes give is a noisy measurement of the true.teff.benchmark

#Likelihood - per spectrum (vectorized over K = num.nodes)

for (j in 1:M) {
  #Fitting the true value
  observed.node.logg.per.spec[j,1:K] ~ dmnorm( (true.logg.benchmark.per.spec[j,1:K] + (norm.bias.vector.logg.node[j,1:K]*scale.logg)) , logg.InvCovMat[1:K,1:K] )
  
  #Based on logg only -> ncoeff = 3 coefficients    
  #norm.bias.vector.logg.node[j,1:K] <- alpha[1,1:K,setup.code[j]] + alpha[2,1:K,setup.code[j]]*only.the.true.logg.benchmark[j,1:K] + alpha[3,1:K,setup.code[j]]*pow(only.the.true.logg.benchmark[j,1:K],2) 
  
  #Based on logg & FEH only -> ncoeff = 6 coefficients
  norm.bias.vector.logg.node[j,1:K] <- alpha[1,1:K,setup.code[j]] + alpha[2,1:K,setup.code[j]]*only.the.true.logg.benchmark[j,1:K] + alpha[3,1:K,setup.code[j]]*pow(only.the.true.logg.benchmark[j,1:K],2) + alpha[4,1:K,setup.code[j]]*only.the.true.feh.benchmark[j,1:K] + alpha[5,1:K,setup.code[j]]*pow(only.the.true.feh.benchmark[j,1:K],2) + alpha[6,1:K,setup.code[j]]*only.the.true.logg.benchmark[j,1:K]*only.the.true.feh.benchmark[j,1:K]

}

# Notice what was done with the bias function above. Let's assume, for example, that I know that everytime I measure the Teff of the Sun, my value comes out with a bias of +50 K. I could write that my measurement minus the bias was drawn from the normal distribution: (my.measurement - 50 K) ~ dnorm(mu = 5777 K, sigma = my.random.error). What I am assuming is that this is equivalent to say that I am instead drawing from a normal distribution with mu = 5827 K (5777 + 50): my.measurement ~ dnorm(mu = 5827 K + my.bias , sigma = my.random.error)

# 3)
# This is extra. Since the quantities I know how to understand are the sigmas and correlations out of the Covariance Matrix:
# Convert logg.invCovMat to node.sd and correlations:

logg.CovMat <- inverse( logg.InvCovMat )
for ( varIdx in 1:K ) {
  node.sd.logg[varIdx] <- sqrt(logg.CovMat[varIdx,varIdx])
}
for ( varIdx1 in 1:K ) {
  for ( varIdx2 in 1:K ) {
    logg.Rho[varIdx1,varIdx2] <- ( logg.CovMat[varIdx1,varIdx2] / (node.sd.logg[varIdx1]*node.sd.logg[varIdx2]) )
  }
}

}"



############################################
# FEH MODEL FOR BENCHMARKS

model.feh.matrix.bias <- "model{

# Prior on the true feh of the benchmarks
for (i in 1:N) { # Running over N which is the num.bench

the.true.feh.bench[i] ~ dnorm(given.feh.benchmarks[i], given.tau.feh.bench[i]) 
given.tau.feh.bench[i] <- pow(given.error.feh.benchmarks[i],-2)                 

}

for (i in 1:n.not.missing) { 

true.feh.benchmark.per.spec[i.not.missing[i],j.not.missing[i]] <- the.true.feh.bench[code.bench.not.missing[i]]

}

#Prior on the missing values

for (l in 1:n.missing) { 

true.feh.benchmark.per.spec[i.of.missing.values[l],j.of.missing.values[l]] ~ dunif(-3.0,0.5) 

}

# Prior on the feh.InvCovMat

feh.InvCovMat ~ dwish( prior.matrix[1:K,1:K] , prior.deg.freedom )

# Other Priors needed for coefficients

for (l in 1:K) {       # On number of nodes
for (m in 1:P) {    # On number of setups
for (j in 1:ncoeff) { # 3 coefficients needed for fitting a quadratic function
alpha[j,l,m] ~ dnorm(0,pow(0.1,-2)) # Remember here that dnorm(mu,tau) where mu = mean and tau = 1/sigma^2, and sigma is the standard deviation
}
}
}


# I will eventually need also a vector that really has only the the.true.feh.bench[i], without taking into account the missing values, for the fitting of the bias function (see below)

for (i in 1:M) {
for (j in 1:K) {
only.the.true.feh.benchmark[i,j] <- (given.feh.benchmarks[star.code[i]]-offset.feh)/scale.feh
}
}

# 2)

#Likelihood - per spectrum (vectorized over K = num.nodes)

for (j in 1:M) {
#Fitting the true value
observed.node.feh.per.spec[j,1:K] ~ dmnorm( (true.feh.benchmark.per.spec[j,1:K] + (norm.bias.vector.feh.node[j,1:K]*scale.feh)) , feh.InvCovMat[1:K,1:K] )

#Bias function as a quadratic function of feh    
norm.bias.vector.feh.node[j,1:K] <- alpha[1,1:K,setup.code[j]] + alpha[2,1:K,setup.code[j]]*only.the.true.feh.benchmark[j,1:K] + alpha[3,1:K,setup.code[j]]*pow(only.the.true.feh.benchmark[j,1:K],2) 

 #Based on FEH & LOGG only -> ncoeff = 6 coefficients
  #norm.bias.vector.feh.node[j,1:K] <- alpha[1,1:K,setup.code[j]] + alpha[2,1:K,setup.code[j]]*only.the.true.feh.benchmark[j,1:K] + alpha[3,1:K,setup.code[j]]*pow(only.the.true.feh.benchmark[j,1:K],2) + alpha[4,1:K,setup.code[j]]*only.the.true.logg.benchmark[j,1:K] + alpha[5,1:K,setup.code[j]]*pow(only.the.true.logg.benchmark[j,1:K],2) + alpha[6,1:K,setup.code[j]]*only.the.true.feh.benchmark[j,1:K]*only.the.true.logg.benchmark[j,1:K]

}


# 3)
# This is extra. Since the quantities I know how to understand are the sigmas and correlations out of the Covariance Matrix:
# Convert feh.invCovMat to node.sd and correlations:

feh.CovMat <- inverse( feh.InvCovMat )
for ( varIdx in 1:K ) {
node.sd.feh[varIdx] <- sqrt(feh.CovMat[varIdx,varIdx])
}
for ( varIdx1 in 1:K ) {
for ( varIdx2 in 1:K ) {
feh.Rho[varIdx1,varIdx2] <- ( feh.CovMat[varIdx1,varIdx2] / (node.sd.feh[varIdx1]*node.sd.feh[varIdx2]) )
}
}

}"

########################################
# FEH with LOGG MODEL FOR BENCHMARKS

model.feh.logg.matrix.bias <- "model{
# We set up the model to use multivariate normal distributions. For a given spectrum of the benchmarks, the vector of logg measured by the nodes feh.nodes.vector = (feh1, feh2,...,fehn) is a random drawing from a multivariate normal distribution with a mean = vector of true.logg of that benchmark (true.logg, true.teff,...,true.teff) - where obviously the true.teff is the same irrespective of node - and there is a covariance matrix (JAGS expects the inverse covariance matrix actually): teff.nodes.vector ~ dmnorm ( mu = true.teff.vector, InvCovMatrix ). So we need to first set up all vectors accordingly.

# 1)
# This first part is to define the vector with TRUE values of the benchmark feh
# We do not use the given feh value directly, but use the given value as a prior for the true value

# Prior on the true feh of the benchmarks
for (i in 1:N) { # Running over N which is the num.bench

  the.true.feh.bench[i] ~ dnorm(given.feh.benchmarks[i], given.tau.feh.bench[i])  # The given benchmark feh value is a prior on the true value of the star's Teff
  given.tau.feh.bench[i] <- pow(given.error.feh.benchmarks[i],-2)                  # dnorm uses the precision, tau, which is the inverse of the variance: tau = 1 / (sigma^2). So, from here we also believe on the error of the Teff value
  
  the.true.logg.bench[i] ~ dnorm(given.logg.benchmarks[i], given.tau.logg.bench[i])  # The given benchmark logg value is a prior on the true value of the star's Teff
  given.tau.logg.bench[i] <- pow(given.error.logg.benchmarks[i],-2)                  # dnorm uses the precision, tau, which is the inverse of the variance: tau = 1 / (sigma^2). So, from here we also believe on the error of the Teff value

}

# Now we will create the vector of true loggs: vector.true.teff = (true.logg, ..., true.logg)
# Here we have to remember that we use directly true.logg only when the node provided a measurement. If the value was missing, we actually use a broad non-informative prior to simulate that the value of the node measurement was unknown. This means that the actual vector.true.teffs is not one per benchmark, but one per spectrum! Because for each spectrum there is a different number of node values missing.

# true.logg.benchmark.per.spec is a matrix with nrow = number of spectrum and ncol = number of nodes, equivalent to the input node.teff matrix

for (i in 1:n.not.missing) { # When the node gave a value (n.not.missing), we include the true.teff of the benchmark estimated above

  true.feh.benchmark.per.spec[i.not.missing[i],j.not.missing[i]] <- the.true.feh.bench[code.bench.not.missing[i]]
  true.logg.benchmark.per.spec[i.not.missing[i],j.not.missing[i]] <- the.true.logg.bench[code.bench.not.missing[i]]

}

#Prior on the missing values

for (l in 1:n.missing) { # When the node was missing a value (n.missing), we include a value drawn from the broad prior

  true.feh.benchmark.per.spec[i.of.missing.values[l],j.of.missing.values[l]] ~ dunif(-3.5,0.5)   #dnorm(5300,pow(800,-2))   #dunif(-1,1)  #dunif(3000,8000)   #dnorm(0,1) #when scaled by meanTbm and sigTbm  #
  true.logg.benchmark.per.spec[i.of.missing.values[l],j.of.missing.values[l]] ~ dunif(0,6)   #dnorm(5300,pow(800,-2))   #dunif(-1,1)  #dunif(3000,8000)   #dnorm(0,1) #when scaled by meanTbm and sigTbm  #

}

# Prior on the teff.InvCovMat

feh.InvCovMat ~ dwish( prior.matrix[1:K,1:K] , prior.deg.freedom )

# Other Priors needed for coefficients
# Priors on the coefficients alpha for the fit of the bias function. We fit one bias function per node. The bias.teff is for now a quadratic function of the true.teff
#for (i in 1:2) { # For missing and non-missing values
for (l in 1:K) {       # On number of nodes
  for (m in 1:P) {    # On number of setups
    for (j in 1:ncoeff) { # 3 coefficients needed for fitting a quadratic function to 1 parameter
      alpha[j,l,m] ~ dnorm(0,pow(0.1,-2)) # Remember here that dnorm(mu,tau) where mu = mean and tau = 1/sigma^2, and sigma is the standard deviation
    }
  }
}
#}

# I will eventually need also a vector that really has only the the.true.teff.bench[i], without taking into account the missing values, for the fitting of the bias function (see below)

for (i in 1:M) {
  for (j in 1:K) {
    only.the.true.feh.benchmark[i,j] <- (given.feh.benchmarks[star.code[i]]-offset.feh)/scale.feh
    only.the.true.logg.benchmark[i,j] <- (given.logg.benchmarks[star.code[i]]-offset.logg)/scale.logg
  }
}

# 2)
# Now, we finally write that what the nodes give is a noisy measurement of the true.teff.benchmark

#Likelihood - per spectrum (vectorized over K = num.nodes)

for (j in 1:M) {
  #Fitting the true value
  observed.node.feh.per.spec[j,1:K] ~ dmnorm( (true.feh.benchmark.per.spec[j,1:K] + (norm.bias.vector.feh.node[j,1:K]*scale.feh)) , feh.InvCovMat[1:K,1:K] )
  
  #Based on logg only -> ncoeff = 3 coefficients    
  #norm.bias.vector.feh.node[j,1:K] <- alpha[1,1:K,setup.code[j]] + alpha[2,1:K,setup.code[j]]*only.the.true.feh.benchmark[j,1:K] + alpha[3,1:K,setup.code[j]]*pow(only.the.true.feh.benchmark[j,1:K],2) 
  
  #Based on logg & FEH only -> ncoeff = 6 coefficients
  norm.bias.vector.feh.node[j,1:K] <- alpha[1,1:K,setup.code[j]] + alpha[2,1:K,setup.code[j]]*only.the.true.feh.benchmark[j,1:K] + alpha[3,1:K,setup.code[j]]*pow(only.the.true.feh.benchmark[j,1:K],2) + alpha[4,1:K,setup.code[j]]*only.the.true.logg.benchmark[j,1:K] + alpha[5,1:K,setup.code[j]]*pow(only.the.true.logg.benchmark[j,1:K],2) + alpha[6,1:K,setup.code[j]]*only.the.true.feh.benchmark[j,1:K]*only.the.true.logg.benchmark[j,1:K]

}

# Notice what was done with the bias function above. Let's assume, for example, that I know that everytime I measure the Teff of the Sun, my value comes out with a bias of +50 K. I could write that my measurement minus the bias was drawn from the normal distribution: (my.measurement - 50 K) ~ dnorm(mu = 5777 K, sigma = my.random.error). What I am assuming is that this is equivalent to say that I am instead drawing from a normal distribution with mu = 5827 K (5777 + 50): my.measurement ~ dnorm(mu = 5827 K + my.bias , sigma = my.random.error)

# 3)
# This is extra. Since the quantities I know how to understand are the sigmas and correlations out of the Covariance Matrix:
# Convert logg.invCovMat to node.sd and correlations:

feh.CovMat <- inverse( feh.InvCovMat )
for ( varIdx in 1:K ) {
  node.sd.feh[varIdx] <- sqrt(feh.CovMat[varIdx,varIdx])
}
for ( varIdx1 in 1:K ) {
  for ( varIdx2 in 1:K ) {
    feh.Rho[varIdx1,varIdx2] <- ( feh.CovMat[varIdx1,varIdx2] / (node.sd.feh[varIdx1]*node.sd.feh[varIdx2]) )
  }
}

}"





############################################
# Initial values for free parameters

create.inits.ncoeff <- function(variable=c('TEFF'),num.chains=num.chains.for.teff,n.setups=num.setups,n.nodes=num.nodes,n.coeff=ncoeff) {
  variable <- str_trim(variable)
  if (length(variable) != 1) { stop(paste('Choose one variable at a time from TEFF or LOGG or FEH')) }
  if (!(variable %in% c("TEFF","LOGG","FEH"))) { stop(paste('Variable has to be one of TEFF, LOGG or FEH - and written in uppercase'))}
  
  list.of.names <- c('teff.InvCovMat','logg.InvCovMat','feh.InvCovMat')
  name.invCovMat <- list.of.names[which(c("TEFF","LOGG","FEH") == variable)]
  
  list.of.inits <- list()
  for (ik in 1:num.chains) {
    list.of.inits[[ik]] <-     list(alpha = array(rnorm(n.coeff*n.setups*n.nodes,0,0.1),dim=c(n.coeff,n.nodes,n.setups)),
                                    InvCovMat = diag(rgamma(n.nodes,0.1,0.1)),
                                    .RNG.name="base::Super-Duper", 
                                    .RNG.seed=((2*ik)+1))
    names(list.of.inits[[ik]]) <- c('alpha',name.invCovMat,'.RNG.name','.RNG.seed')
  }
  return(list.of.inits)
}


#inits.teff <- function(n.setups=num.setups,n.nodes=num.nodes) {
#    list(alpha = array(rnorm(3*n.setups*n.nodes,0,0.1),dim=c(3,n.nodes,n.setups)),
#     teff.InvCovMat = diag(rgamma(n.nodes,0.1,0.1))
#)
#}

#inits.logg <- function(n.setups=num.setups,n.nodes=num.nodes) {
#  list(alpha = array(rnorm(3*n.setups*n.nodes,0,0.1),dim=c(3,n.nodes,n.setups)),
#       logg.InvCovMat = diag(rgamma(n.nodes,0.1,0.1))
#  )
#}

#inits.feh <- function(n.setups=num.setups,n.nodes=num.nodes) {
#  list(alpha = array(rnorm(3*n.setups*n.nodes,0,0.1),dim=c(3,n.nodes,n.setups)),
#       feh.InvCovMat = diag(rgamma(n.nodes,0.1,0.1))
#  )
#}


############################################
# MAIN FUNCTION TO RUN THE MCMC SIMULATION FOR THE BENCHMARK DATA

run.the.reference.model <- function(input.data,the.model,num.chains,num.adapt,num.bur,size.sample,thin.factor,variables.to.monitor,model.inits) {
  
  ptm <- proc.time() # To print out how long it took to run (in seconds)
  
  simple_test <- run.jags(method="parallel", # Modified to run chains in parallel
                          data=input.data,
                          model=the.model,
                          n.chains=num.chains,
                          adapt=num.adapt,
                          monitor=variables.to.monitor,
                          inits=model.inits,
                          burnin=num.bur,
                          sample=size.sample,
                          summarise=F,
                          thin=thin.factor,
                          plots=F
  )
  
  # Automatically print all summaries
  for (ik in 1:length(variables.to.monitor)) {
     print(paste0('Summary of variables = ',variables.to.monitor[[ik]]))
     print(paste(" ",sep=" "))
     print(summary(simple_test,vars=variables.to.monitor[[ik]]))
     print(paste(" ",sep=" "))
  }
  
  
  print(proc.time() - ptm)
  gc()
  
  return(simple_test)
}

#simple_test <- run.jags(method="rjags",
#                        data=teff.data,
#                        mode=model.teff.matrix.bias,
#                        n.chains=3,
#                        adapt=1000,
#                        monitor=c("the.true.teff.bench","alpha","teff.Rho","node.sd.teff"),#"beta1","beta2","beta3"),#"node.sd.teff","teff.Rho","alpha1","alpha2","alpha3","new.tau"),#,"alpha4","alpha5","alpha6"),#,"beta1","beta2","beta3","gamma1","gamma2","gamma3","gamma4"),
#                        inits=list(inits1,inits2,inits3),
#                        burnin=1000,
#                        sample=2000,
#                        summarise=F,
#                        thin=100,
#                        plots=F
#)


#simple_test_logg <- run.jags(method="rjags",
#                             data=logg.data,
#                             mode=model.logg.matrix.bias,
#                             n.chains=3,
#                             adapt=1000,
#                             monitor=c("the.true.logg.bench","alpha","logg.Rho","node.sd.logg"),#"beta1","beta2","beta3"),#"node.sd.teff","teff.Rho","alpha1","alpha2","alpha3","new.tau"),#,"alpha4","alpha5","alpha6"),#,"beta1","beta2","beta3","gamma1","gamma2","gamma3","gamma4"),
#                             inits=list(inits1.logg,inits2.logg,inits3.logg),
#                             burnin=1000,
#                             sample=2000,
#                             summarise=F,
#                             thin=50,
#                             plots=F
#)


#simple_test_feh <- run.jags(method="rjags",
#                            data=feh.data,
#                            mode=model.feh.matrix.bias,
#                            n.chains=7,
#                            adapt=1000,
#                            monitor=c("the.true.feh.bench","alpha","feh.Rho","node.sd.feh"),#"beta1","beta2","beta3"),#"node.sd.teff","teff.Rho","alpha1","alpha2","alpha3","new.tau"),#,"alpha4","alpha5","alpha6"),#,"beta1","beta2","beta3","gamma1","gamma2","gamma3","gamma4"),
#                            inits=list(inits1.feh,inits2.feh,inits3.feh,inits4.feh,inits5.feh,inits6.feh,inits7.feh),
#                            burnin=1000,
#                            sample=1000,
#                            summarise=F,
#                            thin=200,
#                            plots=F
#)



############################################
# Function to print the correlation matrix

show.correlation.matrix <- function(the.model,variable=c('TEFF'),n.nodes=num.nodes,name.nodes=list.nodes) {
  variable <- str_trim(variable)
  if (length(variable) != 1) { stop(paste('Choose one variable at a time from TEFF or LOGG or FEH')) }
  if (!(variable %in% c("TEFF","LOGG","FEH"))) { stop(paste('Variable has to be one of TEFF, LOGG or FEH - and written in uppercase'))}
  
  list.of.names <- c('teff.Rho','logg.Rho','feh.Rho')
  name.rho <- list.of.names[which(c("TEFF","LOGG","FEH") == variable)]
  

  print(paste0('This is the final correlation matrix for ',variable,':'))
  print(' ')
  estimated.cor.matrix <- summary(the.model,vars=c(name.rho))[,4]
  estimated.cor.matrix <- matrix(round(estimated.cor.matrix,digits=2),n.nodes,n.nodes)
  colnames(estimated.cor.matrix) <- name.nodes
  rownames(estimated.cor.matrix) <- name.nodes

  print(estimated.cor.matrix)
}

############################################
# Function to plot all density.plots, in a series of 6 each time

make.density.plots <- function(the.model) {
  fit_mcmc <- as.mcmc(the.model)
  num.cols <- ncol(fit_mcmc)
  k <- 1
  print(paste('num.cols = ',num.cols))
  for (i in 1:trunc(num.cols/6)) {
    j <- 6*i
    if (j > num.cols) { j <- num.cols }
    print(paste0('k = ',k,' and j = ',j))
    #plot(densityplot(fit_mcmc[,k:j]))
    k <- j+1
  }
}

############################################
# Function to print the autocorrelation disgnostics for the alpha coefficients

corr.diag.of.variable <- function(the.model,variable=c('alpha')) {
  fit_mcmc <- as.mcmc(the.model)
  the.cols <- grepl(variable,colnames(fit_mcmc))
  print(autocorr.diag(fit_mcmc[,the.cols]))
}

############################################
# Function to plot the homogenised parameters against the reference parameters

plot.against.reference.bench <- function(the.model,variable=c('TEFF'),bench.data=bench.param,observed.snr=mean.bench.snr) {

  variable <- str_trim(variable)
  if (length(variable) != 1) { stop(paste('Choose one variable at a time from TEFF or LOGG or FEH')) }
  if (!(variable %in% c("TEFF","LOGG","FEH"))) { stop(paste('Variable has to be one of TEFF, LOGG or FEH - and written in uppercase'))}

  test <- summary(the.model,vars=c('the.true'))
  the.col <- which(colnames(bench.data) == variable)
  par(mfrow=c(2,2))
  
  plot(bench.data$TEFF,(bench.data[,the.col]-test[,4]),main=variable,xlab=paste0('Given TEFF of Reference Stars'),ylab=paste0('Delta ',variable,' (given minus homog)'))
  
  plot(bench.data$LOGG,(bench.data[,the.col]-test[,4]),main=variable,xlab='Given LOGG of Reference Stars',ylab=paste0('Delta ',variable,' (given minus homog)'))
  
  plot(bench.data$FEH,(bench.data[,the.col]-test[,4]),main=variable,xlab='Given FEH of Reference Stars',ylab=paste0('Delta ',variable,' (given minus homog)'))

  plot(observed.snr,(bench.data[,the.col]-test[,4]),main=variable,xlab='Observed SNR',ylab=paste0('Delta ',variable,' (given minus homog)'))
    par(mfrow=c(1,1))
  
}


############################################
# Function to plot the node biases

look.at.node.bias.functions <- function(the.model,variable=c('TEFF'),bench.data=bench.param,observed.data=observed.node.teff.spectrum,
                                        the.setups=vector.of.setups,the.stars=star.code,col.of.setups=metadata.of.bench.spectra$SETUP,
                                        nodes=list.nodes,mean.param.bench,sd.param.bench,observed.snr=snr.spec.vec) {
  variable <- str_trim(variable)
  if (length(variable) != 1) { stop(paste('Choose one variable at a time from TEFF or LOGG or FEH')) }
  if (!(variable %in% c("TEFF","LOGG","FEH"))) { stop(paste('Variable has to be one of TEFF, LOGG or FEH - and written in uppercase'))}
  
  this.bench.col <- which(colnames(bench.data) == variable)
  list.of.setups <- levels(as.factor(col.of.setups))
  
  coef.alphas <- summary(the.model,vars=c('alpha'))[,4]
  coef.alphas.1 <- summary(the.model,vars=c('alpha'))[,1] 
  coef.alphas.3 <- summary(the.model,vars=c('alpha'))[,3] 
  
  for (each.node in nodes) {
    num.of.node <- which(nodes == each.node)
    for (each.setup in list.of.setups) {
      num.setup <- which(list.of.setups == each.setup)
      print(each.node)
      print(each.setup)
      plot.x <- bench.data[the.stars[the.setups == num.setup],this.bench.col]
      plot.y <- (observed.data[(the.setups == num.setup),num.of.node]-plot.x)
      plot.s <- observed.snr[(the.setups == num.setup),num.of.node]
      str(plot.s)
      str(plot.y)
      plot.t <- bench.data[the.stars[the.setups == num.setup],'TEFF']
      plot.l <- bench.data[the.stars[the.setups == num.setup],'LOGG']
      plot.f <- bench.data[the.stars[the.setups == num.setup],'FEH']
      
      
      if (sum(is.na(plot.y)) == length(plot.y)) { 
        plot.x <- c(-5,5)
        plot.y <- plot.x
      }
      
      
      par(mfrow=c(2,2))
      
      #X
      #plot(plot.x,plot.y,pch = 16,col = rgb(0,0,0,0.5),
      #     main=paste0(each.node,' - ',each.setup),xlab=paste0('Given ',variable,' of Reference Stars (per spectrum)'),
      #     ylab=paste0('Delta ',variable,' Observed - Given'))
      
      x <- plot.x
      order.x <- order(x)
      x <- x[order.x]
      nor.x <- (x-mean.param.bench)/sd.param.bench

      #print(plot.x[order.x])
      #print(plot.y[order.x])
      
            
      k1 <- ((3*num.of.node)-2)+((num.setup-1)*(length(nodes)*3))
      k2 <- ((3*num.of.node)-1)+((num.setup-1)*(length(nodes)*3))
      k3 <- ((3*num.of.node))+((num.setup-1)*(length(nodes)*3))      
      
      #print(paste('node = ',each.node,' setup = ',each.setup,' k1, k2, k3 =',k1,k2,k3,sep=" "))
      
      nor.y <- coef.alphas[k1]+coef.alphas[k2]*nor.x+coef.alphas[k3]*nor.x^2
      y <- nor.y*sd.param.bench
      #lines(x,y,col='red',lwd=3)
      
      y.1 <- sd.param.bench*(coef.alphas.1[k1]+coef.alphas.1[k2]*nor.x+coef.alphas.1[k3]*nor.x^2)
      #lines(x,y.1,col='blue',lwd=3) 
      
      y.3 <- sd.param.bench*(coef.alphas.3[k1]+coef.alphas.3[k2]*nor.x+coef.alphas.3[k3]*nor.x^2)
      #lines(x,y.3,col='blue',lwd=3)

      #TEFF
      print('teff')
      #str(plot.t)
      #print(length(plot.t))
      #print(length(plot.y))
      plot.tx = plot.t[order.x]  #Puts LOGG in the same order as y (i.e. as x, the independent variable the bias is calculated on)
      order.t <- order(plot.tx)  #Puts LOGG of order x, into LOGG order.
      plot(plot.t,plot.y,pch = 16,col = rgb(0,0,0,0.25),
           main=paste0(each.node,' - ',each.setup),xlab=paste0('Given TEFF of Benchmarks (per spectrum)'),
           ylab=paste0('Delta ',variable,' Observed - Bench'))
      lines(plot.tx[order.t],y[order.t],col='red',lwd=3)
      lines(plot.tx[order.t],y.1[order.t],col='blue',lwd=3)
      lines(plot.tx[order.t],y.3[order.t],col='blue',lwd=3)
      points(plot.t,plot.y,pch = 16,col = rgb(0,0,0,0.25))


      #LOGG
      print('logg')
      plot.lx = plot.l[order.x]  #Puts LOGG in the same order as y (i.e. as x, the independent variable the bias is calculated on)
      order.l <- order(plot.lx)  #Puts LOGG of order x, into LOGG order.
      plot(plot.l,plot.y,pch = 16,col = rgb(0,0,0,0.25),
           main=paste0(each.node,' - ',each.setup),xlab=paste0('Given LOGG of Benchmarks (per spectrum)'),
           ylab=paste0('Delta ',variable,' Observed - Bench'))
      lines(plot.lx[order.l],y[order.l],col='red',lwd=3)
      lines(plot.lx[order.l],y.1[order.l],col='blue',lwd=3)
      lines(plot.lx[order.l],y.3[order.l],col='blue',lwd=3)
      points(plot.l,plot.y,pch = 16,col = rgb(0,0,0,0.25))

      #FEH
      print('feh')
      plot.fx = plot.f[order.x]
      order.f <- order(plot.fx)
      plot(plot.f,plot.y,pch = 16,col = rgb(0,0,0,0.5),
           main=paste0(each.node,' - ',each.setup),xlab=paste0('Given FEH of Benchmarks (per spectrum)'),
           ylab=paste0('Delta ',variable,' Observed - Bench'))

      lines(plot.fx[order.f],y[order.f],col='red',lwd=3)
      lines(plot.fx[order.f],y.1[order.f],col='blue',lwd=3)
      lines(plot.fx[order.f],y.3[order.f],col='blue',lwd=3)
      
      #SNR
      plot(plot.s,plot.y,pch = 16,col = rgb(0,0,0,0.5),
           main=paste0(each.node,' - ',each.setup),xlab=paste0('Measured SNR (per spectrum)'),
           ylab=paste0('Delta ',variable,' Observed - Bench'))
      
      #plot.tx = plot.t[order.x]
      #order.t <- order(plot.tx)
      #lines(plot.tx[order.t],y[order.t],col='red',lwd=3)
      #lines(plot.tx[order.t],y.1[order.t],col='blue',lwd=3) 
      #lines(plot.tx[order.t],y.3[order.t],col='blue',lwd=3)
      
      par(mfrow=c(1,1))
      
    }
  }
  
}

############################################
# Function to plot the node biases

look.at.node.bias.multifunctions <- function(the.model,variable=c('TEFF'),variable2=c('LOGG'),bench.data=bench.param,observed.data=observed.node.teff.spectrum,
                                        observed.data2=observed.node.logg.spectrum,n.coeff=ncoeff,     
                                        the.setups=vector.of.setups,the.stars=star.code,col.of.setups=metadata.of.bench.spectra$SETUP,
                                        nodes=list.nodes,mean.param.bench,sd.param.bench,mean.param2.bench,sd.param2.bench,
                                        observed.snr=snr.spec.vec) {
  variable <- str_trim(variable)
  if (length(variable) != 1) { stop(paste('Choose one variable at a time from TEFF or LOGG or FEH')) }
  if (!(variable %in% c("TEFF","LOGG","FEH"))) { stop(paste('Variable has to be one of TEFF, LOGG or FEH - and written in uppercase'))}
  
  this.bench.col <- which(colnames(bench.data) == variable)
  this.bench.col2 <- which(colnames(bench.data) == variable2)
  list.of.setups <- levels(as.factor(col.of.setups))
  
  coef.alphas <- summary(the.model,vars=c('alpha'))[,4]   #Median
  coef.alphas.1 <- summary(the.model,vars=c('alpha'))[,1]  #5th percentile
  coef.alphas.3 <- summary(the.model,vars=c('alpha'))[,3]  #9th percentile

  for (each.node in nodes) {
    num.of.node <- which(nodes == each.node)
    for (each.setup in list.of.setups) {
      num.setup <- which(list.of.setups == each.setup)
      print(each.node)
      print(each.setup)
      
      #Reference values in two parameters
      plot.x <- bench.data[the.stars[the.setups == num.setup],this.bench.col]
      plot.x2 <- bench.data[the.stars[the.setups == num.setup],this.bench.col2]
      
      #Delta of observed and reference
      plot.y <- (observed.data[(the.setups == num.setup),num.of.node]-plot.x)
      #str(plot.s)
      
      #SNR vector
      plot.s <- observed.snr[(the.setups == num.setup),num.of.node]
      #str(plot.y)
      
      #Parmeter vectors
      plot.t <- bench.data[the.stars[the.setups == num.setup],'TEFF']
      plot.l <- bench.data[the.stars[the.setups == num.setup],'LOGG']
      plot.f <- bench.data[the.stars[the.setups == num.setup],'FEH']
      
      if (sum(is.na(plot.y)) == length(plot.y)) { 
        plot.x <- c(-5,5)
        plot.y <- plot.x
      }
      
      
      par(mfrow=c(2,2))
      
      #X
      #plot(plot.x,plot.y,pch = 16,col = rgb(0,0,0,0.5),
      #     main=paste0(each.node,' - ',each.setup),xlab=paste0('Given ',variable,' of Reference Stars (per spectrum)'),
      #     ylab=paste0('Delta ',variable,' Observed - Given'))
      
      #Scaled main parameter
      x <- plot.x
      order.x <- order(x)
      x <- x[order.x]
      nor.x <- (x-mean.param.bench)/sd.param.bench
      
      #Scaled second parameter in order of main parameter
      x2 <- plot.x2
      x2 <- x2[order.x]
      nor.x2 <- (x2-mean.param2.bench)/sd.param2.bench
      
      #print(plot.x[order.x])
      #print(plot.y[order.x])
      
      if (n.coeff == 3){
        print('here3')
        k1 <- ((n.coeff*num.of.node)-2)+((num.setup-1)*(length(nodes)*n.coeff))
        k2 <- ((n.coeff*num.of.node)-1)+((num.setup-1)*(length(nodes)*n.coeff))
        k3 <- ((n.coeff*num.of.node))+((num.setup-1)*(length(nodes)*n.coeff))
        
        nor.y <- coef.alphas[k1]+coef.alphas[k2]*nor.x+coef.alphas[k3]*nor.x^2
        y <- nor.y*sd.param.bench
        #lines(x,y,col='red',lwd=3)
        
        y.1 <- sd.param.bench*(coef.alphas.1[k1]+coef.alphas.1[k2]*nor.x+coef.alphas.1[k3]*nor.x^2)
        #lines(x,y.1,col='blue',lwd=3) 
        
        y.3 <- sd.param.bench*(coef.alphas.3[k1]+coef.alphas.3[k2]*nor.x+coef.alphas.3[k3]*nor.x^2)
        #lines(x,y.3,col='blue',lwd=3)
      }else if (n.coeff == 6){
        print('here6')
        k1 <- ((n.coeff*num.of.node)-5)+((num.setup-1)*(length(nodes)*n.coeff))
        k2 <- ((n.coeff*num.of.node)-4)+((num.setup-1)*(length(nodes)*n.coeff))
        k3 <- ((n.coeff*num.of.node)-3)+((num.setup-1)*(length(nodes)*n.coeff))
        k4 <- ((n.coeff*num.of.node)-2)+((num.setup-1)*(length(nodes)*n.coeff))
        k5 <- ((n.coeff*num.of.node)-1)+((num.setup-1)*(length(nodes)*n.coeff))
        k6 <- ((n.coeff*num.of.node)-0)+((num.setup-1)*(length(nodes)*n.coeff))
        
        nor.y <- coef.alphas[k1]+coef.alphas[k2]*nor.x+coef.alphas[k3]*nor.x^2+coef.alphas[k4]*nor.x2+coef.alphas[k5]*nor.x2^2+coef.alphas[k6]*nor.x*nor.x2

        y <- nor.y*sd.param.bench
        #lines(x,y,col='red',lwd=3)
        
        y.1 <- sd.param.bench*(coef.alphas.1[k1]+coef.alphas.1[k2]*nor.x+coef.alphas.1[k3]*nor.x^2+coef.alphas.1[k4]*nor.x2+coef.alphas.1[k5]*nor.x2^2+coef.alphas.1[k6]*nor.x*nor.x2)
        #lines(x,y.1,col='blue',lwd=3) 
        
        y.3 <- sd.param.bench*(coef.alphas.3[k1]+coef.alphas.3[k2]*nor.x+coef.alphas.3[k3]*nor.x^2+coef.alphas.3[k4]*nor.x2+coef.alphas.3[k5]*nor.x2^2+coef.alphas.3[k6]*nor.x*nor.x2)
        #lines(x,y.3,col='blue',lwd=3)
      }
      
      #print(paste('node = ',each.node,' setup = ',each.setup,' k1, k2, k3 =',k1,k2,k3,sep=" "))
      #print(coef.alphas)
      #stop()
      
      #TEFF
      print('teff')
      #str(plot.t)
      #print(length(plot.t))
      #print(length(plot.y))
      plot.tx = plot.t[order.x]  #Puts LOGG in the same order as y (i.e. as x, the independent variable the bias is calculated on)
      order.t <- order(plot.tx)  #Puts LOGG of order x, into LOGG order.
      plot(plot.t,plot.y,pch = 16,col = rgb(0,0,0,0.25),
           main=paste0(each.node,' - ',each.setup),xlab=paste0('Given TEFF of Benchmarks (per spectrum)'),
           ylab=paste0('Delta ',variable,' Observed - Bench'))
      lines(plot.tx[order.t],y[order.t],col='red',lwd=3)
      #lines(plot.tx[order.t],y.1[order.t],col='blue',lwd=3)
      #lines(plot.tx[order.t],y.3[order.t],col='blue',lwd=3)
      #points(plot.t,plot.y,pch = 16,col = rgb(0,0,0,0.25))
      
      
      #LOGG
      print('logg')
      plot.lx = plot.l[order.x]  #Puts LOGG in the same order as y (i.e. as x, the independent variable the bias is calculated on)
      order.l <- order(plot.lx)  #Puts LOGG of order x, into LOGG order.
      plot(plot.l,plot.y,pch = 16,col = rgb(0,0,0,0.25),
           main=paste0(each.node,' - ',each.setup),xlab=paste0('Given LOGG of Benchmarks (per spectrum)'),
           ylab=paste0('Delta ',variable,' Observed - Bench'))
      lines(plot.lx[order.l],y[order.l],col='red',lwd=3)
      #lines(plot.lx[order.l],y.1[order.l],col='blue',lwd=3)
      #lines(plot.lx[order.l],y.3[order.l],col='blue',lwd=3)
      #points(plot.l,plot.y,pch = 16,col = rgb(0,0,0,0.25))
      
      #FEH
      print('feh')
      plot.fx = plot.f[order.x]
      order.f <- order(plot.fx)
      plot(plot.f,plot.y,pch = 16,col = rgb(0,0,0,0.5),
           main=paste0(each.node,' - ',each.setup),xlab=paste0('Given FEH of Benchmarks (per spectrum)'),
           ylab=paste0('Delta ',variable,' Observed - Bench'))
      
      lines(plot.fx[order.f],y[order.f],col='red',lwd=3)
      #lines(plot.fx[order.f],y.1[order.f],col='blue',lwd=3)
      #lines(plot.fx[order.f],y.3[order.f],col='blue',lwd=3)
      
      #SNR
      plot(plot.s,plot.y,pch = 16,col = rgb(0,0,0,0.5),
           main=paste0(each.node,' - ',each.setup),xlab=paste0('Measured SNR (per spectrum)'),
           ylab=paste0('Delta ',variable,' Observed - Bench'))
      
      #plot.tx = plot.t[order.x]
      #order.t <- order(plot.tx)
      #lines(plot.tx[order.t],y[order.t],col='red',lwd=3)
      #lines(plot.tx[order.t],y.1[order.t],col='blue',lwd=3) 
      #lines(plot.tx[order.t],y.3[order.t],col='blue',lwd=3)
      
      par(mfrow=c(1,1))
      
    }
  }
  
}



############################################
# NOW TO PREPARE THE HOMOGENISATION OF THE WHOLE SAMPLE

prepare.for.full.sample.homog.mparam <- function(variable='TEFF',variable2='FEH',this.model,full.metadata=nodes.ids,all.param=clean.nodes.param,bench.data=bench.param,n.coeff=ncoeff) {
  # Will order the stars by CNAME 
  ii <- order(full.metadata$CNAME) 
  my.selec.ids <- full.metadata[ii,]
  my.selec.param <- all.param[ii,,]
  
  # Prepare CNAME, GES_FLD, GES_TYPE, SETUP, and RV per CNAME
  my.selec.cnames <- as.vector(as.character(levels(as.factor(my.selec.ids$CNAME))))
  my.selec.gesfld <- vector('character')
  my.selec.gestype <- vector('character')
  my.selec.setup <- vector('character')
  my.selec.vel <- vector('numeric')
  print(length(my.selec.cnames))

  for (i in 1:length(my.selec.cnames)) {
    my.selec.gesfld[i] <- as.character(levels(as.factor(my.selec.ids$GES_FLD[which(my.selec.ids$CNAME == my.selec.cnames[i])])))
    my.selec.setup[i] <- paste(as.character(levels(as.factor(my.selec.ids$SETUP[which(my.selec.ids$CNAME == my.selec.cnames[i])]))),collapse="|")
    my.selec.gestype[i] <- paste(as.character(levels(as.factor(my.selec.ids$GES_TYPE[which(my.selec.ids$CNAME == my.selec.cnames[i])]))),collapse="|")
    ii <- which(my.selec.ids$CNAME == my.selec.cnames[i])
    if (length(ii) == 1) {
      my.selec.vel[i] <- round(as.numeric(my.selec.ids$VEL[ii]),digits=2)
    } else {
      my.selec.vel[i] <- NA # What to do with RV if there is more than 1 spectrum?
    }
  }
  
  tot.num.spectra <- dim(my.selec.param)[1]
  tot.num.stars <- length(my.selec.cnames)
  num.nodes <- dim(my.selec.param)[2]
  setups.used <- levels(as.factor(my.selec.ids$SETUP))
  num.setups <- length(setups.used)
  
  full.sample.param <- matrix(NA,ncol=dim(my.selec.param)[2],nrow=dim(my.selec.param)[1])
  full.sample.param2 <- matrix(NA,ncol=dim(my.selec.param)[2],nrow=dim(my.selec.param)[1])
  
  spectrum.not.missing <- vector("numeric")
  node.not.missing <- vector("numeric")
  spectrum.missing <- vector("numeric")
  node.missing <- vector("numeric")
  star.code.not.missing <- vector("numeric")
  setup.not.missing <- vector("numeric")
  
  tot.num.not.missing <- 0
  tot.missing <- 0
  
  num.col <- which(dimnames(my.selec.param)[[3]] == variable)
  if (variable == 'TEFF') {
    round.to.this <- 0
  } else {
    round.to.this <- 2
  }
  
  num.col2 <- which(dimnames(my.selec.param)[[3]] == variable2)
  if (variable2 == 'TEFF') {
    round.to.this2 <- 0
  } else {
    round.to.this2 <- 2
  }
  
  for (j in 1:tot.num.spectra) {
    for (i in 1:num.nodes) {
      full.sample.param[j,i] <- round(as.numeric(my.selec.param[j,i,num.col]),digits=round.to.this)
      full.sample.param2[j,i] <- round(as.numeric(my.selec.param[j,i,num.col2]),digits=round.to.this2)
      if (is.na(full.sample.param[j,i])) {
        tot.missing <- tot.missing + 1
        full.sample.param[j,i] <- 10
        full.sample.param2[j,i] <- 10
        spectrum.missing[tot.missing] <- j
        node.missing[tot.missing] <- i
      } else {
        if (is.na(full.sample.param2[j,i])) {
          tot.missing <- tot.missing + 1
          full.sample.param[j,i] <- 10
          full.sample.param2[j,i] <- 10
          spectrum.missing[tot.missing] <- j
          node.missing[tot.missing] <- i
        } else {    
          tot.num.not.missing <- tot.num.not.missing+1
          spectrum.not.missing[tot.num.not.missing] <- j
          node.not.missing[tot.num.not.missing] <- i
          star.code.not.missing[tot.num.not.missing] <- which(my.selec.cnames == as.character(my.selec.ids$CNAME[j]))
          setup.not.missing[tot.num.not.missing] <- which(setups.used == my.selec.ids$SETUP[j])
        }
      }
    }
  }
  
  mean.alpha <- array(data=summary(this.model,vars=c('alpha'))[,4],dim=c(3,num.nodes,num.setups))
  sd.alpha   <- array(data=summary(this.model,vars=c('alpha'))[,5],dim=c(3,num.nodes,num.setups))
  
  mean.node.sd <- round(as.vector(summary(this.model,vars=c('node.sd'))[,4]),digits=round.to.this)
  sd.node.sd <- round(as.vector(summary(this.model,vars=c('node.sd'))[,5]),digits=round.to.this)
  
  mean.Rho <- matrix(data=round(summary(this.model,vars=c('Rho'))[,4],digits=2),ncol=num.nodes,nrow=num.nodes)
  sd.Rho <- matrix(data=round(summary(this.model,vars=c('Rho'))[,5],digits=2),ncol=num.nodes,nrow=num.nodes)
  
  mean.biases <- matrix(0,ncol=ncol(full.sample.param),nrow=nrow(full.sample.param))
  sd.biases <- matrix(0,ncol=ncol(full.sample.param),nrow=nrow(full.sample.param))
  renorm.full.sample.param <- matrix(10,ncol=ncol(full.sample.param),nrow=nrow(full.sample.param))
  renorm.full.sample.param2 <- matrix(10,ncol=ncol(full.sample.param2),nrow=nrow(full.sample.param2))
  
  col.in.bench <- which(colnames(bench.data) == variable)
  mean.bench.param <- mean(bench.data[,col.in.bench],na.rm=T)
  sd.bench.param <- sd(bench.data[,col.in.bench],na.rm=T)

  col2.in.bench <- which(colnames(bench.data) == variable2)
  mean.bench.param2 <- mean(bench.data[,col2.in.bench],na.rm=T)
  sd.bench.param2 <- sd(bench.data[,col2.in.bench],na.rm=T)
  
  for (j in 1:tot.num.not.missing) {

    alpha1 <- rnorm(1000,mean.alpha[1,node.not.missing[j],setup.not.missing[j]],sd.alpha[1,node.not.missing[j],setup.not.missing[j]])
    alpha2 <- rnorm(1000,mean.alpha[2,node.not.missing[j],setup.not.missing[j]],sd.alpha[2,node.not.missing[j],setup.not.missing[j]])
    alpha3 <- rnorm(1000,mean.alpha[3,node.not.missing[j],setup.not.missing[j]],sd.alpha[3,node.not.missing[j],setup.not.missing[j]])
    alpha4 <- rnorm(1000,mean.alpha[1,node.not.missing[j],setup.not.missing[j]],sd.alpha[1,node.not.missing[j],setup.not.missing[j]])
    alpha5 <- rnorm(1000,mean.alpha[2,node.not.missing[j],setup.not.missing[j]],sd.alpha[2,node.not.missing[j],setup.not.missing[j]])
    alpha6 <- rnorm(1000,mean.alpha[3,node.not.missing[j],setup.not.missing[j]],sd.alpha[3,node.not.missing[j],setup.not.missing[j]])

    renorm.full.sample.param[spectrum.not.missing[j],node.not.missing[j]] <- (full.sample.param[spectrum.not.missing[j],node.not.missing[j]]-mean.bench.param)/sd.bench.param
    renorm.full.sample.param2[spectrum.not.missing[j],node.not.missing[j]] <- (full.sample.param2[spectrum.not.missing[j],node.not.missing[j]]-mean.bench.param2)/sd.bench.param2
      
    bias.param <- alpha1 + alpha2 * renorm.full.sample.param[spectrum.not.missing[j],node.not.missing[j]] + alpha3 * ((renorm.full.sample.param[spectrum.not.missing[j],node.not.missing[j]])^2) + alpha4 * renorm.full.sample.param2[spectrum.not.missing[j],node.not.missing[j]] + alpha5 * ((renorm.full.sample.param2[spectrum.not.missing[j],node.not.missing[j]])^2) + alpha6 * renorm.full.sample.param[spectrum.not.missing[j],node.not.missing[j]] * renorm.full.sample.param2[spectrum.not.missing[j],node.not.missing[j]]

    mean.biases[spectrum.not.missing[j],node.not.missing[j]] <- round(mean(bias.param*sd.bench.param),digits=round.to.this)
    sd.biases[spectrum.not.missing[j],node.not.missing[j]]   <- round(sd(bias.param*sd.bench.param),digits=round.to.this)

    #print(full.sample.param[spectrum.not.missing[j],node.not.missing[j]])
    #print(full.sample.param2[spectrum.not.missing[j],node.not.missing[j]])
    #print(renorm.full.sample.param[spectrum.not.missing[j],node.not.missing[j]])
    #print(renorm.full.sample.param2[spectrum.not.missing[j],node.not.missing[j]])
    #print(bias.param)
    #print(sd.bench.param)
    #print(sd.biases[spectrum.not.missing[j],node.not.missing[j]])
    if (sd.biases[spectrum.not.missing[j],node.not.missing[j]] == 0) {
      #print(sd(bias.param*sd.bench.param)) - not rounded
      #print(sd.biases[spectrum.not.missing[j],node.not.missing[j]]) - rounded
      sd.biases[spectrum.not.missing[j],node.not.missing[j]]   <- 0.01
    }
  }
 
  print(tot.num.stars)
  
  all.data <- list(full.sample.param,tot.num.stars,tot.num.spectra,star.code.not.missing,spectrum.not.missing,node.not.missing,tot.num.not.missing,
                   tot.missing,spectrum.missing,node.missing,mean.biases,sd.biases,mean.node.sd,mean.Rho,
                   my.selec.cnames,my.selec.gesfld,my.selec.gestype,my.selec.setup,my.selec.vel)
  return(all.data)
}

####################################################
# Model to homogenise all TEFF values of the sample

homog.all.teff <- "model{

# For each star there is one true.teff that I do not know so I start with a uniform prior - star in the sense of individual CNAME
for (i in 1:N) { # for each individual CNAME
  true.teff[i] ~ dunif(3000,8000)
}

# In my formulation of bias, the node does not sample from a distribution centered in true.teff, but one centered in true.teff + bias.teff
#for (j in 1:M) { # for each individual spectrum
#    for (k in 1:L) { # for each node
#    alpha1[j,k] ~ dnorm(mean.alpha[1,k,setup.code[j]],pow(sd.alpha[1,k,setup.code[j]],-2))
#    alpha2[j,k] ~ dnorm(mean.alpha[2,k,setup.code[j]],pow(sd.alpha[2,k,setup.code[j]],-2))
#    alpha3[j,k] ~ dnorm(mean.alpha[3,k,setup.code[j]],pow(sd.alpha[3,k,setup.code[j]],2))
#    renorm.observed.teff[j,k] <- (observed.teff[j,k]-offset.teff)/scale.teff
#    bias.teff[j,k] <- alpha1[j,k] + alpha2[j,k] * renorm.observed.teff[j,k] + alpha3[j,k] * pow(renorm.observed.teff[j,k],2)
#    vector.of.true.teff[j,k] <- true.teff[star.code[j]] + (bias.teff[j,k] * scale.teff)
#    }
#}

# As above but with treatment for missing values
for (l in 1:tot.num.not.missing) {
#    alpha1[spectrum.not.missing[l],node.not.missing[l]] ~ dnorm(mean.alpha[1,node.not.missing[l],setup.not.missing[l]],pow(sd.alpha[1,node.not.missing[l],setup.not.missing[l]],-2))
#    alpha2[spectrum.not.missing[l],node.not.missing[l]] ~ dnorm(mean.alpha[2,node.not.missing[l],setup.not.missing[l]],pow(sd.alpha[2,node.not.missing[l],setup.not.missing[l]],-2))
#    alpha3[spectrum.not.missing[l],node.not.missing[l]] ~ dnorm(mean.alpha[3,node.not.missing[l],setup.not.missing[l]],pow(sd.alpha[3,node.not.missing[l],setup.not.missing[l]],-2))

#    renorm.observed.teff[spectrum.not.missing[l],node.not.missing[l]] <- (observed.teff[spectrum.not.missing[l],node.not.missing[l]]-offset.teff)/scale.teff

#    bias.teff[spectrum.not.missing[l],node.not.missing[l]] <- alpha1[spectrum.not.missing[l],node.not.missing[l]] + 
#alpha2[spectrum.not.missing[l],node.not.missing[l]] * renorm.observed.teff[spectrum.not.missing[l],node.not.missing[l]] + 
#alpha3[spectrum.not.missing[l],node.not.missing[l]] * pow(renorm.observed.teff[spectrum.not.missing[l],node.not.missing[l]],2)

bias.teff[spectrum.not.missing[l],node.not.missing[l]] ~ dnorm(mean.bias.teff[spectrum.not.missing[l],node.not.missing[l]],pow(sd.bias.teff[spectrum.not.missing[l],node.not.missing[l]],-2))

vector.of.true.teff[spectrum.not.missing[l],node.not.missing[l]] <- true.teff[star.code.not.missing[l]] + bias.teff[spectrum.not.missing[l],node.not.missing[l]]
}

for (l in 1:tot.missing) {

vector.of.true.teff[spectrum.missing[l],node.missing[l]] <- 0

}


# Build the teff.InvCovMatrix
#for (k in 1:L) {
#  node.sd[k] ~ dnorm(mean.node.sd[k],pow(sd.node.sd[k],-2))
#}

for (k1 in 1:L) {
for (k2 in 1:L) {
#        teff.Rho[k1,k2] ~ dnorm(mean.teff.Rho[k1,k2],pow(sd.teff.Rho[k1,k2],-2))
teff.CovMat[k1,k2] <- mean.node.sd[k1]*mean.node.sd[k2]*mean.teff.Rho[k1,k2]
}
}
teff.InvCovMat <- inverse(teff.CovMat)


# And what the nodes give is a noisy measurement of the biased vector.of.true.teff
for (j in 1:M) { # for each individual spectrum

observed.teff[j,1:L] ~ dmnorm( vector.of.true.teff[j,1:L], teff.InvCovMat[1:L,1:L] )

}


}"



####################################################
# Model to homogenise all LOGG values of the sample


homog.all.logg <- "model{

# For each star there is one true.logg that I do not know so I start with a uniform prior - star in the sense of individual CNAME
for (i in 1:N) { # for each individual CNAME
true.logg[i] ~ dunif(0.0,5.5)
}

# In my formulation of bias, the node does not sample from a distribution centered in true.logg, but one centered in true.logg + bias.logg
# As above but with treatment for missing values
for (l in 1:tot.num.not.missing) {

bias.logg[spectrum.not.missing[l],node.not.missing[l]] ~ dnorm(mean.bias.logg[spectrum.not.missing[l],node.not.missing[l]],pow(sd.bias.logg[spectrum.not.missing[l],node.not.missing[l]],-2))

vector.of.true.logg[spectrum.not.missing[l],node.not.missing[l]] <- true.logg[star.code.not.missing[l]] + bias.logg[spectrum.not.missing[l],node.not.missing[l]]
}

for (l in 1:tot.missing) {

vector.of.true.logg[spectrum.missing[l],node.missing[l]] <- 10

}


# Build the logg.InvCovMatrix
#for (k in 1:L) {
#  node.sd[k] ~ dnorm(mean.node.sd[k],pow(sd.node.sd[k],-2))
#}

for (k1 in 1:L) {
for (k2 in 1:L) {
#        logg.Rho[k1,k2] ~ dnorm(mean.logg.Rho[k1,k2],pow(sd.logg.Rho[k1,k2],-2))
logg.CovMat[k1,k2] <- mean.node.sd[k1]*mean.node.sd[k2]*mean.logg.Rho[k1,k2]
}
}
logg.InvCovMat <- inverse(logg.CovMat)


# And what the nodes give is a noisy measurement of the biased vector.of.true.logg
for (j in 1:M) { # for each individual spectrum

observed.logg[j,1:L] ~ dmnorm( vector.of.true.logg[j,1:L], logg.InvCovMat[1:L,1:L] )

}


}"


####################################################
# Model to homogenise all FEH values of the sample

homog.all.feh <- "model{

# For each star there is one true.feh that I do not know so I start with a uniform prior - star in the sense of individual CNAME
for (i in 1:N) { # for each individual CNAME
true.feh[i] ~ dunif(-3.5,0.5)
}

# In my formulation of bias, the node does not sample from a distribution centered in true.feh, but one centered in true.feh + bias.feh
# As above but with treatment for missing values
for (l in 1:tot.num.not.missing) {

bias.feh[spectrum.not.missing[l],node.not.missing[l]] ~ dnorm(mean.bias.feh[spectrum.not.missing[l],node.not.missing[l]],pow(sd.bias.feh[spectrum.not.missing[l],node.not.missing[l]],-2))

vector.of.true.feh[spectrum.not.missing[l],node.not.missing[l]] <- true.feh[star.code.not.missing[l]] + bias.feh[spectrum.not.missing[l],node.not.missing[l]]
}

for (l in 1:tot.missing) {

vector.of.true.feh[spectrum.missing[l],node.missing[l]] <- 10

}


# Build the feh.InvCovMatrix
#for (k in 1:L) {
#  node.sd[k] ~ dnorm(mean.node.sd[k],pow(sd.node.sd[k],-2))
#}

for (k1 in 1:L) {
for (k2 in 1:L) {
#        feh.Rho[k1,k2] ~ dnorm(mean.feh.Rho[k1,k2],pow(sd.feh.Rho[k1,k2],-2))
feh.CovMat[k1,k2] <- mean.node.sd[k1]*mean.node.sd[k2]*mean.feh.Rho[k1,k2]
}
}
feh.InvCovMat <- inverse(feh.CovMat)


# And what the nodes give is a noisy measurement of the biased vector.of.true.feh
for (j in 1:M) { # for each individual spectrum

observed.feh[j,1:L] ~ dmnorm( vector.of.true.feh[j,1:L], feh.InvCovMat[1:L,1:L] )

}


}"

########################################
# Function with inits for full homogenisation

inits.for.full.homog <- function(variable=c('TEFF'),num.chains=num.chains.for.teff,tot.num.stars=data.for.full.homog.teff[[2]],min.param=3000,max.param=8000) {
  variable <- str_trim(variable)
  if (length(variable) != 1) { stop(paste('Choose one variable at a time from TEFF or LOGG or FEH')) }
  if (!(variable %in% c("TEFF","LOGG","FEH"))) { stop(paste('Variable has to be one of TEFF, LOGG or FEH - and written in uppercase'))}
  
  list.of.names <- c('true.teff','true.logg','true.feh')
  name.true <- list.of.names[which(c("TEFF","LOGG","FEH") == variable)]
  
  list.of.inits <- list()
  for (ik in 1:num.chains) {
    list.of.inits[[ik]] <-     list(first = runif(tot.num.stars,min.param,max.param),
                                    .RNG.name="base::Super-Duper", 
                                    .RNG.seed=((2*ik)+1))
    names(list.of.inits[[ik]]) <- c(name.true,'.RNG.name','.RNG.seed')
  }
  return(list.of.inits)
}

############################################
# MAIN FUNCTION TO RUN THE MCMC SIMULATION FOR THE WHOLE DATA

run.the.full.homog <- function(input.data,the.model,num.chains,num.adapt,num.bur,size.sample,thin.factor,variables.to.monitor,model.inits) {
  
  ptm <- proc.time() # To print out how long it took to run (in seconds)
  
  simple_test <- run.jags(method="parallel", # Modified to run chains in parallel
                          data=input.data,
                          model=the.model,
                          n.chains=num.chains,
                          adapt=num.adapt,
                          monitor=variables.to.monitor,
                          inits=model.inits,
                          burnin=num.bur,
                          sample=size.sample,
                          summarise=F,
                          thin=thin.factor,
                          plots=F
  )
  
  print(proc.time() - ptm)
  gc()
  
  return(simple_test)
}


#####################################################
# Create a table with the homog data

create.data.table <- function(homog.teff=full.teff.homog,homog.logg=full.logg.homog,homog.feh=full.feh.homog,metadata.here=data.for.full.homog.teff) {
  fit_teff <- as.mcmc(homog.teff)
  mean.homog.teff <- apply(fit_teff,2,mean)
  sd.homog.teff <- apply(fit_teff,2,sd)
  
  fit_logg <- as.mcmc(homog.logg)
  mean.homog.logg <- apply(fit_logg,2,mean)
  sd.homog.logg <- apply(fit_logg,2,sd)
  
  fit_feh <- as.mcmc(homog.feh)
  mean.homog.feh <- apply(fit_feh,2,mean)
  sd.homog.feh <- apply(fit_feh,2,sd)
  
  print(str(mean.homog.teff))
  print(str(mean.homog.logg))
  print(str(mean.homog.feh))
  print(str(metadata.here[[15]]))
  data.homog <- as.data.frame(cbind(metadata.here[[15]],metadata.here[[16]],metadata.here[[17]],metadata.here[[18]],metadata.here[[19]],mean.homog.teff,
                                  sd.homog.teff,mean.homog.logg,sd.homog.logg,mean.homog.feh,sd.homog.feh))
  
colnames(data.homog) <- c('CNAME','GES_FLD','GES_TYPE','SETUP','VEL','TEFF','E_TEFF','LOGG','E_LOGG','FEH','E_FEH')
data.homog$VEL <- round(as.numeric(as.vector(data.homog$VEL)),digits=2)
data.homog$TEFF <- round(as.numeric(as.vector(data.homog$TEFF)),digits=0)
data.homog$E_TEFF <- round(as.numeric(as.vector(data.homog$E_TEFF)),digits=0)
data.homog$LOGG <- round(as.numeric(as.vector(data.homog$LOGG)),digits=2)
data.homog$E_LOGG <- round(as.numeric(as.vector(data.homog$E_LOGG)),digits=2)
data.homog$FEH <- round(as.numeric(as.vector(data.homog$FEH)),digits=2)
data.homog$E_FEH <- round(as.numeric(as.vector(data.homog$E_FEH)),digits=2)
data.homog$CNAME <- as.character(as.vector(data.homog$CNAME))
data.homog$GES_FLD <- as.character(as.vector(data.homog$GES_FLD))
data.homog$GES_TYPE <- as.character(as.vector(data.homog$GES_TYPE))
data.homog$SETUP <- as.character(as.vector(data.homog$SETUP))
filter.no.teff <- (data.homog$E_TEFF > 1400)
data.homog$TEFF[filter.no.teff] <- NA
data.homog$E_TEFF[filter.no.teff] <- NA
filter.no.logg <- (data.homog$E_LOGG >= 1.5)
data.homog$LOGG[filter.no.logg] <- NA
data.homog$E_LOGG[filter.no.logg] <- NA
filter.no.feh <- (data.homog$E_FEH >= 0.9)
data.homog$FEH[filter.no.feh] <- NA
data.homog$E_FEH[filter.no.feh] <- NA

return(data.homog)
}


plot.the.final.benchmarks <- function(data.homog=all.data.homog,bench.data=bench.param) {
  
  filter.benchmarks <- which(data.homog$GES_FLD %in% bench.data$GES_FLD)
  homog.bench <- data.homog[filter.benchmarks,]
  ii.h <- order(homog.bench$GES_FLD)
  homog.bench <- homog.bench[ii.h,]
  
  ii.bp <- order(bench.data$GES_FLD)
  sorted.bench.param <- bench.data[ii.bp,]
  
  for (each.param in c('TEFF','LOGG','FEH')) {
    this.col <- which(colnames(sorted.bench.param) == each.param)
    that.col <- which(colnames(homog.bench) == each.param)
    
    plot.lim <- c(min(sorted.bench.param[,this.col],homog.bench[,that.col],na.rm=TRUE),max(sorted.bench.param[,this.col],homog.bench[,that.col],na.rm=T))
    
    plot(sorted.bench.param[,this.col],homog.bench[,that.col],xlim=plot.lim,ylim=plot.lim,xlab=paste0(each.param,' (Reference)'),ylab=paste0(each.param,' (Homog)'))
    x <- seq(min(plot.lim,na.rm=T),max(plot.lim,na.rm=T),0.1)
    y <- x
    lines(x,y,lwd=3,col='red')
    
    plot(sorted.bench.param[,this.col],(sorted.bench.param[,this.col]-homog.bench[,that.col]),xlab=paste0(each.param,' (Reference)'),ylab=paste0('Delta ',each.param,' (Reference - Homog)'))

    print(' ')
    print(paste0('Standard Deviation ',each.param,' (Reference - Homog)'))
    print(sd((sorted.bench.param[,this.col]-homog.bench[,that.col]),na.rm=T))
    print(' ')
    print(paste0('Summary ',each.param,' (Reference - Homog)'))
    print(summary((sorted.bench.param[,this.col]-homog.bench[,that.col])))
        
  }
}
  

















#####################################################
# Below are pieces of code that are not used - they wer used only in some tests that are no longer relevant

dummy.function <- function() { # Just to avoid running what is below
  ################################################
  #
  # START HOMOG OF WHOLE SAMPLE
  #
  
  # THIS IS TO RUN WHOLE SAMPLE BUT WITH SORTED CNAMES
  
  
  
  #logg
  
  mean.node.sd.logg <- round(as.vector(summary(simple_test_logg,vars=c('node.sd.logg'))[,4]),digits=2)
  sd.node.sd.logg <- round(as.vector(summary(simple_test_logg,vars=c('node.sd.logg'))[,5]),digits=2)
  mean.logg.Rho <- matrix(data=round(summary(simple_test_logg,vars=c('logg.Rho'))[,4],digits=2),ncol=num.nodes,nrow=num.nodes)
  sd.logg.Rho <- matrix(data=round(summary(simple_test_logg,vars=c('logg.Rho'))[,5],digits=2),ncol=num.nodes,nrow=num.nodes)
  
  
  logg.mean.biases <- matrix(0,ncol=ncol(full.sample.logg),nrow=nrow(full.sample.logg))
  logg.sd.biases <- matrix(0,ncol=ncol(full.sample.logg),nrow=nrow(full.sample.logg))
  renorm.full.sample.logg <- matrix(10,ncol=ncol(full.sample.logg),nrow=nrow(full.sample.logg))
  
  for (j in 1:tot.num.not.missing) {
    alpha1 <- rnorm(1000,mean.alpha.logg[1,node.not.missing[j],setup.not.missing[j]],sd.alpha.logg[1,node.not.missing[j],setup.not.missing[j]])
    alpha2 <- rnorm(1000,mean.alpha.logg[2,node.not.missing[j],setup.not.missing[j]],sd.alpha.logg[2,node.not.missing[j],setup.not.missing[j]])
    alpha3 <- rnorm(1000,mean.alpha.logg[3,node.not.missing[j],setup.not.missing[j]],sd.alpha.logg[3,node.not.missing[j],setup.not.missing[j]])
    renorm.full.sample.logg[spectrum.not.missing[j],node.not.missing[j]] <- (full.sample.logg[spectrum.not.missing[j],node.not.missing[j]]-mean.bench.logg)/sd.bench.logg
    bias.logg <- alpha1 + alpha2 * renorm.full.sample.logg[spectrum.not.missing[j],node.not.missing[j]] + alpha3 * ((renorm.full.sample.logg[spectrum.not.missing[j],node.not.missing[j]])^2)
    logg.mean.biases[spectrum.not.missing[j],node.not.missing[j]] <- round(mean(bias.logg*sd.bench.logg),digits=2)
    logg.sd.biases[spectrum.not.missing[j],node.not.missing[j]] <- round(sd(bias.logg*sd.bench.logg),digits=2)
  }
  
  
  #feh
  
  mean.node.sd.feh <- round(as.vector(summary(simple_test_feh,vars=c('node.sd.feh'))[,4]),digits=2)
  sd.node.sd.feh <- round(as.vector(summary(simple_test_feh,vars=c('node.sd.feh'))[,5]),digits=2)
  mean.feh.Rho <- matrix(data=round(summary(simple_test_feh,vars=c('feh.Rho'))[,4],digits=2),ncol=num.nodes,nrow=num.nodes)
  sd.feh.Rho <- matrix(data=round(summary(simple_test_feh,vars=c('feh.Rho'))[,5],digits=2),ncol=num.nodes,nrow=num.nodes)
  
  
  feh.mean.biases <- matrix(0,ncol=ncol(full.sample.feh),nrow=nrow(full.sample.feh))
  feh.sd.biases <- matrix(0,ncol=ncol(full.sample.feh),nrow=nrow(full.sample.feh))
  renorm.full.sample.feh <- matrix(0,ncol=ncol(full.sample.feh),nrow=nrow(full.sample.feh))
  
  for (j in 1:tot.num.not.missing) {
    alpha1 <- rnorm(1000,mean.alpha.feh[1,node.not.missing[j],setup.not.missing[j]],sd.alpha.feh[1,node.not.missing[j],setup.not.missing[j]])
    alpha2 <- rnorm(1000,mean.alpha.feh[2,node.not.missing[j],setup.not.missing[j]],sd.alpha.feh[2,node.not.missing[j],setup.not.missing[j]])
    alpha3 <- rnorm(1000,mean.alpha.feh[3,node.not.missing[j],setup.not.missing[j]],sd.alpha.feh[3,node.not.missing[j],setup.not.missing[j]])
    renorm.full.sample.feh[spectrum.not.missing[j],node.not.missing[j]] <- (full.sample.feh[spectrum.not.missing[j],node.not.missing[j]]-mean.bench.feh)/sd.bench.feh
    bias.feh <- alpha1 + alpha2 * renorm.full.sample.feh[spectrum.not.missing[j],node.not.missing[j]] + alpha3 * ((renorm.full.sample.feh[spectrum.not.missing[j],node.not.missing[j]])^2)
    feh.mean.biases[spectrum.not.missing[j],node.not.missing[j]] <- round(mean(bias.feh*sd.bench.feh),digits=2)
    feh.sd.biases[spectrum.not.missing[j],node.not.missing[j]] <- round(sd(bias.feh*sd.bench.feh),digits=2)
  }
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
################################################
#
# START HOMOG OF WHOLE SAMPLE
#

# THIS IS TO TEST ONLY BENCHMARKS

all.individual.cnames <- as.vector(as.character(levels(as.factor(nodes.ids$CNAME))))
cnames.of.bench <- levels(as.factor(nodes.ids$CNAME[which(nodes.ids$GES_FLD %in% bench.param$GES_FLD)]))

filter.cnames <- (nodes.ids$CNAME %in% cnames.of.bench)

my.selec.ids <- nodes.ids[filter.cnames,]
my.selec.param <- clean.nodes.param[filter.cnames,,]
my.selec.cnames <- as.vector(as.character(levels(as.factor(my.selec.ids$CNAME))))


tot.num.spectra <- dim(my.selec.param)[1]
all.individual.gesfld <- as.vector(as.character(levels(as.factor(my.selec.ids$GES_FLD))))
tot.num.stars <- length(my.selec.cnames)

full.sample.teff <- matrix(NA,ncol=dim(my.selec.param)[2],nrow=dim(my.selec.param)[1])
#full.sample.teff <- matrix(NA,ncol=num.nodes,nrow=100)
spectrum.not.missing <- vector("numeric")
node.not.missing <- vector("numeric")
spectrum.missing <- vector("numeric")
node.missing <- vector("numeric")
star.code.not.missing <- vector("numeric")
setup.not.missing <- vector("numeric")

tot.num.not.missing <- 0
tot.missing <- 0

for (j in 1:tot.num.spectra) {
  for (i in 1:num.nodes) {
     full.sample.teff[j,i] <- trunc(as.numeric(my.selec.param[j,i,1]))
     if (is.na(full.sample.teff[j,i])) {
       tot.missing <- tot.missing + 1
       full.sample.teff[j,i] <- 0
       spectrum.missing[tot.missing] <- j
       node.missing[tot.missing] <- i
     } else {
       tot.num.not.missing <- tot.num.not.missing+1
       spectrum.not.missing[tot.num.not.missing] <- j
       node.not.missing[tot.num.not.missing] <- i
       star.code.not.missing[tot.num.not.missing] <- which(my.selec.cnames == as.character(my.selec.ids$CNAME[j]))
       setup.not.missing[tot.num.not.missing] <- which(setups.used == my.selec.ids$SETUP[j])
     }
  }
}

mean.alpha.teff <- array(data=summary(simple_test,vars=c('alpha'))[,4],dim=c(3,8,2))
sd.alpha.teff <- array(data=summary(simple_test,vars=c('alpha'))[,5],dim=c(3,8,2))

mean.node.sd.teff <- round(as.vector(summary(simple_test,vars=c('node.sd'))[,4]),digits=0)
sd.node.sd.teff <- round(as.vector(summary(simple_test,vars=c('node.sd'))[,5]),digits=0)
mean.teff.Rho <- matrix(data=round(summary(simple_test,vars=c('teff.Rho'))[,4],digits=2),ncol=num.nodes,nrow=num.nodes)
sd.teff.Rho <- matrix(data=round(summary(simple_test,vars=c('teff.Rho'))[,5],digits=2),ncol=num.nodes,nrow=num.nodes)


teff.mean.biases <- matrix(0,ncol=ncol(full.sample.teff),nrow=nrow(full.sample.teff))
teff.sd.biases <- matrix(0,ncol=ncol(full.sample.teff),nrow=nrow(full.sample.teff))
renorm.full.sample.teff <- matrix(0,ncol=ncol(full.sample.teff),nrow=nrow(full.sample.teff))

for (j in 1:tot.num.not.missing) {
  alpha1 <- rnorm(1000,mean.alpha.teff[1,node.not.missing[j],setup.not.missing[j]],sd.alpha.teff[1,node.not.missing[j],setup.not.missing[j]])
  alpha2 <- rnorm(1000,mean.alpha.teff[2,node.not.missing[j],setup.not.missing[j]],sd.alpha.teff[2,node.not.missing[j],setup.not.missing[j]])
  alpha3 <- rnorm(1000,mean.alpha.teff[3,node.not.missing[j],setup.not.missing[j]],sd.alpha.teff[3,node.not.missing[j],setup.not.missing[j]])
  renorm.full.sample.teff[spectrum.not.missing[j],node.not.missing[j]] <- (full.sample.teff[spectrum.not.missing[j],node.not.missing[j]]-mean.bench.teff)/sd.bench.teff
  bias.teff <- alpha1 + alpha2 * renorm.full.sample.teff[spectrum.not.missing[j],node.not.missing[j]] + alpha3 * ((renorm.full.sample.teff[spectrum.not.missing[j],node.not.missing[j]])^2)
  teff.mean.biases[spectrum.not.missing[j],node.not.missing[j]] <- round(mean(bias.teff*sd.bench.teff),digits=0)
  teff.sd.biases[spectrum.not.missing[j],node.not.missing[j]] <- round(sd(bias.teff*sd.bench.teff),digits=0)
}









#logg

homog.logg.data <- list(
  observed.logg = full.sample.logg, # All logg measurements (nrows = spectra, ncol = nodes)
  #  offset.logg = mean.bench.logg, # the same mean.bench.logg used in the Benchmarks step
  #  scale.logg = sd.bench.logg,    # the same sd.bench.logg used in the Benchmarks step
  N = tot.num.stars,   # total number of individual CNAMEs
  M = tot.num.spectra, # total number of individual spectra
  L = num.nodes,       # number of nodes
  star.code.not.missing = star.code.not.missing,
  #  setup.not.missing = setup.not.missing,
  spectrum.not.missing = spectrum.not.missing,
  node.not.missing = node.not.missing,
  tot.num.not.missing = tot.num.not.missing,
  tot.missing = tot.missing,
  spectrum.missing = spectrum.missing,
  node.missing = node.missing,
  mean.bias.logg = logg.mean.biases,
  sd.bias.logg = logg.sd.biases,
  #    star.code = full.sample.star.code,   # this code connects each spectrum to its CNAME. With this, multiple observations will still contribute to the same CNAME
  #  setup.code = full.sample.setup.code, # this code is needed to select which bias to use (U520 or U580)
  #  mean.alpha = mean.alpha.logg, # the mean coefficients for the bias as determined in a previous step
  #  sd.alpha = sd.alpha.logg, # the sd of the mean coefficients for the bias as determined in a previous step
  mean.node.sd = mean.node.sd.logg, # the node typical SD in logg as determined in a previous step
  mean.logg.Rho = mean.logg.Rho#, # the correlation coefficients in logg as determined in a previous step
  #  sd.node.sd = sd.node.sd.logg,
  #  sd.logg.Rho = sd.logg.Rho
)




#feh


homog.feh.data <- list(
  observed.feh = full.sample.feh, # All feh measurements (nrows = spectra, ncol = nodes)
  #  offset.feh = mean.bench.feh, # the same mean.bench.feh used in the Benchmarks step
  #  scale.feh = sd.bench.feh,    # the same sd.bench.feh used in the Benchmarks step
  N = tot.num.stars,   # total number of individual CNAMEs
  M = tot.num.spectra, # total number of individual spectra
  L = num.nodes,       # number of nodes
  star.code.not.missing = star.code.not.missing,
  #  setup.not.missing = setup.not.missing,
  spectrum.not.missing = spectrum.not.missing,
  node.not.missing = node.not.missing,
  tot.num.not.missing = tot.num.not.missing,
  tot.missing = tot.missing,
  spectrum.missing = spectrum.missing,
  node.missing = node.missing,
  mean.bias.feh = feh.mean.biases,
  sd.bias.feh = feh.sd.biases,
  #    star.code = full.sample.star.code,   # this code connects each spectrum to its CNAME. With this, multiple observations will still contribute to the same CNAME
  #  setup.code = full.sample.setup.code, # this code is needed to select which bias to use (U520 or U580)
  #  mean.alpha = mean.alpha.feh, # the mean coefficients for the bias as determined in a previous step
  #  sd.alpha = sd.alpha.feh, # the sd of the mean coefficients for the bias as determined in a previous step
  mean.node.sd = mean.node.sd.feh, # the node typical SD in feh as determined in a previous step
  mean.feh.Rho = mean.feh.Rho#, # the correlation coefficients in feh as determined in a previous step
  #  sd.node.sd = sd.node.sd.feh,
  #  sd.feh.Rho = sd.feh.Rho
)



inits.homog <- function() {
  list(
  true.teff = runif(tot.num.stars,3000,8000)
  )
}

inits1.homog <- inits.homog()
inits2.homog <- inits.homog()
inits3.homog <- inits.homog()
inits4.homog <- inits.homog()

inits.homog.logg <- function() {
  list(
    true.logg = runif(tot.num.stars,0.0,5.5)
  )
}

inits1.homog.logg <- inits.homog.logg()
inits2.homog.logg <- inits.homog.logg()
inits3.homog.logg <- inits.homog.logg()
inits4.homog.logg <- inits.homog.logg()


inits.homog.feh <- function() {
  list(
    true.feh = runif(tot.num.stars,-3.0,0.5)
  )
}

inits1.homog.feh <- inits.homog.feh()
inits2.homog.feh <- inits.homog.feh()
inits3.homog.feh <- inits.homog.feh()
inits4.homog.feh <- inits.homog.feh()


homog_test <- run.jags(method="rjags",
                       data=homog.teff.data,
                       model=homog.all.teff,
                       n.chains=4,
                       adapt=1000,
                       monitor=c("true.teff"),
                       inits=list(inits1.homog,inits2.homog,inits3.homog,inits4.homog),
                       burnin = 1000,
                       sample = 1000,
                       summarise = F,
                       thin=1,
                       plots=F
)


homog_logg_test <- run.jags(method="rjags",
                       data=homog.logg.data,
                       model=homog.all.logg,
                       n.chains=4,
                       adapt=1000,
                       monitor=c("true.logg"),
                       inits=list(inits1.homog.logg,inits2.homog.logg,inits3.homog.logg,inits4.homog.logg),
                       burnin = 1000,
                       sample = 1000,
                       summarise = F,
                       thin=1,
                       plots=F
)

homog_feh_test <- run.jags(method="rjags",
                            data=homog.feh.data,
                            model=homog.all.feh,
                            n.chains=4,
                            adapt=1000,
                            monitor=c("true.feh"),
                            inits=list(inits1.homog.feh,inits2.homog.feh,inits3.homog.feh,inits4.homog.feh),
                            burnin = 1000,
                            sample = 1000,
                            summarise = F,
                            thin=1,
                            plots=F
)

fit_homog <- as.mcmc(homog_test)
mean.homog.teff <- apply(fit_homog,2,mean)
sd.homog.teff <- apply(fit_homog,2,sd)

fit_homog_logg <- as.mcmc(homog_logg_test)
mean.homog.logg <- apply(fit_homog_logg,2,mean)
sd.homog.logg <- apply(fit_homog_logg,2,sd)

fit_homog_feh <- as.mcmc(homog_feh_test)
mean.homog.feh <- apply(fit_homog_feh,2,mean)
sd.homog.feh <- apply(fit_homog_feh,2,sd)

data.homog <- as.data.frame(cbind(my.selec.cnames,my.selec.gesfld,my.selec.gestype,my.selec.setup,my.selec.vel,mean.homog.teff,
                                  sd.homog.teff,mean.homog.logg,sd.homog.logg,mean.homog.feh,sd.homog.feh))
colnames(data.homog) <- c('CNAME','GES_FLD','GES_TYPE','SETUP','VEL','TEFF','E_TEFF','LOGG','E_LOGG','FEH','E_FEH')
data.homog$VEL <- round(as.numeric(as.vector(data.homog$VEL)),digits=2)
data.homog$TEFF <- round(as.numeric(as.vector(data.homog$TEFF)),digits=0)
data.homog$E_TEFF <- round(as.numeric(as.vector(data.homog$E_TEFF)),digits=0)
data.homog$LOGG <- round(as.numeric(as.vector(data.homog$LOGG)),digits=2)
data.homog$E_LOGG <- round(as.numeric(as.vector(data.homog$E_LOGG)),digits=2)
data.homog$FEH <- round(as.numeric(as.vector(data.homog$FEH)),digits=2)
data.homog$E_FEH <- round(as.numeric(as.vector(data.homog$E_FEH)),digits=2)
data.homog$CNAME <- as.character(as.vector(data.homog$CNAME))
data.homog$GES_FLD <- as.character(as.vector(data.homog$GES_FLD))
data.homog$GES_TYPE <- as.character(as.vector(data.homog$GES_TYPE))
data.homog$SETUP <- as.character(as.vector(data.homog$SETUP))
filter.no.teff <- (data.homog$E_TEFF > 1400)
data.homog$TEFF[filter.no.teff] <- NA
data.homog$E_TEFF[filter.no.teff] <- NA
filter.no.logg <- (data.homog$E_LOGG >= 1.5)
data.homog$LOGG[filter.no.logg] <- NA
data.homog$E_LOGG[filter.no.logg] <- NA
filter.no.feh <- (data.homog$E_FEH >= 0.9)
data.homog$FEH[filter.no.feh] <- NA
data.homog$E_FEH[filter.no.feh] <- NA

filter.benchmarks <- which(data.homog$GES_FLD %in% bench.param$GES_FLD)
homog.bench <- data.homog[filter.benchmarks,]
ii.h <- order(homog.bench$GES_FLD)
homog.bench <- homog.bench[ii.h,]

ii.bp <- order(bench.param$GES_FLD)
sorted.bench.param <- bench.param[ii.bp,]

plot(sorted.bench.param$TEFF,homog.bench$TEFF,xlim=c(3000,8000),ylim=c(3000,8000))
x <- seq(3000,8000,10)
y <- x
lines(x,y,lwd=3,col='red')
plot(sorted.bench.param$TEFF,(sorted.bench.param$TEFF-homog.bench$TEFF))
sd((sorted.bench.param$TEFF-homog.bench$TEFF))
summary((sorted.bench.param$TEFF-homog.bench$TEFF))

plot(sorted.bench.param$LOGG,homog.bench$LOGG,xlim=c(0,5),ylim=c(0,5))
x <- seq(-5,5,0.1)
y <- x
lines(x,y,lwd=3,col='red')
plot(sorted.bench.param$LOGG,(sorted.bench.param$LOGG-homog.bench$LOGG))
sd((sorted.bench.param$LOGG-homog.bench$LOGG))
summary((sorted.bench.param$LOGG-homog.bench$LOGG))

plot(sorted.bench.param$FEH,homog.bench$FEH,xlim=c(-3,0.5),ylim=c(-3,0.5))
x <- seq(-5,5,0.1)
y <- x
lines(x,y,lwd=3,col='red')
plot(sorted.bench.param$FEH,(sorted.bench.param$FEH-homog.bench$FEH))
sd((sorted.bench.param$FEH-homog.bench$FEH),na.rm=T)
summary((sorted.bench.param$FEH-homog.bench$FEH))

}

#par(def.par)#- reset to default




# R. Smiljanic
# (start, end) = (Apr 2018, Jul 2018)
# Gaia-ESO iDR6
# Homogenisation code for the atmospheric parameters
# Bayesian MCMC inference of the CovMatrix of the Nodes (i.e., sigmas and correlations of the atmospheric parameters) and the bias functions
# The CovMatrix and bias functions can then be applied on the "other stars" to estimate the best value of their atmospheric parameters

# THIS IS A PRELIMINARY VERSION. IT RUNS, BUT THE RESULTS ARE LESS THAN OPTIMAL. IT IS ASSUMED THAT THE NODE RESULTS ARE AFFECTED BY A BIAS FUNCTION (AN OFFSET) AND A STOCHASTIC ERROR. THE PARAMETRIZATION OF THE BIAS AS A FUNCTION OF TEFF ONLY IS PRELIMINARY. THE STOCHASTIC ERROR SHOULD EVENTUALLY BE PARAMETRIZED AS A FUNCTION OF SOMETHING AS WELL (e.g., ?SNR?), BUT FOR NOW IT IS NOT. NODE FLAGS ARE NOT TAKEN INTO ACCOUNT. MY CHOICE OF PRIORS HAS TO BE REVIEWED AND OPTMIZED.

# I am not renormalizing the variables (i.e., subtracting the mean and deviding by the standard deviation). It is not clear for me which mean and sd to choose (from the benchmarks? from the nodes? all nodes? or per node?).

#################################
# Load the libraries needed to run the MCMC Bayesian analysis and do some plots

library(rjags);library(runjags);library(coda);library(jagstools);library(lattice);library(ggplot2);library(ggthemes)#;library(plot3D)

# Save current parameters for late use 

def.par <- par(no.readonly = TRUE)

# Source the functions needed to read the Node and Benchmark files

source('~/Survey/Resultados/DR6/Codes/part_of_homog_plots.R')

# Path and file with atmospheric parameters of the benchmark stars; columns to be extracted from the benchmarks

bench.path <- paste('/Users/rodolfo/Survey/Resultados/DR5/Support/')
file.bench <- paste('GES_iDR5_FGKMCoolWarm_Benchmarks_AcceptedParams_01082016.fits')
columns.to.read <- c('GES_FLD','GES_TYPE','TEFF','E_TEFF','LOGG','E_LOGG','FEH','DELTA_ION','DELTA_LTE','XI') # Need the real errors of [Fe/H]; in the meantime using will be using the sqrt(delta_ion^2 + delta_lte^2)

# Load the parameters of the benchmark stars

bench.param <- load.bench.fits(fitsname=file.bench,columns=columns.to.read,path.file=bench.path,only.fgk=TRUE)
errors.metallicity <- round(sqrt(bench.param$DELTA_ION**2 + bench.param$DELTA_LTE**2),digits=2) # Proper errors are not provided. Using delta_ion and delta_lte for now
errors.metallicity[errors.metallicity < 0.05] <- 0.05 # Because this is not the proper error, assume a lower limit of 0.05 dex
new.col.names <- c(colnames(bench.param),'E_FEH')
bench.param <- cbind(bench.param,errors.metallicity) # This is the matrix that holds the parameters of the benchmarks
colnames(bench.param) <- new.col.names
rm(errors.metallicity,new.col.names)

# Load the parameters of the Nodes

release <- c('iDR5')
list.nodes <- c('CAUP','EPINARBO','Lumba','Nice','OACT','UCM','Vilnius')
node.files.path <- c('/Users/rodolfo/Survey/Resultados/DR5/Parameters/')

out.list <- load.nodes(list.nodes,release,node.files.path) # Reads the node files and outputs a list that is broken into 3 below:
nodes.param <- out.list[[2]]              # Array (nrow = num.stars, ncol = num.nodes, n.third = parameters [TEFF, E_TEFF, LOGG, E_LOGG, FEH, E_FEH, XI, E_XI] )
nodes.ids <- as.data.frame(out.list[[1]]) # The metadata associated to nodes.param (nrow = numstars, ncol = it varies; includes cnames, setup, snr, etc)
nodes.flags <-  out.list[[3]]             # Array (nrow = num.stars, ncol = num.nodes, n.third = flag) where flag = PECULI, REMARK, TECH in that order
rm(out.list,release)

# To make sure that we use only stars from bench.param that are also included in the nodes.ids

filter.bench <- (bench.param$GES_FLD %in% nodes.ids$GES_FLD)
bench.param <- bench.param[filter.bench,]

# From the Node results, select all and only entries of the benchmarks

filter.bench <- (nodes.ids$GES_FLD %in% bench.param$GES_FLD)
node.measured.param.bench <- nodes.param[filter.bench,,]
metadata.of.bench.spectra <- nodes.ids[filter.bench,]

# Now, let's define all data structures that will be used in the model
# 1) Numbers

num.nodes   <- dim(nodes.param)[2]             # Number of Nodes
num.bench   <- nrow(bench.param)               # Number of individual benchmark stars
num.spectra <- nrow(metadata.of.bench.spectra) # Total number of spectra of benchmark stars

# 2) Series of Vector with the SNR of all benchmark spectra that were analyzed (in each row it repeats the SNR value of that spectrum by n times, where n = num.nodes)

#This is needed if we are going to parametize something against SNR

snr.spec <- as.numeric(as.vector(metadata.of.bench.spectra$SNR))
# There are some weird numbers of SNR > 500 up to ~25000! I am not sure if these should be believed.
# At this moment I am converting this very high values into something more "believable"; anything above 900 gets a random SNR between 900 and 1000 (truncated integer)
filter.high.snr <- (snr.spec > 900)
snr.spec[filter.high.snr] <- trunc(runif(snr.spec[filter.high.snr],900,1000))
snr.spec.vec <- matrix(rep(snr.spec,num.nodes),ncol=num.nodes,nrow=length(snr.spec))

# 3) Create the vectors with the given reference parameters

given.teff.bench <- bench.param$TEFF
given.sigma.teff.bench <- bench.param$E_TEFF

# 4) The measurements from the nodes

observed.node.teff.spectrum <- matrix(NA,ncol=ncol(node.measured.param.bench),nrow=nrow(node.measured.param.bench))
for (ik in seq(1,ncol(observed.node.teff.spectrum))) {
    observed.node.teff.spectrum[,ik] <- trunc(as.numeric(as.vector(node.measured.param.bench[,ik,'TEFF']))) # To be sure they are numbers
}

# Optional test removing all NAs
#observed.node.teff.spectrum <- na.omit(observed.node.teff.spectrum)
#all.row.num <- seq(1,nrow(metadata.of.bench.spectra),1)
#filter.remove <- !(all.row.num %in% na.action(observed.node.teff.spectrum))
#metadata.of.bench.spectra <- metadata.of.bench.spectra[filter.remove,]
#num.spectra <- nrow(metadata.of.bench.spectra)

# Define star.code: it traces back to which benchmark a given line in observed.node.teff.spectrum corresponds to
star.code <- vector('numeric',length=nrow(observed.node.teff.spectrum))
for (ik in seq(1,nrow(observed.node.teff.spectrum))) {
    star.code[ik] <- which(bench.param$GES_FLD == metadata.of.bench.spectra$GES_FLD[ik])
}

# Now, here comes an important trick. There are many missing values (NA or NaN) in the 
#measurements given by the nodes. Missing values need special treatment.

#What we are trying to do boils down to inferring "sigma" given that 
#a node measured "x[i]" and knowing that the real values were "y[i]" for each star "i", 
#or in R jags notation: x[i] ~ dnorm ( mu = y[i] , tau = 1/sigma^2 ); 
#but it happens that sometimes x[i] is missing because the node failed to analyse that star. 
#What I would like to do is, in these cases of missing measurement, simply assume 
#a broad non informative prior for the missing x[i]: x[i] ~ dunif(3000,8000); 
#in other words: ok, the node did not provide a measurement, but I know the value 
#should be between 3000 and 8000 K. **BUT**, I can not write these two lines of code together: 
#x[i] ~ dnorm ( mu = y[i] , tau = 1/sigma^2 ) and x[i] ~ dunif(3000,8000) because JAGS correctly
#thinks that I am trying to define the same "node" x[i] twice and throws an error. 
#On the other hand, JAGS knows how to deal with the situation where it is "y[i]" that has 
#the missing values. If we write x[i] ~ dnorm ( mu = y[i] , tau = 1/sigma^2 ) and 
#y[i] ~ dunif(3000,9000) it automatically knows that it should treat **ONLY** the 
#missing values of y[i] with the non-informative priors (really! I tested it. 
#If y = c(1,2,NA,3,4), with the y[i] ~ dunif(3000,9000) it will substitute only the third element 
#of y for the prior, the rest remains unchanged.

# So, what is the numerical solution for this? 
#I am proposing the following trick where the values to be randomized are inverted. 
#Let's say that the measured values are x = c(1.1,1.9,NA,4.2,4.8,6.3) and the 
#true values are y = c(1,2,3,4,5,6) and I want to infer the sigma as before in 
#x[i] ~ dnorm ( mu = y[i] , tau = 1/sigma^2 ). 
#I would like to treat the missing value with a prior x[i] ~ dunif (0,10), but as I said, 
#this can not be done directly. So, what I do is, I manipulate x and y by taking out 
#the NA from "x" and inputing there the corresponding true value; and removing the 
#corresponding true value from y and inputing there the NA: 
#new.x = c(1.1,1.9,3,4.2,4.8,6.3) and new.y = c(1,2,NA,4,5,6) telling JAGS that when 
#it sees the NA in new.y it should use new.y ~ dunif(0,10). So, instead of randomizing 
#the measured value, I am randomizing the true value. In essence, I think the two formulations 
#are equivalent in what concerns the determination of the "sigma" that we are interested in. 
#Before, I was trying to say that when the node measured y[3], it got a random value x[3] between 0 and 10. 
#Now, what I am saying that it measured x[3] out of a random value of y[3] between 0 and 10.

# This is the new matrix without missing values, where the NAs and NaNs where substituted by 
#the reference benchmak values
manipulated.node.teff.spectrum <- matrix(NA,ncol=ncol(observed.node.teff.spectrum),nrow=nrow(observed.node.teff.spectrum))

# And I also need to save the position (i,j) in the matrix of each element that was changed 
#from NA to a value, those that were not changed, and the star.code of the changed values.
#(Perhaps there is a better and less verbose way to do this, but this was the solution I could 
#come up with)
i.of.missing.teff <- vector()
j.of.missing.teff <- vector()
i.not.missing.teff <- vector()
j.not.missing.teff <- vector()
code.bench.missing <- vector()     # The star.code of the benchmark for which the measurement was missing
code.bench.not.missing <- vector() # The star.code of the benchmark for which measurements were not missing
ik <- 1
jk <- 1
for (i in seq(1,nrow(observed.node.teff.spectrum))) {
    for (j in seq(1,ncol(observed.node.teff.spectrum))) {
        if (is.na(observed.node.teff.spectrum[i,j])) {
            manipulated.node.teff.spectrum[i,j] <- bench.param$TEFF[star.code[i]]
            i.of.missing.teff[ik] <- i
            j.of.missing.teff[ik] <- j
            code.bench.missing[ik] <- star.code[i]
            ik <- ik+1
        } else {
            manipulated.node.teff.spectrum[i,j] <-observed.node.teff.spectrum[i,j]
            i.not.missing.teff[jk] <- i
            j.not.missing.teff[jk] <- j
            code.bench.not.missing[jk] <- star.code[i]
            jk <- jk+1
        }
    }
}

num.missing.values <- ik-1 # save the total number of missing entries
num.not.missing.values <- jk-1 # save the total number of non-missing entries

# 5) The Inverse Covariance Matrix
#  We will model the problem using a multi variate Gaussian where we know the vector of measured 
#values of the node, know the vector with true.teff of the benchmark and want to infer 
#the "Inverse Covariance Matrix". Actually we want the Covariance Matrix, but in the 
#JAGS notation we write the problem with the Inverse Covariance Matrix:
#  observed.node.teff[j,1:K] ~ dmnorm( true.teff.benchmark[j,1:K] , teff.InvCovMat[1:K,1:K] )

# The Inverse Covariance Matrix needs a prior.
# From what I read in the internet, the recommended prior on the inverse covariance matrix 
#is a generic Wishart distribution, which is a multidimensional generalization of a gamma 
#distribution. W( V, n ) takes as input the degrees of freedom "n" and V = Sigma_0 * n^-1 
#where Sigma_0 is a first guess of the covariance matrix. The non-informative choice for 
#the degrees of freedom is to set it equal to the number of rows (in this case the number 
#of nodes).

prior.deg.freedom = num.nodes
prior.matrix = (1/num.nodes)*diag(100^2,ncol=num.nodes,nrow=num.nodes) # Assuming that the typical sigma of Teff = 100 K

# Some links:
#http://doingbayesiandataanalysis.blogspot.com/2017/06/bayesian-estimation-of-correlations-and.html
#http://www.themattsimpson.com/2012/08/20/prior-distributions-for-covariance-matrices-the-scaled-inverse-wishart-prior/
#https://stats.stackexchange.com/questions/66338/how-to-specify-the-wishart-distribution-scale-matrix
#https://en.wikipedia.org/wiki/Wishart_distribution

#6) We are going to fit one bias function per setup (i.e., U520, U580 - this can eventually be 
#extended for Giraffe, but maybe not easilly. In WG11 all nodes analyze all spectra of all 
#setups. That is not the case in WG10. This would break the symmetry of the data structures 
#per setup). We create a vector with a code for each setup:

# Finds out all the setups. Below it converts each setup "name" to a nummeric code
setups.used <- levels(as.factor((metadata.of.bench.spectra$SETUP)))    
vector.of.setups <- as.character(as.vector(metadata.of.bench.spectra$SETUP))
ik <- 0
for (each.setup in setups.used) {
    ik <- ik+1
    vector.of.setups[(vector.of.setups == each.setup)] <- as.character(ik)
}
vector.of.setups <- as.numeric(as.vector(vector.of.setups)) 
num.setups <- length(setups.used) # Total number of setups

# Now, let's create the data list needed to run the model

teff.data <- list(given.teff.benchmarks = given.teff.bench,                       # Benchmark Teff
                  given.error.teff.benchmarks = given.sigma.teff.bench,          # Benchmark e_Teff
                  observed.node.teff.per.spec = manipulated.node.teff.spectrum,  # Node Teff
                  star.code = star.code,   # To trace to which benchmark each line of teff.by.node corresponds to
#                  snr.spec.vec = snr.spec.vec, # SNR
                  K = num.nodes,           #N.nodes,
                  N = num.bench,           #N.benchmarks
                  M = num.spectra,         #N.spectra
                  P = num.setups,          #N.setups
                  setup.code = vector.of.setups, # To trace to which setup the spectrum belongs                  
                  # The missing and not missing values data
                  n.missing = num.missing.values, #N of missing values
                  n.not.missing = num.not.missing.values, # N of NOT missing values
                  i.of.missing.values = i.of.missing.teff,
                  j.of.missing.values = j.of.missing.teff,
                  i.not.missing = i.not.missing.teff,
                  j.not.missing = j.not.missing.teff,
#                  code.bench.missing = code.bench.missing,
                  code.bench.not.missing = code.bench.not.missing,
                  #For wishart (dwish) prior on inverse covariance matrix:
                  prior.matrix = prior.matrix,
                  prior.deg.freedom = prior.deg.freedom
                           )


# And now let's write the model itself

model.teff.matrix.bias <- "model{
# We set up the model to use multivariate normal distributions. 
#For a given spectrum of the benchmarks, the vector of Teff measured by the nodes 
#teff.nodes.vector = (Teff1, Teff2,...,Teffn) (1:n are the nodes) is a random drawing from 
#a multivariate normal distribution with a mean = vector of true.teff of that benchmark 
#(true.teff, true.teff,...,true.teff) - where obviously the true.teff is the same 
#irrespective of node - and there is a covariance matrix (JAGS expects the inverse 
#covariance matrix actually): teff.nodes.vector ~ dmnorm ( mu = true.teff.vector, InvCovMatrix ). 
#So we need to first set up all vectors accordingly.

# 1)
# This first part is to define the vector with TRUE values of the benchmark Teff
# We do not use the given Teff value directly, but use the given value as a prior for the true value

# Prior on the true Teff of the benchmarks
for (i in 1:N) { # Running over N which is the num.bench

    the.true.teff.bench[i] ~ dnorm(given.teff.benchmarks[i], given.tau.teff.bench[i])  # The given benchmark Teff value is a prior on the true value of the star's Teff
    given.tau.teff.bench[i] <- pow(given.error.teff.benchmarks[i],-2)                  # dnorm uses the precision, tau, which is the inverse of the variance: tau = 1 / (sigma^2). So, from here we also believe on the error of the Teff value

}

# Now we will create the vector of true Teffs: vector.true.teff = (true.teff, ..., true.teff)
# Here we have to remember that we use directly true.teff only when the node provided a measurement. 
#If the value was missing, we actually use a broad non-informative prior to simulate that the value 
#of the node measurement was unknown. 
#This means that the actual vector.true.teffs is not one per benchmark, but one per spectrum! 
#Because for each spectrum there is a different number of node values missing.

#true.teff.benchmark.per.spec is a matrix with nrow = number of spectrum and ncol = number of nodes, 
#equivalent to the input node.teff matrix

for (i in 1:n.not.missing) { # When the node gave a value (n.not.missing), we include the true.teff of the benchmark estimated above

    true.teff.benchmark.per.spec[i.not.missing[i],j.not.missing[i]] <- the.true.teff.bench[code.bench.not.missing[i]]

}

#Prior on the missing values

for (l in 1:n.missing) { # When the node was missing a value (n.missing), we include a value drawn from the broad prior

    true.teff.benchmark.per.spec[i.of.missing.values[l],j.of.missing.values[l]] ~ dunif(3000,8000)

}

# Prior on the teff.InvCovMat

teff.InvCovMat ~ dwish( prior.matrix[1:K,1:K] , prior.deg.freedom )

# Other Priors needed for coefficients
# Priors on the coefficients alpha for the fit of the bias function. We fit one bias function per node. The bias.teff is for now a quadratic function of the true.teff
    for (l in 1:K) {
       for (m in 1:P) {

            alpha1[l,m] ~ dnorm(0,1)   # Remember here that dnorm(mu,tau) where mu = mean and tau = 1/sigma^2, and sigma is the standard deviation 
            alpha2[l,m] ~ dnorm(0,1E5)
            alpha3[l,m] ~ dnorm(0,1E6)

        }
    }

# I will eventually need also a vector that really has only the the.true.teff.bench[i], without taking into account the missing values, for the fitting of the bias function (see below)

for (i in 1:M) {
   for (j in 1:K) {
       only.the.true.teff.benchmark[i,j] <- the.true.teff.bench[star.code[i]]
   }
}

# 2)
# Now, we finally write that what the nodes give is a noisy measurement of the true.teff.benchmark

#Likelihood - per spectrum (vectorized over K = num.nodes)

for (j in 1:M) {
    observed.node.teff.per.spec[j,1:K] ~ dmnorm( (true.teff.benchmark.per.spec[j,1:K] + bias.vector.teff.node[j,1:K]) , teff.InvCovMat[1:K,1:K] )
     
     bias.vector.teff.node[j,1:K] <- alpha1[1:K,setup.code[j]] + alpha2[1:K,setup.code[j]]*only.the.true.teff.benchmark[j,1:K] + alpha3[1:K,setup.code[j]]*only.the.true.teff.benchmark[j,1:K]^2

}

# Notice what was done with the bias function above. Let's assume, for example, that I know that everytime I measure 
# the Teff of the Sun, my value comes out with a bias of +50 K. I could write that my measurement minus the bias was drawn from 
# the normal distribution: (my.measurement - 50 K) ~ dnorm(mu = 5777 K, sigma = my.random.error). 
# What I am assuming is that this is equivalent to say that I am instead drawing from a normal distribution 
# with mu = 5827 K (5777 + 50): my.measurement ~ dnorm(mu = 5827 K + my.bias , sigma = my.random.error)

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

# Intial values for free parameters
inits <- function() {
list(alpha1 = matrix(rnorm(num.setups*num.nodes,0,1),ncol=num.setups,nrow=num.nodes),    # Here, outside JAGS, rnorm actually uses mu and sd, the standard deviation.
     alpha2 = matrix(rnorm(num.setups*num.nodes,0,0.01),ncol=num.setups,nrow=num.nodes),
     alpha3 = matrix(rnorm(num.setups*num.nodes,0,0.001),ncol=num.setups,nrow=num.nodes),
     teff.InvCovMat = diag(rgamma(num.nodes,0.1,0.1))
)
}

# Two sets of inital values, because I am running two chains
inits1 <- inits()
inits2 <- inits()

simple_test <- run.jags(method="rjags",
                        data=teff.data,
                        mode=model.teff.matrix.bias,
                        n.chains=2,
                        adapt=1000,
                        monitor=c("the.true.teff.bench","node.sd.teff","teff.Rho","alpha1","alpha2","alpha3"),
                        inits=list(inits1,inits2),
                        burnin=1000,
                        sample=1000,
                        summarise=F,
                        thin=1,
                        plots=F
                        )

# This outputs a summary of the MCMC of these parameters
print(summary(simple_test,vars=c('the.true.teff.bench')))
print(summary(simple_test,vars=c('node.sd')))
print(summary(simple_test,vars=c('teff.Rho')))
print(summary(simple_test,vars=c('alpha1','alpha2','alpha3')))


# Let's look at the covariance matrix

print(paste(' '))
print(paste('This is the final covariance matrix:'))
print(paste(' '))
estimated.cov.matrix <- summary(simple_test,vars=c('teff.Rho'))
estimated.cov.matrix <- matrix(round(estimated.cov.matrix,digits=2),num.nodes,num.nodes)
colnames(estimated.cov.matrix) <- list.nodes
rownames(estimated.cov.matrix) <- list.nodes
print(estimated.cov.matrix)
print(paste(' '))

# This will plot many figures, one at a time. Will need to flip between them to see them all.

# These are examples of plots to see the tracing and posterior
fit_mcmc <- as.mcmc(simple_test)
plot(xyplot(fit_mcmc[,73:78]))
plot(densityplot(fit_mcmc[,73:78]))

plot(xyplot(fit_mcmc[,31:36]))
plot(densityplot(fit_mcmc[,31:36]))

# This will plot the Delta (given Teff of benchmarks minus estimated Teff of benchmarks) as a function Teff, logg and [Fe/H], so we can evaluate the trends that are left

results.teff <- summary(simple_test,vars=c('the.true.teff.bench'))
plot(bench.param$TEFF,(bench.param$TEFF-results.teff[,4]),ylab='Delta',xlab='Given Teff of Benchmarks')

plot(bench.param$LOGG,(bench.param$TEFF-results.teff[,4]),ylab='Delta',xlab='Given Teff of Benchmarks')

plot(bench.param$FEH,(bench.param$TEFF-results.teff[,4]),ylab='Delta',xlab='Given Teff of Benchmarks')

# Now let's look at a few plots to see if the bias functions are doing a good job - which clearly they are not...

coef.alphas <- summary(simple_test,vars=c('alpha1','alpha2','alpha3'))[,4]

plot(bench.param$TEFF[star.code[vector.of.setups == 1]],(observed.node.teff.spectrum[(vector.of.setups == 1),1]-bench.param$TEFF[star.code[vector.of.setups == 1]]),main=c('CAUP - U580'),xlab=c('Given Benchmark Teff (per spectrum)'),ylab=c('Delta'))
x <- bench.param$TEFF[star.code[vector.of.setups == 1]]
y <- coef.alphas[1]+coef.alphas[15]*x+coef.alphas[29]*x^2
lines(x,y,col='red',lwd=3)

plot(bench.param$TEFF[star.code[vector.of.setups == 2]],(observed.node.teff.spectrum[(vector.of.setups == 2),1]-bench.param$TEFF[star.code[vector.of.setups == 2]]),main=c('CAUP - U520'),xlab=c('Given Benchmark Teff (per spectrum)'),ylab=c('Delta'))
x <- bench.param$TEFF[star.code[vector.of.setups == 2]]
y <- coef.alphas[8]+coef.alphas[22]*x+coef.alphas[37]*x^2
lines(x,y,col='red',lwd=3)

plot(bench.param$TEFF[star.code[vector.of.setups == 1]],(observed.node.teff.spectrum[(vector.of.setups == 1),2]-bench.param$TEFF[star.code[vector.of.setups == 1]]),main=c('EPINARBO - U580'),xlab=c('Given Benchmark Teff (per spectrum)'),ylab=c('Delta'))
x <- bench.param$TEFF[star.code[vector.of.setups == 1]]
y <- coef.alphas[2]+coef.alphas[16]*x+coef.alphas[30]*x^2
lines(x,y,col='red',lwd=3)

plot(bench.param$TEFF[star.code[vector.of.setups == 2]],(observed.node.teff.spectrum[(vector.of.setups == 2),2]-bench.param$TEFF[star.code[vector.of.setups == 2]]),main=c('EPINARBO - U520'),xlab=c('Given Benchmark Teff (per spectrum)'),ylab=c('Delta'))
x <- bench.param$TEFF[star.code[vector.of.setups == 2]]
y <- coef.alphas[9]+coef.alphas[23]*x+coef.alphas[38]*x^2
lines(x,y,col='red',lwd=3)

plot(bench.param$TEFF[star.code[vector.of.setups == 1]],(observed.node.teff.spectrum[(vector.of.setups == 1),3]-bench.param$TEFF[star.code[vector.of.setups == 1]]),main=c('Lumba - U580'),xlab=c('Given Benchmark Teff (per spectrum)'),ylab=c('Delta'))
x <- bench.param$TEFF[star.code[vector.of.setups == 1]]
y <- coef.alphas[3]+coef.alphas[17]*x+coef.alphas[31]*x^2
lines(x,y,col='red',lwd=3)

plot(bench.param$TEFF[star.code[vector.of.setups == 2]],(observed.node.teff.spectrum[(vector.of.setups == 2),3]-bench.param$TEFF[star.code[vector.of.setups == 2]]),main=c('Lumba - U520'),xlab=c('Given Benchmark Teff (per spectrum)'),ylab=c('Delta'))
x <- bench.param$TEFF[star.code[vector.of.setups == 2]]
y <- coef.alphas[10]+coef.alphas[24]*x+coef.alphas[39]*x^2
lines(x,y,col='red',lwd=3)

par(def.par)#- reset to default

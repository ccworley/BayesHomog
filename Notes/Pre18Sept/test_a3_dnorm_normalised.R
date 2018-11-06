library(rjags);library(runjags);library(coda);library(jagstools);library(lattice);library(ggplot2);library(ggthemes)#;library(plot3D)
#library(gdata)
node.files.path <- c('/Users/charlotteworley/Documents/GES/WG10/iDR6ParameterHomog/FITSFiles/iDR6/FITSFiles/')

# Source the functions needed to read the Node and Benchmark files
source('load_input.R')

#Define node file paths etc
print(node.files.path)
release <- c('iDR6')
list.nodes <- c('Lumba','OACT','EPINARBO')  #,'IAC')  #,'MaxPlanck')

#Call load node files function and create common matrix (RS code)
#out.list <- load.nodes(list.nodes,release,node.files.path) # Reads the node files and outputs a list that is broken into 3 below:
nodes.param <- out.list[[2]]              # Array (nrow = num.stars, ncol = num.nodes, n.third = parameters [TEFF, E_TEFF, LOGG, E_LOGG, FEH, E_FEH, XI, E_XI] )
nodes.ids <- as.data.frame(out.list[[1]]) # The metadata associated to nodes.param (nrow = numstars, ncol = it varies; includes cnames, setup, snr, etc)
nodes.flags <-  out.list[[3]]             # Array (nrow = num.stars, ncol = num.nodes, n.third = flag) where flag = PECULI, REMARK, TECH in that order
#rm(out.list,release)

#Plot node teff
lteff <- nodes.param[,"Lumba","TEFF"]
oteff <- nodes.param[,"OACT","TEFF"]
eteff <- nodes.param[,"EPINARBO","TEFF"]
#icteff <- nodes.param[,"IAC","TEFF"]
plot(lteff,eteff,type="p")
plot(lteff,oteff,type="p")
#plot(lteff,icteff,type="p")

#Define Benchmark file paths etc
bench.path <- paste('/Users/charlotteworley/Documents/GES/WG10/iDR6ParameterHomog/Inputs/')
file.bench <- paste('GES_iDR6_FGKMCoolWarm_Benchmarks_AcceptedParams_11092018.fits')
columns.to.read <- c('GES_FLD','GES_TYPE','TEFF','E_TEFF','LOGG','E_LOGG','FEH','SIG1_MH0','SIG2_MH0','DELTA_ION','DELTA_LTE','XI') # Need the real errors of [Fe/H]; in the meantime using will be using the sqrt(delta_ion^2 + delta_lte^2)


# Loads benchmarks and adds in metallicity column
bench.param <- load.bench.fits(fitsname=file.bench,columns=columns.to.read,path.file=bench.path,only.fgk=TRUE)
errors.metallicity <- round(sqrt(bench.param$SIG1_MH0**2 +bench.param$SIG2_MH0**2 +bench.param$DELTA_ION**2 + bench.param$DELTA_LTE**2),digits=2) # Proper errors are not provided. Using delta_ion and delta_lte for now
errors.metallicity[errors.metallicity < 0.05] <- 0.05 # Because this is not the proper error, assume a lower limit of 0.05 dex
new.col.names <- c(colnames(bench.param),'E_FEH')
bench.param <- cbind(bench.param,errors.metallicity) # This is the matrix that holds the parameters of the benchmarks
colnames(bench.param) <- new.col.names
rm(errors.metallicity,new.col.names)  #removes these as they are now part of bench.param

#Checks number of benchmarks before filtering by those in iDR
print(dim(bench.param))

# To make sure that we use only stars from bench.param that are also included in the nodes.ids
# Removes those benchmarks that aren't in iDR
filter.bench <- (bench.param$GES_FLD %in% nodes.ids$GES_FLD)
bench.param <- bench.param[filter.bench,]
print(dim(bench.param))

# From the Node results, select all and only entries of the benchmarks
# Restrict to one setup for the moment
anasetup = "HR15N"
print(dim(nodes.ids))
#filter.bench <- (nodes.ids$GES_FLD %in% bench.param$GES_FLD)  #All setups
filter.bench <- (nodes.ids$GES_FLD %in% bench.param$GES_FLD & nodes.ids$SETUP %in% anasetup)
node.measured.param.bench <- nodes.param[filter.bench,,]
metadata.of.bench.spectra <- nodes.ids[filter.bench,]
#print(dim(metadata.of.bench.spectra))
#print(unique(metadata.of.bench.spectra$SETUP))

# Now, let's define all data structures that will be used in the model
# 1) Numbers
num.nodes   <- dim(nodes.param)[2]             # Number of Nodes
num.bench   <- nrow(bench.param)               # Number of individual benchmark stars
num.spectra <- nrow(metadata.of.bench.spectra) # Total number of spectra of benchmark stars

# 2) Series of Vector with the SNR of all benchmark spectra that were analyzed (in each row it repeats the SNR value of that spectrum by n times, where n = num.nodes)
#This is needed if we are going to parametize something against SNR
snr.spec <- as.numeric(as.vector(metadata.of.bench.spectra$SNR))
print(snr.spec)
print(max(snr.spec))  #See RS code if some SNR > 900, currently <695
print(min(snr.spec))
#Create matrix with duplicate snr vectors per node in node.param
snr.spec.vec <- matrix(rep(snr.spec,num.nodes),ncol=num.nodes,nrow=length(snr.spec))

# 3) Create the vectors with the given reference parameters
given.teff.bench <- bench.param$TEFF
given.sigma.teff.bench <- bench.param$E_TEFF

##Normalise the BM?
mean.bm.teff <- mean(given.teff.bench, na.rm=TRUE)
sd.bm.teff <- sd(given.teff.bench, na.rm=TRUE)

unnormalised.given.teff.bench <- given.teff.bench
unnormalised.given.sigma.teff.bench <- given.sigma.teff.bench

given.teff.bench <- (given.teff.bench-mean.bm.teff)/(sd.bm.teff)
given.sigma.teff.bench <- (given.sigma.teff.bench)/(sd.bm.teff)

# 4) The measurements from the nodes
observed.node.teff.spectrum <- matrix(NA,ncol=ncol(node.measured.param.bench),nrow=nrow(node.measured.param.bench))
for (ik in seq(1,ncol(observed.node.teff.spectrum))) {
  #print(as.numeric(as.vector(node.measured.param.bench[,ik,'TEFF'])))
  #print(trunc(as.numeric(as.vector(node.measured.param.bench[,ik,'TEFF']))))
  observed.node.teff.spectrum[,ik] <- trunc(as.numeric(as.vector(node.measured.param.bench[,ik,'TEFF']))) # To be sure they are numbers
}

##Normalise the node TEFFs by the mean and std of all the node TEFFS
#mean.all.teff <- mean(observed.node.teff.spectrum, na.rm=TRUE)
#sd.all.teff <- sd(observed.node.teff.spectrum, na.rm=TRUE)

unnormalised.observed.node.teff.spectrum <- observed.node.teff.spectrum
#observed.node.teff.spectrum <- (observed.node.teff.spectrum-mean.all.teff)/(sd.all.teff)
observed.node.teff.spectrum <- (observed.node.teff.spectrum-mean.bm.teff)/(sd.bm.teff)
# *** Need to make sure what is NaN in the unnormalised is NaN in normalised. Done below in when creating manipulated matrix


# Define star.code: it traces back to which benchmark a given line in observed.node.teff.spectrum 
# corresponds to, where the star.code is integer 1->Nbm correpsonding to list as in bench.param$GES_FLD
star.code <- vector('numeric',length=nrow(observed.node.teff.spectrum))
for (ik in seq(1,nrow(observed.node.teff.spectrum))) {
  star.code[ik] <- which(bench.param$GES_FLD == metadata.of.bench.spectra$GES_FLD[ik])
}

# Remove NA and NaN as per method in RS code
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
    if (is.na(unnormalised.observed.node.teff.spectrum[i,j])) {    #If NaN in observed
    #if (is.na(observed.node.teff.spectrum[i,j])) {    
      manipulated.node.teff.spectrum[i,j] <- (bench.param$TEFF[star.code[i]]-mean.bm.teff)/(sd.bm.teff) 
      observed.node.teff.spectrum[i,j] <- unnormalised.observed.node.teff.spectrum[i,j] #fills in NaN to match
      i.of.missing.teff[ik] <- i
      j.of.missing.teff[ik] <- j
      code.bench.missing[ik] <- star.code[i]
      ik <- ik+1
    } else {
      manipulated.node.teff.spectrum[i,j] <- observed.node.teff.spectrum[i,j]
      i.not.missing.teff[jk] <- i
      j.not.missing.teff[jk] <- j
      code.bench.not.missing[jk] <- star.code[i]
      jk <- jk+1
    }
  }
}

num.missing.values <- ik-1 # save the total number of missing entries
num.not.missing.values <- jk-1 # save the total number of non-missing entries

# 5) The Inverse Covariance Matrix  - see RS code
# From what I read in the internet, the recommended prior on the inverse covariance matrix 
#is a generic Wishart distribution, which is a multidimensional generalization of a gamma 
#distribution. W( V, n ) takes as input the degrees of freedom "n" and V = Sigma_0 * n^-1 
#where Sigma_0 is a first guess of the covariance matrix. The non-informative choice for 
#the degrees of freedom is to set it equal to the number of rows (in this case the number 
#of nodes).

prior.deg.freedom = num.nodes
prior.matrix = (1/num.nodes)*diag(100^2,ncol=num.nodes,nrow=num.nodes) 
#prior.matrix = (1/num.nodes)*diag(dgamma(0.1,0.1),ncol=num.nodes,nrow=num.nodes) 

# Assuming that the typical sigma of Teff = 100 K
# Some links:
#http://doingbayesiandataanalysis.blogspot.com/2017/06/bayesian-estimation-of-correlations-and.html
#http://www.themattsimpson.com/2012/08/20/prior-distributions-for-covariance-matrices-the-scaled-inverse-wishart-prior/
#https://stats.stackexchange.com/questions/66338/how-to-specify-the-wishart-distribution-scale-matrix
#https://en.wikipedia.org/wiki/Wishart_distribution

#6) We are going to fit one bias function per setup  - jsut HR15N at this point
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

sd.of.alpha1 = 0.1         #0.1        #400  1.0E-6  #
sd.of.alpha2 = 1     #1   1.0E-6  #
sd.of.alpha3 = 0.1     #1E-4


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
                  prior.deg.freedom = prior.deg.freedom,
                  sd.of.alpha1 = sd.of.alpha1,  #1      #400
                  sd.of.alpha2 = sd.of.alpha2,  #0.0003  #1
                  sd.of.alpha3 = sd.of.alpha3  #0.001  #1E-4
)

#Below defines a model that will be put into JAGS, hence the coding is not R and needs to be within double quotes.
#In JAGS, within the model, each relation (either ~ or <-) is called a node. When a node is defined in 
#terms of other nodes on the right it is a child node, where the nodes that define it are it's parent nodes.
#Nodes not defined in terms of other nodes are constant nodes.
#A node defined using ~ is stochastic (having a ran dom probability distribution or pattern that may be analysed 
#statistically but may not be predicted precisely) representing a random variable in the model.
#A node defined using <- is deterministic, the value of which is determined exactly by the values of its parents.
#(Mathematical model in which outcomes are precisely determined through known relationships among states 
#and events, without any room for random variation. In such models, a given input will always produce the 
#same output, such as in a known chemical reaction.)
#Nodes are embedded in named arrays e.g. MU=[mu[1],mu[2]...mu[N]]. Each element in the array (mu[n]) is a node in itself.

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

for (l in 1:n.missing) { # When the node was missing a value (n.missing), we include a value drawn from normal about BM mean and tau=1/sd^2  #the broad prior

true.teff.benchmark.per.spec[i.of.missing.values[l],j.of.missing.values[l]] ~ dnorm(5300,pow(800,-2))   #dunif(-1,1)  #dunif(3000,8000)   #dnorm(0,1) #when scaled by meanTbm and sigTbm  #

}

# Prior on the teff.InvCovMat
teff.InvCovMat ~ dwish( prior.matrix[1:K,1:K] , prior.deg.freedom )

# Other Priors needed for coefficients
# Priors on the coefficients alpha for the fit of the bias function. We fit one bias function per node. 
# The bias.teff is for now a quadratic function of the true.teff
for (l in 1:K) {
  for (m in 1:P) {

    # Remember here that dnorm(mu,tau) where mu = mean and tau = 1/sigma^2,
    # and sigma is the standard deviation 
    alpha1[l,m] ~ dnorm(0,pow(sd.of.alpha1,-2))   
    alpha2[l,m] ~ dnorm(0,pow(sd.of.alpha2,-2))
    alpha3[l,m] ~ dnorm(0,pow(sd.of.alpha3,-2))  #dunif(-10000, -0.01)    #
  
  }
}

# I will eventually need also a vector that really has only the the.true.teff.bench[i], without taking into account the missing values, 
# for the fitting of the bias function (see below)

for (i in 1:M) {
for (j in 1:K) {
only.the.true.teff.benchmark[i,j] <- the.true.teff.bench[star.code[i]]
}
}

# 2)
# Now, we finally write that what the nodes give is a noisy measurement of the true.teff.benchmark

#Likelihood - per spectrum (vectorized over K = num.nodes)
#Multivariate normal distribution where mu=Cor_Teff_(node)=Teff_(node)+bias_to_BM
#The bias is a polynomial order 2, where alpha are the coefficients

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
  list(alpha1 = matrix(rnorm(num.setups*num.nodes,0,sd.of.alpha1),ncol=num.setups,nrow=num.nodes),    # Here, outside JAGS, rnorm actually uses mu and sd, the standard deviation.
       alpha2 = matrix(rnorm(num.setups*num.nodes,0,sd.of.alpha2),ncol=num.setups,nrow=num.nodes),
       alpha3 = matrix(rnorm(num.setups*num.nodes,0,sd.of.alpha3),ncol=num.setups,nrow=num.nodes),     #matrix(runif(num.setups*num.nodes,-10000, -0.01),ncol=num.setups,nrow=num.nodes),   #
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


# Let's look at the correlation matrix

print(paste(' '))
print(paste('This is the final correlation matrix:'))
print(paste(' '))
estimated.cov.matrix <- summary(simple_test,vars=c('teff.Rho'))[,4]
estimated.cov.matrix <- matrix(round(estimated.cov.matrix,digits=2),num.nodes,num.nodes)
colnames(estimated.cov.matrix) <- list.nodes
rownames(estimated.cov.matrix) <- list.nodes
print(estimated.cov.matrix)
print(paste(' '))

# This will plot many figures, one at a time. Will need to flip between them to see them all.

# These are examples of plots to see the tracing and posterior
fit_mcmc <- as.mcmc(simple_test)
plot(xyplot(fit_mcmc[,0:5]))
plot(densityplot(fit_mcmc[,0:5]))

plot(xyplot(fit_mcmc[,5:10]))
plot(densityplot(fit_mcmc[,5:10]))

#Alphas
plot(densityplot(fit_mcmc[,37:45]))

# This will plot the Delta (given Teff of benchmarks minus estimated Teff of benchmarks) as a function Teff, logg and [Fe/H], so we can evaluate the trends that are left

results.teff <- summary(simple_test,vars=c('the.true.teff.bench'))
plot(bench.param$TEFF,(given.teff.bench-results.teff[,4]),ylab='Delta (Given-Inferred)',xlab='Given Teff of Benchmarks')

plot(bench.param$LOGG,(given.teff.bench-results.teff[,4]),ylab='Delta (Given-Inferred)',xlab='Given log g of Benchmarks')

plot(bench.param$FEH,(given.teff.bench-results.teff[,4]),ylab='Delta (Given-Inferred)',xlab='Given [Fe/H] of Benchmarks')

# Now let's look at a few plots to see if the bias functions are doing a good job - which clearly they are not...

#Extract coefficients
coef.alphas <- summary(simple_test,vars=c('alpha1','alpha2','alpha3'))[,4]

#Define vectors to do independent quadratic fit for comparison
#Basic order 2 fit to get a feel for coefficients
ii <- order(given.teff.bench[star.code[vector.of.setups == 1]])
xx <- given.teff.bench[star.code[vector.of.setups == 1]]
xx <- xx[ii]

plot(given.teff.bench[star.code[vector.of.setups == 1]],(observed.node.teff.spectrum[(vector.of.setups == 1),1]-given.teff.bench[star.code[vector.of.setups == 1]]),main=c('Lumba - HR15N'),xlab=c('Given Benchmark Teff (per spectrum)'),ylab=c('Delta'))
x <- sort(given.teff.bench[star.code[vector.of.setups == 1]])
y <- coef.alphas[1]+coef.alphas[4]*x+coef.alphas[7]*x^2  #update date which ones depending on how many nodes
#lines(x,y,col='red',lwd=3)
points(x,y,col='red',lwd=3)
lines(x,y,col='red',lwd=3)

#Fit for Lumba
yy <- (observed.node.teff.spectrum[(vector.of.setups == 1),1]-given.teff.bench[star.code[vector.of.setups == 1]])
yyL <- yy[ii]
fitpL <- lm( yyL~xx+I(xx^2))
xxx <- x
lines(xxx, predict(fitpL, data.frame(x=xxx)), col='green')
points(xxx, predict(fitpL, data.frame(x=xxx)), col='green')

plot(given.teff.bench[star.code[vector.of.setups == 1]],(observed.node.teff.spectrum[(vector.of.setups == 1),2]-given.teff.bench[star.code[vector.of.setups == 1]]),main=c('OACT - HR15N'),xlab=c('Given Benchmark Teff (per spectrum)'),ylab=c('Delta'))
x <- sort(given.teff.bench[star.code[vector.of.setups == 1]])
y <- coef.alphas[2]+coef.alphas[5]*x+coef.alphas[8]*x^2
points(x,y,col='red',lwd=3)
lines(x,y,col='red',lwd=3)
lines(x,-1*y,col='blue',lwd=3)

#Fit for OACT
yy <- (observed.node.teff.spectrum[(vector.of.setups == 1),2]-given.teff.bench[star.code[vector.of.setups == 1]])
yyO <- yy[ii]
fitpO <- lm( yyO~xx+I(xx^2))
xxx <- x
lines(xxx, predict(fitpO, data.frame(x=xxx)), col='green')
points(xxx, predict(fitpO, data.frame(x=xxx)), col='green')


plot(given.teff.bench[star.code[vector.of.setups == 1]],(observed.node.teff.spectrum[(vector.of.setups == 1),3]-given.teff.bench[star.code[vector.of.setups == 1]]),main=c('EPINARBO - HR15N'),xlab=c('Given Benchmark Teff (per spectrum)'),ylab=c('Delta'))
x <- sort(given.teff.bench[star.code[vector.of.setups == 1]])
y <- coef.alphas[3]+coef.alphas[6]*x+coef.alphas[9]*x^2
points(x,y,col='red',lwd=3)
lines(x,y,col='red',lwd=3)
#lines(x,-1*y,col='blue',lwd=3)

#Fit for EPINARBO
yy <- (observed.node.teff.spectrum[(vector.of.setups == 1),3]-given.teff.bench[star.code[vector.of.setups == 1]])
yyE <- yy[ii]
fitpE <- lm( yyE~xx+I(xx^2))
xxx <- x
lines(xxx, predict(fitpE, data.frame(x=xxx)), col='green')
points(xxx, predict(fitpE, data.frame(x=xxx)), col='green')

#plot(bench.param$TEFF[star.code[vector.of.setups == 1]],(observed.node.teff.spectrum[(vector.of.setups == 1),4]-given.teff.bench[star.code[vector.of.setups == 1]]),main=c('IAC - HR15N'),xlab=c('Given Benchmark Teff (per spectrum)'),ylab=c('Delta'))
#x <- sort(bench.param$TEFF[star.code[vector.of.setups == 1]])
#y <- coef.alphas[4]+coef.alphas[8]*x+coef.alphas[12]*x^2
#lines(x,y,col='red',lwd=3)

par(def.par)#- reset to default

library(rjags);library(runjags);library(coda);library(jagstools);library(lattice);library(ggplot2);library(ggthemes)#;library(plot3D)
#library(gdata)

# Original Code form RS updated to read in WG10
# 12/09/2018 Normalisation of values implemented
# 19/09/2018 New solution defining Bias in terms of accepted values only (See test_a3... version and readmes)



node.files.path <- c('/Users/charlotteworley/Documents/GES/WG10/iDR6ParameterHomog/FITSFiles/iDR6/FITSFiles/')

# Source the functions needed to read the Node and Benchmark files
source('load_input.R')

#Define node file paths etc
print(node.files.path)
release <- c('iDR6')
list.nodes.all <- c('Lumba','OACT','EPINARBO','IAC','MaxPlanck')

#Call load node files function and create common matrix (RS code)
#out.list <- load.nodes(list.nodes.all,release,node.files.path) # Reads the node files and outputs a list that is broken into 3 below:
nodes.param <- out.list[[2]]              # Array (nrow = num.stars, ncol = num.nodes, n.third = parameters [TEFF, E_TEFF, LOGG, E_LOGG, FEH, E_FEH, XI, E_XI] )
nodes.ids <- as.data.frame(out.list[[1]]) # The metadata associated to nodes.param (nrow = numstars, ncol = it varies; includes cnames, setup, snr, etc)
nodes.flags <-  out.list[[3]]             # Array (nrow = num.stars, ncol = num.nodes, n.third = flag) where flag = PECULI, REMARK, TECH in that order

# Restrict to one setup for the moment
anasetup = "HR15N"
if (length(grep("HR15N",anasetup)) > 0) {
  list.nodes <- c('Lumba','OACT','EPINARBO','MaxPlanck')
} else if (length(grep("HR9",anasetup)) > 0) {
  list.nodes <- c('OACT','EPINARBO')
} else if (length(grep("HR10|HR21",anasetup)) > 0) {
  list.nodes <- c('Lumba','IAC')
} else if (length(grep("HR21",anasetup)) > 0) {
  list.nodes <- c('Lumba','IAC','MaxPlanck')
} else if (length(grep("HR10",anasetup)) > 0) {
  list.nodes <- c('IAC','MaxPlanck')
}

#Reduce to nodes that anlaysed this setup
nodes.param <- nodes.param[,list.nodes,]             # Array (nrow = num.stars, ncol = num.nodes, n.third = parameters [TEFF, E_TEFF, LOGG, E_LOGG, FEH, E_FEH, XI, E_XI] )
nodes.flags <- nodes.flags[,list.nodes,]             # Array (nrow = num.stars, ncol = num.nodes, n.third = flag) where flag = PECULI, REMARK, TECH in that order

#Extract the rows that are benchmarks only
#filter.bench <- (nodes.param$GES_FLD %in% nodes.ids$GES_FLD)
#bench.param <- bench.param[filter.bench,]

#rm(out.list,release)

#Plot node teff
#lteff <- nodes.param[,"Lumba","TEFF"]
#oteff <- nodes.param[,"OACT","TEFF"]
#eteff <- nodes.param[,"EPINARBO","TEFF"]
#icteff <- nodes.param[,"IAC","TEFF"]
#plot(lteff,eteff,type="p")
#plot(lteff,oteff,type="p")
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
given.logg.bench <- bench.param$LOGG
given.sigma.logg.bench <- bench.param$E_LOGG
given.feh.bench <- bench.param$FEH
given.sigma.feh.bench <- bench.param$E_FEH

#Normalise the BM?
mean.bench.teff <- mean(given.teff.bench, na.rm=TRUE)
sd.bench.teff <- sd(given.teff.bench, na.rm=TRUE)
mean.bench.logg <- mean(given.logg.bench, na.rm=TRUE)
sd.bench.logg <- sd(given.logg.bench, na.rm=TRUE)
mean.bench.feh <- mean(given.feh.bench, na.rm=TRUE)
sd.bench.feh <- sd(given.feh.bench, na.rm=TRUE)

#unnormalised.given.teff.bench <- given.teff.bench
#unnormalised.given.sigma.teff.bench <- given.sigma.teff.bench
#given.teff.bench <- (given.teff.bench-mean.bm.teff)/(sd.bm.teff)
#given.sigma.teff.bench <- (given.sigma.teff.bench)/(sd.bm.teff)

# 4) The measurements from the nodes
observed.node.teff.spectrum <- matrix(NA,ncol=ncol(node.measured.param.bench),nrow=nrow(node.measured.param.bench))
observed.node.logg.spectrum <- matrix(NA,ncol=ncol(node.measured.param.bench),nrow=nrow(node.measured.param.bench))
observed.node.feh.spectrum <- matrix(NA,ncol=ncol(node.measured.param.bench),nrow=nrow(node.measured.param.bench))
for (ik in seq(1,ncol(observed.node.teff.spectrum))) {
  #print(as.numeric(as.vector(node.measured.param.bench[,ik,'TEFF'])))
  #print(trunc(as.numeric(as.vector(node.measured.param.bench[,ik,'TEFF']))))
  observed.node.teff.spectrum[,ik] <- trunc(as.numeric(as.vector(node.measured.param.bench[,ik,'TEFF']))) # To be sure they are numbers
  observed.node.logg.spectrum[,ik] <- trunc(as.numeric(as.vector(node.measured.param.bench[,ik,'LOGG']))*100)/100 # To be sure they are numbers
  #For Nice and IAC feh<-mh - code up!!!
  observed.node.feh.spectrum[,ik] <- trunc(as.numeric(as.vector(node.measured.param.bench[,ik,'FEH']))*100)/100 # To be sure they are numbers
}

##Normalise the node TEFFs by the mean and std of all the node TEFFS
#mean.all.teff <- mean(observed.node.teff.spectrum, na.rm=TRUE)
#sd.all.teff <- sd(observed.node.teff.spectrum, na.rm=TRUE)

#unnormalised.observed.node.teff.spectrum <- observed.node.teff.spectrum
#observed.node.teff.spectrum <- (observed.node.teff.spectrum-mean.all.teff)/(sd.all.teff)
#observed.node.teff.spectrum <- (observed.node.teff.spectrum-mean.bm.teff)/(sd.bm.teff)
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
manipulated.node.logg.spectrum <- matrix(NA,ncol=ncol(observed.node.teff.spectrum),nrow=nrow(observed.node.teff.spectrum))
manipulated.node.feh.spectrum <- matrix(NA,ncol=ncol(observed.node.teff.spectrum),nrow=nrow(observed.node.teff.spectrum))

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
    #if (is.na(unnormalised.observed.node.teff.spectrum[i,j])) {    #If NaN in observed
    if (is.na(observed.node.teff.spectrum[i,j])) {    #If NaN in observed
      manipulated.node.teff.spectrum[i,j] <- bench.param$TEFF[star.code[i]]#-mean.bm.teff)/(sd.bm.teff) 
      manipulated.node.logg.spectrum[i,j] <- bench.param$LOGG[star.code[i]]#-mean.bm.teff)/(sd.bm.teff) 
      manipulated.node.feh.spectrum[i,j] <- bench.param$FEH[star.code[i]]#-mean.bm.teff)/(sd.bm.teff) 
      #observed.node.teff.spectrum[i,j] <- unnormalised.observed.node.teff.spectrum[i,j] #fills in NaN to match
      i.of.missing.teff[ik] <- i
      j.of.missing.teff[ik] <- j
      code.bench.missing[ik] <- star.code[i]
      ik <- ik+1
      
    } else if (is.na(observed.node.logg.spectrum[i,j])) {    #If NaN in observed
        manipulated.node.teff.spectrum[i,j] <- bench.param$TEFF[star.code[i]]#-mean.bm.teff)/(sd.bm.teff) 
        manipulated.node.logg.spectrum[i,j] <- bench.param$LOGG[star.code[i]]#-mean.bm.teff)/(sd.bm.teff) 
        manipulated.node.feh.spectrum[i,j] <- bench.param$FEH[star.code[i]]#-mean.bm.teff)/(sd.bm.teff) 
        #observed.node.teff.spectrum[i,j] <- unnormalised.observed.node.teff.spectrum[i,j] #fills in NaN to match
        i.of.missing.teff[ik] <- i
        j.of.missing.teff[ik] <- j
        code.bench.missing[ik] <- star.code[i]
        ik <- ik+1

    } else if (is.na(observed.node.feh.spectrum[i,j])) {    #If NaN in observed
      manipulated.node.teff.spectrum[i,j] <- bench.param$TEFF[star.code[i]]#-mean.bm.teff)/(sd.bm.teff) 
      manipulated.node.logg.spectrum[i,j] <- bench.param$LOGG[star.code[i]]#-mean.bm.teff)/(sd.bm.teff) 
      manipulated.node.feh.spectrum[i,j] <- bench.param$FEH[star.code[i]]#-mean.bm.teff)/(sd.bm.teff) 
      #observed.node.teff.spectrum[i,j] <- unnormalised.observed.node.teff.spectrum[i,j] #fills in NaN to match
      i.of.missing.teff[ik] <- i
      j.of.missing.teff[ik] <- j
      code.bench.missing[ik] <- star.code[i]
      ik <- ik+1
      
    } else {
      manipulated.node.teff.spectrum[i,j] <- observed.node.teff.spectrum[i,j]
      manipulated.node.logg.spectrum[i,j] <- observed.node.logg.spectrum[i,j]
      manipulated.node.feh.spectrum[i,j] <- observed.node.feh.spectrum[i,j]
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

#sd.of.alpha1 = 0.1  #1      #400  1.0E-6  #
#sd.of.alpha2 = 0.1  #0.0003  #1   1.0E-6  #
#sd.of.alpha3 = 0.008  #0.001  #1E-4
ncoeff <- 3  #6


# Now, let's create the data list needed to run the model
teff.data <- list(given.teff.benchmarks = given.teff.bench,                       # Benchmark Teff
                  given.error.teff.benchmarks = given.sigma.teff.bench,          # Benchmark e_Teff
                  observed.node.teff.per.spec = manipulated.node.teff.spectrum,  # Node Teff
                  given.logg.benchmarks = given.logg.bench,                       # Benchmark Teff
                  given.error.logg.benchmarks = given.sigma.logg.bench,          # Benchmark e_Teff
                  observed.node.logg.per.spec = manipulated.node.logg.spectrum,  # Node Teff
                  given.feh.benchmarks = given.teff.bench,                       # Benchmark Teff
                  given.error.feh.benchmarks = given.sigma.feh.bench,          # Benchmark e_Teff
                  observed.node.feh.per.spec = manipulated.node.feh.spectrum,  # Node Teff
                  star.code = star.code,   # To trace to which benchmark each line of teff.by.node corresponds to
                  snr.spec.vec = snr.spec.vec,   # SNR
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
                  #For normalizing Teff values and estimating biases
                  offset.teff = mean.bench.teff,
                  scale.teff = sd.bench.teff,
                  offset.logg = mean.bench.logg,
                  scale.logg = sd.bench.logg,
                  offset.feh = mean.bench.feh,
                  scale.feh = sd.bench.feh,
                  ncoeff = ncoeff#,
                  #For wishart (dwish) prior on inverse covariance matrix:
                  #prior.matrix = prior.matrix,
                  #prior.deg.freedom = prior.deg.freedom  #,
                  #sd.of.alpha1 = sd.of.alpha1,  #1      #400
                  #sd.of.alpha2 = sd.of.alpha2#,  #0.0003  #1
                  #sd.of.alpha3 = sd.of.alpha3  #0.001  #1E-4
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

#teff.InvCovMat ~ dwish( prior.matrix[1:K,1:K] , prior.deg.freedom )

# Other Priors needed for coefficients
# Priors on the coefficients alpha for the fit of the bias function. We fit one bias function per node. The bias.teff is for now a quadratic function of the true.teff
#for (i in 1:2) { # For missing and non-missing values
for (l in 1:K) {       # On number of nodes
for (m in 1:P) {    # On number of setups
for (j in 1:3) { # 3 coefficients needed for fitting a quadratic function
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
}
}


# Prior for the cholesky factor matrix

cholesky.mat[1,1] <- 1                     # the c_11 term is always = 1
for (i in 2:K) {
  cholesky.mat[i,1] ~ dunif(-0.5,0.5)    # the c_i1 terms are random correlation coefficients
  cholesky.mat[1,i] <- 0
}

cholesky.mat[2,2] <- pow((1-pow(cholesky.mat[2,1],2)),1/2)  #c_22 is defined as a function of c_21

cholesky.mat[2,3] <- 0 # c_23 is always 0

upper.bound[3,2] <- cholesky.mat[3,1] * cholesky.mat[2,1] + pow(((pow(cholesky.mat[3,1],2)-1)*(pow(cholesky.mat[2,1],2)-1)),1/2)
lower.bound[3,2] <- cholesky.mat[3,1] * cholesky.mat[2,1] - pow(((pow(cholesky.mat[3,1],2)-1)*(pow(cholesky.mat[2,1],2)-1)),1/2)
prior.rho[3,2] ~ dunif(lower.bound[3,2],upper.bound[3,2])

cholesky.mat[3,2] <- pow(cholesky.mat[2,2],-1) * (prior.rho[3,2] - (cholesky.mat[3,1]*cholesky.mat[2,1])) # c_32 includes the rho_32 which is random but with bounds that depend on c_21 and c_31, plus it also depends on c_22

cholesky.mat[3,3] <- pow((1-(pow(cholesky.mat[3,1],2)+pow(cholesky.mat[3,2],2))),1/2) # c_33 is a function of c_31 and c_32


#Prior on the beta coefficients
for (i in 1:K) {
beta1[i] ~ dgamma(1,1E-6)
beta2[i] ~ dnorm(0,pow(100,-2)) I(0,)
beta3[i] ~ dnorm(0,pow(100,-2)) I(0,)
}


# 2)
# Now, we finally write that what the nodes give is a noisy measurement of the true.teff.benchmark

#Likelihood - per spectrum (vectorized over K = num.nodes)

for (j in 1:M) {
observed.node.teff.per.spec[j,1:K] ~ dmnorm( (true.teff.benchmark.per.spec[j,1:K] + (norm.bias.vector.teff.node[j,1:K]*scale.teff)) , teff.InvCovMat[j,1:K,1:K] )

norm.bias.vector.teff.node[j,1:K] <- alpha[1,1:K,setup.code[j]] + alpha[2,1:K,setup.code[j]]*only.the.true.teff.benchmark[j,1:K] + alpha[3,1:K,setup.code[j]]*pow(only.the.true.teff.benchmark[j,1:K],2)

node.var[j,1:K] <- beta1[1:K]*pow(exp(1),(pow(snr.spec.vec[j,1:K],-1)*beta2[1:K])) + beta3[1:K]
#node.var[j,1:K] <- beta1[1:K]*pow(snr.spec.vec[j,1:K],-1) + pow(beta2[1:K],2)
node.random[j,1:K] <- pow(node.var[j,1:K],1/2)

for (i in 1:K) {
for (l in 1:K) {
New.Mat[j,i,l] <- node.random[j,i]*cholesky.mat[i,l]
}
}

teff.CovMat[j,1:K,1:K] <- New.Mat[j,1:K,1:K] %*% t(New.Mat[j,1:K,1:K])
teff.InvCovMat[j,1:K,1:K] <- inverse(teff.CovMat[j,1:K,1:K])

}

# Notice what was done with the bias function above. Let's assume, for example, that I know that everytime I measure the Teff of the Sun, my value comes out with a bias of +50 K. I could write that my measurement minus the bias was drawn from the normal distribution: (my.measurement - 50 K) ~ dnorm(mu = 5777 K, sigma = my.random.error). What I am assuming is that this is equivalent to say that I am instead drawing from a normal distribution with mu = 5827 K (5777 + 50): my.measurement ~ dnorm(mu = 5827 K + my.bias , sigma = my.random.error)

# 3)
# This is extra. Since the quantities I know how to understand are the sigmas and correlations out of the Covariance Matrix:
# Convert teff.invCovMat to node.sd and correlations:

#  teff.CovMat <- inverse( teff.InvCovMat )
#  for ( varIdx in 1:K ) {
#      node.sd.teff[varIdx] <- sqrt(teff.CovMat[varIdx,varIdx])
#   }
#  for ( varIdx1 in 1:K ) {
#       for ( varIdx2 in 1:K ) {
#           teff.Rho[varIdx1,varIdx2] <- ( teff.CovMat[varIdx1,varIdx2] / (node.sd.teff[varIdx1]*node.sd.teff[varIdx2]) )
#  }
#}

}"

# Intial values for free parameters
# Intial values for free parameters
inits <- function() {
  list(alpha = array(rnorm(ncoeff*num.setups*num.nodes,0,0.1),dim=c(ncoeff,num.nodes,num.setups)),
       beta1 = rgamma(num.nodes,0.1,0.1),
       beta2 = abs(rnorm(num.nodes,0,100)),
       beta3 = abs(rnorm(num.nodes,0,100))
       #         beta2 = abs(rnorm(num.nodes,0,100))
       #teff.InvCovMat = diag(rgamma(num.nodes,0.1,0.1))
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
                        monitor=c("the.true.teff.bench","alpha","beta1","beta2","beta3"),
                        inits=list(inits1,inits2),
                        burnin=1000,
                        sample=2000,
                        summarise=F,
                        thin=1,
                        plots=F
)

# This outputs a summary of the MCMC of these parameters
print(summary(simple_test,vars=c('the.true.teff.bench')))
#print(summary(simple_test,vars=c('node.sd')))
#print(summary(simple_test,vars=c('teff.Rho')))
print(summary(simple_test,vars=c('alpha')))
print(summary(simple_test,vars=c('beta1')))
print(summary(simple_test,vars=c('beta2')))
print(summary(simple_test,vars=c('beta3')))


# Let's look at the correlation matrix

#print(paste(' '))
#print(paste('This is the final correlation matrix:'))
#print(paste(' '))
#estimated.cov.matrix <- summary(simple_test,vars=c('teff.Rho'))[,4]
#colnames(estimated.cov.matrix) <- list.nodes
#estimated.cov.matrix <- matrix(round(estimated.cov.matrix,digits=2),num.nodes,num.nodes)
#rownames(estimated.cov.matrix) <- list.nodes
#print(estimated.cov.matrix)
#print(paste(' '))

# This will plot many figures, one at a time. Will need to flip between them to see them all.

# These are examples of plots to see the tracing and posterior
fit_mcmc <- as.mcmc(simple_test)
plot(xyplot(fit_mcmc[,0:5]))  #plots per BM i=0-5
plot(densityplot(fit_mcmc[,0:5])) #plots per BM i=0-5

plot(xyplot(fit_mcmc[,5:10])) #plots per BM i=5-10
plot(densityplot(fit_mcmc[,5:10])) #plots per BM i=5-10

#Alphas
plot(densityplot(fit_mcmc[,length(fit_mcmc[1,])-8:length(fit_mcmc[1,])]))

# This will plot the Delta (given Teff of benchmarks minus estimated Teff of benchmarks) as a function Teff, logg and [Fe/H], so we can evaluate the trends that are left
par(mfrow=c(3,1))
results.teff <- summary(simple_test,vars=c('the.true.teff.bench'))
plot(bench.param$TEFF,(given.teff.bench-results.teff[,4]),ylab='Delta (Given-Inferred)',xlab='Given Teff of Benchmarks')

plot(bench.param$LOGG,(given.teff.bench-results.teff[,4]),ylab='Delta (Given-Inferred)',xlab='Given log g of Benchmarks')

plot(bench.param$FEH,(given.teff.bench-results.teff[,4]),ylab='Delta (Given-Inferred)',xlab='Given [Fe/H] of Benchmarks')

# Now let's look at a few plots to see if the bias functions are doing a good job - which clearly they are not...

#Extract coefficients
coef.alphas <- summary(simple_test,vars=c('alpha'))[,4]

#Define vectors to do independent quadratic fit for comparison
#Basic order 2 fit to get a feel for coefficients
it <- order(given.teff.bench[star.code[vector.of.setups == 1]])
x <- given.teff.bench[star.code[vector.of.setups == 1]]
xt <- x[it]  #teff in teff order

ig <- order(given.logg.bench[star.code[vector.of.setups == 1]])
xg <- given.logg.bench[star.code[vector.of.setups == 1]]
xg <- xg[ig]  #logg in logg order

#constant vectors for plotting
xts <- (given.teff.bench[star.code[vector.of.setups == 1]]-mean.bench.teff)/sd.bench.teff
xgs <- (given.logg.bench[star.code[vector.of.setups == 1]]-mean.bench.logg)/sd.bench.logg

xtst<- xts[it]  #scaled BM teff in teff order, if xgs is priority
xtsg<- xts[ig]  #scaled BM teff in logg order, if xgs is priority
xgst <- xgs[it]  #scaled BM logg in teff order, if xts is priority
xgsg <- xgs[ig]  #scaled BM logg in logg order, if xts is priority

#Iterate over nodes per setup
par(mfrow=c(1,2))
ind <- 0
for (each.node in list.nodes) {
  ind <- ind+1
  
  if (ncoeff == 3) {
    y <- (coef.alphas[1+(ind-1)*3]+coef.alphas[2+(ind-1)*3]*xts+coef.alphas[3+(ind-1)*3]*xts^2)*sd.bench.teff #update date which ones depending on how many nodes
  } else if (ncoeff == 6) {
    y <- (coef.alphas[1+(ind-1)*6]+coef.alphas[2+(ind-1)*6]*xts+coef.alphas[3+(ind-1)*6]*xts^2+coef.alphas[4+(ind-1)*6]*xgs+coef.alphas[5+(ind-1)*6]*xgs^2+coef.alphas[6+(ind-1)*6]*xts*xgs)*sd.bench.teff #update date which ones depending on how many nodes
  }

  #Plot delta teff against node teff
  plot(given.teff.bench[star.code[vector.of.setups == 1]],(observed.node.teff.spectrum[(vector.of.setups == 1),ind]-given.teff.bench[star.code[vector.of.setups == 1]]),main=c(paste(list.nodes[ind], ': HR15N', sep = "")),xlab=c('Given Benchmark Teff (per spectrum)'),ylab=c('Delta Teff'))
  #lines(x,y,col='red',lwd=3)
  points(xt,y[it],col='red',lwd=3)
  lines(xt,y[it],col='red',lwd=3)
  
  #Do quadratic fit with TEFF separately against TEFF
  yy <- (observed.node.teff.spectrum[(vector.of.setups == 1),ind]-given.teff.bench[star.code[vector.of.setups == 1]])
  yyn <- yy[it]
  fitpn <- lm( yyn~xt+I(xt^2))
  lines(xt, predict(fitpn, data.frame(x=xt)), col='green')
  points(xt, predict(fitpn, data.frame(x=xt)), col='green')
  
  #Plot delta teff against node logg
  plot(given.logg.bench[star.code[vector.of.setups == 1]],(observed.node.teff.spectrum[(vector.of.setups == 1),ind]-given.teff.bench[star.code[vector.of.setups == 1]]),main=c(paste(list.nodes[ind], ': HR15N', sep = "")),xlab=c('Given Benchmark log g (per spectrum)'),ylab=c('Delta Teff'))
  #y <- (coef.alphas[1]+coef.alphas[2]*x+coef.alphas[3]*x^2)*sd.bench.teff #update date which ones depending on how many nodes
  #y <- y[ig]  #(coef.alphas[1+(ind-1)*6]+coef.alphas[2+(ind-1)*6]*xtsg+coef.alphas[3+(ind-1)*6]*xtsg^2+coef.alphas[4+(ind-1)*6]*xgs+coef.alphas[5+(ind-1)*6]*xgs^2+coef.alphas[6+(ind-1)*6]*xtsg*xgs)*sd.bench.teff #update date which ones depending on how many nodes
  #lines(x,y,col='red',lwd=3)
  points(xg,y[ig],col='red',lwd=3)
  lines(xg,y[ig],col='red',lwd=3)
  
  #Do quadratic fit with TEFF separately against LOGG
  yy <- (observed.node.teff.spectrum[(vector.of.setups == 1),ind]-given.teff.bench[star.code[vector.of.setups == 1]])
  yyn <- yy[ig]
  fitpn <- lm( yyn~xg+I(xg^2))
  lines(xg, predict(fitpn, data.frame(x=xg)), col='green')
  points(xg, predict(fitpn, data.frame(x=xg)), col='green')
  

}

par(mfrow=c(1,1))

#par(def.par)#- reset to default

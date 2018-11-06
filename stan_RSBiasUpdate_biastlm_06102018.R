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
# Restrict to one setup for the moment
anasetup = "HR15N"
if (length(grep("HR15N",anasetup)) > 0) {
  list.nodes <- c('Lumba','OACT','EPINARBO','MaxPlanck')
} else if (length(grep("HR9",anasetup)) > 0) {
  list.nodes <- c('OACT','EPINARBO')
} else if (length(grep("HR10|HR21",anasetup)) > 0) {
  list.nodes <- c('Lumba','IAC')
} else if (length(grep("HR10",anasetup)) > 0) {
  list.nodes <- c('IAC','MaxPlanck')
}

#Call load node files function and create common matrix (RS code)
#out.list <- load.nodes(list.nodes,release,node.files.path) # Reads the node files and outputs a list that is broken into 3 below:
nodes.param <- out.list[[2]]              # Array (nrow = num.stars, ncol = num.nodes, n.third = parameters [TEFF, E_TEFF, LOGG, E_LOGG, FEH, E_FEH, XI, E_XI] )
nodes.ids <- as.data.frame(out.list[[1]]) # The metadata associated to nodes.param (nrow = numstars, ncol = it varies; includes cnames, setup, snr, etc)
nodes.flags <-  out.list[[3]]             # Array (nrow = num.stars, ncol = num.nodes, n.third = flag) where flag = PECULI, REMARK, TECH in that order

#Reduce to nodes that anlaysed this setup for TEFF only
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
print(unique(bench.param$GES_FLD))
bench.gesfld = unique(bench.param$GES_FLD)

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

spe.bench <- vector('numeric',length=num.bench)
for (ik in seq(1,length(spe.bench))) {
  print(length(which(metadata.of.bench.spectra$GES_FLD == bench.param$GES_FLD[ik])))
  spe.bench[ik] <- length(which(metadata.of.bench.spectra$GES_FLD == bench.param$GES_FLD[ik]))[1]
}
max.num.spe = max(spe.bench)

#Calibrator , Node, spectrum, >0 = missing is_missing[C, N, V]
#spectra per node per BM without results
#array(0, dim = c(N, D, L))  = N=list of matrices with D columns and L rows
mat.res.miss <- array(0,dim=c(num.bench,num.nodes,max.num.spe)) #is_missing[C, N, V], array indicating missing values
mat.res.val <- array(0,dim=c(num.bench,num.nodes,max.num.spe))  #matrix[N, V] estimates[C];array of values

mat_spe_isnr <- array(0,dim=c(max.num.spe,num.bench))  #matrix[V, C] spectrum_isnr

ik <- 1
jk <- 1
for (i in seq(1,num.bench)) {
  print(bench.param$GES_FLD[i])
  icr <- which(metadata.of.bench.spectra$GES_FLD == bench.param$GES_FLD[i])
  snr.bm = as.numeric(as.vector(metadata.of.bench.spectra$SNR[icr]))
  for (j in seq(1,num.nodes)) {
    print(list.nodes[j])
    node.bench.res = trunc(as.numeric(as.vector(node.measured.param.bench[icr,j,'TEFF'])))
    num.spe = length(node.bench.res)
    for (k in seq(1,num.spe)) {
      print(node.bench.res[k])
      if (is.na(node.bench.res[k])) {   #result is missing
        mat.res.miss[i,j,k] = 7
        ik = ik + 1 
        browser()
      } else {   #result is present
        mat.res.val[i,j,k] = node.bench.res[k] 
        jk = jk + 1
      }
      if (j==1) {
        mat_spe_isnr[k,i] = 1/snr.bm[k]
      } 
    #print(ik)
    #print(jk)
    }
  }
}

num.missing.values <- ik-1 # save the total number of missing entries
num.not.missing.values <- jk-1 # save the total number of non-missing entries

#all_mu_calibrator[C, S]  - BM=columns, parameters=rows
all.bm.param = cbind(given.teff.bench,given.logg.bench)  #,given.feh.bench)

#sd.of.alpha1 = 0.1  #1      #400  1.0E-6  #
#sd.of.alpha2 = 0.1  #0.0003  #1   1.0E-6  #
#sd.of.alpha3 = 0.008  #0.001  #1E-4
ncoeff <- 6  #6
nparam <- 2

#alpha_bounds = dict(teff=(100, 1000), logg=(0.1, 1.0), feh=(0.1, 1.0))

# Now, let's create the data list needed to run the model
stan.data <- list(N = num.nodes,  #N.nodes,
                  C = num.bench,  #N.benchmarks
                  V = max.num.spe,  #maximum number of visits to any calibrator (visits=spectra per benchmark)
                  visits = spe.bench, #number of visits to each calibrator
                  S = nparam, #the number of astrophysical parameters that can contribute to the node systematic variance (rank-ordered per param)                  

                  # The missing and not missing values data
                  TM = num.missing.values,  # total number of missing data points.
                  is_missing = mat.res.miss, # whether a node estimate is missing (>0 = missing)

                  #Node reported values (Nnodes,Nspectra)
                  estimates = mat.res.val,  # Node Teff
                  spectrum_isnr = mat_spe_isnr,  # inverse-snr of the spectrum (1/SNR)

                  #Accepted BM values - single parameter
                  mu_calibrator = given.teff.bench,                    # Benchmark param
                  sigma_calibrator = given.sigma.teff.bench,           # Benchmark e_param

                  #All accepted BM values for all parameters
                  all_mu_calibrator = all.bm.param,  # all non-spectroscopic calibrator values  available (for modeling systematic variance)

                  lower_alpha_sq = 100 ,  # lower bound on alpha parameter for f(snr)
                  upper_alpha_sq = 1000 ,  # upper bound on alpha parameter for f(snr)

                  lower_bound = 3000 ,  # lower univariate bound on the missing data values
                  upper_bound = 8000   # upper univariate bound on the missing data values

                  #For coefficients normalizing Teff values and estimating biases
                  #ncoeff = ncoeff,
                  
)


library(rstan)
fit1 <- stan(
  file = "ac_ensemble.stan",  # Stan program
  data = stan.data,    # named list of data
  chains = 2,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 1000,            # total number of iterations per chain
  cores = 1,              # number of cores (using 2 just for the vignette)
  refresh = 1          # show progress every 'refresh' iterations
)


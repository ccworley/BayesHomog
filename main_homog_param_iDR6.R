# R Smiljanic
# GES-WG11 iDR4 Homogenization of Atmospheric Parameters - modified for iDR5 - modified for iDR6
# Feb 2015 
# Modified Sep 2016 
# Jul-Oct 2018 - modified for iDR6
#
# Oct2018 - Modified to send to Clare WG10

# This is the main control of the homogenisation, including a lot of functions to make plots

def.par <- par(no.readonly = TRUE)

library(FITSio);library(stringr);library(gplots);library(robustbase);library(MASS)

# DEFINITIONS
release <- c('iDR6')
list.nodes <- c('Lumba','EPINARBO','IAC','OACT','MaxPlanck')
node.template.file <- c('/Users/charlotteworley/Documents/GES/WG10/iDR6ParameterHomog/FITSFiles/iDR6/FITSFiles/NodeTemplate/GES_iDR6_WG10_NodeTemplate_14062018_plus2ndExt.fits')
node.files.path <- c('/Users/charlotteworley/Documents/GES/WG10/iDR6ParameterHomog/FITSFiles/iDR6/FITSFiles/')

# Source the file with definitions of functions for plots and other tasks
source('homog_plots.R')

# Read the Node FITS files
# out.list <- load.nodes(list.nodes,release,node.files.path) # The list of columns to be read from the files has to be chnaged inside the function
# nodes.param <- out.list[[2]]
# nodes.ids <- as.data.frame(out.list[[1]])
# nodes.flags <-  out.list[[3]]
#rm(out.list)

#Data set of all setups and nodes
# nodes.param.all <- nodes.param
# nodes.flags.all <- nodes.flags
# nodes.ids.all <- nodes.ids

#Data set of all setups and nodes - save back to main variables as needed.
#Uncomment this BUT COMMENT the above!!!
nodes.param <- nodes.param.all
nodes.flags <- nodes.flags.all
nodes.ids <- nodes.ids.all


# Restrict to one setup for the moment
anasetup = "HR15N"   #"HR9B"
if ("HR15N" %in% anasetup) {
  list.nodes <- c('Lumba','OACT','EPINARBO')  #,'MaxPlanck')
} else if ("HR9B" %in% anasetup) {
  list.nodes <- c('OACT','EPINARBO')
} else if ("HR10|HR21" %in% anasetup) {
  list.nodes <- c('Lumba','IAC')
} else if ("HR21" %in%anasetup) {
  list.nodes <- c('Lumba','IAC','MaxPlanck')
} else if ("HR10" %in% anasetup) {
  list.nodes <- c('IAC','MaxPlanck')
}

#Filter sample to just the nodes and rows for this setup
#**** need to restrict this further for HR21 *****
if ("HR21" %in% anasetup) {
  #GES_TYPEs for HR21-only Milky Way Bulge plus Standards
  hr21gestypes <-c('AR_MW_BL','AR_SD_BM','AR_SD_GC','GE_MW_BL','GE_SD_BC','GE_SD_BM','GE_SD_BW','GE_SD_CR','GE_SD_GC','GE_SD_K2','GE_SD_MC','GE_SD_OC','GE_SD_PC','GE_SD_RV','GE_SD_TL')
  filter.setup <- (nodes.ids$SETUP %in% anasetup & nodes.ids$GES_TYPE %in% hr21gestypes)
  #print('here')
} else {
  filter.setup <- (nodes.ids$SETUP %in% anasetup)
  #print('heretoo')
}

nodes.param <- nodes.param[filter.setup,list.nodes,]             # Array (nrow = num.stars, ncol = num.nodes, n.third = parameters [TEFF, E_TEFF, LOGG, E_LOGG, FEH, E_FEH, XI, E_XI] )
nodes.flags <- nodes.flags[filter.setup,list.nodes,]             # Array (nrow = num.stars, ncol = num.nodes, n.third = flag) where flag = PECULI, REMARK, TECH in that order
nodes.ids <- nodes.ids[filter.setup,]

# Final list of all different flags used by the Nodes
# 
nodes.only.flags <- list.flags(nodes.flags)
final.list.flags <- sort(levels(factor(substr(unlist(nodes.only.flags),1,5))))
rm(nodes.only.flags)

# Based on the list produced above, and looking at the flags meaning from WG14 dictionary, below are the flags that mark results that should be ignored:
# List of flags to ignore for iDR6 - if these appear, all node results are excluded
flags.to.ignore.allnodes <- c("10005", "10106", "10108", "20020", "20030", "20040", "20070", "25000")

# 10005 SNR<5; 10106 broken spectrum (picket fence...); 10108 leaking from simcal causing spurious emission features; 
# 10320 suspected multiple system; 11100,11150,11200,11250 vsini>100,150,200,250, 13027 suspicious parameter multiple system, 200N0 SBN; 
# 25000 Halpha in emission
# NOTE: Confidence level is ignored

#
# For IDR6
# These flags, if found, exclude only the results of node that used the flag

flags.to.ignore.specific <- c("10015","10303","10304","10305","11100","11200","13022", # EPINARBO
                              "10015","10050","10302","10308","13020","13021","13022", # LUMBA
                              "10015","10020","10303","11100","11150","11200","13020","13022", # OACT
                              "10302", #IAC
                              "10020","10303","10305") # MaxPlanck


# EPINARBO: 10015 inaccurate result because of SNR; 10303, 10304, 10305, code convergence issues; 11100,11200 vsini issues; 13022, 13027 suspicious parameter
#      IAC: 10302 no convergence - no parameters given 
#    Lumba: 10015 10050 inaccurate result because of SNR; 10302 code convergence; 13020 13021 13022 suspicious stellar parameters
#     OACT: 10015,10020 inaccurate result because of SNR; 10303 code convergence; 11100,11150,11200 vsini>100,150,200; 13020,13022 suspicious parameter grid edge
#MaxPlanck: 10005,10010 SNR<5,10; 10302,10303 code convergence

# Remove the flagged results
list.flag.corrected <- eliminate.flagged.dr5(nodes.param,nodes.flags,flags.to.ignore.allnodes,flags.to.ignore.specific)
correc.flags.nodes.param <- list.flag.corrected[[1]]
print(sum(is.na(match(nodes.param[,,"TEFF"],"NaN"))))
#Check numbers when Flags are NaNned match AGs report - yes for HR15N
print(sum(is.na(match(correc.flags.nodes.param[,,"TEFF"],"NaN"))))
print(sum(is.na(match(correc.flags.nodes.param[,'Lumba',"TEFF"],"NaN"))))
print(sum(is.na(match(correc.flags.nodes.param[,'EPINARBO',"TEFF"],"NaN"))))
print(sum(is.na(match(correc.flags.nodes.param[,'OACT',"TEFF"],"NaN"))))

recom.flags <- list.flag.corrected[[2]]

# Remove weird results of IAC and MaxPlanck
correc.nodes.param <- correct.iacmp.grid(correc.flags.nodes.param)
#correc.nodes.param <- correct.iacaip.grid(correc.flags.nodes.param)

# Find and remove outliers
filter.outliers <- find.outliers(correc.nodes.param)
clean.nodes.param <- apply.filter.outliers(correc.nodes.param,filter.outliers)

# clean.nodes.param is the data that I want, so remove intermediate steps:
#rm(flags.to.ignore.allnodes,flags.to.ignore.specific,list.flag.corrected,correc.flags.nodes.param,correc.nodes.param,filter.outliers)

# Read the Benchmarks and the reference parameters
bench.path <- paste('/Users/charlotteworley/Documents/GES/WG10/iDR6ParameterHomog/Inputs/')
file.bench <- paste('GES_iDR6_FGKMCoolWarm_Benchmarks_AcceptedParams_11092018.fits')
#file.bench <- paste('GES_WG11xmatch_HR15N.fits')
columns.to.read <- c('GES_FLD','GES_TYPE','TEFF','E_TEFF','LOGG','E_LOGG','FEH','SIG1_MH0','SIG2_MH0','DELTA_ION','DELTA_LTE','XI') # Need the real errors of [Fe/H]; in the meantime using will be using the sqrt(delta_ion^2 + delta_lte^2)

# Load the parameters of the benchmark stars
bench.param <- load.bench.fits(fitsname=file.bench,columns=columns.to.read,path.file=bench.path,only.fgk=FALSE)
rm(bench.path,file.bench,columns.to.read)

# Do some plots
#sun.nodes(clean.nodes.param,nodes.ids,bench.param)                # Compare results for the Sun by all nodes
#bench.nodes(clean.nodes.param,nodes.ids,bench.param)              # Compare results for all benchmarks by all nodes
#compare.580and520.toref(clean.nodes.param,nodes.ids,bench.param)  # Compare results from U580 and U520
#plot.teff.logg.nodes(clean.nodes.param,nodes.ids)                 # Teff-logg diagram per node in 6 bins of FEH
#plot.clusters.nodes(clean.nodes.param,nodes.ids)                  # Teff-logg diagram of each cluster
#plot.metal.clusters.nodes(clean.nodes.param,nodes.ids)            # Check for correlations between FEH and TEFF or LOGG in clusters


##########################################
# FOR iDR6 - May 2018
# Compare the Node results from iDR6 with the recommended iDR5 results

#recom.param <- load.ges.dr5.data(list.columns.abun,list.columns.metadata,path.to.files='/work/rsmiljanic/Survey/DR5/Support/',abun.file='gesiDR5RecommendedAstroAnalysis.fits',metadata.file='results13_11_57_46_531.fits')

# Reorder the recom.param from iDR5 according to the order of CNAMEs in the iDR6 data (clean.nodes.param and nodes.ids)
# Saves this into test.recom.param

#test.recom.param <- as.data.frame(matrix(NA,ncol=5,nrow=length(nodes.ids$CNAME)))
#colnames(test.recom.param) <- c('CNAME','TEFF','LOGG','FEH','SNR')
#test.recom.param$CNAME <- nodes.ids$CNAME
#test.recom.param$SNR <- as.numeric(as.vector(nodes.ids$SNR))
#for (each.cname in test.recom.param$CNAME) {
#    ik <- which(test.recom.param$CNAME == each.cname)
#    jk <- which(recom.param$CNAME == each.cname)
#    if (length(jk) > 0) {
#        test.recom.param$TEFF[ik] <- recom.param$TEFF[jk]
#        test.recom.param$LOGG[ik] <- recom.param$LOGG[jk]
#        test.recom.param$FEH[ik] <- recom.param$FEH[jk]
#    }
#}

# Compute deltas and produce plots
#deltas.with.dr5 <- node.trends.with.recom(clean.nodes.param,test.recom.param,nodes.flags,release=substr(this.release,2,4),dir.to.save.plots='/work/rsmiljanic/Survey/DR6/Codes/Plots/')

# Rename data because from now on, the functions expect "nodes.param" but I want to use the "clean.nodes.param"

old.nodes.param <- nodes.param
nodes.param <- clean.nodes.param

########################################
#
# Now we start the bayesian approach to homogenise the results
#

# File that defines the functions to do the modelling
source('models_rjags_homog_param.R')

# Adds the column E_FEH, it is possible to choose more than one column
bench.param <- adds.feh.error(data.of.bench=bench.param,error.columns=c('SIG1_MH0')) 

# To make sure that we use only stars from bench.param that are also included in the nodes.ids
filter.bench <- (bench.param$GES_FLD %in% nodes.ids$GES_FLD)
bench.param <- bench.param[filter.bench,]

# From the Node results, select all and only entries of the benchmarks (for specific SETUP : already selected above)
filter.bench <- (nodes.ids$GES_FLD %in% bench.param$GES_FLD)
node.measured.param.bench <- nodes.param[filter.bench,,]
metadata.of.bench.spectra <- nodes.ids[filter.bench,]
rm(filter.bench)

# Now, let's define all data structures that will be used in the model
# 1) Numbers

num.nodes   <- dim(nodes.param)[2]             # Number of Nodes
num.bench   <- nrow(bench.param)               # Number of individual benchmark stars
num.spectra <- nrow(metadata.of.bench.spectra) # Total number of spectra of benchmark stars

# 2) Series of Vector with the SNR of all benchmark spectra that were analyzed 
# (in each row it repeats the SNR value of that spectrum by n times, where n = num.nodes)
# This is needed if we are going to parametize something against SNR
# There are some weird numbers of SNR > 500 up to ~25000! I am not sure if these should be believed.
# At this moment I am converting this very high values into something more "believable"; anything above 900 gets a random SNR between 900 and 1000 (truncated integer)

# Commented because parametrisation with SNR did not work so this is not needed
#Qant to check delta with snr
snr.spec.vec <- create.snr.vector(snr.data=metadata.of.bench.spectra$SNR)

# 3) Create the vectors with the given reference parameters

given.teff.bench <- bench.param$TEFF
given.sigma.teff.bench <- bench.param$E_TEFF
given.logg.bench <- bench.param$LOGG
given.sigma.logg.bench <-  bench.param$E_LOGG
given.feh.bench <- bench.param$FEH
given.sigma.feh.bench <-  bench.param$E_FEH

# 4) Separate the measurements from the nodes by TEFF, LOGG, FEH
all.observed <- node.measurements.of.reference.stars(data.of.nodes=node.measured.param.bench)
observed.node.teff.spectrum <- all.observed[[1]]
observed.node.logg.spectrum <- all.observed[[2]]
observed.node.feh.spectrum <- all.observed[[3]]
rm(all.observed)

# Define star.code: it traces back to which benchmark a given line in observed.node.teff/logg.spectrum corresponds to
star.code <- vector('numeric',length=nrow(observed.node.teff.spectrum))
for (ik in seq(1,nrow(observed.node.teff.spectrum))) {
  star.code[ik] <- which(bench.param$GES_FLD == metadata.of.bench.spectra$GES_FLD[ik])
}

# For FEH this is all different, because some "benchmarks" do not have FEH - they can be used for TEFF and LOGG but not for FEH
# most variables for FEH I am calling "corr.variables"

# WARNING: Inside the function below, I chose to assign Sigma_FEH = \pm 0.2 for those benchmarks without error in FEH

new.data.for.feh <- correct.data.for.feh(bench.data=bench.param,orig.star.code=star.code,orig.metadata=metadata.of.bench.spectra,observed.feh=observed.node.feh.spectrum,observed.snr=snr.spec.vec)
corr.bench.param <- new.data.for.feh[[1]]
corr.metadata.of.bench.spectra <- new.data.for.feh[[2]]
corr.observed.node.feh.spectrum <- new.data.for.feh[[3]]
corr.star.code <- new.data.for.feh[[4]]
corr.given.feh.bench <- new.data.for.feh[[5]]
corr.given.sigma.feh.bench <- new.data.for.feh[[6]]
positions.with.feh <- new.data.for.feh[[7]]
corr.snr.spec.vec  <- new.data.for.feh[[8]]
rm(new.data.for.feh)

# Create the manipulated data vectors that take care of missing values - First for Teff and logg

manipulated.data.teff.logg <- create.manipulated.teff.logg(observed.teff=observed.node.teff.spectrum,observed.logg=observed.node.logg.spectrum,bench.data=bench.param,the.stars=star.code)
manipulated.node.teff.spectrum <- manipulated.data.teff.logg[[1]]
manipulated.node.logg.spectrum <- manipulated.data.teff.logg[[2]]
i.of.missing.teff <- manipulated.data.teff.logg[[3]]
j.of.missing.teff <- manipulated.data.teff.logg[[4]]
i.not.missing.teff <- manipulated.data.teff.logg[[5]]
j.not.missing.teff <- manipulated.data.teff.logg[[6]]
code.bench.missing <- manipulated.data.teff.logg[[7]]
code.bench.not.missing <- manipulated.data.teff.logg[[8]]
keep.track.of.missing <- manipulated.data.teff.logg[[9]]
num.missing.values <- manipulated.data.teff.logg[[10]]
num.not.missing.values <- manipulated.data.teff.logg[[11]]
rm(manipulated.data.teff.logg)

# Create the manipulated data vectors that take care of missing values - Now for FEH - This is separated because some benchmarks do not have FEH values

manipulated.data.feh <- create.manipulated.feh(observed.feh=observed.node.feh.spectrum,bench.data=corr.bench.param,the.stars=corr.star.code,positions=positions.with.feh)
corr.manipulated.node.feh.spectrum <- manipulated.data.feh[[1]]
i.of.missing.feh <- manipulated.data.feh[[2]]
j.of.missing.feh <- manipulated.data.feh[[3]]
i.not.missing.feh <- manipulated.data.feh[[4]]
j.not.missing.feh <- manipulated.data.feh[[5]]
corr.code.bench.missing <- manipulated.data.feh[[6]]
corr.code.bench.not.missing <- manipulated.data.feh[[7]]
corr.keep.track.of.missing <- manipulated.data.feh[[8]]
corr.num.missing.values <- manipulated.data.feh[[9]]
corr.num.not.missing.values <- manipulated.data.feh[[10]]
rm(manipulated.data.feh)

# Prior for the convariance matrix

prior.deg.freedom = num.nodes
prior.matrix.teff <- create.prior.covariance(n.nodes=num.nodes,typical.sigma=100)
prior.matrix.logg <- create.prior.covariance(n.nodes=num.nodes,typical.sigma=0.2)
prior.matrix.feh <- create.prior.covariance(n.nodes=num.nodes,typical.sigma=0.2)

# Data on setups needed to separate the bias functions

data.on.setups <- create.setup.vector(column.setups=metadata.of.bench.spectra$SETUP)
vector.of.setups <- data.on.setups[[1]]
num.setups <- data.on.setups[[2]]
rm(data.on.setups)

# Re-work other variables for [Fe/H], because some benchmarks do not have reference metallicity values

corr.num.bench <- length(bench.param$FEH[!is.na(bench.param$FEH)])
corr.num.spectra <- nrow(corr.metadata.of.bench.spectra)
#corr.snr.spec.vec <- snr.spec.vec[positions.with.feh,]
corr.vector.of.setups <- vector.of.setups[positions.with.feh]

# Now, let's create the data list needed to run the model

teff.data <- list(given.teff.benchmarks = given.teff.bench,                      # Benchmark Teff
                  given.error.teff.benchmarks = given.sigma.teff.bench,          # Benchmark e_Teff
                  observed.node.teff.per.spec = manipulated.node.teff.spectrum,  # Node Teff
                  star.code = star.code,         # To trace to which benchmark each line of teff.by.node corresponds to
                  #snr.spec.vec = snr.spec.vec,   # SNR
                  K = num.nodes,                 #N.nodes,
                  N = num.bench,                 #N.benchmarks
                  M = num.spectra,               #N.spectra
                  P = num.setups,                #N.setups
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
                  prior.matrix = prior.matrix.teff,
                  prior.deg.freedom = prior.deg.freedom,
                  #For normalizing Teff values and estimating biases
                  offset.teff = mean(given.teff.bench),
                  scale.teff = sd(given.teff.bench)#,
                  #                  missing.code.for.bias = keep.track.of.missing # 1 is for not missing; 2 is for missing
)


logg.data <- list(given.logg.benchmarks = given.logg.bench,                      # Benchmark logg
                  given.error.logg.benchmarks = given.sigma.logg.bench,          # Benchmark e_logg
                  observed.node.logg.per.spec = manipulated.node.logg.spectrum,  # Node logg
                  star.code = star.code,         # To trace to which benchmark each line of logg.by.node corresponds to
                  #snr.spec.vec = snr.spec.vec,   # SNR
                  K = num.nodes,                 #N.nodes,
                  N = num.bench,                 #N.benchmarks
                  M = num.spectra,               #N.spectra
                  P = num.setups,                #N.setups
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
                  prior.matrix = prior.matrix.logg,
                  prior.deg.freedom = prior.deg.freedom,
                  #For normalizing logg values and estimating biases
                  offset.logg = mean(given.logg.bench),
                  scale.logg = sd(given.logg.bench)#,
                  #                  missing.code.for.bias = keep.track.of.missing # 1 is for not missing; 2 is for missing
)


feh.data <- list(given.feh.benchmarks = corr.given.feh.bench,                      # Benchmark feh
                 given.error.feh.benchmarks = corr.given.sigma.feh.bench,          # Benchmark e_feh
                 observed.node.feh.per.spec = corr.manipulated.node.feh.spectrum,  # Node feh
                 star.code = corr.star.code,         # To trace to which benchmark each line of feh.by.node corresponds to
                 #snr.spec.vec = corr.snr.spec.vec,   # SNR
                 K = num.nodes,                 #N.nodes,
                 N = corr.num.bench,                 #N.benchmarks
                 M = corr.num.spectra,               #N.spectra
                 P = num.setups,                #N.setups
                 setup.code = corr.vector.of.setups, # To trace to which setup the spectrum belongs                  
                 # The missing and not missing values data
                 n.missing = corr.num.missing.values, #N of missing values
                 n.not.missing = corr.num.not.missing.values, # N of NOT missing values
                 i.of.missing.values = i.of.missing.feh,
                 j.of.missing.values = j.of.missing.feh,
                 i.not.missing = i.not.missing.feh,
                 j.not.missing = j.not.missing.feh,
                 #                  code.bench.missing = code.bench.missing,
                 code.bench.not.missing = corr.code.bench.not.missing,
                 #For wishart (dwish) prior on inverse covariance matrix:
                 prior.matrix = prior.matrix.feh,
                 prior.deg.freedom = prior.deg.freedom,
                 #For normalizing feh values and estimating biases
                 offset.feh = mean(given.feh.bench,na.rm=T),
                 scale.feh = sd(given.feh.bench,na.rm=T)#,
                 #                  missing.code.for.bias = keep.track.of.missing # 1 is for not missing; 2 is for missing
)

# Sets of inital values, has to be the same number as the number of chains
num.chains.for.teff <- 8
num.chains.for.logg <- 8
num.chains.for.feh <- 8

list.inits.teff <- create.inits(variable=c('TEFF'),num.chains=num.chains.for.teff,n.setups=num.setups,n.nodes=num.nodes)
list.inits.logg <- create.inits(variable=c('LOGG'),num.chains=num.chains.for.logg,n.setups=num.setups,n.nodes=num.nodes)
list.inits.feh <- create.inits(variable=c('FEH'),num.chains=num.chains.for.feh,n.setups=num.setups,n.nodes=num.nodes)


# Here are the functions that actually run the MCMC chains to create the reference model based on the benchmarks

model.for.teff <- run.the.reference.model(input.data=teff.data,the.model=model.teff.matrix.bias,num.chains=num.chains.for.teff,num.adapt=1000,num.bur=1000,
                                          size.sample=1000,thin.factor=1,variables.to.monitor=c("the.true.teff.bench","alpha","teff.Rho","node.sd.teff"),
                                          model.inits=list.inits.teff)

model.for.logg <- run.the.reference.model(input.data=logg.data,the.model=model.logg.matrix.bias,num.chains=num.chains.for.logg,num.adapt=1000,num.bur=1000,
                                          size.sample=1000,thin.factor=1,variables.to.monitor=c("the.true.logg.bench","alpha","logg.Rho","node.sd.logg"),
                                          model.inits=list.inits.logg)

model.for.feh <- run.the.reference.model(input.data=feh.data,the.model=model.feh.matrix.bias,num.chains=num.chains.for.feh,num.adapt=1000,num.bur=1000,
                                          size.sample=1000,thin.factor=1,variables.to.monitor=c("the.true.feh.bench","alpha","feh.Rho","node.sd.feh"),
                                          model.inits=list.inits.feh)


#
# Series of plots and other diagnostics of the model
#

#show.correlation.matrix(model.for.teff,c('TEFF'))
#show.correlation.matrix(model.for.logg,c('LOGG'))
show.correlation.matrix(model.for.feh,c('FEH'))

#corr.diag.of.variable(model.for.teff,variable=c('alpha'))
#corr.diag.of.variable(model.for.logg,variable=c('alpha'))
corr.diag.of.variable(model.for.feh,variable=c('alpha'))

plot.against.reference.bench(model.for.teff,c('TEFF'),bench.param)
plot.against.reference.bench(model.for.logg,c('LOGG'),bench.param)
plot.against.reference.bench(model.for.feh,c('FEH'),corr.bench.param)


look.at.node.bias.functions(the.model=model.for.teff,variable=c('TEFF'),bench.data=bench.param,observed.data=observed.node.teff.spectrum,
                            the.setups=vector.of.setups,the.stars=star.code,col.of.setups=metadata.of.bench.spectra$SETUP,
                            nodes=list.nodes,mean.param.bench=mean(given.teff.bench,na.rm=T),sd.param.bench=sd(given.teff.bench,na.rm=T),
                            observed.snr=snr.spec.vec)

look.at.node.bias.functions(the.model=model.for.logg,variable=c('LOGG'),bench.data=bench.param,observed.data=observed.node.logg.spectrum,
                            the.setups=vector.of.setups,the.stars=star.code,col.of.setups=metadata.of.bench.spectra$SETUP,
                            nodes=list.nodes,mean.param.bench=mean(given.logg.bench,na.rm=T),sd.param.bench=sd(given.logg.bench,na.rm=T),
                            observed.snr=snr.spec.vec)

look.at.node.bias.functions(the.model=model.for.feh,variable=c('FEH'),bench.data=corr.bench.param,observed.data=corr.observed.node.feh.spectrum,
                            the.setups=corr.vector.of.setups,the.stars=corr.star.code,col.of.setups=corr.metadata.of.bench.spectra$SETUP,
                            nodes=list.nodes,mean.param.bench=mean(given.feh.bench,na.rm=T),sd.param.bench=sd(given.feh.bench,na.rm=T),
                            observed.snr=corr.snr.spec.vec)
  
stop('Reference Models done')
# Save the models for future use

#save(model.for.teff,file='my_reference_teff_model.Rdata')
#save(model.for.logg,file='my_reference_logg_model.Rdata')
#save(model.for.feh,file='my_reference_feh_model.Rdata')

# Load back the needed model

#load(file = 'my_reference_teff_model.Rdata')
#load(file = 'my_reference_logg_model.Rdata')
#load(file = 'my_reference_feh_model.Rdata')

#########################
# Now to homogenise the whole sample - CNAME by CNAME

data.for.full.homog.teff <- prepare.for.full.sample.homog(variable='TEFF',this.model=model.for.teff,full.metadata=nodes.ids,all.param=clean.nodes.param,bench.data=bench.param)
data.for.full.homog.logg <- prepare.for.full.sample.homog(variable='LOGG',this.model=model.for.logg,full.metadata=nodes.ids,all.param=clean.nodes.param,bench.data=bench.param)
data.for.full.homog.feh <- prepare.for.full.sample.homog(variable='FEH',this.model=model.for.feh,full.metadata=nodes.ids,all.param=clean.nodes.param,bench.data=bench.param)

homog.teff.data <- list(
  observed.teff = data.for.full.homog.teff[[1]],         # All Teff measurements (nrows = spectra, ncol = nodes)
  N = data.for.full.homog.teff[[2]],                     # Total number of individual CNAMEs
  M = data.for.full.homog.teff[[3]],                     # Total number of individual spectra
  L = length(list.nodes),                                # Number of nodes
  star.code.not.missing = data.for.full.homog.teff[[4]], # star.code of non-missing values
  spectrum.not.missing = data.for.full.homog.teff[[5]],  # row number of observed.teff where there are real measurements (used together with node.not.missing)
  node.not.missing = data.for.full.homog.teff[[6]],      # col number of ovserved.teff where there are real measurements (used together with spectrum.not.missing)
  tot.num.not.missing = data.for.full.homog.teff[[7]],   # Total number of real measurements
  tot.missing = data.for.full.homog.teff[[8]],           # Total number of missing measurements
  spectrum.missing = data.for.full.homog.teff[[9]],      # row number of observed.teff where there are missing measurements (used together with node.missing)
  node.missing = data.for.full.homog.teff[[10]],         # col number of ovserved.teff where there are missing measurements (used together with spectrum.missing)
  mean.bias.teff = data.for.full.homog.teff[[11]],       # Mean value of the bias at the given Teff of the measurement
  sd.bias.teff = data.for.full.homog.teff[[12]],         # SD value  of the bias at the given Teff of the measurement
  mean.node.sd = data.for.full.homog.teff[[13]],         # The node typical SD in Teff 
  mean.teff.Rho = data.for.full.homog.teff[[14]]         # The correlation coefficients in Teff 
)

homog.logg.data <- list(
  observed.logg = data.for.full.homog.logg[[1]],         # All logg measurements (nrows = spectra, ncol = nodes)
  N = data.for.full.homog.logg[[2]],                     # Total number of individual CNAMEs
  M = data.for.full.homog.logg[[3]],                     # Total number of individual spectra
  L = length(list.nodes),                                # Number of nodes
  star.code.not.missing = data.for.full.homog.logg[[4]], # star.code of non-missing values
  spectrum.not.missing = data.for.full.homog.logg[[5]],  # row number of observed.logg where there are real measurements (used together with node.not.missing)
  node.not.missing = data.for.full.homog.logg[[6]],      # col number of ovserved.logg where there are real measurements (used together with spectrum.not.missing)
  tot.num.not.missing = data.for.full.homog.logg[[7]],   # Total number of real measurements
  tot.missing = data.for.full.homog.logg[[8]],           # Total number of missing measurements
  spectrum.missing = data.for.full.homog.logg[[9]],      # row number of observed.logg where there are missing measurements (used together with node.missing)
  node.missing = data.for.full.homog.logg[[10]],         # col number of ovserved.logg where there are missing measurements (used together with spectrum.missing)
  mean.bias.logg = data.for.full.homog.logg[[11]],       # Mean value of the bias at the given logg of the measurement
  sd.bias.logg = data.for.full.homog.logg[[12]],         # SD value  of the bias at the given logg of the measurement
  mean.node.sd = data.for.full.homog.logg[[13]],         # The node typical SD in logg 
  mean.logg.Rho = data.for.full.homog.logg[[14]]         # The correlation coefficients in logg 
)

homog.feh.data <- list(
  observed.feh = data.for.full.homog.feh[[1]],          # All feh measurements (nrows = spectra, ncol = nodes)
  N = data.for.full.homog.feh[[2]],                     # Total number of individual CNAMEs
  M = data.for.full.homog.feh[[3]],                     # Total number of individual spectra
  L = length(list.nodes),                               # Number of nodes
  star.code.not.missing = data.for.full.homog.feh[[4]], # star.code of non-missing values
  spectrum.not.missing = data.for.full.homog.feh[[5]],  # row number of observed.feh where there are real measurements (used together with node.not.missing)
  node.not.missing = data.for.full.homog.feh[[6]],      # col number of ovserved.feh where there are real measurements (used together with spectrum.not.missing)
  tot.num.not.missing = data.for.full.homog.feh[[7]],   # Total number of real measurements
  tot.missing = data.for.full.homog.feh[[8]],           # Total number of missing measurements
  spectrum.missing = data.for.full.homog.feh[[9]],      # row number of observed.feh where there are missing measurements (used together with node.missing)
  node.missing = data.for.full.homog.feh[[10]],         # col number of ovserved.feh where there are missing measurements (used together with spectrum.missing)
  mean.bias.feh = data.for.full.homog.feh[[11]],        # Mean value of the bias at the given feh of the measurement
  sd.bias.feh = data.for.full.homog.feh[[12]],          # SD value  of the bias at the given feh of the measurement
  mean.node.sd = data.for.full.homog.feh[[13]],         # The node typical SD in feh 
  mean.feh.Rho = data.for.full.homog.feh[[14]]          # The correlation coefficients in feh 
)

# Number of chains and initial values 

num.chains.for.teff <- 4
num.chains.for.logg <- 4
num.chains.for.feh  <- 4

list.inits.full.teff <- inits.for.full.homog(variable=c('TEFF'),num.chains=num.chains.for.teff,tot.num.stars=data.for.full.homog.teff[[2]],
                                             min.param=3000,max.param=8000)

list.inits.full.logg <- inits.for.full.homog(variable=c('LOGG'),num.chains=num.chains.for.logg,tot.num.stars=data.for.full.homog.logg[[2]],
                                             min.param=0.0,max.param=5.5)

list.inits.full.feh <- inits.for.full.homog(variable=c('FEH'),num.chains=num.chains.for.feh,tot.num.stars=data.for.full.homog.feh[[2]],
                                             min.param=-3.5,max.param=0.5)


# Run the full homogenisation per parameter

full.teff.homog <- run.the.full.homog(input.data=homog.teff.data,the.model=homog.all.teff,num.chains=num.chains.for.teff,num.adapt=1000,num.bur=1000,
                                         size.sample=1000,thin.factor=1,variables.to.monitor=c("true.teff"),
                                         model.inits=list.inits.full.teff)

full.logg.homog <- run.the.full.homog(input.data=homog.logg.data,the.model=homog.all.logg,num.chains=num.chains.for.logg,num.adapt=1000,num.bur=1000,
                                      size.sample=1000,thin.factor=1,variables.to.monitor=c("true.logg"),
                                      model.inits=list.inits.full.logg)

full.feh.homog <- run.the.full.homog(input.data=homog.feh.data,the.model=homog.all.feh,num.chains=num.chains.for.feh,num.adapt=100,num.bur=100,
                                      size.sample=100,thin.factor=1,variables.to.monitor=c("true.feh"),
                                      model.inits=list.inits.full.feh)


# Save the results for use later
#save(full.teff.homog,file='results_homog_teff_full_sample.Rdata')
#save(full.teff.homog,file='results_homog_logg_full_sample.Rdata')
#save(full.teff.homog,file='results_homog_feh_full_sample.Rdata')

# Load previous results
#load(file='results_homog_teff_full_sample.Rdata')
#load(file='results_homog_logg_full_sample.Rdata')
#load(file='results_homog_feh_full_sample.Rdata')

# Needs all parameters (TEFF,LOGG,FEH)
all.data.homog <- create.data.table(homog.teff=full.teff.homog,homog.logg=full.logg.homog,homog.feh=full.feh.homog,metadata.here=data.for.full.homog.feh)
plot.the.final.benchmarks(data.homog=all.data.homog,bench.data=bench.param)



par(def.par)#- reset to default
#



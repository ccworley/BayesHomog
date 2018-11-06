library(rjags);library(runjags);library(coda);library(jagstools);library(lattice);library(ggplot2);library(ggthemes)#;library(plot3D)
#Looking at dependence on BM parameters part of ACs code

# Source the functions needed to read the Node and Benchmark files
source('load_input.R')

#ACs code
# C=number of BMs
# S=number of parameters
# N= number of nodes
# transformed data {
#   int TCV; // total number of calibrator visits
#   real AMCS[C, S]; // all_mu_calibrator values scaled between [0, 1]
#   
#   TCV = sum(visits);
#   for (s in 1:S) {
#     real offset;
#     real scale;
#     
#     offset = min(all_mu_calibrator[:, s]);
#     scale = max(all_mu_calibrator[:, s]) - offset;
#     
#     for (c in 1:C)
#       AMCS[c, s] = (all_mu_calibrator[c, s] - offset)/scale;
#   }
# }
# 
# for (s in 1:S)
#   talk[s] = vs_a[s, n] * pow(1 - AMCS[c, s], vs_b[s, n]);
# 
# sigma[n] = sqrt(alpha_sq[n] * spectrum_isnr[v, c] +
#                   vs_c[n] * exp(sum(talk)));


#Per parameter find the offset and the scale
teffoff = min(bench.param$TEFF)
teffscale = max(bench.param$TEFF) - teffoff
amcsteff = (bench.param$TEFF - teffoff)/teffscale
loggoff = min(bench.param$LOGG)
loggscale = max(bench.param$LOGG) - loggoff
amcslogg = (bench.param$LOGG - loggoff)/loggscale
ifehfin = is.finite(bench.param$FEH)
fehoff = min(bench.param$FEH[ifehfin])
fehscale = max(bench.param$FEH[ifehfin]) - fehoff
amcsfeh = (bench.param$FEH - fehoff)/fehscale

pwr = 2
par(mfrow = c(4, 4))
plot(bench.param$TEFF,bench.param$TEFF - teffoff,ylab='BM-min(BM)',xlab='Given Teff of Benchmarks')
plot(bench.param$TEFF,amcsteff,ylab='AMCS = (BM-minBM)/(maxBM-minBM)',xlab='Given Teff of Benchmarks')
plot(bench.param$TEFF,1-amcsteff,ylab='1-AMCS',xlab='Given Teff of Benchmarks')
plot(bench.param$TEFF,(1-amcsteff)^pwr,ylab='(1-AMCS)^2',xlab='Given Teff of Benchmarks')
plot(bench.param$LOGG,bench.param$LOGG - loggoff,ylab='BM-min(BM)',xlab='Given LOGG of Benchmarks')
plot(bench.param$LOGG,amcslogg,ylab='AMCS = (BM-minBM)/(maxBM-minBM)',xlab='Given LOGG of Benchmarks')
plot(bench.param$LOGG,1-amcslogg,ylab='1-AMCS',xlab='Given LOGG of Benchmarks')
plot(bench.param$LOGG,(1-amcslogg)^pwr,ylab='(1-AMCS)^2',xlab='Given LOGG of Benchmarks')
plot(bench.param$FEH,bench.param$FEH - fehoff,ylab='BM-min(BM)',xlab='Given FEH of Benchmarks')
plot(bench.param$FEH,amcsfeh,ylab='AMCS = (BM-minBM)/(maxBM-minBM)',xlab='Given FEH of Benchmarks')
plot(bench.param$FEH,1-amcsfeh,ylab='1-AMCS',xlab='Given FEH of Benchmarks')
plot(bench.param$FEH,(1-amcsfeh)^pwr,ylab='(1-AMCS)^2',xlab='Given FEH of Benchmarks')

pwteff = (1-amcsteff)^pwr
pwlogg = (1-amcslogg)^pwr
pwfeh = (1-amcsfeh)^pwr

talkbm = rbind(pwteff,pwlogg,pwfeh)
talkbmsum = colSums(talkbm,na.rm=TRUE)
talkbmsumexp = exp(colSums(talkbm,na.rm=TRUE))

plot(bench.param$TEFF,pwteff,ylab='(1-AMCS)^2',xlab='Given Teff of Benchmarks')
points(bench.param$TEFF,pwlogg,col='red',lwd=3)
points(bench.param$TEFF,pwfeh,col='green',lwd=3)

plot(bench.param$TEFF,talkbmsum,ylab='sum((1-AMCS)^2)',xlab='Given Teff of Benchmarks')
plot(bench.param$TEFF,talkbmsumexp,ylab='exp(sum((1-AMCS)^2))',xlab='Given Teff of Benchmarks')

scatter3D(x = bench.param$TEFF, y = bench.param$LOGG, z = amcs, ticktype = "detailed", pch = 16, bty = "f",theta=45, phi=45,
          colvar=amcs, xlab = "TEFF", ylab = "LOGG", zlab = "Weight")

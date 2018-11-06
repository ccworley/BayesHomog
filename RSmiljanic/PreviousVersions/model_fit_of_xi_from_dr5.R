library(graphics)
library(gplots)
library(FITSio)
library(stringr)
#library(mixtools)
#library(mclust)
#library(celestial)
#
library(rjags);library(runjags);library(coda);library(jagstools);library(lattice)#;library(ggplot2);library(ggthemes);library(plot3D)

#
def.par <- par(no.readonly = TRUE)
#
wg15.dr5 <- readFITS(file="/Users/rodolfo/Survey/Homogen/iDR5/gesiDR5RecommendedAstroAnalysis.fits", hdu=1, maxLines=5500)
num.stars <- length(wg15.dr5$col[[2]])
list.of.columns <- c("RECWG","WG","WGSOURCE","ISWGPARAMS","GESTYPE","INSTRUMENT","RECGRATINGS","GRATINGS","GESFIELD","GESOBJECT",
                     "CNAME","TEFF","TEFFERR","NNTEFF","ENNTEFF","NNETEFF","LOGG","LOGGERR","NNLOGG","ENNLOGG","NNELOGG",
                     "LIMLOGG","FEH","FEHERR","NNFEH","ENNFEH","NNEFEH","XI","XIERR","NNXI","ENNXI","NNEXI","OBJRA","OBJDEC","SNR",
                     "RADVEL","RADVELERR","ROTVEL","ROTVELERR","VRAD","VRADERR","VRADPROV","VRADOFFSET","VRADFILENAME","VSINI",
                     "VSINIERR","VSINILIM","PECULI","REMARK","TECH")
#
extr.wg15.dr5 <- matrix(' ',num.stars,length(list.of.columns))
colnames(extr.wg15.dr5) <- list.of.columns
extr.wg15.dr5 <- as.data.frame(extr.wg15.dr5)
#
for (ik in seq(1,length(list.of.columns),1)) {
    numb.col <- which(wg15.dr5$colNames == list.of.columns[ik])
    extr.wg15.dr5[,ik] <- str_trim(wg15.dr5$col[[numb.col]])
}
numeric.columns <- c("ISWGPARAMS","TEFF","TEFFERR","NNTEFF","ENNTEFF","NNETEFF","LOGG","LOGGERR","NNLOGG","ENNLOGG","NNELOGG",
                     "LIMLOGG","FEH","FEHERR","NNFEH","ENNFEH","NNEFEH","XI","XIERR","NNXI","ENNXI","NNEXI","OBJRA","OBJDEC","SNR",
                     "RADVEL","RADVELERR","ROTVEL","ROTVELERR","VRAD","VRADERR","VRADOFFSET","VSINI",
                     "VSINIERR","VSINILIM")
for (each.col in numeric.columns) {
    this.col <- which(colnames(extr.wg15.dr5) == each.col)
    extr.wg15.dr5[,this.col] <- as.vector(as.numeric(extr.wg15.dr5[,this.col]))
    extr.wg15.dr5[(extr.wg15.dr5[,this.col] < -90),this.col] <- NA
}
ges.dr5 <- cbind(extr.wg15.dr5)
ges.dr5$TECH <- as.vector(as.character(ges.dr5$TECH))
ges.dr5$PECULI <- as.vector(as.character(ges.dr5$PECULI))
ges.dr5$REMARK <- as.vector(as.character(ges.dr5$REMARK))
ges.dr5$RECWG <- as.vector(as.character(ges.dr5$RECWG))
rm(wg15.dr5,extr.wg15.dr5)
#
filter.param <- (ges.dr5$RECWG %in% c('WG11','WG12')) & (ges.dr5$INSTRUMENT == 'UVES') & !is.na(ges.dr5$TEFF) & !is.na(ges.dr5$LOGG) & !is.na(ges.dr5$FEH) & !is.na(ges.dr5$XI) & !is.na(ges.dr5$VSINI) & (((ges.dr5$SNR >= 40) & (ges.dr5$VSINI <= 5) & (ges.dr5$LOGG < 3.5)) | ((ges.dr5$SNR >= 40) & (ges.dr5$VSINI <= 5) & (ges.dr5$TEFF > 5250) & (ges.dr5$LOGG >= 3.5)) | ((ges.dr5$SNR >= 20) & (ges.dr5$VSINI <= 10) & (ges.dr5$TEFF <= 5250) & (ges.dr5$LOGG >= 3.5))) & (ges.dr5$PECULI == '') & !grepl("10306",ges.dr5$TECH) & !grepl("10315",ges.dr5$TECH) & !grepl("10391",ges.dr5$TECH) & !grepl("12007",ges.dr5$TECH) & !grepl("13000",ges.dr5$TECH) & !grepl("13023",ges.dr5$TECH) & !grepl("10106",ges.dr5$TECH) & !grepl("10110",ges.dr5$TECH) & !grepl("10103",ges.dr5$TECH) & !grepl("10104",ges.dr5$TECH) & !grepl("10309",ges.dr5$TECH) & !grepl("12004",ges.dr5$TECH) & !grepl("13020",ges.dr5$TECH) & !grepl("13021",ges.dr5$TECH) & !grepl("13022",ges.dr5$TECH) 

#filter.param2 <- (ges.dr5$RECWG %in% c('WG11','WG12')) & (ges.dr5$INSTRUMENT == 'UVES') & !is.na(ges.dr5$TEFF) & !is.na(ges.dr5$LOGG) & !is.na(ges.dr5$FEH) & !is.na(ges.dr5$XI) & (ges.dr5$SNR >= 40) & !is.na(ges.dr5$VSINI) & (ges.dr5$VSINI <= 5) & (ges.dr5$PECULI == '') & !grepl("10306",ges.dr5$TECH) & !grepl("10315",ges.dr5$TECH) & !grepl("10391",ges.dr5$TECH) & !grepl("12007",ges.dr5$TECH) & !grepl("13000",ges.dr5$TECH) & !grepl("13023",ges.dr5$TECH) & !grepl("10106",ges.dr5$TECH) & !grepl("10110",ges.dr5$TECH) & !grepl("10103",ges.dr5$TECH) & !grepl("10104",ges.dr5$TECH) & !grepl("10309",ges.dr5$TECH) & !grepl("12004",ges.dr5$TECH) & !grepl("13020",ges.dr5$TECH) & !grepl("13021",ges.dr5$TECH) & !grepl("13022",ges.dr5$TECH)
#ges.reserve <- ges.dr5

ges.dr5 <- ges.dr5[filter.param,]


bensby.param <- read.table('~/Survey/Resultados/DR6/Support/bensby_etal_2014_param_sorted.dat',header=T)      
roederer.param <- read.table('~/Survey/Resultados/DR6/Support/roederer_etal_2014_param_sorted.dat',header=T) 
#allende.param <- read.table('~/Survey/Resultados/DR6/Support/allendeprieto_etal_2004_param.dat',header=T)
#bruntt.param <- read.table('~/Survey/Resultados/DR6/Support/bruntt_etal_2012_param.dat',header=T)
#takeda.param <- read.table('~/Survey/Resultados/DR6/Support/takeda_etal_2005_param.dat',header=T)
#barklem.param <- read.table('~/Survey/Resultados/DR6/Support/barklem_etal_2005_param.dat',header=T)
#francois.param <- read.table('~/Survey/Resultados/DR6/Support/francois_eal_param.dat',header=T)      
#johnson.param <- read.table('~/Survey/Resultados/DR6/Support/johnson_2002_param.dat',header=T)
#yong.param <- read.table('~/Survey/Resultados/DR6/Support/yong_etal_2008_param.dat',header=T)
#bonifacio.param <- read.table('~/Survey/Resultados/DR6/Support/bonifacio_etal_2009_param.dat',header=T)      
#
#complete.sample <- as.data.frame(cbind(c(ges.dr5$TEFF,bensby.param$TEFF,roederer.param$TEFF,allende.param$TEFF,bruntt.param$TEFF,takeda.param$TEFF,barklem.param$TEFF,francois.param$TEFF,johnson.param$TEFF,yong.param$TEFF,bonifacio.param$TEFF),c(ges.dr5$LOGG,bensby.param$LOGG,roederer.param$LOGG,allende.param$LOGG,bruntt.param$LOGG,takeda.param$LOGG,barklem.param$LOGG,francois.param$LOGG,johnson.param$LOGG,yong.param$LOGG,bonifacio.param$LOGG),c(ges.dr5$FEH,bensby.param$FEH,roederer.param$FEH,(allende.param$FE-7.5),bruntt.param$FEH,takeda.param$FEH,barklem.param$FEH,francois.param$FEH,johnson.param$FEH,yong.param$FEH,bonifacio.param$FEH),c(ges.dr5$XI,bensby.param$XI,roederer.param$XI,allende.param$XI,bruntt.param$XI,takeda.param$XI,barklem.param$XI,francois.param$XI,johnson.param$XI,yong.param$XI,bonifacio.param$XI)))

complete.sample <- as.data.frame(cbind(c(ges.dr5$TEFF,bensby.param$TEFF,roederer.param$TEFF),c(ges.dr5$LOGG,bensby.param$LOGG,roederer.param$LOGG),c(ges.dr5$FEH,bensby.param$FEH,roederer.param$FEH),c(ges.dr5$XI,bensby.param$XI,roederer.param$XI)))

colnames(complete.sample) <- c('TEFF','LOGG','FEH','XI')


filter.main.sequence <- (complete.sample$LOGG >= 3.5) & (complete.sample$TEFF > 5250) 
filter.cool.stars <- (complete.sample$LOGG >= 3.5) & (complete.sample$TEFF <= 5250)
filter.giants <- (complete.sample$LOGG < 3.5) 
#
main.sequence <- complete.sample[filter.main.sequence,]
cool.stars <- complete.sample[filter.cool.stars,]
giant.stars <- complete.sample[filter.giants,]


N.stars.ms = nrow(main.sequence)
N.stars.cs = nrow(cool.stars)
N.stars.gs = nrow(giant.stars)

# renormalize
#
#
mean.teff <- mean(main.sequence$TEFF)
sd.teff <- sd(main.sequence$TEFF)
mean.logg <- mean(main.sequence$LOGG)
sd.logg <- sd(main.sequence$LOGG)
mean.feh <- mean(main.sequence$FEH)
sd.feh <- sd(main.sequence$FEH)
mean.xi <- mean(main.sequence$XI)
sd.xi <- sd(main.sequence$XI)
#
mean.teff.cs <- mean(cool.stars$TEFF)
sd.teff.cs <- sd(cool.stars$TEFF)
mean.logg.cs <- mean(cool.stars$LOGG)
sd.logg.cs <- sd(cool.stars$LOGG)
mean.feh.cs <- mean(cool.stars$FEH)
sd.feh.cs <- sd(cool.stars$FEH)
mean.xi.cs <- mean(cool.stars$XI)
sd.xi.cs <- sd(cool.stars$XI)
#
mean.teff.gs <- mean(giant.stars$TEFF)
sd.teff.gs <- sd(giant.stars$TEFF)
mean.logg.gs <- mean(giant.stars$LOGG)
sd.logg.gs <- sd(giant.stars$LOGG)
mean.feh.gs <- mean(giant.stars$FEH)
sd.feh.gs <- sd(giant.stars$FEH)
mean.xi.gs <- mean(giant.stars$XI)
sd.xi.gs <- sd(giant.stars$XI)
#
main.sequence$TEFF <- (main.sequence$TEFF-mean.teff)/(sd.teff)
main.sequence$LOGG <- (main.sequence$LOGG-mean.logg)/(sd.logg)
main.sequence$FEH  <- (main.sequence$FEH-(mean.feh))/(sd.feh)
main.sequence$XI  <- (main.sequence$XI-mean.xi)/(sd.xi)
#
cool.stars$TEFF <- (cool.stars$TEFF-mean.teff.cs)/(sd.teff.cs)
cool.stars$LOGG <- (cool.stars$LOGG-mean.logg.cs)/(sd.logg.cs)
cool.stars$FEH  <- (cool.stars$FEH-(mean.feh.cs))/(sd.feh.cs)
cool.stars$XI  <- (cool.stars$XI-mean.xi.cs)/(sd.xi.cs)
#
giant.stars$TEFF <- (giant.stars$TEFF-mean.teff.gs)/(sd.teff.gs)
giant.stars$LOGG <- (giant.stars$LOGG-mean.logg.gs)/(sd.logg.gs)
giant.stars$FEH  <- (giant.stars$FEH-(mean.feh.gs))/(sd.feh.gs)
giant.stars$XI  <- (giant.stars$XI-mean.xi.gs)/(sd.xi.gs)

#main.sequence$XIERR <- main.sequence$XIERR+0.05
#
main.sequence.data <- list(xi = main.sequence$XI,
#                           xi.sd = main.sequence$XIERR,
                           teff = main.sequence$TEFF,
#                           teff.sd = main.sequence$TEFFERR,
                           logg = main.sequence$LOGG,
#                           logg.sd = main.sequence$LOGGERR,
                           feh = main.sequence$FEH,
#                           feh.sd = main.sequence$FEHERR,
                           N = N.stars.ms
                           )
#
giant.stars.data <- list(xi = giant.stars$XI,
                           teff = giant.stars$TEFF,
                           logg = giant.stars$LOGG,
                           feh = giant.stars$FEH,
                           N = N.stars.gs
                           )
#
cool.stars.data <- list(xi = cool.stars$XI,
                           teff = cool.stars$TEFF,
                           logg = cool.stars$LOGG,
                           feh = cool.stars$FEH,
                           N = N.stars.cs
                           )

calib.model <- "model{
#Likelihood
for (i in 1:N) {
    xi[i]~dnorm(true.xi[i],typical.tau)
    true.xi[i] <- alpha[1] + alpha[2]*teff[i] + alpha[3]*pow(teff[i],2) + alpha[4]*logg[i] + alpha[5]*pow(logg[i],2) + alpha[6]*feh[i] + alpha[7]*pow(feh[i],2)
#+ alpha[8]*(teff[i]*logg[i]) + alpha[9]*(teff[i]*feh[i]) + alpha[10]*(logg[i]*feh[i])
}


#Priors for the coefficients
for (j in 1:7) {
alpha[j]~dnorm(0, 0.0001)
}
typical.tau <- pow(true.sigma,-2)
#true.sigma~dunif(0.01,1.0)
true.sigma~dgamma(1,1)
}"

my.inits <- function() {
    list(
        alpha = rnorm(7, 0, 0.0001)
    )
}
my.inits1 <- my.inits()
my.inits2 <- my.inits()
my.inits3 <- my.inits()
my.inits4 <- my.inits()
#my.inits5 <- my.inits()

my.params <- c("alpha","typical.tau")


simple_test.ms <- run.jags(method="rjags",
                        data=main.sequence.data,
                        inits=list(my.inits1,my.inits2,my.inits3,my.inits4),#,my.inits5),
                        mode=calib.model,
                        n.chains=4,
                        adapt=2000,
                        monitor=c(my.params),
                        burnin=2000,
                        sample=5000,
                        summarise=F,
                        thin=10,
                        plots=F
                        )

simple_test.cs <- run.jags(method="rjags",
                        data=cool.stars.data,
                        inits=list(my.inits1,my.inits2,my.inits3,my.inits4),#,my.inits5),
                        mode=calib.model,
                        n.chains=4,
                        adapt=2000,
                        monitor=c(my.params),
                        burnin=2000,
                        sample=5000,
                        summarise=F,
                        thin=10,
                        plots=F
                        )


simple_test.gs <- run.jags(method="rjags",
                        data=giant.stars.data,
                        inits=list(my.inits1,my.inits2,my.inits3,my.inits4),#,my.inits5),
                        mode=calib.model,
                        n.chains=4,
                        adapt=2000,
                        monitor=c(my.params),
                        burnin=2000,
                        sample=5000,
                        summarise=F,
                        thin=5,
                        plots=F
                        )

print(summary(simple_test.ms,vars=c('alpha',"typical.tau")))
fit.mcmc.ms <- as.mcmc(simple_test.ms)
print(xyplot(fit.mcmc.ms))
print(densityplot(fit.mcmc.ms))
print(heidel.diag(fit.mcmc.ms))
print(autocorr.diag(fit.mcmc.ms))
print(raftery.diag(fit.mcmc.ms))


print(summary(simple_test.cs,vars=c('alpha',"typical.tau")))
fit.mcmc.cs <- as.mcmc(simple_test.cs)
print(xyplot(fit.mcmc.cs))
print(densityplot(fit.mcmc.cs))
print(heidel.diag(fit.mcmc.cs))
print(autocorr.diag(fit.mcmc.cs))
print(raftery.diag(fit.mcmc.cs))


print(summary(simple_test.gs,vars=c('alpha',"typical.tau")))
fit.mcmc.gs <- as.mcmc(simple_test.gs)
print(xyplot(fit.mcmc.gs))
print(densityplot(fit.mcmc.gs))
print(heidel.diag(fit.mcmc.gs))
print(autocorr.diag(fit.mcmc.gs))
print(raftery.diag(fit.mcmc.gs))

coef.main.sequence.4 <- summary(simple_test.ms,vars=c('alpha'))[,4]
coef.cool.stars.4 <- summary(simple_test.cs,vars=c('alpha'))[,4]
coef.giant.stars.4 <- summary(simple_test.gs,vars=c('alpha'))[,4]


#coef.main.sequence <- summary(simple_test,vars=c('alpha'))[,4]

old.ges.main <- function(teff,logg,feh) {
    cor.teff <- (teff-5500)
    cor.logg <- (logg-4.0)
    #
    comp.xi <- 1.05 + (2.51E-4)*(cor.teff) + (1.5E-7)*(cor.teff^2) + (-0.14)*(cor.logg) + (-0.05E-1)*(cor.logg^2) + 0.05*feh + (0.01)*(feh^2)
    return(comp.xi)
}

old.ges.cool <- function(teff,logg,feh) {
    cor.teff <- (5000-5500)
    cor.logg <- (logg-4.0)
    #
    comp.xi <- 1.05 + (2.51E-4)*(cor.teff) + (1.5E-7)*(cor.teff^2) + (-0.14)*(cor.logg) + (-0.05E-1)*(cor.logg^2) + 0.05*feh + (0.01)*(feh^2)
    return(comp.xi)
}

old.ges.giants <- function(teff,logg,feh) {
    cor.teff <- (teff-5500)
    cor.logg <- (logg-4.0)
    #
    comp.xi <- 1.25 + (4.01E-4)*(cor.teff) + (3.1E-7)*(cor.teff^2) + (-0.14)*(cor.logg) + (-0.05)*(cor.logg^2) + 0.05*feh + (0.01)*(feh^2)
    return(comp.xi)
}



full.calib.main.sequence.1 <- function(teff,logg,feh) {
    #
    comp.xi <- 1.08+((6.30E-4)*(teff-5778))+((1.23E-7)*((teff-5778)^2))+((-3.28E-1)*(logg-4.18))+((1.47E-1)*((logg-4.18)^2))+((6.11E-2)*(feh+0.30))+((2.43E-2)*((feh+0.30)^2))
    #
    return(comp.xi)
}

full.calib.cool.stars.1 <- function(teff,logg,feh) {
    #
    comp.xi <- 0.48+((8.07E-4)*(teff-4981))+((1.02E-7)*((teff-4981)^2))+((-1.42E-1)*(logg-4.27))+((7.21E-1)*((logg-4.27)^2))+((2.38E-1)*(feh+0.28))+((1.89E-1)*((feh+0.28)^2))
    #
    return(comp.xi)
}

full.calib.giant.stars.1 <- function(teff,logg,feh) {
    #
    comp.xi <- 1.48+((5.00E-4)*(teff-4790))+((1.97E-7)*((teff-4790)^2))+((-5.16E-1)*(logg-2.35))+((-7.24E-2)*((logg-2.35)^2))+((2.26E-1)*(feh+0.76))+((5.52E-2)*((feh+0.76)^2))
    #
    return(comp.xi)
}

full.calib.main.sequence.4 <- function(teff,logg,feh) {
    cor.teff <- (teff-mean.teff)/sd.teff
    cor.logg <- (logg-mean.logg)/sd.logg
    cor.feh <- (feh-(mean.feh))/sd.feh
    #
    cor.xi <- coef.main.sequence.4[1]+(coef.main.sequence.4[2]*(cor.teff))+(coef.main.sequence.4[3]*(cor.teff^2))+(coef.main.sequence.4[4]*cor.logg)+(coef.main.sequence.4[5]*(cor.logg^2))+(coef.main.sequence.4[6]*cor.feh)+(coef.main.sequence.4[7]*(cor.feh^2))
    #
    comp.xi <- (cor.xi*sd.xi)+mean.xi
    return(comp.xi)
}

full.calib.cool.stars.4 <- function(teff,logg,feh) {
    cor.teff <- (teff-mean.teff.cs)/sd.teff.cs
    cor.logg <- (logg-mean.logg.cs)/sd.logg.cs
    cor.feh <- (feh-(mean.feh.cs))/sd.feh.cs
    #
    cor.xi <- coef.cool.stars.4[1]+(coef.cool.stars.4[2]*(cor.teff))+(coef.cool.stars.4[3]*(cor.teff^2))+(coef.cool.stars.4[4]*cor.logg)+(coef.cool.stars.4[5]*(cor.logg^2))+(coef.cool.stars.4[6]*cor.feh)+(coef.cool.stars.4[7]*(cor.feh^2))
    #
    comp.xi <- (cor.xi*sd.xi.cs)+mean.xi.cs
    return(comp.xi)
}

full.calib.giant.stars.4 <- function(teff,logg,feh) {
    cor.teff <- (teff-mean.teff.gs)/sd.teff.gs
    cor.logg <- (logg-mean.logg.gs)/sd.logg.gs
    cor.feh <- (feh-(mean.feh.gs))/sd.feh.gs
    #
    cor.xi <- coef.giant.stars.4[1]+(coef.giant.stars.4[2]*(cor.teff))+(coef.giant.stars.4[3]*(cor.teff^2))+(coef.giant.stars.4[4]*cor.logg)+(coef.giant.stars.4[5]*(cor.logg^2))+(coef.giant.stars.4[6]*cor.feh)+(coef.giant.stars.4[7]*(cor.feh^2))
    #
    comp.xi <- (cor.xi*sd.xi.gs)+mean.xi.gs
    return(comp.xi)
}
################
test.4.main <- full.calib.main.sequence.4(complete.sample$TEFF[filter.main.sequence],complete.sample$LOGG[filter.main.sequence],complete.sample$FEH[filter.main.sequence])
ges.4.main <- old.ges.main(complete.sample$TEFF[filter.main.sequence],complete.sample$LOGG[filter.main.sequence],complete.sample$FEH[filter.main.sequence])

test.4.cool <- full.calib.cool.stars.4(complete.sample$TEFF[filter.cool.stars],complete.sample$LOGG[filter.cool.stars],complete.sample$FEH[filter.cool.stars])
ges.4.cool <- old.ges.cool(complete.sample$TEFF[filter.cool.stars],complete.sample$LOGG[filter.cool.stars],complete.sample$FEH[filter.cool.stars])

test.4.giants <- full.calib.giant.stars.4(complete.sample$TEFF[filter.giants],complete.sample$LOGG[filter.giants],complete.sample$FEH[filter.giants])
ges.4.giants <- old.ges.giants(complete.sample$TEFF[filter.giants],complete.sample$LOGG[filter.giants],complete.sample$FEH[filter.giants])


postscript(file="xi_calib_teff_warm_main.eps",horizontal=FALSE, onefile=FALSE, height=6, width=6, pointsize=10)
par(oma=c(0.5,0.5,0.8,0.5))
par(mar=c(5,6,4,2))
par(mgp=c(4,1.2,0))
#
plotCI(complete.sample$TEFF[filter.main.sequence],complete.sample$XI[filter.main.sequence],xlab="Teff (K)",ylab="xi (km/s)",main='Warm main-sequence stars',xaxt='n',yaxt='n',cex.lab=1.4,cex.main=1.2)
plotCI(complete.sample$TEFF[filter.main.sequence],test.4.main,pch=19,col='red',add=T)
plotCI(complete.sample$TEFF[filter.main.sequence],ges.4.main,pch=19,col='blue',add=T)
leg.text <- c('GES iDR5','New calib.', 'Old. calib')
legend("topleft",legend=leg.text,pch=c(1,19,19),col=c('black','red','blue'),inset=0.02)

box(lwd=3)
axis(side=1, at=seq(1000,8000,500), labels=seq(1000,8000,500), cex.axis=1.5, lwd=2)
axis(side=1, at=seq(1000,8000,100), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
axis(side=2, at=seq(-5.0,5.0,0.5), labels=seq(-5.0,5.0,0.5), cex.axis=1.5, lwd=2)
axis(side=2, at=seq(-5.0,5.0,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
#
dev.off()


postscript(file="xi_calib_logg_warm_main.eps",horizontal=FALSE, onefile=FALSE, height=6, width=6, pointsize=10)
par(oma=c(0.5,0.5,0.8,0.5))
par(mar=c(5,6,4,2))
par(mgp=c(4,1.2,0))
#
plotCI(complete.sample$LOGG[filter.main.sequence],complete.sample$XI[filter.main.sequence],xlab="logg (dex)",ylab="xi (km/s)",main='Warm main-sequence stars',xaxt='n',yaxt='n',cex.lab=1.4,cex.main=1.2)
plotCI(complete.sample$LOGG[filter.main.sequence],test.4.main,pch=19,col='red',add=T)
plotCI(complete.sample$LOGG[filter.main.sequence],ges.4.main,pch=19,col='blue',add=T)
leg.text <- c('GES iDR5','New calib.', 'Old. calib')
legend("topleft",legend=leg.text,pch=c(1,19,19),col=c('black','red','blue'),inset=0.02)

box(lwd=3)
axis(side=1, at=seq(-1,5,0.5), labels=seq(-1,5,0.5), cex.axis=1.5, lwd=2)
axis(side=1, at=seq(-1,5,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
axis(side=2, at=seq(-5.0,5.0,0.5), labels=seq(-5.0,5.0,0.5), cex.axis=1.5, lwd=2)
axis(side=2, at=seq(-5.0,5.0,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
#
dev.off()


postscript(file="xi_calib_feh_warm_main.eps",horizontal=FALSE, onefile=FALSE, height=6, width=6, pointsize=10)
par(oma=c(0.5,0.5,0.8,0.5))
par(mar=c(5,6,4,2))
par(mgp=c(4,1.2,0))
#
plotCI(complete.sample$FEH[filter.main.sequence],complete.sample$XI[filter.main.sequence],xlab="[Fe/H] (dex)",ylab="xi (km/s)",main='Warm main-sequence stars',xaxt='n',yaxt='n',cex.lab=1.4,cex.main=1.2)
plotCI(complete.sample$FEH[filter.main.sequence],test.4.main,pch=19,col='red',add=T)
plotCI(complete.sample$FEH[filter.main.sequence],ges.4.main,pch=19,col='blue',add=T)
leg.text <- c('GES iDR5','New calib.', 'Old. calib')
legend("topleft",legend=leg.text,pch=c(1,19,19),col=c('black','red','blue'),inset=0.02)

box(lwd=3)
axis(side=1, at=seq(-5,5,0.5), labels=seq(-5,5,0.5), cex.axis=1.5, lwd=2)
axis(side=1, at=seq(-5,5,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
axis(side=2, at=seq(-5.0,5.0,0.5), labels=seq(-5.0,5.0,0.5), cex.axis=1.5, lwd=2)
axis(side=2, at=seq(-5.0,5.0,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
#
dev.off()





postscript(file="xi_calib_teff_cool_main.eps",horizontal=FALSE, onefile=FALSE, height=6, width=6, pointsize=10)
par(oma=c(0.5,0.5,0.8,0.5))
par(mar=c(5,6,4,2))
par(mgp=c(4,1.2,0))
#
plotCI(complete.sample$TEFF[filter.cool.stars],complete.sample$XI[filter.cool.stars],xlab="Teff (K)",ylab="xi (km/s)",main='Cool main-sequence stars',xaxt='n',yaxt='n',cex.lab=1.4,cex.main=1.2)
plotCI(complete.sample$TEFF[filter.cool.stars],test.4.cool,pch=19,col='red',add=T)
plotCI(complete.sample$TEFF[filter.cool.stars],ges.4.cool,pch=19,col='blue',add=T)
leg.text <- c('GES iDR5','New calib.', 'Old. calib')
legend("topleft",legend=leg.text,pch=c(1,19,19),col=c('black','red','blue'),inset=0.02)

box(lwd=3)
axis(side=1, at=seq(1000,8000,100), labels=seq(1000,8000,100), cex.axis=1.5, lwd=2)
axis(side=1, at=seq(1000,8000,10), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
axis(side=2, at=seq(-5.0,5.0,0.5), labels=seq(-5.0,5.0,0.5), cex.axis=1.5, lwd=2)
axis(side=2, at=seq(-5.0,5.0,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
#
dev.off()


postscript(file="xi_calib_logg_cool_main.eps",horizontal=FALSE, onefile=FALSE, height=6, width=6, pointsize=10)
par(oma=c(0.5,0.5,0.8,0.5))
par(mar=c(5,6,4,2))
par(mgp=c(4,1.2,0))
#
plotCI(complete.sample$LOGG[filter.cool.stars],complete.sample$XI[filter.cool.stars],xlab="logg (dex)",ylab="xi (km/s)",main='Cool main-sequence stars',xaxt='n',yaxt='n',cex.lab=1.4,cex.main=1.2)
plotCI(complete.sample$LOGG[filter.cool.stars],test.4.cool,pch=19,col='red',add=T)
plotCI(complete.sample$LOGG[filter.cool.stars],ges.4.cool,pch=19,col='blue',add=T)
leg.text <- c('GES iDR5','New calib.', 'Old. calib')
legend("topleft",legend=leg.text,pch=c(1,19,19),col=c('black','red','blue'),inset=0.02)

box(lwd=3)
axis(side=1, at=seq(-1,5,0.5), labels=seq(-1,5,0.5), cex.axis=1.5, lwd=2)
axis(side=1, at=seq(-1,5,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
axis(side=2, at=seq(-5.0,5.0,0.5), labels=seq(-5.0,5.0,0.5), cex.axis=1.5, lwd=2)
axis(side=2, at=seq(-5.0,5.0,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
#
dev.off()


postscript(file="xi_calib_feh_cool_main.eps",horizontal=FALSE, onefile=FALSE, height=6, width=6, pointsize=10)
par(oma=c(0.5,0.5,0.8,0.5))
par(mar=c(5,6,4,2))
par(mgp=c(4,1.2,0))
#
plotCI(complete.sample$FEH[filter.cool.stars],complete.sample$XI[filter.cool.stars],xlab="[Fe/H] (dex)",ylab="xi (km/s)",main='Cool main-sequence stars',xaxt='n',yaxt='n',cex.lab=1.4,cex.main=1.2)
plotCI(complete.sample$FEH[filter.cool.stars],test.4.cool,pch=19,col='red',add=T)
plotCI(complete.sample$FEH[filter.cool.stars],ges.4.cool,pch=19,col='blue',add=T)
leg.text <- c('GES iDR5','New calib.', 'Old. calib')
legend("topleft",legend=leg.text,pch=c(1,19,19),col=c('black','red','blue'),inset=0.02)

box(lwd=3)
axis(side=1, at=seq(-5,5,0.5), labels=seq(-5,5,0.5), cex.axis=1.5, lwd=2)
axis(side=1, at=seq(-5,5,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
axis(side=2, at=seq(-5.0,5.0,0.5), labels=seq(-5.0,5.0,0.5), cex.axis=1.5, lwd=2)
axis(side=2, at=seq(-5.0,5.0,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
#
dev.off()




postscript(file="xi_calib_teff_giants.eps",horizontal=FALSE, onefile=FALSE, height=6, width=6, pointsize=10)
par(oma=c(0.5,0.5,0.8,0.5))
par(mar=c(5,6,4,2))
par(mgp=c(4,1.2,0))
#
plotCI(complete.sample$TEFF[filter.giants],complete.sample$XI[filter.giants],xlab="Teff (K)",ylab="xi (km/s)",main='Giant stars',xaxt='n',yaxt='n',cex.lab=1.4,cex.main=1.2)
plotCI(complete.sample$TEFF[filter.giants],test.4.giants,pch=19,col='red',add=T)
plotCI(complete.sample$TEFF[filter.giants],ges.4.giants,pch=19,col='blue',add=T)
leg.text <- c('GES iDR5','New calib.', 'Old. calib')
legend("topleft",legend=leg.text,pch=c(1,19,19),col=c('black','red','blue'),inset=0.02)

box(lwd=3)
axis(side=1, at=seq(1000,8000,500), labels=seq(1000,8000,500), cex.axis=1.5, lwd=2)
axis(side=1, at=seq(1000,8000,100), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
axis(side=2, at=seq(-5.0,5.0,0.5), labels=seq(-5.0,5.0,0.5), cex.axis=1.5, lwd=2)
axis(side=2, at=seq(-5.0,5.0,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
#
dev.off()


postscript(file="xi_calib_logg_giants.eps",horizontal=FALSE, onefile=FALSE, height=6, width=6, pointsize=10)
par(oma=c(0.5,0.5,0.8,0.5))
par(mar=c(5,6,4,2))
par(mgp=c(4,1.2,0))
#
plotCI(complete.sample$LOGG[filter.giants],complete.sample$XI[filter.giants],xlab="logg (dex)",ylab="xi (km/s)",main='Giants stars',xaxt='n',yaxt='n',cex.lab=1.4,cex.main=1.2)
plotCI(complete.sample$LOGG[filter.giants],test.4.giants,pch=19,col='red',add=T)
plotCI(complete.sample$LOGG[filter.giants],ges.4.giants,pch=19,col='blue',add=T)
leg.text <- c('GES iDR5','New calib.', 'Old. calib')
legend("topleft",legend=leg.text,pch=c(1,19,19),col=c('black','red','blue'),inset=0.02)

box(lwd=3)
axis(side=1, at=seq(-1,5,0.5), labels=seq(-1,5,0.5), cex.axis=1.5, lwd=2)
axis(side=1, at=seq(-1,5,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
axis(side=2, at=seq(-5.0,5.0,0.5), labels=seq(-5.0,5.0,0.5), cex.axis=1.5, lwd=2)
axis(side=2, at=seq(-5.0,5.0,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
#
dev.off()


postscript(file="xi_calib_feh_giants.eps",horizontal=FALSE, onefile=FALSE, height=6, width=6, pointsize=10)
par(oma=c(0.5,0.5,0.8,0.5))
par(mar=c(5,6,4,2))
par(mgp=c(4,1.2,0))
#
plotCI(complete.sample$FEH[filter.giants],complete.sample$XI[filter.giants],xlab="[Fe/H] (dex)",ylab="xi (km/s)",main='Giant stars',xaxt='n',yaxt='n',cex.lab=1.4,cex.main=1.2)
plotCI(complete.sample$FEH[filter.giants],test.4.giants,pch=19,col='red',add=T)
plotCI(complete.sample$FEH[filter.giants],ges.4.giants,pch=19,col='blue',add=T)
leg.text <- c('GES iDR5','New calib.', 'Old. calib')
legend("topleft",legend=leg.text,pch=c(1,19,19),col=c('black','red','blue'),inset=0.02)

box(lwd=3)
axis(side=1, at=seq(-5,5,0.5), labels=seq(-5,5,0.5), cex.axis=1.5, lwd=2)
axis(side=1, at=seq(-5,5,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
axis(side=2, at=seq(-5.0,5.0,0.5), labels=seq(-5.0,5.0,0.5), cex.axis=1.5, lwd=2)
axis(side=2, at=seq(-5.0,5.0,0.1), labels=FALSE, cex.axis=1.5, lwd=1, tcl=-0.25)
#
dev.off()

################
par(def.par)#- reset to default
#

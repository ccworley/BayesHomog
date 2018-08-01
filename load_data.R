require(astro)

wg10file = '/Users/charlotteworley/Documents/GES/WG10/ParameterHomog/FITSFiles/20180220/GES_iDR5_WG10_Params_CombMat_20180220.fits'

# read FITS image
#x = read.fits(wg10file)
summary(x)

# show keywords in primary header
phdrkeywrd = x$hdr[[1]][,"key"]  #"value" , "comment"
ext1col = x$dat[[2]]$meta   #coluns names for binary table
ext1data = x$dat[[2]]$table  #binary table
cname = ext1data[,"CNAME"]
teff = as.numeric(ext1data[,"TEFF"])
logg = as.numeric(ext1data[,"LOGG"])

lmteff = as.numeric(ext1data[,"LM_TEFF"])
epteff = as.numeric(ext1data[,"EP_TEFF"])
icteff = as.numeric(ext1data[,"IC_TEFF"])
otteff = as.numeric(ext1data[,"OT_TEFF"])
lmlogg = as.numeric(ext1data[,"LM_LOGG"])
eplogg = as.numeric(ext1data[,"EP_LOGG"])
iclogg = as.numeric(ext1data[,"IC_LOGG"])
otlogg = as.numeric(ext1data[,"OT_LOGG"])


ireal <- which(is.finite(logg))

# create a plot
#layout(matrix(c(1,2,3,4), 2, 2), widths=c(1, 1), heights=c(1, 1))
par(mfrow = c(2, 4))
par(mar = c(4, 5, 1, 1))  #bottom, left, top, right
#par("mar"=c(3.1,3.1,1.1,1.1))

#Teff
plot(teff, lmteff-teff, pch = ".", xlab="Rec Teff", ylab="LM Teff")  #,ylim=rev(range(logg[ireal])),xlim=rev(range(teff[ireal])))
#label("topleft", txt="WG10", cex=1, lwd=0)
box()

plot(teff, epteff-teff, pch = ".", xlab="Rec Teff", ylab="EP Teff")
#label("topleft", txt="HDU 2", cex=2, lwd=0)
box()

plot(teff, icteff-teff, pch = ".", xlab="Rec Teff", ylab="IC Teff")  #,ylim=rev(range(logg[ireal])),xlim=rev(range(teff[ireal])))
#label("topleft", txt="WG10", cex=1, lwd=0)
box()

plot(teff, otteff-teff, pch = ".", xlab="Rec Teff", ylab="OT Teff")
#label("topleft", txt="HDU 2", cex=2, lwd=0)
box()

#Logg
plot(logg, lmlogg-logg, pch = ".", xlab="Rec logg", ylab="LM logg")  #,ylim=rev(range(logg[ireal])),xlim=rev(range(logg[ireal])))
#label("topleft", txt="WG10", cex=1, lwd=0)
box()

plot(logg, eplogg-logg, pch = ".", xlab="Rec logg", ylab="EP logg")
#label("topleft", txt="HDU 2", cex=2, lwd=0)
box()

plot(logg, iclogg-logg, pch = ".", xlab="Rec logg", ylab="IC logg")  #,ylim=rev(range(logg[ireal])),xlim=rev(range(logg[ireal])))
#label("topleft", txt="WG10", cex=1, lwd=0)
box()

plot(logg, otlogg-logg, pch = ".", xlab="Rec logg", ylab="OT logg")
#label("topleft", txt="HDU 2", cex=2, lwd=0)
box()

#par("mar"=c(3.1,0,1.1,1.1))
#aplot(sin, type="n", axes=FALSE, xlab="", ylab="")
#hdr = x$hdr[[2]]
#ktxt = paste("** astro.fits: HDU 2 Header **\n\nKey\n-----\n",
#             paste(hdr[,"key"],collapse="\n",sep=""),collapse="",sep="")
#vtxt = paste("\n\nValue\n-----\n",paste(hdr[,"value"],collapse="\n",sep=""),
#             collapse="",sep="")
#mtxt = paste("\n\nComment\n-----\n",paste(hdr[,"comment"],collapse="\n",sep=""),
#             collapse="",sep="")
#label("topleft", txt=ktxt, align="left", bty="n")
#label("topleft", txt=vtxt, align="left", bty="n", inset=c(1,0.08))
#label("topleft", txt=mtxt, align="left", bty="n", inset=c(2,0.08))
#label("bottom", txt="note: 'astro.fits' has been automatically deleted",
#      bty="n", col="blue", cex=1.5)
unlink(wg10file)
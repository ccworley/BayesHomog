
library("plot3D", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library");
# 3D Plot of Half of a Torus
par(mar = c(2, 2, 2, 2))
par(mfrow = c(1, 1))
R <- 3
r <- 2
x1 <- xtst
y1 <- xgsg
M <- mesh(x1, y1)
#zm <- (coef.alphas[1+(ind-1)*6]+coef.alphas[2+(ind-1)*6]*x1+coef.alphas[3+(ind-1)*6]*x1^2+coef.alphas[4+(ind-1)*6]*y1+coef.alphas[5+(ind-1)*6]*y1^2+coef.alphas[6+(ind-1)*6]*x1*y1)*sd.bench.teff
zm <- (coef.alphas[1+(ind-1)*6]+coef.alphas[2+(ind-1)*6]*M$x+coef.alphas[3+(ind-1)*6]*M$x^2+coef.alphas[4+(ind-1)*6]*M$y+coef.alphas[5+(ind-1)*6]*M$y^2+coef.alphas[6+(ind-1)*6]*M$x*M$y)*sd.bench.teff

x2 <- xt
y2 <- xg
M2 <- mesh(x2, y2)

surf3D(x = M2$x,
       y = M2$y,
       z = zm,
       xlab = "TEFF", ylab = "LOGG", zlab = "DeltaTEFF",
       ticktype = "detailed",
       colkey=FALSE,
       bty="b2",
       main="Bias Function ")


x <- given.teff.bench[star.code[vector.of.setups == 1]]
y <- given.logg.bench[star.code[vector.of.setups == 1]]
z <- (observed.node.teff.spectrum[(vector.of.setups == 1),ind]-given.teff.bench[star.code[vector.of.setups == 1]])


par(mfrow = c(2, 2))
scatter3D(x = x, y = y, z = z, ticktype = "detailed", pch = 16, bty = "f",theta=45, phi=45,
          colvar=z, xlab = "TEFF", ylab = "LOGG", zlab = "DeltaTEFF",
          surf=list(x = M2$x,
                    y = M2$y,
                    z = zm))
                    #colkey=FALSE,
                    #colvar=zm,
                    #bty="b2",
                    #main="Bias Function "))

scatter3D(x = x, y = y, z = z, ticktype = "detailed", pch = 16, bty = "f",theta=90, phi=45,
          colvar=z, xlab = "TEFF", ylab = "LOGG", zlab = "DeltaTEFF",
          surf=list(x = M2$x,
                    y = M2$y,
                    z = zm))

scatter3D(x = x, y = y, z = z, ticktype = "detailed", pch = 16, bty = "f",theta=0, phi=30,
          colvar=z, xlab = "TEFF", ylab = "LOGG", zlab = "DeltaTEFF",
          surf=list(x = M2$x,
                    y = M2$y,
                    z = zm))

scatter3D(x = x, y = y, z = z, ticktype = "detailed", pch = 16, bty = "f",theta=180, phi=25,
          colvar=z, xlab = "TEFF", ylab = "LOGG", zlab = "DeltaTEFF",
          surf=list(x = M2$x,
                    y = M2$y,
                    z = zm))

#colkey=FALSE,
#colvar=zm,
#bty="b2",
#main="Bias Function "))

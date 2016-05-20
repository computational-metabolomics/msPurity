iwNormGauss <- function(sdlim=3){

  # get a gaussian curve
  x <- seq(0-sdlim, 0+sdlim, 0.05)
  y <- dnorm(x, mean = 0)

  # linear scaling to range 0 to 1
  y <- (y-min(y))/(max(y)-min(y))
  x <- (x-min(x))/(max(x)-min(x))

  # Create function that outputs a 0 to 1 value depending on
  # the position of the x (i.e. mz)
  f <- approxfun(x, y)
}


iwNormRcosine <- function(){
  s <- sapa::taper(type="raised cosine")

  y <- as.vector(s)
  x <- seq(y)

  # linear scaling to range 0 to 1
  y <- (y-min(y))/(max(y)-min(y))
  x <- (x-min(x))/(max(x)-min(x))
  f <- approxfun(x, y)

}
#
#
#
# library(sapa)
# s <- taper(type="raised cosine", n.sample = 101)
#
# y <- as.vector(s)
# x <- seq(-0.5, +0.5,0.01)
#
# # linear scaling to range 0 to 1
# y <- (y-min(y))/(max(y)-min(y))
# f <- approxfun(x, y)
# par(mar = c(5,5,2,5))
# plot(f, xlim=c(-0.5, 0.5), xlab="Isolation window (Da)", ylab="Contribution")
# par(new = T)
#
#
# y2 <- rep(0, length(y))
# y2[7] <- 10000
# plot(x, y2, pch=16, axes=F, xlab=NA, ylab="", type="h", lwd=3)
# axis(side = 4)
# mtext(side = 4, line = 3, 'Intensity')
#
# lines(x=0, y=2300, type="h", lwd=3, col="red")





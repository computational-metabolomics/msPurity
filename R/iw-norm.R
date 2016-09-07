iwNormGauss <- function(sdlim=3, minOff=-0.5, maxOff=+0.5){

  # get a gaussian curve
  x <- seq(0-sdlim, 0+sdlim, 0.05)
  y <- dnorm(x, mean = 0)

  # linear scaling to range 0 to 1
  y <- (y-min(y))/(max(y)-min(y))

  # scaling to the min and max of the window
  x <- seq(minOff, maxOff, length.out = length(y))

  # Create function that outputs a 0 to 1 value depending on
  # the position of the x (i.e. mz)
  f <- approxfun(x, y)
}


iwNormQE.5 <- function(){
    y <- c(0.0000, 0.0000, 0.0000, 0.0550, 0.2336, 0.4437, 0.6509, 0.8210,
           0.9339, 0.9915, 0.9975, 0.9555, 0.8694, 0.7428, 0.5805, 0.3986,
           0.2208, 0.0710, 0.0000, 0.0000, 0.0000)
    x <- seq(-1, 1, 0.1)
    f <- approxfun(x, y)
}


iwNormRcosine <- function(minOff = -0.5, maxOff = +0.5){
   s <- sapa::taper(type="raised cosine")

   y <- as.vector(s)
   x <- seq(minOff, maxOff, length.out = length(y))

   # linear scaling to range 0 to 1
   y <- (y-min(y))/(max(y)-min(y))
   f <- approxfun(x, y)

}





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
   # s <- sapa::taper(type="raised cosine")
   # y <- as.vector(s)
   # Taken from above function from package sapa. Generates a rasied cosine curve
   #
   y <- c(3e-04, 0.001, 0.0023, 0.0041, 0.0064, 0.0091, 0.0123, 0.016, 0.02, 0.0244, 0.0291,
          0.0341, 0.0393, 0.0448, 0.0504, 0.0562, 0.062, 0.0679, 0.0738, 0.0796, 0.0853, 0.0909,
          0.0962, 0.1014, 0.1062, 0.1108, 0.115, 0.1188, 0.1222, 0.1252, 0.1277, 0.1297, 0.1313, 0.1323,
          0.1328, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329,
          0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329,
          0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1328, 0.1323, 0.1313, 0.1297, 0.1277,
          0.1252, 0.1222, 0.1188, 0.115, 0.1108, 0.1062, 0.1014, 0.0962, 0.0909, 0.0853, 0.0796, 0.0738,
          0.0679, 0.062, 0.0562, 0.0504, 0.0448, 0.0393, 0.0341, 0.0291, 0.0244, 0.02, 0.016, 0.0123, 0.0091,
          0.0064, 0.0041, 0.0023, 0.001, 3e-04)

   x <- seq(minOff, maxOff, length.out = length(y))

   # linear scaling to range 0 to 1
   y <- (y-min(y))/(max(y)-min(y))
   f <- approxfun(x, y)

}





context ("checking pcalc functions")

test_that("checking pcalc functions", {
  print ("\n")
  print("########################################################")
  print("## Checking pcalc functions                           ##")
  print("########################################################")

  # peak matrix with c13 single bond
  pm1 <- rbind(c(100, 1000),
              c(101.003, 10))

  # Check when the precursor is the M peak
  p1 <- pcalc(pm1, mzmin = 98, mzmax = 102, mztarget=100, ppm=5)
  p2 <- pcalc(pm1, mzmin = 98, mzmax = 102, mztarget=100, ppm=5, isotopes = TRUE)
  expect_equal(p1, c(0.990099, 2.000000))
  expect_equal(p2, c(1, 1))

  # Check when the precursor is the M+1 peak
  p3 <- pcalc(pm1, mzmin = 98, mzmax = 102, mztarget=101.003, ppm=5)
  p4 <- pcalc(pm1, mzmin = 98, mzmax = 102, mztarget=101.003, ppm=5, isotopes = TRUE)

  expect_equal(p3, c(0.00990099, 2.000000))
  expect_equal(p4, c(1, 1))

  # peak matrix with c13 double bond
  pm2 <- rbind(c(100, 1000),
              c(100.502, 10))

  # Check when the precursor is the M peak
  p5 <- pcalc(pm2, mzmin = 98, mzmax = 102, mztarget=100, ppm=5)
  p6 <- pcalc(pm2, mzmin = 98, mzmax = 102, mztarget=100, ppm=5, isotopes = TRUE)
  expect_equal(p5, c(0.990099, 2.000000))
  expect_equal(p6, c(1, 1))

})

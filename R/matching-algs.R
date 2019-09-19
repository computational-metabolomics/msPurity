cdpc <- function(L, Q){

  # Should only used the matching peaks for profile comparison
  remv_idx = c(which(L==0), which(Q==0))

  Lr <- L[-remv_idx]
  Qr <- Q[-remv_idx]

    # total length of matching peaks
  LnQ = length(Lr)

  if (LnQ<=1){
    # this can only be performed when there are adjacent values (i.e. there is hardcut off at 2)
    fr = 0
  }else{
    # loop through and get comparison of neighbour profiles
    profile_diff = rep(0, length(Qr)-1)
    c = 1
    for (i in 2:LnQ){
      abundance_ratio1 <- Lr[i]/Lr[i-1]
      abundance_ratio2 <- Qr[i]/Qr[i-1]
      if (abundance_ratio1<abundance_ratio2){
        n= 1
      }else{
        n = -1
      }

      prof <- abundance_ratio1^n * abundance_ratio2^-n

      profile_diff[c] <- prof

      c=c+1
    }

    # get an overall profile diff for the number of hits (notice that is one less than )
    fr = 1/(LnQ-1) * sum(profile_diff)

  }


  return((1/(length(Q)+LnQ)) * (length(Q)  * dpc(L, Q) + LnQ * fr))

}


dpc <- function(A,B) {
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}

overlap <- function(start1, end1, start2, end2){
  #Does the range (start1, end1) overlap with (start2, end2)? (returns boolean)
  #based on "De Morgan's laws", see http://nedbatchelder.com/blog/201310/range_overlap_in_two_compares.html
  return((end1 >= start2) & (end2 >= start1))

}


# recode.R
#
# Given a pair of vectors representing haplotypes coded with
# binary markers, recode as additive, recessive or dominant.
#
# Minor allele assumed to be the reference allele and
#
#       add  rec  dom  lof  2df
#  aa     1   1    1    1   0,0
#  aA     0   0    1   -2   0,1
#  Aa     0   0    1   -2   0,1
#  AA    -1   0    0    1   1,0
#
# David R. McWilliams <dmcwilli@wfubmc.edu>
#
# 14-Dec-2010 Start
# 05-Jan-2011 Generalize to any allele representation
# 22-Aug-2011 Add coding for lof test
# 19-Oct-2011 Handle the case where all are missing (return vector of NA)

recode <- function(h1, h2, model, missing=0)  {
  models <- c("additive", "recessive", "dominant", "lof", "2df")
  len <- length(h1)
  if(len != length(h2) ){
    stop("recode: Vectors h1 and h2 are different lengths.")
  } else if( !(grep(model, models, ignore.case=TRUE))) {
    stop("recode: Supplied model ", model, " not recognized.")
  }

  alleles <- (unique(c(h1,h2))[which(unique(c(h1,h2)) != missing)])
  
  if( length(alleles) ==2 ) {
    a1 <- alleles[1]
    a2 <- alleles[2]
  } else if(length(alleles) == 1) {
    a1 <- alleles[1]
    a2 <- alleles[1]
  } else if (length(alleles) == 0) {
    # Do nothing, return after allocating vector of NA below
  } else {
    stop("recode: Problem with allele number.")
  }

  # pre-allocate results container for efficiency
  if(length(grep(model, "2df", ignore.case=TRUE))) {
    r <- matrix(nrow=len, ncol=2)
  } else { 
    r <- as.vector(rep(NA, len))
  }

  if(length(alleles) == 0) {
    return(r)
  }

  num1 <- length(which(h1==a1)) + length(which(h2==a1))
  num2 <- length(which(h1==a2)) + length(which(h2==a2))
    
  if( length(grep(model, "additive", ignore.case=TRUE))) {
    for (i in 1:len) {
      if( h1[i] == a2 & h2[i] == a2 ) {
        r[i] <- ifelse( (num1 >= num2), 1, -1)
      }
      else if( (h1[i] == a1 & h2[i] == a2) | (h1[i] == a2 & h2[i] == a1) ) {
        r[i] <- 0
      }
      else if( h1[i] == a1 & h2[i] == a1 ) {
        r[i] <- ifelse( (num1 >= num2), -1, 1)
      }
    }
  } else if(length(grep(model, "recessive", ignore.case=TRUE))) {
    for (i in 1:len) {
      if( h1[i] == a2 & h2[i] == a2 ) {
        r[i] <- ifelse( (num1 >= num2), 1, 0)
      }
      else if( (h1[i] == a1 & h2[i] == a2) | (h1[i] == a2 & h2[i] == a1) ) {
        r[i] <- 0
      }
      else if( h1[i] == a1 & h2[i] == a1 ) {
        r[i] <- ifelse( (num1 >= num2), 0, 1)
      }
    }
  } else if(length(grep(model, "dominant", ignore.case=TRUE))) {
    for (i in 1:len) {
      if( h1[i] == a2 & h2[i] == a2 ) {
        r[i] <- ifelse( (num1 >= num2), 1, 0)
      }
      else if( (h1[i] == a1 & h2[i] == a2) | (h1[i] == a2 & h2[i] == a1) ) {
        r[i] = 1
      }
      else if( h1[i] == a1 & h2[i] == a1 ) {
        r[i] <- ifelse( (num1 >= num2), 0, 1)
      }
    }
  } else if(length(grep(model, "lof", ignore.case=TRUE))) {
    for (i in 1:len) {
      r[i] <- ifelse( (h1[i] == h2[i]), 1, -2) 
    }
  } else if(length(grep(model, "2df", ignore.case=TRUE))) {
    for(i in 1:len) {
      if( h1[i] == a2 & h2[i] == a2 ) {
        if(num1 >= num2) {
          r[i,] <- c(0,0)
        } else {
          r[i,] <- c(1,0)
        }
      } else if( (h1[i] == a1 & h2[i] == a2) | (h1[i] == a2 & h2[i] == a1) ) {
        r[i,] <- c(0,1)
      } else if(  h1[i] == a1 & h2[i] == a1 ) {
        if( num1 >= num2 ) {
          r[i,] <- c(1,0)
        } else {
          r[i,] <- c(0,0)
        }
      }
    }
  }
  r
} # end main

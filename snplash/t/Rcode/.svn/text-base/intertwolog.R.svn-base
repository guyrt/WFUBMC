# intertwolog.R
#
# Loop over columns of SNP data, recoding the genetic model and
# predicting the response using the resulting vector.  This version is
# meant to reproduce intertwolog and hence expects 2 covariates.
#
# Depends on the function 'recode' defined in recode.R
#
# David R. McWilliams <dmcwilli@wfubmc.edu>
#
# 20-May-2011 Adapt multi.lm.R for glm and to test intertwolog

# intertwolog.glm
#
# Accepts: A dataframe containing a response and snp data in biallelic
#          format, the response column index, and the index of the
#          start column for the snp data, which must be in one
#          contiguous block.
#
# Returns: A dataframe with the coefficient information provided by glm()
#

intertwolog <- function (data.df, qvar, snp.start, gen.model, cov1, cov2) {

  # call:

  # subset.i2l.centered.df <- intertwolog(sim2000.geno.df,
  #            sim2000.response-1,
  #            snp.start=2,
  #            gen.model="additive",
  #            sim2000.phen.df$cov1,
  #            sim2000.phen.df$cov2))

  debug <- 0
  
  orig.col <- ncol(data.df)
  d <- data.df[, snp.start:orig.col]
  n.col <- ncol(d)
  if(n.col %% 2) {
    print("Odd number of columns: aborting.")
    return(0)
  }
  n.snp <- n.col/2

  # Pre-allocate a matrix to hold the results.
  #
  # For 2 variables and two covariates there are 6 coefficients and 4
  # statistics.  Since this is logistic regression, might as well save
  # the AIC and deviance as well. i.e. 6*4+2=26 columns per
  # regression.  Add 2 columns to hold the snp ids, 26+2 = 28.
  #
  # The number of rows is n!/(n-r)!r!, which simplifies to
  # n*(n-1)/2 for r=2 as here.
  results.m <- matrix(nrow=n.snp*(n.snp-1)/2, ncol=28)
  row.num <- 1

  if(debug) {
    cat(c("dim results.m: ", dim(results.m), fill=T))
    print("\n")
  }
  
  for (i in seq(1, n.col-3, 2)) {
    snp1.num <- (i+1)/2

    for (j in seq(i+2, n.col-1, 2)) {
      snp2.num <- (j+1)/2

      if(debug) {
        cat(c(i, snp1.num, j, snp2.num), fill=T)
      }

      # Subset the dataframe, removing lines with zeros.  Note that
      # this is by snps.  We assume proper QC on the genotype and snp
      # data has been performed.  Must 'and' all 4 columns to prevent
      # unequal length vectors.

      not.zero <- (d[,i] !=0 & d[,i+1] !=0 & d[,j] !=0 & d[,j+1] !=0)
      sub.df   <- d[not.zero, c(i,i+1,j,j+1)]

      snp1     <- recode(sub.df[,1], sub.df[,2], gen.model)
      cross1   <- snp1 - mean(snp1)

      snp2     <- recode(sub.df[,3], sub.df[,4], gen.model)
      cross2   <- snp2 - mean(snp2)
      
      c1       <- cov1[not.zero]
      c2       <- cov2[not.zero]

      y        <- qvar[not.zero]

      this.glm <- glm(y ~ snp1 + snp2 + cross1:cross2 + c1 + c2,
                      family="binomial")

      if(debug) {
        print(this.glm)
        print("")
      }
      
      glm.sum <- summary(this.glm)$coefficients
      if( (dim(glm.sum))[1] > 1) {
        new.row <- c(
                     snp1.num,
                     snp2.num,
                     this.glm$deviance,
                     this.glm$aic,
                     t(glm.sum[,1]),
                     t(glm.sum[,2]),
                     t(glm.sum[,3]),
                     t(glm.sum[,4]))
        if( length(new.row) == 28) {
          results.m[row.num, ] <- new.row
        } else {
          results.m[row.num, ] <- c(snp1.num, snp2.num, rep(NA,26))
        }
        
      } else {
        results.m[row.num, ] <- c(snp1.num, snp2.num, rep(NA,26))
      }
      row.num <- row.num + 1

      not.zero <- NULL
      sub.df   <- NULL
      snp1     <- NULL
      snp2     <- NULL
      c1       <- NULL
      c2       <- NULL
      y        <- NULL  
      this.glm <- NULL  
    } # end for(j ...)
  } # end for(i ...)

  results.df <- data.frame(results.m)

  names(results.df) <- c("snp1", "snp2", "dev", "aic",
                         "int",    "b1",    "b2",    "cov1",    "cov2",    "b1.b2",
                         "int.se", "b1.se", "b2.se", "cov1.se", "cov2.se", "b1.b2.se",
                         "int.z",  "b1.z",  "b2.z",  "cov1.z",  "cov2.z",  "b1.b2.z",
                         "int.p",  "b1.p",  "b2.p",  "cov1.p",  "cov2.p",  "b1.b2.p"
                         )
  results.df
} # end main

# par.glm.cov.R
#
# Loop over columns of SNP data, recoding the genetic model
# and predicting the response using the resulting vector.
#
# Depends on the function 'recode' defined in recode.R
#
# David R. McWilliams <dmcwilli@wfubmc.edu>
#
# 15-Dec-2010 Initiate project
# 20-Dec-2010 Pass in model as a parameter
# 09-Jun-2011 Make parallel; adapted from multi.lm

# par.glm.cov
#
# Accepts: A dataframe containing a response and snp data in biallelic
#          format, the response column index, and the index of the
#          start column for the snp data, which must be in one
#          contiguous block.
#
# Returns: A dataframe with the coefficient information provided by (g)lm
#

# Need test to see if package doMC is already loaded

par.glm.cov <- function (data.df, qvar, snp.start, gen.model, cov1, cov2) {
  d <- data.df
  n.col <- ncol(d)
  n.snp <- (n.col-snp.start+1)/2

  snp.num <- 1

  mcols <- 4*4+2
  results.m <- matrix(nrow=n.snp, ncol=mcols)
  
  results.m <-foreach(i=seq(snp.start, n.col-1, 2), .combine=rbind) %dopar% {
    nz <- d[,i]!=0 & d[,i+1]!=0
    snp <- recode(d[nz,i], d[nz,i+1], gen.model)

    if (all(is.na(snp))) {
      rep(NA, mcols)
    } else {
      int.res <- rep(NA,4)
      snp.res <- rep(NA,4)
      cv1.res <- rep(NA,4)
      cv2.res <- rep(NA,4)
      aic.res <- NA
      dev.res <- NA

      this.glm <- glm(qvar[nz] ~ snp + cov1[nz] + cov2[nz], family="binomial")
      this.coef <- summary(this.glm)$coefficients
      snp.num <- snp.num + 1

      nms <- row.names(this.coef)
      for (i in 1:length(nms)) {
        if(length(grep("Int", nms[i], value=T)) > 0) {
          int.res <- this.coef[i,]
        }
        if(length(grep("snp", nms[i], value=T)) > 0) {
          snp.res <- this.coef[i,]
        }
        if(length(grep("cov1", nms[i], value=T)) > 0) {
          cv1.res <- this.coef[i,]
        }
        if(length(grep("cov2", nms[i], value=T)) > 0) {
          cv2.res <- this.coef[i,]
        }
      }
      aic.res <- this.glm$aic
      dev.res <- this.glm$deviance

      c(int.res, snp.res, cv1.res, cv2.res, aic.res, dev.res)
    }
  } # end foreach(i in snp)
#  snp.num <- snp.num -1

  results.df <- data.frame(results.m)
  row.names(results.df) <- NULL

  names(results.df) <- c("b0", "b0.se", "b0.z", "b0.p",
                         "b1", "b1.se", "b1.z", "b1.p",
                         "c1", "c1.se", "c1.z", "c1.p",
                         "c2", "c2.se", "c2.z", "c2.p",
                         "aic", "null.dev")
  results.df
} # end main

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
# 22-Jan-2013 Include the additive model for a snp as a covariate
#             in the lack-of-fit calculation.

# Accepts: A dataframe containing a response and snp data in biallelic
#          format, the response column index, and the index of the
#          start column for the snp data, which must be in one
#          contiguous block.
#
# Returns: A dataframe with the coefficient information provided by glm
#
# Requires: The doMC package must be already loaded.

# ToDo: Test for package doMC and fail gracefully if absent

par.glm.cov <- function (data.df, qvar, snp.start, gen.model, cov1, cov2) {
  d <- data.df
  n.col <- ncol(d)
  n.snp <- (n.col-snp.start+1)/2

  is.lof <- (gen.model=="lof")
  
  snp.num <- 1

  mcols <- 4*4+2
  if(is.lof) {
    mcols <- mcols+4
  }
  results.m <- matrix(nrow=n.snp, ncol=mcols)
  
  results.m <-foreach(i=seq(snp.start, n.col-1, 2), .combine=rbind) %dopar% {
    nz <- d[,i]!=0 & d[,i+1]!=0
    snp <- recode(d[nz,i], d[nz,i+1], gen.model)

    if(is.lof) {
      snp.add <- recode(d[nz,i], d[nz,i+1], "add")
    }

    if (all(is.na(snp))) {
      rep(NA, mcols)
    } else {
      int.res     <- rep(NA,4)
      snp.res     <- rep(NA,4)
      cv1.res     <- rep(NA,4)
      cv2.res     <- rep(NA,4)
      snp.add.res <- rep(NA,4)
      aic.res     <- NA
      dev.res     <- NA

      if(is.lof) {
        this.glm <- glm(qvar[nz] ~ snp + cov1[nz] + cov2[nz] + snp.add, family="binomial")
      } else {
        this.glm <- glm(qvar[nz] ~ snp + cov1[nz] + cov2[nz], family="binomial")
      }
      this.coef <- summary(this.glm)$coefficients
      snp.num <- snp.num + 1

      nms <- row.names(this.coef)
      for (i in 1:length(nms)) {
        if(length(grep("Int", nms[i], value=T)) > 0) {
          int.res <- this.coef[i,]
        }
        if(length(grep("^snp$", nms[i], value=T)) > 0) {
          snp.res <- this.coef[i,]
        }
        if(length(grep("cov1", nms[i], value=T)) > 0) {
          cv1.res <- this.coef[i,]
        }
        if(length(grep("cov2", nms[i], value=T)) > 0) {
          cv2.res <- this.coef[i,]
        }
        if(length(grep("snp.add", nms[i], value=T)) > 0) {
          snp.add.res <- this.coef[i,]
        }
      }

      out.row <- c(int.res, snp.res, cv1.res, cv2.res)
      if(is.lof) {
        out.row <- append(out.row, snp.add.res)
      }
      aic.res <- this.glm$aic
      dev.res <- this.glm$deviance

      out.row <- append(out.row, aic.res)
      out.row <- append(out.row, dev.res)
    }
  } # end foreach(i in snp)

  results.df <- data.frame(results.m)
  row.names(results.df) <- NULL

  var.names <- c("b0", "b0.se", "b0.z", "b0.p",
                 "b1", "b1.se", "b1.z", "b1.p",
                 "c1", "c1.se", "c1.z", "c1.p",
                 "c2", "c2.se", "c2.z", "c2.p")

  if(is.lof) {
    var.names <- append(var.names,
                        c("snp.add", "snp.add.se", "snp.add.z", "snp.add.p"))
  }
  var.names <- append(var.names, c("aic", "null.dev"))
  names(results.df) <- var.names
  
  results.df
} # end main

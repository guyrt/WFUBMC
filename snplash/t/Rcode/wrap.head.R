wrap.head <- function(x) {
  for(i in 1:length(x)) {
    if(nchar(x[i]) == 0) {
      cat("\n")
      next
    }
    if(length(grep("\\*+", x[i]))>0) {
      cat("\n")
      next
    }
    writeLines(strwrap(x[i]))
  }
} # end wrap.head()
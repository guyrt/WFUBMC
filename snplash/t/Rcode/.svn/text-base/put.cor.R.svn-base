# put.cor.R
#
# Function to place the correlation for the test_<engine>.Snw
# scripts.
#
# David R. McWilliams <dmcwilli@wfubmc.edu>
# 01-Sep-2011

put.cor <- function(x,y) {
  # put the correlation value in the lower right corner
  plot.lims <- par("usr")
  t.x  <- plot.lims[1] + 0.75*(plot.lims[2] - plot.lims[1])
  t.y  <- plot.lims[3] + 0.20*(plot.lims[4] - plot.lims[3])

  text(t.x, t.y,
     paste("r = ", round(cor(x,y, use="complete"), 6)))

} # end put.cor()

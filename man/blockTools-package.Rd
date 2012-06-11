\name{blockTools-package}
\alias{blockTools-package}
\alias{blockTools}
\docType{package}
\title{Block, randomly assign, and diagnose potential interference in randomized experiments}
\description{Block units into experimental blocks, with one unit per
treatment condition, by creating a measure of multivariate distance
between all possible pairs of units.  Maximum, minimum, or an
allowable range of differences between units on one variable can be
set.  Randomly assign units to treatment conditions.  Diagnose
potential interference problems between units assigned to different
treatment conditions. Write outputs to .tex and .csv files.
}
\details{
\tabular{ll}{
Package: \tab blockTools\cr
Type: \tab Package\cr
Version: \tab 0.5-5\cr
Date: \tab 2012-06-11\cr
License: \tab  GPL (>=2)\cr
}
Given raw data, \code{block} creates experimental blocks,
\code{assignment} assigns units to treatment conditions, \code{diagnose}
detects possible interference problems, and \code{outTeX} and
\code{outCSV} write block or assignment output objects to a set of .tex
and .csv files, respectively. 
}
\author{
Ryan T. Moore

Maintainer: Ryan T. Moore \email{rtm@wustl.edu}
}
\references{
\url{http://rtm.wustl.edu/software.blockTools.htm}
}

\keyword{ package }

\examples{
data(x100)

## block
out <- block(x100, groups = "g", n.tr = 2, id.vars = c("id"), block.vars
             = c("b1", "b2"), algorithm="optGreedy", distance =
             "mahalanobis", level.two = FALSE, valid.var = "b1",
             valid.range = c(0,500), verbose = TRUE)
## assign
assg <- assignment(out, seed = 123)
## diagnose
diag <- diagnose(object = assg, data = x100, id.vars = "id",
                suspect.var = "b2", suspect.range = c(0,50))
## create .tex files of assigned blocks
outTeX(assg)
## create .csv files of unassigned blocks
outCSV(out)
}

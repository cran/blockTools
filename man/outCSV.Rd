\name{outCSV}
\alias{outCSV}
\title{Export blocked or assigned data to .csv format files}
\description{Exports output from \code{block} or \code{assignment} to a
  set of .csv files using \code{write.csv}.
}
\usage{
outCSV(block.obj, namesCol = NULL, digits = 2, ...)
}
\arguments{
  \item{block.obj}{a list of dataframes, such as output from
    \code{block} or \code{assignment}.}
  \item{namesCol}{an optional vector of column names to be used in
    output files.}
  \item{digits}{number of digits to which to round multivariate
    distances in output files.}
  \item{\dots}{additional arguments passed to \code{write.csv}.}
}

\value{
A set of .csv files, one for each element of the input list of blocked
or assigned units.  Each file is named ``GroupXXX.csv", where ``XXX" is
the group name taken from the input object.
}

\author{Ryan T. Moore}

\seealso{\code{\link{outTeX}}, \code{\link{write.csv}},
  \code{\link{block}}, \code{\link{assignment}}
}
\examples{
data(x100)

## First, block
out <- block(x100, groups = "g", n.tr = 2, id.vars = c("id"), block.vars
             = c("b1", "b2"), algorithm="optGreedy", distance =
             "mahalanobis", level.two = FALSE, valid.var = "b1",
             valid.range = c(0,500), verbose = TRUE)
## Second, assign
assg <- assignment(out, seed = 123)

## create three .csv files of blocks
outCSV(out)
## create three .csv files of assigned blocks
##   (note: overwrites blocked .csv files)
outCSV(assg)
}

\keyword{design}
\keyword{IO}
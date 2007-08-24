\name{outTeX}
\alias{outTeX}
\title{Export blocked or assigned data to .tex format files} 
\description{Exports output from \code{block} or \code{assignment} to a
  set of .tex files using \code{xtable}.
}
\usage{
outTeX(block.obj, namesCol = NULL, digits = 2, ...)
}
\arguments{
  \item{block.obj}{a list of dataframes, such as output from
    \code{block} or \code{assignment}.}
  \item{namesCol}{an optional vector of column names to be used in
    output files.}
  \item{digits}{the number of digits to which to round multivariate
    distances in output files.}
  \item{\dots}{additional arguments passed to \code{write.tex}.}
}

\value{
A set of .tex files, one for each element of the input list of blocked
or assigned units.  Each file is named ``GroupXXX.tex", where ``XXX" is
the group name taken from the input object.  The tables in these .tex
files can easily be integrated into an existing .tex document using
\LaTeX code \code{\include{GroupXXX}}.
}

\author{Ryan T. Moore}

\seealso{\code{\link{outCSV}}, \code{\link{block}}, \code{\link{assignment}}
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

## create three .tex files of blocks
outTeX(out)
## create three .tex files of assigned blocks
##   (note: overwrites blocked .tex files)
outTeX(assg)
}
\keyword{design}
\keyword{IO}
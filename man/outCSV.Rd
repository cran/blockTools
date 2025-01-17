% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outCSV.R
\name{outCSV}
\alias{outCSV}
\title{Export blocked or assigned data to .csv format files}
\usage{
outCSV(block.obj, namesCol = NULL, file.names = NULL, digits = 2, ...)
}
\arguments{
\item{block.obj}{A list of dataframes, such as output from \code{block} or \code{assignment}.}

\item{namesCol}{An optional character vector of column names to be used in output files.}

\item{file.names}{An optional character vector of file names specifying the output file names.}

\item{digits}{An integer representing the number of decimal places to which to round multivariate distances in output files, passed to \code{round()}.}

\item{...}{Additional arguments passed to \code{write.csv}.}
}
\value{
A set of .csv files, one for each element of the input list of blocked or assigned units, written by \code{write.csv()}.
}
\description{
Exports output from \code{block} or \code{assignment} to a set of .csv files using \code{write.csv}.
}
\details{
Under the default (\code{file.names = NULL}), each file is named \dQuote{GroupXXX.csv}, where \dQuote{XXX} is the group name taken from the input object.
}
\examples{
data(x100)

# Block and assign:
out <- block(x100, groups = "g", n.tr = 2, id.vars = "id", block.vars = c("b1", "b2"))
assg <- assignment(out, seed = 123)

# create three .csv files of blocks
\dontrun{outCSV(out)}
# create three .csv files of assigned blocks
#   (note: overwrites blocked .csv files)
\dontrun{outCSV(assg)}
# create three .csv files with custom file names
\dontrun{outCSV(assg, file.names = c("file1", "file2", "file3"))}

}
\seealso{
\code{\link{outTeX}}, \code{\link{write.csv}}, \code{\link{block}}, \code{\link{assignment}}
}
\author{
Ryan T. Moore
}
\keyword{IO}
\keyword{design}

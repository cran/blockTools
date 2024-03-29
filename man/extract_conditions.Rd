\name{extract_conditions}
\alias{extract_conditions}
\title{Create vector of integers containing treatment condition identifiers}
\description{
  Creates a vector of integers which represent unique treatment conditions 
  in an object output from \code{assignment}.
}
\usage{
extract_conditions(assg.obj, data, id.var)
}

\arguments{
  \item{assg.obj}{an output object from \code{assignment}.}
  \item{data}{the data frame that was input into \code{block} for blocking.}
  \item{id.var}{a string specifying which column of \code{data} contains identifying information.}
}

\details{
Under the current implementation, \code{level.two} in \code{block}
should be set to \code{FALSE}.
}

\value{ 
A numeric vector of integers with \code{nrow(data)} elements with
lowest value equal to 1, corresponding to the treatment condition column 
from the assignment object each unit is in. For example, if the columns of the assignment object are \code{Treatment} and \code{Control} (in that order), then \code{Treatment} will be represented by a \code{1} and \code{Control} will be represented by a \code{2}. 

For units in \code{data} that are not in \code{assg.obj}, the value of NA
is assigned.
}

\author{Ryan T. Moore}

\seealso{\code{\link{assignment}}}

\examples{
data(x100)
out <- block(x100, groups = "g", n.tr = 2, id.vars = c("id"), 
             block.vars = c("b1", "b2"))

assg <- assignment(out)

extract_conditions(assg, x100, id.var = "id")

# Treatment conditions represented by integers.
}

\keyword{design}


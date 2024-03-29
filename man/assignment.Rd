\name{assignment}
\alias{assignment}

\title{Randomly assign blocked units to treatment conditions}

\description{Using an output object from \code{block}, assign elements
  of each row to treatment condition columns.  Each element is equally 
  likely to be assigned to each column.}

\usage{ assignment(block.obj, seed = NULL, namesCol = NULL) }

\arguments{
  \item{block.obj}{an output object from \code{block}, or a user-specified
    block object.}
  \item{seed}{a user-specified random seed.}
  \item{namesCol}{an optional vector of column names for the output
    table.}
}
\details{

\code{block.obj} can be specified directly by the user.  It can be a
single dataframe or matrix with blocks as rows and treatment
conditions as columns.  \code{assignment} is designed to take a list
with two elements.  The first element should be named \code{$blocks},
and should be a list of dataframes.  Each dataframe should have blocks
as rows and treatment conditions as columns.  The second element
should be a logical named \code{$level.two}.  A third element, such as
\code{$call} in a \code{block} output object, is currently ignored.
  
Specifying the random seed yields constant assignment, and thus allows
for easy replication of experimental protocols.

If \code{namesCol = NULL}, then \dQuote{Treatment 1},
\dQuote{Treatment 2}, \ldots are used.  If \code{namesCol} is supplied by the user 
and is of length \code{n.tr} (or 2*\code{n.tr}, where \code{level.two = TRUE}), 
then either \code{"Distance"} or \code{"Max Distance"} is appended to it as 
appropriate (consistent with \code{namesCol} usage in \code{block}).  If 
\code{namesCol} is supplied and is of length \code{n.tr}+1 (or 2*\code{n.tr}+1, 
where \code{level.two = TRUE}), then the last user-supplied name is used for the 
last column of each dataframe.
}

\value{A list with elements
  \item{assg}{a list of dataframes, each containing a group's blocked
    units assigned to treatment conditions.  If there are two treatment
    conditions, then the last column of each dataframe displays the
    multivariate distance between the two units.  If there are more than
    two treatment conditions, then the last column of each dataframe
    displays the largest of the multivariate distances between all
    possible pairs in the block.}
  \item{call}{the orginal call to \code{assignment}.}
}

\author{Ryan T. Moore}

\seealso{\code{\link{block}}, \code{\link{diagnose}}}
\examples{
data(x100)

# First, block
out <- block(x100, groups = "g", n.tr = 2, id.vars = c("id"), block.vars
             = c("b1", "b2"), algorithm="optGreedy", distance =
             "mahalanobis", level.two = FALSE, valid.var = "b1",
             valid.range = c(0,500), verbose = TRUE)

# Second, assign
assigned <- assignment(out, seed = 123)

# assigned$assg contains 3 data frames
}

\keyword{design}


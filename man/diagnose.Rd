\name{diagnose}
\alias{diagnose}

\title{Diagnose whether units assigned to different treatment
  conditions may be subject to interference or pairwise imbalance}
\description{List all pairs of units assigned to different treatment
  conditions whose difference on a specified variable falls within a
  specified range.
}
\usage{
diagnose(object, data, id.vars, suspect.var, suspect.range = NULL)
}

\arguments{
  \item{object}{a dataframe or list of dataframes of assigned units,
    such as output from \code{assignment}.}
  \item{data}{a dataframe with auxiliary information on assigned units,
    including the specified variable \code{suspect.var}.}
  \item{id.vars}{a required string or vector of two strings specifying which
    column(s) of \code{data} contain identifying information.}
  \item{suspect.var}{a string specifying which column of \code{data}
    contains the variable suspected of interference or imbalance.}
  \item{suspect.range}{a vector defining the range of
    \code{suspect.var} within which units in different treatment
    conditions must fall to be considered suspect.}
}
\details{
\code{object} requires rows to correspond to blocks and columns to
correspond to treatment conditions, such as output from \code{assignment}.

\code{data} should include identifying variables and variable suspected
of interference or imbalance.  Typically, \code{data} may be the same
dataframe input into \code{block}.

An example of specified identifying variables is \code{id.vars = c("id",
  "id2")}.  Unlike \code{block}, \code{diagnose} requires that the length of
\code{id.vars} correspond to the level of the original blocking.  See
\code{block} documentation for details.

An example of specified suspect range is \code{suspect.var = "b2"},
\code{suspect.range = c(0,50)} identifies all units assigned to
different treatment conditions no more than 50 units apart on variable
\code{"b2"}. 
}

\value{
A list of dataframes, each containing a group's pairs of units assigned
to different treatments falling within \code{suspect.range} on the
variable \code{suspect.var}.  The last column of each dataframe displays
the observed difference between the two units. 
}

\author{Ryan T. Moore}

\seealso{\code{\link{assignment}}, \code{\link{block}}}
\examples{
data(x100)

## First, block
out <- block(x100, groups = "g", n.tr = 2, id.vars = c("id"), block.vars
             = c("b1", "b2"), algorithm="optGreedy", distance =
             "mahalanobis", level.two = FALSE, valid.var = "b1",
             valid.range = c(0,500), verbose = TRUE)
## Second, assign
assg <- assignment(out, seed = 123)
## Third, diagnose
diag <- diagnose(object = assg, data = x100, id.vars = "id",
                suspect.var = "b2", suspect.range = c(0,50))
}

\keyword{design}


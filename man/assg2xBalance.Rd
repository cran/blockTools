\name{assg2xBalance}
\alias{assg2xBalance}
\title{Calculate balance statistics from an assignment object}
\description{Calculate several balance statistics for experimental units assigned to treatment conditions.  Naturally accepts output from the \code{assignment} function, and passes it to \code{xBalance} from \code{library(RItools)}.  Provides balance summaries for the entire experiment and by group.    
}
\usage{
assg2xBalance(assg.obj, data, id.var, bal.vars, to.report = "all")
}


\arguments{
  \item{assg.obj}{an output object from \code{assignment}.}
  \item{data}{the data frame that was input into \code{block} for blocking.}
  \item{id.var}{a string specifying the column of \code{data} containing identifying information.}
  \item{bal.vars}{a string or vector of strings specifying which column(s) of \code{data} contain the variables on which balance is to be checked.}
  \item{to.report}{a string or vector of strings passed to \code{xBalance} listing the measures to report for each group.  See Details for more information.}
}

\details{
As of \code{RItools} version 0.1-11, \code{to.report} must be a subset of \code{c("std.diffs","z.scores","adj.means", "adj.mean.diffs","adj.mean.diffs.null.sd", "chisquare.test","p.values", "all")}.  The default, \code{all}, returns all measures.
}

\value{
A list of output objects from \code{xBalance}.  For each group defined in the \code{assignment} object, one list element is assigned the name of that group and summarizes the balance in that group according to \code{to.report}.  The last element of the list is named \code{"Overall"} and summarizes balance across all groups.  The elements of this list are themselves objects of class \code{c("xbal", "list")}.  

If \code{assg.obj} has only one group, the first element of the output list is named \code{"Group1"}, and the second is named \code{"Output"}.  In this case, these two elements will be identical.
}

\references{
Hansen, Ben B. and Jake Bowers. 2008. "Covariate balance in simple,
stratified and clustered comparative studies". \emph{Statistical
Science} 23(2):219--236.

Bowers, Jake and Mark Fredrickson and Ben
Hansen. 2010. "RItools:Randomization Inference Tools". R package
version 0.1-11.

Moore, Ryan T. 2012. "Multivariate Continuous Blocking to Improve Political Science
Experiments". \emph{Political Analysis}, 20(4):460--479, Autumn.
}

\author{Ryan T. Moore}

\seealso{
  \code{\link{assignment}}.
}

\examples{
data(x100)
b <- block(x100, groups = "g", id.vars = "id", block.vars = c("b1", "b2"))
a <- assignment(b)
axb <- assg2xBalance(a, x100, id.var = "id", bal.vars = c("b1", "b2"))
axb

# axb is a list with 4 elements (one for each of 3 groups, plus one for 'Overall')
}

\keyword{design}
\keyword{multivariate}

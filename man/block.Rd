\name{block}
\alias{block}
\title{Block units into homogeneous experimental blocks}
\description{
  Block units into experimental blocks, with one unit per treatment
  condition.  Blocking begins by creating a measure of
  multivariate distance between all possible pairs of units.  Maximum,
  minimum, or an allowable range of differences between units on one
  variable can be set.
}
\usage{
block(data, vcov.data = NULL, groups = NULL, n.tr = 2, id.vars,
  block.vars = NULL, algorithm = "optGreedy", distance = "mahalanobis",
  weight = NULL, row.sort = NULL, level.two = FALSE, valid.var = NULL,
  valid.range = NULL, seed, verbose = FALSE, ...)
}

\arguments{
  \item{data}{a dataframe or matrix, with units in rows and variables in
    columns.}
  \item{vcov.data}{an optional matrix of data used to estimate the
    variance-covariance matrix for calculating multivariate distance.}
  \item{groups}{an optional column name from \code{data}, specifying subgroups
    within which blocking occurs.}
  \item{n.tr}{the number of treatment conditions per block.}
  \item{id.vars}{a required string or vector of two strings specifying which
    column(s) of \code{data} contain identifying information.}
  \item{block.vars}{an optional string or vector of strings specifying which
    column(s) of \code{data} contain the blocking variables.}
  \item{algorithm}{a string specifying the blocking algorithm.
    \code{"optGreedy"}, \code{"naiveGreedy"}, \code{"randGreedy"}, and
    \code{"sortGreedy"} algorithms are currently available.  See Details
    for more information.}
  \item{distance}{either a) a string defining how the multivariate
    distance used for blocking is calculated (options include
    \code{"mahalanobis"}, \code{"mcd"}, and \code{"mve"}), or b) a
    user-defined $k$-by-$k$ matrix, where $k$ is the number of rows 
in
    \code{data}.}
  \item{weight}{either a vector of length equal to the number of blocking 
    variables or a square matrix with dimensions equal to the number of 
    blocking variables used to explicitly weight blocking variables.}
  \item{row.sort}{an optional vector of integers from 1 to
    \code{nrow(data)} used to sort the rows of data when
    \code{algorithm = sortGreedy}.}
  \item{level.two}{a logical defining the level of blocking.}
  \item{valid.var}{an optional string defining a variable on which
    units in the same block must fall within the range defined by
    \code{valid.range}.} 
  \item{valid.range}{an optional vector defining the range of
    \code{valid.var} within which units in the same block must fall.}
  \item{seed}{an optional integer value for the random seed set in
    \code{cov.rob}, used to calculate measures of the
    variance-covariance matrix robust to outliers.}
  \item{verbose}{a logical specifying whether \code{groups} names and 
    block numbers are printed as blocks are created.}
  \item{\dots}{additional arguments passed to \code{cov.rob}.}
}
\details{
If \code{vcov.data = NULL}, then \code{block} calculates the
variance-covariance matrix using the \code{block.vars} from
\code{data}. 

If \code{groups} is not user-specified, \code{block} temporarily creates
a variable in \code{data} called \code{"groups"}, which takes the value
1 for every unit.

Where possible, one unit is assigned to each condition in each block.
If there are fewer available units than treatment conditions, available
units are used.

If \code{n.tr} $> 2$, then the \code{optGreedy} algorithm finds the best
possible pair match, then the best match to either member of the pair,
then the best match to any member of the triple, \ldots.  Other
algorithms proceed similarly.

An example of \code{id.vars} is \code{id.vars = c("id", "id2")}.  If
two-level blocking is selected, \code{id.vars} should be ordered 
(unit id, subunit id).  See details for \code{level.two} below for more 
information.

If \code{block.vars = NULL}, then all variables in \code{data} except
the \code{id.vars} are taken as blocking variables.  E.g.,
\code{block.vars = c("b1", "b2")}. 

\code{"optGreedy"} calls an optimal-greedy algorithm, sequentially
finding the best match in the entire dataset; \code{"naiveGreedy"} finds
the best match proceeding down the dataset from the first unit to the
last; \code{"randGreedy"} randomly selects a unit, finds its best match,
and repeats; \code{"sortGreedy"} resorts the dataset according to
\code{row.sort}, then implements a \code{naiveGreedy} algorithm.
  
The \code{optGreedy} algorithm breaks ties by randomly selecting one of
the minimum-distance pairs.  The \code{naiveGreedy}, \code{sortGreedy},
and \code{randGreedy} algorithms break ties by randomly selecting one of
the minimum-distance matches to the particular unit in question.

The \code{distance = "mcd"} and \code{distance = "mve"} options call
\code{cov.rob} to calculate measures of multivariate spread robust to
outliers.  The \code{distance = "mcd"} option calculates the Minimum
Covariance Determinant estimate; the \code{distance = "mve"} option
calculates the Minimum Volume Ellipsoid estimate.

A user-specified distance matrix must have diagonals equal to 0,
indicating zero distance between a unit and itself.  Only the lower
triangle of the matrix is used.

If \code{weight} is a vector, then it is used as the diagonal of a square 
weighting matrix with non-diagonal elements equal to zero.  The
weighting is done by using as the Mahalanobis distance scaling matrix
$((((chol(Sigma))')^\{-1\})'W((chol(Sigma))')^\{-1\})^\{-1\}$, where 
$chol(Sigma)$ is the Cholesky decomposition of the usual variance-covariance 
matrix and $W$ is the weighting matrix.  Differences should be smaller 
on covariates given higher weights.

If \code{level.two = TRUE}, then the best subunit block-matches in
different units are found.  E.g., provinces could be matched based on
the most similar cities within them.  All subunits in the data should
have unique names.  Thus, if subunits are numbered 1 to (number of 
subunits in unit) within each unit, then they should
be renumbered, e.g., 1 to (total number of subunits in all units). 

An example of a variable restriction is \code{valid.var="b2"},
\code{valid.range = c(10,50)}, which requires that units in the same
block be at least 10 units apart, but no more than 50 units apart, on
variable \code{"b2"}.
}

\value{ A list with elements
  \item{blocks}{a list of dataframes, each containing a group's blocked
    units.  If there are two treatment conditions, then the last column
    of each dataframe displays the multivariate distance between the two
    units.  If there are more than two treatment conditions, then the
    last column of each dataframe displays the largest of the
     multivariate distances between all possible pairs in the block.}
  \item{level.two}{a logical indicating whether \code{level.two =
      TRUE}.}
  \item{call}{the orginal call to \code{block}.}
}

\references{
King, Gary, Emmanuela Gakidou, Nirmala Ravishankar, Ryan T. Moore, Jason
Lakin, Manett Vargas, Martha Mar\'ia T\'ellez-Rojo and Juan Eugenio
Hern\'andez \'Avila and Mauricio Hern\'andez \'Avila and H\'ector
Hern\'andez Llamas.  2007.  "A 'Politically Robust' Experimental Design
for Public Policy Evaluation, with Application to the Mexican Universal
Health Insurance Program". \emph{Journal of Policy Analysis and
Management} 26(3):  479-509.
}

\author{Ryan T. Moore}

\seealso{\code{\link{assignment}}, \code{\link{diagnose}}}

\examples{
data(x100)
out <- block(x100, groups = "g", n.tr = 2, id.vars = c("id"), block.vars
             = c("b1", "b2"), algorithm="optGreedy", distance =
             "mahalanobis", level.two = FALSE, valid.var = "b1",
             valid.range = c(0,500), verbose = TRUE)
## out$blocks contains 3 data frames

## To illustrate two-level blocking, with multiple level two units per
##  level one unit:
for(i in (1:nrow(x100))){if(even(i)){x100$id[i] <- x100$id[i-1]}}

out <- block(x100, groups = "g", n.tr = 2, id.vars = c("id", "id2"),
              block.vars = c("b1", "b2"), algorithm="optGreedy",
              distance = "mahalanobis", level.two = TRUE, valid.var =
              "b1", valid.range = c(0,500), verbose = TRUE) 
}

\keyword{design}
\keyword{multivariate}

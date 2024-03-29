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
  weight = NULL, optfactor = 10^7, row.sort = NULL, level.two = FALSE, 
  valid.var = NULL, valid.range = NULL, seed.dist, namesCol = NULL, verbose = FALSE, ...)
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
    column(s) of \code{data} contain the numeric blocking variables.}
  \item{algorithm}{a string specifying the blocking algorithm.
    \code{"optGreedy"}, \code{"optimal"}, \code{"naiveGreedy"}, \code{"randGreedy"}, and
    \code{"sortGreedy"} algorithms are currently available.  See Details
    for more information.}
  \item{distance}{either a) a string defining how the multivariate
    distance used for blocking is calculated (options include
    \code{"mahalanobis"}, \code{"mcd"}, \code{"mve"}, and \code{"euclidean"}), or b) a
    user-defined $k$-by-$k$ matrix of distances, where $k$ is the number of rows in 
    \code{data}.}
  \item{weight}{either a vector of length equal to the number of blocking 
    variables or a square matrix with dimensions equal to the number of 
    blocking variables used to explicitly weight blocking variables.}
  \item{optfactor}{a number by which distances are multiplied then divided when \code{algorithm = "optimal"}.}
  \item{row.sort}{an optional vector of integers from 1 to
    \code{nrow(data)} used to sort the rows of data when
    \code{algorithm = "sortGreedy"}.}
  \item{level.two}{a logical defining the level of blocking.}
  \item{valid.var}{an optional string defining a variable on which
    units in the same block must fall within the range defined by
    \code{valid.range}.} 
  \item{valid.range}{an optional vector defining the range of
    \code{valid.var} within which units in the same block must fall.}
  \item{seed.dist}{an optional integer value for the random seed set in
    \code{cov.rob}, used to calculate measures of the
    variance-covariance matrix robust to outliers.}
  \item{namesCol}{an optional vector of column names for the output
    table.}
  \item{verbose}{a logical specifying whether \code{groups} names and 
    block numbers are printed as blocks are created.}
  \item{\dots}{additional arguments passed to \code{cov.rob}.}
}
\details{
If \code{vcov.data = NULL}, then \code{block} calculates the
variance-covariance matrix using the \code{block.vars} from
\code{data}. 

If \code{groups} is not user-specified, \code{block} temporarily creates
a variable in \code{data} called \code{groups}, which takes the value
1 for every unit.

Where possible, one unit is assigned to each condition in each block.
If there are fewer available units than treatment conditions, available
units are used.

If \code{n.tr} $> 2$, then the \code{optGreedy} algorithm finds the best
possible pair match, then the best match to either member of the pair,
then the best match to any member of the triple, etc.  After finding the best pair match to a given unit, the other greedy algorithms proceed by finding the third, fourth, etc. best match to that given unit.

An example of \code{id.vars} is \code{id.vars = c("id", "id2")}.  If
two-level blocking is selected, \code{id.vars} should be ordered 
(unit id, subunit id).  See details for \code{level.two} below for more 
information.

If \code{block.vars = NULL}, then all variables in \code{data} except
the \code{id.vars} are taken as blocking variables.  E.g.,
\code{block.vars = c("b1", "b2")}.  

The algorithm \code{optGreedy} calls an optimal-greedy algorithm, repeatedly
finding the best remaining match in the entire dataset;
\code{optimal} finds the set of blocks that minimizes the sum of the
distances in all blocks; \code{naiveGreedy} finds the best match
proceeding down the dataset from the first unit to the last;
\code{randGreedy} randomly selects a unit, finds its best match, and
repeats; \code{sortGreedy} resorts the dataset according to
\code{row.sort}, then implements the \code{naiveGreedy} algorithm.
  
The \code{optGreedy} algorithm breaks ties by randomly selecting one of
the minimum-distance pairs.  The \code{naiveGreedy}, \code{sortGreedy},
and \code{randGreedy} algorithms break ties by randomly selecting one of
the minimum-distance matches to the particular unit in question.

As of version 0.5-1, blocking is done in C for all algorithms except
\code{optimal} (see following paragraphs for more details on the
\code{optimal} algorithm implementation).

The \code{optimal} algorithm uses two functions from the
\pkg{nbpMatching} package: \code{distancematrix} prepares a distance
matrix for optimal blocking, and \code{nonbimatch} performs the
optimal blocking by minimizing the sum of distances in blocks.
\code{nonbimatch}, and thus the \code{block} algorithm
\code{optimal}, requires that \code{n.tr = 2}.

Because \code{distancematrix} takes the integer \code{floor} of the
distances, and one may want much finer precision, the multivariate
distances calculated within \code{block} are multiplied by
\code{optfactor} prior to optimal blocking.  Then
\code{distancematrix} prepares the resulting distance matrix, and
\code{nonbimatch} is called on the output.  The distances are then
untransformed by dividing by \code{optfactor} before being returned by
\code{block}.

The choice of \code{optfactor} can determine whether the Fortran code
can allocate enough memory to solve the optimization problem.  For
example, blocking the first 14 units of \code{\link{x100}} by executing
\code{block(x100[1:14, ], id.vars = "id", block.vars = c("b1", "b2"),
algorithm = "optimal", optfactor = 10^8)} fails for Fortran memory
reasons, while the same code with \code{optfactor = 10^5} runs
successfully.  Smaller values of \code{optfactor} imply easier
computation, but less precision.

Most of the algorithms in \code{block} make prohibited blockings by
using a distance of \code{Inf}.  However, the optimal algorithm calls
\code{Fortran} code from \pkg{nbpMatching} and requires integers.
Thus, a distance of \code{99999*max(dist.mat)} is used to effectively
prohibit blockings.  This follows the procedure demonstrated in the
example of \code{help(nonbimatch)}.

In order to enable comparisons of block-quality across groups, when
\code{distance} is a string, $Sigma$ is calculated using units from
all groups.

The \code{distance = "mcd"} and \code{distance = "mve"} options call
\code{cov.rob} to calculate measures of multivariate spread robust to
outliers.  The \code{distance = "mcd"} option calculates the Minimum
Covariance Determinant estimate (Rousseeuw 1985); the \code{distance = "mve"} option
calculates the Minimum Volume Ellipsoid estimate (Rousseeuw and van Zomeren 1990).  When \code{distance = "mcd"}, the interquartile range on blocking variables should not be zero.

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
subunits in unit) within each unit, then they should be renumbered,
e.g., 1 to (total number of subunits in all units).  \code{level.two}
blocking is not currently implemented for \code{algorithm =
"optimal"}.  Units with no blocked subunit are put into their own
blocks.  However, unblocked subunits within a unit that does have a
blocked subunit are not put into their own blocks.

An example of a variable restriction is \code{valid.var = "b2"},
\code{valid.range = c(10,50)}, which requires that units in the same
block be at least 10 units apart, but no more than 50 units apart, on
variable \code{"b2"}.  As of version 0.5-3, variable restrictions are
implemented in all algorithms except \code{optimal}.  Note that
employing a variable restriction may result in fewer than the maximum possible number of blocks.  See 
\url{https://www.ryantmoore.org/html/software.blockTools.html} for details.

If \code{namesCol = NULL}, then \dQuote{Unit 1},
\dQuote{Unit 2}, \ldots are used.  If \code{level.two = FALSE}, then \code{namesCol} should be of length \code{n.tr}; if \code{level.two = TRUE}, then \code{namesCol} should be of length 2*\code{n.tr}, and in the order shown in the example below.
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

Moore, Ryan T. 2012. "Multivariate Continuous Blocking to Improve Political Science
Experiments." \emph{Political Analysis} 20(4):460-479.

Rousseeuw, Peter J. 1985. "Multivariate Estimation with High Breakdown Point". \emph{Mathematical Statistics and Applications} 8:283-297.

Rousseeuw, Peter J. and Bert C. van Zomeren. 1990. "Unmasking Multivariate Outliers and Leverage Points". \emph{Journal of the American Statistical Association} 85(411):633-639.
}

\author{Ryan T. Moore \email{rtm@american.edu} and Keith Schnakenberg \email{keith.schnakenberg@gmail.com}}

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
for(i in (1:nrow(x100))){if((i \%\% 2) == 0){x100$id[i] <- x100$id[i-1]}}

out2 <- block(x100, groups = "g", n.tr = 2, id.vars = c("id", "id2"),
              block.vars = c("b1", "b2"), algorithm="optGreedy",
              distance = "mahalanobis", level.two = TRUE, valid.var =
              "b1", valid.range = c(0,500), namesCol = c("State 1", "City 1", 
	      "State 2", "City 2"), verbose = TRUE) 
}

\keyword{design}
\keyword{multivariate}

\name{seqblock}
\alias{seqblock}
\title{Sequential assignment of unit(s) into experimental conditions using covariates}
\description{
Sequentially assign units into experimental conditions. Blocking begins by creating a measure of multivariate distance between a \emph{current} unit and one or multiple \emph{prior}, already-assigned unit(s). Then, average distance between current unit and each treatment condition is calculated, and random assignment is biased toward conditions more dissimilar to current unit. Argument values can be specified either as argument to the function, or via a query. The query directly asks the user to identify the blocking variables and to input, one-by-one, each unit's variable values.
}

\usage{
seqblock(object = NULL, id.vars, id.vals, exact.vars = NULL, exact.vals = NULL, 
  exact.restr = NULL, exact.alg = "single", covar.vars = NULL, covar.vals = NULL, 
  covar.restr = NULL, covars.ord = NULL, n.tr = 2, tr.names = NULL, assg.prob = NULL,
  seed = NULL, seed.dist, assg.prob.stat = NULL, trim = NULL, assg.prob.method = NULL,
  assg.prob.kfac = NULL, distance = NULL, file.name = NULL, query = FALSE, 
  verbose = TRUE, ...)
}

\arguments{
  \item{object}{a character string giving the file name of a \code{.RData} file containing a list output from the \code{seqblock} function which contains at least one previously assigned unit.}
  \item{id.vars}{a string or vector of strings specifying the name of the identifying variable(s); if \code{query = FALSE} and the object argument is not given, then the \code{id.vars} argument is required.}
  \item{id.vals}{a vector of ID values for every unit being assigned to a treatment group; those are corresponding to the \code{id.vars}.}
  \item{exact.vars}{a string or vector of strings containing the names of each of the exact blocking variables.}
  \item{exact.vals}{a vector containing the unit's values on each of the exact blocking variables.}
  \item{exact.restr}{a list object containing the restricted values that the exact blocking variables can take on. Thus the first element of \code{exact.restr} is a vector containing all of the possible values that the first exact blocking variable (see \code{exact.vars} above) can take on; the second element is a vector containing all of the possible values for the second exact blocking variable; and so on.}
  \item{exact.alg}{a string specifying the blocking algorithm. Currently the only acceptable value is \code{"single"}. This algorithm creates a variable with a unique level for every possible combination of the values in all of the exact variables. See Details section below.}
  \item{covar.vars}{a string or vector of strings containing the names of each of the non-exact blocking variables.}
  \item{covar.vals}{a vector containing the unit's values on each of the non-exact blocking variables.}
  \item{covar.restr}{a list object containing the restricted values that the non-exact blocking variables can take on. Thus the first element of \code{covar.restr} is a vector containing all of the possible values that the first non-exact blocking variable (see \code{covar.vars} above) can take on; the second element is a vector containing all of the possible values for the second non-exact blocking variable; and so on.}
  \item{covars.ord}{a string or vector of strings containing the name of the non-exact blocking variables ordered so that the highest priority covariate comes first, followed by the second highest priority covariate, then the third, etc.}
  \item{n.tr}{the number of treatment groups. If not specified, this defaults to \code{n.tr = 2}.}
  \item{tr.names}{a string or vector of strings containing the names of the different treatment groups.}
  \item{assg.prob}{a numeric vector containing the probabilities that a unit will be assigned to the treatment groups; this vector should sum to 1.}
  \item{seed}{an optional integer value for the random seed, which is used when assigning units to treatment groups.}
  \item{seed.dist}{an optional integer value for the random seed set in \code{cov.rob}, used to calculate measures of the variance-covariance matrix robust to outliers.}
  \item{assg.prob.stat}{a string specifying which assignment probability summary statistic to use; valid values are \code{mean}, \code{median}, and \code{trimmean}. If not specified, this defaults to \code{assg.prob.stat = "mean"}.}
  \item{trim}{a numeric value specifying what proportion of the observations are to be dropped from each tail when the assignment probability summary statistic (\code{assg.prob.stat}) is set equal to \code{trimmean}. Blocks on each tail of the distribution are dropped before the mean is calculated. If not specified, this defaults to \code{trim = 0.1}.}
  \item{assg.prob.method}{a string specifying which algorithm should be used when assigning treatment probabilities. Acceptable values are \code{ktimes}, \code{fixed}, \code{prop}, \code{prop2}, and \code{wprop}. If not specified, this defaults to \code{assg.prob.method = "ktimes"}.}
  \item{assg.prob.kfac}{a numeric value for \code{k}, the factor by which the most likely experimental condition will be multiplied relative to the other conditions. If not specified, this defaults to \code{assg.prob.kfac = 2}.}
  \item{distance}{a string specifying how the multivariate distance used for blocking covariates are calculated. If not specified, this defaults to \code{distance = "mahalanobis"}.}
  \item{file.name}{a string containing the name of the file that one would like the output to be written to. Ideally this file name should have the extension .RData.}
  \item{query}{a logical stating whether the console should ask the user questions to input the data and assign a treatment condition. If not specified, this defaults to \code{query = FALSE}.}
  \item{verbose}{a logical stating whether the function should print the name of the output file, the current working directory, the treatment group that the most recent unit was assigned to, and the dataframe \code{x} returned by the function as part of the \code{bdata} list. If not specified, this defaults to \code{verbose = TRUE}.}
  \item{\dots}{additional arguments.}
}

\details{
The \code{seqblock} function's code is primarily divided into two parts: the first half deals with instances, in which the unit being assigned is the first unit in a given study to receive an assignment; the second half addresses subsequent units that are assigned after at least one first assignment has already been made. If the \code{object} argument is left as \code{NULL}, the function will run the first half; if the \code{object} argument is specified, the second part will be executed. When \code{object = NULL}, the researcher has no past file from which to pull variable names and past data; this corresponds to the case when the unit being assigned is the first one. If the researcher does specify \code{object}, it implies the user is drawing data from a past file, which means this is not the first unit in the study to be assigned to a treatment.

However, the function can be called for subsequent units even when \code{object} is not specified. By setting \code{query = TRUE}, the console will ask the researcher whether this is the first unit to be assigned in the study. Based on the researcher's response, it will decide which part of the code to run.

If the \code{object} and \code{file.name} arguments are set to the same value, then \code{seqblock} overwrites the specified file with a new file, which now contains both the previously-assigned units and the newly-assigned unit. To create a new file when a new unit is assigned, use a new \code{file.name}. 

The \code{single} algorithm (see \code{exact.alg} in the Arguments section above) creates a variable that has a unique level for every possible combination of the exact variables. As an example, say that there were 3 exact blocking variables: \emph{party} (Democrat, Republican); \emph{region} (North, South); and \emph{education} (HS, NHS). The \code{single} algorithm creates one level for units with the following values: Democrat-North-HS. It would create another level for Democrat-North-NHS; a third level for Republican-North-HS; and so forth, until every possible combination of these 3 variables has its own level. Thus if there are \eqn{k} exact blocking variables and each exact blocking variable has \eqn{q_{i}} values it can take on, then there are a total of \eqn{\prod_{1}^{k} q_{i}} levels created.

The \code{distance = "mcd"} and \code{distance = "mve"} options call \code{cov.rob} to calculate measures of multivariate spread robust to outliers. The \code{distance = "mcd"} option calculates the Minimum Covariance Determinant estimate (Rousseeuw 1985); the \code{distance = "mve"} option calculates the Minimum Volume Ellipsoid estimate (Rousseeuw and van Zomeren 1990). When \code{distance = "mcd"}, the interquartile range on blocking variables should not be zero. The \code{distance = "euclidean"} option calculates the Euclidean distance between the new unit and the previously-assigned units. The default \code{distance = "mahalanobis"} option calculates the Mahalanobis distance.
}

\value{ A list (called \code{bdata}) with elements
  \item{x}{a dataframe containing the names and values for the different ID and blocking variables, as well as each unit's initial treatment assignment.}
  \item{nid}{a string or vector of strings containing the name(s) of the ID variable(s).}
  \item{nex}{a string or vector of strings containing the name(s) of the exact blocking variable(s).}
  \item{ncv}{a string or vector of strings containing the name(s) of the non-exact blocking variable(s).}
  \item{rex}{a list of the restricted values of the exact blocking variables.}
  \item{rcv}{a list of the restricted values of the non-exact blocking variables.}
  \item{ocv}{a vector of the order of the non-exact blocking variables.}
  \item{trn}{a string or vector of strings containing the name(s) of the different treatment groups.}
  \item{apstat}{a string specifying the assignment probability summary statistic that was used.}
  \item{mtrim}{a numeric value specifying the proportion of observations to be dropped when the assignment probability statistic takes on the value \code{"trimmean"}.}
  \item{apmeth}{a string specifying the assignment probability algorithm that was used.}
  \item{kfac}{the assignment probability \emph{kfactor}; see \emph{assg.prob.kfac} in the Arguments section above.}
  \item{assgpr}{a vector of assignment probabilities to each treatment group.}
  \item{distance}{a string specifying how the multivarite distance used for blocking was calculated. }
  \item{trd}{a list with the length equal to the number of previously assigned treatment conditions; each object in the list contains a vector of the distance between each unit in one treatment group and the new unit. This will be \code{NULL} when there are no non-exact blocking variables.}
  \item{tr.sort}{a string vector of treatment conditions, sorted from the largest to the smallest. Set to \code{NULL} when there are no non-exact blocking variables.}
  \item{p}{a vector of assignment probabilities to each treatment group used in assigning a treatment condition to the new unit.}
  \item{distance}{a string specifying how the multivarite distance used for blocking is calculated}
  \item{trcount}{a table containing the counts for each experimental/treatment conditions.}
  \item{datetime}{the date and time at which each unit was assigned their treatment group.}
  \item{orig}{a dataframe containing the names and values for the different id and blocking variables, as well as each unit's treatment assignment.}
}

\references{
Moore, Ryan T. and Sally A. Moore. 2013. "Blocking for Sequential Political Experiments." \emph{Political Analysis} 21(4):507-523. 

Moore, Ryan T. 2012. "Multivariate Continuous Blocking to Improve Political Science
Experiments." \emph{Political Analysis} 20(4):460-479.

Rousseeuw, Peter J. 1985. "Multivariate Estimation with High Breakdown Point". \emph{Mathematical Statistics and Applications} 8:283-297.

Rousseeuw, Peter J. and Bert C. van Zomeren. 1990. "Unmasking Multivariate Outliers and Leverage Points". \emph{Journal of the American Statistical Association} 85(411):633-639.
}

\author{Ryan T. Moore \email{rtm@wustl.edu}, Tommy Carroll \email{tcarroll22@wustl.edu}, Jonathan Homola \email{homola@wustl.edu} and Jeong Hyun Kim \email{jeonghyun.kim@wustl.edu}}

\seealso{\code{\link{assignment}}, \code{\link{block}}}

\examples{
## Assign first unit (assume a 25 year old member of the Republican Party) to a treatment group.
## Save the results in file "sdata.RData":
## seqblock(query = FALSE, id.vars = "ID", id.vals = 001, exact.vars = "party", 
##   exact.vals = "Republican", covar.vars = "age", covar.vals = 25, file.name = "sdata.RData")

## Assign next unit (age 30, Democratic Party):
## seqblock(query = FALSE, object = "sdata.RData", id.vals = 002, exact.vals = "Democrat", 
##   covar.vars = "age", covar.vals = 30, file.name = "sdata.RData")
}
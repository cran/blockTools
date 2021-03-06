\name{invertRIconfInt}
\alias{invertRIconfInt}

\title{Calculate treatment effect confidence intervals by inverting the randomization test}

\description{
Using an output object from \code{seqblock} or any other matrix or dataframe that includes a treatment and an outcome variable for multiple units, as well as blocking and non-blocking variables for the respective unit(s), \code{invertRIconfInt} calculates treatment effect confidence intervals by inverting the randomization inference test.}

\usage{
invertRIconfInt(dat, outcome.var, tr.var, tau.abs.min = -1, tau.abs.max = 1, 
  tau.length = 10, n.sb.p = 100, id.vars, id.vals, exact.vars = NULL, exact.vals = NULL,
  exact.restr = NULL, exact.alg = "single", covar.vars = NULL, covar.vals = NULL, 
  covar.restr = NULL, covars.ord = NULL, n.tr = 2, tr.names = NULL, assg.prob = NULL,
  seed = NULL, seed.dist, assg.prob.stat = NULL, trim = NULL, assg.prob.method = NULL, 
  assg.prob.kfac = NULL, distance = "mahalanobis", file.name = "sbout.RData", 
  query = FALSE, verbose = TRUE)
}

\arguments{
  \item{dat}{a matrix or dataframe containing the names and values of the different blocking and non-blocking variables, as well as each unit's treatment assignment and outcome}
  \item{outcome.var}{a string specifying the name of the outcome variable}
  \item{tr.var}{a string specifying the name of the treatment variable}
  \item{tau.abs.min}{lower bound of the range across which the confidence intervals will be computed}
  \item{tau.abs.max}{upper bound of the range across which the confidence intervals will be computed}
  \item{tau.length}{the number of (evenly spaced) possible treatment effects across the range specified by \code{tau.abs.min} and \code{tau.abs.max} for which location inside or outside the confidence intervals will be computed}
  \item{n.sb.p}{the number of times that sequential blocking will be performed on the dataset}
  \item{id.vars}{see the \code{seqblock} documentation}
  \item{id.vals}{see the \code{seqblock} documentation}
  \item{exact.vars}{see the \code{seqblock} documentation}
  \item{exact.vals}{see the \code{seqblock} documentation}
  \item{exact.restr}{see the \code{seqblock} documentation}
  \item{exact.alg}{see the \code{seqblock} documentation}
  \item{covar.vars}{see the \code{seqblock} documentation}
  \item{covar.vals}{see the \code{seqblock} documentation}
  \item{covar.restr}{see the \code{seqblock} documentation}
  \item{covars.ord}{see the \code{seqblock} documentation}
  \item{n.tr}{see the \code{seqblock} documentation}
  \item{tr.names}{see the \code{seqblock} documentation}
  \item{assg.prob}{see the \code{seqblock} documentation}
  \item{seed}{see the \code{seqblock} documentation}
  \item{seed.dist}{see the \code{seqblock} documentation}
  \item{assg.prob.stat}{see the \code{seqblock} documentation}
  \item{trim}{see the \code{seqblock} documentation}
  \item{assg.prob.method}{see the \code{seqblock} documentation}
  \item{assg.prob.kfac}{see the \code{seqblock} documentation}
  \item{distance}{see the \code{seqblock} documentation}
  \item{file.name}{see the \code{seqblock} documentation}
  \item{query}{see the \code{seqblock} documentation}
  \item{verbose}{see the \code{seqblock} documentation}
}

\details{
\code{invertRIconfInt} takes a data matrix (or data frame) containing names and values of different blocking and non-blocking variables, as well as each unit's treatment assignment and outcome as input and returns a list of treatment effect confidence intervals.

Apart from specifying the treatment and outcome variable, the user can set all other arguments to \code{seqblock} when running \code{invertRIconfInt}. The function will then calculate the confidence intervals by employing a method described in Ho and Imai (2006), which inverts Fisher's exact test. The resulting confidence intervals are distribution-free, nonparametric and have accurate coverage probabilities.
}

\value{
A list with elements
  \item{ci95}{vector of treatment effects within the 95\% confidence interval}
  \item{ci90}{vector of treatment effects within the 90\% confidence interval}
  \item{ci80}{vector of treatment effects within the 80\% confidence interval}
}

\references{
Moore, Ryan T. and Sally A. Moore. 2013. "Blocking for Sequential Political Experiments." \emph{Political Analysis} 21(4): 507-523. 

Ho, Daniel E., and Kosuke Imai. 2006. "Randomization inference with natural experiments: An analysis of ballot effects in the 2003 California recall election." \emph{Journal of the American Statistical Association} 101(475): 888-900.
}

\author{Ryan T. Moore \email{rtm@wustl.edu} and Jonathan Homola \email{homola@wustl.edu}}

\seealso{\code{\link{seqblock}}}
\examples{
## Create an example data matrix with 50 observations that contains an ID variable, 
##   a dummy variable indicating gender, an age variable (between 18 and 55), a 
##   treatment variable and an outcome variable (between 15 and 20). 
## id <- seq(1,50,1)
## gender <- sample(c(1,2),50,replace=T)
## age <- sample(seq(18,55,1),50,replace=T)
## treat <- sample(c(1,2),50,replace=T)
## out <- treat + sample(seq(15,20,1),50,replace=T)
## data <- cbind(id, gender, age, out, treat)

## Check summary statistics for the created data
## aggregate(out~treat, data, mean)

## Run invertRIconfInt()
## invertRIconfInt(data, outcome.var="out", tr.var="treat", tau.abs.min = -3, 
##  tau.abs.max = 3, id.vars = "id", id.vals = "id", exact.vars = c("gender", "age"), 
##  exact.vals = c("gender", "age"))
}

\keyword{design}
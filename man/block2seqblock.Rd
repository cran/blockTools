\name{block2seqblock}
\alias{block2seqblock}

\title{Prepare prior nonsequential assignments for subsequent sequential assignments}

\description{Converts output objects from the \code{block} and \code{assignment} functions
into an object in the format of one output by the \code{seqblock} function. This allows the
user to block and assign multiple units at the beginning of an experiment (using \code{block} and
\code{assignment}) and then sequentially block and assign more units to the experiment over time
(using \code{seqblock}).}

\usage{
block2seqblock(block.obj, assg.obj, data, exact.restr = NULL, covar.restr = NULL, 
  covar.order = NULL, trn = NULL, apstat = "mean", mtrim = 0.1, apmeth = "ktimes", 
  kfac = 2, assgpr = c(0.5, 0.5), distance = NULL, datetime = NULL, orig, seed = NULL,
  file.name = "sbout.RData", verbose = FALSE)
}

\arguments{
  \item{block.obj}{an output object from \code{block}, or a user-specified block object}
  \item{assg.obj}{an output object from \code{assignment}, or a user-specified assignment object}
  \item{data}{a matrix or dataframe containing the original data used to block the units in the study}
  \item{exact.restr}{a list object containing the restricted values that the exact blocking variables can take on. Thus the first element of \code{exact.restr} is a vector containing all of the possible values that the first exact blocking variable can take on; the second element is a vector containing all of the possible values for the second exact blocking variable; and so on}
  \item{covar.restr}{a list object containing the restricted values that the non-exact blocking variables can take on. Thus the first element of \code{covar.restr} is a vector containing all of the possible values that the first non-exact blocking variable can take on; the second element is a vector containing all of the possible values for the second non-exact blocking variable; and so on}
  \item{covar.order}{a string or vector of strings containing the name of the non-exact blocking variables ordered so that the highest priority covariate comes first, followed by the second highest priority covariate, then the third, etc.}
  \item{trn}{a string or vector of strings containing the names of the different treatment groups}
  \item{apstat}{a string specifying the assignment probability summary statistic that was used}
  \item{mtrim}{a numeric value specifying the proportion of observations to be dropped when the assignment probability statistic takes on the value \code{"trimmean"}.}
  \item{apmeth}{a string specifying the assignment probability algorithm that was used.}
  \item{kfac}{the assignment probability \emph{kfactor}; see \emph{assg.prob.kfac} in the Arguments section above} 
  \item{assgpr}{a vector of assignment probabilities to each treatment group}
  \item{distance}{a string specifying how the multivarite distance used for blocking is calculated} 
  \item{datetime}{the date and time that the units were assigned to the treatment group; by default this is set to be a vector of NA; however the user could also specify a specific datetime and all of the units from the block object will be given the same datetime stamp} 
  \item{orig}{a dataframe containing the names and values for the different id and blocking variables, as well as each unit's initial treatment assignment}
  \item{seed}{an optional integer value for the random seed set which is used when assigning units to treatment groups} 
  \item{file.name}{a string containing the name of the file that one would like the output to be written to. Ideally this file name should have the extension .RData}
  \item{verbose}{a logical stating whether the function should print the name of the output file, the current working directory, and the dataframe \code{x} returned by the function as part of the \code{bdata} list}
}
\details{
The function converts data from a blocked experiment into a form allowing subsequent sequential blocking.  Minimally, the user sets only the arguments \code{block.obj}, \code{assg.obj} and \code{data}.  Then, \code{block2seqblock} uses the call to \code{block}, the assignment object, and the original data to create an object that is ready to be input into \code{seqblock}.

If the user explicitly specifies \code{groups}, \code{id.vars} and \code{block.vars} in the initial \code{block} function that is used to create the \code{block.obj}, \code{block2seqblock} will order the variables in the output it produces according to the order specified in the initial \code{block} function call. If the user does not explicitly specify the blocking variables in the \code{block} function call, \code{block2seqblock} will order the variables according to the order in the initial matrix or dataframe that was used to run the original \code{block} function.

As part of the function, variables that are of class \code{factor} in the original matrix or dataframe specified in \code{data}, will be converted into class \code{character}. 
 
The \code{trn} argument uses the \code{n.tr} argument from \code{block} to extract the names of the treatment variables. Most other arguments are set to default values that mirror those in the \code{seqblock} function.  One exception is the \code{datetime} argument, which defaults to a vector of \code{NA}'s instead of the current datetime.}

\value{A list (called \code{bdata}) with elements
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
  \item{distance}{a string specifying how the multivarite distance used for blocking is calculated}
  \item{trd}{a list with the length equal to the number of previously assigned treatment conditions; each object in the list contains a vector of the distance between each unit in one treatment group and the new unit.  Set to \code{NULL} when there are no non-exact blocking variables.}
  \item{tr.sort}{a string vector of treatment conditions, sorted from the largest to the smallest}
  \item{p}{a vector of assignment probabilities to each treatment group used in assigning a treatment condition to the new unit.}
  \item{trcount}{a table containing the counts for each experimental/treatment conditions.}
  \item{datetime}{the date and time that the user was assigned a treatment group.}
  \item{orig}{a dataframe containing the names and values for the different id and blocking variables, as well as each unit's initial treatment assignment.}
}

\author{Tommy Carroll \email{tcarroll22@wustl.edu}, Jonathan Homola \email{homola@wustl.edu}, and Ryan T. Moore \email{rtm@wustl.edu}}

\seealso{\code{\link{block}}, \code{\link{assignment}}, \code{\link{seqblock}}}
\examples{
# data(x100)
# out <- block(x100, n.tr = 2, id.vars = c("id"), block.vars = c("b1", "b2"), 
#   algorithm = "optGreedy", distance = "mahalanobis", valid.var = "b1", valid.range = c(0,500))
# assg.out <- assignment(out, seed = 123)
# b2sb <- block2seqblock(block.obj = out, assg.obj = assg.out, data = x100)
# sb <- seqblock("sbout.RData", id.vals = 1101, covar.vals = c(100, 200), file.name = "sb101.RData")
}

\keyword{design}
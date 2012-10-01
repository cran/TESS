\name{globalBiDe.likelihood}
\alias{globalBiDe.likelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{globalBiDe.likelihood: Probability density of a tree under a global, time-dependent birth-death process}
\description{globalBiDe.likelihood computes the probability of a reconstructed phylogenetic tree under time-dependent. The rates may be any positive function of time or a constant. Additionally, mass-extinction event can be provided and a random taxon sampling probability.}
\usage{
globalBiDe.likelihood(tree, lambda, mu, massExtinctionTimes = c(), massExtinctionSurvivalProbabilities = c(), samplingProbability = 1, MRCA = TRUE, log = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{The tree in 'phylo' format.}
  \item{lambda}{The speciation rate function or constant.}
  \item{mu}{The extinction rate function or constant.}
  \item{massExtinctionTimes}{The set of mass-extinction times after the start of the process.}
  \item{massExtinctionSurvivalProbabilities}{The set of survival probabilities for each speciation event. The set must have the same length as the set of mass-extinction times.}
  \item{samplingProbability}{The probability for a species to be included in the sample.}
  \item{MRCA}{Does the process start with the most recent common ancestor?}
  \item{log}{log-likelihood?}
}
\value{Returns the (log) probability of the tree.}
\references{Hoehna, Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2012, submitted }
\author{Sebastian Hoehna}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
newick <- "((A:0.25,B:0.25):0.75,(C:0.5,D:0.5):0.5);"
tree <- read.tree(text=newick)

l <- function(x) { if (x > 0.5 || x < 0.3) { return (1) } else { return (2) } }
e <- function(x) { if (x > 0.5 || x < 0.3) { return (0.95) } else { return (0.5) } }

globalBiDe.likelihood(tree,l,e,MRCA=TRUE,log=TRUE)

# constant rate likelihood function
globalBiDe.likelihood(tree,2.0,1.0,MRCA=TRUE,log=TRUE)

# constant rate pure birth likelihood function
globalBiDe.likelihood(tree,2.0,0.0,MRCA=TRUE,log=TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ likelihood }
%\keyword{ tree }
%\keyword{ birth-death }% __ONLY ONE__ keyword per line

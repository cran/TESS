\name{globalBiDe.likelihood}
\alias{globalBiDe.likelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{globalBiDe.likelihood: Probability density of a tree under a global, time-dependent birth-death process}
\description{globalBiDe.likelihood computes the probability of a reconstructed phylogenetic tree under time-dependent. The rates may be any positive function of time or a constant. Additionally, mass-extinction event can be provided and a random taxon sampling probability. You have several options for the start of the process (origin vs MRCA) and the condition of the process (time, survival or taxa; note that survival and taxa implicetily condition on the time too!). See equation (5) in the manuscript for more information. Note that constant rates lead to much faster computations.}
\usage{
globalBiDe.likelihood(tree, lambda, mu, massExtinctionTimes = c(), 
   massExtinctionSurvivalProbabilities = c(), samplingProbability = 1, 
   MRCA = TRUE, CONDITION = "survival", log = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{The tree in 'phylo' format.}
  \item{lambda}{The speciation rate function or constant.}
  \item{mu}{The extinction rate function or constant.}
  \item{massExtinctionTimes}{The set of mass-extinction times after the start of the process.}
  \item{massExtinctionSurvivalProbabilities}{The set of survival probabilities for each speciation event. The set must have the same length as the set of mass-extinction times.}
  \item{samplingProbability}{The probability for a species to be included in the sample.}
  \item{MRCA}{Does the process start with the most recent common ancestor? If not, the tree must have a root edge!}
  \item{CONDITION}{do we condition the process on time|survival|taxa??}
  \item{log}{log-likelihood?}
}
\value{Returns the (log) probability of the tree.}
\references{Hoehna, S.: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2012, submitted }
\author{Sebastian Hoehna}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(cettiidae)

l <- function(x) { if (x > 0.5 || x < 0.3) { return (1) } else { return (2) } }
e <- function(x) { if (x > 0.5 || x < 0.3) { return (0.95) } else { return (0.5) } }

globalBiDe.likelihood(cettiidae,l,e,MRCA=TRUE,log=TRUE)

# constant rate likelihood function
globalBiDe.likelihood(cettiidae,2.0,1.0,MRCA=TRUE,log=TRUE)

# constant rate pure birth likelihood function
globalBiDe.likelihood(cettiidae,2.0,0.0,MRCA=TRUE,log=TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }
%\keyword{ birth-death }% __ONLY ONE__ keyword per line

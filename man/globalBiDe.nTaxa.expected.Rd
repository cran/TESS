\name{globalBiDe.nTaxa.expected}
\alias{globalBiDe.nTaxa.expected}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{globalBiDe.nTaxa.expected: The expected number of taxa at present of a tree under a global, time-dependent birth-death process (E[N(T)])}
\description{globalBiDe.nTaxa.expected computes the expected number of taxa at the present time T (the process start at time s and times increases until the present) under time-dependent. The rates may be any positive function of time or a constant. Additionally, mass-extinction event can be provided and a random taxon sampling probability. You have several options for the start of the process (origin vs MRCA). See equation (??) in the manuscript for more information. Note that constant rates (scalars) lead to much faster computations.}
\usage{
globalBiDe.nTaxa.expected(t_low,t_high,lambda,mu,massExtinctionTimes=c(),
   massExtinctionSurvivalProbabilities=c(),samplingProbability=1.0,MRCA=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t_low}{The time when the process starts.}
  \item{t_high}{The time when the process end (e.g. the present).}
  \item{lambda}{The speciation rate function or constant.}
  \item{mu}{The extinction rate function or constant.}
  \item{massExtinctionTimes}{The set of mass-extinction times after the start of the process.}
  \item{massExtinctionSurvivalProbabilities}{The set of survival probabilities for each speciation event. The set must have the same length as the set of mass-extinction times.}
  \item{samplingProbability}{The probability for a species to be included in the sample.}
  \item{MRCA}{Does the process start with the most recent common ancestor? If not, the tree must have a root edge!}
}
\value{Returns the expected number of taxa.}
\references{Hoehna, S.: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2012, submitted }
\author{Sebastian Hoehna}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

l <- function(x) { if (x > 0.5 || x < 0.3) { return (1) } else { return (2) } }
e <- function(x) { if (x > 0.5 || x < 0.3) { return (0.95) } else { return (0.5) } }

globalBiDe.nTaxa.expected(0,5,l,e,MRCA=TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
%\keyword{ htest }
%\keyword{ birth-death }% __ONLY ONE__ keyword per line

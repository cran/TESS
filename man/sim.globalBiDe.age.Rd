\name{sim.globalBiDe.age}
\alias{sim.globalBiDe.age}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sim.globalBiDe.age: Simulate a reconstructed tree for a given age under a global, time-dependent birth-death process.}
\description{sim.globalBiDe.age simlates a reconstructed phylogenetic tree under a global, time-dependent birth-death process conditioned on the age of the tree. The rates may be any positive function of time or a constant. Additionally, mass-extinction event can be provided and a random taxon sampling probability. It is possible to start either with the origin (1 species) or with the most recent common ancestor (2 species).}
\usage{
sim.globalBiDe.age(lambda, mu, massExtinctionTimes = c(), massExtinctionSurvivalProbabilities = c(), samplingProbability = 1, age, MRCA = TRUE, approx = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{The speciation rate function or constant.}
  \item{mu}{The extinction rate function or constant.}
  \item{massExtinctionTimes}{The set of mass-extinction times after the start of the process.}
  \item{massExtinctionSurvivalProbabilities}{The set of survival probabilities for each speciation event. The set must have the same length as the set of mass-extinction times.}
  \item{samplingProbability}{The probability for a species to be included in the sample.}
  \item{age}{The age of the tree, i.e. the time to simulate.}
  \item{MRCA}{Does the process start with the most recent common ancestor?}
  \item{approx}{Use approximatian (faster) of the integrals?}
}
\value{Returns a tree in 'phylo' format.
}
\references{S. Hoehna, Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2012, submitted }
\author{Sebastian Hoehna}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

l <- function(x) { if (x > 0.5 || x < 0.3) { return (1) } else { return (2) } }
e <- function(x) { if (x > 0.5 || x < 0.3) { return (0.95) } else { return (0.5) } }

sim.globalBiDe.age(l,e,age=1,MRCA=TRUE,approx=TRUE)

# simulation under global rates
sim.globalBiDe.age(2.0,1.0,age=1,MRCA=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ likelihood }
%\keyword{ tree }
%\keyword{ birth-death }% __ONLY ONE__ keyword per line

\name{sim.globalBiDe.taxa}
\alias{sim.globalBiDe.taxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sim.globalBiDe.taxa.taxa: Simulate a reconstructed tree for a given number of taxa under a global, time-dependent birth-death process.}
\description{sim.globalBiDe.taxa simlates a reconstructed phylogenetic tree under a global, time-dependent birth-death process conditioned on the number of taxa sampled. The rates may be any positive function of time or a constant. Additionally, mass-extinction event can be provided and a random taxon sampling probability. It is possible to start either with the origin (1 species) or with the most recent common ancestor (2 species).}
\usage{
sim.globalBiDe.taxa(lambda, mu, massExtinctionTimes = c(), massExtinctionSurvivalProbabilities = c(), samplingProbability = 1, nTaxa, max, MRCA = TRUE, approx = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{The speciation rate function or constant.}
  \item{mu}{The extinction rate function or constant.}
  \item{massExtinctionTimes}{The set of mass-extinction times after the start of the process.}
  \item{massExtinctionSurvivalProbabilities}{The set of survival probabilities for each speciation event. The set must have the same length as the set of mass-extinction times.}
  \item{samplingProbability}{The probability for a species to be included in the sample.}
  \item{nTaxa}{Number of species sampled.}
  \item{max}{Maximum time of the tree.}
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

sim.globalBiDe.taxa(l,e,nTaxa=10,max=10,MRCA=TRUE,approx=TRUE)

# simulation under constant rates
sim.globalBiDe.taxa(2.0,1.0,nTaxa=10,max=10,MRCA=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ likelihood }
%\keyword{ tree }
%\keyword{ birth-death }% __ONLY ONE__ keyword per line

\name{tess.likelihood}
\alias{tess.likelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tess.likelihood: Probability density of a tree under a tree-wide time-dependent birth-death process}
\description{tess.likelihood computes the probability of a reconstructed phylogenetic tree under time-dependent diversification rates. The rates may be any positive function of time or a constant. Additionally, mass-extinction event can be provided and a uniform taxon sampling probability. You have several options for the start of the process (origin vs MRCA) and the condition of the process (time, survival or taxa; note that survival and taxa implicitly condition on the time too!). See equation (5) in Hoehna (2013) for more information. Note that constant rates lead to much faster computations.
The likelihood can be computed for incompletely sampled trees if you give a sampling probability != 1.0. You have two options for the sampling strategy: uniform|diversified. The detailed description of these can be found in the references. }
\usage{
tess.likelihood(times,
                lambda,
                mu,
                massExtinctionTimes=c(),
                massExtinctionSurvivalProbabilities=c(),
                missingSpecies = c(),
                timesMissingSpecies = c(),
                samplingProbability=1.0,
                samplingStrategy="uniform",
                MRCA=TRUE,
                CONDITION="survival",
                log=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{The branching times of the phylogeny.}
  \item{lambda}{The speciation rate function or constant.}
  \item{mu}{The extinction rate function or constant.}
  \item{massExtinctionTimes}{The set of mass-extinction times after the start of the process.}
  \item{massExtinctionSurvivalProbabilities}{The set of survival probabilities for each speciation event. The set must have the same length as the set of mass-extinction times.}
  \item{missingSpecies}{The number of species missed which originated in a given time interval (empirical taxon sampling).}
  \item{timesMissingSpecies}{The times intervals of the missing species (empirical taxon sampling).}
  \item{samplingProbability}{The probability for a species to be included in the sample.}
  \item{samplingStrategy}{The strategy how samples were obtained. Options are: uniform|diversified|age.}
  \item{MRCA}{Does the process start with the most recent common ancestor? If not, the tree must have a root edge!}
  \item{CONDITION}{do we condition the process on time|survival|taxa?}
  \item{log}{Should we log-transform the likelihood?}
}
\value{Returns the (log) probability of the tree, i.e. the likelihood of the parameters given the tree.}
\references{S. Hoehna: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2013, Bioinformatics, 29:1367-1374
Hoehna, S. Likelihood Inference of Non-Constant Diversification Rates with Incomplete Taxon Sampling. 2014, PLoS one, Public Library of Science, 9, e84184}
\author{Sebastian Hoehna}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(cettiidae)
times <- as.numeric( branching.times(cettiidae) )

l <- function(x) { if (x > 0.5 || x < 0.3) { return (1) } else { return (2) } }
e <- function(x) { if (x > 0.5 || x < 0.3) { return (0.95) } else { return (0.5) } }

tess.likelihood(times,l,e,MRCA=TRUE,log=TRUE)

# constant rate likelihood function
tess.likelihood(times,2.0,1.0,MRCA=TRUE,log=TRUE)

# constant rate pure birth likelihood function
tess.likelihood(times,2.0,0.0,MRCA=TRUE,log=TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }
%\keyword{ birth-death }% __ONLY ONE__ keyword per line

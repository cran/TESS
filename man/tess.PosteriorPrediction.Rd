\name{tess.PosteriorPrediction}
\alias{tess.PosteriorPrediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tess.PosteriorPrediction: Approximation of the posterior predictive distribution.}
\description{tess.PosteriorPrediction calls the simulation function exactly once for each sampled parameter combination. In that way, posterior predictive simulations can be obtained.}
\usage{
tess.PosteriorPrediction(simulationFunction,parameters)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simulationFunction}{The simulation function which will be called internally by simulationFunction(parameters).}
  \item{parameters}{A matrix of parameters where the rows represent samples of parameters and the column the different parameters.}
}
\value{Returns samples simulated from the posterior predictive distribution.}
\references{Hoehna, S.: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2012, submitted }
\author{Sebastian Hoehna}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# We first run an MCMC to obtain samples from the posterior distribution and then simulate the posterior predictive distribution.

# The bird phylogeny as the test data set
data(cettiidae)

# The log-likelihood function
likelihood <- function(params) {
  # We use the parameters as diversification rate and turnover rate.
  # Thus we need to transform first
  b <- params[1] + params[2]
  d <- params[2]
  
  lnl <- globalBiDe.likelihood(cettiidae,b,d,samplingProbability=1.0,log=TRUE)
  return (lnl)
}

prior_diversification <- function(x) { dexp(x,rate=0.1,log=TRUE) }
prior_turnover <- function(x) { dexp(x,rate=0.1,log=TRUE) }
priors <- c(prior_diversification,prior_turnover)

# Note, the number of iterations and the burnin is too small here and should be adapted for real analyses
samples <- tess.mcmc(likelihood,priors,c(1,0.1),c(TRUE,TRUE),c(0.1,0.1),10,10)

tmrca <- max(branching.times(cettiidae))
# The simulation function
sim <- function(params) {
  # We use the parameters as diversification rate and turnover rate.
  # Thus we need to transform first
  b <- params[1] + params[2]
  d <- params[2]
  
  tree <- sim.globalBiDe.age(n=1,age=tmrca,b,d,samplingProbability=1.0)[[1]]
  return (tree)
}

trees <- tess.PosteriorPrediction(sim,samples)
#plot(trees)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
%\keyword{ tree }
%\keyword{ birth-death }% __ONLY ONE__ keyword per line

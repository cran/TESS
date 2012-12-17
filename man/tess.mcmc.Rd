\name{tess.mcmc}
\alias{tess.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tess.mcmc: Metropolis-Hastings algorithm.}
\description{tess.mcmc constructs a Markov chain Monte Carlo algorithm (MCMC) by implementing a general Metropolis-Hastings algorithm.Any model can be used where the likelihood is known and thus can be passed in as an argument. The parameters have to be continuous. Proposals are taken from a normal distribution centered around the current value. The algorithm creates sampels from the posterior probility distribution and returns these a CODA mcmc object.}
\usage{
tess.mcmc(likelihoodFunction,priors,parameters,logTransforms,iterations,
             burnin=round(iterations/3),thining=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{likelihoodFunction}{The log-likelihood function which will be called internally by likelihoodFunction(parameters).}
  \item{priors}{A list of functions of the log-prior-densities of each parameter.}
  \item{parameters}{The initial parameter value list.}
  \item{logTransforms}{A vector of booleans telling if log-transform for the parameters should be used (e.g. for rates).}
  \item{iterations}{The number of iterations for the MCMC.}
  \item{burnin}{The number of iterations to burn before starting the MCMC.}
  \item{thining}{The frequency of taking a sample of the parameters.}
}
\value{Returns the posterior samples for the parameters.}
\references{Hoehna, S.: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2012, submitted }
\author{Sebastian Hoehna}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

data(cettiidae)

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
samples <- tess.mcmc(likelihood,priors,runif(2,0,1),c(TRUE,TRUE),10,10)
#plot(samples)
summary(samples)
colMeans(samples)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
%\keyword{ tree }
%\keyword{ birth-death }% __ONLY ONE__ keyword per line

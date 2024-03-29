\name{tess.process.output}
\alias{tess.process.output}
\alias{globalBiDe.output.summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tess.process.output: Summarizing the output of a diversification rate estimation including mass-extinction events. See the tess.analysis function for more information on how such output is generated and the tess.plot.output how the output can be visualized. Also have a look at the vignette for more in detail description and examples.}
\description{tess.process.output summarizes the output generated by a tess.analysis(...) run.}
\usage{
tess.process.output(dir,
                    tree=NULL,
                    numExpectedRateChanges=2,
                    numExpectedMassExtinctions=2,
                    burnin=0.25,
                    numIntervals=100,
                    criticalBayesFactors=c(2,6,10))

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{The directory from which the CoMET output will be read.}
  \item{tree}{The tree analyzed with CoMET in phylo format. By default, looks for a tree in the target directory.}
  \item{numExpectedRateChanges}{The number of expected diversification-rate changes.}
  \item{numExpectedMassExtinctions}{The number of expected mass-extinction events.}
  \item{burnin}{The fraction of samples that will be discarded as burnin.}
  \item{numIntervals}{The number of discrete intervals in which to break the tree.}
  \item{criticalBayesFactors}{The Bayes factor thresholds to use to assess significance of events.}
}

\details{
The output of a CoMET analysis is stored in a directory with different files containing the MCMC samples from the posterior distribution.
For example, the tess.analysis function stores the times and survival probabilities of the mass-extinction events in a file.
This function, converts the output by counting the number of events that fall into a given time-bin.
This pre-processing of the output simplifies the plotting.
}

\value{
This function returns a list with the following elements:

  \item{posterior}{An object of class 'mcmc' that contains the trace of the model's posterior probability.}
  \item{numSpeciationCategories}{An object of class 'mcmc' that contains samples from the posterior distribution of the number of speciation categories (minimum 1, since this includes the initial speciation rate).}
  \item{numExtinctionCategories}{An object of class 'mcmc' that contains samples from the posterior distribution of the number of extinction categories (minimum 1, since this includes the initial extinction rate).}
  \item{numMassExtinctions}{An object of class 'mcmc' that contains samples from the posterior distribution of the number of mass-extinction events.}

  \item{speciation rates}{An object of class 'mcmc' that contains speciation rates sampled from the posterior distribution for each of numIntervals discrete time intervals. Rows correspond to samples from the posterior distribution, columns correspond to intervals.}
  \item{speciation change times}{An object of class 'mcmc' that contains speciation-rate-change events sampled from the posterior distribution for each of numIntervals discrete time intervals. A value of 1 indicates an event was contained in the interval, 0 that no event was contained in the interval. Rows correspond to samples from the posterior distribution, columns correspond to intervals.}
  \item{speciation Bayes factors}{A vector of class 'numeric' that contains the Bayes factor support for there being a speciation-rate-change event for each of numIntervals discrete time intervals. The ith element corresponds to the Bayes factor support for an event in the ith interval.}
  \item{speciationRateChangeCriticalPosteriorProbabilities}{A vector of posterior probabilities that correspond to critical Bayes factor thresholds (specified by the argument criticalBayesFactors). Element i is the posterior probability of a speciation-rate-change event in an interval needed to produce Bayes factor support of criticalBayesFactors[i].}

  \item{extinction rates}{An object of class 'mcmc' that contains extinction rates sampled from the posterior distribution for each of numIntervals discrete time intervals. Rows correspond to samples from the posterior distribution, columns correspond to intervals.}
  \item{extinction change times}{An object of class 'mcmc' that contains extinction-rate-change events sampled from the posterior distribution for each of numIntervals discrete time intervals. A value of 1 indicates an event was contained in the interval, 0 that no event was contained in the interval. Rows correspond to samples from the posterior distribution, columns correspond to intervals.}
  \item{extinction Bayes factors}{A vector of class 'numeric' that contains the Bayes factor support for there being a extinction-rate-change event for each of numIntervals discrete time intervals. The ith element corresponds to the Bayes factor support for an event in the ith interval.}
  \item{extinctionRateChangeCriticalPosteriorProbabilities}{A vector of posterior probabilities that correspond to critical Bayes factor thresholds (specified by the argument criticalBayesFactors). Element i is the posterior probability of a extinction-rate-change event in an interval needed to produce Bayes factor support of criticalBayesFactors[i].}

  \item{net-diversification rates}{An object of class 'mcmc' that contains net-diversification (speciation - extinction) rates sampled from the posterior distribution for each of numIntervals discrete time intervals. Rows correspond to samples from the posterior distribution, columns correspond to intervals.}
  \item{relative-extinction rates}{An object of class 'mcmc' that contains relative-extinction (extinction / speciation) rates sampled from the posterior distribution for each of numIntervals discrete time intervals. Rows correspond to samples from the posterior distribution, columns correspond to intervals.}


  \item{mass extinction times}{An object of class 'mcmc' that contains mass-extinction events sampled from the posterior distribution for each of numIntervals discrete time intervals. A value of 1 indicates an event was contained in the interval, 0 that no event was contained in the interval. Rows correspond to samples from the posterior distribution, columns correspond to intervals.}
  \item{mass extinction Bayes factors}{A vector of class 'numeric' that contains the Bayes factor support for there being a mass-extinction event for each of numIntervals discrete time intervals. The ith element corresponds to the Bayes factor support for an event in the ith interval.}
  \item{massExtinctionCriticalPosteriorProbabilities}{A vector of posterior probabilities that correspond to critical Bayes factor thresholds (specified by the argument criticalBayesFactors). Element i is the posterior probability of a mass-extinction event in an interval needed to produce Bayes factor support of criticalBayesFactors[i].}

  \item{criticalBayesFactors}{The critical Bayes factor values used for the Bayes factor tests (default 2 ln BF = \{2,6,10\}).}
  \item{tree}{The tree analyzed with CoMET (just in case).}
  \item{intervals}{The discrete intervals used to compute the interval-specific parameters.}

}

\author{Michael R. May}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# Load the data, compute the sampling fraction rho
data(conifers)
totalConiferSpecies <- 630
sampledConiferSpecies <- conifers$Nnode+1
rho <- sampledConiferSpecies / totalConiferSpecies

# Run a tess analysis
tess.analysis(tree = conifers,
              initialSpeciationRate=c(1.0),
              initialExtinctionRate=c(0.5),
              empiricalHyperPriors = FALSE,
              numExpectedRateChanges = 2,
              numExpectedMassExtinctions = 2,
              samplingProbability = rho,
              MAX_ITERATIONS = 200,
              BURNIN=100)

# Process the output
coniferOutput <- tess.process.output(dir=getwd(),
                                     numExpectedRateChanges=2,
                                     numExpectedMassExtinctions=2)

# Plot the output
tess.plot.output(coniferOutput)

\dontshow{
	unlink("*.txt")
	unlink("*.pdf")
	unlink("*.tre")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }
%\keyword{ birth-death }% __ONLY ONE__ keyword per line

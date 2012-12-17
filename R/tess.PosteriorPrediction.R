require(coda)


################################################################################
# 
# @brief General model adequacy test using posterior predictive testing. Prior
#        predictive testing can be achieved by providing samples from the prior
#        instead.
#
# @date Last modified: 2012-12-17
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-11-19, version 1.0
#
# @param    simulationFunction     function      the simulation function
# @param    parameters             matrix        set of parameter samples
# @return                          scalar        the upper quantile of observing such a value
#
################################################################################


tess.PosteriorPrediction <- function(simulationFunction,parameters) {


  samples <- list()
  for ( i in 1:length(parameters[,1])) {

    # get the current set of parameter values
    theta <- parameters[i,]

    # simulate a new observation under the current parameter values
    samples[[i]] <- simulationFunction(theta)

  }

  return (samples)
}

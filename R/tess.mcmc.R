require(coda)


################################################################################
# 
# @brief General Markov chain Monte Carlo algorithm using Metropolis-Hastings.
#
# @date Last modified: 2012-12-17
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-11-06, version 1.1
#
# @param    likelihoodFunction     function      the log-likelihood function
# @param    priors                 list          a list of functions of the log-prior density per parameter
# @param    parameters             vector        initial set of paramete values
# @param    logTransform           vector        should be a log-transform be used for parameter proposals
# @param    iterations             scalar        the number of iterations
# @param    burnin                 scalar        number of burnin iterations
# @param    thining                scalar        number of iterations between samples
# @return                          list          samples from the posterior
#
################################################################################


tess.mcmc <- function(likelihoodFunction,priors,parameters,logTransforms,iterations,burnin=round(iterations/3),thining=1) {

  # create a list for the samples
  chain = array(dim = c(floor(iterations/thining)+1,length(parameters))) #reserve memory for the chain, for large chains we might consider writing to a file instead of storing in memory
   

  # pre-compute current posterior probability
  pp <- likelihoodFunction(parameters)
  for ( j in 1:length(parameters) ) {
    pp <- pp + priors[[j]](parameters[j])
  }

  for (i in 1:(burnin+iterations)) {

    # propose new values
    for ( j in 1:length(parameters) ) {
      if ( logTransforms[j] == TRUE ) {
        if (parameters[j] == 0) {
          stop("Cannot propose new value for a parameter with value 0.0.")
        }
        eta           <- log(parameters[j]) ### propose a new value for parameter[j]
        new_eta       <- eta + rnorm(1,0,1)
        new_val       <- exp(new_eta)
        hr            <- log(new_val / parameters[j]) # calculate the Hastings ratio
        parameters[j] <- new_val
        new_pp        <- likelihoodFunction(parameters)
        for ( k in 1:length(parameters) ) {
          new_pp <- new_pp + priors[[k]](parameters[k])
        }
        # accept / reject
        if ( is.finite(new_pp) && is.finite(hr) &&  new_pp-pp+hr > log(runif(1,0,1)) ) {
          pp <- new_pp
        } else {
          parameters[j] <- exp(eta)
        }
      } else {
        eta           <- parameters[j] ### propose a new value for parameter[j]
        new_val       <- eta + rnorm(1,0,1)
        hr            <- 0.0 # calculate the Hastings ratio
        parameters[j] <- new_val
        new_pp        <- likelihoodFunction(parameters)
        for ( k in 1:length(parameters) ) {
          new_pp <- new_pp + priors[[k]](parameters[k])
        }
        # accept / reject
        if ( is.finite(new_pp) && is.finite(hr) &&  new_pp-pp+hr > log(runif(1,0,1)) ) {
          pp <- new_pp
        } else {
          parameters[j] <- eta
        }
      }

    }

    # sample the parameter
    if (i >= burnin) {
      if ( (i-burnin) %% thining == 0 ) {
        chain[(i-burnin)/thining+1,] <- parameters
      }
    }
    
  }

  return(as.mcmc(chain)) #return a mcmc object, used by coda to plot
}

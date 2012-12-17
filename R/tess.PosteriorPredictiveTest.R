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
# @since 2012-12-17, version 1.0
#
# @param    samples                list          Samples from the posterior predictive distribution
# @param    observation            any           the observed value
# @param    statistic              function      the function computing the statistic
# @return                          scalar        the upper quantile of observing such a value
#
################################################################################


tess.PosteriorPredictiveTest <- function(samples,observation,statistic) {

  obs <- statistic(observation)

  sampled_statistics <- c()
  count <- 0
  for ( i in 1:length(samples)) {

    # compute the statistic for the i-th sample
    sampled_statistics[i] <- statistic(samples[[i]])

    if ( sampled_statistics[i] < obs ) {
      count <- count + 1
    }
  }

  p <- count / length(samples)

  return (list(samples=sampled_statistics,pvalue=p))
}

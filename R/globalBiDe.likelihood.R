################################################################################
#
# globalBiDe.likelihood.R
#
# Copyright (c) 2012- Sebastian Hoehna
#
# This file is part of TESS.
# See the NOTICE file distributed with this work for additional
# information regarding copyright ownership and licensing.
#
# TESS is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
#  TESS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with TESS; if not, write to the
# Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
# Boston, MA  02110-1301  USA
#
################################################################################



################################################################################
# 
# @brief Computation of the likelihood for a given tree.
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.2
# @since 2012-09-22, version 1.0
#
# @param    tree                                          phylo         the tree
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    samplingStrategy                              string        Which strategy was used to obtain the samples (taxa). Options are: random|diversified|age
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    CONDITITON                                    string        do we condition the process on nothing|survival|taxa?
# @param    log                                           boolean       likelhood in log-scale?

# @return                                                 scalar        probability of the speciation times
#
################################################################################

globalBiDe.likelihood <- function(tree,lambda,mu,massExtinctionTimes=c(),massExtinctionSurvivalProbabilities=c(),t.crit=c(),samplingProbability=1.0,samplingStrategy="random",MRCA=TRUE,CONDITION="survival",log=TRUE) {

  if ( length(massExtinctionTimes) != length(massExtinctionSurvivalProbabilities) ) {
    stop("Number of mass-extinction times needs to equals the number of mass-extinction survival probabilities!")
  }

  if ( CONDITION != "time" && CONDITION != "survival" && CONDITION != "taxa" ) {
    stop("Wrong choice of argument for \"CONDITION\". Possible option are time|survival|taxa.")
  }

  if ( samplingStrategy != "random" && samplingStrategy != "diversified" && samplingStrategy != "age") {
    stop("Wrong choice of argument for \"samplingStrategy\". Possible option are random|diversified.")
  } 

  if ( (!is.numeric(lambda) && !inherits(lambda, "function")) || (!is.numeric(mu) && !inherits(mu, "function"))) {
    stop("Unexpected parameter types for lambda and mu!")
  }
  
  # test if we got constant values for the speciation and extinction rates
  if ( is.numeric(lambda) && is.numeric(mu) && samplingStrategy != "age" ) {
    # call computations for constant rates (much faster)
    p <- globalBiDe.likelihood.constant(tree,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA,CONDITION,log)
    return (p)
  } else  {
    
    # convert the speciation rate into a function if necessary
    if ( is.numeric(lambda) ) {
      speciation <- function (x) rep(lambda,length(x))
    } else {
      speciation <- lambda
    }
    # convert the extinction rate into a function if necessary
    if ( is.numeric(mu) ) {
      extinction <- function (x) rep(mu,length(x))
    } else {
      extinction <- mu
    }

    if ( samplingStrategy == "age" ) {
      p <- globalBiDe.likelihood.age(tree,speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,t.crit,samplingProbability,MRCA,CONDITION,log)
    } else {
      p <- globalBiDe.likelihood.function(tree,speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,t.crit,samplingProbability,samplingStrategy,MRCA,CONDITION,log)
    }
    
    return (p)
  }

}



################################################################################
# 
# @brief Computation of the likelihood for a given tree.
#
# Here we use equation (6) from Hoehna, S., 2013, Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes
# For the diversified taxon sampling, see Hoehna, S., 2013, A Birth-Death Process with Decreasing Diversification Rate and Diversified Taxon Sampling
#
# @date Last modified: 2012-09-22
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-22, version 1.0
#
# @param    tree                                          phylo         the tree
# @param    lambda                                        scalar        speciation rate function
# @param    mu                                            scalar        extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    samplingStrategy                              string        Which strategy was used to obtain the samples (taxa). Options are: random|diversified
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    CONDITITON                                    string        do we condition the process on nothing|survival|taxa?
# @return                                                 scalar        probability of the speciation times
#
################################################################################

globalBiDe.likelihood.constant <- function(tree,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA,CONDITION,log) {

  times <- as.numeric(branching.times(tree))
  PRESENT <- max(times)
  nTaxa <- length(times) + 1
  
  # if we have a root edge, we need to add this time
  if ( MRCA == FALSE ) {
    if ( !is.null(tree$root) ) {
      PRESENT <- PRESENT + tree$root
    } else if ( !is.null(tree$rootEdge) ) {
      PRESENT <- PRESENT + tree$rootEdge
    } else {
      MRCA <- TRUE
    }
  }
  times <- PRESENT - sort(times,decreasing=TRUE)

  # if we condition on the MRCA, then we need to remove the root speciation event
  if ( MRCA == TRUE ) {
    times <- times[-1]
  }

  # set the random taxon sampling probability
  if (samplingStrategy == "random") {
    rho <- samplingProbability
  } else {
    rho <- 1.0
  }
  
  # initialize the log likelihood
  lnl <- 0

  # what do we condition on?
  # did we condition on survival?
  if ( CONDITION == "survival" || CONDITION == "taxa" )    lnl <- - globalBiDe.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,rho,0,PRESENT,PRESENT,log=TRUE)
  
  # multiply the probability of a descendant of the initial species
  lnl <- lnl + globalBiDe.equations.p1.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,rho,0,PRESENT,log=TRUE)

  # add the survival of a second species if we condition on the MRCA
  if ( MRCA == TRUE ) {
    lnl <- 2*lnl
  } 

  # did we condition on observing n species today
  if ( CONDITION == "taxa" )    lnl <- lnl - globalBiDe.equations.pN.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,rho,nTaxa,0,PRESENT,MRCA,log=TRUE)

  # if we assume diversified sampling, we need to multiply with the probability that all missing species happened after the last speciation event
  if ( samplingStrategy == "diversified" ) {
    # We use equation (5) of Hoehna et al. "Inferring Speciation and Extinction Rates under Different Sampling Schemes"
    lastEvent <- times[length(times)]
    p_0_T <- 1.0 - globalBiDe.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,1.0,0,PRESENT,PRESENT,log=FALSE) * exp((mu-lambda)*PRESENT)
    p_0_t <- 1.0 - globalBiDe.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,1.0,lastEvent,PRESENT,PRESENT,log=FALSE)*exp((mu-lambda)*(PRESENT-lastEvent))
    F_t <- p_0_t / p_0_T
    # get an estimate of the actual number of taxa
    m <- round(nTaxa / samplingProbability)
    lnl <- lnl + (m-nTaxa) * log(F_t) + log(choose(m,nTaxa))
  }

  # multiply the probability for each speciation time
  if ( length(times) > 0 ) {
    lnl <- lnl + length(times)*log(lambda) + sum(globalBiDe.equations.p1.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,rho,times,PRESENT,log=TRUE))
  }

  if ( log == FALSE ) {
    lnl <- exp(lnl)
  }
  
  if (is.nan(lnl)) lnl <- -Inf
  
  return (lnl)
}



################################################################################
# 
# @brief Computation of the likelihood for a given tree. 
#
# Here we use equation (6) from Hoehna, S., 2013, Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes
# For the diversified taxon sampling, see Hoehna, S., 2013, A Birth-Death Process with Decreasing Diversification Rate and Diversified Taxon Sampling
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-22, version 1.0
#
# @param    tree                                          phylo         the tree
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    samplingStrategy                              string        Which strategy was used to obtain the samples (taxa). Options are: random|diversified
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    CONDITITON                                    string        do we condition the process on nothing|survival|taxa?
# @return                                                 scalar        probability of the speciation times
#
################################################################################

globalBiDe.likelihood.function <- function(tree,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,t.crit,samplingProbability,samplingStrategy,MRCA,CONDITION,log) {

  lnl <- -Inf
   
  times <- as.numeric(branching.times(tree))
  PRESENT <- max(times)
  nTaxa <- length(times) + 1
  
  # if we have a root edge, we need to add this time
  if (MRCA == FALSE) {
    if ( !is.null(tree$root) ) {
      PRESENT <- PRESENT + tree$root
    } else if ( !is.null(tree$rootEdge) ) {
      PRESENT <- PRESENT + tree$rootEdge
    } 
  }
  times <- PRESENT - sort(times,decreasing=TRUE)

  # if we condition on the MRCA, then we need to remove the root speciation event
  if (MRCA == TRUE) {
    times <- times[-1]
  }

  # set the random taxon sampling probability
  if (samplingStrategy == "random") {
    rho <- samplingProbability
  } else {
    rho <- 1.0
  }

  # prepare the integrals
  tryCatch({
    approxFuncs <- tess.prepare.pdf(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,PRESENT,c(t.crit,times))
    
  # initialize the log likelihood
  lnl <- 0

  # what do we condition on?
  # did we condition on survival?
  if ( CONDITION == "survival" || CONDITION == "taxa" )    lnl <- - globalBiDe.equations.pSurvival.fastApprox(approxFuncs$r,approxFuncs$s,rho,0,PRESENT,PRESENT,log=TRUE)


  # multiply the probability of a descendant of the initial species
  lnl <- lnl + globalBiDe.equations.p1.fastApprox(approxFuncs$r,approxFuncs$s,rho,0,PRESENT,log=TRUE)

  # add the survival of a second species if we condition on the MRCA
  if ( MRCA == TRUE ) {
    lnl <- 2*lnl 
  } 

  # did we condition on observing n species today
  if ( CONDITION == "taxa" )    lnl <- lnl - globalBiDe.equations.pN.fastApprox(approxFuncs$r,approxFuncs$s,rho,nTaxa,0,PRESENT,MRCA,log=TRUE)

  # if we assume diversified sampling, we need to multiply with the probability that all missing species happened after the last speciation event
  if ( samplingStrategy == "diversified" ) {
    # We use equation (5) of Hoehna et al. "Inferring Speciation and Extinction Rates under Different Sampling Schemes"
    lastEvent <- times[length(times)]
    
    p_0_T <- 1.0 - globalBiDe.equations.pSurvival.fastApprox(approxFuncs$r,approxFuncs$s,1.0,0,PRESENT,PRESENT,log=FALSE) * exp(approxFuncs$r(PRESENT))
    p_0_t <- (1.0 - globalBiDe.equations.pSurvival.fastApprox(approxFuncs$r,approxFuncs$s,1.0,lastEvent,PRESENT,PRESENT,log=FALSE) * exp(approxFuncs$r(PRESENT) - approxFuncs$r(lastEvent)))
    F_t <- p_0_t / p_0_T

    # get an estimate of the actual number of taxa
    m <- round(nTaxa / samplingProbability)
    lnl <- lnl + (m-nTaxa) * log(F_t) + log(choose(m,nTaxa))
  }

  # multiply the probability for each speciation time
  if ( length(times) > 0 ) {
    lnl <- lnl + sum(log(lambda(times)) + globalBiDe.equations.p1.fastApprox(approxFuncs$r,approxFuncs$s,rho,times,PRESENT,log=TRUE))
  }

  }, warning = function(w) {  }, error = function(e) { })
  
  if ( log == FALSE ) {
    lnl <- exp(lnl)
  }
  
  if (is.nan(lnl)) lnl <- -Inf
  
  return (lnl)
    

}



################################################################################
# 
# @brief Computation of the likelihood for a given tree. 
#
# Here we use equation (xxx) from Hoehna, S. and Catalan, A., 2013, Estimating Diversification Rates and Patterns in Mammals.
#
# @date Last modified: 2013-03-06
# @author Sebastian Hoehna
# @version 1.2
# @since 2013-02-06, version 1.2
#
# @param    tree                                          phylo         the tree
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    samplingStrategy                              string        Which strategy was used to obtain the samples (taxa). Options are: random|diversified
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    CONDITITON                                    string        do we condition the process on nothing|survival|taxa?
# @return                                                 scalar        probability of the speciation times
#
################################################################################

globalBiDe.likelihood.age <- function(tree,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,t.crit,samplingProbability,MRCA,CONDITION,log) {

  lnl <- -Inf
  
  times <- as.numeric(branching.times(tree))
  PRESENT <- max(times)
  nTaxa <- length(times) + 1
  
  if ( inherits(samplingProbability, "numeric") ) {
    rho <- function(t) rep(samplingProbability,length(t))
  } else {
    rho <- samplingProbability
  }

  
  # if we have a root edge, we need to add this time
  if (MRCA == FALSE) {
    if ( !is.null(tree$root) ) {
      PRESENT <- PRESENT + tree$root
    } else if ( !is.null(tree$rootEdge) ) {
      PRESENT <- PRESENT + tree$rootEdge
    } 
  }
  times <- PRESENT - sort(times,decreasing=TRUE)

  # if we condition on the MRCA, then we need to remove the root speciation event
  if (MRCA == TRUE) {
    nInitial <- 2
  } else {
    nInitial <- 1
  }

  # prepare the integrals
  tryCatch({
    approxFuncs <- tess.prepare.pdf(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,PRESENT,t.crit)
  
  # initialize the log likelihood
  lnl <- 0

  # what do we condition on?
  # did we condition on survival?
  if ( CONDITION != "survival" )    lnl <- globalBiDe.equations.pSurvival.fastApprox(approxFuncs$r,approxFuncs$s,rho(0),0,PRESENT,PRESENT,log=TRUE)

  # add the survival of a second species if we condition on the MRCA
  if ( MRCA == TRUE ) {
    lnl <- 2*lnl
  }

  # for every speciation event
  for (i in seq_len(length(times)-1)) {
    # Multiply with the probability of no speciation between t_i and t_(i+1) if i species are alive in this interval and
    # multiply the probability density of observing a speciation event exactly at time t_(i+1)
    lnl <- lnl + globalBiDe.equations.pWaiting.ageSampling(lambda,approxFuncs$r,approxFuncs$s,rho,times[i],times[i+1],PRESENT,i+nInitial-1,log=TRUE) + globalBiDe.equations.pSurvival.fastApprox(approxFuncs$r,approxFuncs$s,rho(times[i+1]),times[i+1],PRESENT,PRESENT,log=TRUE) + log(lambda(times[i+1]))
  }

  # Multiply with the probability that we didn't observe a speciation event between the last speciation event t_(N-1) and the present t_N
  lnl <- lnl + globalBiDe.equations.pWaiting.ageSampling(lambda,approxFuncs$r,approxFuncs$s,rho,times[length(times)],PRESENT,PRESENT,nTaxa,log=TRUE)
  

  }, warning = function(w) {  }, error = function(e) {  })

  if ( log == FALSE ) {
    lnl <- exp(lnl)
  }
  
  if (is.nan(lnl)) lnl <- -Inf
  
  return (lnl)
}





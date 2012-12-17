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
# @version 1.1
# @since 2012-09-22, version 1.0
#
# @param    tree                                          phylo         the tree
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    CONDITITON                                    string        do we condition the process on nothing|survival|taxa?
# @param    log                                           boolean       likelhood in log-scale?

# @return                                                 scalar        probability of the speciation times
#
################################################################################

globalBiDe.likelihood <- function(tree,lambda,mu,massExtinctionTimes=c(),massExtinctionSurvivalProbabilities=c(),samplingProbability=1.0,MRCA=TRUE,CONDITION="survival",log=TRUE) {

  if ( length(massExtinctionTimes) != length(massExtinctionSurvivalProbabilities) ) {
    stop("Number of mass-extinction times needs to equals the number of mass-extinction survival probabilities!")
  }

  if ( CONDITION != "time" && CONDITION != "survival" && CONDITION != "taxa" ) {
    stop("Wrong choice of argument for \"CONDITION\". Possible option are time|survival|taxa.")
  }
  
  # test if we got constant values for the speciation and extinction rates
  if ( class(lambda) == "numeric" && class(mu) == "numeric" ) {
    # call computations for constant rates (much faster)
    p <- globalBiDe.likelihood.constant(tree,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA,CONDITION,log)
    return (p)
  } else if ( class(lambda) == "function" && class(mu) == "numeric" ) {
    extinction <- function (x) mu
    p <- globalBiDe.likelihood.function(tree,lambda,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA,CONDITION,log)
    return (p)
  } else if ( class(lambda) == "numeric" && class(mu) == "function" ) {
    speciation <- function (x) lambda
    p <- globalBiDe.likelihood.function(tree,speciation,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA,CONDITION,log)
    return (p)
  } else if ( class(lambda) == "function" && class(mu) == "function" ) {
    p <- globalBiDe.likelihood.function(tree,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA,CONDITION,log)
    return (p)
  } else {
    stop("Unexpected parameter types for lambda and mu!")
  }

}



################################################################################
# 
# @brief Computation of the likelihood for a given tree.
#
# Here we use equation (6) from Hoehna, S., 2013, Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes
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
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    CONDITITON                                    string        do we condition the process on nothing|survival|taxa?
# @return                                                 scalar        probability of the speciation times
#
################################################################################

globalBiDe.likelihood.constant <- function(tree,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA,CONDITION,log) {

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

  # initialize the log likelihood
  lnl <- 0

  # what do we condition on?
  # did we condition on survival?
  if ( CONDITION == "survival" || CONDITION == "taxa" )    lnl <- - globalBiDe.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,0,PRESENT,PRESENT,log=TRUE)
  
  # multiply the probability of a descendant of the initial species
  lnl <- lnl + globalBiDe.equations.p1.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,0,PRESENT,log=TRUE)

  # add the survival of a second species if we condition on the MRCA
  if ( MRCA == TRUE ) {
    lnl <- 2*lnl 
  } 

  # did we condition on observing n species today
  if ( CONDITION == "taxa" )    lnl <- lnl - globalBiDe.equations.pN.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,0,PRESENT,MRCA,log=TRUE)


  # multiply the probability for each speciation time
  for ( i in seq_len(length(times)) ) {
    lnl <- lnl + log(lambda) + globalBiDe.equations.p1.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,times[i],PRESENT,log=TRUE)
  }


  if ( log == FALSE ) {
    lnl <- exp(lnl)
  }
  
  return (lnl)
}



################################################################################
# 
# @brief Computation of the likelihood for a given tree. 
#
# Here we use equation (6) from Hoehna, S., 2013, Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes
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
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    CONDITITON                                    string        do we condition the process on nothing|survival|taxa?
# @return                                                 scalar        probability of the speciation times
#
################################################################################

globalBiDe.likelihood.function <- function(tree,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA,CONDITION,log) {
  
  times <- as.numeric(branching.times(tree))
  PRESENT <- max(times)
  nTaxa <- length(times) + 1
  
  approxFuncs <- tess.prepare.pdf(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,PRESENT,c())
  
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

  # initialize the log likelihood
  lnl <- 0

  # what do we condition on?
  # did we condition on survival?
  if ( CONDITION == "survival" || CONDITION == "taxa" )    lnl <- - globalBiDe.equations.pSurvival.fastApprox(approxFuncs$r,approxFuncs$s,samplingProbability,0,PRESENT,PRESENT,log=TRUE)


  # multiply the probability of a descendant of the initial species
  lnl <- lnl + globalBiDe.equations.p1.fastApprox(approxFuncs$r,approxFuncs$s,samplingProbability,0,PRESENT,log=TRUE)

  # add the survival of a second species if we condition on the MRCA
  if ( MRCA == TRUE ) {
    lnl <- 2*lnl 
  } 

  # did we condition on observing n species today
  if ( CONDITION == "taxa" )    lnl <- lnl - globalBiDe.equations.pN.fastApprox(approxFuncs$r,approxFuncs$s,samplingProbability,nTaxa,0,PRESENT,MRCA,log=TRUE)

  # multiply the probability for each speciation time
  for (i in seq_len(length(times)) ) {
    lnl <- lnl + log(lambda(times[i])) + globalBiDe.equations.p1.fastApprox(approxFuncs$r,approxFuncs$s,samplingProbability,times[i],PRESENT,log=TRUE)
  }

  if ( log == FALSE ) {
    lnl <- exp(lnl)
  }
  
  return (lnl)
}



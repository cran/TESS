################################################################################
#
# sim.globalBiDe.age.R
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
# @brief Simulate a tree for a given age (either using numerical integration
#        or the analytical solution for constant rates).
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-23, version 1.0
#
# @param    n                                             scalar        number of simulations
# @param    age                                           scalar        time of the process
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    approx                                        boolean       should linear function approximation be used
# @return                                                 phylo         a random tree
#
################################################################################
sim.globalBiDe.age <- function(n,age,lambda,mu,massExtinctionTimes=c(),massExtinctionSurvivalProbabilities=c(),samplingProbability=1.0,samplingStrategy="random",MRCA=TRUE) {

  if ( length(massExtinctionTimes) != length(massExtinctionSurvivalProbabilities) ) {
    stop("Number of mass-extinction times needs to equals the number of mass-extinction survival probabilities!")
  }

  if ( samplingStrategy != "random" && samplingStrategy != "diversified" ) {
    stop("Wrong choice of argument for \"samplingStrategy\". Possible option are random|diversified.")
  }

  if ( (!is.numeric(lambda) && !inherits(lambda, "function")) || (!is.numeric(mu) && !inherits(mu, "function"))) {
    stop("Unexpected parameter types for lambda and mu!")
  }
  
  # test if we got constant values for the speciation and extinction rates
  if ( is.numeric(lambda) && is.numeric(mu) ) {
    # call simulation for constant rates (much faster)
    trees <- sim.globalBiDe.age.constant(n,age,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA)
    return (trees)
  } else {
    
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

    trees <- sim.globalBiDe.age.function(n,age,speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA)
    return (trees)
  }

}
 


################################################################################
# 
# @brief Simulate a tree for a given age.
#
# 1) Draw n times the number of taxa at the present time (Equation (10)).
# 2) For each nTaxa, draw simulate a tree using the function sim.globalBiDe.taxa.age.constant
#
# @date Last modified: 2013-01-14
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    n                                             scalar        number of simulations
# @param    age                                           scalar        time of the process
# @param    lambda                                        scalar        speciation rate function
# @param    mu                                            scalar        extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @return                                                 phylo         a random tree
#
################################################################################
sim.globalBiDe.age.constant <- function(n,age,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA) {

  # set the random taxon sampling probability
  if (samplingStrategy == "random") {
    rho <- samplingProbability
  } else {
    rho <- 1.0
  }

  
  # precompute the probabilities
  p_s <- globalBiDe.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,rho,0,age,age,log=FALSE)
  r   <- (mu-lambda)*age - log(rho)
  for (j in seq_len(length(massExtinctionTimes)) ) {
    cond <-  (0 < massExtinctionTimes[j]) & (age >= massExtinctionTimes[j])
    r  <- r - ifelse(cond, log(massExtinctionSurvivalProbabilities[j]), 0.0)
  }

  # randomly draw a number of taxa
  m <- rgeom(n, p_s * exp(r)) + 1
  if ( MRCA == TRUE ) m <- m + rgeom(n, p_s * exp(r)) + 1
  
  trees <- list()
  # for each simulation
  for ( i in 1:n ) {
    nTaxa <- m[i]
#    nTaxa <- rbinom(1,m[i],rho)

    # check if we actually have a tree
    if ( nTaxa < 2 ) {
      tree <- 1
    } else {

      # delegate the call to the simulation condition on both, age and nTaxa
      tree <- sim.globalBiDe.taxa.age.constant(1,nTaxa,age,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA)[[1]]
    }

    trees[[i]] <- tree
  } # end for each simulation

  return (trees)
}



################################################################################
# 
# @brief Simulate a tree for a given age. This is the fast approximation
#        procedure using precomputed integral functions.
#
# 1) Draw n times the number of taxa at the present time (Equation (10)).
# 2) For each nTaxa, draw simulate a tree using the function sim.globalBiDe.taxa.age.function
#
# @date Last modified: 2013-01-14
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-18, version 1.0
#
# @param    n                                             scalar        number of simulations
# @param    age                                           scalar        time of the process
# @param    lambda                                        scalar        speciation rate function
# @param    mu                                            scalar        extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @return                                                 phylo         a random tree
#
################################################################################
sim.globalBiDe.age.function <- function(n,age,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA) {

  # set the random taxon sampling probability
  if (samplingStrategy == "random") {
    rho <- samplingProbability
  } else {
    rho <- 1.0
  }

  # approximate the rate integral and the survival probability integral for fas computations
  approxFuncs <- tess.prepare.pdf(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,age,c())

  # precompute the probabilities
  p_s <- globalBiDe.equations.pSurvival.fastApprox(approxFuncs$r,approxFuncs$s,rho,0,age,age,log=FALSE)
  r   <- approxFuncs$r(age) - log(rho)

  # randomly draw a number of taxa
  nTaxa <- rgeom(n, p_s * exp(r)) + 1
  if ( MRCA == TRUE ) nTaxa <- nTaxa + rgeom(n, p_s * exp(r)) + 1
  
  trees <- list()
  # for each simulation
  for ( i in 1:n ) {

    # check if we actually have a tree
    if ( nTaxa[i] < 2 ) {
      tree <- 1
    } else {

      # delegate the call to the simulation condition on both, age and nTaxa
      tree <- sim.globalBiDe.taxa.age.function(1,nTaxa[i],age,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA=MRCA,approxFuncs$r,approxFuncs$s)[[1]]
    }

    trees[[i]] <- tree
  } # end for each simulation
  
  return (trees)
    
}

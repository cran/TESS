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
#        or linear function approximation).
#
# @date Last modified: 2012-09-23
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-23, version 1.0
#
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    age                                           scalar        time of the process
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    approx                                        boolean       should linear function approximation be used
# @return                                                 phylo         a random tree
#
################################################################################
sim.globalBiDe.age <- function(lambda,mu,massExtinctionTimes=c(),massExtinctionSurvivalProbabilities=c(),samplingProbability=1.0,age,MRCA=TRUE,approx=TRUE) {
  # test if we got constant values for the speciation and extinction rates
  if ( class(lambda) == "numeric" && class(mu) == "numeric" ) {
    # call simulation for constant rates (much faster)
    tree <- sim.globalBiDe.age.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,age,MRCA)
    return (tree)
  } else if ( class(lambda) == "function" && class(mu) == "numeric" ) {
    # should we use linear function approximation?
    if ( approx == TRUE ) {
      extinction <- function (x) mu
      tree <- sim.globalBiDe.age.function.fastApprox(lambda,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,age,MRCA,1000,1000)
      return (tree)
    } else {
      extinction <- function (x) mu
      tree <- sim.globalBiDe.age.function(lambda,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,age,MRCA,1000,1000)
      return (tree)
    }
  } else if ( class(lambda) == "numeric" && class(mu) == "function" ) {
    # should we use linear function approximation?
    if ( approx == TRUE ) {
      speciation <- function (x) lambda
      tree <- sim.globalBiDe.age.function.fastApprox(speciation,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,age,MRCA,1000,1000)
      return (tree)
    } else {
      speciation <- function (x) lambda
      tree <- sim.globalBiDe.age.function(speciation,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,age,MRCA,1000,1000)
      return (tree)
    }
  } else if ( class(lambda) == "function" && class(mu) == "function" ) {
    # should we use linear function approximation?
    if ( approx == TRUE ) {
      tree <- sim.globalBiDe.age.function.fastApprox(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,age,MRCA,1000,1000)
      return (tree)
    } else {
      tree <- sim.globalBiDe.age.function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,age,MRCA,1000,1000)
      return (tree)
    }
  } else {
    stop("Unexpected parameter types for lambda and mu!")
  }

}
 


################################################################################
# 
# @brief Simulate a tree for a given age.
#
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    lambda        scalar        speciation rate function
# @param    mu            scalar        extinction rate function
# @param    age           scalar        time of the process
# @return                 phylo         a random tree
#
################################################################################
sim.globalBiDe.age.constant <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,age,MRCA) {
  # randomly draw a number of taxa
  u <- runif(1,0,1)

  nTaxa <- 0

  # precompute the probabilities
  p_s <- globalBiDe.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,0,age,age,log=FALSE)
  r   <- (mu-lambda)*age
  
  if (MRCA == TRUE) {
    nTaxa <- 1
    while (u > 0) {
      nTaxa <- nTaxa + 1
      p <- (nTaxa-1) * (p_s * exp(r))^2 * ( 1 - p_s * exp(r))^(nTaxa-2)
      u <- u - p
    }
  } else { 
    while (u > 0) {
      nTaxa <- nTaxa + 1
      p <- p_s * exp(r) * ( 1 - p_s * exp(r))^(nTaxa-1)
      u <- u - p
    }
  }

  # check if we actually have a tree
  if ( nTaxa < 2 ) {
    return (1)
  } else {

    # delegate the call to the simulation condition on both, age and nTaxa
    tree <- sim.globalBiDe.taxa.age.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,age,MRCA)

    return (tree)
  }
}


      
################################################################################
# 
# @brief Simulate a tree for a given age.
#
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    lambda        function      speciation rate function
# @param    mu            function      extinction rate function
# @param    age           scalar        time of the process
# @return                 phylo         a random tree
#
################################################################################
sim.globalBiDe.age.function <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,age,MRCA,nTimeSteps,nBlocks) {
  # randomly draw a number of taxa
  u <- runif(1,0,1)

  nTaxa <- 0

  # precompute the probabilities
  p_s <- globalBiDe.equations.pSurvival(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,0,age,age,log=FALSE)
  r   <- globalBiDe.equations.rate(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,0,age,age)

  if (MRCA == TRUE) {
    nTaxa <- 1
    while (u > 0) {
      nTaxa <- nTaxa + 1
      p <- (nTaxa-1) * (p_s * exp(r))^2 * ( 1 - p_s * exp(r))^(nTaxa-2)
      u <- u - p
    }
  } else { 
    while (u > 0) {
      nTaxa <- nTaxa + 1
      p <- p_s * exp(r) * ( 1 - p_s * exp(r))^(nTaxa-1)
      u <- u - p
    }
  }

  # check if we actually have a tree
  if ( nTaxa < 2 ) {
    return (1)
  } else {

    # delegate the call to the simulation condition on both, age and nTaxa
    tree <- sim.globalBiDe.taxa.age.function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,age,MRCA,nTimeSteps,nBlocks)
  
    return (tree)
  }
}



################################################################################
# 
# @brief Simulate a tree for a given age. This is the fast approximation
#        procedure using precomputed integral functions.
#
# @date Last modified: 2012-09-18
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-18, version 1.0
#
# @param    lambda        function      speciation rate function
# @param    mu            function      extinction rate function
# @param    age           scalar        time of the process
# @return                 phylo         a random tree
#
################################################################################
sim.globalBiDe.age.function.fastApprox <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,age,MRCA,nTimeSteps,nBlocks) {

  #env <- parent.env(environment(NULL))
  assign("maxTime", age, envir = .TessEnv) 

  N <- get("N_DISCRETIZATION_NBLOCKS",envir=.TessEnv)

  # adjust mass-extinction times
  if (length(massExtinctionTimes) > 0) {
    for (i in 1:length(massExtinctionTimes)) {
      tmp <- N * massExtinctionTimes[i] / age
      massExtinctionTimes[i] <- (round(tmp) / N) * age
    }
  }

  # precompute the rate integral
  riv <- c()
  riv[1] <- 0
  for (i in 2:(N+1)) {
    riv[i] <- riv[i-1] + globalBiDe.equations.rate(lambda, mu, massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability, (i-2)*age/N, (i-1)*age/N, age)
  }
  assign("rateIntegralValues", riv, envir = .TessEnv)
  
  # precompute the survival integral
  siv <- c()
  siv[1] <- 0
  survival <- function(t) {
    s <- mu(t) * exp(rateIntegral(t))
    return (s)
  }
  for (i in 2:(N+1)) {
    lower <- (i-2)*age/N
    upper <- (i-1)*age/N
    siv[i] <- siv[i-1] + tess.integrate(survival,lower, upper)
    if ( length(massExtinctionTimes) > 0 ) {
      for (j in 1:length(massExtinctionTimes) ) {
        if ( lower < massExtinctionTimes[j] && upper >= massExtinctionTimes[j] ) {
          siv[i] <- siv[i] - (massExtinctionSurvivalProbabilities[j]-1)*rateIntegral(massExtinctionTimes[j])
        }
      }
    }
  }
  assign("survivalIntegralValues", siv, envir = .TessEnv)
  
  # randomly draw a number of taxa
  u <- runif(1,0,1)

  nTaxa <- 0

  # precompute the probabilities
  p_s <- globalBiDe.equations.pSurvival.fastApprox(0,age,log=FALSE)
  r   <- rateIntegral(age)

  if (MRCA == TRUE) {
    nTaxa <- 1
    while (u > 0) {
      nTaxa <- nTaxa + 1
      p <- (nTaxa-1) * (p_s * exp(r))^2 * ( 1 - p_s * exp(r))^(nTaxa-2)
      u <- u - p
    }
  } else { 
    while (u > 0) {
      nTaxa <- nTaxa + 1
      p <- p_s * exp(r) * ( 1 - p_s * exp(r))^(nTaxa-1)
      u <- u - p
    }
  }

  # check if we actually have a tree
  if ( nTaxa < 2 ) {
    return (1)
  } else {
    # delegate the call to the simulation condition on both, age and nTaxa
    tree <- sim.globalBiDe.taxa.age.function.fastApprox(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,age,MRCA=MRCA,nTimeSteps=nTimeSteps,nBlocks=nBlocks,precomputed=TRUE)
  
    return (tree)
  }
}

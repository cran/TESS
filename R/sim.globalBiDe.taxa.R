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
# @brief Simulate a tree for a given number of taxa.
#
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    nTaxa                                         scalar        number of taxa at present
# @param    max                                           scalar        the maximal time for the age
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    approx                                        boolean       should linear function approximation be used
# @return                                                 phylo         a random tree
#
################################################################################
sim.globalBiDe.taxa <- function(lambda,mu,massExtinctionTimes=c(),massExtinctionSurvivalProbabilities=c(),samplingProbability=1.0,nTaxa,max,MRCA=TRUE,approx=TRUE) {
  # test if we got constant values for the speciation and extinction rates
  if ( class(lambda) == "numeric" && class(mu) == "numeric" ) {
    # call simulation for constant rates (much faster)
    tree <- sim.globalBiDe.taxa.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,max,MRCA)
    return (tree)
  } else if ( class(lambda) == "function" && class(mu) == "numeric" ) {
    # should we use linear function approximation?
    if ( approx == TRUE ) {
      extinction <- function (x) mu
      tree <- sim.globalBiDe.taxa.function.fastApprox(lambda,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,max,MRCA,1000,1000)
      return (tree)
    } else {
      extinction <- function (x) mu
      tree <- sim.globalBiDe.taxa.function(lambda,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,max,MRCA,1000,1000)
      return (tree)
    }
  } else if ( class(lambda) == "numeric" && class(mu) == "function" ) {
    # should we use linear function approximation?
    if ( approx == TRUE ) {
      speciation <- function (x) lambda
      tree <- sim.globalBiDe.taxa.function.fastApprox(speciation,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,max,MRCA,1000,1000)
      return (tree)
    } else {
      speciation <- function (x) lambda
      tree <- sim.globalBiDe.taxa.function(speciation,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,max,MRCA,1000,1000)
      return (tree)
    }
  } else if ( class(lambda) == "function" && class(mu) == "function" ) {
    # should we use linear function approximation?
    if ( approx == TRUE ) {
      tree <- sim.globalBiDe.taxa.function.fastApprox(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,max,MRCA,1000,1000)
      return (tree)
    } else {
      tree <- sim.globalBiDe.taxa.function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,max,MRCA,1000,1000)
      return (tree)
    }
  } else {
    stop("Unexpected parameter types for lambda and mu!")
  }

}



################################################################################
# 
# @brief Simulate a tree for a given number of taxa.
#
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    lambda        scalar      speciation rate function
# @param    mu            scalar      extinction rate function
# @param    nTaxa         scalar        number of taxa in the tree at the present time
# @return                 phylo         a random tree
#
################################################################################
sim.globalBiDe.taxa.constant <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,max,MRCA) {

  div <- lambda - mu

  # randomly draw the age of the tree
  u <- runif(1,0,1)

  T <- log( (-lambda + mu * u^(1.0/nTaxa) ) / (lambda * (-1 + u^(1.0/nTaxa)))) / div;


  # delegate the call to the simulation condition on both, age and nTaxa
  tree <- sim.globalBiDe.taxa.age.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,T,MRCA)

  return (tree)
}



################################################################################
# 
# @brief Simulate a tree for a given number of taxa.
#
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    lambda        function      speciation rate function
# @param    mu            function      extinction rate function
# @param    nTaxa         scalar        number of taxa in the tree at the present time
# @return                 phylo         a random tree
#
################################################################################
sim.globalBiDe.taxa.function <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,max,MRCA,nTimeSteps,nBlocks) {
  # compute the normalizing constant first
  integrand <- function(x) {
    prob <- globalBiDe.equations.pN(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,0,x,log=FALSE)
    return (prob)
  }
  normalizingConstant <- tess.integrate(integrand,lower=0,upper=max)

  # randomly draw the age of the tree
  u <- runif(1,0,1)
  high <- max
  low <- 0
  f_low <- 0
  f_high <- 1
  repeat {
    # stoping condition
    if ( (f_high - f_low) < 1E-4 || (high - low) < 1E-4) {
      break
    }
    middle <- low + (high - low) / 2

    f_middle <- f_low + tess.integrate(integrand,lower=low,upper=middle) / normalizingConstant
    
    if ( f_middle > u ) {
      high <- middle
      f_high <- f_middle
    } else {
      low <- middle
      f_low <- f_middle
    } 
  }

  T <- low + (u - f_low) * (high - low) / (f_high - f_low)

  # delegate the call to the simulation condition on both, age and nTaxa
  tree <- sim.globalBiDe.taxa.age.function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,T,nTimeSteps,nBlocks)
  
  return (tree)
}




################################################################################
# 
# @brief Simulate a tree for a given number of taxa. This is the fast approximation
#        procedure using precomputed integral functions.
#
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    lambda        function      speciation rate function
# @param    mu            function      extinction rate function
# @param    nTaxa         scalar        number of taxa in the tree at the present time
# @return                 phylo         a random tree
#
################################################################################
sim.globalBiDe.taxa.function.fastApprox <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,max,MRCA,nTimeSteps,nBlocks) {

  assign("maxTime", max, envir = .TessEnv)
  N <- get("N_DISCRETIZATION_NBLOCKS",envir=.TessEnv)

  # adjust mass-extinction times
  if (length(massExtinctionTimes) > 0) {
    for (i in 1:length(massExtinctionTimes)) {
      tmp <- N * massExtinctionTimes[i] / max
      massExtinctionTimes[i] <- (round(tmp) / N) * max
    }
  }


  # precompute the rate integral
  riv <- c()
  riv[1] <- 0
  for (i in 2:(N+1)) {
    riv[i] <- riv[i-1] + globalBiDe.equations.rate(lambda, mu, massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability, (i-2)*max/N, (i-1)*max/N, max)
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
    lower <- (i-2)*max/N
    upper <- (i-1)*max/N
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

  # compute the cumulative distribution function for the time of the process
  integrand <- function(x) {
    prob <- globalBiDe.equations.pN.fastApprox(nTaxa,0,x,log=FALSE)
    return (prob)
  }
  originIntegralValues <- c()
  originIntegralValues[1] <- 0
  for (i in 2:(N+1)) {
    lower <- (i-2)*max/N
    upper <- (i-1)*max/N
    originIntegralValues[i] <- originIntegralValues[i-1] + tess.integrate(integrand,lower,upper)
  }

  originIntegral <- function(x) {
    if ( x == 0 ) return ( 0.0 )
  
    index <- ceiling(x*N/max)
    y <- 0
    if ( index == x*N/max) {
      y <- originIntegralValues[index + 1]
    } else {
      y <- originIntegralValues[index] + (originIntegralValues[index+1] - originIntegralValues[index])*( x/max - ((index-1)/N))
    }
  
    return (y)
  }

  normalizingConstant <- originIntegral(max)
  # randomly draw the age of the tree
  u <- runif(1,0,1)
  high <- max
  low <- 0
  f_low <- 0
  f_high <- 1
  repeat {
    # stoping condition
    if ( (f_high - f_low) < 1E-4 || (high - low) < 1E-4) {
      break
    }
    middle <- low + (high - low) / 2

    f_middle <- originIntegral(middle) / normalizingConstant
    
    if ( f_middle > u ) {
      high <- middle
      f_high <- f_middle
    } else {
      low <- middle
      f_low <- f_middle
    } 
  }

  T <- low + (u - f_low) * (high - low) / (f_high - f_low)

  # delegate the call to the simulation condition on both, age and nTaxa
  tree <- sim.globalBiDe.taxa.age.function.fastApprox(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,age=T,MRCA,nTimeSteps,nBlocks,precomputed=FALSE)
  
  return (tree)
}

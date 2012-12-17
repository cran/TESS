################################################################################
#
# globalBiDe.equations.R
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
# @date Last modified: 2012-11-05
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    n                                             scalar        number of simulations
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    nTaxa                                         scalar        number of taxa at present
# @param    age                                           scalar        the age
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    approx                                        boolean       should linear function approximation be used
# @return                                                 phylo         a random tree
#
################################################################################
sim.globalBiDe.taxa.age <- function(n,nTaxa,age,lambda,mu,massExtinctionTimes=c(),massExtinctionSurvivalProbabilities=c(),samplingProbability=1.0,MRCA=TRUE) {

  if ( length(massExtinctionTimes) != length(massExtinctionSurvivalProbabilities) ) {
    stop("Number of mass-extinction times needs to equals the number of mass-extinction survival probabilities!")
  }

  if ( (!inherits(lambda, "numeric") && !inherits(lambda, "function")) || (!inherits(mu, "numeric") && !inherits(mu, "function"))) {
    stop("Unexpected parameter types for lambda and mu!")
  }
  
  # test if we got constant values for the speciation and extinction rates
  if ( inherits(lambda, "numeric") && inherits(mu, "numeric") ) {
    # call simulation for constant rates (much faster)
    trees <- sim.globalBiDe.taxa.age.constant(n,nTaxa,age,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA)
    return (trees)
  } else {
    # convert the speciation rate into a function if necessary
    if ( inherits(lambda, "numeric") ) {
      speciation <- function (x) lambda
    } else {
      speciation <- lambda
    }
    # convert the extinction rate into a function if necessary
    if ( inherits(mu, "numeric") ) {
      extinction <- function (x) mu
    } else {
      extinction <- mu
    }

    approxFuncs <- tess.prepare.pdf(speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,age,c())
    
    trees <- sim.globalBiDe.taxa.age.function(n,nTaxa,age,speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA,approxFuncs$r,approxFuncs$s)
    
    return (trees)
  } 

}


################################################################################
# 
# @brief Simulate a tree conditioned on the number of taxa and
#        the age of the tree.
#
#
# Simulate the n speciation times using the inverse cdf given in Equation (9).
#
# @date Last modified: 2012-12-17
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    n                                             scalar        number of simulations
# @param    nTaxa                                         scalar        number of taxa at present
# @param    age                                           scalar        the age
# @param    lambda                                        scalar        speciation rate function
# @param    mu                                            scalar        extinction rate function
# @return                                                 phylo         a random tree
#
################################################################################
sim.globalBiDe.taxa.age.constant <- function(n,nTaxa,age,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA) {

  # if we have mass-extinction events, then we use the function-integration approach
  # because the closed form solutions only apply to models without mass-extinctions.
  if ( length(massExtinctionTimes) > 0 ) {
    speciation <- function(x) lambda
    extinction <- function(x) mu
    approxFuncs <- tess.prepare.pdf(speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,age,c())    
    trees <- sim.globalBiDe.taxa.age.function(n,nTaxa,age,speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA,approxFuncs$r,approxFuncs$s)
    return (trees)
  } else {

    div    <- lambda - mu
    origin <- age 
    b <- samplingProbability * lambda
    d <- mu - lambda *(1-samplingProbability)

    # start to simulate n trees
    trees <- list()
    for ( j in 1:n ) {

      # now draw the random time speciation times
      times         <- c()
      x <- 1
      if (MRCA == TRUE) x <- 2
      if ( (nTaxa-x) > 0) {
        
        # for each speciation time
        u           <- runif(nTaxa-x,0,1)
        times       <- 1/(b-d)*log((b - d*exp((-b+d)*age) -d*(1-exp((-b+d)*age)) *u )/(b - d*exp((-b+d)*age) -b*(1.0-exp((-b+d)*age)) *u )   )  

        # now we have the vector of speciation times, which just need to be converted to a phylogeny
        times         <- c(sort(times,decreasing=FALSE),age)
      } else {
        times         <- age
      }
      tree          <- tess.create.phylo(times, root = !MRCA)

      trees[[j]]    <- tree
    } # end for each simulation

    return (trees)
  }
}




################################################################################
# 
# @brief Simulate a tree conditioned on the number of taxa and
#        the age of the tree. This is the fast approximation
#        procedure using precomputed integral functions.
#
# Simulate the n speciation times using the inverse cdf given in Equation (9).
#
# @date Last modified: 2012-12-17
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-18, version 1.0
#
# @param    n                                             scalar        number of simulations
# @param    nTaxa                                         scalar        number of taxa at present
# @param    age                                           scalar        the age
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    approx                                        boolean       should linear function approximation be used
# @return                                                 phylo         a random tree
#
################################################################################
sim.globalBiDe.taxa.age.function <- function(n,nTaxa,age,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA,r,s) {

  tess.pdf <- function(x) lambda(x) * globalBiDe.equations.p1.fastApprox(r, s, samplingProbability, x, age, log=FALSE)

#  inverse <- n * nTaxa > 10000
  inverse <- TRUE

  if ( inverse ) {
    n2 <- 1001
    obj.pdf <- function(t, state, pars) list(tess.pdf(t))
    times <- seq(0, age, length=n2)
    ## This step is slow for large n - perhaps 1/2s for 1000 points
    zz <- lsoda(0, times, obj.pdf, tcrit=age)[,2]
    const <- last(zz)
    zz <- zz / const ## Normalise
    icdf <- approxfun(zz, times) ## Interpolate
    speciationEventSim <- function(n)  icdf(runif(n))
  } else {
    sup <- find.max(tess.pdf, age, 101)
    speciationEventSim <- function(n)  rejection.sample.simple(n, tess.pdf, c(0, age), sup)
  }

  # start to simulate n trees
  trees <- list()
  for ( j in 1:n ) {
    
    # now draw the random time speciation times
    times         <- c()
    x <- 1
    if (MRCA == TRUE) x <- 2
    if ( (nTaxa-x) > 0) {
      times <- speciationEventSim(nTaxa-x)

      # now we have the vector of speciation times, which just need to be converted to a phylogeny
      times         <- c(age - sort(times,decreasing=TRUE),age)
    } else {
      times         <- age
    }

    # create the tree
    tree          <- tess.create.phylo( times, root = !MRCA )

    trees[[j]]    <- tree
  } #end for each simulation
  
  return (trees)
}

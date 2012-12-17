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
# @date Last modified: 2012-11-05
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    n                                             scalar        number of simulations
# @param    nTaxa                                         scalar        number of taxa at present
# @param    max                                           scalar        the maximal time for the age
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    MRCA                                          boolean       do we start the tree with the MRCA (two species)?
# @param    t_crit                                        vector        critical times when jumps in the rate functions occur
# @return                                                 list          list of random trees (type phylo)
#
################################################################################
sim.globalBiDe.taxa <- function(n,nTaxa,max,lambda,mu,massExtinctionTimes=c(),massExtinctionSurvivalProbabilities=c(),samplingProbability=1.0,MRCA=TRUE,t_crit=c()) {

  if ( length(massExtinctionTimes) != length(massExtinctionSurvivalProbabilities) ) {
    stop("Number of mass-extinction times needs to equals the number of mass-extinction survival probabilities!")
  }

  if ( (!inherits(lambda, "numeric") && !inherits(lambda, "function")) || (!inherits(mu, "numeric") && !inherits(mu, "function"))) {
    stop("Unexpected parameter types for lambda and mu!")
  }
  
  # test if we got constant values for the speciation and extinction rates
  if ( class(lambda) == "numeric" && class(mu) == "numeric" ) {
    # call simulation for constant rates (much faster)
    trees <- sim.globalBiDe.taxa.constant(n,nTaxa,max,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA)
    return (trees)
  } else  {
    
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
    
    trees <- sim.globalBiDe.taxa.function(n,nTaxa,max,speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA,t_crit)
    return (trees)
  }

}



################################################################################
# 
# @brief Simulate a tree for a given number of taxa.
#
# 1) Simulate the time of the process using Monte Carlo sampling, see Equation (11).
# 2) Simulate the tree by calling sim.globalBiDe.taxa.age.constant
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    n                                             scalar        number of simulations
# @param    nTaxa                                         scalar        number of taxa at present
# @param    max                                           scalar        the maximal time for the age
# @param    lambda                                        scalar        speciation rate function
# @param    mu                                            scalar        extinction rate function
# @param    samplingProbability                           scalar        probability of random sampling at present
# @return                                                 list          list of random trees (type phylo)
#
################################################################################
sim.globalBiDe.taxa.constant <- function(n,nTaxa,max,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA) {

  if ( length(massExtinctionTimes) > 0 ) {
    # compute the cumulative distribution function for the time of the process
    pdf <- function(x) globalBiDe.equations.pN.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,0,x,MRCA,log=FALSE)

    # preparations for the rejection sampling
    xx <- seq(0, max, length=101)
    yy <- pdf(xx)
    i <- which.max(yy)
    sup <- optimize(pdf, range(na.omit(xx[(i-1):(i+1)])), maximum=TRUE)$max

    # randomly draw the age of the tree
    # use rejection sampling
    T <- rejection.sample.simple(n,pdf,c(0,max),sup)
  } else {
    u <- runif(n,0,1)
    T <- log((-lambda * samplingProbability - lambda * u^(1/nTaxa) + mu * u^(1/nTaxa) + lambda * samplingProbability * u^(1/nTaxa))/(lambda * samplingProbability * (-1 + u^(1/nTaxa)))) / (lambda - mu)
  }
    
  trees <- list()
  start <- Sys.time()
  # for each simulation
  for ( i in 1:n ) {

    # delegate the call to the simulation condition on both, age and nTaxa
    tree <- sim.globalBiDe.taxa.age.constant(1,nTaxa,T[i],lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA)
  
    # add the new tree to the list
    trees[[i]] <- tree[[1]]
  }

  return (trees)
}



################################################################################
# 
# @brief Simulate a tree for a given number of taxa.
#
# 1) Simulate the time of the process using Monte Carlo sampling, see Equation (11).
# 2) Simulate the tree by calling sim.globalBiDe.taxa.age.function
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    lambda        function      speciation rate function
# @param    mu            function      extinction rate function
# @param    nTaxa         scalar        number of taxa in the tree at the present time
# @return                 phylo         a random tree
#
################################################################################
sim.globalBiDe.taxa.function <- function(n,nTaxa,max,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA,t_crit=c()) {

  # approximate the rate integral and the survival probability integral for fas computations
  approxFuncs <- tess.prepare.pdf(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,max,t_crit)

  # compute the cumulative distribution function for the time of the process
  pdf <- function(x) globalBiDe.equations.pN.fastApprox(approxFuncs$r,approxFuncs$s,samplingProbability,nTaxa,0,x,MRCA,log=FALSE)

  # preparations for the rejection sampling
  xx <- seq(0, max, length=101)
  yy <- pdf(xx)
  i <- which.max(yy)
  sup <- optimize(pdf, range(na.omit(xx[(i-1):(i+1)])), maximum=TRUE)$max

  # randomly draw the age of the tree
  # use rejection sampling
  T <- rejection.sample.simple(n,pdf,c(0,max),sup)
  
  trees <- list()
  start <- Sys.time()
  # for each simulation
  for ( i in 1:n ) {

    # delegate the call to the simulation condition on both, age and nTaxa
    tree <- sim.globalBiDe.taxa.age.function(1,nTaxa,age=T[i],lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA,approxFuncs$r,approxFuncs$s)

    # add the new tree to the list
    trees[[i]] <- tree[[1]]
  }
  
  return (trees)
}

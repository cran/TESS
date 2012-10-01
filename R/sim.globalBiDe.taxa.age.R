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
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
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
sim.globalBiDe.taxa.age <- function(lambda,mu,massExtinctionTimes=c(),massExtinctionSurvivalProbabilities=c(),samplingProbability=1.0,nTaxa,age,MRCA=TRUE,approx=TRUE) {
  # test if we got constant values for the speciation and extinction rates
  if ( class(lambda) == "numeric" && class(mu) == "numeric" ) {
    # call simulation for constant rates (much faster)
    tree <- sim.globalBiDe.taxa.age.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,age,MRCA)
    return (tree)
  } else if ( class(lambda) == "function" && class(mu) == "numeric" ) {
    # should we use linear function approximation?
    if ( approx == TRUE ) {
      extinction <- function (x) mu
      tree <- sim.globalBiDe.taxa.age.function.fastApprox(lambda,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,age,MRCA,1000,1000)
      return (tree)
    } else {
      extinction <- function (x) mu
      tree <- sim.globalBiDe.taxa.age.function(lambda,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,age,MRCA,1000,1000)
      return (tree)
    }
  } else if ( class(lambda) == "numeric" && class(mu) == "function" ) {
    # should we use linear function approximation?
    if ( approx == TRUE ) {
      speciation <- function (x) lambda
      tree <- sim.globalBiDe.taxa.age.function.fastApprox(speciation,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,age,MRCA,1000,1000)
      return (tree)
    } else {
      speciation <- function (x) lambda
      tree <- sim.globalBiDe.taxa.age.function(speciation,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,age,MRCA,1000,1000)
      return (tree)
    }
  } else if ( class(lambda) == "function" && class(mu) == "function" ) {
    # should we use linear function approximation?
    if ( approx == TRUE ) {
      tree <- sim.globalBiDe.taxa.age.function.fastApprox(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,age,MRCA,1000,1000)
      return (tree)
    } else {
      tree <- sim.globalBiDe.taxa.age.function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,age,MRCA,1000,1000)
      return (tree)
    }
  } else {
    stop("Unexpected parameter types for lambda and mu!")
  }

}


################################################################################
# 
# @brief Simulate a tree conditioned on the number of taxa and
#        the age of the tree.
#
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    lambda                                        scalar        speciation rate function
# @param    mu                                            scalar        extinction rate function
# @param    nTaxa                                         scalar        number of taxa at present
# @param    T                                             scalar        the age
# @return                                                 phylo         a random tree
#
################################################################################
sim.globalBiDe.taxa.age.constant <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,T,MRCA,nTimeSteps=1000,nBlocks=1000) {

  div    <- lambda - mu
  origin <- T 

  # now draw the random time speciation times
  times         <- c()
  x <- 1
  if (MRCA == TRUE) x <- 2
  if ( (nTaxa-x) > 0) {
    for (i in 1:(nTaxa-x)) { # for each speciation time
      u           <- runif(1,0,1)
      times[i]    <- 1.0/div * log((lambda - mu * exp((-div)*origin) - mu * (1.0 - exp((-div) * origin)) * u )/(lambda - mu * exp((-div) * origin) - lambda * (1.0 - exp(( -div ) * origin)) * u ) )
    } # end loop over each speciation event
  }

  # now we have the vector of speciation times, which just need to be converted to a phylogeny
  times         <- c(T,T - sort(times))

  # create a list with all taxa
  taxa          <- list()
  for (i in 1:nTaxa) {
    taxon       <- list(name=sprintf("Tip_%i",i),time=0.0)
    taxa[[i]]   <- taxon
  }

  # build tree by merging taxa
  i            <- length(times)
  while (length(taxa) >= 2 ) {
    # get left child
    index       <- ceiling(length(taxa)*runif(1,0,1))
    left        <- taxa[[index]]
    taxa[index] <- NULL
    # get right child
    index       <- ceiling(length(taxa)*runif(1,0,1))
    right       <- taxa[[index]]
    taxa[index] <- NULL
    # build the new parent
    t           <- times[i]
    parent      <- list(name=sprintf("(%s:%f,%s:%f)",left$name,t-left$time,right$name,t-right$time),time=t)
    taxa[[length(taxa)+1]]    <- parent
    i           <- i - 1
  }

  newick <- ""
  if (MRCA == FALSE) {
    newick      <- sprintf("%s:%f;",taxa[[1]]$name,T-taxa[[1]]$time)
  } else { 
    newick      <- sprintf("%s;",taxa[[1]]$name)
  }

  # construct a tree from the phylo class
  tree          <- read.tree(text=newick)

  return (tree)
}



################################################################################
# 
# @brief Simulate a tree conditioned on the number of taxa and
#        the age of the tree.
#
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
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
sim.globalBiDe.taxa.age.function <- function(lambda,mu,massExtinctionTimes=c(),massExtinctionSurvivalProbabilities=c(),samplingProbability=1.0,nTaxa,T,MRCA,nTimeSteps=1000,nBlocks=1000) {

  age <- T
  # we need to compute the inverse of the cumulative distribution function first
  divisor       <- 1 - globalBiDe.equations.pSurvival(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,0,T,T,log=FALSE)*exp(globalBiDe.equations.rate(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,0,T,T))
  normalizingConstant <- 1.0
  integrand     <- function(x) {
    a <- lambda(x) * globalBiDe.equations.p1(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,x,T,log=FALSE) / divisor / normalizingConstant
    return (a)
  }
  normalizingConstant <- tess.integrate(integrand,lower=0,upper=age)
  
  index         <- 2
  lastTime      <- 0
  lastUsedTime  <- 0
  # we break the inverse cumulative distribution function into 'nBlocks' equally spaced discrete intervals
  values        <- c()

  # the first value is obviously 0
  values[1]     <- 0
  intervalSize  <- T/nTimeSteps
  F_remainder   <- 0
  for (i in 1:nTimeSteps) { # loop over each time step
    # compute the time at step i
    time           <- i*intervalSize

    # compute the integral of the density function between the last time and the new time
    F_N            <- F_remainder + tess.integrate(integrand,lower=lastTime,upper=time)
    blocks         <- floor(F_N*nBlocks)
    F_remainder    <- F_N - (blocks/nBlocks)
    if (blocks > 0) {
      # compute the time increase for each small discrete step of the inverse cumulative distribution function
      increasePerBlock <- (time - lastUsedTime) / blocks
      for (j in 1:blocks) { # loop over each block
        # store the time of the inverse cumulative distribution function f(x) = t, where x = index/nBlocks
        values[index] <- values[index-1] + increasePerBlock

        # increment the index
        index         <- index + 1
      } # end loop over each new block
      lastUsedTime <- time
    }
    lastTime       <- time
  } # end loop over each time step
  # add the final stoppping time
  values[nBlocks+1] <- T

  # now draw the random time speciation times
  times         <- c()
  x <- 1
  if (MRCA == TRUE) x <- 2
  if ( (nTaxa-x) > 0) {
    for (i in 1:(nTaxa-x)) { # for each speciation time
      u           <- runif(1,0,1)
      index       <- ceiling(u*nBlocks)
      times[i]    <- values[index] + (values[index+1] - values[index])*( u - ((index-1)/nBlocks))
    } # end loop over each speciation event
  }

  # now we have the vector of speciation times, which just need to be converted to a phylogeny
  times         <- c(T,T - sort(times))
  
  # create a list with all taxa
  taxa          <- list()
  for (i in 1:nTaxa) {
    taxon       <- list(name=sprintf("Tip_%i",i),time=0.0)
    taxa[[i]]   <- taxon
  }

  # build tree by merging taxa
  i            <- length(times)
  while (length(taxa) >= 2 ) {
    # get left child
    index       <- ceiling(length(taxa)*runif(1,0,1))
    left        <- taxa[[index]]
    taxa[index] <- NULL
    # get right child
    index       <- ceiling(length(taxa)*runif(1,0,1))
    right       <- taxa[[index]]
    taxa[index] <- NULL
    # build the new parent
    t           <- times[i]
    parent      <- list(name=sprintf("(%s:%f,%s:%f)",left$name,t-left$time,right$name,t-right$time),time=t)
    taxa[[length(taxa)+1]]    <- parent
    i           <- i - 1
  }

  newick <- ""
  if (MRCA == FALSE) {
    newick      <- sprintf("%s:%f;",taxa[[1]]$name,age-taxa[[1]]$time)
  } else { 
    newick      <- sprintf("%s;",taxa[[1]]$name)
  }

  # construct a tree from the phylo class
  tree          <- read.tree(text=newick)
  
  return (tree)
}



################################################################################
# 
# @brief Simulate a tree conditioned on the number of taxa and
#        the age of the tree. This is the fast approximation
#        procedure using precomputed integral functions.
#
# @date Last modified: 2012-09-18
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-18, version 1.0
#
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
sim.globalBiDe.taxa.age.function.fastApprox <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,age,MRCA,nTimeSteps,nBlocks,precomputed=FALSE) {

  T <- age
  if ( precomputed == FALSE ) {
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
  }
  
  # we need to compute the inverse of the cumulative distribution function first
  divisor       <- 1 - globalBiDe.equations.pSurvival.fastApprox(0,age,age,log=FALSE)*exp(rateIntegral(age))
  normalizingConstant <- 1.0
  integrand     <- function(x) {
    a <- lambda(x) * globalBiDe.equations.p1.fastApprox(x,age,log=FALSE) / divisor / normalizingConstant
    return (a)
  }
  normalizingConstant <- tess.integrate(integrand,lower=0,upper=age)
  
  index         <- 2
  lastTime      <- 0
  lastUsedTime  <- 0
  # we break the inverse cumulative distribution function into 'nBlocks' equally spaced discrete intervals
  values        <- c()

  # the first value is obviously 0
  values[1]     <- 0
  intervalSize  <- age/nTimeSteps
  F_remainder   <- 0
  for (i in 1:nTimeSteps) { # loop over each time step
    # compute the time at step i
    time           <- i*intervalSize

    # compute the integral of the density function between the last time and the new time
    F_N            <- F_remainder + tess.integrate(integrand,lower=lastTime,upper=time)
    blocks         <- floor(F_N*nBlocks)
    F_remainder    <- F_N - (blocks/nBlocks)
    if (blocks > 0) {
      # compute the time increase for each small discrete step of the inverse cumulative distribution function
      increasePerBlock <- (time - lastUsedTime) / blocks
      for (j in 1:blocks) { # loop over each block
        # store the time of the inverse cumulative distribution function f(x) = t, where x = index/nBlocks
        values[index] <- values[index-1] + increasePerBlock

        # increment the index
        index         <- index + 1
      } # end loop over each new block
      lastUsedTime <- time
    }
    lastTime       <- time
  } # end loop over each time step
  # add the final stoppping time
  values[nBlocks+1] <- age

  # now draw the random time speciation times
  times         <- c()
  x <- 1
  if (MRCA == TRUE) x <- 2
  if ( (nTaxa-x) > 0) {
    for (i in 1:(nTaxa-x)) { # for each speciation time
      u           <- runif(1,0,1)
      index       <- ceiling(u*nBlocks)
      times[i]    <- values[index] + (values[index+1] - values[index])*( u - ((index-1)/nBlocks))
    } # end loop over each speciation event
  }
  
  # now we have the vector of speciation times, which just need to be converted to a phylogeny
  times         <- c(age,age - sort(times))
  
  # create a list with all taxa
  taxa          <- list()
  for (i in 1:nTaxa) {
    taxon       <- list(name=sprintf("Tip_%i",i),time=0.0)
    taxa[[i]]   <- taxon
  }

  # build tree by merging taxa
  i            <- length(times)
  while (length(taxa) >= 2 ) {
    # get left child
    index       <- ceiling(length(taxa)*runif(1,0,1))
    left        <- taxa[[index]]
    taxa[index] <- NULL
    # get right child
    index       <- ceiling(length(taxa)*runif(1,0,1))
    right       <- taxa[[index]]
    taxa[index] <- NULL
    # build the new parent
    t           <- times[i]
    parent      <- list(name=sprintf("(%s:%f,%s:%f)",left$name,t-left$time,right$name,t-right$time),time=t)
    taxa[[length(taxa)+1]]    <- parent
    i           <- i - 1
  }

  newick <- ""
  if (MRCA == FALSE) {
    newick      <- sprintf("%s:%f;",taxa[[1]]$name,age-taxa[[1]]$time)
  } else { 
    newick      <- sprintf("%s;",taxa[[1]]$name)
  }

  # construct a tree from the phylo class
  tree          <- read.tree(text=newick)
  
  return (tree)
}

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
# @brief  Calculate the probability of survival in the interval [t_low,t_high].
#
# See Equation (2) from Hoehna, S., 2013, Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes
#
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    t_low                                         scalar        starting time
# @param    t_high                                        scalar        end time
# @param    T                                             scalar        present time (time goes forward and the origin/MRCA might be 0)
# @param    log                                           bool          if in log-scale
# @return                                                 scalar        probability of survival in [t,tau]
#
################################################################################
globalBiDe.equations.pSurvival.constant <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t_low,t_high,T,log=FALSE) {

  # compute the rate
  rate <- mu - lambda
   
  # add mass-extinction
  part1 <- 0.0
  part2 <- 0.0
  accumulatedMassExtinction <- 1.0
  prev_time <- t_low
  if ( length(massExtinctionTimes) > 0 ) {
    for (j in 1:length(massExtinctionTimes) ) {
      cond <-  (t_low < massExtinctionTimes[j]) & (t_high >= massExtinctionTimes[j])
      part1  <- part1 + ifelse(cond,  accumulatedMassExtinction * ( exp(rate*massExtinctionTimes[j]) - exp(rate*prev_time)), 0.0)
      prev_time <- ifelse(cond, massExtinctionTimes[j], prev_time)
      accumulatedMassExtinction <- accumulatedMassExtinction / ifelse(cond, massExtinctionSurvivalProbabilities[j], 1.0)
      part2  <- part2 + ifelse(cond, (massExtinctionSurvivalProbabilities[j]-1)*exp( rate*(massExtinctionTimes[j]) ) * accumulatedMassExtinction, 0.0) 
    }
  }

  part1 <- part1 + accumulatedMassExtinction * ( exp(rate*t_high) - exp(rate*prev_time))

  # add sampling
  cond <- (t_low < T) & (t_high >= T)
  accumulatedMassExtinction <- accumulatedMassExtinction / ifelse(cond, samplingProbability, 1.0)
  part2 <- part2 + ifelse(cond, (samplingProbability-1)*exp( rate*(T) ) * accumulatedMassExtinction, 0.0)

  den <- 1.0 + part1 * (mu/rate) / exp(rate*t_low) - part2 / exp(rate*t_low)

  res <- 1.0 / den

  if ( log == TRUE ) {
    res <- log( res )
  }

  return (res)

}


################################################################################
# 
# @brief  Calculate the probability of survival in the interval [t_low,t_high]
#         using a faster approximation.
#
# See Equation (2) from Hoehna, S., 2013, Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes
#
# @date Last modified: 2013-03-01
# @author Sebastian Hoehna
# @version 1.2
# @since 2012-09-18, version 1.0
#
# @param    r             function      rate integral function
# @param    s             function      survival probability function
# @param    t_low         scalar        starting time
# @param    t_high        scalar        end time
# @param    T             scalar        present time (time goes forward and the origin/MRCA might be 0)
# @param    log           bool          if in log-scale
# @return                 scalar        probability of survival in [t,tau]
#
################################################################################
globalBiDe.equations.pSurvival.fastApprox <- function(r, s, rho, t_low, t_high, T, log=FALSE) {
#  den <- (1 + ( s(t_high) - s(t_low) ) / exp(r(t_low)))
  b <- s(t_low)
  c <- exp(r(t_high) - r(t_low))
  d <- s(t_high)

  den <- (1 + b - c*d)

  # add sampling
  cond <- (t_low < T) & (t_high >= T)
  den[cond] <- den[cond] - (rho-1)*exp( r(T) - r(t_low[cond]) - log(rho) )

#  res <- ifelse( is.finite(den) & den > 0, 1/den, 0 )
  res <- rep(0,length(den))
  tmp <- is.finite(den) & den > 0
  res[tmp] <- 1/den
  
  if ( log == TRUE )
    res <- log( res )

  return (res)
}


################################################################################
# 
# @brief  Calculate the probability of n lineage existing at time t
#         if we started with 1 (or 2) lineage at time s.
#
# 
# @see    Equation (3), (5) and (12) in Hoehna, S.: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    i                                             scalar        number of lineages
# @param    s                                             scalar        start time
# @param    t                                             scalar        stop time
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    log                                           bool          if in log-scale
# @return                                                 scalar        log-probability
#
################################################################################
globalBiDe.equations.pN.constant <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,i,s,t,MRCA=FALSE,log=FALSE) {
  if ( i < 1 ) { # we assume conditioning on survival
    p <- 0
  } else if (i == 1) {
    if ( MRCA == TRUE ) { # we assume conditioning on survival of the two species
      p <- 0
    } else {
      p   <-  globalBiDe.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,s,t,t,log=TRUE) + (mu-lambda)*(t-s) - log(samplingProbability)
    }
  } else {
    p_s <- globalBiDe.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,s,t,t,log=FALSE)
    r   <- (mu-lambda)*(t-s) - log(samplingProbability)
    for (j in seq_len(length(massExtinctionTimes)) ) {
      cond <-  (s < massExtinctionTimes[j]) & (t >= massExtinctionTimes[j])
      r  <- r - ifelse(cond, log(massExtinctionSurvivalProbabilities[j]), 0.0)
    }
    e   <- p_s * exp(r)
    e[e > 1] <- 1

    if ( MRCA == FALSE ) {
      p   <- log(p_s) + r + log( 1 - e) * (i-1)
    } else {
      p   <- log(i-1) + 2*log(p_s) + 2*r + log( 1 - e) * (i-2)
    }
  }

   if ( log == FALSE ) {
     p <- exp( p )
   }

  return (p)
}


################################################################################
# 
# @brief  Calculate the probability of n lineage existing at time t
#         if we started with 1 (or 2) lineage(s) at time s. This is the fast
#         approximation.
#
# @see    Equation (3), (5) and (12) in Hoehna, S.: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-18, version 1.0
#
# @param    rate                                          function      rate integral function
# @param    surv                                          function      survival integral function
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    i                                             scalar        number of lineages
# @param    s                                             scalar        start time
# @param    t                                             scalar        stop time
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    log                                           bool          if in log-scale
# @return                                                 scalar        log-probability
#
################################################################################
globalBiDe.equations.pN.fastApprox <- function(rate,surv,samplingProbability,i,s,t,MRCA=FALSE,log=FALSE) {

  if ( i < 1 ) { # we assume conditioning on survival
    p <- 0
  } else if (i == 1) {
    if ( MRCA == TRUE ) { # we assume conditioning on survival of the two species
      p <- 0
    } else {
      p   <- globalBiDe.equations.pSurvival.fastApprox(rate,surv,samplingProbability,s,t,t,log=TRUE) + rate(t) - rate(s) - log(samplingProbability)
    }
  }
  else {
    p_s <- globalBiDe.equations.pSurvival.fastApprox(rate,surv,samplingProbability,s,t,t,log=FALSE)
    r   <- rate(t) - rate(s) - log(samplingProbability)
    e   <- p_s * exp(r)
    e[e > 1] <- 1

    if ( MRCA == FALSE ) {
      p   <- log(p_s) + r + log( 1 - e) * (i-1)
    } else {
      p   <- log(i-1) + 2*log(p_s) + 2*r + log( 1 - e) * (i-2)
    }

  }

  if ( log == FALSE ) {
    p <- exp( p )
  }


  return (p)
}


################################################################################
# 
# @brief  Calculate the expected number of taxa when the process start at time
#         s with 1 (or two) species and ends at time t.
#
# 
# @see    Equation (3), (5) and (12) in Hoehna, S.: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    s                                             scalar        start time
# @param    t                                             scalar        stop time
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @return                                                 scalar        expected number of taxa
#
################################################################################
globalBiDe.equations.nTaxa.expected.constant <- function(s,t,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA=FALSE) {

  p_s <- globalBiDe.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,s,t,t,log=FALSE)
  r   <- (mu-lambda)*(t-s) - log(samplingProbability)
  for (j in seq_len(length(massExtinctionTimes)) ) {
    cond <-  (s < massExtinctionTimes[j]) & (t >= massExtinctionTimes[j])
    r  <- r - ifelse(cond, log(massExtinctionSurvivalProbabilities[j]), 0.0)
  }
  e   <- p_s * exp(r)
  e[e > 1] <- 1

  if ( MRCA == FALSE ) {
    n  <- 1.0 / e
  } else {
    n   <- 2.0 / e
  }

  return (n)
}


################################################################################
# 
# @brief  Calculate the expected number of taxa when the process start at time
#         s with 1 (or two) species and ends at time t.
#
# @see    Equation (3), (5) and (12) in Hoehna, S.: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-18, version 1.0
#
# @param    rate                                          function      rate integral function
# @param    surv                                          function      survival integral function
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    s                                             scalar        start time
# @param    t                                             scalar        stop time
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @return                                                 scalar        the expected number of species
#
################################################################################
globalBiDe.equations.nTaxa.expected.fastApprox <- function(s,t,rate,surv,samplingProbability,MRCA=FALSE) {

  
  p_s <- globalBiDe.equations.pSurvival.fastApprox(rate,surv,samplingProbability,s,t,t,log=FALSE)
  r   <- rate(t) - rate(s) - log(samplingProbability)
  e   <- p_s * exp(r)
  e[e > 1] <- 1

  if ( MRCA == FALSE ) {
    n   <- 1.0 / e
  } else {
    n   <- 2.0 / e
  }

  return (n)
}


################################################################################
# 
# @brief  Calculate the probability of exactly 1 lineage surviving until time T
#         if we started with 1 lineage at time t.
#
# @see    Equation (4) in Hoehna, S.: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    t                                             scalar        time
# @param    T                                             scalar        present time
# @param    log                                           bool          if in log-scale
# @return                                                 scalar        ln-probability
#
################################################################################
globalBiDe.equations.p1.constant <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,T,log=FALSE) {
  
  # compute the survival probability
  a <- globalBiDe.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,T,T,log=TRUE)

  # compute the rate
  rate <- (mu - lambda)*(T-t)
  # add mass-extinction
  for (j in seq_len(length(massExtinctionTimes)) ) {
    rate <- rate - ifelse( t < massExtinctionTimes[j] & T >= massExtinctionTimes[j], log(massExtinctionSurvivalProbabilities[j]), 0 )
  }

  # add sampling
  rate <- rate - log(samplingProbability)
  
  p <- 2*a + rate

  if ( log == FALSE ) {
    p <- exp( p )
  }

  return (p)
}


################################################################################
# 
# @brief  Calculate the probability of exactly 1 lineage surviving until time T
#         if we started with 1 lineage at time t.
#
# @see    Equation (4) in Hoehna, S.: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    r                                             function      rate integral function
# @param    s                                             function      survival integral function
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    t                                             scalar        time
# @param    T                                             scalar        present time
# @param    log                                           bool          if in log-scale
# @return                                                 scalar        ln-probability
#
################################################################################
globalBiDe.equations.p1.fastApprox <- function(r,s,samplingProbability,t,T,log=FALSE) {
  a <- globalBiDe.equations.pSurvival.fastApprox(r,s,samplingProbability,t,T,T,log=TRUE)
  b <- r(T) - r(t) - log(samplingProbability)
  p <- 2*a + b

  if ( log == FALSE ) {
    p <- exp( p )
  }

  return (p)
}


################################################################################
# 
# @brief  Calculate the probability of no speciation event on the reconstructed
#         process between time t and t'.
#
# @date Last modified: 2013-02-06
# @author Sebastian Hoehna
# @version 1.2
# @since 2013-02-06, version 1.2
#
# @param    lambda        scalar        speciation rate
# @param    mu            scalar        extinction rate
# @param    t             scalar        time
# @param    T             scalar        present time
# @param    log           bool          if in log-scale
# @return                 scalar        ln-probability
#
################################################################################
globalBiDe.equations.pWaiting.constant <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,t_prime,T,n,log=FALSE) {

  u <- 1.0 - globalBiDe.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,t_prime,T,log=FALSE) * exp( (mu-lambda)*(t_prime-t) - log(samplingProbability) ) 
  tmp <- u * globalBiDe.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,T,T,log=FALSE) / globalBiDe.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,t_prime,T,log=FALSE)

  # some modification so that we can compute the log
  tmp[tmp > 1.0] <- 1.0
  p <- log(1.0 - tmp ) * n
  
  if (log == FALSE) {
    p <- exp(p)
  }

  return(p)

}


################################################################################
# 
# @brief  Calculate the probability of no speciation event on the reconstructed
#         process between time t and t'.
#
# @date Last modified: 2013-02-06
# @author Sebastian Hoehna
# @version 1.2
# @since 2013-02-06, version 1.2
#
# @param    lambda        scalar        speciation rate
# @param    mu            scalar        extinction rate
# @param    t             scalar        time
# @param    T             scalar        present time
# @param    log           bool          if in log-scale
# @return                 scalar        ln-probability
#
################################################################################
globalBiDe.equations.pWaiting.ageSampling <- function(lambda,r,s,samplingProbability,t,t_prime,T,n,log=FALSE) {

  pdf <- function(x) lambda(x)*globalBiDe.equations.pSurvival.fastApprox(r, s, samplingProbability(x), x, T, T, log=FALSE)
  obj.pdf <- function(t, state, pars) list(pdf(t))
  ## This step is slow for large n - perhaps 1/2s for 1000 points
  n2 <- 11
  times <- seq(t, t_prime, length=n2)
  zz <- ode(0, times, obj.pdf, method = "lsoda", tcrit=t_prime)[,2]

  p <- -n*zz[n2]
  
#  tryCatch( {zz <- integrate(pdf,lower=t,upper=t_prime)$value
#    p <- -n*zz},warning = function(w) {}, error = function(e) {
#      cat("Caught an error in integrate:\n")
#      cat("t =",t,"\n")
#      cat("t_prime =",t_prime,"\n")
#      cat("pdf(t) =",pdf(t),"\n")
#      cat("pdf(t_prime) =",pdf(t_prime),"\n")
#      zz <- integrate(pdf,lower=t,upper=t_prime)$value
#      stop("error in integrate")})

  if (log == FALSE) {
    p <- exp(p)
  }

  return(p)

}




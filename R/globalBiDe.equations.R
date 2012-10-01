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
# @brief Computation of the diversification rate integral for a given time.
#
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    lambda        function      speciation rate function
# @param    mu            function      extinction rate function
# @param    t             scalar        starting time
# @param    tau           scalar        end time
# @param    T             scalar        present time
# @return                 scalar        integral of diversification rates in [t,tau]
#
################################################################################

globalBiDe.equations.rate <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,tau,T) {
 
  r <- tess.integrate(mu,lower=t,upper=tau) - tess.integrate(lambda,lower=t,upper=tau)

  # add mass-extinction rates
  if ( length(massExtinctionTimes) > 0 ) {
    for (i in 1:length(massExtinctionTimes)) {
      if ( massExtinctionTimes[i] > t && massExtinctionTimes[i] <= tau ) {
        r <- r - log( massExtinctionSurvivalProbabilities[i] )
#        r <- r - massExtinctionSurvivalProbabilities[i] + 1
      }
    }
  }

  # add the sampling probability
  if ( t < T && tau >= T ) {
    r <- r - log( samplingProbability )
#    r <- r - samplingProbability + 1
  }

  return (r)
}


################################################################################
# 
# @brief  Calculate the probability of survival in the interval [t_low,t_high].
#
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    lambda        scalar        speciation rate
# @param    mu            scalar        extinction rate
# @param    t_low         scalar        starting time
# @param    t_high        scalar        end time
# @param    T             scalar        present time
# @param    log           bool          if in log-scale
# @return                 scalar        probability of survival in [t,tau]
#
################################################################################
globalBiDe.equations.pSurvival.constant <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t_low,t_high,T,log=FALSE) {
	
   nom <- 1
   den <- 1

   # compute the rate
   rate <- mu - lambda
   
   # add mass-extinction
   accumulatedMassExtinction <- 0
   if ( length(massExtinctionTimes) > 0 ) {
      for (j in 1:length(massExtinctionTimes) ) {
         if ( t_low < massExtinctionTimes[j] && t_high >= massExtinctionTimes[j] ) {
            accumulatedMassExtinction <- accumulatedMassExtinction + log(massExtinctionSurvivalProbabilities[j])
            den  <- den  - (massExtinctionSurvivalProbabilities[j]-1)*exp( rate*(massExtinctionTimes[j]-t_low) - accumulatedMassExtinction)
         }
      }
   }

   # add sampling
   if ( t_low < T && t_high >= T ) {
      accumulatedMassExtinction <- accumulatedMassExtinction + log(samplingProbability)
      den  <- den  - (samplingProbability-1)*exp( rate*(T-t_low) - accumulatedMassExtinction )
   }

   # Only calculate the integral if the interval has at least a certain range, otherwise we treat the integral as 0.
   if ( (t_high-t_low) > 1E-10 ) {
      den <- 1 + (mu/rate)*(exp(rate*(t_high-t_low) - accumulatedMassExtinction) - 1.0)
   }
    
   res <- nom / den

   if ( log == TRUE ) {
     res <- log( res )
   }

   return (res)
}


################################################################################
# 
# @brief  Calculate the probability of survival in the interval [t_low,t_high].
#
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    lambda        function      speciation rate function
# @param    mu            function      extinction rate function
# @param    t_low         scalar        starting time
# @param    t_high        scalar        end time
# @param    T             scalar        present time
# @param    log           bool          if in log-scale
# @return                 scalar        probability of survival in [t,tau]
#
################################################################################
globalBiDe.equations.pSurvival <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t_low,t_high,T,log=FALSE) {

   integrand <- function(x) {
      area <- mu(x) * exp( globalBiDe.equations.rate(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t_low,x,T) )
      
      return(area)
   }
	
   nom <- 1
   den <- 1

   # Only calculate the integral if the interval has at least a certain range, otherwise we treat the integral as 0.
   if ( (t_high-t_low) > 1E-10 ) {
      den <- 1 + tess.integrate(integrand,lower=t_low,upper=t_high) 
   }

   # add mass-extinction
   if ( length(massExtinctionTimes) > 0 ) {
      for (j in 1:length(massExtinctionTimes) ) {
         if ( t_low < massExtinctionTimes[j] && t_high >= massExtinctionTimes[j] ) {
            den <- den - (massExtinctionSurvivalProbabilities[j]-1)*exp( globalBiDe.equations.rate(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t_low,massExtinctionTimes[j],T) )
         }
      }
   }

    # add sampling
   if ( t_low < T && t_high >= T ) {
      den <- den - (samplingProbability-1)*exp( globalBiDe.equations.rate(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t_low,T,T) )
   }

    
   res <- nom / den

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
# @date Last modified: 2012-09-18
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-18, version 1.0
#
# @param    lambda        function      speciation rate function
# @param    mu            function      extinction rate function
# @param    t_low         scalar        starting time
# @param    t_high        scalar        end time
# @param    T             scalar        present time
# @param    log           bool          if in log-scale
# @return                 scalar        probability of survival in [t,tau]
#
################################################################################
globalBiDe.equations.pSurvival.fastApprox <- function(t_low,t_high,T,log=FALSE) {
	
	
   nom <- 1
   den <- 1

   # Only calculate the integral if the interval has at least a certain range, otherwise we treat the integral as 0.
   if ( (t_high-t_low) > 1E-10 ) {
      den <- 1 + ( survivalIntegral(t_high) - survivalIntegral(t_low) ) / exp(rateIntegral(t_low))
   }
    
   res <- nom / den

   if ( log == TRUE ) {
     res <- log( res )
   }

   return (res)
}




################################################################################
# 
# @brief  Calculate the probability of at least i lineages alive at time t and
#         exactly i lineages survive until T.
#
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    lambda        function      speciation rate function
# @param    mu            function      extinction rate function
# @param    i             scalar        number of lineages
# @param    t             scalar        time
# @param    T             scalar        present time
# @param    log           bool          if in log-scale
# @return                 scalar        probability
#
################################################################################
globalBiDe.equations.pReconstructed <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,i,t,T,log=FALSE) {
  prob <- 0
  if (i == 0) {
    prob <- 0
  } else {
    z <- (1 - globalBiDe.equations.pSurvival(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,0,t,T,log=FALSE)*exp(globalBiDe.equations.rate(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,0,t,T))) * globalBiDe.equations.pSurvival(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,0,T,T,log=FALSE)/globalBiDe.equations.pSurvival(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,0,t,T,log=FALSE)
    prob <- (1-z)*(z)^(i-1)
  }

   if ( log == TRUE ) {
     prob <- log( prob )
   }
  
  return (prob)
}


################################################################################
# 
# @brief  Calculate the probability of n lineage existing at time t
#         if we started with 1 lineage at time s.
#
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    lambda        function      speciation rate function
# @param    mu            function      extinction rate function
# @param    i             scalar        number of lineages
# @param    s             scalar        start time
# @param    t             scalar        stop time
# @param    log           bool          if in log-scale
# @return                 scalar        log-probability
#
################################################################################
globalBiDe.equations.pN <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,i,s,t,log=FALSE) {
  p <- 0
  if (i == 1) {
    p   <- globalBiDe.equations.pSurvival(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,s,t,t,log=TRUE) + globalBiDe.equations.rate(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,s,t,t)
  }
  else {
    p_s <- globalBiDe.equations.pSurvival(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,s,t,t,log=FALSE)
    r   <- globalBiDe.equations.rate(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,s,t,t)
    p   <- log(p_s) + r + log( 1 - p_s * exp(r)) * (i-1)
  }

   if ( log == FALSE ) {
     p <- exp( p )
   }

  return (p)
}


################################################################################
# 
# @brief  Calculate the probability of n lineage existing at time t
#         if we started with 1 lineage at time s. This is the fast
#         approximation.
#
# @date Last modified: 2012-09-18
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-18, version 1.0
#
# @param    lambda        function      speciation rate function
# @param    mu            function      extinction rate function
# @param    i             scalar        number of lineages
# @param    s             scalar        start time
# @param    t             scalar        stop time
# @param    log           bool          if in log-scale
# @return                 scalar        log-probability
#
################################################################################
globalBiDe.equations.pN.fastApprox <- function(i,s,t,log=FALSE) {
  p <- 0
  if (i == 1) {
    p   <- globalBiDe.equations.pSurvival.fastApprox(s,t,t,log=TRUE) + rateIntegral(t) - rateIntegral(s)
  }
  else {
    p_s <- globalBiDe.equations.pSurvival.fastApprox(s,t,t,log=FALSE)
    r   <- rateIntegral(t) - rateIntegral(s)
    e   <- min(1, exp(r))
    p   <- log(p_s) + r + log( 1 - p_s * e ) * (i-1)
  }

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
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    lambda        scalar        speciation rate
# @param    mu            scalar        extinction rate
# @param    t             scalar        time
# @param    T             scalar        present time
# @param    log           bool          if in log-scale
# @return                 scalar        ln-probability
#
################################################################################
globalBiDe.equations.p1.constant <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,T,log=FALSE) {
  a <- globalBiDe.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,T,T,log=TRUE)

  # compute the rate
  rate <- (mu - lambda)*(T-t)
  # add mass-extinction
  if ( length(massExtinctionTimes) > 0 ) {
    for (j in 1:length(massExtinctionTimes) ) {
      if ( t < massExtinctionTimes[j] && T >= massExtinctionTimes[j] ) {
        rate <- rate - log(massExtinctionSurvivalProbabilities[j])
      }
    }
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
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    lambda        function      speciation rate function
# @param    mu            function      extinction rate function
# @param    t             scalar        time
# @param    T             scalar        present time
# @param    log           bool          if in log-scale
# @return                 scalar        ln-probability
#
################################################################################
globalBiDe.equations.p1 <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,T,log=FALSE) {
  a <- globalBiDe.equations.pSurvival(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,T,T,log=TRUE)
  b <- globalBiDe.equations.rate(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,T,T)
  p <- 2*a + b

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
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    lambda        function      speciation rate function
# @param    mu            function      extinction rate function
# @param    t             scalar        time
# @param    T             scalar        present time
# @param    log           bool          if in log-scale
# @return                 scalar        ln-probability
#
################################################################################
globalBiDe.equations.p1.fastApprox <- function(t,T,log=FALSE) {
  a <- globalBiDe.equations.pSurvival.fastApprox(t,T,T,log=TRUE)
  b <- rateIntegral(T) - rateIntegral(t)
  p <- 2*a + b

  if ( log == FALSE ) {
    p <- exp( p )
  }

  return (p)
}

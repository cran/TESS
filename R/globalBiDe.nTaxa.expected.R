################################################################################
#
# globalBiDe.nTaxa.expected.R
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
# @brief Computation of the expected number of taxa after a given time.
#
# @date Last modified: 2013-02-01
# @author Sebastian Hoehna
# @version 1.1
# @since 2013-02-01, version 1.1.1
#
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of random sampling at present
# @param    MRCA                                          boolean       does the tree start at the mrca?

# @return                                                 scalar        probability of the speciation times
#
################################################################################

globalBiDe.nTaxa.expected <- function(t_low,t_high,lambda,mu,massExtinctionTimes=c(),massExtinctionSurvivalProbabilities=c(),samplingProbability=1.0,MRCA=TRUE) {

  if ( length(massExtinctionTimes) != length(massExtinctionSurvivalProbabilities) ) {
    stop("Number of mass-extinction times needs to equals the number of mass-extinction survival probabilities!")
  }

  if ( (!inherits(lambda, "numeric") && !inherits(lambda, "function")) || (!inherits(mu, "numeric") && !inherits(mu, "function"))) {
    stop("Unexpected parameter types for lambda and mu!")
  }
  
  # test if we got constant values for the speciation and extinction rates
  if ( class(lambda) == "numeric" && class(mu) == "numeric" ) {
    # call computations for constant rates (much faster)
    p <- globalBiDe.equations.nTaxa.expected.constant(t_low,t_high,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA)
    return (p)
  } else {
    
    # convert the speciation rate into a function if necessary
    if ( inherits(lambda, "numeric") ) {
      speciation <- function (x) rep(lambda,length(x))
    } else {
      speciation <- lambda
    }
    # convert the extinction rate into a function if necessary
    if ( inherits(mu, "numeric") ) {
      extinction <- function (x) rep(mu,length(x))
    } else {
      extinction <- mu
    }
    
    approxFuncs <- tess.prepare.pdf(speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,t_high,c())
    
    n <- globalBiDe.equations.nTaxa.expected.fastApprox(t_low,t_high,approxFuncs$r,approxFuncs$s,samplingProbability,MRCA)
    return (n)
  }

}



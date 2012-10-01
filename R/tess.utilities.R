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


# Global Variables
.TessEnv <- new.env()
assign("N_DISCRETIZATION_NBLOCKS",  10000, envir = .TessEnv)
assign("rateIntegralValues",  c(), envir = .TessEnv)
assign("survivalIntegralValues",  c(), envir = .TessEnv)
assign("maxTime",  1, envir = .TessEnv)



################################################################################
# 
# @brief  This function returns the precomputed rate integral for the
#         interval [0,x]
#
# @date Last modified: 2012-09-19
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-19, version 1.0
#
# @param    x             scalar        time
# @return                 scalar        integrate(rate,0,x)
#
################################################################################
rateIntegral <- function(x) {
  # get the global variables
  N <- get("N_DISCRETIZATION_NBLOCKS",envir=.TessEnv)
  T <- get("maxTime",envir=.TessEnv)
  values <- get("rateIntegralValues",envir=.TessEnv)

  if ( x == 0 ) return ( 0.0 )
  
  index <- ceiling(x*N/T)
  y <- 0
  if ( index == x*N/T) {
    y <- values[index + 1]
  } else {
    y <- values[index] + (values[index+1] - values[index])*( x*N/T - (index-1))
  }

  return (y)
}



################################################################################
# 
# @brief  This function returns the precomputed survival probability integral 
#         for the interval [0,x]
#
# @date Last modified: 2012-09-19
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-19, version 1.0
#
# @param    x             scalar        time
# @return                 scalar        integrate(mu(t)*exp(rate(t,x)),0,x)
#
################################################################################
survivalIntegral <- function(x) {
  # get the global variables
  N <- get("N_DISCRETIZATION_NBLOCKS",envir=.TessEnv)
  T <- get("maxTime",envir=.TessEnv)
  values <- get("survivalIntegralValues",envir=.TessEnv)

  if ( x == 0 ) return ( 0.0 )
  
  index <- ceiling(x*N/T)
  y <- 0
  if ( index == x*N/T) {
    y <- values[index + 1]
  } else {
    y <- values[index] + (values[index+1] - values[index])*(x*N/T - (index-1))
  }
  
  return (y)
}

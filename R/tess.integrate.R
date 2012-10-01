################################################################################
#
# tess.integrate.R
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
# @brief Compute the integral in the inveral [a,b]
#
# @date Last modified: 2012-09-11
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    func          function      function to integrate
# @param    lower         scalar        lower end of interval
# @param    upper         scalar        upper end of interval
# @param    tolerance     scalar        precision tolerance
# @return                 scalar        integral of function in [lower,upper]
#
################################################################################
tess.integrate <- function(func,lower,upper,tolerance=1E-06) {

#  tmp <- Vectorize(func)
#  return (integrate(tmp,lower,upper)$value)
  
  # check if the length of the interval is larger than 0
  if (lower == upper) {
    return (0.0)
  }
    
  # compute the 4 interval borders
  midpoint         <- (lower+upper) / 2
    
  # compute the function values at each interval endpoint
  f_a              <- func( lower )
  f_m              <- func( midpoint )
  f_b              <- func( upper )
    
  if ( !is.finite(f_a) || !is.finite(f_m) || !is.finite(f_b) ) {
    stop("Infinite function value in integration.")
  }
    
  integral_left    <- (midpoint-lower) / 2.0 * ( f_a + f_m )
  integral_right   <- (upper-midpoint) / 2.0 * ( f_m + f_b )
    
  integral         <- tess.integrate.recursive( func, lower, midpoint, f_a, f_m, integral_left, tolerance) + tess.integrate.recursive( func, midpoint, upper, f_m, f_b, integral_right, tolerance)
    
  return (integral)
    
}


tess.integrate.recursive <- function(func,lower,upper,f_a,f_b,rough_integral,tolerance=1E-06) {
  # check if the length of the interval is larger than 0
  if (lower == upper) {
    return (0.0)
  }
    
  # compute the 4 interval borders
  midpoint         <- (lower+upper) / 2
    
  # compute the function values at each interval endpoint
  f_m              <- func( midpoint )
     
  if ( !is.finite(f_m) ) {
    stop("Infinite function value in integration.")
  }

  integral_left    <- (midpoint-lower) / 2.0 * ( f_a + f_m )
  integral_right   <- (upper-midpoint) / 2.0 * ( f_m + f_b )
    
  integral         <- integral_left + integral_right;
  if ( abs( rough_integral - integral ) > tolerance && abs( 1.0 - rough_integral / integral ) > tolerance && abs(lower-upper) > tolerance ) {
    integral       <- tess.integrate.recursive( func, lower, midpoint, f_a, f_m, integral_left, tolerance) + tess.integrate.recursive( func, midpoint, upper, f_m, f_b, integral_right, tolerance)
  }
        
  return (integral)
    
}



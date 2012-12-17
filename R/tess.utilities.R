################################################################################
#
# tess.utilities.R
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



tess.create.phylo <- function(times,root=FALSE,tip.label=NULL) {
  n <- as.integer(length(times))+1
  if ( root ) {
    n <- n-1
  }
  nbr <- 2*n - 2 

  # create the data types for edges and edge-lengths
  edge <- matrix(NA, nbr, 2)
  edge.length <- numeric(nbr)
  
  h <- numeric(2*n - 1) # initialized with 0's
  pool <- 1:n
  # VERY VERY IMPORTANT: the root MUST have index n+1 !!!
  nextnode <- 2L*n - 1L
  if ( n > 1) {
    for (i in 1:(n - 1)) {
      # sample two nodes that have no parent yet
      y <- sample(pool, size = 2)
      # compute the edge indices (we just order the edges from 1 to 2n-2)
      ind <- (i - 1)*2 + 1:2
      # set the source node of the new edges (i.e. the new internal node)
      edge[ind, 1] <- nextnode
      # set the destination of the new edges (i.e. the two sampled nodes)
      edge[ind, 2] <- y
      # compute the edge length from the difference between the node heights (child <-> parent)
      edge.length[ind] <- times[i] - h[y]
      # store the node height of this new internal node
      # we cannot use times because then we would get into trouble with the indices and we would need to check for tip nodes ...
      h[nextnode] <- times[i]
      # reset the pool of available nodes to merge
      pool <- c(pool[! pool %in% y], nextnode)
      # increase the node index counter
      nextnode <- nextnode - 1L
    }
  }

  phy <- list(edge = edge, edge.length = edge.length)
  if (is.null(tip.label))
    tip.label <- paste("t", 1:n, sep = "")
  phy$tip.label <- sample(tip.label)
  phy$Nnode <- n - 1L

  if ( root ) {
    phy$root.edge <- times[n] - times[n-1]
    phy$root <- times[n] - times[n-1]
  }

  class(phy) <- "phylo"
  
  phy <- reorder(phy)
  ## to avoid crossings when converting with as.hclust:
  phy$edge[phy$edge[, 2] <= n, 2] <- 1:n

  phy
}




tess.prepare.pdf <- function(lambda, mu, massExtinctionTimes,
                         massExtinctionSurvivalProbabilities,
                         samplingProbability,
                         age, t.crit.f) {
  
  ## Constants for now, but this is the resolution of reporting
  ## (different to calculation) for the two stages of integration.
  n1 <- 101

  
  if ( length(massExtinctionTimes) == 0 ) {
    t.crit <- sort(unique(c(0, t.crit.f,age)))
    t.crit <- t.crit[t.crit <= age]
    t.crit.i <- c()
    t.crit.p <- c(rep(NA,length(t.crit)-1),1.0)
  } else {
    t.crit <- sort(unique(c(0, t.crit.f, massExtinctionTimes)))
    t.crit <- t.crit[t.crit < age]
    t.crit.i <- match(t.crit, massExtinctionTimes)
    t.crit.p <- massExtinctionSurvivalProbabilities[t.crit.i]
    t.crit <- c(t.crit, age)
    t.crit.p <- c(t.crit.p, 1.0)
  }

  times <- seq(0, age, length=n1)
  rs <- tess.ode.piecewise(lambda, mu, times, t.crit, t.crit.p)
    
  return (rs)
}

## Utility functions.
first <- function(x) x[[1]]
last <- function(x) x[[length(x)]]


## This carries out integration for the time-varying
## speciation/extinction model with functions lambda and mu,
## outputting at the times in the vector 'times'.  The vectors
## 't.crit'  and 't.crit.p' contain break points or mass-extinction events.
tess.ode.piecewise <- function(lambda, mu, times, t.crit, t.crit.p) {
  obj <- function(t, state, parms) {
    r <- mu(t) - lambda(t)
    s <- mu(t) * exp(state[[1]])
    list(c(r, s))
  }


  tmax <- times[length(times)]
  n <- length(t.crit) - 1L
  t.split <- r.split <- s.split <- vector("list", n)
  bin <- findInterval(times, t.crit, TRUE)
  y <- c(0, 0)

  for ( i in seq_len(n) ) {
    j <- i + 1L
    t.split[[i]] <- ti <- c(t.crit[[i]], times[bin == i], t.crit[[j]])
    yi <- lsoda(y, ti, obj, tcrit=t.crit[[j]])

    r.split[[i]] <- yi[,2]
    s.split[[i]] <- yi[,3]
    y <- yi[nrow(yi),-1]

    if ( !is.na(t.crit.p[[j]]) ) {
      y[1] <- y[1] - log(t.crit.p[[j]])
      y[2] <- y[2] - (t.crit.p[[j]]-1) * exp(y[1])
    }
  }

  ## There should be a more general way of dealing with single values
  ## within an interval -- this will come up elsewhere in the general
  ## case.  For now, fudge a second point.  You'll need three points
  ## here if you want to replace approxfun with splinefun.
  t.split <- c(t.split, list(c(tmax, tmax + 1e-8)))
  r.split <- c(r.split, list(rep.int(y[[1]], 2)))
  s.split <- c(s.split, list(rep.int(last(s.split[[n]]), 2)))

  list(r=make.piecewise(t.split, r.split),
       s=make.piecewise(t.split, s.split))
}


## Constructs a piecewise function, allowing jumps in a function.  Let
## xx and yy be lists, with each element containing a vector of x or y
## values respectively.  That is x[[i]] and y[[i]] will each contain n
## points, defining a smooth line.  This uses R's linear interpolation
## function to construct functions over each pair of x/y coordinates,
## returning a vectorised function.
make.piecewise <- function(xx, yy) {
  if ( length(xx) == 1 )
    return(approxfun(xx, yy))
  ff <- mapply(approxfun, xx, yy)
  t.crit <- c(sapply(xx, first), last(last(xx)))
  function(x) {
    i <- findInterval(x, t.crit, TRUE)
    ret <- numeric(length(x))
    for ( j in unique(i) ) {
      k <- i == j
      ret[k] <- ff[[j]](x[k])
    }
    ret
  }
}



find.max <- function(pdf,m,n) {
  
  xx <- seq(0, m, length=n)
  yy <- pdf(xx)
  i <- which.max(yy)
  sup <- optimize(pdf, range(na.omit(xx[(i-1):(i+1)])), maximum=TRUE)$max

  return (sup)
}


## Sample n random variates from a distribution with a PDF
## proportional to 'f', which has domain 'r', knowing that max(f) over
## this domain is always less than or 'sup'
rejection.sample.simple <- function(n, f, r, sup) {
  ok <- numeric(0)
  repeat {
    u <- runif(n, r[1], r[2])
    ok <- c(ok, u[f(u) / sup > runif(n)])
    if ( length(ok) >= n )
      break
  }
  ok[seq_len(n)]
}


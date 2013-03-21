\name{mammalia}
\alias{mammalia}
\title{Dated family level mammalian phylogeny from Meredith et al. (2011): Impacts of the cretaceous terrestrial revolution and kpg extinction on mammal diversification.}
\description{
  This phylogeny describes the species relationship and divergence times of the class Mammalia with 1-3 species included per family, published in Meredith et al. (2011).
}
\usage{
data(mammalia)
}
\format{
  The phylogeny is stored as an object of class \code{"phylo"} which
  structure is described in the help page of the function
  \code{\link{read.tree}} of the package \code{\link{ape}}.
}
\source{
   Meredith, R. et al. (2011): Impacts of the cretaceous terrestrial revolution and kpg extinction on mammal diversification. Science, 334:521-524
}

\examples{
data(mammalia)
op <- par()
par(cex = 0.3)
plot(mammalia)
par(op)
}
\keyword{datasets}

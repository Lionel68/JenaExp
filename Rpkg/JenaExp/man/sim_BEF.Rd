\name{sim_BEF}
\alias{sim_BEF}
\title{Simulate a Biodiversity-Ecosystem Function design}
\usage{ sim_BEF(n_species=10,n_rep=5,gradient=1:10) }
\arguments{
\item{n_species}{The number of species used in the experiment, default is 10 species, for now the function only work for n_species<676.}
\item{n_rep}{The number of replicates per species richness levels, default is 5, can be an integer or a numeric vector of length equal to the gradient argument.}
\item{gradient}{The levels of species richness.}
}
\value{
A site x species community matrix filled with 0/1 for species absence/presence in particular sites.
}
\description{}
\author{Lionel Hertzog with code borrowed from Cameron Wagg}
\examples{
#using default values
sim_BEF()
#with more species and a longer gradient
sim_BEF(n_species=20,n_rep=4,gradient=c(1,2,4,8,16))
#with varying number of replicates per species richness level
sim_BEF(n_species=20,n_rep=c(20,rep(5,4),4),gradient=c(1,2,4,8,16,20)
}
\keyword{simulation,BEF,design}

\name{calcFD}
\alias{calcFD}
\title{Compute FD community level Functional Diversity index}
\usage{ calcFD(species_matrix,trait_table,weighted=FALSE,distance = "euclidean",standardize = "None",calc.FRic=FALSE,calc.FDiv=FALSE,calc.FEve=FALSE,calc.FDis=FALSE,corr="sqrt") }
\arguments{
\item{species_matrix}{ a site x species community matrix either as 0/1 or as abundance. }
\item{trait_table}{a species x traits data frame, traits can be numeric, factor, or character.}
\item{weighted}{should the traits index be weighted by the abundance of the species, default is FALSE.}
\item{distance}{the distance metric to use to compute the distance between the species in the MDS space, default is "eucledian", can be any distance metric defined in the function vegdist in the package vegan.}
\item{standardize}{should the trait be standardize, default is "None" can be any character string accepted by the function decostand in the package vegan.}
\item{calc.FXx}{should the various index from the dbFD function in the FD package be computed, default is FALSE.}
\item{corr}{the correction to apply to the initial species distance matrix if it is not eucledian, default is "sqrt" but can be any character string accepted by the dbFD function, like "none".}
}
\value{
a data frame with the FD value for the different sites in the community
}
\description{}
\author{Cameron Wagg and Lionel Hertzog}
\examples{
data(BCI)
traits<- data.frame(T1=rnorm(dim(BCI)[2],10,2),T2=runif(dim(BCI)[2],0,10))
rownames(traits)<-colnames(BCI)
fd_exp <- calcFD(BCI,traits,weighted=TRUE)
}
\keyword{FD}

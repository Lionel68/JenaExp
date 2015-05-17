\name{helper_calcFD}
\alias{helper_calcFD}
\title{An helper function}
\usage{}
\arguments{}
\value{}
\description{}
\author{Lionel Hertzog}
\examples{
data(BCI)
traits<- data.frame(T1=rnorm(dim(BCI)[2],10,2),T2=runif(dim(BCI)[2],0,10))
rownames(traits)<-colnames(BCI)
fd_exp <- calcFD(BCI,traits,weighted=TRUE)
}
\keyword{FD}

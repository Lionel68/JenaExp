\name{price}
\alias{price}
\title{Decompose BEF effect using the extended Price equation}
\usage{ price(data,control,avg=TRUE) }
\arguments{
\item{data}{A data frame of sites x species contribution to an ecosystem function at comparison sites.}
\item{control}{A data frame of sites x species contribution to an ecosystem function at control sites. The number of row can be bigger than 1.}
\item{avg}{If the number of control sites is bigger than 1 should the results be averaged for each comparison sites over the control sites.}
}
\details{
The column order in data and control must be matched, this means that species X must be in the same column in both data and control (ie in column 1). 
For now results are automatically averaged if there is more than one control sites.
If there are no species in common between one comparison site and the baseline site then the function return a NA value.
The equations are based on Fox and Kerr 2012 Oikos paper.
}
\value{
A data frame where every sites is one row with the 5 terms in the price equations: SRL, SRG, SCL, SCG, CDE. See the reference to understand the meaning of these different columns.
}
\references{
Fox, Jeremy W., and Benjamin Kerr. "Analyzing the effects of species gain and loss on ecosystem function using the extended Price equation partition." Oikos 121.2 (2012): 290-298.
}
\author{Lionel Hertzog}
\examples{
#ten species in ten comparison sites
data<-matrix(rnorm(100,10),ncol=10,nrow=10,dimnames=list(paste0("Site",1:10),paste0("Sp",1:10)))
#two control sites
control<-matrix(rnorm(20,10),ncol=10,dimnames=list(paste0("Control",1:2),paste0("Sp",1:10)))
#the function call
price(data,control)
}
\keyword{Price equation,BEF}

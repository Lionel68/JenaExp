\name{biomass}
\alias{biomass}
\title{Biomass data from the Jena experiment between 2003 and 2008}
\usage{data(biomass)}
\arguments{}
\value{A site x species community matrix with plot code as rownames and species abbreviations as column names, format is: Gen.spe}
\description{Biomass data collected twice per year on the Jena field site. Only the biomass of the sown species was kept, the values were averaged per species over the month and the year.}
\author{Lionel Hertzog}
\source{http://www.esajournals.org/doi/abs/10.1890/09-0863.1}
\examples{
data(biomass)
summary(as.numeric(biomass))
}
\keyword{biomass}

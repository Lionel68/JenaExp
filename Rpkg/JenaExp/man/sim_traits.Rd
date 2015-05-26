\name{sim_traits}
\alias{sim_traits}
\title{Simulate traits data}
\usage{ sim_traits(n_traits=5,n_species=10,distr="rnorm",sp_names=NULL,params=NULL) }
\arguments{
\item{n_species}{The number of species used in the experiment, default is 10 species.}
\item{n_traits}{The number of traits simulated default is 5, must be 1 or bigger.}
\item{distr}{Either a character string specifying a distribution function, or a vector of character strings of length equal to n_traits specifying the distribution for the j's traits. See Details for a list of accepted functions.}
\item{sp_names}{A vector of character strings of length equal to n_species giving the row names in the outputed table, default is set to NULL.}
\item{params}{A list of parameters for the various distribution function specified in the distr argument, must match in length and order the distr vector.}
}
\details{
The following functions are currently accepted by the distr argument: "rnorm", "rlnorm", "runif", "rcauchy" for the continuous distributions and "factor" for factorial traits.
When no parameters are specified the function use the default values for the continuous distribution, see ?rnorm for example. For factorial traits 4 levels with equal probabilities are taken per default, the function "base::sample" is used to draw the different levels.
The parameters of the various function can also be given as argument to the function, see the examples below. The given parameter are then interpreted in the same order are they are used in the functions. For example in rnorm the first parameter will be interpreted as the mean and the second one as the standard deviation.
For factorial traits the number of levels and the probabilities associated with each level can be specified. The params argument is then a bit more complex, it should look like this: params = list(list(n_levels,c(p_1,p_2 ... p_i))), for i = 1 ... n_levels, examples are given below. Also note that as the function "base::sample" is used the probabilities do not need to sum to one. 
If no sp_names argument is given the row names in the returned data frame are like this: Sp_1, Sp_2 ... Sp_j, for j = 1 ... n_species.
}

\value{
A species x traits data frame, factorial traits are returned as factor while continuous traits are returned as numeric variables.
}

\author{Lionel Hertzog}
\examples{
#using default values
sim_traits()
#with more species and specifying the distribution of the traits
sim_traits(n_species=20,n_traits=4,distr=c("rnorm","runif","factor","rlnorm"))
#specifying parameter values for the distribution functions
param1 <- list(c(10,5),c(-2,2),c(3,1))
sim_traits(n_species=20,n_traits=3,distr=c("rnorm","runif","rlnorm"),params=param1)
#with factorial traits variables the params must be specified differently
param2 <- list(c(10,5),list(5,c(0.1,0.4,0.3,0.1,0.1)))
jena<-sim_traits(n_traits=2,distr=c("rnorm","factor"),params=param2)
}
\keyword{simulation,BEF,traits}

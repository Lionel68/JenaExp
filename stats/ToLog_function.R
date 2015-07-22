#########################
#Function files investiguating log-transformation vs GLM
#author Lionel Hertzog, date 21.07.2015
#code heavily based on the supplementary material of
#Ives, A. R. (2015), For testing the significance of regression coefficients, go ahead and log-transform count data. Methods in Ecology and Evolution, 6: 828–835. doi: 10.1111/2041-210X.12386




###################################################################
### Univariate, negative binomial functions (Figs. 1 and 2A-C)
EstPvalue = function(y, x, add = 0,transformation) {
  if (transformation == "negbin") {
    try({
      z<-glm.nb(y ~ x, control = glm.control(maxit = 30))
      conv<-z$converged
      if(conv==1){
        val <- summary(z)[[11]][8]
        b1.est <- coef(z)[2]
      }
      else{
        val <- NA
        b1.est <- NA
      }
     })    
  }
  if (transformation == "qpois") {
    z <- glm(y ~ x, family = quasipoisson, control = list(maxit = 1000))
    conv<-z$converged
    if(conv==1){
      val <- summary(z)[[12]][8]
      b1.est <- coef(z)[2]
    }
    else{
      val<-NA
      b1.est<-NA
    }
  }
  if (transformation == "glmm") {
    try({
      z <- glmer(y ~ x + (1 | as.factor(1:length(x))), family = "poisson")
      conv <- (attributes(z)$optinfo$conv$opt == 0)
      if (conv == 1) {
        val <- summary(z)[[10]][2, 4]
        b1.est <- fixef(z)[2]
      } else {
        val <- NA
        b1.est <- NA
      }
    })    
  }
  if (transformation == "log") {
    z <- lm(log(y + add) ~ x)
    val <- summary(z)[[4]][8]
    b1.est <- coef(z)[2]
  }
  if(exists("val")){
    ret <- c(val, b1.est)
    return(ret)
  }
  else{
    return(c(NA,NA))
  }
}


# Function to estimate p-values
EstStats.pvalue <- function(Data, transformation, b1,Add = 0, alpha = 0.05) {
  pvalues <- apply(Data$y,2,function(y) EstPvalue(y = y,x = Data$x,transformation = transformation,add = Add))
  ret <- c(mean(pvalues[1, ] < alpha,na.rm=TRUE), mean(pvalues[2, ]-b1,na.rm=TRUE))  
  return(ret)
}

# Function to fit all models
GetAnalyses = function(Data, alpha = 0.05,b1,GLMM=FALSE) {
  
  NB = EstStats.pvalue(Data, transformation = "negbin",b1)
  QPois = EstStats.pvalue(Data, transformation = "qpois",b1)
  DataHalf <- Data
  DataHalf$y[DataHalf$y == 0] <- 0.5
  LogPlusHalf = EstStats.pvalue(DataHalf, transformation = "log",b1)
  Log1 = EstStats.pvalue(Data, transformation = "log",b1, Add = 1)
  Log0001 = EstStats.pvalue(Data, transformation = "log", b1,Add = 1e-04)
  
  
  if(GLMM){
    GLMM<-EstStats.pvalue(Data,transformation = "glmm",b1)
    d<-data.frame(Model=c("NB","QuasiP","GLMM","LogHalf","Log1","Log0001"),Reject=c(NB[1],QPois[1],GLMM[1],LogPlusHalf[1],Log1[1],Log0001[1]),Bias=c(NB[2],QPois[2],GLMM[2],LogPlusHalf[2],Log1[2],Log0001[2]))
  }
  else{
    d<-data.frame(Model=c("NB","QuasiP","LogHalf","Log1","Log0001"),Reject=c(NB[1],QPois[1],LogPlusHalf[1],Log1[1],Log0001[1]),Bias=c(NB[2],QPois[2],LogPlusHalf[2],Log1[2],Log0001[2]))
    
  }
  
  return(d)
}

# Function to simulate and fit data for the univariate negative binomial model
compute.stats <- function(NRep = 50, b1.range = 0, n.range = 100, dispersion.range = 1, b0.range = log(1), alpha = 0.05,seed=20150721) {
  val<-expand.grid(b1=b1.range,n=n.range,disper=dispersion.range,b0=b0.range)
  tmp<-mapply(function(b1,n,disper,b0){
    set.seed(seed)
    x <- runif(n,-2,2)
    mean.y <- exp(b0 + b1 * x)
    y <- replicate(NRep, rnbinom(n, disper, mu = mean.y))
    Data <- list(x = x, y = y)
    g <- GetAnalyses(Data, alpha,b1)
    g$b1<-b1
    g$n<-n
    g$disper<-disper
    g$b0<-b0
    return(g)
  },b1=val$b1,b0=val$b0,n=val$n,disper=val$disper,SIMPLIFY=FALSE)
  output<-rbind.fill(tmp)
  return(output)
}

#########################################################################
### Univariate, lognormal-Poisson hierarchical model functions

# Function to simulate and fit data for the univariate negative binomial model
compute.statsGLMM <- function(NRep = 50, b1.range = 0, n.range = 100, sd.eps.range = 1, b0.range = log(1), alpha = 0.05,seed=20150723) {
  val<-expand.grid(b1=b1.range,n=n.range,b0=b0.range,sd.eps=sd.eps.range)
  tmp<-mapply(function(b1,n,b0,sd.eps){
    set.seed(seed)
    x <- runif(n,-2,2)
    eps<-rnorm(n,0,sd.eps)
    mean.y <- exp(b0 + b1 * x+eps)
    y <- replicate(NRep, rpois(n, mean.y))
    Data <- list(x = x, y = y)
    g <- GetAnalyses(Data, alpha,b1,GLMM=TRUE)
    g$b1<-b1
    g$n<-n
    g$b0<-b0
    g$sd.eps<-sd.eps
    return(g)
  },b1=val$b1,b0=val$b0,n=val$n,sd.eps=val$sd.eps,SIMPLIFY=FALSE)
  output<-rbind.fill(tmp)
  return(output)
}
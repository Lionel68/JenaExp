###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
# FD calculations used in Wagg et al. manuscript in prep. for analysis of traits and biodiversity effects in the Jena TBE data from 2012
#modified on 15.05 by Lionel Hertzog to include mixed trait data and new FD index computation
# Literature: 
#  Clark, C.M., Flynn, D.F.B., Butterfield, B.J. & Reich, P.B. (2012). Testing the link between functional diversity and ecosystem functioning in a Minnesota grassland experiment. PLoS ONE 7(12): e52821. doi:10.1371/journal.pone.0052821
# Laliberté, E. and P. Legendre (2010) A distance-based framework for measuring functional diversity from multiple traits. Ecology 91:299-305.
# Petchey, O.L. & Gaston, K.J. (2002) Functional diversity (FD), species richness and community composition. Ecology Letters 5: 402-411. 
# Petchey, O.L. & Gaston, K.J. (2006) Functional diversity: back to basics and looking forward. Ecology Letters 9: 741-758.


library(FD)
library(psych)

#load the function
source_github <- function(u) {
  # load package
  require(RCurl)
  
  # read script lines from website
  script <- getURL(u, ssl.verifypeer = FALSE)
  
  # parase lines and evealuate in the global environement
  eval(parse(text = script))
}
source("https://github.com/Lionel68/JenaExp/blob/master/FD/calcFD_LH20150514.R")


###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
## GENERATE EXAMPLE DATA
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### Generate fake biodiversity experiment data
set.seed(13)
sp = LETTERS[1:20]
sp2 = combn(sp, 2)[,sample(dim(combn(sp, 2))[2], 10)] # 2 sp mix
sp3 = combn(sp, 3)[,sample(dim(combn(sp, 3))[2], 10)] # 3 sp mix
sp4 = combn(sp, 4)[,sample(dim(combn(sp, 4))[2], 10)] # 4 sp mix
sp6 = combn(sp, 6)[,sample(dim(combn(sp, 6))[2], 10)] # 6 sp mix
sp9 = combn(sp, 9)[,sample(dim(combn(sp, 9))[2], 10)] # 9 sp mix

sp.comb = sp

for (i in 1:10){
	sp.comb = c(sp.comb, paste(sp2[,i], collapse=""))
}
for (i in 1:10){
	sp.comb = c(sp.comb, paste(sp3[,i], collapse=""))
}
for (i in 1:10){
	sp.comb = c(sp.comb, paste(sp4[,i], collapse=""))
}
for (i in 1:10){
	sp.comb = c(sp.comb, paste(sp6[,i], collapse=""))
}
for (i in 1:10){
	sp.comb = c(sp.comb, paste(sp9[,i], collapse=""))
}

sp.mat = matrix(nrow = length(sp.comb), ncol = length(sp), dimnames=list(sp.comb, sp))
sp.list = strsplit(sp.comb, "")

sr = rep(NA, dim(sp.mat)[1])
for (i in 1:dim(sp.mat)[1]){
	sr[i] = length(sp.list[[i]])
	sp.mat[i, sp.list[[i]]] = sample(20:2000, length(sp.list[[i]]), replace=TRUE)*(1/sample(1:sr[i]))
}



##### Generate trait table for each species
traits = data.frame(Trt1=rnorm(length(sp),10,3),Trt2=runif(length(sp),0,10),Trt3=sample(c("A","B","C","D"),length(sp),replace=TRUE),Trt4=factor(sample(1:4,length(sp),replace=TRUE),ordered = TRUE),row.names=sp)

# check for negative values and change them to positive proportions
traits[traits< 0] = -1*traits[traits< 0] / 10

# remove row names
#rownames(sp.mat) = as.character(1:dim(sp.mat)[1]) (?)

#sp.mat # species matrix with absolute abundances
#sr # species richness levels
#traits # trait table
#sp.comb # species combination

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
# Example analysis using the function:
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 

library(vegan)


# Fake data looks like typical biodiversity experiment?
plot(sr, rowSums(sp.mat, na.rm=TRUE), las=1, ylab="Total biomass", xlab="Species richness")
anova(lm(rowSums(sp.mat, na.rm=TRUE) ~ sr + as.factor(sr)))
summary(lm(rowSums(sp.mat, na.rm=TRUE) ~ sr + I(sr^2)))
cf = summary(lm(rowSums(sp.mat, na.rm=TRUE) ~ sr + I(sr^2)))$coef[,1]
x = seq(0, max(sr), 0.1)
y = cf[1] + cf[2]*x + cf[3]*x^2
lines(x,y)

plot(log(sr), rowSums(sp.mat, na.rm=TRUE), las=1, ylab="Total biomass", xlab="Species richness (Log)")
summary(lm(rowSums(sp.mat, na.rm=TRUE) ~ log(sr)))
abline(lm(rowSums(sp.mat, na.rm=TRUE) ~ log(sr)))
## ## ## ## ## 


Clust = hclust(gowdis(traits), method="complete")
plot(as.dendrogram(Clust), las=1, cex=0.75, type="triangle", xlab="Species", ylab="Height", main="Functional trait dendrogram") # for nice dendrograms see http://gastonsanchez.com/blog/how-to/2012/10/03/Dendrograms.html

# Run function at end of script below before using calc.FD()

###  # using defult settings
test = calc.FD(sp.mat, traits)
#test

# MDS ordination of traits
#plot(test$mds.points, pch=rownames(traits))

###  # using standardized trait values (?decostand, for options)
#test = calc.FD(sp.mat, traits, standardize="standardize") hard to do with mixed trait data, in a later implementation we could target continuous traits and apply the correction only to them
#test

# MDS ordination of traits
#plot(test$mds.points, pch=rownames(traits))


###  # using standardized trait values and weighting by observed proportions
test = calc.FD(sp.mat, traits, weighted = TRUE)
#test



##
# compare the two FD methods
pairs.panels(cbind(test$FD, test$FDmpd, sr), labels=c("FD", "FDmpd", "SR"))


# note correlation can be due to 0 values in monocultures / monocultures should be removed from downstream analyses using these measures
cor.test(test$FD, test$FDmpd) # monocultures included
cor.test(test$FD[which(sr!=1)], test$FDmpd[which(sr!=1)]) # monocultures excluded

##
# compare the two FD methods with species richness
par(mfrow=c(1,2))
plot(sr, test$FD, las=1, ylab="Petchey-Gaston FD", xlab="Species richness")
plot(sr, test$FDmpd, las=1, xlab="Species richness", ylab="FD (mean MDS paiwise distance)")

## 
# relation to total biomass?
par(mfrow=c(1,2))
plot(test$FD, rowSums(sp.mat, na.rm=TRUE), las=1, ylab="Total BM", xlab="FD (Petchey-Gaston)")
plot(test$FDmpd, rowSums(sp.mat, na.rm=TRUE), las=1, ylab="Total BM", xlab="FD (mean MDS paiwise distance)")

# omit monocultures
plot(test$FD[which(sr!=1)], rowSums(sp.mat, na.rm=TRUE)[which(sr!=1)], las=1, ylab="Total BM", xlab="FD (Petchey-Gaston)")
cor.test(test$FD[which(sr!=1)], rowSums(sp.mat, na.rm=TRUE)[which(sr!=1)])
plot(test$FDmpd[which(sr!=1)], rowSums(sp.mat, na.rm=TRUE)[which(sr!=1)], las=1, ylab="Total BM", xlab="FD (mean MDS paiwise distance)")
cor.test(test$FDmpd[which(sr!=1)], rowSums(sp.mat, na.rm=TRUE)[which(sr!=1)])

# compare the two FD methods weighted and unweighted
###  # using defult settings
test = calc.FD(sp.mat, traits)
#test

test.wt = calc.FD(sp.mat, traits, weighted = TRUE)
#test.wt

#pairs.panels(data.frame(rowSums(sp.mat, na.rm=TRUE), test$FD, test.wt$FD, test$FDmpd, test.wt$FDmpd, sr))

#compare with other FD index 1) unweighted
test<-calc.FD(sp.mat,traits,calc.FDis=TRUE,calc.FEve=TRUE,calc.FDiv=TRUE,calc.FRic=TRUE,corr="none")
test$BM<-rowSums(sp.mat,na.rm=TRUE)
pairs.panels(test[,c(3,8,1,2,4:7)])

#compare with other FD index 2) weighted
test.wt<-calc.FD(sp.mat,traits,weighted=TRUE,calc.FDis=TRUE,calc.FEve=TRUE,calc.FDiv=TRUE,calc.FRic=TRUE,corr="none")
test.wt$BM<-rowSums(sp.mat,na.rm=TRUE)
pairs.panels(test.wt[,c(3,8,1,2,4:7)])


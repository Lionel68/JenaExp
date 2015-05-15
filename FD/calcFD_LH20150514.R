###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
# FD calculations used in Wagg et al. manuscript in prep. for analysis of traits and biodiversity effects in the Jena TBE data from 2012
#modified by Lionel Hertzog on 14.05.2015 to include non-numeric trait data
# Literature: 
#  Clark, C.M., Flynn, D.F.B., Butterfield, B.J. & Reich, P.B. (2012). Testing the link between functional diversity and ecosystem functioning in a Minnesota grassland experiment. PLoS ONE 7(12): e52821. doi:10.1371/journal.pone.0052821
# Laliberté, E. and P. Legendre (2010) A distance-based framework for measuring functional diversity from multiple traits. Ecology 91:299-305.
# Petchey, O.L. & Gaston, K.J. (2002) Functional diversity (FD), species richness and community composition. Ecology Letters 5: 402-411. 
# Petchey, O.L. & Gaston, K.J. (2006) Functional diversity: back to basics and looking forward. Ecology Letters 9: 741-758.



###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### Explanation of input options / conditions
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
# package vegan and FD is required for standardizing, clusterig, distances and calculating banch lengths

# species_matrix = species matrix with colums as individual species and rows as sites / plots. detected species can be absolute values or proportions. Non-detected species should be 0's, while species not present / sown at all should be NA's.

# trait_table = table of traits of each species with matrix rows as individual species and row names as species names. I have not tried this with factoral traits, so will only work if trats are numeric.

# weighted = should the distance matricies be weighted by the proportional difference in species abundance? where weighed distances are: dist * (1-abs(pi-pj)) # the absolute difference in species proportions, such that distance weightings are 0 if one species is 100% and distance weightings are 1 if species are equal in proportional abundance.

# distance = the distance metric used to calculate the distance matrix among species in their trait values. Should be one of methods in function vegdist()

# clustering = clustering method for hclust() function

# standardize = standardize method for trait values (standardize columns in trait_table?). Should be either "None" for no standardization or one of the methods in function decostand() in the package vegan

# Note branch lengths can be better manipulated (ie for standarization etc...) . see https://github.com/ibartomeus/fundiv/blob/master/R/Xtree.R

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### The function
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 

calc.FD<-function(species_matrix,trait_table,weighted=FALSE,distance = "euclidean",standardize = "None") {

# # # # # # # Generate trait distance matrix # # # # # # # # #
  #apply standardization if asked for
  if(standardize!="None"){
    trait_table<-decostand(trait_table, method=standardize)
  }
  #if all traits variables are numeric compute the euclidean distance
  if(all(apply(trait_table,2,is.numeric))){
    D <- dist(trait_table)
  }
  #otherwise use the gowdis distance from package FD
  else{
    D <- gowdis(trait_table)
  }
  

# # # # # # # Make species matrix binary if distances are to be unweighted # # # # # # # # #
	if (weighted == FALSE)
		species_matrix[species_matrix>0] <- 1

# set up NMDS distance matrix
	mds <- metaMDS(D, k=2,)
	mds.points <- mds$points # save mds points for output

	# calculate distances between points in the MDS space 
	#(I prefer "euclidean" for ease and options here, but could also be calculated by hand using sqrt((x1-x2)^2 + (y1-y2)^2)  )
	MDSdist = as.matrix(vegdist(mds$points, method=distance)) 	
	
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # Go through each sampled sites and compute FD and FDmpd
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 output<-apply(species_matrix,1,function(x) helper_fun(x))
 return(as.data.frame(t(output)))
}

helper_fun<-function(site){
  #remove non present species
  site<-site[site!=0 & !is.na(site)]
  if(sum(site,na.rm=TRUE)<=1){
    return(c(FD=0,FDmpd=0))
  }
  else{
    ps<-site/sum(site)
    #compute FDmpd
    smat<-1-as.matrix(dist(ps,method="manhattan"))
    d_mds<-MDSdist[names(ps),names(ps)]*smat 
    FDmpd<-mean(as.numeric(d_mds[lower.tri(d_mds)]))
    #compute FD
    Dsub<-dist(trait_table[rownames(trait_table)%in%names(ps),])
    clust<-hclust(Dsub,method="complete")
    FD<- treedive(t(as.matrix(site)), clust, match.force=FALSE)
    #return
    return(c(FD=FD,FDmpd=FDmpd))
  }
}


species_matrix<-matrix(floor(rlnorm(500,1,1)),ncol=100,nrow=50,dimnames=list(paste0("Site",1:50),paste0("Species",1:100)))
trait_table<-data.frame(T1=rnorm(100,0,2),T2=runif(100,-2,2),T3=factor(sample(c("A","B","C","D"),100,replace=TRUE)),row.names=paste0("Species",1:100))

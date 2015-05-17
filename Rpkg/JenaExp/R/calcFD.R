###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
# FD calculations used in Wagg et al. manuscript in prep. for analysis of traits and biodiversity effects in the Jena TBE data from 2012
#modified by Lionel Hertzog on 14.05.2015 to include non-numeric trait data and FD index computation
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

# clac.FXxx = should the function return one or several of the FD index computed by the function dbFD in the FD package. Note if calc.FDiv=TRUE, calc.FRic must also be TRUE

# Note branch lengths can be better manipulated (ie for standarization etc...) . see https://github.com/ibartomeus/fundiv/blob/master/R/Xtree.R

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### The function
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 

calcFD <- function(species_matrix,trait_table,weighted=FALSE,distance = "euclidean",standardize = "None",calc.FRic=FALSE,calc.FDiv=FALSE,calc.FEve=FALSE,calc.FDis=FALSE,corr="sqrt") {
  
# # # # # # # Generate trait distance matrix # # # # # # # # #
  #apply standardization if asked for
  if(standardize!="None"){
    trait_table <- decostand(trait_table, method=standardize)
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
	if (!weighted){
	  species_matrix[species_matrix>0] <- 1
	}
		

# set up NMDS distance matrix
	mds <- metaMDS(D, k=2,)
	mds.points <- mds$points # save mds points for output

	# calculate distances between points in the MDS space 
	#(I prefer "euclidean" for ease and options here, but could also be calculated by hand using sqrt((x1-x2)^2 + (y1-y2)^2)  )
	MDSdist = as.matrix(vegdist(mds$points, method=distance)) 	
	
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # Go through each sampled sites and compute FD and FDmpd
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 output <- apply(species_matrix,1,function(x) helper_calcFD(x,MDSdist,trait_table))
 output <- as.data.frame(t(output))


####################
#if asked for compute other FD index from the dbFD function
##########

 if(calc.FRic | calc.FDiv | calc.FEve | calc.FDis){
   #for large number of species can take a long time to run or even crash....
  fds<-dbFD(x = trait_table,a = species_matrix,w.abun = weighted,calc.FRic = calc.FRic,calc.FDiv = calc.FDiv,corr=corr)
  output$Rich<-fds$nbsp
  
  if(calc.FDis){
    output$FDis<-fds$FDis
  }
  if(calc.FEve){
    output$FEve<-fds$FEve
  }    
  if(calc.FRic){
    output$FRic<-fds$FRic
  }
  if(calc.FDiv){
    output$FDiv<-fds$FDiv
  }
 }

 return(output)
}




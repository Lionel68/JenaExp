###################################
##Script playing with FD indices##
 ################################

library(vegan)
library(FD)
#generate some random trait data with different types and sd
traits<-data.frame(Phenology=factor(sample(1:5,50,replace=TRUE,prob = c(0.5,0.3,0.2,0.1,0.1))),Height=rlnorm(50,3,1),SLA=rnorm(50,20,5),Seed_mass=rlnorm(50,1,0.5),stringsAsFactors = FALSE)
rownames(traits)<-paste0("Sp",1:50)
#generate a site x species community matrix
comm<-matrix(sample(0:1,20*50,replace=TRUE),ncol=50,nrow=20)
#add some names
rownames(comm)<-paste("Site",1:20,sep="")
colnames(comm)<-paste("Sp",1:50,sep="")

#take the Gower distance between each species to feed into the FD calculation using podani extension to ordinal variables
dis<-gowdis(traits,ord="podani")
attr(dis,"Labels")<-paste0("Sp",1:50)
#compute the FDis, FRic, FEeve and FDiv
fd<-dbFD(x = traits,a=comm,w.abun=TRUE,stand.x=TRUE) #takes long time??

#function to compute the slope of the mimimum spanning tree
a_mst<-function(traits,comm){
  #compute mst for each sample
  mst<-apply(comm,1,function(x){
    trait<-traits[rownames(traits)%in%names(x)[x!=0],]
    dis<-gowdis(trait,ord="podani")
    mst<-spantree(dis)
  })
  #get for each sample the distance between the species
  l<-lapply(mst,function(x) x$dist)
  #order the distance
  l<-lapply(l,function(x) x[order(x)])
  #get some rank
  rank<-lapply(l,function(x) 1:length(x))
  #compute cumulative distance
  l<-lapply(l,function(x) cumsum(x))
  #estimate the slopes
  a<-mapply(FUN=function(x,y) coef(lm(y~log(x)))[2],x=rank,y=l)
  names(a)<-rownames(comm)
  return(a)
}

as<-a_mst(traits,comm)

#put all FD together
all<-data.frame(S=fd$nbsp,FRic=fd$FRic,FDiv=fd$FDiv,FDis=fd$FDis,FEve=fd$FEve,a_MST=as)


##########################
#next step would be to have only two traits
#and to generate illustration on how varying abundances, traits distribution affect
#the FD indices, one aim for this would be to have then a Shiny App
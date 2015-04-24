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

#generate two traits for 5 species
traits<-data.frame(Size=rlnorm(5,1,1),Mass=rnorm(5,10,1))
traits<-apply(traits,2,scale)
rownames(traits)<-paste0("Sp",1:5)
#generate community abundance data (only for one site)
abd<-matrix(rpois(5,3),ncol=5,nrow=1)
colnames(abd)<-paste0("Sp",1:5)
#compute the FDs
fd<-dbFD(traits,abd)
#compute a_MST
fd$a_mst<-a_mst(traits,abd)
#compute the weighted centroids in trait space
p<-abd/sum(abd)
ck<-apply(traits,2,function(x) sum(p*x))

#plot FDis
plot(Mass~Size,traits,pch=16,cex=15*(abd/sum(abd)),main=paste("FDis:",round(fd$FDis,2),sep=" "))
points(ck[1],ck[2],pch=13,cex=2,col="red")
apply(traits,1,function(x) lines(c(ck[1],x[1]),c(ck[2],x[2])))

#plot FRic
chv<-t(convhulln(traits))
poly<-cbind(traits[chv,1],traits[chv,2])
poly<-poly[!duplicated(poly),]
plot(Mass~Size,traits,pch=16,cex=15*(abd/sum(abd)),main=paste("FRic:",round(fd$FRic,2),sep=" "))
polygon(poly[,1],poly[,2])

#plot FDiv
edge<-traits[traits[,1]%in%poly[,1] & traits[,2]%in%poly[,2],]
g<-apply(edge,2,function(x) sum(x)/length(x))
#average square rooted distance from each species to this centroid
dg<-sqrt(rowSums((traits-g)**2))
dg_avg<-sum(dg)/dim(traits)[1]
#abundance weighted difference for each species from this average centroids
d<-sum((dg-dg_avg)*p)
d_abs<-sum(abs(dg-dg_avg)*p)
lims<-c(g[1]-dg_avg-0.5,g[1]+dg_avg+0.5,g[2]-dg_avg-0.5,g[2]+dg_avg+0.5)
plot(Mass~Size,traits,pch=16,cex=15*(abd/sum(abd)),main=paste("FDiv:",round(fd$FDiv,2),sep=" "),xlim=c(lims[1],lims[2]),ylim=c(lims[3],lims[4]))
#average distance from centroids
draw.circle(g[1],g[2],radius=dg_avg,lwd = 3)
#centroids
points(g[1],g[2],pch=13,col="red",cex=2)
#distance between each species and the average distance
#still not resolved ....

#FEve
mst<-spantree(gowdis(traits))
p<-plot(mst,main=paste("FEve:",round(fd$FEve,2),sep=" "))
points(p$sites[,1],p$sites[,2],pch=16,cex=15*(abd/sum(abd)))



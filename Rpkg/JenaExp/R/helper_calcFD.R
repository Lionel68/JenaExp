#this is an helper function to use in the apply context to go through all sites
helper_calcFD<-function(site,MDSdist,trait_table){
  #remove non present species
  site<-site[site!=0 & !is.na(site)]
  #if monoculture or no species in sites return 0
  if(length(site)==1){
    return(c(FD=0,FDmpd=0))
  }
  else{
    ps<-site/sum(site)
    #compute FDmpd
    smat<-1-as.matrix(dist(ps,method="manhattan"))
    d_mds<-MDSdist[names(ps),names(ps)]*smat 
    FDmpd<-mean(as.numeric(d_mds[lower.tri(d_mds)]))
    #compute FD
    trait_sub<-trait_table[rownames(trait_table)%in%names(ps),]
    #if all traits are numeric use euclidian distance
    if(all(apply(trait_sub,2,is.numeric))){
      Dsub<-dist(trait_sub)
    }
    #otherwise use Gower distance
    else{
      Dsub<-gowdis(trait_sub)
    }    
    clust<-hclust(Dsub,method="complete")
    FD<- treedive(t(as.matrix(site)), clust, match.force=FALSE)
    
    return(c(FD=FD,FDmpd=FDmpd))
  }
}

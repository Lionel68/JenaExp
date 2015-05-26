sim_BEF<-function(n_species=10,n_rep=5,gradient=1:10){
  #construct the richness gradient
  if(length(n_rep)!=1 & length(n_rep)!=length(gradient)){
    print("The number of parameter argument must either be of length 1 or of the same length as the richness gradient")
    return(invisible(NULL))
  }
  else{
    #each element of this vector is a plot with the defined richness level
    rich<-rep(gradient,n_rep)
    n_rep<-as.data.frame(table(rich),stringsAsFactors=FALSE)
    n_rep$rich<-as.numeric(n_rep$rich)
  }
  #create a species list
  if(n_species<=26){
    sp<-LETTERS[1:n_species]
  }
  else{
    sp<-unlist(lapply(1:26,function(x) paste(LETTERS[x],letters[1:26],sep="")))[1:n_species]
  }
  #randomly choose from all possible combination at each level of diversity
  #when the number of possible combination is lower than the number of replication, single species
  #combination may be appearing more than once
  Exp<-list()
  for(i in 1:dim(n_rep)[1]){
    #when the number of possible combination is bigger than 1000 rely on self-made function
    if(choose(length(sp),n_rep[i,1])<=1000){
      #all possible combination
      comb<-combn(sp,n_rep[i,1],FUN=function(x) paste(x,collapse=" "))
      if(length(comb)<n_rep[i,2]){
        Exp[[i]]<-comb[sample(1:length(comb),n_rep[i,2],replace=TRUE)]
      }
      else{
        #sample the number of asked replication from the possible combination
        Exp[[i]]<-comb[sample(1:length(comb),n_rep[i,2],replace=FALSE)]
      }
    }
    #if the number of possible combination is too big use a sampling alternatives to speed up things
    else{
      comb<-sapply(1:n_rep[i,2],function(x) paste(sample(sp,n_rep[i,1]),collapse=" "))
      #check if two combinations are not repeated
      if(length(unique(comb))<n_rep[i,2]){
        comb<-unique(comb)
        diff<-n_rep[i,2]-length(comb)
        comb<-c(comb,sapply(1:diff,function(x) paste(sample(sp,n_rep[i,1]),collapse=" ")))
      }
      Exp[[i]]<-comb
    }
  }
  
  #generate the site x species matrix
  comm<-matrix(0,nrow=n_species,ncol=length(rich),dimnames=list(sp,paste("Site",1:length(rich))))
  #add ones in the sites where a species is occuring
  is_present<-lapply(strsplit(unlist(Exp)," "),function(x) sp%in%x)
  comm[unlist(is_present)]<-1  
  comm<-t(comm)
  return(comm)
}
  
  
 







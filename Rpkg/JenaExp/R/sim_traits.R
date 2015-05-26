sim_traits<-function(n_traits=5,n_species=10,distr="rnorm",params=NULL,sp_names=NULL){
  
  #if no species names are provided create them
  if(is.null(sp_names)){
    sp_names<-paste0("Sp",1:n_species)
  }
  
  #if the distr argument is wrongly set
  if(length(distr)!=1 & length(distr)!=n_traits){
    print("The parameter distr should either be of length one or of length equal to the n_traits parameters")
    return(invisible(NULL))
  }
  
  #get the distribution provided by the users
  fun_distr <- tryCatch(lapply(distr,match.fun), error=function(w) {print("Wrong trait distribution");return(NULL)})
  if(is.null(fun_distr)){
    return(invisible(fun_distr))
  }
  
  #if no parameters are given the default distribution values are used for the continuous distribution
  #for the factorial traits per default 4 levels are used
  if(is.null(params)){
    #if all traits distribution are continuous
    if(length(grep("factor",distr))==0){
      if(length(fun_distr)==1){
        df<-matrix(sapply(fun_distr,FUN=function(x) x(n_species*n_traits)),ncol=n_traits,nrow=n_species,dimnames=list(sp_names,paste0("Traits",1:n_traits)))
      }
      else{
        df<-matrix(sapply(fun_distr,FUN=function(x) x(n_species)),ncol=n_traits,nrow=n_species,dimnames=list(sp_names,paste0("Traits",1:n_traits)))
      }
    }
    #if some traits are factorial
    else{
      if(length(fun_distr)==1){
        df<-matrix(sapply(fun_distr,FUN=function(x) x(sample(1:4,n_species*n_traits,replace=TRUE),levels=1:4)),ncol=n_traits,nrow=n_species,dimnames=list(sp_names,paste0("Traits",1:n_traits)))
      }
      else{
        df<-matrix(0,ncol=n_traits,nrow=n_species,dimnames=list(sp_names,paste0("Trait",1:n_traits)))
        for(i in 1:n_traits){
          if(distr[i]=="factor"){
            df[,i]<-fun_distr[[i]](sample(1:4,n_species,replace=TRUE),levels=1:4)
          }
          else{
            df[,i]<-fun_distr[[i]](n_species)
          }
        }
      }
    }
  }
  #if parameters are given
  else{
    #if all traits distribution are continuous
    if(length(grep("factor",distr))==0){
      if(length(fun_distr)==1){
        df<-matrix(mapply(FUN=function(x,y) x(n_species*n_traits,y[1],y[2]),x=fun_distr,y=params),ncol=n_traits,nrow=n_species,dimnames=list(sp_names,paste0("Traits",1:n_traits)))
      }
      else{
        df<-matrix(mapply(FUN=function(x,y) x(n_species,y[1],y[2]),x=fun_distr,y=params),ncol=n_traits,nrow=n_species,dimnames=list(sp_names,paste0("Traits",1:n_traits)))
      }
    }
    #if some traits are factorial
    else{
      if(length(fun_distr)==1){
        df<-matrix(mapply(FUN=function(x,y) x(sample(1:y[[1]],n_species*n_traits,replace=TRUE,prob=y[[2]]),levels=1:y[[1]]),x=fun_distr,y=params),ncol=n_traits,nrow=n_species,dimnames=list(sp_names,paste0("Traits",1:n_traits)))
      }
      else{
        df<-matrix(0,ncol=n_traits,nrow=n_species,dimnames=list(sp_names,paste0("Trait",1:n_traits)))
        for(i in 1:n_traits){
          if(distr[i]=="factor"){
            df[,i]<-fun_distr[[i]](sample(1:params[[i]][[1]],n_species,replace=TRUE,prob=params[[i]][[2]]),levels=1:params[[i]][[1]])
          }
          else{
            df[,i]<-fun_distr[[i]](n_species,params[[i]][1],params[[i]][[2]])
          }
        }
      }
    }
  }  
  return(df)    
}
  
  
  

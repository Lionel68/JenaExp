##############################
#function to decompose the BEF effect
#following Fox and Kerr 2012 Oikos Paper
###########################################


price<-function(data,control,avg=TRUE){
  
  #some basic values
  s <- apply(data,1,function(x) length(x[x!=0])) #species richness in comparison sites
  s_dot <- apply(control,1,function(x) length(x[x!=0])) #species richness in control sites
  s_c <- apply(control,1,function(x) apply(data,1,function(y) length(which(x!=0 & y!=0)))) #a table with the number of species shared between each control plot and the comparison plot
  s_c[s_c==0]<-NA #remove plot when there is no species in common from the analysis
  
  z <- apply(data,1,mean) #the average species contribution at the comparison sites
  z_dot <- apply(control,1,mean) #the average species contribution at the control sites
  z_c <- apply(control,1,function(x) apply(data,1,function(y) mean(y[x!=0 & y!=0]))) #the average contribution of species present at both sites in the comparison sites
  z_c_dot <- apply(control,1,function(x) apply(data,1,function(y) mean(x[x!=0 & y!=0]))) #the average contribution of species present at both sites in the control sites
  
  
  
  
  srl <- (s_c-s)*z #the random species loss term
  if(avg){
    srl<-apply(srl,1,mean,na.rm=TRUE)
  }
  
  srg <- (s_dot-s_c)*z_dot #the random species gain term
  if(avg){
    srg<-apply(srg,1,mean,na.rm=TRUE)
  }
  
  scl <- (z_c-z)*s_c #the species composition lost effect
  if(avg){
    scl<-apply(scl,1,mean,na.rm=TRUE)
  }
  
  scg <- -(z_c_dot-z_dot)*s_c #the species composition gain effect
  if(avg){
    scg <- apply(scg,1,mean,na.rm=TRUE)
  }
  
  cde <- (z_c_dot-z_c)*s_c #the context dependance effect
  if(avg){
    cde <- apply(cde,1,mean,na.rm=TRUE)
  }
  
  df <- data.frame(Plot=rownames(data),SRL=srl,SRG=srg,SCL=scl,SCG=scg,CDE=cde)
  rownames(df)<-1:dim(df)[1]
  return(df)  
}
#server

#load libraries
library(FD)
library(MASS)
library(shiny)
library(geometry)
library(plotrix)
#functions
sim_traits<-function(S,var_t1=1,var_t2=1){
  tmp<-data.frame(T1=rlnorm(S,1,var_t1),T2=rnorm(S,10,var_t2))
  rownames(tmp)<-paste0("Sp",1:S)
  return(tmp)  
}

sim_comm<-function(S,var=1){
  tmp<-matrix(rnegbin(S,5,1/var)+1,ncol=S,nrow=1)
  colnames(tmp)<-paste0("Sp",1:S)
  return(tmp)
}

ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
}  

#function to compute the slope of the mimimum spanning tree
a_mst<-function(traits,comm){
  #compute mst for each sample
  dis<-gowdis(traits,ord="podani")
  mst<-spantree(dis)  
  #get for each sample the distance between the species
  l<-mst$dist
  #order the distance
  l<-l[order(l)]
  #get some rank
  rank<-1:length(l)
  #compute cumulative distance
  l<-cumsum(l)
  dat<-data.frame(Rank=rank,l=l)
  a<-coef(lm(log(l)~log(Rank),dat))[2]  
  return(list(coeff=a,dat=dat))
}

#from users input compute FD index 
shinyServer(function(input,output){
  r_comm<-reactive({   
    comm<-sim_comm(input$S,input$var)
  })  
  r_traits<-reactive({
    traits<-sim_traits(input$S,input$var_t1,input$var_t2)
    
  })
  r_fd<-reactive({
    fd<-dbFD(r_traits(),r_comm(),stand.x=TRUE)
  })
 
  output$Plot1<-renderPlot({
    #the trait data
    dat<-r_traits()
    names(dat)<-c("T1","T2")
    #scale to mean 0 and unit standard deviation
    dat<-apply(dat,2,scale)
    #the community data
    com<-r_comm()
    ps<-(com/sum(com))
    #the weighted centroid coordinates
    cs<-apply(dat,2,function(x) sum(x*ps))
    #the vertices from the convex hull
    hpts <- chull(x = dat[,1], y = dat[,2])
    edge<-dat[hpts,]
    #the centroid based only on the vertices position
    g<-apply(edge,2,function(x) sum(x)/length(x))
    #average square rooted distance from each species to this centroid
    dg<-sqrt(rowSums((dat-g)**2))
    dg_avg<-sum(dg)/dim(dat)[1]
    #abundance weighted difference for each species from this average centroids
    d<-sum((dg-dg_avg)*ps)
    d_abs<-sum(abs(dg-dg_avg)*ps)
    #the minimum spanning tree for FEve
    mst<-spantree(gowdis(dat))
    
    par(mfrow=c(1,2))
    #FDis
    plot(T2~T1,dat,pch=16,cex=15*ps,main=paste("FDis = ",round(r_fd()$FDis,2),sep=""))    
    points(cs[1],cs[2],pch=13,col="red",cex=2)
    apply(dat,1,function(x) lines(c(cs[1],x[1]),c(cs[2],x[2])))
    #FRic
    plot(T2~T1,dat,pch=16,cex=15*ps,main=paste("FRic = ",round(r_fd()$FRic,2),sep=""))
    ConvexHull(dat[,1],dat[,2],"black")
   
  })
  
  output$Plot2<-renderPlot({
    dat<-r_traits()
    names(dat)<-c("T1","T2")
    #scale to mean 0 and unit standard deviation
    dat<-apply(dat,2,scale)
    #the community data
    com<-r_comm()
    ps<-(com/sum(com))
    #the weighted centroid coordinates
    cs<-apply(dat,2,function(x) sum(x*ps))
    #the vertices from the convex hull
    hpts <- chull(x = dat[,1], y = dat[,2])
    edge<-dat[hpts,]
    #the centroid based only on the vertices position
    g<-apply(edge,2,function(x) sum(x)/length(x))
    #average square rooted distance from each species to this centroid
    dg<-sqrt(rowSums((dat-g)**2))
    dg_avg<-sum(dg)/dim(dat)[1]
    #abundance weighted difference for each species from this average centroids
    d<-sum((dg-dg_avg)*ps)
    d_abs<-sum(abs(dg-dg_avg)*ps)
    #the minimum spanning tree for FEve
    mst<-spantree(gowdis(dat))
    
    par(mfrow=c(1,2))
    #FDiv
    plot(T2~T1,dat,pch=16,cex=15*ps,main=paste("FDiv = ",round(r_fd()$FDiv,2),sep=""))
    #average distance from centroids
    draw.circle(g[1],g[2],radius=dg_avg,lwd = 3)
    #centroids
    points(g[1],g[2],pch=13,col="red",cex=2)
    #FEve
    p<-plot(mst,main=paste("FEve = ",round(r_fd()$FEve,2),sep=""))
    points(p$sites[,1],p$sites[,2],pch=16,cex=15*ps)
  })
  
  output$Plot3<-renderPlot({
    dat<-r_traits()
    names(dat)<-c("T1","T2")
    dat<-apply(dat,2,scale)
    com<-r_comm()
    #compute the slope of the cumulative distance
    a<-a_mst(dat,com)
    plot(log(l)~log(Rank),a$dat,xlab="Rank of the branches (log)",ylab="Cumulative branch length (log)",main=paste("Slope of MST = ",round(a$coeff,2),sep=""),pch=16)
    
  })
  
  
  
  #output$text<-renderText(r_fd()$FDis)  
})

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
    dat<-r_traits()
    names(dat)<-c("T1","T2")
    dat<-apply(dat,2,scale)
    com<-r_comm()
    ps<-(com/sum(com))
    plot(T2~T1,dat,pch=16,cex=15*ps,main=paste("FDis = ",round(r_fd()$FDis,3),sep=""))
    cs<-apply(dat,2,function(x) sum(x*ps))
    points(cs[1],cs[2],pch=13,col="red",cex=2)
    apply(dat,1,function(x) lines(c(cs[1],x[1]),c(cs[2],x[2])))
  })
  
  output$Plot2<-renderPlot({
    dat<-r_traits()
    names(dat)<-c("T1","T2")
    dat<-apply(dat,2,scale)
    com<-r_comm()
    ps<-(com/sum(com))
    plot(T2~T1,dat,pch=16,cex=15*ps,main=paste("FRic = ",round(r_fd()$FRic,3),sep=""))
    ConvexHull(dat[,1],dat[,2],"black")
  })
  
  output$Plot3<-renderPlot({
    dat<-r_traits()
    names(dat)<-c("T1","T2")
    dat<-apply(dat,2,scale)
    com<-r_comm()
    ps<-(com/sum(com))
    hpts <- chull(x = dat[,1], y = dat[,2])
    edge<-dat[hpts,]
    g<-apply(edge,2,function(x) sum(x)/length(x))
    #average square rooted distance from each species to this centroid
    dg<-sqrt(rowSums((dat-g)**2))
    dg_avg<-sum(dg)/dim(dat)[1]
    #abundance weighted difference for each species from this average centroids
    d<-sum((dg-dg_avg)*ps)
    d_abs<-sum(abs(dg-dg_avg)*ps)
    lims<-c(g[1]-dg_avg-0.5,g[1]+dg_avg+0.5,g[2]-dg_avg-0.5,g[2]+dg_avg+0.5)
    plot(T2~T1,dat,pch=16,cex=15*ps,main=paste("FDiv = ",round(r_fd()$FDiv,3),sep=""),xlim=c(lims[1],lims[2]),ylim=c(lims[3],lims[4]))
    #average distance from centroids
    draw.circle(g[1],g[2],radius=dg_avg,lwd = 3)
    #centroids
    points(g[1],g[2],pch=13,col="red",cex=2)
  })
  
  output$Plot4<-renderPlot({
    dat<-r_traits()
    names(dat)<-c("T1","T2")
    com<-r_comm()
    ps<-com/sum(com)
    mst<-spantree(gowdis(dat))
    p<-plot(mst,main=paste("FEve = ",round(r_fd()$FEve,3),sep=""))
    points(p$sites[,1],p$sites[,2],pch=16,cex=15*ps)    
  })
  
  #output$text<-renderText(r_fd()$FDis)  
})

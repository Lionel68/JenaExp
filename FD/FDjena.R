###############################
#Script looking at the relations
#between FDjena and other FD index
####################################

#load the libraris
library(FD)
library(psych)
library(ggplot2)
library(reshape2)

#load traits data, excel ID: 11259
traits<-read.table("/home/lionel/Documents/PhD/Data/FD/TBE_traits.csv",sep=",",head=TRUE,stringsAsFactors = FALSE)
#load plot infos
plot<-read.table("/home/lionel/Documents/PhD/Data/tbe_info.csv",sep=",",head=TRUE,stringsAsFactors = FALSE)
#load some cover data
#could not find cover data from the TBE yet ...

#compute FD
comm_bin<-plot[,19:38] #only keep the plot x species matrix
traits_d<-traits[traits$Species%in%names(comm_bin),c(1:6)] #only keep species sown in the TBE
rownames(traits_d)<-traits[traits$Species%in%names(comm_bin),9]

fd<-dbFD(traits_d,comm_bin)
fd<-data.frame(Plot=plot$plot,Rich=fd$nbsp,FDjena=plot$FD,FRic=fd$FRic,FDis=fd$FDis,FDiv=fd$FDiv,FEve=fd$FEve)
pairs.panels(fd[,-1])
fd_m<-melt(fd,id.vars = 1:3)
ggplot(fd_m,aes(x=FDjena,y=value))+geom_point(position=position_jitter(w=0.1,h=0))+facet_wrap(~variable,scales = "free")+stat_summary(fun.data=mean_cl_normal,geom="errorbar",color="red",width=.1,na.rm=TRUE)+stat_summary(fun.y=mean,na.rm=TRUE,geom="point",col="red")+stat_smooth()

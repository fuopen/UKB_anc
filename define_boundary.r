library(sp)
library(fastcluster)
library(raster)
library(colorRamps)
library(parallel)

if(!exists('d1')){
    d1<-readRDS('make_boundary.rds')
    county.spdf.list<-d1[[1]]
    v2.cl.srb.spdf<-d1[[2]]
    Gb.counties<-d1[[3]]
}
make.boundary<-function(){
    res<-sapply(county.spdf.list,function(x)attr(x,'county.mean.ancestry'))

    res.large<-apply(res,2,function(x)rownames(res)[which.max(x)])

    temp=v2.cl.srb.spdf@data

    ourinds=which(rowSums(temp[,1:23]>0.5)>0)

    whichbest=1:nrow(temp)*0
    maxes=whichbest
    for(i in 1:23){
        whichbest[temp[,i]>maxes]=i
        maxes[temp[,i]>maxes]=temp[temp[,i]>maxes,i]
    }

    spincounties<-v2.cl.srb.spdf%over%Gb.counties

    newcounties=spincounties[ourinds,"NAME_2"]

    infcounties=table(newcounties,colnames(temp[ourinds,1:23])[whichbest[ourinds]])

    whichbest2=1:nrow(infcounties)*0
    maxes=whichbest2
    for(i in 1:22){
        whichbest2[infcounties[,i]>maxes]=i
        maxes[infcounties[,i]>maxes]=infcounties[infcounties[,i]>maxes,i]
    }
    newlabels=colnames(infcounties)[whichbest2]
    names(newlabels)=rownames(infcounties)
    oldlabels=res.large[names(newlabels)]

    return(list(newlabels,oldlabels))
}

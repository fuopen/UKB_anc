library(parallel)
library(fastcluster)
library(sp)
library(raster)
library(colorRamps)
library(igraph)
library(gplots)

load.map<-function(){
    gbirl<<-readRDS('data/GB_IRELAND.rds')
    gb.map<<-readRDS('data/new_boundary.rds')
    irl<<-gbirl[5,]
    my.box.gb<<-bbox(gbirl)
}

if(!exists('spdf')){
    spdf<-readRDS('v2_self_BI_487409.rds')
    counties<-readRDS('self_BI_counties_487409.rds')
}

if(!exists('gb.bound')){
    gb.bound<-readRDS('data/New_GB_boundaries.rds')
    load.map()
}

if(!exists('Sp.gb.data')){
    SpInCounties<-spdf%over%gb.map
    Sp.data<-spdf@data
    rownames(Sp.data)<-rownames(counties)
    sna<-is.na(SpInCounties$assign)
    SpInCounties<-SpInCounties[!sna,]
    Sp.data<-Sp.data[!sna,]
    Sp.gb.data<-Sp.data[,1:23]
    Sp.gb.data$G_Cornwall_Tip<-Sp.gb.data$G_Cornwall_Tip+Sp.gb.data$G_Cornwall
    Sp.gb.data$G_Cornwall<-NULL
}


raster.boundary<-function(n=200){
    panel.box<-rbind(c(-11,2),c(50,61.5))
    eval(parse(text=paste0('base.raster',n,'<-raster(nr=n,ncol=n,xmn=panel.box[1,1],xmx=panel.box[1,2],ymn=panel.box[2,1],ymx=panel.box[2,2],crs=CRS("+init=epsg:4326"))')))
    eval(parse(text=paste0('base.raster',n,'[]<-1:ncell(base.raster',n,')')))
    raster.ls<-list()
    for(k in names(gb.bound)){
        eval(parse(text=paste0(k,'.raster',n,'<<-rasterize(gb.bound[["',k,'"]],base.raster',n,')')))
        eval(parse(text=paste0(k,'.raster',n,'<<-mask(base.raster',n,',',k,'.raster',n,')')))
        eval(parse(text=paste0('raster.ls[["',k,'"]]<-',k,'.raster',n)))
        #eval(parse(text=paste0(k,'.raster',n,'.idx<<-which(!is.na(getValues(',k,'.raster',n,')))')))
    }
    return(raster.ls)
}

if(!exists('gb22.ras')){
    gb22.ras<-readRDS('gb_22_counties_1000.rds')
}

if(!exists('cornwall.ras')){
    cornwall.ras<-readRDS('cornwall_split_1000.rds')
}

merge.raster<-function(){
    gs<-gb22.ras
    gs$G_Cornwall<-cornwall.ras$G_Cornwall
    gs$G_Cornwall_Tip<-cornwall.ras$G_Cornwall_Tip
    return(gs)
}

if(!exists('mr')){
    mr<-merge.raster()
}

color=grDevices::colors()[grep('gr(a|e)y',grDevices::colors(),invert=T)]
if(!exists('gb.color')){
    gb.color<-readRDS('23_regions_color1.rds')
    names(gb.color)<-names(mr)
}

official.names<-c('Anglia','South Central England','Central England','North Central England','North Yorkshire','Northumberland','South East England','Merseyside','South England','South East Wales','North East Scotland','South West Wales','Devon','North West Wales','South Yorkshire','Cornwall Tip','Orkney','Lincolnshire','Cumbria','North West Scotland','North Ireland and South Scotland','Ireland','Cornwall')

names(official.names)<-names(gb.color)

plot.boundary.raster<-function(){
    pdf('raster_gb_boundary23.pdf')
    par(oma=c(0,0,0,0),mar=c(0,0,0,0))
    plot(0,0,type='n',axes=F,ann=F,xlim=c(-11,2),ylim=c(50,61.5))
    for(i in names(mr)){
        plot(mr[[i]],col=gb.color[[i]],add=T,legend=F,xaxs='i',yaxs='i')
    }
    
    legend('top',legend=official.names,ncol=4,fill=gb.color[names(official.names)],cex=0.5,bty='n')

    dev.off()
}

plot.boundary.raster.seprate<-function(){
    pdf('raster_gb_boundary23_separate.pdf')
    par(oma=c(0,0,0,0),mar=c(0,0,0,0))
    plot(0,0,type='n',axes=F,ann=F,xlim=c(-11,2),ylim=c(50,61.5))
    for(i in names(mr)){
        plot(mr[[i]],col=gb.color[[i]],add=T,legend=F,xaxs='i',yaxs='i')
    }
    dev.off()

    #legend('top',legend=official.names,ncol=4,fill=gb.color[names(official.names)],cex=0.5,bty='n')

    #dev.off()
}

#plot.boundary.raster<-function(n){
#    rb<-raster.boundary(n)
#    pdf(paste0('raster_gb_boundar',n,'.pdf'))
#    par(oma=c(0,0,0,0),mar=c(0,0,0,0))
#    plot(0,0,type='n',axes=F,ann=F,xlim=c(-11,2),ylim=c(50,61.5))
#    for(i in names(gb.bound)){
#        plot(get(paste0(i,'.raster',n)),col=gb.color[[i]],add=T,legend=F,xaxs='i',yaxs='i')
#    }
#    
#    legend('top',legend=official.names,ncol=4,fill=gb.color[names(official.names)],cex=0.5,bty='n')
#
#    dev.off()
#}

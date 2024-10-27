library(sp)
library(fastcluster)
library(raster)
library(colorRamps)
library(parallel)

color=grDevices::colors()[grep('gr(a|e)y',grDevices::colors(),invert=T)]

#gb.color<-readRDS('data/23_regions_color1.rds')

if(!exists('England_Wales')){
    gbirl<-readRDS('data/GB_IRELAND.rds')
    panel.box<-bbox(gbirl)
}

if(!exists('GB.counties')){
    gb.counties<-readRDS('data/GBR_adm2.rds')
    gb.counties<-spTransform(gb.counties,CRS('+init=epsg:4326'))
}


if(!exists('v2.cl.srb.spdf')){
    long_lat.dt<-read.table('long_lat_ukbiobank.tsv',as.is=T,header=T,sep="\t")[,1:2]

    v2.eth<-readRDS('biobank_v2_eth.background.rds')
    v2<-readRDS('biobank_v2_results.rds')
    v2.cl<-v2[[1]]
    
    uk.countries<-c('England','Scotland','Wales')
    norirl<-'Northern_Ireland'
    irl<-'Republic_of_Ireland'

    BI<-c('British','Irish')
    nonNA<-!is.na(long_lat.dt[[1]])
    ESW.ids<-nonNA & (!is.na(v2.eth[[2]])) & v2.eth[[2]]%in%uk.countries & v2.eth[[4]]%in%BI
    NI.ids<-(!is.na(v2.eth[[2]])) & v2.eth[[2]]==norirl & v2.eth[[4]]%in%BI
    IRL.ids<-(!is.na(v2.eth[[2]])) & v2.eth[[2]]==irl & v2.eth[[4]]%in%BI
    
    ESW.cord<-long_lat.dt[ESW.ids,1:2]
    esw.cl<-v2.cl[ESW.ids,]
    noi.cl<-v2.cl[NI.ids,]
    irl.cl<-v2.cl[IRL.ids,]

    colnames(ESW.cord)<-c('long','lat')
    coordinates(ESW.cord)<- ~ long+lat
    proj4string(ESW.cord)<-proj4string(gbirl)
    
    noi.mean.ancestry<-colMeans(noi.cl)
    irl.mean.ancestry<-colMeans(irl.cl)

    v2.cl.srb.spdf<-SpatialPointsDataFrame(ESW.cord,data=data.frame(esw.cl))
}

if(!exists('county.spdf.list')){
    v2InCounties<-v2.cl.srb.spdf%over%gb.counties

    v2.counties.spdf.nonNA<-!is.na(v2InCounties$NAME_2)

    v2InCounties.nonNA<-v2InCounties[v2.counties.spdf.nonNA,]
    v2.counties.spdf<-v2.cl.srb.spdf[v2.counties.spdf.nonNA,]
    county.sp.sizes<-table(v2InCounties$NAME_2)
    county.spdf.list<-mclapply(names(county.sp.sizes),function(county.name){
        if(county.sp.sizes[county.name]==0){
            return(NULL)
        }
        county.spdf<-v2.counties.spdf[v2InCounties.nonNA$NAME_2==county.name,]
        attr(county.spdf,'county.name')<-county.name
        attr(county.spdf,'county.sample.size')<-county.sp.sizes[county.name]
        attr(county.spdf,'county.mean.ancestry')<-colMeans(county.spdf@data[,1:23])
        attr(county.spdf,'county.meangb.percent')<-c(mean(rowSums(county.spdf@data[,1:23])),mean(rowSums(county.spdf@data[,24:127])))
        return(county.spdf)
    },mc.cores=12)
    county.spdf.list<-county.spdf.list[sapply(county.spdf.list,function(x)!is.null(x))]
    names(county.spdf.list)<-sapply(county.spdf.list,function(x)attr(x,'county.name'))
}

if(!exists('cornwall.spdf')){
    cornwall.spdf<-v2.cl.srb.spdf[,c('G_Cornwall','G_Cornwall_Tip')]
}

make.raster.gbirl<-function(n=250){
    panel.box<-rbind(c(-11,2),c(50,61.5))
    eval(parse(text=paste0('base.raster',n,'<-raster(nr=n,ncol=n,xmn=panel.box[1,1],xmx=panel.box[1,2],ymn=panel.box[2,1],ymx=panel.box[2,2],crs=CRS("+init=epsg:4326"))')))
    eval(parse(text=paste0('base.raster',n,'[]<-1:ncell(base.raster',n,')')))
    eval(parse(text=paste0('gb.raster',n,'<<-rasterize(gbirl[c(1,3,4),],base.raster',n,')')))
    eval(parse(text=paste0('noi.raster',n,'<<-rasterize(gbirl[2,],base.raster',n,')')))
    eval(parse(text=paste0('irl.raster',n,'<<-rasterize(gbirl[5,],base.raster',n,')')))
    eval(parse(text=paste0('gb.raster',n,'<<-mask(base.raster',n,',gb.raster',n,')')))
    eval(parse(text=paste0('noi.raster',n,'<<-mask(base.raster',n,',noi.raster',n,')')))
    eval(parse(text=paste0('irl.raster',n,'<<-mask(base.raster',n,',irl.raster',n,')')))
    eval(parse(text=paste0('gb.raster',n,'.coords<<-coordinates(gb.raster',n,')')))
    eval(parse(text=paste0('gb.raster',n,'.idx<<-which(!is.na(getValues(gb.raster',n,')))')))
    eval(parse(text=paste0('noi.raster',n,'.coords<<-coordinates(noi.raster',n,')')))
    eval(parse(text=paste0('noi.raster',n,'.idx<<-which(!is.na(getValues(noi.raster',n,')))')))
    eval(parse(text=paste0('irl.raster',n,'.coords<<-coordinates(irl.raster',n,')')))
    eval(parse(text=paste0('irl.raster',n,'.idx<<-which(!is.na(getValues(irl.raster',n,')))')))
    cornwall<-gb.counties[gb.counties$NAME_2=='Cornwall',]
    eval(parse(text=paste0('cornwall.raster',n,'<<-rasterize(cornwall,base.raster',n,')')))
    eval(parse(text=paste0('cornwall.raster',n,'<<-mask(base.raster',n,',cornwall.raster',n,')')))
    eval(parse(text=paste0('cornwall.raster',n,'.idx<<-which(!is.na(getValues(cornwall.raster',n,')))')))
    eval(parse(text=paste0('cornwall.raster',n,'.coords<<-coordinates(cornwall.raster',n,')')))
}

get.weight.gbirl<-function(n,q,spdf){#at the moment,try constant
    if(!exists(paste0('gb.raster',n))){
        md<-make.raster.gbirl(n)
    }
    ras<-get(paste0('gb.raster',n))
    ras.c<-get(paste0('gb.raster',n,'.coords'))
    pixel.index<-get(paste0('gb.raster',n,'.idx'))

    ras.c.ukb<-ras.c[pixel.index,]
    bm.dt<-as.matrix(spdf@data)
    bm<-spdf@coords
    ras.c.ukb.mean<-do.call(cbind,mclapply(1:nrow(ras.c.ukb),function(i){
        X0<-ras.c.ukb[i,]
        bm.c<-sweep(bm,2,X0)
        bm.c<-rowSums(bm.c^2,na.rm=T)
        bm.w<-exp(-q*bm.c)
        d1<-sum(bm.w,na.rm=T)
        while(d1<200){
            q<-q/2
            bm.w<-exp(-q*bm.c)
            d1<-sum(bm.w,na.rm=T)
        }
        
        bm.w<-matrix(bm.w/d1,nrow=1)
        b<-bm.w%*%bm.dt
        return(b[1,])
    },mc.cores=12))
    return(ras.c.ukb.mean)
}

get.cornwall.pix<-function(n,q=50){
    cr<-get(paste0('cornwall.raster',n,'.idx'))
    gbr<-get(paste0('gb.raster',n,'.idx'))
    cornwall.weight<-get.weight.gbirl(n,q,cornwall.spdf)
    cornwall.anc.pix<-cornwall.weight[,match(cr,gbr)]
    cornwall.tip.pix<-get(paste0('cornwall.raster',n))
    cornwall.body.pix<-get(paste0('cornwall.raster',n))
    cornwall.tip.exc.idx<-cr[cornwall.anc.pix[2,]<0.21]
    cornwall.tip.inc.idx<-cr[cornwall.anc.pix[2,]>=0.21]
    cornwall.tip.pix[cornwall.tip.exc.idx]=NA
    cornwall.body.pix[cornwall.tip.inc.idx]=NA
    cornwall.both<-list(cornwall.tip.pix,cornwall.body.pix)
}
plot.ma.gbirl<-function(n,q,dir){
    if(!dir.exists(dir)){
        dir.create(dir)
    }
    px.mean.gb<-get.weight.gbirl(n,q,ancestry.spdf)
    plot.px<-get(paste0('gb.raster',n))
    plot.px.gbidx<-get(paste0('gb.raster',n,'.idx'))
    plot.px.noiidx<-get(paste0('noi.raster',n,'.idx'))
    plot.px.irlidx<-get(paste0('irl.raster',n,'.idx'))
   
    for(i in 1:nrow(px.mean.gb)){
        pdf(paste0(dir,'/GB_regions_',colnames(ancestry.spdf@data)[i],'_',n,'_',q,'_entropy_v2.pdf'))
        plot.px[plot.px.gbidx]<-px.mean.gb[i,]
        plot.px[plot.px.noiidx]<-noi.mean.ancestry[i]
        plot.px[plot.px.irlidx]<-irl.mean.ancestry[i]
        plot(plot.px,main=colnames(ancestry.spdf@data)[i],col=blue2green2red(100))
        plot(gbirl,add=T)
        dev.off()
    }
}

library(sp)
library(fastcluster)
library(raster)
library(colorRamps)
library(parallel)

#color=grDevices::colors()[grep('gr(a|e)y',grDevices::colors(),invert=T)]

#gb.color<-readRDS('23_regions_color1.rds')

if(!exists('England_Wales')){
    gbirl<-readRDS('GB_IRELAND.rds')
    panel.box<-bbox(gbirl)
    England_Wales<-readRDS('UK_main_cities_England_Wales.rds')
    England_Wales<-spTransform(England_Wales,proj4string(gbirl))
    Glasgow<-readRDS('Scotland/Glasgow.rds')
    Dundee<-readRDS('Scotland/Dundee.rds')
    Edinburgh<-readRDS('Scotland/Edinburgh.rds')
    Aberdeen<-readRDS('Scotland/Aberdeen.rds')
}

if(!exists('GB.counties')){
    gb.counties<-readRDS('../map/GBR_adm2.rds')
    gb.counties<-spTransform(gb.counties,CRS('+init=epsg:4326'))
    
    london.area<-c('London Borough','London Borough (city)','London Borough (royal)')

    GB.counties.ids<-!gb.counties$TYPE_2 %in% london.area
    
    Gb.counties<-gb.counties[GB.counties.ids,]
}

make.entropy<-function(vec,eps=1e-20){ ###for each sample, we can calculate its' entropy to measure how admixed this sample is
    vec[vec<eps]<-eps
    vec.entropy<--sum(vec*log(vec),na.rm=T)
    return(vec.entropy)
}

make.entropy.sub<-function(vec,sub.idx=1:23,eps=1e-20){
    vec.new<-vec[sub.idx]
    vec.new<-vec.new/sum(vec.new)
    vec.new.entropy<--sum(vec.new*log(vec.new),na.rm=T)
    return(vec.new.entropy)
}

my.bootstrap<-function(vec,rep=1000){
    n<-length(vec)
    mean.bootstrap<-sapply(1:rep,function(i){
        vec.res<-sample(vec,replace=T)
        return(mean(vec.res))
    })
    return(quantile(mean.bootstrap,c(0.05,0.95)))
}

if(!exists('v2.cl.srb.spdf')){
    long_lat.dt<-read.table('../long_lat_ukbiobank.tsv',as.is=T,header=T,sep="\t")[,1:2]

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

    esw.cl.entropy<-apply(esw.cl,1,make.entropy)
    esw.cl.gbirl.entropy<-apply(esw.cl,1,make.entropy.sub)
    noi.cl.entropy<-apply(noi.cl,1,make.entropy)
    noi.cl.gbirl.entropy<-apply(noi.cl,1,make.entropy.sub)
    noi.cl.mean.entropy<-mean(noi.cl.entropy)
    noi.cl.mean.gbirl.entropy<-mean(noi.cl.gbirl.entropy)

    irl.cl.entropy<-apply(irl.cl,1,make.entropy)
    irl.cl.gbirl.entropy<-apply(irl.cl,1,make.entropy.sub)
    irl.cl.mean.entropy<-mean(irl.cl.entropy)
    irl.cl.mean.gbirl.entropy<-mean(irl.cl.gbirl.entropy)
    irl.mean.entropy<-c(irl.cl.mean.entropy,irl.cl.mean.gbirl.entropy)
    noi.mean.entropy<-c(noi.cl.mean.entropy,noi.cl.mean.gbirl.entropy)

    colnames(ESW.cord)<-c('long','lat')
    coordinates(ESW.cord)<- ~ long+lat
    proj4string(ESW.cord)<-proj4string(gbirl)
    
    noi.mean.ancestry<-colMeans(noi.cl)
    irl.mean.ancestry<-colMeans(irl.cl)

    v2.cl.srb.spdf<-SpatialPointsDataFrame(ESW.cord,data=data.frame(esw.cl))
    v2.cl.entropy.spdf<-SpatialPointsDataFrame(ESW.cord,data=data.frame(all.entropy=esw.cl.entropy,gb.entropy=esw.cl.gbirl.entropy))
}

if(!exists('county.spdf.list')){
    v2InCounties<-v2.cl.srb.spdf%over%gb.counties

    v2.counties.spdf.nonNA<-!is.na(v2InCounties$NAME_2)

    v2InCounties.nonNA<-v2InCounties[v2.counties.spdf.nonNA,]
    v2.counties.spdf<-v2.cl.srb.spdf[v2.counties.spdf.nonNA,]
    v2.counties.entropy.spdf<-v2.cl.entropy.spdf[v2.counties.spdf.nonNA,]
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
        attr(county.spdf,'county.sample.entropy')<-apply(county.spdf@data,1,make.entropy)
        attr(county.spdf,'county.sample.gb.entropy')<-apply(county.spdf@data,1,make.entropy.sub)
        attr(county.spdf,'county.mean.entropy')<-mean(attr(county.spdf,'county.sample.entropy'))
        attr(county.spdf,'county.mean.gb.entropy')<-mean(attr(county.spdf,'county.sample.gb.entropy'))
        attr(county.spdf,'county.mean.entropy.CI')<-my.bootstrap(attr(county.spdf,'county.sample.entropy'))
        attr(county.spdf,'county.mean.gb.entropy.CI')<-my.bootstrap(attr(county.spdf,'county.sample.gb.entropy'))
        return(county.spdf)
    },mc.cores=12)
    county.spdf.list<-county.spdf.list[sapply(county.spdf.list,function(x)!is.null(x))]
    names(county.spdf.list)<-sapply(county.spdf.list,function(x)attr(x,'county.name'))
}

if(!exists('city.spdf.list')){
    attr(v2.cl.srb.spdf,'all.entropy')<-apply(v2.cl.srb.spdf@data,1,make.entropy)
    attr(v2.cl.srb.spdf,'gb.entropy')<-apply(v2.cl.srb.spdf@data,1,make.entropy.sub)
    
    v2InCities<-v2.cl.srb.spdf%over% England_Wales
    v2InGlasgow<-v2.cl.srb.spdf%over% Glasgow
    v2InDundee<-v2.cl.srb.spdf%over% Dundee
    v2InEdinburgh<-v2.cl.srb.spdf%over% Edinburgh
    v2InAberdeen<-v2.cl.srb.spdf%over% Aberdeen 

    v2.cities.spdf.nonNA<-!is.na(v2InCities$tcity15nm)
    v2.Glasgow.spdf.nonNA<-!is.na(v2InGlasgow[[1]])
    v2.Dundee.spdf.nonNA<-!is.na(v2InDundee[[1]])
    v2.Edinburgh.spdf.nonNA<-!is.na(v2InEdinburgh[[1]])
    v2.Aberdeen.spdf.nonNA<-!is.na(v2InAberdeen[[1]])
   
    Scotland.cities<-c('Glasgow','Dundee','Edinburgh','Aberdeen')
    Scotland.cities.spdf.list<-lapply(Scotland.cities,function(scity){
        scity.nonNA<-get(paste0('v2.',scity,'.spdf.nonNA'))
        scity.sp.size<-sum(scity.nonNA)
        scity.spdf<-v2.cl.srb.spdf[scity.nonNA,]
        attr(scity.spdf,'city.name')<-scity
        attr(scity.spdf,'city.sample.size')<-scity.sp.size
        attr(scity.spdf,'city.mean.ancestry')<-colMeans(scity.spdf@data[,1:23])
        attr(scity.spdf,'city.meangb.percent')<-c(mean(rowSums(scity.spdf@data[,1:23])),mean(rowSums(scity.spdf@data[,24:127])))
        attr(scity.spdf,'city.sample.entropy')<-apply(scity.spdf@data,1,make.entropy)
        attr(scity.spdf,'city.sample.gb.entropy')<-apply(scity.spdf@data,1,make.entropy.sub)
        attr(scity.spdf,'city.mean.entropy')<-mean(attr(scity.spdf,'city.sample.entropy'))
        attr(scity.spdf,'city.mean.gb.entropy')<-mean(attr(scity.spdf,'city.sample.gb.entropy'))
        attr(scity.spdf,'city.mean.entropy.CI')<-my.bootstrap(attr(scity.spdf,'city.sample.entropy'))
        attr(scity.spdf,'city.mean.gb.entropy.CI')<-my.bootstrap(attr(scity.spdf,'city.sample.gb.entropy'))
        return(scity.spdf)
    })
    names(Scotland.cities.spdf.list)<-Scotland.cities
    v2InCities.nonNA<-v2InCities[v2.cities.spdf.nonNA,]
    v2.cl.srb.city.spdf<-v2.cl.srb.spdf[v2.cities.spdf.nonNA,]

    city.sp.sizes<-table(v2InCities.nonNA[[3]])
    city.spdf.list<-mclapply(names(city.sp.sizes),function(city.name){
        if(city.sp.sizes[city.name]==0){
            return(NULL)
        }
        city.spdf<-v2.cl.srb.city.spdf[v2InCities.nonNA[[3]]==city.name,]
        attr(city.spdf,'city.name')<-city.name
        attr(city.spdf,'city.sample.size')<-city.sp.sizes[city.name]
        attr(city.spdf,'city.mean.ancestry')<-colMeans(city.spdf@data[,1:23])
        attr(city.spdf,'city.meangb.percent')<-c(mean(rowSums(city.spdf@data[,1:23])),mean(rowSums(city.spdf@data[,24:127])))
        attr(city.spdf,'city.sample.entropy')<-apply(city.spdf@data,1,make.entropy)
        attr(city.spdf,'city.sample.gb.entropy')<-apply(city.spdf@data,1,make.entropy.sub)
        attr(city.spdf,'city.mean.entropy')<-mean(attr(city.spdf,'city.sample.entropy'))
        attr(city.spdf,'city.mean.gb.entropy')<-mean(attr(city.spdf,'city.sample.gb.entropy'))
        attr(city.spdf,'city.mean.entropy.CI')<-my.bootstrap(attr(city.spdf,'city.sample.entropy'))
        attr(city.spdf,'city.mean.gb.entropy.CI')<-my.bootstrap(attr(city.spdf,'city.sample.gb.entropy'))
        return(city.spdf)
    },mc.cores=12)
    city.spdf.list<-city.spdf.list[sapply(city.spdf.list,function(x)!is.null(x))]
    names(city.spdf.list)<-sapply(city.spdf.list,function(x)attr(x,'city.name'))
    Scotland.attrs<-lapply(Scotland.cities.spdf.list,attributes)
    EW.attrs<-lapply(city.spdf.list,attributes)
    all.attrs<-EW.attrs
    for(i in 1:length(Scotland.attrs)){
        all.attrs[[length(EW.attrs)+i]]<-Scotland.attrs[[i]]
    }
    all.attrs<-all.attrs[sapply(all.attrs,function(x)!is.null(x))]
    all.attrs.big.cities<-all.attrs[which(sapply(all.attrs,function(x)x$city.sample.size>5000))]

    names(all.attrs.big.cities)<-sapply(all.attrs.big.cities,function(x)x$city.name)
}

sample.normalise<-function(mt){
    nmt<-mt-matrix(rep(colMeans(mt),each=nrow(mt)),nrow=nrow(mt))
    return(nmt)
}

get.sim<-function(x){
    nx<-sample.normalise(x)
    dist<-1-nx%*%t(nx)
    hc<-hclust(as.dist(dist))
    od<-hc$order
    permut.x<-x[od,]
    return(permut.x)
}

get.sim.idx<-function(x){
    nx<-sample.normalise(x)
    dist<-1-nx%*%t(nx)
    hc<-hclust(as.dist(dist))
    od<-hc$order
    return(return(od))
}

my.barplot<-function(mat,rec.cord.bot,rec.cord.top,ifGB=T,bar.width=1.0){
    xb0<-rec.cord.bot[1]
    yb0<-rec.cord.bot[2]
    xt1<-rec.cord.top[1]
    yt1<-rec.cord.top[2]

    length<-xt1-xb0
    width<-yt1-yb0
    if(ifGB){
        sub<-1:23
        mat<-mat[,1:23]
        mat.rowsum<-rowSums(mat)
        mat<-mat/mat.rowsum
    }
    else{
        sub<-1:ncol(mat)
    }
    mat<-t(apply(mat,1,cumsum))
    n.sub<-length(sub)
    n.sap<-nrow(mat)
    x.seqb<-rep(0:(n.sap-1),each=n.sub)
    x.seqt<-rep(1:n.sap,each=n.sub)
    y.seqb<-as.vector(t(cbind(rep(0,n.sap),mat[,-n.sub])))
    y.seqt<-as.vector(t(mat))
    
    xb<-xb0+(xt1-xb0)*x.seqb/(bar.width*n.sap)
    xt<-xb0+(xt1-xb0)*x.seqt/(bar.width*n.sap)
    yb<-yb0+(yt1-yb0)*y.seqb
    yt<-yb0+(yt1-yb0)*y.seqt
    
    mt<-cbind(xb,yb,xt,yt)
    return(mt)
}

small.str.plot<-function(rec.cord.bot,rec.cord.top,anc.mat,col=gb.color,ifGB){
    anc.mat<-as.matrix(anc.mat)
    mt<-my.barplot(anc.mat,rec.cord.bot,rec.cord.top,ifGB=ifGB)
    rect(mt[,1],mt[,2],mt[,3],mt[,4],col=col,border=F)
}

plot.each.place<-function(spdf.list,nm,place='county',dir){
    spdf.tmp<-spdf.list[[nm]]
    anc.mat<-as.matrix(spdf.tmp@data)
    anc.od<-get.sim.idx(anc.mat)
    anc.mat<-anc.mat[anc.od,]
    pdf(paste0(dir,'/',nm,'.pdf'))
    mt<-my.barplot(anc.mat,c(0,0),c(3,1))
    tmp.sample.entropy<-attr(spdf.tmp,paste0(place,'.sample.gb.entropy'))
    tmp.sample.entropy<-tmp.sample.entropy[anc.od]
    par(oma=rep(0,4),mar=rep(0,4))
    layout(matrix(1:2,ncol=1),height=c(1,5))
    x0.entro.p<-matrix(mt[,1],nrow=23)[1,]
    x1.entro.p<-matrix(mt[,3],nrow=23)[1,]
    x.entro.p<-(x0.entro.p+x1.entro.p)/2
    y.entro.p<-tmp.sample.entropy
    plot(x.entro.p,y.entro.p,col='grey',pch='.',xlim=c(-0.01,3.01),xaxs='i',yaxs='i',ann=F,axes=F)
    plot(0,0,type='n',xlim=c(-0.01,3.01),ylim=c(-0.01,1.01),xaxs='i',yaxs='i',ann=F,axes=F) 
    rect(mt[,1],mt[,2],mt[,3],mt[,4],col=gb.color,border=F)
    dev.off()
}

plot.combine.place<-function(spdf.list,nms,place='county',dir){
    #spdf.tmp<-spdf.list[[nm]]
    anc.mat<-as.matrix(do.call(rbind,lapply(spdf.list[nms],function(x)x@data)))
    #anc.mat<-as.matrix(spdf.tmp@data)
    anc.od<-get.sim.idx(anc.mat)
    anc.mat<-anc.mat[anc.od,]
    pdf(paste0(dir,'/',paste(nms,collapse='_'),'.pdf'))
    mt<-my.barplot(anc.mat,c(0,0),c(3,1))
    #tmp.sample.entropy<-attr(spdf.tmp,paste0(place,'.sample.gb.entropy'))
    tmp.sample.entropy<-do.call(c,lapply(spdf.list[nms],function(x)attr(x,paste0(place,'.sample.gb.entropy'))))
    tmp.sample.entropy<-tmp.sample.entropy[anc.od]
    par(oma=rep(0,4),mar=rep(0,4))
    layout(matrix(1:2,ncol=1),height=c(1,5))
    x0.entro.p<-matrix(mt[,1],nrow=23)[1,]
    x1.entro.p<-matrix(mt[,3],nrow=23)[1,]
    x.entro.p<-(x0.entro.p+x1.entro.p)/2
    y.entro.p<-tmp.sample.entropy
    plot(x.entro.p,y.entro.p,col='grey',pch='.',xlim=c(-0.01,3.01),xaxs='i',yaxs='i',ann=F,axes=F)
    plot(0,0,type='n',xlim=c(-0.01,3.01),ylim=c(-0.01,1.01),xaxs='i',yaxs='i',ann=F,axes=F) 
    rect(mt[,1],mt[,2],mt[,3],mt[,4],col=gb.color,border=F)
    dev.off()
}

run.each.place<-function(dir1,dir2){
    if(!dir.exists(dir1)){
        dir.create(dir1)
    }
    if(!dir.exists(dir2)){
        dir.create(dir2)
    }
    #mclapply(names(Scotland.cities.spdf.list),plot.each.place,spdf.list=Scotland.cities.spdf.list,place='city',dir=dir1,mc.cores=4)
    #mclapply(names(city.spdf.list),plot.each.place,spdf.list=city.spdf.list,place='city',dir=dir1,mc.cores=12)
    mclapply(names(county.spdf.list),plot.each.place,spdf.list=county.spdf.list,place='county',dir=dir2,mc.cores=12)

}

make.raster.gbirl<-function(n=250){
    panel.box<-rbind(c(-14,5),c(47,63))
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
}

get.weight.gbirl<-function(n,q,spdf){#at the moment,try constant
    md<-make.raster.gbirl(n)
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

plot.ma.single.entropy<-function(n,q,entropy.spdf,dir){
    if(!dir.exists(dir)){
        dir.create(dir)
    }
    px.mean.gb<-get.weight.gbirl(n,q,entropy.spdf)
    plot.px<-get(paste0('gb.raster',n))
    plot.px.gbidx<-get(paste0('gb.raster',n,'.idx'))
    plot.px.noiidx<-get(paste0('noi.raster',n,'.idx'))
    plot.px.irlidx<-get(paste0('irl.raster',n,'.idx'))
    pdf(paste0(dir,'/UK_entropy_average_entropy_',n,'_v2.pdf'))
    par(bty='n',mar=c(2,4,2,5))
    plot.px[plot.px.gbidx]<-px.mean.gb[2,]
    plot.px[plot.px.noiidx]<-noi.mean.entropy[2]
    plot.px[plot.px.irlidx]<-irl.mean.entropy[2]
    plot(plot.px,main=NA,axes=F,ann=F,col=blue2green2red(100))
    #plot(gbirl,add=T)
    plot(England_Wales[England_Wales@data$tcity15nm%in%c('London','Liverpool','Birmingham','Bristol','Manchester'),],add=T)
    plot(Edinburgh,add=T)
    plot(Glasgow,add=T)
    plot(Aberdeen,add=T)

    plot(gb.counties[gb.counties@data$NAME_2%in%c('Isle of Wight','Norfolk','Suffolk','Gwynedd','Anglesey','Carmarthenshire','West Yorkshire','Hartlepool','Stockton-on-Tees','Middlesbrough','Plymouth','Oxfordshire','East Riding of Yorkshire','North East Lincolnshire','Fife','Angus','Perthshire and Kinross','Ceredigion','Northumberland','Cumbria','Cardiff'),],add=T)
    #West Yorkshire low entropy
    #Liverpool high entropy
    #Isle of Wight
    #London 
    #Gloucester
    #Norfolk
    #Edinburgh
    #Angus
    #Aberdeen city
    #Gwynedd
    #Carmarthenshire
    #Hartlepool/Stockton-on-Tees/Middlesbrough/
    #Plymouth
    #Oxfordshire
    #East Riding of Yorkshire
    #North East Lincolnshire
    dev.off()
}

plot.ma.entropy<-function(n,q,entropy.spdf,dir){
    if(!dir.exists(dir)){
        dir.create(dir)
    }
    px.mean.gb<-get.weight.gbirl(n,q,entropy.spdf)
    plot.px<-get(paste0('gb.raster',n))
    plot.px.gbidx<-get(paste0('gb.raster',n,'.idx'))
    plot.px.noiidx<-get(paste0('noi.raster',n,'.idx'))
    plot.px.irlidx<-get(paste0('irl.raster',n,'.idx'))

    Birmingham.coords<-c(-15,-11,52,53)
    Bristol.coords<-c(-5.5,-1.5,49,50)
    Cardiff.coords<-c(-10,-6,50,51)
    Leeds.coords<-c(1.5,5.5,54.5,55.5)
    Liverpool.coords<-c(-14,-10,54,55)
    London.coords<-c(2,6,49.5,50.5)
    Manchester.coords<-c(2,6,53,54)
    Edinburgh.coords<-c(0,4,56,57)
    Glasgow.coords<-c(-11,-7,55.5,56.5)

    for(i in nrow(px.mean.gb):nrow(px.mean.gb)){
        tiff(paste0(dir,'/',colnames(entropy.spdf@data)[i],'_average_entropy_',n,'_v2.tiff'),height=2400,width=2400)
        #pdf(paste0(dir,'/',colnames(entropy.spdf@data)[i],'_average_entropy_',n,'_v2.pdf'))
        par(bty='n',mar=c(2,4,2,5))
        plot.px[plot.px.gbidx]<-px.mean.gb[i,]
        plot.px[plot.px.noiidx]<-noi.mean.entropy[i]
        plot.px[plot.px.irlidx]<-irl.mean.entropy[i]
        #plot(plot.px,main=colnames(entropy.spdf@data)[i],col=blue2green2red(100))
        plot(plot.px,main=NA,axes=F,ann=F,col=blue2green2red(100))
        plot(gbirl,add=T)
        #plot(England_Wales[England_Wales@data$tcity15nm%in%names(all.attrs.big.cities),],add=T,col='darkgrey')
        #plot(Glasgow,col='darkgrey',add=T)
        #plot(Edinburgh,col='darkgrey',add=T)
        my.gb.cols<-rainbow(1000)[round(seq(1,1000,length=23))]
        cm.reorder.mats<-mclapply(names(all.attrs.big.cities),function(nm){
            if(!exists(paste0(nm,'.coords'))){
                return(NULL)
            }
            if(nm %in% England_Wales@data$tcity15nm){
                cm.reorder<-get.sim(as.matrix(city.spdf.list[[nm]]@data)[,1:23])
                return(cm.reorder)
            }
            if(nm=='Edinburgh'||nm=='Glasgow'){
                cm.reorder<-get.sim(as.matrix(Scotland.cities.spdf.list[[nm]]@data)[,1:23])
                return(cm.reorder)
            }
        },mc.cores=3)
        cm.null<-sapply(cm.reorder.mats,function(x)!is.null(x))
        plot.cities<-names(all.attrs.big.cities)[cm.null]
        cm.reorder.mats<-cm.reorder.mats[cm.null]
        names(cm.reorder.mats)<-plot.cities
        for(nm in plot.cities){
            coords<-get(paste0(nm,'.coords'))
            small.str.plot(coords[c(1,3)],coords[c(2,4)],cm.reorder.mats[[nm]],gb.color,ifGB=T)
            text.coords<-c(mean(coords[1:2]),coords[4]+0.2)
            text(text.coords[1],text.coords[2],nm,cex=0.8)
            if(nm%in%England_Wales$tcity15nm){
                city.ct<-coordinates(England_Wales[England_Wales$tcity15nm==nm,])
            }
            if(nm=='Edinburgh'||nm=='Glasgow'){
                city.ct<-colMeans(coordinates(get(nm))) 
            }
            if(nm %in% c('Glasgow','Liverpool','Birmingham','Cardiff','Bristol')){
                lines(c(coords[2],city.ct[1]),c(coords[4],city.ct[2]))
            }
            else{
                lines(c(coords[1],city.ct[1]),c(coords[3],city.ct[2]))
            }
        }
        dev.off()
    }
}

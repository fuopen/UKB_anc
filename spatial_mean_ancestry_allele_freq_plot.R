library(parallel)
library(fastcluster)
library(sp)
library(raster)
library(colorRamps)
library(igraph)
library(gplots)
library(Rfast)

load.map<-function(){
    gbirl<<-readRDS('data/GB_IRELAND.rds')#GB+IRELAND map data
    gb.map<<-readRDS('data/new_boundary.rds')#GB county boundary map data
    irl<<-gbirl[5,]
    my.box.gb<<-bbox(gbirl)
}

make.entropy.sub<-function(vec,sub.idx=1:23,eps=1e-20){
    vec.new<-vec[sub.idx]
    vec.new<-vec.new/sum(vec.new)
    vec.new[vec.new<eps]<-eps
    vec.entropy<--sum(vec.new*log(vec.new),na.rm=T)
    return(vec.entropy)
}
if(!exists('mp')){
    mp<-readRDS('mapping.rds')#individual level data
    withdraw.ids<-readRDS('withdraw_list.rds')#individual level data
    mp$withdraw<-mp[[2]]%in%withdraw.ids[[1]]
    v2.eth<-readRDS('biobank_v2_eth.background.rds')#individual level data
    BI<-c('British','Irish')
    norirl<-'Northern_Ireland'
    irl<-'Republic_of_Ireland'

    mp.ids<-match(mp[[1]],v2.eth[[1]])
    v2.eth<-v2.eth[mp.ids,]
    NI.ids<-(!is.na(v2.eth[[2]])) & v2.eth[[2]]==norirl & v2.eth[[4]]%in%BI
    IRL.ids<-(!is.na(v2.eth[[2]])) & v2.eth[[2]]==irl & v2.eth[[4]]%in%BI

}

if(!exists('ac.new')){
    ac<-readRDS('v2_487409.rds')#ACs data, individual level data

    ids<-readRDS('self_BI_487409.rds')#Self-reported birth country, individual level data

    ac.new<-ac[ids,]
    noi.cl<-ac[NI.ids,]
    irl.cl<-ac[IRL.ids,]

    ac.new$E_italians<-ac.new$E_italians+ac.new$E_Italy_3+ac.new$E_North_Italian
    ac.new$E_Italy_3<-NULL
    ac.new$E_North_Italian<-NULL
    ac.new$E_French<-ac.new$E_French+ac.new$E_French_Basque
    ac.new$E_French_Basque<-NULL
    ac.new$E_german<-ac.new$E_german+ac.new$E_DutchlikeGerman+ac.new$E_BelgianlikeGerman
    ac.new$E_DutchlikeGerman<-NULL
    ac.new$E_BelgianlikeGerman<-NULL
    ac.new$E_Spain_3<-ac.new$E_Spain_3+ac.new$E_PortlikeSpain
    ac.new$E_PortlikeSpain<-NULL

    noi.cl$E_italians<-noi.cl$E_italians+noi.cl$E_Italy_3+noi.cl$E_North_Italian
    noi.cl$E_Italy_3<-NULL
    noi.cl$E_North_Italian<-NULL
    noi.cl$E_French<-noi.cl$E_French+noi.cl$E_French_Basque
    noi.cl$E_French_Basque<-NULL
    noi.cl$E_german<-noi.cl$E_german+noi.cl$E_DutchlikeGerman+noi.cl$E_BelgianlikeGerman
    noi.cl$E_DutchlikeGerman<-NULL
    noi.cl$E_BelgianlikeGerman<-NULL
    noi.cl$E_Spain_3<-noi.cl$E_Spain_3+noi.cl$E_PortlikeSpain
    noi.cl$E_PortlikeSpain<-NULL

    irl.cl$E_italians<-irl.cl$E_italians+irl.cl$E_Italy_3+irl.cl$E_North_Italian
    irl.cl$E_Italy_3<-NULL
    irl.cl$E_North_Italian<-NULL
    irl.cl$E_French<-irl.cl$E_French+irl.cl$E_French_Basque
    irl.cl$E_French_Basque<-NULL
    irl.cl$E_german<-irl.cl$E_german+irl.cl$E_DutchlikeGerman+irl.cl$E_BelgianlikeGerman
    irl.cl$E_DutchlikeGerman<-NULL
    irl.cl$E_BelgianlikeGerman<-NULL
    irl.cl$E_Spain_3<-irl.cl$E_Spain_3+irl.cl$E_PortlikeSpain
    irl.cl$E_PortlikeSpain<-NULL
    
    ac.meananc<-sapply(ac.new,function(x)mean(x>0.01))

    ac.new<-ac.new[,ac.meananc>0.01]

    noi.cl<-noi.cl[,ac.meananc>0.01]
    irl.cl<-irl.cl[,ac.meananc>0.01]

    noi.mean.ancestry<-colMeans(noi.cl)
    irl.mean.ancestry<-colMeans(irl.cl)

    ac.entropy.gb<-apply(as.matrix(ac.new[,1:23]),1,make.entropy.sub)

    spdf<-readRDS('v2_self_BI_487409.rds')#individual level data: spatial dataframe for ACs
    spdf@data<-ac.new
	#spdf<-readRDS('gwas_snp_regional_plot_af.rds')
}


if(!exists('gb.bound')){
    gb.bound<-readRDS('data/New_GB_boundaries.rds')#GB boundary information
    load.map()
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
    px.mean.gb<-get.weight.gbirl(n,q,spdf)
    plot.px<-get(paste0('gb.raster',n))
    plot.px.gbidx<-get(paste0('gb.raster',n,'.idx'))
    plot.px.noiidx<-get(paste0('noi.raster',n,'.idx'))
    plot.px.irlidx<-get(paste0('irl.raster',n,'.idx'))

    for(i in 1:nrow(px.mean.gb)){
        pdf(paste0(dir,'/GB_regions_',colnames(spdf@data)[i],'_',n,'_',q,'.pdf'))
        par(bty='n')
        plot.px[plot.px.gbidx]<-px.mean.gb[i,]
        plot.px[plot.px.noiidx]<-noi.mean.ancestry[i]
        plot.px[plot.px.irlidx]<-irl.mean.ancestry[i]
        plot(plot.px,main=colnames(spdf@data)[i],col=blue2green2red(100),ann=F,axes=F)
        #plot(gbirl,add=T,border=NA)
        dev.off()
    }
}

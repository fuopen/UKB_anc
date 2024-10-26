library(parallel)
library(fastcluster)
library(sp)
library(raster)
library(colorRamps)
library(igraph)
library(gplots)
library(Rfast)

if(!exists('gb.raster')){
    gb.raster<-readRDS('data/rastered_23GB_regions.rds')
}

gb.raster.mat<-function(){
    grm<-sapply(gb.raster,function(x)x@data@values)
}

if(!exists('grm')){
    grm<-gb.raster.mat()
    grm.nv<-apply(grm,2,function(x)which(!is.na(x)))
}

if(!exists('IRL.ids')){
    mp<-readRDS('mapping.rds')
    withdraw.ids<-readRDS('withdraw_list.rds')
    mp$withdraw<-mp[[2]]%in%withdraw.ids[[1]]
    v2.eth<-readRDS('biobank_v2_eth.background.rds')
    BI<-c('British','Irish')
    norirl<-'Northern_Ireland'
    irl<-'Republic_of_Ireland'

    mp.ids<-match(mp[[1]],v2.eth[[1]])
    v2.eth<-v2.eth[mp.ids,]
    NI.ids<-(!is.na(v2.eth[[2]])) & v2.eth[[2]]==norirl & v2.eth[[4]]%in%BI
    IRL.ids<-(!is.na(v2.eth[[2]])) & v2.eth[[2]]==irl & v2.eth[[4]]%in%BI

}

if(!exists('ac.new')){
    ac<-readRDS('v2_487409.rds')

    ids<-readRDS('self_BI_487409.rds')

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


    spdf<-readRDS('v2_self_BI_487409.rds')
    spdf@data<-ac.new
    counties<-readRDS('self_BI_counties_487409.rds')
}

if(!exists('spdf.birth.region')){
    spdf.birth.region<-readRDS('426879_ind_birth_place.rds')
}

test.correct<-function(test.dir='test_correct/'){
    if(!dir.exists(test.dir)){
        dir.create(test.dir)
    }
    spdf.gb<-names(spdf)[1:23]
    plot.rg<-function(i){
        spdf.tmp<-spdf[spdf.birth.region==spdf.gb[i],]
        pdf(paste0(test.dir,spdf.gb[i],'.pdf'))
        plot(gb.raster[[spdf.gb[i]]])
        plot(spdf.tmp,pch='.',add=T)
        dev.off()
    }
    mclapply(1:23,plot.rg,mc.cores=10)
}

if(!exists('region.anc')){
    region.anc<-by(ac.new,spdf.birth.region,function(x)x)
    region.anc$G_Noi<-rbind(region.anc$G_Noi,noi.cl)
    region.anc$G_Ireland<-irl.cl
}

make.raster.gbirl<-function(n=250,correct.anglia=F){
    panel.box<-rbind(c(-11,2),c(50,61.5))
    cord.spdf<-data.frame(coordinates(spdf))
    coordinates(cord.spdf)<- ~birth.lon+birth.lat
    projection(cord.spdf)<-projection(spdf)
    eval(parse(text=paste0('base.raster',n,'<-raster(nr=n,ncol=n,xmn=panel.box[1,1],xmx=panel.box[1,2],ymn=panel.box[2,1],ymx=panel.box[2,2],crs=CRS("+init=epsg:4326"))')))
    eval(parse(text=paste0('base.raster',n,'[]<-1:ncell(base.raster',n,')')))
    spot.each<-function(i){
        eval(parse(text=paste0('spdf.raster',n,'<-rasterize(cord.spdf[',i,'],base.raster',n,')')))
        eval(parse(text=paste0('spdf.raster',n,'<-mask(base.raster',n,',spdf.raster',n,')')))
        eval(parse(text=paste0('spdf.raster',n,'.idx<-which(!is.na(getValues(spdf.raster',n,')))')))
        eval(parse(text=paste0('sap<-sapply(grm.nv,function(x)sum(spdf.raster',n,'.idx%in%x))')))
        eval(parse(text=paste0('lk<-names(grm.nv)[which.max(sap)]')))
    }
    if(!correct.anglia){
        me<-mclapply(1:length(cord.spdf),spot.each,mc.cores=16)
    }
    else{
        spot.each.ag<-function(i){
            eval(parse(text=paste0('spdf.raster',n,'<-rasterize(cord.spdf[',i,'],base.raster',n,')')))
            eval(parse(text=paste0('spdf.raster',n,'<-mask(base.raster',n,',spdf.raster',n,')')))
            eval(parse(text=paste0('spdf.raster',n,'.idx<-which(!is.na(getValues(spdf.raster',n,')))')))
            eval(parse(text=paste0('sap<-spdf.raster',n,'.idx%in%grm.nv[[1]]')))
            eval(parse(text=paste0('lk<-any(sap)')))
        }
        me<-mclapply(1:length(cord.spdf),spot.each.ag,mc.cores=16)
    }
}

if(!exists('reg.anc')){
    irish<-readRDS('Irish.rds')
    north.irish<-readRDS('North_Irish.rds')
    #reg.ancr<-by(Sp.data,SpInCounties$assign,function(x)x)
    #reg.anc<-lapply(reg.ancr,function(x)x[,1:22])
}

color=grDevices::colors()[grep('gr(a|e)y',grDevices::colors(),invert=T)]
if(!exists('new.colors')){
    gb.color<-readRDS('23_regions_color1.rds')
    names(gb.color)<-names(gb.raster)

    set.seed(20181223)
    other.colors<-color[!color%in% gb.color]
    other.colors<-sample(other.colors,ncol(spdf)-length(gb.color),replace=T)
    names(other.colors)<-names(spdf)[!names(spdf)%in%names(gb.color)]
    new.colors<-c(gb.color,other.colors)
    new.colors<-new.colors[names(spdf)]

    new.colors[['E_norwegian']]<-'chartreuse3'
    new.colors[['E_Switzerland_1']]<-'chartreuse4'
    new.colors[['E_Portugal']]<-'darkolivegreen2'
    new.colors[['E_TSI']]<-'gold'
    new.colors[['W_syrian']]<-'lightpink3'
    new.colors[['E_Denmark']]<-'springgreen'
}

official.names<-c('Anglia','South Central England','Central England','North Central England','North Yorkshire','Northumberland','South East England','Merseyside','South England','South East Wales','North East Scotland','South West Wales','Devon','North West Wales','South Yorkshire','Cornwall Tip','Orkney','Lincolnshire','Cumbria','North West Scotland','North Ireland and South Scotland','Ireland','Cornwall')

names(official.names)<-names(gb.color)

official.others<-c('France','French Basque','Sardinian','North Italian','Italy1','Adygei','Russia','Ukraine','Bulgaria','Croatia','Cyprus','DutchlikeGerman','Germany', 'Greece', 'Lithuania', 'Mordovia', 'Norway', 'Poland', 'Sicily','Spain', 'Switzerland','Yugoslavia','Czech Republic','Romania', 'Hungary', 'Austria', 'BelgianlikeGerman', 'Portugal', 'PortlikeSpain', 'Belgium', 'Netherlands', 'Finland', 'Italy2', 'Tuscany', 'Biaka Pygmies', 'Mandenka', 'Mbuti Pygmies', 'LWK', 'Bedouin', 'Palestine', 'Hazara', 'Mozabite', 'Balochi', 'San', 'India', 'SriLanka', 'Pathan', 'Kalash', 'Burusho', 'Sindhi', 'Yizu', 'Miaozu', 'Daur', 'Hezhen', 'Xibo', 'Uygur', 'Dai', 'She', 'Tu', 'Yakut', 'Papuan', 'Melanesian', 'Lahu', 'Naxi', 'Cambodia', 'Pima', 'Colombia', 'Karitiana', 'Nama', 'GIH', 'MKK', 'Chechen', 'Chuvash', 'Syria', 'East Turkey', 'South Turkey', 'AMHARA', 'ANUAK', 'ARI', 'Somali', 'GUMUZ', 'amaXhosa', 'ColouredWellington', 'Karretjie', 'GuiGhanaKgal', 'Khomani', 'Khwe', 'Yoruba', 'Juhoan', 'Maya', 'Japan', 'China', 'Denmark', 'Sweden', 'Algeria', 'Egypt', 'Morocco', 'Libya', 'Nepal', 'Myanmar', 'Bangladesh', 'Thailand', 'Malaysia', 'Philippines')

names(official.others)<-colnames(ac)[24:127]

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

if(!exists('neighbor')){
    neighbour<-readRDS("neighbour_region.rds")
}

barplot.eachgbreg2<-function(mat,reg.name,dir='GB_region_ancestry_fig_update') {#without normalisation
    if(!dir.exists(dir)){
        dir.create(dir)
    }
    reg.nn<-neighbour[[reg.name]]
    mat<-as.matrix(mat)
    mat.which<-apply(mat,1,which.max)
    mat.reg<-colnames(mat)[mat.which]
    mro<-names(sort(table(mat.which),decreasing=T))
    mat.reg.ord<-colnames(mat)[as.integer(mro)]
    #if(reg.name!=mat.reg.ord[1]){
    #    rg.idx<-mat.reg.ord==reg.name
    #    mat.reg.ord[rg.idx]=mat.reg.ord[1]
    #    mat.reg.ord[1]=reg.name
    #}
    rg.nidx<-which(!mat.reg.ord%in%reg.nn)
    mat.reg.ord<-c(reg.nn,mat.reg.ord[rg.nidx])
    me<-colnames(mat)
    mes<-me[!me%in%mat.reg.ord]
    mat<-mat[,c(mat.reg.ord,mes)]
    mat.reg.list<-by(mat,mat.reg,function(x){
        x<-as.matrix(x)
        if(nrow(x)>2){
            x.id<-hclust.vector(x)
            x<-x[rev(x.id$order),]
        }
        return(x)
    })
    mat<-do.call(rbind,mat.reg.list[mat.reg.ord])
    png(paste0(dir,'/',reg.name,'_strplot_',nrow(mat),'_order.png'),width=8.25,height=2.25,units='in',res=1200,pointsize=4)
    par(mar=c(0,0,0,0),oma=c(0,0,0,0))
    barplot(t(mat),col=new.colors[colnames(mat)],border=NA,space=0,ann=F,axes=F,xaxt='n',yaxt='n')
    dev.off()
}

barplot.gb<-function(dir2='GB_region_ancestry_fig_update'){
    if(!dir.exists(dir2)){
        dir.create(dir2)
    }
    mclapply(names(region.anc),function(regn){
        barplot.eachgbreg2(region.anc[[regn]],regn,dir=dir2)},mc.cores=6)
}


barplot.eachgbreg.meanmain<-function(mat,reg.name,ifbootstrap=F){
    mat<-as.matrix(mat)
    mat<-mat/rowSums(mat)
    mat.mean<-mean(mat[,reg.name],na.rm=T)
    if(ifbootstrap){
        mk<-mat[,reg.name]
        ns<-length(mk)
        bstr.means<-sapply(1:10000,function(i)mean(mk[sample(ns,replace=T)],na.rm=T))
        return(c(mat.mean,quantile(bstr.means,c(0.05,0.95))))
    }
    else{
        return(mat.mean)
    }
}

get.meanmain<-function(bs=T){
    gmm<-mclapply(names(reg.anc),function(x)barplot.eachgbreg.meanmain(reg.anc[[x]],x,bs),mc.cores=12)
    gmdt<-do.call(rbind,gmm)
}

raster.boundary<-function(n=200){
    panel.box<-rbind(c(-11,2),c(50,61.5))
    eval(parse(text=paste0('base.raster',n,'<-raster(nr=n,ncol=n,xmn=panel.box[1,1],xmx=panel.box[1,2],ymn=panel.box[2,1],ymx=panel.box[2,2],crs=CRS("+init=epsg:4326"))')))
    eval(parse(text=paste0('base.raster',n,'[]<-1:ncell(base.raster',n,')')))
    for(k in names(gb.bound)){
        eval(parse(text=paste0(k,'.raster',n,'<<-rasterize(gb.bound[["',k,'"]],base.raster',n,')')))
        eval(parse(text=paste0(k,'.raster',n,'<<-mask(base.raster',n,',',k,'.raster',n,')')))
        #eval(parse(text=paste0(k,'.raster',n,'.idx<<-which(!is.na(getValues(',k,'.raster',n,')))')))
    }
}

make.legend<-function(){
    par(mar=rep(0,4),oma=rep(0,4),xpd=T)
    plot(0,0,type='n',ann=F,axes=F)
    legend('left',legend=c(official.names,official.others[names(new.colors)[24:38]])[names(new.colors)],fill=new.colors,ncol=7,cex=0.1,bty='n')
}
make.entropy.sub<-function(vec,sub.idx=1:23,eps=1e-20){
    vec.new<-vec[sub.idx]
    vec.new<-vec.new/sum(vec.new)
    vec.new[vec.new<eps]<-eps
    vec.entropy<--sum(vec.new*log(vec.new),na.rm=T)
    return(vec.entropy)
}

find.mean.eachregion<-function(){
    rr.mean<-lapply(region.anc,function(x){v<-rowSums(x[,1:23]);x[v>0.95,]})
    rr.rpop.mean<-sapply(names(rr.mean),function(x)mean(rr.mean[[x]][,x]))
    rr.rpop.entropy<-sapply(rr.mean,function(x)mean(apply(x[,1:23],1,make.entropy.sub)))
    rout<-cbind(rr.rpop.mean,rr.rpop.entropy)
    colnames(rout)<-c('main_pop_mean','entropy')
    return(rout)
}

find.ancestry.from.birthandneibour<-function(){
	nb<-readRDS("neighbour_region.rds")
	rg<-region.anc
	fb<-function(pop){
		joint.pop<-nb[[pop]]
		rs<-rowSums(rg[[pop]][,joint.pop])
		ez<-sum(rs>=0.50)
	}
	ld<-lapply(names(rg),fb)
	names(ld)<-names(rg)
	lds<-sum(ld)
	ld.ratio<-sum(lds)/sum(sapply(rg,nrow))	
}
find.ancestry.from.birthandneibour2<-function(){
	nb<-readRDS("neighbour_region.rds")
	rg<-region.anc
	fb<-function(pop){
		joint.pop<-nb[[pop]]
		pop.anc<-rg[[pop]]
		pop.anc.single<-apply(pop.anc,1,function(x)any(x>=0.5))
		rs<-rowSums(pop.anc[pop.anc.single,joint.pop])
		ez<-c(sum(rs>=0.50),sum(pop.anc.single))
	}
	ld<-lapply(names(rg),fb)
	names(ld)<-names(rg)
	lds.r<-do.call(rbind,ld)
}

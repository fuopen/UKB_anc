library(parallel)

read.genmap<-function(chr){
    file<-paste0('gen_maps/genetic_map_chr',chr,'_combined_b37.txt')
    gen.map<-read.table(file,as.is=T,header=T)
    return(gen.map)
}

read.site<-function(chr){
    file<-paste0('imp_1kg_hap/snp_chr',chr,'.rds')
    site<-readRDS(file)
    site<-data.frame(chr=as.integer(site[,1]),pos=as.integer(site[,2]),rsid=site[,3],ref=site[,4],alt=site[,5],stringsAsFactors=F)
    return(site)
}

read.chr<-function(chr){
    gen.map<-read.genmap(chr)
    site<-read.site(chr)
    interp.x<-gen.map[,1]
    interp.y<-gen.map[,3]
    interp.xout<-site$pos
    interp.res<-approx(x=interp.x,y=interp.y,xout=interp.xout,rule=2)$y
    site$cm<-interp.res/100
    snps<-site[,c('rsid','chr','cm','pos','ref','alt')]
    write.table(snps,paste0('imp_1kg_hap/chr',chr,'_snpfile'),row.names=F,col.names=F,quote=F)
    #return(snps)
}

read.chr.rate<-function(chr){
    rfile<-paste0('imp_1kg_hap/chr',chr,'_ratefile')
    gen.map<-read.genmap(chr)
    site<-read.site(chr)
    interp.x<-gen.map[,1]
    interp.y<-gen.map[,3]
    interp.xout<-site$pos
    interp.res<-approx(x=interp.x,y=interp.y,xout=interp.xout,rule=2)
    interp.str<-lapply(interp.res,function(x)paste(x,collapse=' '))
    cat(paste0(':sites:',length(interp.res[[1]]),'\n'),file=rfile)
    cat(paste0(interp.str[[1]],'\n'),file=rfile,append=T)
    cat(paste0(interp.str[[2]],'\n'),file=rfile,append=T)
}

run.gm<-function(){
    mclapply(1:22,read.chr,mc.cores=10)
}
run.rf<-function(){
    mclapply(1:22,read.chr.rate,mc.cores=10)
}

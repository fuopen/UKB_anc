library(parallel)

base.dir<-'/data/muscovy/not-backed-up/shu/hapmix/'

get.exam<-function(){
    rk<-system('cat /data/muscovy/not-backed-up/shu/hapmix/HapmixRelease/example.par',intern=T)

    rk.strs<-strsplit(rk,':')
    rk.dt<-do.call(rbind,rk.strs)
    return(rk.dt)
}

if(!exists('par.dt')){
    par.dt<-get.exam()
}

get.chr.par<-function(chr,alt.dir){
    par.dt[9,2]<-paste0(base.dir,'imp_1kg_hap/ceu_hap_chr',chr)
    par.dt[10,2]<-paste0(base.dir,'imp_1kg_hap/yri_hap_chr',chr)
    par.dt[11,2]<-paste0(base.dir,'imp_1kg_hap/chr',chr,'_snpfile')
    par.dt[12,2]<-paste0(base.dir,'imp_1kg_hap/chr',chr,'_snpfile')
    geno.snp.file<-paste0(alt.dir,'/data/chr',chr,'_snpfile')
    if(!file.exists(geno.snp.file)){
        file.copy(paste0(base.dir,'imp_1kg_hap/chr',chr,'_snpfile'),geno.snp.file)
    }
    par.dt[13,2]<-paste0(alt.dir,'/data/chr',chr,'_snpfile')

    par.dt[14,2]<-paste0(alt.dir,'/data/chr',chr,'_geno')
    par.dt[15,2]<-paste0(alt.dir,'/data/sample_ind')

    par.dt[18,2]<-paste0(base.dir,'imp_1kg_hap/chr',chr,'_ratefile')
    par.dt[19,2]<-'UKB_AFR'
    if(!dir.exists(paste0(alt.dir,'/out'))){
        dir.create(paste0(alt.dir,'/out'))
    }
    par.dt[20,2]<-chr
    par.dt[21,2]<-paste0(alt.dir,'/out/')
    par.dt[22,2]<-'DIPLOID'
    par.file<-paste0(alt.dir,'/chr',chr,'.par')
    par.str<-paste0(par.dt[,1],':',par.dt[,2])
    write.table(par.str,par.file,row.names=F,col.names=F,quote=F)
}

run.par<-function(){
    sub.dirs<-list.dirs('/data/muscovy/not-backed-up/shu/hapmix/imp_geno',recursive=F)
    fz<-function(xdir){
       lapply(1:22,get.chr.par,alt.dir=xdir)
    }
    mclapply(sub.dirs,fz,mc.cores=15)
}

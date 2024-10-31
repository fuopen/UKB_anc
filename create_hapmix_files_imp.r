library(parallel)

geno.dir<-'tmp_ukb_hap'

hap.dir<-'tmp_1kg_hap'

out.dir<-'hapmix_all_files'

if(!exists('af.gender')){
    af.gender<-readRDS('african_gender.rds')
}

read.chr<-function(chr){
    ceu.file<-paste0(hap.dir,'/ceu_chr',chr,'.rds')
    yri.file<-paste0(hap.dir,'/yri_chr',chr,'.rds')

    afr.file<-paste0(geno.dir,'/afr_chr',chr,'_update.rds')

    ceu<-readRDS(ceu.file)
    yri<-readRDS(yri.file)
    afr<-readRDS(afr.file)

    colnames(ceu)<-ceu[1,]
    colnames(yri)<-yri[1,]
    ceu<-ceu[-1,]
    yri<-yri[-1,]
    ref.anno<-ceu[,1:5]
    gen.anno<-afr[,1:5]
    ref.str<-paste(chr,ref.anno[,2],ref.anno[,3],ref.anno[,4],ref.anno[,5],sep=':')
    afr.str<-paste(chr,gen.anno[,2],gen.anno[,3],gen.anno[,4],gen.anno[,5],sep=':')
    if(all(ref.str %in% afr.str) && all(afr.str %in% ref.str)){
        afr<-afr[match(ref.str,afr.str),6:ncol(afr)]
        ceu<-ceu[,6:ncol(ceu)]
        yri<-yri[,6:ncol(yri)]
        write.table(ceu,paste0(out.dir,'/ceu_hapfile_chr',chr),row.names=F,col.names=F,quote=F,sep='')
        write.table(yri,paste0(out.dir,'/yri_hapfile_chr',chr),row.names=F,col.names=F,quote=F,sep='')
        write.table(afr,paste0(out.dir,'/afr_geno_chr',chr),row.names=F,col.names=F,quote=F,sep='')
        saveRDS(ref.anno,paste0(out.dir,'/anno_chr',chr,'.rds'))
        #if(!file.exists(paste0(out.dir,'/african_ind'))){
        #    af.gender<-af.gender[match(colnames(afr)[-(1:5)],as.character(af.gender$id)),]
        #    af.gender$pop<-rep('african_ukb',nrow(af.gender))
        #    write.table(af.gender,paste0(out.dir,'/african_ind'),row.names=F,col.names=F,quote=F)
        #}
    }
    else{
        stop('ref and afr doesn\'t match!')
    }
}

run.all<-function(){
    mclapply(1:22,read.chr,mc.cores=10)
}

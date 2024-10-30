library(parallel)

get.1kg.chr<-function(chr){
    ceu.strs<-system(paste0('bcftools view /data/muscovy/not-backed-up/shu/1kg/raw_data/CEU_1kg_chr',chr,".bcf|grep -v '##'|cut -f 1-5,10-|sed 's/|//g'"),intern=T)
    yri.strs<-system(paste0('bcftools view /data/muscovy/not-backed-up/shu/1kg/raw_data/YRI_1kg_chr',chr,".bcf|grep -v '##'|cut -f 1-5,10-|sed 's/|//g'"),intern=T)
    line.num<-read.table(paste0('tmp_1kg_hap/chr',chr,'_line_num.txt'),as.is=T,header=F)[[1]]
    ceu.use.str<-ceu.strs[line.num]
    yri.use.str<-yri.strs[line.num]
    ceu.mat<-do.call(rbind, strsplit(ceu.use.str,'\t'))
    yri.mat<-do.call(rbind, strsplit(yri.use.str,'\t'))
    saveRDS(ceu.mat,paste0('tmp_1kg_hap/ceu_chr',chr,'.rds'))
    saveRDS(yri.mat,paste0('tmp_1kg_hap/yri_chr',chr,'.rds'))
}


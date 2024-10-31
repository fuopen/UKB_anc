library(parallel)
source('prepare_par_file.r')

res.dir<-'/data/muscovy/not-backed-up/shu/hapmix/hapmix_res'

file.dir<-'hapmix_all_files'

if(!dir.exists(res.dir)){
    dir.create(res.dir)
    #system(paste0('ln -s /data/muscovy/not-backed-up/shu/hapmix/HapmixRelease/bin '
}

if(!exists('file.seq')){
    file.seq<-seq(1,8600,100)
}

split.file<-function(batch.num){
    b.start<-file.seq[batch.num]
    b.end<-b.start+99
    
    b.dir<-paste0(res.dir,'/batch_',b.start,'_',b.end)
    if(!dir.exists(b.dir)){
        dir.create(b.dir)
    }
        
    b.chrs.dir<-paste0(b.dir,'/data/')
    if(!dir.exists(b.chrs.dir)){
        dir.create(b.chrs.dir)
    }
    sample.cmd<-paste0('sed -n ',b.start,',',b.end,'p ',file.dir,'/african_ind > ',b.chrs.dir,'/batch_sample_ind')
    system(sample.cmd)

    cp.fun<-function(i){
        file.copy(paste0(file.dir,'/chr',i,'_snpfile'),b.chrs.dir)
        get.chr.par(i,b.dir)
        split.cmd<-paste0("awk -v start=",b.start," -v end=",b.end," 'BEGIN{FS=\"\"}{for(k=start;k<=end;k++){if(k!=end){printf(\"%s\",$k)}else{printf(\"%s\\n\",$k)}}}' ",file.dir,"/afr_geno_chr",i," > ",b.chrs.dir,"/batch_geno_chr",i)
        record<-system(split.cmd,intern=T)
        return(record)
    }
    mclapply(1:22,cp.fun,mc.cores=16)
}

run.split<-function(){
    lapply(1:length(file.seq),split.file)
}

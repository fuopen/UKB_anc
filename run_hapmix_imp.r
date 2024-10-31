library(parallel)

res.dir<-'imp_geno'

if(!exists('file.seq')){
    file.seq<-seq(1,8500,50)
}

run.batch<-function(batch.num){
    b.start<-file.seq[batch.num]
    b.end<-b.start+49
    b.dir<-paste0(res.dir,'/batch_',b.start,'_',b.end,'/')
    lapply(1:22,function(i){
        run.cmd<-paste0('./run_hapmix_imp.sh ',b.dir,' ',i)
        rc<-system(run.cmd,intern=T)
    })
}

run<-function(start,end,nc=16){
    mclapply(start:end,run.batch,mc.cores=nc)
}

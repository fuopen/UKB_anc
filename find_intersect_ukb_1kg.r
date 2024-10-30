library(parallel)

find.inse.chr<-function(chr){
    dir.tmp<-'tmp_dir/'
    if(!dir.exists(dir.tmp)){
        dir.create(dir.tmp)
    }
    no_flip.cmd<-paste0('grep -xf ukb_snp_list/chr',chr,'_ukb.txt 1kg_snp_list/chr',chr,'_1kg.txt >',dir.tmp,'/chr',chr,'_match.txt')
    flip.cmd<-paste0('./flip.sh ',chr)
    system(no_flip.cmd)
    system(flip.cmd)
}

run.chrs<-function(){
    mclapply(1:22,find.inse.chr,mc.cores=12)
}


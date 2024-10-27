library(parallel)
beta.dir<-'new_gwas_results_bg3_merge_whradjbmi/'

se.dir<-beta.dir

tmp.dir<-'tmp_merge_whradjbmi/'

if(!dir.exists(tmp.dir)){
    dir.create(tmp.dir)
}

anno.file<-'new_gwas_snps_anno_withaf.rds'
anno<-''

if(!exists('anno')){
    anno<-readRDS(anno.file)
    anno.str<-strsplit(anno$vid,':')
    anno$chr<-sapply(anno.str,function(x)x[[1]])
    anno$pos<-sapply(anno.str,function(x)x[[2]])
    anno$a1<-sapply(anno.str,function(x)x[[3]])
    anno$a2<-sapply(anno.str,function(x)x[[4]])
}

if(!exists('sp.size')){
	sp.size<-342784
	names(sp.size)<-'whradjbmi'
}

read.files<-function(f_code){
    beta.file<-paste0(beta.dir,f_code,'.ac.rds')
    beta<-readRDS(beta.file)
    rf.ac<-cbind.data.frame(anno[,-1],beta[,c(1:2,4)])
	rf.list<-list(ac=rf.ac)
    return(rf.list)
}

run.ldsc<-function(f_code){
    rf<-read.files(f_code)
    rf.ac<-rf$ac[,c(1,4:7,2:3,8:10)]
    rf.ac[,10]<-10^(-rf.ac[,10])
    spsize<-sp.size[[f_code]]

    colnames(rf.ac)[8:10]<-c('beta','se','pval')
	if(!file.exists(paste0(tmp.dir,f_code,'.actmp.txt'))){
    	write.table(rf.ac,paste0(tmp.dir,f_code,'.actmp.txt'),row.names=F,col.names=T,quote=F)
    }
	munge.ac.cmd<-paste0('/data/muscovy/shu/ldsc/munge_sumstats.py --out ',tmp.dir,f_code,'.ac --merge-alleles /data/muscovy/shu/ldsc/eur_w_ld_chr/w_hm3.snplist --N ',spsize,' --sumstats ',tmp.dir,f_code,'.actmp.txt')
    ldsc.ac.cmd<-paste0('/data/muscovy/shu/ldsc/ldsc.py --h2 ',tmp.dir,f_code,'.ac.sumstats.gz --ref-ld-chr /data/muscovy/shu/ldsc/eur_w_ld_chr/ --out ',tmp.dir,f_code,'.AC --w-ld-chr /data/muscovy/shu/ldsc/eur_w_ld_chr/')
    extract.ac.cmd<-paste0('grep -E "h2:|^Lambda GC:|Mean Chi|^Intercept" ',tmp.dir,f_code,'.AC.log')
    system(munge.ac.cmd)
    system(ldsc.ac.cmd)
    f.ac<-system(extract.ac.cmd,intern=T)
    convert.table<-function(fstr,type){
        ks<-sapply(strsplit(fstr,':'),function(x)x[[2]])
        ks.data<-sapply(strsplit(ks,' '),function(x)x[[2]])
        ks.data<-as.numeric(ks.data)
        names(ks.data)<-paste0(type,'.',c('h2','lambda_GC','MeanChi2','Intercept'))
        return(ks.data)
    }
    ac.res<-convert.table(f.ac,'AC')
    rm.cmd.ac<-paste0('rm ',tmp.dir,f_code,'.actmp.txt ',tmp.dir,f_code,'.ac.sumstats.gz ')
	ret<-ac.res
}
run.ldsc2<-function(f_code){
    spsize<-sp.size[[f_code]]
	munge.ac.cmd<-paste0('/data/muscovy/shu/ldsc/munge_sumstats.py --out ',tmp.dir,f_code,'.ac --merge-alleles /data/muscovy/shu/ldsc/eur_w_ld_chr/w_hm3.snplist --N ',spsize,' --sumstats ',tmp.dir,f_code,'.actmp.txt')
    ldsc.ac.cmd<-paste0('/data/muscovy/shu/ldsc/ldsc.py --h2 ',tmp.dir,f_code,'.ac.sumstats.gz --ref-ld-chr /data/muscovy/shu/ldsc/eur_w_ld_chr/ --out ',tmp.dir,f_code,'.AC --w-ld-chr /data/muscovy/shu/ldsc/eur_w_ld_chr/')
    
    extract.ac.cmd<-paste0('grep -E "h2:|^Lambda GC:|Mean Chi|^Intercept" ',tmp.dir,f_code,'.AC.log')
    system(munge.ac.cmd)
    system(ldsc.ac.cmd)
    f.ac<-system(extract.ac.cmd,intern=T)
    convert.table<-function(fstr,type){
        ks<-sapply(strsplit(fstr,':'),function(x)x[[2]])
        ks.data<-sapply(strsplit(ks,' '),function(x)x[[2]])
        ks.data<-as.numeric(ks.data)
        names(ks.data)<-paste0(type,'.',c('h2','lambda_GC','MeanChi2','Intercept'))
        return(ks.data)
    }
    ac.res<-convert.table(f.ac,'AC')
    rm.cmd.ac<-paste0('rm ',tmp.dir,f_code,'.actmp.txt ',tmp.dir,f_code,'.ac.sumstats.gz ')
	ret<-ac.res
}

run.all.ldsc<-function(){
	codes<-names(sp.size)
    me<-mclapply(codes,function(code)try(run.ldsc(code)),mc.cores=12)
    names(me)<-codes
    return(me)
}

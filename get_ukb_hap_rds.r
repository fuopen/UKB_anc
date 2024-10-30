library(parallel)

get.ukb.chr<-function(chr){
    afr.strs<-system(paste0("zgrep -v '##' /data/muscovy/not-backed-up/shu/ukb_27960/geno/african_chr",chr,".vcf.gz|cut -f 1-5,10-"),intern=T)
    if(!chr %in% 21:22){
        line.num<-read.table(paste0('tmp_ukb_hap/chr',chr,'_line_num.txt'),as.is=T,header=F)[[1]]
        afr.use.str<-afr.strs[line.num]
        afr.mat<-do.call(rbind,strsplit(afr.use.str,'\t'))
        colnames(afr.mat)<-afr.mat[1,]
        afr.mat<-afr.mat[-1,]
        afr.mat[afr.mat=='1,0,0,1'|afr.mat=='0,1,1,0']='1'
        afr.mat[afr.mat=='1,0,1,0']='0'
        afr.mat[afr.mat=='0,1,0,1']='2'
        match.str<-strsplit(afr.mat[,3],',')
        afr.mat[,1]<-sapply(match.str,function(x)x[[2]])
        afr.mat[,3]<-sapply(match.str,function(x)x[[1]])
        saveRDS(afr.mat,paste0('tmp_ukb_hap/afr_chr',chr,'_update.rds'))
    }
    else{
        line.num.match<-read.table(paste0('tmp_ukb_hap/chr',chr,'_line_num_match.txt'),as.is=T,header=F)[[1]]
        line.num.flip.match<-read.table(paste0('tmp_ukb_hap/chr',chr,'_line_num_flip_match.txt'),as.is=T,header=F)[[1]]
        afr.use.match.str<-afr.strs[line.num.match]
        afr.use.flip.str<-afr.strs[line.num.flip.match]
        afr.mat.match<-do.call(rbind,strsplit(afr.use.match.str,'\t'))
        afr.mat.flip.match<-do.call(rbind,strsplit(afr.use.flip.str,'\t'))
        colnames(afr.mat.match)<-afr.mat.match[1,]
        colnames(afr.mat.flip.match)<-colnames(afr.mat.match)
        afr.mat.match<-afr.mat.match[-1,]
        afr.mat.flip.match[,c(4:5)]<-afr.mat.flip.match[,c(5:4)]

        afr.mat.match[afr.mat.match=='1,0,0,1'|afr.mat.match=='0,1,1,0']='1'
        afr.mat.match[afr.mat.match=='1,0,1,0']='0'
        afr.mat.match[afr.mat.match=='0,1,0,1']='2'
        
        afr.mat.flip.match[afr.mat.flip.match=='1,0,0,1'|afr.mat.flip.match=='0,1,1,0']='1'
        afr.mat.flip.match[afr.mat.flip.match=='1,0,1,0']='2'
        afr.mat.flip.match[afr.mat.flip.match=='0,1,0,1']='0'

        afr.geno.mat<-rbind(afr.mat.match,afr.mat.flip.match)
        saveRDS(afr.geno.mat,paste0('tmp_ukb_hap/afr_chr',chr,'_update.rds'))
    }
}

run<-function(seq){
    mclapply(seq,get.ukb.chr,mc.cores=11)
}

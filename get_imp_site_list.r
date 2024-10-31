library(parallel)
ukb.sites.dir<-'/data/garganey/not-backed-up/shu/tmp_impute/hapmix_8.5K'
GP3.sites.dir<-'1kg_snp_list'

out.dir<-'imp_site_list'
out.geno.dir<-'imp_geno'

if(!dir.exists(out.geno.dir)){
    dir.create(out.geno.dir)
}

if(!dir.exists(out.dir)){
    dir.create(out.dir)
}

read.chr<-function(chr){
    ukb.anno<-readRDS(paste0(ukb.sites.dir,'/chr',chr,'_anno.rds'))
    GP3.anno<-read.table(paste0(GP3.sites.dir,'/chr',chr,'_1kg.txt'),as.is=T,header=T,comment.char='',check.names=F)
    ukb.anno<-ukb.anno[,c(1,4,2,6,5)]
    ukb.str<-paste0(ukb.anno[,1],ukb.anno[,2],ukb.anno[,3],ukb.anno[,4],ukb.anno[,5])
    GP3.str<-paste0(GP3.anno[,1],GP3.anno[,2],GP3.anno[,3],GP3.anno[,4],GP3.anno[,5])
    inse.str<-intersect(ukb.str,GP3.str)
    ukb.match.ids<-match(inse.str,ukb.str)
    GP3.match.ids<-match(inse.str,GP3.str)
    saveRDS(ukb.match.ids,paste0(out.dir,'/ukb_match_ids_chr',chr,'.rds'))
    write.table(GP3.anno[GP3.match.ids,],paste0(out.dir,'/GP3_match_chr',chr,'_snps.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}

run.read.chr<-function(){
    mclapply(1:22,read.chr,mc.cores=12)
}

if(!exists('sp.ids')){
    sp.ids<-seq(1,8500,500)
    sp.seq.mat<-cbind(sp.ids,sp.ids+499)
}

get.ukb.geno.chr<-function(chr){
    ids<-readRDS(paste0(out.dir,'/ukb_match_ids_chr',chr,'.rds'))
    get.each<-function(x){
        geno.file<-paste0(ukb.sites.dir,'/chr',chr,'_range_',x[1],'_',x[2],'.rds')
        geno<-readRDS(geno.file)
        geno<-geno[ids,]
        print(c(x[1],x[2]))
        sub.seq<-seq(x[1],x[2],50)
        sapply(sub.seq,function(i){
            print(c(i,i+49))
            j<-i%%500
            geno.sub<-geno[,j:(j+49)]
            sub.dir<-paste0(out.geno.dir,'/batch_',i,'_',i+49,'/')
            if(!dir.exists(sub.dir)){
                dir.create(sub.dir)
                dir.create(paste0(sub.dir,'data/'))
            }
            write.table(geno.sub,paste0(sub.dir,'data/chr',chr,'_geno'),row.names=F,col.names=F,quote=F,sep='')
        })
    }
    apply(sp.seq.mat,1,get.each)
}
run.ukb.geno.chr<-function(){
    mclapply(1:22,get.ukb.geno.chr,mc.cores=12)
}

get.1kg.geno.chr<-function(chr){
    cmd<-paste0('./get_1kg_line_number_imp.sh ',chr)
    system(cmd)
    ceu.strs<-system(paste0('bcftools view /data/muscovy/not-backed-up/shu/1kg/raw_data/CEU_1kg_chr',chr,".bcf|grep -v '##'|cut -f 1-5,10-|sed 's/|//g'"),intern=T) 
    yri.strs<-system(paste0('bcftools view /data/muscovy/not-backed-up/shu/1kg/raw_data/YRI_1kg_chr',chr,".bcf|grep -v '##'|cut -f 1-5,10-|sed 's/|//g'"),intern=T) 
    line.num<-read.table(paste0('imp_1kg_hap/chr',chr,'_line_num.txt'),as.is=T,header=F)[[1]]
    ceu.use.str<-ceu.strs[line.num]
    yri.use.str<-yri.strs[line.num]
    ceu.mat<-do.call(rbind, strsplit(ceu.use.str,'\t'))
    yri.mat<-do.call(rbind, strsplit(yri.use.str,'\t'))
    saveRDS(ceu.mat,paste0('imp_1kg_hap/ceu_chr',chr,'.rds'))
    saveRDS(yri.mat,paste0('imp_1kg_hap/yri_chr',chr,'.rds'))
}

run.1kg.geno<-function(){
    mclapply(1:22,get.1kg.geno.chr,mc.cores=12)
}

get.1kg.hap.chr<-function(chr){
    ceu<-readRDS(paste0('imp_1kg_hap/ceu_chr',chr,'.rds'))
    yri<-readRDS(paste0('imp_1kg_hap/yri_chr',chr,'.rds'))
    ceu.anno<-ceu[-1,1:5]
    ceu.hap<-ceu[-1,-(1:5)]
    yri.hap<-yri[-1,-(1:5)]
    saveRDS(ceu.anno,paste0('imp_1kg_hap/snp_chr',chr,'.rds'))
    write.table(ceu.hap,paste0('imp_1kg_hap/ceu_hap_chr',chr),row.names=F,col.names=F,quote=F,sep='')
    write.table(yri.hap,paste0('imp_1kg_hap/yri_hap_chr',chr),row.names=F,col.names=F,quote=F,sep='')
}
run.1kg.hap<-function(){
    mclapply(1:22,get.1kg.hap.chr,mc.cores=12)
}

###run recomb_rate_imp.r here
if(!exists('af.sample')){
    af.sample<-readRDS('african_ids_clean.rds')
}
get.ukb.sp.chr<-function(){
    get.each<-function(x){
        geno.file<-paste0(ukb.sites.dir,'/chr',22,'_range_',x[1],'_',x[2],'.rds')
        geno<-readRDS(geno.file)
        sps<-colnames(geno)
        sps<-as.integer(sapply(strsplit(sps,'_'),function(s)s[[1]]))
        sub.seq<-seq(x[1],x[2],50)
        sapply(sub.seq,function(i){
            j<-i%%500
            sps.sub<-sps[j:(j+49)]
            sub.dir<-paste0(out.geno.dir,'/batch_',i,'_',i+49,'/')
			if(!dir.exists(sub.dir)){
				dir.create(sub.dir)
			}
			if(!dir.exists(paste0(sub.dir,'/data/'))){
				dir.create(paste0(sub.dir,'/data'))
			}
            af.sub<-af.sample[match(sps.sub,af.sample[,1]),]
            write.table(af.sub,paste0(sub.dir,'data/sample_ind'),row.names=F,col.names=F,quote=F,sep=' ')
        })
    }
    apply(sp.seq.mat,1,get.each)
}


library(parallel)
if(!exists('phased.eu.af.geno')){
    phased.eu.af.geno<-readRDS('trio_geno/trio_eu_af_pgs_new.rds')
    phased.gg<-cbind(rowSums(phased.eu.af.geno[,1:2]),rowSums(phased.eu.af.geno[,3:4]),rowSums(phased.eu.af.geno[,5:6]),rowSums(phased.eu.af.geno[,7:8]),rowSums(phased.eu.af.geno[,9:10]),rowSums(phased.eu.af.geno[,11:12]))
}

if(!exists('imp.anno')){
    imp.anno<-readRDS('imp_anno_withmaf.rds')
    imp.order<-1:nrow(imp.anno)
    imp.chr<-mclapply(1:22,function(i)imp.order[!is.na(imp.anno$maf) & imp.anno[,2]==i & imp.anno$maf>=0.01 ],mc.cores=6)
}

if(!exists('phased.anc')){
    phased.anc<-readRDS('trio_geno/trio_eu_af_AF_anc_new.rds')
}

if(!exists('unphased.bh')){
    unphased.bh<-readRDS('trio_geno/trio_eu_af_AF_anc_bothhet_new.rds')
}

if(!exists('ra')){
    ra<-readRDS('uncertainty_updated.rds')
    ra.ids<-ra<0.1
}

if(!exists('hap.anc')){
    hap.anc<-readRDS('trio_geno/trios_hap_anc_unphased_phased.rds')
}

if(!exists('t7')){
    t7<-readRDS('traits7_results/Standing_height_EU_pgs.50.rds')

    colnames(phased.eu.af.geno)[seq(1,12,2)]->ca
    match(ca,colnames(t7))->mat.id
}
if(!exists('file.st')){
    anc.dir<-'ancestry_geno/'
    file.st<-mat.id%/%50*50+1
    files<-sapply(file.st,function(i)paste0(anc.dir,'sample_',i,'_',i+49,'.rds'))
    anc.dir2<-'ancestry_geno_check2/'
    anc.files<-sapply(file.st,function(i)paste0(anc.dir2,'sample_',i,'_',i+49,'.rds'))
}

if(!exists('unphased.geno')){
    unphased.geno<-readRDS('trio_geno/trio_6_ind_original_genotype.rds')
}

read.geno<-function(){
    g.col<-mat.id%%50
    g.cols<-(g.col-1)*2
    sf<-mclapply(files,readRDS,mc.cores=6)
    sf.col<-lapply(1:length(sf),function(i)sf[[i]][,g.cols[i]+1:2])
    return(sf.col)
}
read.geno.anc<-function(){
    g.col<-mat.id%%50
    sf<-mclapply(anc.files,readRDS,mc.cores=6)
    sf.col<-lapply(1:length(sf),function(i)sf[[i]][,g.col[i]])
    sf.dt<-do.call(cbind,sf.col)
    return(sf.dt)
}

if(!exists('unphased.eu.af.geno')){
    unphased.eu.af.geno<-readRDS('trio_geno/trio_eu_af_pgs_unphased.rds')
}
if(!exists('unphased.anc')){
    unphased.anc<-readRDS('trio_geno/trio_eu_af_pgs_unphased_af_anc.rds')
}

get.eu.anc<-function(){
    ug.eu.anc<-1-unphased.anc
    hap.eu.anc.up.1<-sapply(hap.anc[[1]],function(x)x[,1])
    hap.eu.anc.up.2<-sapply(hap.anc[[1]],function(x)x[,3])
    hap.eu.anc.ph.1<-sapply(hap.anc[[2]],function(x)x[,1])
    hap.eu.anc.ph.2<-sapply(hap.anc[[2]],function(x)x[,3])
    lz<-list(ug.eu.anc,hap.eu.anc.up.1,hap.eu.anc.up.2,hap.eu.anc.ph.1,hap.eu.anc.ph.2)
}
make.geno<-function(trio.geno){
    ap<-apply(trio.geno,1,function(x){
        y<-x
        if(all(x=='1')|all(x=='0')|all(x=='2')){
            return(y)
        }
        if(x[1]=='0'){
            if(x[2]=='1' && x[3]=='1'){
                y[2]<-'3'
                y[3]<-'3'
            }
            else if(x[2]=='0' && x[3]=='1'){
                y[3]<-'3'
            }
            else if(x[2]=='1' && x[3]=='0'){
                y[2]<-'3'
            }
            else{#bad cases, child 0 and one/both of parents are 2,e
                y<-y
                #y<-rep(NA,3)
            }    
            return(y)
        }
        else if(x[1]=='2'){
            if(x[2]=='1' && x[3]=='1'){
                y[2]<-'4'
                y[3]<-'4'
            }
            else if(x[2]=='2' && x[3]=='1'){
                y[3]<-'4'
            }
            else if(x[2]=='1' && x[3]=='2'){
                y[2]<-'4'
            }
            else{
                #y<-rep(NA,3)
                y<-y
            }
            return(y)
        }
        else{#child is het 1
            if(x[2]=='0' && x[3]=='1'){
                y[1]='4'
                y[3]='4'
            }
            else if(x[2]=='0' && x[3]=='2'){
                y[1]='4'
            }
            else if(x[2]=='1' && x[3]=='0'){
                y[1]='3'
                y[2]='4'
            }
            else if(x[2]=='1' && x[3]=='2'){
                y[1]='4'
                y[2]='3'
            }
            else if(x[2]=='2' && x[3]=='0'){
                y[1]='3'
            }
            else if(x[2]=='2' && x[3]=='1'){
                y[1]='3'
                y[3]='3'
            }
            else{
                #y<-rep(NA,3)
                y<-y
            }
            return(y)
        }
    })
    return(t(ap))
}

if(!exists('phased.geno')){
    phased.geno<-readRDS('trio_geno/genotype_after_conversion.rds')
}

if(!exists('unphased.anc.copy')){
    unphased.anc.copy<-readRDS('trio_geno/trio_copy_num_euro_african.rds')
}

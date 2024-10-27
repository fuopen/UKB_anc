#***********************************************************************
#***********************************************************************
library(parallel)
pthr<--log10(5e-8)

res.dir<-'/data/muscovy/not-backed-up/shu/'
var.anno.file<-paste0(res.dir,'gwas/new_gwas_snps_anno.rds')

if(!exists('c.anno')){
	c.anno<-readRDS('/data/muscovy/not-backed-up/shu/gwas/cp_dict/cp_anno.rds')
}

read.anno<-function(){
    anno<-readRDS(var.anno.file)
    colnames(anno)<-c('variant','rsid')
    anno.info<-strsplit(anno[,1],':')
    anno$chr<-sapply(anno.info,function(x)x[[1]])
    anno$rsid<-anno[,2]
    anno$pos<-sapply(anno.info,function(x)x[[2]])
    
    anno$ichr<-as.integer(anno$chr)
    anno$pos<-as.integer(anno$pos)

    return(anno)
}

if(!exists('ac.anno')){
    ac.anno<-read.anno()
}

if(!exists('match.ord')){
	match.ord.file<-'/datas/garganey/not-backed-up/shu/gwas_new/merge_match.rds'
	match.ord<-readRDS(match.ord.file)
}

read.trait<-function(file.ac,file.pc){
    pv<-ac.anno[,c('ichr','variant','rsid','pos')]
    if(file.exists(file.ac)&& file.exists(file.pc)){
        rds.ac<-readRDS(file.ac)
        pv$ac.beta<-rds.ac$ac.beta
        pv$ac.pv<-rds.ac$ac.pv
        rds.pc<-readRDS(file.pc)[match.ord,]
		pv$pc.beta<-rds.pc$pc.beta
        pv$pc.pv<-rds.pc$pc.pv
        return(pv)
    }
    else{
        return(NULL)
    }
}

chrome.filter0<-function(dt,which.pv='pc.pv',wd=1e5){
    dt<-dt[!is.na(dt[[which.pv]]),]
    dt<-dt[order(dt$pos),]
    dt<-dt[dt[[which.pv]]>=pthr,]
	if(nrow(dt)==0){
		return(NULL)
	}
	if(nrow(dt)==1){
		return(dt)
	}
	peak.dt<-dt[1,]
	for(i in 2:nrow(dt)){
		dtz<-dt[1:i,]
		dtz<-dtz[dt$pos[i]-dtz$pos<=wd,]
		#print(dtz)
		dtz.peak.idx<-which.max(dtz[[which.pv]])
		dtz.peak<-dtz[dtz.peak.idx,]
		#print(dtz.peak)
		if(nrow(dtz)>dtz.peak.idx){
			next
		}
		else{
			if(peak.dt[nrow(peak.dt),'variant']%in%dtz$variant){
				peak.dt[nrow(peak.dt),]<-dtz.peak
			}
			else{
				peak.dt<-rbind(peak.dt,dtz.peak)
			}
		}
		#print(peak.dt)
	}
	peak.dt
}


get.peak<-function(pv.dt,which.pv='pc.pv'){#from read.trait
    #my.peak<-by(pv.dt,pv.dt$ichr,chrome.filter)
    #my.peak<-do.call(rbind.data.frame,my.peak)
    my.peak<-mclapply(1:22,function(i)chrome.filter0(pv.dt[pv.dt$ichr==i,],which.pv=which.pv),mc.cores=12)
    return(my.peak)
}

output.peak<-function(file.ac,file.pc){
    pdt<-read.trait(file.ac,file.pc)
    gp.pc<-get.peak(pdt,'pc.pv')
    gp.ac<-get.peak(pdt,'ac.pv')
	return(list(dt=pdt,pc.pv=gp.pc,ac.pv=gp.ac))
}

make.acpc<-function(op.list){
	dall<-op.list[[1]]
	dpc<-do.call(rbind,op.list[[2]])
	dac<-do.call(rbind,op.list[[3]])
	dallz<-dall[dall$variant %in%c(dpc$variant,dac$variant),]
}

if(!exists('ukb.sp.size')){
	ukb.sp.size<-readRDS('/data/muscovy/not-backed-up/shu/gwas/ukb_sample_size.rds')
}

#generate_new_result<-function(fid){
#	filename.ac<-paste0('/data/muscovy/not-backed-up/shu/gwas/new_gwas_results_38traits_bg3_merge/',fid,'.rds')
#	filename.pc<-paste0('/datas/garganey/not-backed-up/shu/gwas_new/pc40_38traits_bg3_merge/',fid,'.rds')
#	wc<-readRDS(paste0('/data/muscovy/not-backed-up/shu/gwas/cp_dict/f.50.0.0.rds'))
#	op<-output.peak(filename.ac,filename.pc)
#	mp<-make.acpc(op)
#	c.anno$alternate_id<-gsub('_',':',c.anno$alternate_id)
#	m1<-match(mp$variant,c.anno$alternate_id)
#	mp.cn<-c.anno[m1,c('maf','info_score'),]	
#	mp<-cbind(mp,mp.cn)
#	mp$wb_sample_copy<-wc[m1]
#	mp	
#}

test.plot<-function(op.list,chr,which.pv='pc.pv',outfile=paste0('test_capture/chr',chr,'_',which.pv,'.png')){
    whole.chr<-op.list[[1]]
    peak.chr<-op.list[[which.pv]][[chr]]
    png(outfile,width=4.75,height=3.45,units='in', res=1200,pointsize=4)
    plot(whole.chr[whole.chr$ichr==chr,'pos'],whole.chr[whole.chr$ichr==chr,which.pv],col='lightgreen',pch='.',xlab=paste0('chr',chr),ylab=expression("-log"[10]*"(p-value)"),main=NA)
    if(!is.null(peak.chr)){
		points(peak.chr$pos,peak.chr[[which.pv]],col='red',pch=1)
    }
	dev.off()
}

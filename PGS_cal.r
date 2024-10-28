library(parallel)

snp.dir<-'new_clump_snps_bg3_r_0.1/'

if(!exists('pvs')){
	pvs<-c('5e-08','1e-05','0.01','0.05')
}

if(!exists('pheno.list')){
	pheno.list<-read.table('gwas_files_path.txt',as.is=T,sep='\t')[[1]]
}
read.info<-function(pv){
	each.chr<-function(chr){
		aa<-readRDS(paste0(snp.dir,'chr',chr,'_pv_',pv,'.rds'))
		phenos<-lapply(pheno.list,function(pheno){
			aa.snps<-rownames(aa)[aa[,paste0(pheno,'.AC')]]
		})
		names(phenos)<-pheno.list
		return(phenos)
	}
	lz<-mclapply(1:22,each.chr,mc.cores=12)
	rliz<-lapply(pheno.list,function(pp){
		lzk<-lapply(1:22,function(i)lz[[i]][[pp]])
	})
	names(rliz)<-pheno.list
	rliz	
}

#do sanity check

sanity.check<-function(r1,r2){
	p1<-names(r1)
	p2<-names(r2)
	if(!identical(p1,p2)){
		stop('Error:not the same names!')
	}
	p12<-sapply(p1,function(x)sapply(1:22,function(i)mean(r1[[x]][[i]]%in%r2[[x]][[i]])))
	p21<-sapply(p2,function(x)sapply(1:22,function(i)mean(r2[[x]][[i]]%in%r1[[x]][[i]])))
	py<-list(p12=p12,p21=p21)
}

#calculate the PGS

geno.dir<-'geno_35k'

clump.dir<-'PGS_snps_betas'

if(!exists('pv.ll')){
	pv.ll<-lapply(pvs,function(pv){
		rr<-readRDS(paste0(clump.dir,'/pv',pv,'.rds'))
	})
	names(pv.ll)<-pvs
}

out.dir<-'PGS_35K/'

cal.pgs.snp<-function(chr){
	chr.anno<-readRDS(paste0(geno.dir,'/chr',chr,'_anno.rds'))
	chr.geno<-readRDS(paste0(geno.dir,'/chr',chr,'_sample35K.rds'))
	chr.anno$vid<-paste0(chr.anno$CHR,':',chr.anno$POS,':',chr.anno$ALT,':',chr.anno$COUNTED)
	each.pv<-function(pv){
		ll<-pv.ll[[pv]]
		each.pheno<-function(lld){
			if(nrow(lld[[chr]])==0){
				ret<-matrix(0,nrow=ncol(chr.geno))
				rownames(ret)<-colnames(chr.geno)
				return(ret)
			}
			else{
				chr.ids<-match(rownames(lld[[chr]]),chr.anno$vid)
				chr.geno<-as.matrix(chr.geno)[chr.ids,]
				chr.pgs<-t(lld[[chr]])%*%chr.geno
				return(t(chr.pgs))
			}
		}
		ll.z<-mclapply(ll,each.pheno,mc.cores=12)
		ll.zt<-do.call(cbind,ll.z)
		colnames(ll.zt)<-names(ll)
		ll.zt
	}
	epp<-lapply(pvs,each.pv)
	names(epp)<-pvs
	saveRDS(epp,paste0(out.dir,'chr',chr,'.pgs.rds'))
	return(epp)
}


run.all.pgs<-function(){
	pgs.all<-lapply(1:22,cal.pgs.snp)
	
	pgs.pvs<-lapply(pvs,function(x)lapply(pgs.all,function(pgs.chr){
		pgs.pv.chr<-pgs.chr[[x]]
	}))
	names(pgs.pvs)<-pvs
	pgs.new<-lapply(pgs.pvs,function(pv)Reduce('+',pv))
	return(pgs.new)
	#pgs.new<-Reduce('+',pgs.all)
}

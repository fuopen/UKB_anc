library(parallel)

if(!exists('gwas.wb')){
	gwas.wb<-readRDS('gwas_wb_ids.rds')
}

if(!exists('chr15.snps')){
	chr15.snps<-read.table('plink_odds_chr15.raw',as.is=T,header=T,sep='\t',check.names=F)
}

if(!exists('pc40')){
	pc<-readRDS('v2_40PCs.rds')
	colnames(pc)<-paste0('PC',1:40)
}

pc.pred.geno<-function(chr){
	cc<-get(paste0('chr',chr,'.snps'))
	cc<-cc[gwas.wb,-(1:6)]
	pc<-pc[gwas.wb,1:20]
	pc.rg<-lapply(cc,function(x){
		lk<-lm(x ~ .,data=pc)
		lk.sum<-summary(lk)
	})
	pc.rg.coeff<-lapply(pc.rg,function(x)x$coefficients)
	pc.rg.rsq<-sapply(pc.rg,function(x)x$adj.r.squared)
	ll<-list(coeff=pc.rg.coeff,rsq=pc.rg.rsq)
	lk<-cor(cc,pc,use='complete.obs')
	nz<-list(reg=ll,cor=lk)
}

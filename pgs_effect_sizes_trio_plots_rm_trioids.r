library(parallel)

if(!exists('pgs.anno')){
	pgs.anno<-readRDS('data/imp_anno_afmaf_1kg_EAF.rds')
	pgs.anno$id<-paste0(pgs.anno[,2],':',pgs.anno[,4],':',pgs.anno[,5],':',pgs.anno[,6])
	pgs.chr<-mclapply(1:22,function(i){
		pchr<-pgs.anno[pgs.anno[,2]==i,]
		pchrt<-!is.na(pchr$eu.maf)  & pchr$eu.maf>=0.01
	},mc.cores=6)
	pgs.chrids<-mclapply(1:22,function(i){
		pchr<-pgs.anno[pgs.anno[,2]==i,]
		pu<-pchr[!is.na(pchr$eu.maf)  & pchr$eu.maf>=0.01,]
		pu.ids<-paste0(pu[,2],':',pu[,4],':',pu[,5],':',pu[,6])
	},mc.cores=6)
}

if(!exists('covar')){
	covar<-readRDS('eu_af_cov.rds')
}

if(!exists('eu.af.pheno')){
	eu.af.pheno<-readRDS('eu_af_59_traits_pheno.rds')
}

if(!exists('af.anc')){
	af.anc<-readRDS('African_ancestry_8003.rds')
}

test.all.dir<-'trio_geno_newtraits_fxboot/'

if(!dir.exists(test.all.dir)){
	dir.create(test.all.dir)
}

if(!exists('pheno.ids')){
	#pheno.ids<-names(snp.anno.29)
	pheno.ids<-colnames(eu.af.pheno)
}
# as requested by UKB we removed the trio-ids here 
trio.ids<-c()

if(!exists('eu.af.allele.frq')){
	eu.af.allele.frq<-readRDS('data/pgs_8003_all_chr_af_eu_EAF_estimate.rds')
}

out.geno.dir<-'trio_geno_P4'
filt.dir<-'imp_mask_region_5mb'

if(!exists('snp.betas')){
	snp.betas<-readRDS('59_traits_snps_beta_pv0.05.rds')
}

cal.pgs<-function(trio.ids){
    trio.ids<-as.character(trio.ids)
	each.chr<-function(chr){
		chr.all<-readRDS(paste0(out.geno.dir,'/triophased_geno_chr',chr,'.rds'))
		chr.ee<-chr.all[pgs.chr[[chr]],seq(1,24,by=4)]
		chr.ea<-chr.all[pgs.chr[[chr]],seq(2,24,by=4)]
		chr.aa<-chr.all[pgs.chr[[chr]],seq(3,24,by=4)]
		chr.ae<-chr.all[pgs.chr[[chr]],seq(4,24,by=4)]
		anc.all<-readRDS(paste0(out.geno.dir,'/triophased_anc_chr',chr,'.rds'))
		anc.ee<-anc.all[,seq(1,18,by=3)]
		anc.ea<-anc.all[,seq(2,18,by=3)]
		anc.aa<-anc.all[,seq(3,18,by=3)]
		
		filter.chr<-readRDS(paste0(filt.dir,'/chr',chr,'_5mb.rds'))[,trio.ids]
		eu.af.frq.chr<-eu.af.allele.frq[[chr]][,c('af.eaf','eu.eaf')]
		each.pheno<-function(pheno){
			#ep.snp<-snp.anno.29[[pheno]]
			#ep.beta<-snp.beta.29[[pheno]]
			ep.beta<-snp.betas[[pheno]][[chr]]
			if(nrow(ep.beta)==0){
				zz<-rep(0,ncol(chr.ee))
				chr.pgs<-data.frame(MC.EU=zz,MC.AF=zz,NMC.EU=zz,NMC.AF=zz)
				return(chr.pgs)
			}
			ep.snp<-rownames(ep.beta)
			ep.snp.chr.ids<-ep.snp%in%rownames(eu.af.frq.chr)
			ep.beta.chr<-ep.beta[ep.snp.chr.ids,1]
			ep.snp.chr<-ep.snp[ep.snp.chr.ids]
			ep.chr.ids<-match(ep.snp.chr,rownames(eu.af.frq.chr))
			ee.geno<-chr.ee[ep.chr.ids,]
			ea.geno<-chr.ea[ep.chr.ids,]
			aa.geno<-chr.aa[ep.chr.ids,]
			ae.geno<-chr.ae[ep.chr.ids,]
			anc.ee<-anc.ee[ep.chr.ids,]
			anc.ea<-anc.ea[ep.chr.ids,]
			anc.aa<-anc.aa[ep.chr.ids,]
			f.chr<-filter.chr[ep.chr.ids,]
			pheno.frq.chr<-eu.af.frq.chr[ep.chr.ids,]
			ee.geno.chr<-ee.geno-2*pheno.frq.chr$eu.eaf*anc.ee
			ea.geno.chr<-ea.geno-pheno.frq.chr$eu.eaf*anc.ea
			aa.geno.chr<-aa.geno-2*pheno.frq.chr$af.eaf*anc.aa
			ae.geno.chr<-ae.geno-pheno.frq.chr$af.eaf*anc.ea
			ee.chr.pgs<-ep.beta.chr*ee.geno.chr
			ea.chr.pgs<-ep.beta.chr*ea.geno.chr
			aa.chr.pgs<-ep.beta.chr*aa.geno.chr
			ae.chr.pgs<-ep.beta.chr*ae.geno.chr
			ee.chr.pgs[f.chr==0]=0
			ea.chr.pgs[f.chr==0]=0
			aa.chr.pgs[f.chr==0]=0
			ae.chr.pgs[f.chr==0]=0
			chr.pgs<-data.frame(MC.EU=colSums(ee.chr.pgs)+colSums(ea.chr.pgs),MC.AF=colSums(aa.chr.pgs)+colSums(ae.chr.pgs))
		}
		ep.all<-mclapply(pheno.ids,each.pheno,mc.cores=12)
		names(ep.all)<-pheno.ids
		ep.all
	}
	mk<-lapply(1:22,each.chr)
}

cal.pgs53<-function(){
	if(!exists('cp')){
		cp<-cal.pgs()
	}
	pgs.allchr.each<-function(pheno){
		pae.eu<-do.call(cbind,lapply(cp,function(x)x[[pheno]]$MC.EU))
		pae.af<-do.call(cbind,lapply(cp,function(x)x[[pheno]]$MC.AF))
		pae.eu<-rowSums(pae.eu)
		pae.af<-rowSums(pae.af)
		#pae.dt<-data.frame(EE=pae.ee,EA=pae.ea,AA=pae.aa,AE=pae.ae)
		pgs.dt<-data.frame(EU=pae.eu,AF=pae.af)
	}
	ptest<-mclapply(pheno.ids,pgs.allchr.each,mc.cores=12)
	names(ptest)<-pheno.ids
	ptest
}

m.test.dt<-function(code){
    gtp<-get.trait.pgs(code)
    mc<-match(trio.ids,eu.af.ids)
    t.pheno<-eu.af.pheno[mc,paste0('f.',code,'.0.0')]
    t.dt<-data.frame(pheno=t.pheno,age=eu.af.agesex[mc,1],sex=eu.af.agesex[mc,2],eu.up.pgs=gtp[[1]][,1],af.up.pgs=gtp[[1]][,2],eu.p.pgs=gtp[[2]][,1],af.p.pgs=gtp[[2]][,2])
    t.dt<-cbind.data.frame(t.dt,eu.af.anc[mc,main.pops])
    t.dt$hap.af<-af.anc[match(trio.ids,af.anc.ids)]
    return(t.dt)
}


if(!exists('pgs.59')){
	pgs.59<-readRDS('59_traits_eu_af_pgs_mean0filter_pv0.05.rds')
}

if(!exists('pgs.up.59')){
	pgs.up.59<-readRDS('59_traits_trio_pv_0.05_replace_with_old29.rds')
}

get.trait.all.eu.af.pgs<-function(code){
	eu.af.t<-pgs.59[,paste0(code,c('.eu.pgs','.af.pgs')),drop=F]		
	my.dt<-cbind(eu.af.t,covar)
}

get.trait.trio.pgs<-function(code){
	my.dt<-pgs.up.59[[code]]
	my.dt<-lapply(my.dt,function(x){y<-x;colnames(y)<-c('eu.t','af.t');y})	
}

adjust.pgs<-function(code){
    all.pgs.dt<-get.trait.all.eu.af.pgs(code)
	colnames(all.pgs.dt)[1:2]<-c('eu.t','af.t')
    trio.pgs<-get.trait.trio.pgs(code)
	#rownames(trio.pgs[[2]])<-rownames(trio.pgs[[1]])
    rownames(trio.pgs$phased)<-rownames(trio.pgs$unphased)
	ggt.eu.coef<-coefficients(lm(eu.t~.,data=all.pgs.dt[,-2]))
    ggt.af.coef<-coefficients(lm(af.t~.,data=all.pgs.dt[,-1]))
    #trio.dt<-all.pgs.dt[rownames(trio.pgs[[1]]),]
    trio.dt<-all.pgs.dt[as.character(trio.ids),]
    reg.dt<-as.matrix(cbind(rep(1,nrow(all.pgs.dt)),all.pgs.dt[,-(1:2)]))
    coreg.eu.score<-reg.dt%*%matrix(ggt.eu.coef,ncol=1)
    coreg.af.score<-reg.dt%*%matrix(ggt.af.coef,ncol=1)
    
    res.eu.score<-all.pgs.dt[,1]-coreg.eu.score[,1]
    res.af.score<-all.pgs.dt[,2]-coreg.af.score[,1]

    coreg.eu.trio.score<-coreg.eu.score[rownames(trio.pgs[[1]]),1]
    coreg.af.trio.score<-coreg.af.score[rownames(trio.pgs[[1]]),1]
    
    trio.adjust.unphased.pgs<-trio.pgs[[1]]-cbind(coreg.eu.trio.score,coreg.af.trio.score)
    trio.adjust.phased.pgs<-trio.pgs[[2]]-cbind(coreg.eu.trio.score,coreg.af.trio.score)
    #got partial scores for unphased/phased score
    eu.anc<-1-af.anc

    trio.anc.match<-match(rownames(trio.pgs[[1]]),names(af.anc))
    af.anc.trio<-af.anc[trio.anc.match]
    eu.anc.trio<-eu.anc[trio.anc.match]
    
    mean.eu.score<-sum(res.eu.score*eu.anc)/sum(eu.anc)
    mean.af.score<-sum(res.af.score*af.anc)/sum(af.anc)
    res.eu.adj.score<-res.eu.score-mean.eu.score*eu.anc
    res.af.adj.score<-res.af.score-mean.af.score*af.anc
    var.eu.score<-sum(res.eu.adj.score^2)/sum(eu.anc)
    var.af.score<-sum(res.af.adj.score^2)/sum(af.anc)
    
    sd.eu.score<-sqrt(var.eu.score)
    sd.af.score<-sqrt(var.af.score)
    
    trio.norm.unphased.pgs<-cbind((trio.adjust.unphased.pgs[,1]-mean.eu.score*eu.anc.trio)/sqrt(var.eu.score),(trio.adjust.unphased.pgs[,2]-mean.af.score*af.anc.trio)/sqrt(var.af.score))
    trio.norm.phased.pgs<-cbind((trio.adjust.phased.pgs[,1]-mean.eu.score*eu.anc.trio)/sqrt(var.eu.score),(trio.adjust.phased.pgs[,2]-mean.af.score*af.anc.trio)/sqrt(var.af.score))
    trio.norm.pgs<-list(trio.norm.unphased.pgs,trio.norm.phased.pgs)
    names(trio.norm.pgs)<-c('unphased','phased')
    return(trio.norm.pgs)
}

norm.pgs.all<-function(){
	ap<-mclapply(pheno.ids,adjust.pgs,mc.cores=8)
	names(ap)<-pheno.ids
	return(ap)
}
plot.npg<-function(){
    if(!exists('npa')){
        npa<<-norm.pgs.all()
    }
	ph<-readRDS('29_traits_for_plot_annotation.rds')
	n.title<-ph[match(names(npa),ph[,1]),2]
    cp<-c('trio1.child','trio1.father','trio1.mother','trio2.child','trio2.father','trio2.mother')
    #my.cols<-c('red','green','blue','purple','orange','darkblue')
	my.cols<-rainbow(length(npa))
    n.pheno<-names(npa)
    pdf(paste0(test.all.dir,'/Trio_newtraits_normalised_unphased_phased_EU.pdf'))
    plot(0,0,type='n',xlab='unphased normalised PGS',ylab='phased normalised PGS',xlim=c(-1.25,1.25),ylim=c(-1.25,1.25),main='European')
    abline(a=0,b=1,lwd=0.1)
    for(i in 1:length(npa)){
        points(npa[[i]][[1]][,1],npa[[i]][[2]][,1],col=my.cols[i],pch=as.character(1:6),cex=0.5)
    }
    legend('topleft',col=my.cols,legend=n.title,pch=19,cex=0.4,ncol=3)
    legend('bottomright',legend=paste0(1:6,':',cp),cex=0.5)
   	dev.off()
    pdf(paste0(test.all.dir,'/Trio_newtraits_normalised_unphased_phased_AF.pdf'))
    plot(0,0,type='n',xlab='unphased normalised PGS',ylab='phased normalised PGS',xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),main='African')
    abline(a=0,b=1,lwd=0.1)
    for(i in 1:length(npa)){
        points(npa[[i]][[1]][,2],npa[[i]][[2]][,2],col=my.cols[i],pch=as.character(1:6),cex=0.5)
    }
    legend('topleft',col=my.cols,legend=n.title,pch=19,cex=0.4,ncol=3)
    #legend('topleft',col=my.cols,legend=n.title,pch=19,cex=0.4)
    legend('bottomright',legend=paste0(1:6,':',cp),cex=0.5)
    dev.off()
	nqz<-do.call(rbind,lapply(npa,function(x){
		qqq<-cbind(x$unphased,x$phased)
		colnames(qqq)<-c('up.eu','up.af','ph.eu','ph.af')
		qqq
	}))
}

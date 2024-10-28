library(parallel)

if(!exists('af.anc')){
    af.anc<-readRDS('African_ancestry_8003.rds')
}

if(!exists('covar')){
	covar<-readRDS('eu_af_cov.rds')
}

if(!exists('eu.af.pheno')){
    eu.af.pheno<-readRDS('eu_af_53_traits_pheno.rds')
}

if(!exists('gds.mean0.0.05')){
	gds.mean0.0.05<-readRDS('eu_af_53_traits_pgs_pv0.05.rds')
}
if(!exists('gds.mean0.1e4')){
	gds.mean0.1e4<-readRDS('eu_af_53_traits_pgs_pv1e4.rds')
}

if(!exists('trait.info.53')){
	trait.info.53<-readRDS('/homes/shu/workdir/ng_revision/PGS_calculation/53_traits_ldsc_for_figure.rds')
}

if(!exists('fxboot')){
	fxboot<-readRDS('bootstrap_8003.mat.rds')
}

if(!exists('acboot')){
	acboot<-readRDS('bootstrap_cut_0.35_0.55_0.78_0.95.rds')
}

if(!exists('obs.eu.cov')){
	obs.eu.cov<-readRDS('ObsEu_covar.rds')
}

if(!exists('obs.eu.pheno')){
	obs.eu.pheno<-readRDS('ObsEu_53_traits_pheno.rds')
}

if(!exists('obs.eu.pgs.0.05')){
	obs.eu.pgs.0.05<-readRDS('ObsEu_53_traits_pgs_pv0.05.rds')
}
if(!exists('obs.eu.pgs.1e4')){
	obs.eu.pgs.1e4<-readRDS('ObsEu_53_traits_pgs_pv1e4.rds')
}

if(!exists('obs.eu.boot')){
	obs.eu.boot<-readRDS('ObsEu_bootstrap.rds')
}

get.trait.matrix.raw.euaf.binanc<-function(code,gds.mean0,breaks=c(0,0.35,0.55,0.78,0.95,1)){
	gds.euaf<-gds.mean0[[code]]
	if(!identical(rownames(gds.euaf),names(af.anc))){
		stop('ids dont match')
	}
	af.anc.cut<-cut(af.anc,breaks=breaks)
	gtp.list<-by(gds.euaf,af.anc.cut,function(x)x)
	pheno<-eu.af.pheno[,code,drop=F]
	tz<-lapply(gtp.list,function(x){
		pheno.x<-pheno[rownames(x),1]
		afanc.x<-af.anc[rownames(x)]
		covar.x<-covar[rownames(x),]
		t.dt<-data.frame(pheno=pheno.x,eu.pgs=x[,1],af.pgs=x[,2],af.anc=afanc.x)
		t.dt<-cbind.data.frame(t.dt,covar.x)	
	})
	names(tz)<-names(gtp.list)
	tz
}


cut.samples.byanc.raw.euaf.binanc<-function(code,gds.mean0,breaks=c(0,0.35,0.55,0.78,0.95,1),bootstrap=NA){
	dt.code<-get.trait.matrix.raw.euaf.binanc(code,gds.mean0,breaks)
	if(!is.na(bootstrap)){
		for(x in names(dt.code)){
			dt.code[[x]]<-dt.code[[x]][acboot[[x]][,bootstrap],]
		}
	}
	cut.reg<-function(cut.dt){
		t.sum<-summary(lm(pheno ~.,data=cut.dt))
		t.sq<-as.vector(t.sum$coefficients[c('eu.pgs','af.pgs'),1:2])
		names(t.sq)<-c('eu.pgs.beta','af.pgs.beta','eu.pgs.sd','af.pgs.sd')	
		t.sq
	}
	dcut<-lapply(dt.code,cut.reg)
	dcut.dt<-do.call(rbind,dcut)
	colnames(dcut.dt)<-paste0(code,'.',colnames(dcut.dt))
	dcut.dt
}

get.trait.matrix.central.AF<-function(code,gds.mean0){
	gtp<-gds.mean0[[code]]
	pheno<-eu.af.pheno[,code]
	t.dt<-data.frame(pheno=pheno,eu.pgs=gtp[,1],af.pgs=gtp[,2],pgs=gtp[,1]+gtp[,2],af.anc=af.anc)
	t.dt<-cbind.data.frame(t.dt,covar)
	gtp.nm<-rownames(gtp)
	eu.cent.lm<-lm(eu.pgs ~ .,data=t.dt[,-c(1,3,4,5)],na.action=na.exclude)
	eu.cent<-resid(eu.cent.lm)
	af.cent.lm<-lm(af.pgs ~ .,data=t.dt[,-c(1,2,4,5)],na.action=na.exclude)
	af.cent<-resid(af.cent.lm)
	pgs.cent.lm<-lm(pgs ~.,data=t.dt[,-c(1:3,5)],na.action=na.exclude)
	all.cent<-resid(pgs.cent.lm)
	
	t.dt$eu.pgs<-eu.cent
	t.dt$af.pgs<-af.cent
	t.dt$pgs<-all.cent
	t.dt
}

get.trait.matrix.central.EU<-function(code,obs.eu.pgs){
	gb.pgs<-obs.eu.pgs[,code]
	pheno<-obs.eu.pheno[,code]
	t.dt<-data.frame(pheno=pheno,pgs=gb.pgs)
	t.dt<-cbind.data.frame(t.dt,obs.eu.cov)
	gtp.nm<-rownames(gb.pgs)
	
	pgs.cent.lm<-lm(pgs ~.,data=t.dt[,-1],na.action=na.exclude)
	all.cent<-resid(pgs.cent.lm)
	t.dt$pgs<-all.cent
	t.dt
}

cut.samples.byanc.AF<-function(code,gds.mean0,bootstrap=NA){
	dt.code<-get.trait.matrix.central.AF(code,gds.mean0)
	cc<-dt.code[dt.code$af.anc>0.95,]
	if(!is.na(bootstrap)){
		cc<-cc[acboot[[5]][,bootstrap],]
	}
	cut.reg.euaf<-function(cut.dt){
		t.sum.euaf<-summary(lm(pheno ~.,data=cut.dt[,-4]))
		t.sum.all<-summary(lm(pheno ~.,data=cut.dt[,-(2:3)]))
		t.sq.euaf<-as.vector(t.sum.euaf$coefficients[c('eu.pgs','af.pgs'),1:2])
		t.sq.all<-t.sum.all$coefficients['pgs',1:2]
		t.sq<-c(t.sq.euaf,t.sq.all)
		names(t.sq)<-c('eu.pgs.beta','af.pgs.beta','eu.pgs.sd','af.pgs.sd','pgs.beta','pgs.sd')
		t.sq
	}
	dcut.dt<-cut.reg.euaf(cc)
	dcut.dt
}

cut.samples.byanc.EU<-function(code,obs.eu.pgs,bootstrap=NA){
	dt.code<-get.trait.matrix.central.EU(code,obs.eu.pgs)
	if(!is.na(bootstrap)){
		dt.code<-dt.code[obs.eu.boot[,bootstrap],]
	}
	cc<-dt.code
	cut.reg.euaf<-function(cut.dt){
		t.sum.all<-summary(lm(pheno ~.,data=cut.dt))
		t.sq<-t.sum.all$coefficients['pgs',1:2]
		names(t.sq)<-c('pgs.beta','pgs.sd')
		t.sq
	}
	dcut.dt<-cut.reg.euaf(cc)
	dcut.dt
}

if(!exists('pheno.ids')){
	pheno.ids<-trait.info.53$trait
}

get.coefficients.raw.euaf<-function(gds.mean0,breaks=c(0,0.35,0.55,0.78,0.95,1),iter.num=1000){
	csb.all<-do.call(cbind,mclapply(pheno.ids,cut.samples.byanc.raw.euaf.binanc,gds.mean0=gds.mean0,breaks=breaks,bootstrap=NA,mc.cores=12))
	csb.real<-csb.all
	csb.btsr.all<-mclapply(1:iter.num,function(i){
		c.all<-do.call(cbind,lapply(pheno.ids,cut.samples.byanc.raw.euaf.binanc,gds.mean0=gds.mean0,breaks=breaks,bootstrap=i))
		lb<-c.all
	},mc.cores=16)
	rrv<-list(real=csb.real,bootstrap=csb.btsr.all)
}

get.coefficients.AF<-function(gds.mean0,iter.num=1000){
	csb.all.list<-mclapply(pheno.ids,cut.samples.byanc.AF,gds.mean0=gds.mean0,bootstrap=NA,mc.cores=12)
	names(csb.all.list)<-pheno.ids
	csb.all<-do.call(rbind,csb.all.list)
	csb.btsr.all<-mclapply(1:iter.num,function(i){
		c.alls<-lapply(pheno.ids,cut.samples.byanc.AF,gds.mean0=gds.mean0,bootstrap=i)
		names(c.alls)<-pheno.ids
		ca<-do.call(rbind,c.alls)
	},mc.cores=12)
	rrv<-list(real=csb.all,bootstrap=csb.btsr.all)
}

get.coefficients.EU<-function(obs.eu.pgs,iter.num=1000){
	csb.all.list<-mclapply(pheno.ids,cut.samples.byanc.EU,obs.eu.pgs=obs.eu.pgs,bootstrap=NA,mc.cores=12)
	names(csb.all.list)<-pheno.ids
	csb.all<-do.call(rbind,csb.all.list)
	csb.btsr.all<-mclapply(1:iter.num,function(i){
		c.alls<-lapply(pheno.ids,cut.samples.byanc.EU,obs.eu.pgs=obs.eu.pgs,bootstrap=i)
		names(c.alls)<-pheno.ids
		ca<-do.call(rbind,c.alls)
	},mc.cores=12)
	rrv<-list(real=csb.all,bootstrap=csb.btsr.all)
}

gclist.af.file.0.05<-'53_AF_binned_fxboot_pv0.05.rds'
gclist.eu.file.0.05<-'53_EU_binned_fxboot_pv0.05.rds'
gclist.af.file.1e4<-'53_AF_binned_fxboot_pv1e4.rds'
gclist.eu.file.1e4<-'53_EU_binned_fxboot_pv1e4.rds'

gc.filter<-function(gc.list){
	gc.list$real<-gc.list$real[pheno.ids,]
	gc.list$bootstrap<-lapply(gc.list$bootstrap,function(x)x[pheno.ids,])
	gc.list
}

if(!exists('gclist.af.0.05')){
	gclist.af.0.05.53<-readRDS(gclist.af.file.0.05)
	gclist.af.0.05<-gc.filter(gclist.af.0.05.53)
}

if(!exists('gclist.eu.0.05')){
	gclist.eu.0.05.53<-readRDS(gclist.eu.file.0.05)
	gclist.eu.0.05<-gc.filter(gclist.eu.0.05.53)
}

if(!exists('gclist.af.1e4')){
	gclist.af.1e4.53<-readRDS(gclist.af.file.1e4)
	gclist.af.1e4<-gc.filter(gclist.af.1e4.53)
}

if(!exists('gclist.eu.1e4')){
	gclist.eu.1e4.53<-readRDS(gclist.eu.file.1e4)
	gclist.eu.1e4<-gc.filter(gclist.eu.1e4.53)
}
	
generate.plot.af.dt<-function(gc.list){
	cols.beta<-'pgs.beta'
	cols.sd<-'pgs.sd'
	ccs.pgs.real.beta<-gc.list[[1]][,grep(cols.beta,colnames(gc.list[[1]]))]
	ccs.pgs.real.sd<-gc.list[[1]][,grep(cols.sd,colnames(gc.list[[1]]))]
	ccs.pgs.btr.betas<-lapply(gc.list[[2]],function(x)x[,grep(cols.beta,colnames(x))])
	ccs.pgs.btr.sd<-lapply(gc.list[[2]],function(x)x[,grep(cols.sd,colnames(x))])
	pgs.btsr.beta<-simplify2array(ccs.pgs.btr.betas)
	pgs.btsr.sd<-simplify2array(ccs.pgs.btr.sd)
	v.beta<-list(ccs.pgs.real.beta,pgs.btsr.beta)
	v.sd<-list(ccs.pgs.real.sd,pgs.btsr.sd)
	vr<-list(beta=v.beta,sd=v.sd)
}
generate.plot.eu.dt<-function(gc.list){
	cols.beta<-'pgs.beta'
	cols.sd<-'pgs.sd'
	ccs.pgs.real.beta<-gc.list[[1]][,grep(cols.beta,colnames(gc.list[[1]])),drop=F]
	ccs.pgs.real.sd<-gc.list[[1]][,grep(cols.sd,colnames(gc.list[[1]])),drop=F]
	ccs.pgs.btr.betas<-lapply(gc.list[[2]],function(x)x[,grep(cols.beta,colnames(x)),drop=F])
	ccs.pgs.btr.sd<-lapply(gc.list[[2]],function(x)x[,grep(cols.sd,colnames(x)),drop=F])
	pgs.btsr.beta<-simplify2array(ccs.pgs.btr.betas)
	pgs.btsr.sd<-simplify2array(ccs.pgs.btr.sd)
	v.beta<-list(ccs.pgs.real.beta,pgs.btsr.beta)
	v.sd<-list(ccs.pgs.real.sd,pgs.btsr.sd)
	vr<-list(beta=v.beta,sd=v.sd)
}

generate.plot.CI.af<-function(gc.list){
	gpd<-generate.plot.af.dt(gc.list)[[1]]
	gpd.org<-gpd[[1]]
	gpd.ci.025<-apply(gpd[[2]],c(1,2),quantile,probs=0.025)
	gpd.ci.975<-apply(gpd[[2]],c(1,2),quantile,probs=0.975)
	gpCI.each.all<-lapply(rownames(gpd.org),function(i){
		rbind(gpd.org[i,],gpd.ci.025[i,],gpd.ci.975[i,])
	})
	names(gpCI.each.all)<-rownames(gpd.org)
	gpCI<-gpCI.each.all
}
generate.plot.CI.eu<-function(gc.list){
	gpd<-generate.plot.eu.dt(gc.list)[[1]]
	gpd.org<-gpd[[1]]
	gpd.ci.025<-apply(gpd[[2]],c(1,2),quantile,probs=0.025)
	gpd.ci.975<-apply(gpd[[2]],c(1,2),quantile,probs=0.975)
	gpCI.each.all<-lapply(rownames(gpd.org),function(i){
		rbind(gpd.org[i,],gpd.ci.025[i,],gpd.ci.975[i,])
	})
	names(gpCI.each.all)<-rownames(gpd.org)
	gpCI<-gpCI.each.all
}

get.inverse.weight.p4<-function(mat){
	inv.mat<-(1/mat)^2
	inv.rowsum<-rowSums(inv.mat)
	inv.mat<-inv.mat/inv.rowsum
}

get.inverse.weight<-function(mat){
	inv.mat<-(1/mat)^2
	inv.rowsum<-colSums(inv.mat)
	inv.mat<-t(inv.mat)/inv.rowsum
	inv.mat<-t(inv.mat)
}

estimated.weight.mean.euaf<-function(gc.list){
	gpd<-generate.plot.af.dt(gc.list)
	
	gpd.real.blk.beta.all<-gpd[[1]][[1]]
	gpd.real.blk.sd.all<-gpd[[2]][[1]]

	gpd.real.blk.beta.all.eu<-gpd.real.blk.beta.all[,seq(1,ncol(gpd.real.blk.beta.all),by=2)]
	gpd.real.blk.beta.all.af<-gpd.real.blk.beta.all[,seq(2,ncol(gpd.real.blk.beta.all),by=2)]
	
	gpd.real.blk.sd.all.eu<-gpd.real.blk.sd.all[,seq(1,ncol(gpd.real.blk.sd.all),by=2)]
	gpd.real.blk.sd.all.af<-gpd.real.blk.sd.all[,seq(2,ncol(gpd.real.blk.sd.all),by=2)]
	
	gwd.all.eu<-gpd.real.blk.beta.all.eu*get.inverse.weight.p4(gpd.real.blk.sd.all.eu)
	gwd.mean.all.eu<-rowSums(gwd.all.eu)
	
	gwd.all.af<-gpd.real.blk.beta.all.af*get.inverse.weight.p4(gpd.real.blk.sd.all.af)
	gwd.mean.all.af<-rowSums(gwd.all.af)

	gwd.btr.all.eu<-do.call(rbind,mclapply(1:1000,function(i){
		m1<-gpd[[1]][[2]][,,i][,seq(1,ncol(gpd[[1]][[2]][,,i]),by=2)]
		m2<-get.inverse.weight.p4(gpd[[2]][[2]][,,i][,seq(1,ncol(gpd[[2]][[2]][,,i]),by=2)])
		mz<-rowSums(m1*m2)
	},mc.cores=10))
	gwd.btr.all.af<-do.call(rbind,mclapply(1:1000,function(i){
		m1<-gpd[[1]][[2]][,,i][,seq(2,ncol(gpd[[1]][[2]][,,i]),by=2)]
		m2<-get.inverse.weight.p4(gpd[[2]][[2]][,,i][,seq(2,ncol(gpd[[2]][[2]][,,i]),by=2)])
		mz<-rowSums(m1*m2)
	},mc.cores=10))

	gwd.all.cis.eu<-apply(gwd.btr.all.eu,2,quantile,probs=c(0.025,0.975))
	gwd.all.cis.af<-apply(gwd.btr.all.af,2,quantile,probs=c(0.025,0.975))
	
	emw.all.eu<-rbind(gwd.mean.all.eu,gwd.all.cis.eu)
	emw.all.af<-rbind(gwd.mean.all.af,gwd.all.cis.af)
	
	emw.list<-list(eu=emw.all.eu,af=emw.all.af)
}

get.all.dt<-function(pv,breaks=c(0,0.35,0.55,0.78,0.95,1)){
	pdir<-'plot_bin_anc_filter_53traits_pv_0.05_1e4/'
	if(!dir.exists(pdir)){
		dir.create(pdir)
	}
	cols.53<-c('.eu.pgs.beta','.af.pgs.beta','.eu.pgs.sd','.af.pgs.sd')
	pheno.cols<-as.vector(sapply(pheno.ids,function(x)paste0(x,cols.53)))

	filter.53<-function(euaf.list){
		euaf.list$real<-euaf.list$real[,pheno.cols]
		euaf.list$bootstrap<-lapply(euaf.list$bootstrap,function(x)x[,pheno.cols])	
		euaf.list
	}
	euaf.file<-paste0(pdir,'euaf_anc_by8003_pv',pv,'.rds')
	if(!file.exists(euaf.file)){
		euaf<-get.coefficients.raw.euaf(breaks)
		saveRDS(euaf,euaf.file)
	}
	else{
		euaf<-readRDS(euaf.file)
		euaf<-filter.53(euaf)
	}
	euaf.plot.file<-paste0(pdir,'euaf_anc_by8003_plot_pv',pv,'.rds')
	if(!file.exists(euaf.plot.file)){
		euaf.plot<-estimated.weight.mean.euaf(euaf)
		saveRDS(euaf.plot,euaf.plot.file)
	}
	else{
		euaf.plot<-readRDS(euaf.plot.file)
	}
	lv<-list(euaf=euaf.plot)
}

estimated.weight.mean.af<-function(gc.list){
	gpd<-generate.plot.af.dt(gc.list)

	gpd.real.blk.beta.all<-gpd[[1]][[1]]
	gpd.real.blk.sd.all<-gpd[[2]][[1]]
	
	gwd.all<-gpd.real.blk.beta.all*get.inverse.weight(gpd.real.blk.sd.all)
	gwd.mean.all<-colSums(gwd.all)

	gwd.btr.all<-do.call(rbind,mclapply(1:1000,function(i){
		m1<-gpd[[1]][[2]][,,i]
		m2<-get.inverse.weight(gpd[[2]][[2]][,,i])
		mz<-colSums(m1*m2)
	},mc.cores=10))

	gwd.all.cis<-apply(gwd.btr.all,2,quantile,probs=c(0.025,0.975))
	emw.all<-rbind(gwd.mean.all,gwd.all.cis)
}

estimated.weight.mean.eu<-function(gc.list){
	gpd<-generate.plot.eu.dt(gc.list)
	gpd.real.blk.beta.all<-gpd[[1]][[1]]
	gpd.real.blk.sd.all<-gpd[[2]][[1]]

	gwd.all<-gpd.real.blk.beta.all*get.inverse.weight(gpd.real.blk.sd.all)
	gwd.mean.all<-colSums(gwd.all)

	gwd.btr.all<-do.call(rbind,mclapply(1:1000,function(i){
		m1<-matrix(gpd[[1]][[2]][,,i],ncol=1)
		m2<-get.inverse.weight(matrix(gpd[[2]][[2]][,,i],ncol=1))
		mz<-colSums(m1*m2)
	},mc.cores=10))

	gwd.all.cis<-apply(gwd.btr.all,2,quantile,probs=c(0.025,0.975))
	emw.all<-rbind(gwd.mean.all,gwd.all.cis)
}

bin.ratio<-function(beta.list,Eu.obs.list){
	b.real<-beta.list$real[,grep('pgs\\.beta',colnames(beta.list$real))]
	b.real<-b.real[,grep('pgs\\.beta',colnames(b.real))]
	b.real.eu<-b.real[,grep('\\.eu',colnames(b.real))]
	b.real.af<-b.real[,grep('\\.af',colnames(b.real))]
	Eu.obs.real<-Eu.obs.list$real
	Eu.obs.boot<-Eu.obs.list$bootstrap
	colnames(b.real.eu)<-gsub('\\.eu.*','',colnames(b.real.eu))
	colnames(b.real.af)<-gsub('\\.af.*','',colnames(b.real.af))
	b.real.eu<-t(b.real.eu[,rownames(Eu.obs.real)])
	b.real.af<-t(b.real.af[,rownames(Eu.obs.real)])
	b.boot<-lapply(beta.list$bootstrap,function(x)x[,grep('pgs\\.beta',colnames(x))])
	b.boot.eu<-lapply(b.boot,function(x){
		y<-x[,grep('\\.eu',colnames(x))]
		colnames(y)<-gsub('\\.eu.*','',colnames(y))
		yz<-t(y[,rownames(Eu.obs.real)])
	})
	b.boot.af<-lapply(b.boot,function(x){
		y<-x[,grep('\\.af',colnames(x))]
		colnames(y)<-gsub('\\.af.*','',colnames(y))
		yz<-t(y[,rownames(Eu.obs.real)])
	})
	b.real.eu.ratio<-b.real.eu/Eu.obs.real[,1]
	b.real.af.ratio<-b.real.af/Eu.obs.real[,1]

	b.boot.eu.ratio<-lapply(1:length(b.boot.eu),function(i)b.boot.eu[[i]]/Eu.obs.boot[[i]][,1])
	b.boot.af.ratio<-lapply(1:length(b.boot.af),function(i)b.boot.af[[i]]/Eu.obs.boot[[i]][,1])
	b.boot.eu.ratios<-simplify2array(b.boot.eu.ratio)
	b.boot.af.ratios<-simplify2array(b.boot.af.ratio)
	b.boot.eu.ratio.weight<-apply(b.boot.eu.ratios,c(1,2),function(x)1/var(x))
	b.boot.af.ratio.weight<-apply(b.boot.af.ratios,c(1,2),function(x)1/var(x))
	b.boot.eu.ratio.norm<-apply(b.boot.eu.ratio.weight,2,function(x)x/sum(x))
	b.boot.af.ratio.norm<-apply(b.boot.af.ratio.weight,2,function(x)x/sum(x))
	b.boot.eu.ratio.new<-do.call(rbind,lapply(b.boot.eu.ratio,function(x)colSums(x*b.boot.eu.ratio.norm)))
	b.boot.af.ratio.new<-do.call(rbind,lapply(b.boot.af.ratio,function(x)colSums(x*b.boot.af.ratio.norm)))
	eu.wm.real<-colSums(b.real.eu.ratio*b.boot.eu.ratio.norm)
	af.wm.real<-colSums(b.real.af.ratio*b.boot.af.ratio.norm)
	eu.wm.new.ci<-apply(b.boot.eu.ratio.new,2,function(x)quantile(x,probs=c(0.025,0.975)))
	af.wm.new.ci<-apply(b.boot.af.ratio.new,2,function(x)quantile(x,probs=c(0.025,0.975)))
	eu.mt<-rbind(eu.wm.real,eu.wm.new.ci)
	af.mt<-rbind(af.wm.real,af.wm.new.ci)
	euaf<-list(EU=eu.mt,AF=af.mt)
}

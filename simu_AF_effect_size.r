library(parallel)

file.dir<-'/homes/shu/workdir/hapmix/simu_AF_by_EU_effect/sim-coeff'
out.dir<-'/homes/shu/workdir/hapmix/simu_AF_by_EU_effect/sim_af_new'

if(!dir.exists(out.dir)){
	dir.create(out.dir)
}

single.beta.simu<-function(beta_E,seed.id,rho,sigma=1){
	set.seed(seed.id)
	beta_A=rho*beta_E+sigma*sqrt(1-rho^2)*rnorm(length(beta_E),0,1)
}

if(!exists('file.list')){
	file.list<-list.files(file.dir)
}

if(!exists('rhos')){
	rhos<-c(seq(0.1,0.9,0.2),1)
}

if(!exists('sigmas')){
	sigmas<-rbind(c(0.9533,0.9997,1.0,0.9977,0.993,1.0),c(2.645,2.164,2.493,2.437,2.352,2.52))
	rownames(sigmas)<-c('beta.a0','beta.a5')
	colnames(sigmas)<-c('sim-coeff-100_cl','sim-coeff-100_rg','sim-coeff-10k_cl','sim-coeff-10k_rg','sim-coeff-1k_cl','sim-coeff-1k_rg')
}

each.file<-function(id){
	file<-file.list[id]
	zz<-read.table(paste0(file.dir,'/',file),as.is=T,header=T)
	sig.id<-gsub('-chr.*','',file)
	this.sigma<-sigmas[,sig.id]	
	each.rho<-function(rho,cid,Sigma)single.beta.simu(zz[[cid]],2023+id,rho=rho,sigma=Sigma)
	beta.1<-do.call(cbind,lapply(rhos,each.rho,cid='BETA',Sigma=this.sigma[1]))
	beta.2<-do.call(cbind,lapply(rhos,each.rho,cid='BETA.a5',Sigma=this.sigma[2]))
	
	colnames(beta.1)<-paste0('BETA','.AF.',rhos)
	colnames(beta.2)<-paste0('BETA.a5','.AF.',rhos)

	cb<-cbind(zz,beta.1,beta.2)

}

get.all.files<-function(){
	l.pp<-unique(gsub('.*coeff-(.*)-chr.*','\\1',file.list))
	each.uniq<-function(eun){
		ecc<-lapply(1:22,function(chr){
			#file.name<-paste0(file.dir,'/sim-coeff-',eun,'-chr',chr,'.tab')
			#efn<-read.table(file.name,as.is=T,header=T)
			file.name<-paste0('sim-coeff-',eun,'-chr',chr,'.tab')
			fid<-match(file.name,file.list)
			ef<-each.file(fid)
			write.table(ef,paste0(out.dir,'/AF-',eun,'-chr',chr,'.txt'),row.names=F,col.names=T,quote=F)
		})
	}
	eu.all<-mclapply(l.pp,each.uniq,mc.cores=length(l.pp))
}

geno.dir<-'/homes/shu/workdir/hapmix/simu_AF_by_EU_effect/tmp_pgen'
if(!exists('o.pp')){
	o.pp<-unique(gsub('AF-(.*)-chr.*','\\1',list.files(out.dir)))
}	
read.geno.chr<-function(chr){
	geno<-read.table(paste0(geno.dir,'/simugeno_chr',chr,'.raw'),as.is=T,header=T,check.names=F)
	geno.mat<-as.matrix(geno[,-(1:6)])
	geno.info<-colnames(geno.mat)
	geno.info.mat<-do.call(rbind,strsplit(geno.info,'_'))
	colnames(geno.info.mat)<-c('ID','A1')
	
	rownames(geno.mat)<-geno$IID
	o.pp.files<-paste0('AF-',o.pp,'-chr',chr,'.txt')
	each.file<-function(ff){
		eff<-read.table(paste0(out.dir,'/',ff),header=T,as.is=T)
		if(!all(eff$ID%in%geno.info.mat[,'ID'])){
			stop(paste0('error: not all snps in ',ff,' are in the geno file!'))
		}
		geno.info.ff<-geno.info.mat[match(eff$ID,geno.info.mat[,1]),,drop=F]
		eff$A1<-as.character(eff$A1)
		eff$A1[eff$A1=='TRUE']<-'T'	
		eff.mat<-as.matrix(eff[,-(1:2)])
		if(any(eff$A1!=geno.info.ff[,2])){
			flip.col<-which(eff$A1!=geno.info.ff[,2])
			eff$A1 <- geno.info.ff[,'A1']
		}
		else{
			flip.col<-NULL
		}
		rownames(eff.mat)<-paste0(eff$ID,'_',eff$A1)
		return(list(mat=eff.mat,flipcol=flip.col))
	}
	o.pp.dt<-lapply(o.pp.files,each.file)
	names(o.pp.dt)<-o.pp
	
	ind.pgs<-function(oppe){
		print(oppe)
		ef.dt<-o.pp.dt[[oppe]]$mat
		ef.geno<-geno.mat[,rownames(ef.dt),drop=F]
		if(!is.null(o.pp.dt[[oppe]]$flipcol)){
			ef.geno[,o.pp.dt[[oppe]]$flipcol]<-2-ef.geno[,o.pp.dt[[oppe]]$flipcol]
		}
		ef.pgs<-ef.geno%*%ef.dt
		return(ef.pgs)
	}
	all.chr.pgs<-lapply(o.pp,ind.pgs)
	names(all.chr.pgs)<-o.pp
	return(all.chr.pgs)
}

get.all.pgs<-function(){
	all.chrs.pgs<-mclapply(1:22,read.geno.chr,mc.cores=10)
	pp.names<-names(all.chrs.pgs[[1]])
	lq<-lapply(pp.names,function(x)lapply(all.chrs.pgs,function(y)y[[x]]))
	lq.pgs<-lapply(lq,function(x)Reduce('+',x))
	names(lq.pgs)<-pp.names
	lq.pgs
}

if(!exists('aa.pgs')){
	if(!file.exists('/homes/shu/workdir/hapmix/simu_AF_by_EU_effect/AA_PGS_new.rds')){
		aa.pgs<-get.all.pgs()
	}
	else{
		aa.pgs<-readRDS('/homes/shu/workdir/hapmix/simu_AF_by_EU_effect/AA_PGS_new.rds')
	}
}

generate.pheno<-function(i){
	apgs<-aa.pgs[[i]]
	norm.pgs<-apply(apgs,2,function(x){
		x.mean<-mean(x)
		x.sd<-sd(x)
		x.norm<-(x-x.mean)/x.sd
	})
	h2<-c(0.3,0.6)
	nc<-ncol(norm.pgs)
	set.seed(20230402+i)
	pheno.h2.0.3<-norm.pgs+rnorm(nrow(norm.pgs)*ncol(norm.pgs),0,sqrt((1-h2[1])/h2[1]))
	pheno.h2.0.6<-norm.pgs+rnorm(nrow(norm.pgs)*ncol(norm.pgs),0,sqrt((1-h2[2])/h2[2]))
	colnames(pheno.h2.0.3)<-paste0(colnames(pheno.h2.0.3),'_h0.3')
	colnames(pheno.h2.0.6)<-paste0(colnames(pheno.h2.0.6),'_h0.6')
	pp<-cbind(pheno.h2.0.3,pheno.h2.0.6)
	return(pp)
}

all.pheno<-function(){
	all.p<-lapply(1:length(aa.pgs),function(i){
		ni<-names(aa.pgs)[i]
		gp<-generate.pheno(i)
		colnames(gp)<-paste0(colnames(gp),'.',ni)
		gp
	})
	all.p.dt<-do.call(cbind,all.p)
}

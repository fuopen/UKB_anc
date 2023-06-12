library(parallel)

column.map<-function(p1,h1,p2,h2){
    #p:{0,1} h:{0,1}
    col<-sum(c(p1,h1,p2,h2)*2^(3:0))+1
    return(col)
}

w<-function(p1,p2,h1,h2) column.map(p1,h1,p2,h2)

pop1.geno<-function(prob16){#Don't use this function before mean-center genotype
    pop1.1.col<-c(w(0,0,1,0),w(0,1,1,0),w(0,1,1,1),w(0,0,0,1),w(1,0,0,1),w(1,0,1,1))
    pop1.2.col<-w(0,0,1,1)
    pop1<-rowSums(prob16[,pop1.1.col])+2*prob16[,pop1.2.col]
}

pop2.geno<-function(prob16){#Don't use this function before mean-center genotype
    pop2.1.col<-c(w(1,0,1,0),w(1,1,1,0),w(1,0,1,1),w(0,1,0,1),w(1,1,0,1),w(0,1,1,1))
    pop2.2.col<-w(1,1,1,1)
    pop2<-rowSums(prob16[,pop2.1.col])+2*prob16[,pop2.2.col]
}

pop12.geno<-function(prob16){
	pop11.col<-c(w(0,0,1,0),w(0,0,0,1),w(0,0,1,1),w(0,0,1,1))
	pop12.col<-c(w(0,1,1,0),w(0,1,1,1),w(1,0,0,1),w(1,0,1,1))
	pop11.g<-rowSums(prob16[,pop11.col])
	pop12.g<-rowSums(prob16[,pop12.col])
	pop.g12<-data.frame(P11=pop11.g,P12=pop12.g)
}

pop21.geno<-function(prob16){
	pop22.col<-c(w(1,1,1,0),w(1,1,0,1),w(1,1,1,1),w(1,1,1,1))
	pop21.col<-c(w(1,0,1,0),w(1,0,1,1),w(0,1,0,1),w(0,1,1,1))
	pop22.g<-rowSums(prob16[,pop22.col])
	pop21.g<-rowSums(prob16[,pop21.col])
	pop.g21<-data.frame(P22=pop22.g,P21=pop21.g)
}

anc.prob<-function(prob16){
	prob.0.cols<-c(w(0,0,0,0),w(0,0,0,1),w(0,0,1,0),w(0,0,1,1))
	prob.2.cols<-c(w(1,1,0,0),w(1,1,0,1),w(1,1,1,0),w(1,1,1,1))
	prob.0<-rowSums(prob16[,prob.0.cols])
	prob.2<-rowSums(prob16[,prob.2.cols])
	prob.1<-1-prob.0-prob.2
	prob.1[prob.1<0]<-0
	pz<-data.frame(A11=prob.0,A12=prob.1,A22=prob.2)
	return(pz)
}

read.imp<-function(ADMIXINDFILE,OUTDIR,ADMIXPOP,HAPMIX_MODE,output.dir,mcc=16){#don't use this function without mean-centering
	#input arguments are exactly the same as input configuration file for Hapmix (.par file) except the output.dir, which specified by user 
	if(!dir.exists(output.dir)){
		dir.create(output.dir)
	}
	sample.ind<-read.table(ADMIXINDFILE,as.is=T)
	n.sample<-nrow(sample.ind)-1
 	if(!dir.exists(OUTDIR)){
		stop(paste0('error: directory:',OUTDIR,' does not exist!'))
	}
    get.batch<-function(chr){
        get.anc.geno<-do.call(cbind,lapply(0:n.sample,function(i){
            fd<-read.table(paste0(OUTDIR,'/',ADMIXPOP,'.',HAPMIX_MODE,'.',i,'.',chr),as.is=T,header=F)
            fd.p1<-pop1.geno(fd)
            fd.p2<-pop2.geno(fd)
            gn<-cbind(fd.p1,fd.p2)
        }))
    }
    bgb<-do.call(rbind,mclapply(1:22,get.batch,mc.cores=mcc))
    sp.ids<-sample.ind[,1]
	colnames(bgb)<-paste0(c('pop1-','pop2-'),rep(sp.ids,each=2))
    saveRDS(bgb,paste0(output.dir,'/',basename(ADMIXINDFILE),'.pop1_pop2.rds'))
}
read.imp4<-function(ADMIXINDFILE,OUTDIR,ADMIXPOP,HAPMIX_MODE,output.dir,mcc=16){
	#input arguments are exactly the same as input configuration file for Hapmix (.par file) except the output.dir, which specified by user 
	if(!dir.exists(output.dir)){
		dir.create(output.dir)
	}
	sample.ind<-read.table(ADMIXINDFILE,as.is=T)
	n.sample<-nrow(sample.ind)-1
 	if(!dir.exists(OUTDIR)){
		stop(paste0('error: directory:',OUTDIR,' does not exist!'))
	}
    
	get.batch<-function(chr){
        get.anc.geno<-do.call(cbind,lapply(0:n.sample,function(i){
            fd<-read.table(paste0(OUTDIR,'/',ADMIXPOP,'.',HAPMIX_MODE,'.',i,'.',chr),as.is=T,header=F)
            fd.pop12<-pop12.geno(fd)
            fd.pop21<-pop21.geno(fd)
			raw.geno<-rowSums(fd.pop12)+rowSums(fd.pop21)
			fd.anc<-anc.prob(fd)
            gn<-cbind(fd.pop12,fd.pop21,data.frame(geno=raw.geno),fd.anc)
        }))
    }
    bgb<-do.call(rbind,mclapply(1:22,get.batch,mc.cores=mcc))
    sp.ids<-sample.ind[,1]
	pp.names<-c('pop11','pop12','pop22','pop21')
	for(pp in 1:4){
		pp.sp<-bgb[,8*(0:n.sample)+pp]
		colnames(pp.sp)<-sp.ids
    	saveRDS(pp.sp,paste0(output.dir,'/',basename(ADMIXINDFILE),'.',pp.names[pp],'.rds'))
	}
	rg<-bgb[,8*(0:n.sample)+5]
	colnames(rg)<-sp.ids
	saveRDS(rg,paste0(output.dir,'/',basename(ADMIXINDFILE),'.rawgeno.rds'))
	aa.names<-c('A11','A12','A22')
	for(aa in 6:7){
		aa.sp<-bgb[,8*(0:n.sample)+aa]
		colnames(aa.sp)<-sp.ids
		saveRDS(aa.sp,paste0(output.dir,'/',basename(ADMIXINDFILE),'.',aa.names[aa],'.rds'))
	}
}

read.imp.anno<-function(HAPMIX_DATADIR,mcc=10){
    read.imp.chr.anno<-function(chr){
        chr.anno<-read.table(paste0(HAPMIX_DATADIR,'/chr',chr,'_snpfile'),as.is=T,header=F)
    }
    all.anno<-do.call(rbind,mclapply(1:22,read.imp.chr.anno,mc.cores=mcc))
}

out.geno.dir<-''

if(!exists('s.anno')){
	s.anno<-read.imp.anno()
}

estimate_f<-function(mcc=12){
	AFN<-basename(ADMIXINDFILE)
	geno.info<-readRDS(paste0(out.geno.dir,'/',AFN,'.rawgeno.rds'))
	anc.A11<-readRDS(paste0(out.geno.dir,'/',AFN,'.A11.rds'))
	anc.A12<-readRDS(paste0(out.geno.dir,'/',AFN,'.A12.rds'))
	geno.info<-as.matrix(geno.info)
	anc.A11<-as.matrix(anc.A11)
	anc.A12<-as.matrix(anc.A12)
	geno.info[geno.info>2]=2
	nn<-nrow(anc.A11)
	anc.sum<-anc.A11+anc.A12
	anc.sum[anc.sum<1]=1
	anc.A11<-anc.A11/anc.sum
	anc.A12<-anc.A12/anc.sum
	anc.x<-2*anc.A11+anc.A12
	np<-mclapply(1:nn,function(i)coefficients(lm(geno.info[i,]~anc.x[i,])),mc.cores=mcc)
	np.dt<-do.call(rbind,np)
	f.pop1<-np.dt[,2]+np.dt[,1]/2
	f.pop2<-np.dt[,1]/2
	f.pop1[f.pop1<0]=0
	f.pop1[f.pop1>1]=1
	f.pop2[f.pop2<0]=0
	f.pop2[f.pop2>1]=1
	freq.dt<-s.anno
	freq.dt$pop1.freq<-f.pop1
	freq.dt$pop2.freq<-f.pop2
	freq.dt
}

if(!exists('fq')){
	if(!file.exists(paste0(out.geno.dir,'/pop1_pop2_Effect_allele_freq.rds'))){
		fq<-estimate_f()
		saveRDS(fq,paste0(out.geno.dir,'/pop1_pop2_Effect_allele_freq.rds'))
	}
	else{
		fq<-readRDS(paste0(out.geno.dir,'/pop1_pop2_Effect_allele_freq.rds'))
	}
}

mean.center.geno<-function(){
	AFN<-basename(ADMIXINDFILE)
	geno.info<-readRDS(paste0(out.geno.dir,'/',AFN,'.rawgeno.rds'))
	anc.A11<-readRDS(paste0(out.geno.dir,'/',AFN,'.A11.rds'))
	anc.A12<-readRDS(paste0(out.geno.dir,'/',AFN,'.A12.rds'))
	geno.info<-as.matrix(geno.info)
	anc.A11<-as.matrix(anc.A11)
	anc.A12<-as.matrix(anc.A12)
	geno.info[geno.info>2]=2
	nn<-nrow(anc.A11)
	anc.sum<-anc.A11+anc.A12
	anc.sum[anc.sum<1]=1
	anc.A11<-anc.A11/anc.sum
	anc.A12<-anc.A12/anc.sum
	anc.A22<-1-anc.A11-anc.A12
	geno.P11<-as.matrix(readRDS(paste0(out.geno.dir,'/',AFN,'.pop11.rds')))
	geno.P12<-as.matrix(readRDS(paste0(out.geno.dir,'/',AFN,'.pop12.rds')))
	geno.P21<-as.matrix(readRDS(paste0(out.geno.dir,'/',AFN,'.pop21.rds')))
	geno.P22<-as.matrix(readRDS(paste0(out.geno.dir,'/',AFN,'.pop22.rds')))
		
	geno.P1<-geno.P11+geno.P12-(2*fq$pop1.freq*anc.A11+fq$pop1.freq*anc.A12)
	geno.P2<-geno.P22+geno.P21-(2*fq$pop2.freq*anc.A22+fq$pop2.freq*anc.A12)
	
	saveRDS(geno.P1,paste0(out.geno.dir,'/',AFN,'.pop1.rds'))
	saveRDS(geno.P2,paste0(out.geno.dir,'/',AFN,'.pop2.rds'))
}

cal.pgs<-function(betafile){
	betas<-readRDS(betafile)
	betas.in<-intersect(betas$SNP,fq$SNP)
	beta<-betas[match(betas.in,betas$SNP),'BETA']
	fq.ids<-match(betas.in,fq$SNP)
	g.p1<-readRDS(paste0(out.geno.dir,'/',AFN,'.pop1.rds'))[fq.ids,]
	g.p2<-readRDS(paste0(out.geno.dir,'/',AFN,'.pop2.rds'))[fq.ids,]
	pgs.p1<-t(g.p1)*matrix(beta,ncol=1)
	pgs.p2<-t(g.p2)*matrix(beta,ncol=1)
	pgs<-data.frame(pgs.pop1=pgs.p1,pgs.pop2=pgs.p2)
}

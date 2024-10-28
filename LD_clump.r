##***********************************************************************
#***********************************************************************
#Using roundrect function in package shape
library(shape)
library(colorRamps)
library(parallel)

pthr<--log10(5e-8)

res.dir<-'/data/muscovy/not-backed-up/shu/'
var.anno.file<-'new_gwas_snps_anno_withaf.rds'

clump.dir<-'PS_hapmap3_r2_0.1/'

clump.dir.eur<-paste0(clump.dir,'/eur_real/')

if(!dir.exists(clump.dir)){
    dir.create(clump.dir)
}

#clump.dir.eur<-paste0(clump.dir,'simu_bgen')

if(!dir.exists(clump.dir.eur)){
    dir.create(paste0(clump.dir.eur))
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

if(!exists('hapmap3.snps')){
    hapmap3.snps<-readRDS('hapmap3_subset_SNPs_PS_filter_mean0.rds')
}

if(!exists('hap.var.ids')){
    hap.var.ids<-ac.anno$variant %in% hapmap3.snps$vid
}
gwas.dir<-'/homes/shu/workdir/ng_revision/GWAS_bgenie_all_need/'

if(!exists('fp')){
	fp<-read.table('gwas_files_path.txt',as.is=T,sep='\t')
}

clump.by.beta<-function(tmp.dir=clump.dir.eur){
    if(!dir.exists(tmp.dir)){
        dir.create(tmp.dir)
    }
	codes<-fp[[1]]
	files<-fp[[3]]
	each.code<-function(i){
		code<-codes[i]
		file<-files[i]
    	code.dir<-paste0(tmp.dir,'/',code)
    	if(!dir.exists(code.dir)){
        	dir.create(code.dir)
    	}
		gs<-readRDS(file)
		gs<-data.frame(gs)
		gs.hapmap3<-gs[hap.var.ids,]
		gs.hapmap3[[4]]<-10^-(gs.hapmap3[[4]])
		gs.hapmap3$ids<-ac.anno$variant[hap.var.ids]
		gs.hapmap3$chr<-ac.anno$ichr[hap.var.ids]
		gs.hapmap3<-gs.hapmap3[,c('ids','chr','ac.beta','ac.pv')]
		colnames(gs.hapmap3)<-c('VID','CHR','BAC','PAC')		

    	ltff<-function(i){
			chr.gs<-subset(gs.hapmap3,CHR==i)
			write.table(chr.gs,paste0(code.dir,'/chr',i,'.h3.pv'),row.names=F,col.names=T,quote=F)
    	}
		lt<-lapply(1:22,ltff)
	}
	m1<-mclapply(1:length(codes),each.code,mc.cores=12)
}

run.clump.cmd<-function(code,p1=0.05,p2=1,r2=0.1,kb=500,pop='euro'){
    run.cmd<-function(chr,ap='PAC'){
        cmd<-paste0('/data/mallard/shu/apps/bin/plink --bfile /data-tmp/muscovy/not-backed-up/shu/1kg/',pop,'_ld_bed_bychr/',pop,'_ld_chr',chr,' --clump ',clump.dir.eur,'/"',code,'"/chr',chr,'.h3.pv --clump-p1 ',p1,' --clump-p2 ',p2,' --clump-kb ',kb,' --clump-r2 ',r2,' --clump-snp-field VID --clump-field ',ap,' --out ',clump.dir.eur,'/"',code,'"/chr',chr,'.',ap,'_p1_',p1,'.clump')
		st<-system(cmd,intern=T)
		print(cmd)
    }
    lapply(1:22,run.cmd,ap='PAC')
}

run.clump<-function(code,pv=0.05){
    #clump.by.beta()
    run.clump.cmd(code,p1=pv,pop='euro')
}

clump.all.gwas<-function(pv=0.05){
	#codes<-basename(list.dirs(gwas.dir))
	#codes<-codes[grep('^f\\.',codes)]
	codes<-fp[[1]]
	mclapply(codes,run.clump,pv=pv,mc.cores=12)
}	

assemble.chr.clump<-function(code,chr,p,comp='AC'){
    ac.file<-paste0(clump.dir.eur,'/',code,'/chr',chr,'.P',comp,'_p1_',p,'.clump.clumped')
    if(!file.exists(ac.file)){
        ta<-NULL
    }
    else{
        ta<-read.table(ac.file,as.is=T,header=T)$SNP
    }
}

assemble.chr<-function(chr,p,dir='new_clump_snps_bg3_r_0.1/'){
    if(!dir.exists(dir)){
        dir.create(dir)
    }
	codes<-fp[[1]]
    chr.mc.ac<-mclapply(codes,function(x)try(assemble.chr.clump(x,chr=chr,p=p)),mc.cores=12)
	snps<-unique(do.call(c,chr.mc.ac))
    
    snp.match.ac<-sapply(chr.mc.ac,function(x)!is.na(match(snps,x)))    
    colnames(snp.match.ac)<-paste0(codes,'.AC')
    snp.match<-snp.match.ac
    rownames(snp.match)<-snps
    saveRDS(snp.match,paste0(dir,'/chr',chr,'_pv_',p,'.rds'))
}

merge_and_beta<-function(p){
    clump.snp.dir<-'new_clump_snps_bg3_r_0.1/'
	codes<-fp[[1]]
	each.chr<-function(chr){
		snp.chr<-readRDS(paste0(clump.snp.dir,'chr',chr,'_pv_',p,'.rds'))
		code.chr<-function(code){
			#gs<-readRDS(paste0(gwas.dir,code,'.rds'))
			gsf<-fp[fp[[1]]==code,3]
			gs<-readRDS(gsf)
			gs<-data.frame(gs)
			gsk<-gs[match(rownames(snp.chr),ac.anno$variant),'ac.beta']
		}
		all.codes<-gsub('\\.AC$','',colnames(snp.chr))
		beta.all<-do.call(cbind,mclapply(all.codes,code.chr,mc.cores=12))
		chr.beta<-snp.chr*beta.all
		saveRDS(chr.beta,paste0(clump.snp.dir,'betas_chr',chr,'.rds'))
	}
	betas<-lapply(1:22,each.chr)
}

read.chrs<-function(pv){
    clump.snp.dir<-'new_clump_snps_bg3_r_0.1/'
	pv.list<-lapply(1:22,function(chr){
		rr<-readRDS(paste0(clump.snp.dir,'chr',chr,'_pv_',pv,'.rds'))
	})
}
pvs<-c(0.05,0.01,1e-5,5e-8)

if(!exists('pv.traits')){
	pv.traits<-lapply(pvs,read.chrs)
	names(pv.traits)<-pvs
}

get.beta<-function(code){
    clump.snp.dir<-'new_clump_snps_bg3_r_0.1/'
	gsf<-fp[fp[[1]]==code,3]
	gs<-readRDS(gsf)
	each.pv<-function(p){
		each.chr<-function(chr){
			snp.chr<-pv.traits[[as.character(p)]][[chr]]
			snps<-rownames(snp.chr)[snp.chr[,paste0(code,'.AC')]]
			gsk.chr<-gs[match(snps,ac.anno$variant),'ac.beta',drop=F]
			rownames(gsk.chr)<-snps
			gsk.chr<-as.matrix(gsk.chr)
		}
		all.betas<-lapply(1:22,each.chr)
	}
	all.pvs<-lapply(pvs,each.pv)
	names(all.pvs)<-pvs
	all.pvs	
}

all.betas<-function(){
	codes<-fp[[1]]
	gab<-mclapply(codes,get.beta,mc.cores=12)
	names(gab)<-codes
	gab
}

all.betas.pv<-function(gab){
	pab<-lapply(as.character(pvs),function(x)lapply(gab,function(gg)gg[[x]]))
	names(pab)<-pvs
	pab
}	

run.all<-function(p=0.05){
    clump.all.gwas(pv=p)
    #lapply(1:22,assemble.chr,codes=codes,p=p)
    lapply(1:22,assemble.chr,p=p)
    merge_and_beta(p=p)
}

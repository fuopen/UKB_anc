library(parallel)
args<-commandArgs(T)
#args<-'test2.txt'
#snp.file<-paste0(dir,args[1])
#geno.file<-paste0(snp.file,'.geno')

ancestry.file<-args[1]#NxA Matrix with row samples and columns ancestry inferred by ancestry pipeline 
geno.file<-args[2]#NxM Matrix with row samples and columns hard-call genotypes (coded as 0,1,2 and the row of ancestry and genotype should match 
snp.file<-args[3] #SNP annotation file, SNP ID columns should be called "SNP", and the SNP ID columns should match the columns of genotype matrix 
out.dir<-args[4] #Output directory
#This software doesn't need ids for the samples
mcc<-as.integer(args[5])

if(!exists('ancestry')){
    sample.markers<-as.matrix(read.table(geno.file,as.is=T,header=F)) #(ns by nm genotype matrix)
    marker.ids<-read.table(snp.file,as.is=T,header=F)$SNP#colnames(sample.markers)
    sample.ids<-paste0('ukb',1:nrow(sample.markers))
    rownames(sample.markers)<-sample.ids

    ns<-length(sample.ids) ### number of samples
    nm<-length(marker.ids) ### number of markers
    new.markers=t(cbind(floor(sample.markers/2),floor((sample.markers+1)/2)))
    ancestry<-as.matrix(read.table(ancestry.file,as.is=T,header=T))
    new.ancestry=as.data.frame(rbind(ancestry,ancestry))
    na<-ncol(ancestry) ### number of ancestries
    prob_m<-matrix(rep(0,2*ns*na),ncol=na)
}

set.initial.freq<-function(i,del=0.001){### ith site
    freq.lm<-lm(new.markers[,i] ~ .,data=new.ancestry)
    freq.ini<-coef(freq.lm)
    intcpt<-freq.ini[1]
    freq.ini<-freq.ini[-1]+intcpt
    #freq.ini<-ifelse(freq.ini<0 | freq.ini>1,0.001,freq.ini)
    freq.ini[is.na(freq.ini)]<-0
    freq.ini[freq.ini<0]<-0+del
    freq.ini[freq.ini>1]<-1-del
    return(freq.ini)
}

em<-function(m,tol=0.001,minits=3){############### mth marker
    time.record<-matrix(rep(0,3000),ncol=3)
    step.record<-rep(0,1000)
    #print(m)
    #f0<-rep(0.5,na) #initial frequency 0.5 for each of na groups
    tl<-proc.time()
    f0<-set.initial.freq(m)
    tr<-proc.time()
    init.time<-(tr-tl)[1:3]
    f1<-f0
    geno<-new.markers[,m]
    delta<-10 ## initialize difference
    k=0
    ######
    # probablity for each haplotype SNP f^g*(1-f)^g can be rewritten as
    # 1-f+g*(2*f-1)
    # Ancestry*(1-f+g*(2*f-1))=(A*(1-f)+A*(2*f-1)*f)=C1+C2*f
    ######
    C1<-new.ancestry*(1-geno)
    tC2<-t(new.ancestry*(2*geno-1)) ##precalcutaed
    setting.time<-(proc.time()-tr)[1:3]
    while(delta>tol | k<minits){
        tl<-proc.time()
        k=k+1
        f0<-f1
        prob_m=C1+t(tC2*f0)
        prob_m=prob_m/rowSums(prob_m)
        term1=colSums(prob_m)
        term2=colSums(prob_m[geno==1,])
        f1=term2/term1
        delta<-sqrt(sum((f1-f0)^2))
        tr<-proc.time()
        time.record[k,]<-(tr-tl)[1:3]
        step.record[k]<-delta
    }
    rd<-list(inti.t=init.time,set.t=setting.time,time.r=time.record,step.r=step.record)
    return(list(f1,rd))
}

f.em<-mclapply(1:nm,em,mc.cores=mcc)

f.freq<-do.call(rbind,mclapply(f.em,function(x)x[[1]],mc.cores=mcc))
f.rds<-mclapply(f.em,function(x)x[[2]],mc.cores=mcc)

freq.gz<-gzfile(paste0(out.dir,'/',snp.file,'.freq.tab.gz'),'w')
write.table(f.freq,freq.gz,row.names=F,col.names=F,quote=F,sep='\t')
close(freq.gz)
saveRDS(f.rds,paste0(out.dir,'/',snp.file,'.rds'))
q(save='no')

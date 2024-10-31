library(parallel)
run<-function(){
    file.seq<-seq(1,8500,50)
    dir<-'imp_geno/'
    my.files<-paste0(dir,'batch_',file.seq,'_',file.seq+49,'/data/sample_ind')
    me<-do.call(rbind,mclapply(my.files,function(x)read.table(x,as.is=T,header=F),mc.cores=17))
    trio<-read.table('ind_geno_relatedness/identified_trios.txt',as.is=T,header=T)
    trio.ids<-as.vector(t(trio[,-1]))
    file.seq.2<-file.seq+49
    sse<-sapply(match(trio.ids,me[,1]),function(x)which(apply(cbind(file.seq,file.seq.2),1,function(y)y[1] <=x && y[2]>=x)))
    sse.files<-my.files[sse]
    sse.pos<-match(trio.ids,me[,1])%%50
    return(list(sse.files,sse.pos))
}

get.geno<-function(){
    rg<-run()
    rg.files<-rg[[1]]
    rg.pos<-rg[[2]]
    n<-length(rg.files)
    read.each<-function(i){
        geno.dir<-gsub('sample_ind','',rg.files[i])
        geno.files<-paste0(geno.dir,'chr',1:22,'_geno.ORIG')
        geno.data<-do.call('c',mclapply(1:22,function(chr){
            print(paste0('reading ',geno.files[chr]))
            print(paste0('pos ',rg.pos[i]))
            geno.char<-readLines(geno.files[chr])
            geno.str<-sapply(strsplit(geno.char,''),function(x)x[[rg.pos[i]]])
            return(geno.str)
        },mc.cores=14))
        return(geno.data)
    }
    sr<-sapply(1:n,read.each)
    return(sr)
}

if(!exists('gg')){
    gg<-readRDS('trio_geno/trio_6_ind_original_genotype.rds')
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

generate.check<-function(){#order child,father,mother
    my.dt<-rbind(
        c(0,0,0),
        c(0,0,1),
        c(0,0,2),
        c(0,1,0),
        c(0,1,1),
        c(0,1,2),
        c(0,2,0),
        c(0,2,1),
        c(0,2,2),
        c(1,0,0),
        c(1,0,1),
        c(1,0,2),
        c(1,1,0),
        c(1,1,1),
        c(1,1,2),
        c(1,2,0),
        c(1,2,1),
        c(1,2,2),
        c(2,0,0),
        c(2,0,1),
        c(2,0,2),
        c(2,1,0),
        c(2,1,1),
        c(2,1,2),
        c(2,2,0),
        c(2,2,1),
        c(2,2,2))
    return(my.dt)
}

if(!exists('error_geno')){
    error_geno<-c('002','012','020','021','100','122','200','201','202','210','220')
}

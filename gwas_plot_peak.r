library(parallel)

if(!exists('tp')){
	tp<-readRDS('data/98_traits_for_ldsc_EA.rds')
	tp.ldsc<-attr(tp,'ldsc')
}

pthr=-log10(5e-8)

plot.phenos<-function(dir.out='new_traits_peak/'){
    fids<-names(tp)[15:16]
    dir<-dir.out
    if(!dir.exists(dir)){
        dir.create(dir)
    }
    plot.peak<-function(dt.sig,main=NA){
        par(plt=c(0.15,0.85,0.15,0.85))
        dt.sig<-dt.sig[!is.na(dt.sig$wb_sample_copy) & dt.sig$wb_sample_copy>=100,]
        if(nrow(dt.sig)==0){
            return(NULL)
        }
        dt.pch<-rep(19,nrow(dt.sig))
        dt.col<-rep('cyan',nrow(dt.sig))
        plot(dt.sig$ac.pv,dt.sig$pc.pv,xlab=expression("-log"[10]*"(pvalue"[AC]*")"),ylab=expression("-log"[10]*"(pvalue"[PC]*")"),xlim=c(0,35),ylim=c(0,35),xaxs='i',yaxs='i',pch=dt.pch,cex=0.8,col=dt.col,main=main)
        abline(a=0,b=1,lwd=0.8,col='lightblue')
        abline(h=pthr,lwd=0.8,col='darkblue',lty=3)
        abline(v=pthr,lwd=0.8,col='darkblue',lty=3)
    }
    mclapply(1:length(fids),function(i){
        peak.file<-paste0(dir,fids[i],'.peak.pdf')
        pdf(peak.file)
		
        plot.peak(tp[[i]],fids[i])
        dev.off()
    },mc.cores=8)
}

library(colorRamps)

if(!exists('eu.mats')){
    eu.mats<-readRDS('data/combinned_bin_simu_and_real_traits53_pv0.05_pv1e4.rds')
}

plot.XY<-function(pfile='figure4b_update.pdf'){
    eq<-eu.mats[nrow(eu.mats):1,]
    ex<-seq(0.1,1.5,0.2)
    #ex<-c(seq(0.1,0.9,0.2),1,1.3,1.4)
    dot.pos<-seq(0.07,0.13,length=4)-0.1
    ex.pos<-rbind(ex+dot.pos[1],ex+dot.pos[2],ex+dot.pos[3],ex+dot.pos[4])
    rec.cols<-blue2red(30)[seq(1,30,5)][-c(1,6)]
    pdf(pfile)
    par(mar=c(5,5,3,3),xaxs='i',yaxs='i')
    plot(0,0,type='n',xlim=c(0,1.75),ylim=c(0,1.25),xlab=NA,ylab=NA,main=NA,axes=F)
    #clip(0,1.7,0,1.25)
    segments(0,0.1,0.2,0.1,lwd=1.2,lty=2)
    segments(0.2,0.3,0.4,0.3,lwd=1.2,lty=2)
    segments(0.4,0.5,0.6,0.5,lwd=1.2,lty=2)
    segments(0.6,0.7,0.8,0.7,lwd=1.2,lty=2)
    segments(0.8,0.9,1.0,0.9,lwd=1.2,lty=2)
    segments(1.0,1.0,1.6,1.0,lwd=1.2,lty=2)
    #abline(v=c(seq(0.1,0.9,0.2),1,1.3,1.4),lwd=0.1,col='#00000080')
    #abline(a=0,b=1,lwd=0.8,col='darkgreen')
    for(i in 1:4){
        points(ex.pos[i,],eq[,3*(i-1)+1],pch=19,col=paste0(rec.cols[i],'C0'),cex=0.9)
        arrows(ex.pos[i,],eq[,3*(i-1)+2],ex.pos[i,],eq[,3*(i-1)+3],angle=90,code=3,length=0.02,lwd=1.0,col=paste0(rec.cols[i],'C0'))
    }
    axis(1,at=seq(0.1,1.5,0.2),labels=c(ex[1:5],1,NA,NA),cex.lab=0.5,pos=0)
    axis(2,at=seq(0,1.2,0.1),cex.lab=0.5,las=2,line=1,pos=0)
    #legend('topleft',legend=c('[0.1,0.35]','[0.35,0.55]','[0.55,0.78]','[0.78,95]'),pch=19,col=paste0(rec.cols,'C0'),cex=1.5,bty='n')
    dev.off()
}


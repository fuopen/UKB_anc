library(colorRamps)

if(!exists('mean.af')){
    mean.af<-readRDS('data/mean_af_anc_bins.rds')
}

if(!exists('plot.euaf.list')){
    plot.euaf.list<-readRDS('data/simu_24_traits_ratio_bins.rds')
}

if(!exists('bound.beta')){
    bound.beta<-cbind(c(1,0.9995,1.0005),do.call(cbind,lapply(plot.euaf.list[6:1],function(x)x$AF[,5])))
    bound.beta<-as.data.frame(bound.beta)
    names(bound.beta)[-1]<-paste0('rho_',names(bound.beta)[-1])
}

plot.eu.af.panel<-function(rho=0.5,af.anc=0.09,shift.y=0,ifreal=F){
    plot.euaf<-plot.euaf.list[[as.character(rho)]]
    ee.beta<-plot.euaf$EU[,-5]+shift.y
    aa.beta<-plot.euaf$AF[,-5]+shift.y

    eu.bd.beta<-bound.beta[,1]+shift.y
    eu.bd.beta.exp<-bound.beta[,1]*as.numeric(rho)+shift.y
    af.bd.beta<-bound.beta[[paste0('rho_',rho)]]+shift.y

    eu.beta.pos<-mean.af[-6]
    af.beta.pos<-mean.af[-1]

    ee.betas<-cbind(eu.bd.beta,ee.beta)
    aa.betas<-cbind(aa.beta,af.bd.beta)

    ee.betas[ee.betas<0]=0
    aa.betas[aa.betas<0]=0

    #plot(0,0,type='n',xlim=c(-0.03,1.03),ylim=c(-0.1,1.2),xlab=NA,ylab=NA,xaxt='n',pch=19,col='blue',bty='n',yaxt='n')
    #axis(2,at=seq(0,0.4,0.08)+shift.y,labels=seq(0,0.4,0.08),las=2,cex.lab=0.3,line=-1)
    axis(2,at=seq(0,1.2,0.2)+shift.y,labels=c(seq(0,1.0,0.2),NA),las=2,cex.lab=0.3,line=-1)
    rect(0,ee.betas[2,1],1,ee.betas[3,1],border=NA,col='#0000FFA0')
    segments(0,eu.bd.beta.exp[1],1,eu.bd.beta.exp[1],lwd=0.8,lty=3,col='darkblue')
    blk.bx<-c(0,0.1,0.35,0.55,0.78,0.95)
    blk.tx<-c(0.1,0.35,0.55,0.78,0.95,1)
    blk.by<-rep(0,6)+shift.y
    blk.ty<-rep(1.00,6)+shift.y
    rec.cols<-blue2red(30)[seq(1,30,5)]
    for(i in 1:6){
        rect(blk.bx[i],blk.by[i],blk.tx[i],blk.ty[i],col=paste0(rec.cols,'80')[i],border=NA)
    }
    points(eu.beta.pos,ee.betas[1,],pch=19,cex=1,col='blue',bty='n',yaxt='n')
    arrows(eu.beta.pos,ee.betas[2,],eu.beta.pos,ee.betas[3,],lwd=0.8,angle=90,code=3,length=0.02)
    points(af.beta.pos,aa.betas[1,],pch=19,cex=1,col='red')
    arrows(af.beta.pos,aa.betas[2,],af.beta.pos,aa.betas[3,],lwd=0.8,angle=90,code=3,length=0.02)
    #axis(1,at=round(mean.af,2))
    #axis(3,at=c(0,0.1,0.35,0.55,0.78,0.95,1),line=-1)
    #legend('topleft',legend=c('EU','AF'),col=c('blue','red'),pch=19,cex=0.8,bty='n')
}

plot.all.panels<-function(pfile='EDF9.pdf',af.anc=0.09){
    pdf(pfile,height=10)
    plot(0,0,type='n',xlim=c(-0.03,1.03),ylim=c(-0.1,9.0),xlab=NA,ylab=NA,xaxt='n',pch=19,col='blue',bty='n',yaxt='n')
    segments(mean.af,0,mean.af,2.5,lty=2,lwd=0.2)
    rhos<-c('0.1','0.3','0.5','0.7','0.9','1.0')
    shifts<-cumsum(c(0,rep(1.35,6)))
    for(i in 1:length(rhos)){
        plot.eu.af.panel(rhos[i],shift.y=shifts[i])
    }
    axis(1,at=round(mean.af,2),line=-1)
    axis(3,at=c(0,0.1,0.35,0.55,0.78,0.95,1),line=-5)
    dev.off()
}

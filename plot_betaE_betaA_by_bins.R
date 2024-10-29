library(colorRamps)

if(!exists('mean.af')){
    mean.af<-readRDS('mean_af_anc_bins.rds')
}


if(!exists('plot.euaf.list')){
    plot.euaf.list<-readRDS('bin.ratio.rds')
}

if(!exists('bound.beta')){
    bound.beta<-cbind(c(1,0.9995,1.0005),do.call(cbind,lapply(plot.euaf.list,function(x)x$AF[,5])))
    bound.beta<-as.data.frame(bound.beta)
    names(bound.beta)[-1]<-c('withmask.mean.center','withmask.non.mean.center','withoutmask.mean.center','withoutmask.non.mean.center')
}

plot.eu.af.panel<-function(cv='mask.mc',shift.y,center.type=NULL){
    plot.euaf<-plot.euaf.list[[cv]]
    ee.beta<-plot.euaf$EU[,-5]+shift.y
    aa.beta<-plot.euaf$AF+shift.y
    #ee.beta<-plot.euaf$EU+shift.y
    #aa.beta<-plot.euaf$AF+shift.y

    eu.bd.beta<-bound.beta[,1]+shift.y
    af.bd.beta<-bound.beta[[cv]]+shift.y

    eu.beta.pos<-mean.af[-6]
    #af.beta.pos<-mean.af[-1]
    #eu.beta.pos<-mean.af
    af.beta.pos<-mean.af[-1]

    ee.betas<-cbind(eu.bd.beta,ee.beta)
    #aa.betas<-cbind(aa.beta,af.bd.beta)
    aa.betas<-aa.beta
    ee.betas[ee.betas<0]=0
    aa.betas[aa.betas<0]=0

    axis(2,at=seq(0,1.2,0.2)+shift.y,labels=c(seq(0,1.0,0.2),NA),las=2,cex.lab=0.3,line=-1)
    if(!is.null(center.type)){
        mtext(text=center.type,side=2,at=0.5+shift.y,line=2)
    }
    rect(0,ee.betas[2,1],1,ee.betas[3,1],border=NA,col='#0000FFA0')
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
}

plot.all.panels<-function(){
    pdf('ratio.bin.pdf',height=10,width=18)
    layout(matrix(1:2,ncol=2))
    plot(0,0,type='n',xlim=c(-0.03,1.03),ylim=c(-0.1,2.7),xlab=NA,ylab=NA,xaxt='n',pch=19,col='blue',bty='n',yaxt='n')
    title('With mask')
    segments(mean.af,0,mean.af,2.5,lty=2,lwd=0.2)
    shifts<-cumsum(c(0,rep(1.45,2)))
    side.titles=c('Without mean center','With mean center')
    plot.euaf.list<-plot.euaf.list[c(2,1,4,3)]
    for(i in 1:2){
        plot.eu.af.panel(names(plot.euaf.list)[i],shift.y=shifts[i],side.titles[i])
    }
    axis(1,at=round(mean.af,2),line=-1)
    axis(3,at=c(0,0.1,0.35,0.55,0.78,0.95,1),line=-4)
    plot(0,0,type='n',xlim=c(-0.03,1.03),ylim=c(-0.1,2.7),xlab=NA,ylab=NA,xaxt='n',pch=19,col='blue',bty='n',yaxt='n')
    title('Without mask')
    segments(mean.af,0,mean.af,2.5,lty=2,lwd=0.2)
    shifts<-cumsum(c(0,rep(1.45,2)))
    for(i in 1:2){
        plot.eu.af.panel(names(plot.euaf.list)[i+2],shift.y=shifts[i])
    }
    axis(1,at=round(mean.af,2),line=-1)
    axis(3,at=c(0,0.1,0.35,0.55,0.78,0.95,1),line=-4)
    dev.off()
}

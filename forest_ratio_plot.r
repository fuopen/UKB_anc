library(parallel)


if(!exists('dt.plot')){
    dt.plot<-readRDS('data/figure5_plot_ratio_dt.rds')
    dt.plot[1:2,8]<-dt.plot[1:2,8]-0.015
    dt.plot[1:2,9]<-dt.plot[1:2,9]+0.015
}
if(!exists('anno.plot')){
    anno.plot<-readRDS('data/figure5_plot_ratio_anno.rds')
}

plot.multiple.ratio.mainfig.final<-function(pfile="figure5_ratio_plot.pdf"){
    pdf(pfile,width=15,height=9)
    plot.single.panel<-function(shift,ci.mat,col.p='purple'){
        cimat<-as.matrix(ci.mat)+shift
        ni<-nrow(cimat)
        yy<-seq(1,5*ni,5)
        yy[1]<-yy[1]-1
        yy[1:2]<-yy[1:2]-2
        yy[3]<-yy[3]+1
        yy[4:ni]<-yy[4:ni]+2

        yy[14]<-yy[14]+2
        yy[15:ni]<-yy[15:ni]+3
        clip(-2,13,-15.5,5*(ni+1))
        abline(v=0.5+shift,lwd=1.2,lty=3,col='darkgrey')
        clip(-2,13,-15.5,5*(ni+1))
        abline(v=1+shift,lwd=1.2,lty=3,col='red4')
        clip(-2,13,-15.5,5*(ni+1))
        abline(v=1.5+shift,lwd=1.2,lty=3,col='darkgrey')
        points(cimat[,1],yy,cex=0.8,col=col.p,pch=19)
        arrows(cimat[,2],yy,cimat[,3],yy,angle=90,code=3,length=0.01,lwd=0.7,col=col.p)
        axis(1,at=seq(0,2,0.5)+shift,labels=seq(0,2,0.5),line=1,cex.lab=2.5)
    }
    
    plot(0,0,type='n',xlim=c(-3.5,9.0),ylim=c(-1,(nrow(dt.plot))*5),axes=F,ann=F)
    plot.single.panel(0,dt.plot[,1:3],'blue')
    plot.single.panel(3.2,dt.plot[,4:6])
    plot.single.panel(6.4,dt.plot[,7:9],'red')
    ndt<-nrow(dt.plot)
    yt<-seq(1,5*ndt,5)
    yt.m<-yt
    yt.m[1]<-yt.m[1]-1
    yt.m[1:2]<-yt.m[1:2]-2
    yt.m[3]<-yt.m[3]+1
    yt.m[4:ndt]<-yt.m[4:ndt]+2
    yt.m[14]<-yt.m[14]+1
    yt.m[15:ndt]<-yt.m[15:ndt]+3
    text(0,yt.m,anno.plot[,2],cex=1.2,pos=2)
    dev.off()
}

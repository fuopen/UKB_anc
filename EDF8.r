if(!exists('mask.dt.mc')){
    mask.dt.mc<-matrix(c(1.02809606862402,0.988406648944784,1.06880568126796,1.05031182587167,0.982351528138807,1.13441833461341,0.575862443786318,0.5587621812503,0.592533138563997),nrow=1)
    colnames(mask.dt.mc)<-c('beta_eu_all','beta_eu_all.025','beta_eu_all.975','beta_eu_af','beta_eu_af.025','beta_eu_af.975','beta_af_all','beta_af_all.025','beta_af_all.975')
}

if(!exists('mask.dt.nmc')){
    mask.dt.nmc<-matrix(c(0.708563999252722,0.673396142466918,0.744380634515106,0.652147067595492,0.606960374901992,0.71092079431198,0.377384065081166,0.363915746814771,0.390441567717442),nrow=1)
    colnames(mask.dt.nmc)<-c('beta_eu_all','beta_eu_all.025','beta_eu_all.975','beta_eu_af','beta_eu_af.025','beta_eu_af.975','beta_af_all','beta_af_all.025','beta_af_all.975')
}

if(!exists('unmask.dt.mc')){
    unmask.dt.mc<-matrix(c(1.0163439,0.9803494,1.0520222,1.0281554,0.9722360,1.0925679,0.5743359,0.5586898,0.5895154),nrow=1)
    colnames(unmask.dt.mc)<-c('beta_eu_all','beta_eu_all.025','beta_eu_all.975','beta_eu_af','beta_eu_af.025','beta_eu_af.975','beta_af_all','beta_af_all.025','beta_af_all.975')
}
if(!exists('unmask.dt.nmc')){
    unmask.dt.nmc<-matrix(c(0.9028014,0.8716939,0.9327526,0.8273400,0.7894253,0.8690388,0.5959996,0.5810461,0.6097352),nrow=1)
    colnames(unmask.dt.nmc)<-c('beta_eu_all','beta_eu_all.025','beta_eu_all.975','beta_eu_af','beta_eu_af.025','beta_eu_af.975','beta_af_all','beta_af_all.025','beta_af_all.975')
}

plot.res<-function(pfile='ed8.pdf'){
	pdf(pfile)
    plot(0,0,type='n',xlab=NA,ylab=NA,xlim=c(0,7),ylim=c(0.45,1.15),xaxt='n')
    abline(h=1,lwd=0.3,col='darkgreen')
    segments(3.5,0.55,3.5,1.20,lwd=0.8,col='darkgreen')
    legend('bottom',c('mean-centering,with masking','mean-centering,without masking','non-mean-centering, with masking','non-mean-centering, without masking'),col=c('blue','blue','cyan','cyan'),pch=c(19,18,19,18),cex=1.0,bty='n')
    points(c(1,5),mask.dt.mc[1,c(1,4)],pch=19,col='blue',cex=2.4)
    arrows(c(1,5),mask.dt.mc[1,c(2,5)],c(1,5),mask.dt.mc[1,c(3,6)],lwd=1.2,angle=90,code=3,length=0.05)
    points(c(2,6),mask.dt.nmc[1,c(1,4)],pch=19,col='cyan',cex=2.4)
    arrows(c(2,6),mask.dt.nmc[1,c(2,5)],c(2,6),mask.dt.nmc[1,c(3,6)],lwd=1.2,angle=90,code=3,length=0.05)
    points(c(1.5,5.5),unmask.dt.mc[1,c(1,4)],pch=18,col='blue',cex=2.4)
    arrows(c(1.5,5.5),unmask.dt.mc[1,c(2,5)],c(1.5,5.5),unmask.dt.mc[1,c(3,6)],lwd=1.2,angle=90,code=3,length=0.05)
    points(c(2.5,6.5),unmask.dt.nmc[1,c(1,4)],pch=18,col='cyan',cex=2.4)
    arrows(c(2.5,6.5),unmask.dt.nmc[1,c(2,5)],c(2.5,6.5),unmask.dt.nmc[1,c(3,6)],lwd=1.2,angle=90,code=3,length=0.05)
	dev.off()
}

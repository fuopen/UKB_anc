source('compare_res_trio3.r')

gen.bar.hap.mat<-function(chr,id,hap.id=1,shift=0){
    imp.chr.anno<-imp.anno[imp.chr[[chr]],]
    n.s<-nrow(imp.chr.anno)
    imp.width<-diff(imp.chr.anno[,4])
    imp.eu.xb<-c(imp.chr.anno[1,4],imp.chr.anno[2:n.s,4]-imp.width/2)
    imp.eu.xt<-c(imp.chr.anno[1:(n.s-1),4]+imp.width/2,imp.chr.anno[n.s,4])
    imp.eu.yb<-rep(0,length(imp.eu.xb))+shift
    imp.eu.yt<-hap.anc[[2]][[id]][imp.chr[[chr]],hap.id*2-1]+shift
    imp.af.xb<-imp.eu.xb
    imp.af.xt<-imp.eu.xt
    imp.af.yb<-imp.eu.yt
    imp.af.yt<-rep(1,length(imp.eu.yb))+shift

    hap.eu<-cbind(imp.eu.xb,imp.eu.yb,imp.eu.xt,imp.eu.yt)
    hap.af<-cbind(imp.af.xb,imp.af.yb,imp.af.xt,imp.af.yt)
    return(cbind(hap.eu,hap.af))
}
gen.bar.cp.mat<-function(chr,id,shift=0){
    imp.chr.anno<-imp.anno[imp.chr[[chr]],]
    n.s<-nrow(imp.chr.anno)
    imp.width<-diff(imp.chr.anno[,4])
    imp.eu.xb<-c(imp.chr.anno[1,4],imp.chr.anno[2:n.s,4]-imp.width/2)
    imp.eu.xt<-c(imp.chr.anno[1:(n.s-1),4]+imp.width/2,imp.chr.anno[n.s,4])
    imp.eu.yb<-rep(0,length(imp.eu.xb))+shift
    imp.eu.yt<-unphased.anc.copy[[id]][imp.chr[[chr]],1]+shift
    imp.af.xb<-imp.eu.xb
    imp.af.xt<-imp.eu.xt
    imp.af.yb<-imp.eu.yt
    imp.af.yt<-imp.eu.yt+unphased.anc.copy[[id]][imp.chr[[chr]],2]
    imp.uc.xb<-imp.eu.xb
    imp.uc.xt<-imp.eu.xt
    imp.uc.yb<-imp.af.yt
    imp.uc.yt<-rep(1,length(imp.eu.yb))+shift

    hom.eu<-cbind(imp.eu.xb,imp.eu.yb,imp.eu.xt,imp.eu.yt)
    hom.af<-cbind(imp.af.xb,imp.af.yb,imp.af.xt,imp.af.yt)
    hom.uc<-cbind(imp.uc.xb,imp.uc.yb,imp.uc.xt,imp.uc.yt)
    return(cbind(hom.eu,hom.af,hom.uc))
}

id.names<-c('trio1.child','trio1.father','trio1.mother','trio2.child','trio2.father','trio2.mother')

get.chr.ind<-function(chr,id,outdir){
    if(!dir.exists(outdir)){
        dir.create(outdir)
    }
    hap1<-gen.bar.hap.mat(chr,id,1,1.1)
    hap2<-gen.bar.hap.mat(chr,id,2,2.2)
    uc<-gen.bar.cp.mat(chr,id)
    png(paste0(outdir,'/',id.names[id],'_chr',chr,'.png'),width=4.875,height=3.25,units='in',res=1200,pointsize=4)
    plot(0,0,type='n',ann=F,xaxt='n',yaxt='n',main=paste0(id.names[id],'-chr',chr),xlim=range(uc[,1])+10000000,ylim=c(0,3.5),bty='n')
    rect(uc[,1],uc[,2],uc[,3],uc[,4],col='red',border=NA)
    rect(uc[,5],uc[,6],uc[,7],uc[,8],col='yellow',border=NA)
    rect(uc[,9],uc[,10],uc[,11],uc[,12],col='orange',border=NA)
    rect(hap1[,1],hap1[,2],hap1[,3],hap1[,4],col='red',border=NA)
    rect(hap1[,5],hap1[,6],hap1[,7],hap1[,8],col='yellow',border=NA)
    rect(hap2[,1],hap2[,2],hap2[,3],hap2[,4],col='red',border=NA)
    rect(hap2[,5],hap2[,6],hap2[,7],hap2[,8],col='yellow',border=NA)
    text(max(uc[,1])+100000,0.5,'unphased',pos=4)
    text(max(uc[,1])+100000,1.6,'haplotype1',pos=4)
    text(max(uc[,1])+100000,2.7,'haplotype2',pos=4)
    axis(2,at=c(seq(0,1,0.2),seq(0,1,0.2)+1.1,seq(0,1,0.2)+2.2),labels=rep(seq(0,1,0.2),3),cex=0.1)
    legend('top',bty='n',ncol=3,fill=c('red','yellow','orange'),legend=c('European','African','Uncertain'))
    dev.off()
}

get.chr.ind.special<-function(chr,child.id,outdir){
    if(!dir.exists(outdir)){
        dir.create(outdir)
    }
   
    child.uc<-gen.bar.cp.mat(chr,child.id)
    child.hap1<-gen.bar.hap.mat(chr,child.id,1,1.1)
    child.hap2<-gen.bar.hap.mat(chr,child.id,2,2.2)
    father.uc<-gen.bar.cp.mat(chr,child.id+1,3.3)
    father.hap1<-gen.bar.hap.mat(chr,child.id+1,1,4.4)
    father.hap2<-gen.bar.hap.mat(chr,child.id+1,2,5.5)
    mother.uc<-gen.bar.cp.mat(chr,child.id+2,6.6)
    mother.hap1<-gen.bar.hap.mat(chr,child.id+2,1,7.7)
    mother.hap2<-gen.bar.hap.mat(chr,child.id+2,2,8.8)
    
    frame.cord<-function(idx){
        f.b.x0<-min(child.uc[,1])
        f.b.y0<-(idx-1)*1.1-0.03
        f.b.x1<-max(child.uc[,1])+50000
        f.b.y1<-f.b.y0
        
        f.r.x0<-f.b.x1
        f.r.y0<-f.b.y0
        f.r.x1<-f.b.x1
        f.r.y1<-f.b.y1+1+0.06
        
        f.t.x0<-f.b.x0
        f.t.y0<-f.r.y1
        f.t.x1<-f.b.x1
        f.t.y1<-f.r.y1
        fd<-rbind(c(f.b.x0,f.b.y0,f.b.x1,f.b.y1),c(f.r.x0,f.r.y0,f.r.x1,f.r.y1),c(f.t.x0,f.t.y0,f.t.x1,f.t.y1))
    }

    fc<-lapply(c(2,3,6,9),frame.cord)

    png(paste0(outdir,'/',id.names[child.id],'_chr',chr,'.png'),width=4.875,height=3.25,units='in',res=1200,pointsize=4)
    plot(0,0,type='n',ann=F,xaxt='n',yaxt='n',main=paste0(id.names[child.id],'-chr',chr),xlim=range(child.uc[,1])+25000000,ylim=c(0,10.2),bty='n')
    rect(child.uc[,1],child.uc[,2],child.uc[,3],child.uc[,4],col='red',border=NA)
    rect(child.uc[,5],child.uc[,6],child.uc[,7],child.uc[,8],col='yellow',border=NA)
    rect(child.uc[,9],child.uc[,10],child.uc[,11],child.uc[,12],col='orange',border=NA)
    rect(child.hap1[,1],child.hap1[,2],child.hap1[,3],child.hap1[,4],col='red',border=NA)
    rect(child.hap1[,5],child.hap1[,6],child.hap1[,7],child.hap1[,8],col='yellow',border=NA)
    rect(child.hap2[,1],child.hap2[,2],child.hap2[,3],child.hap2[,4],col='red',border=NA)
    rect(child.hap2[,5],child.hap2[,6],child.hap2[,7],child.hap2[,8],col='yellow',border=NA)
    
    rect(father.uc[,1],father.uc[,2],father.uc[,3],father.uc[,4],col='red',border=NA)
    rect(father.uc[,5],father.uc[,6],father.uc[,7],father.uc[,8],col='yellow',border=NA)
    rect(father.uc[,9],father.uc[,10],father.uc[,11],father.uc[,12],col='orange',border=NA)
    rect(father.hap1[,1],father.hap1[,2],father.hap1[,3],father.hap1[,4],col='red',border=NA)
    rect(father.hap1[,5],father.hap1[,6],father.hap1[,7],father.hap1[,8],col='yellow',border=NA)
    rect(father.hap2[,1],father.hap2[,2],father.hap2[,3],father.hap2[,4],col='red',border=NA)
    rect(father.hap2[,5],father.hap2[,6],father.hap2[,7],father.hap2[,8],col='yellow',border=NA)
   
    rect(mother.uc[,1],mother.uc[,2],mother.uc[,3],mother.uc[,4],col='red',border=NA)
    rect(mother.uc[,5],mother.uc[,6],mother.uc[,7],mother.uc[,8],col='yellow',border=NA)
    rect(mother.uc[,9],mother.uc[,10],mother.uc[,11],mother.uc[,12],col='orange',border=NA)
    rect(mother.hap1[,1],mother.hap1[,2],mother.hap1[,3],mother.hap1[,4],col='red',border=NA)
    rect(mother.hap1[,5],mother.hap1[,6],mother.hap1[,7],mother.hap1[,8],col='yellow',border=NA)
    rect(mother.hap2[,1],mother.hap2[,2],mother.hap2[,3],mother.hap2[,4],col='red',border=NA)
    rect(mother.hap2[,5],mother.hap2[,6],mother.hap2[,7],mother.hap2[,8],col='yellow',border=NA)
   
    segments(fc[[1]][,1],fc[[1]][,2],fc[[1]][,3],fc[[1]][,4],lwd=0.3,col='blue')
    segments(fc[[2]][,1],fc[[2]][,2],fc[[2]][,3],fc[[2]][,4],lwd=0.3,col='green')
    segments(fc[[3]][,1],fc[[3]][,2],fc[[3]][,3],fc[[3]][,4],lwd=0.3,col='blue')
    segments(fc[[4]][,1],fc[[4]][,2],fc[[4]][,3],fc[[4]][,4],lwd=0.3,col='green')

    text(max(child.uc[,1])+100000,0.5,'child:unphased',pos=4)
    text(max(child.uc[,1])+100000,1.6,'child:paternal',pos=4,col='blue')
    text(max(child.uc[,1])+100000,2.7,'child:maternal',pos=4,col='green')
    text(max(child.uc[,1])+100000,3.8,'father:unphased',pos=4)
    text(max(child.uc[,1])+100000,4.9,'father:untransmitted',pos=4)
    text(max(child.uc[,1])+100000,6.0,'father:transmitted',pos=4,col='blue')
    text(max(child.uc[,1])+100000,7.1,'mother:unphased',pos=4)
    text(max(child.uc[,1])+100000,8.2,'mother:untransmitted',pos=4)
    text(max(child.uc[,1])+100000,9.3,'mother:transmitted',pos=4,col='green')
    #axis(2,at=c(seq(0,1,0.2),seq(0,1,0.2)+1.1,seq(0,1,0.2)+2.2,seq(0,1,0.2)+3.3,seq(0,1,0.2)+4.4,seq(0,1,0.2)+5.5,seq(0,1,0.2)+6.6),labels=rep(seq(0,1,0.2),7),cex=0.1)
    axis(2,at=rep(seq(0,1,0.5),9)+rep(seq(0,8.8,1.1),each=3),labels=rep(seq(0,1,0.5),9),cex=0.1,las=2)
    legend('top',bty='n',ncol=3,fill=c('red','yellow','orange'),legend=c('European','African','Heterozgous ancestry'))
    legend('topright',legend=c('paternal','maternal'),col=c('blue','green'),lty=1,lwd=0.3)
    dev.off()
}


run.all.ids<-function(chr,out.dir='ancestry_aware_phasing'){
    lapply(c(1,4),get.chr.ind,chr=chr,outdir=out.dir)
}

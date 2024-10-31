library(parallel)

if(!exists('grp29')){
    grp29<-readRDS('data/Groups7_pgs_effect_size_with_CI.rds')
    grp29<-lapply(grp29,function(x)x[[2]])
}

if(!exists('euaf29')){
    euaf29<-readRDS('data/Admixed_sp_pgs_effect_size_with_CI.rds')
    euaf29<-lapply(euaf29,function(x)lapply(x,function(y)y[,1:2]))
}


plot.confintv<-function(plot.dir='Barplots_for_traits'){
	if(!dir.exists(plot.dir)){
		dir.create(plot.dir)
	}
    og<-c('G_SC_England','G_Northumberland','G_Ireland','E_polish','W_India_2','W_CHD_2')
    fk<-colnames(grp29[[1]])
    me<-match(og,fk)
    pheno.titles<-pheno.ids<-names(grp29)
    bp<-function(pheno){
        t8<-grp29[[pheno]]
        tv<-t8[1,me]
		central.intv<-euaf29[[pheno]]
        qz.tv<-central.intv[[1]][5,]
        tv<-c(tv,central.intv[[1]][1,1],central.intv[[1]][2,],central.intv[[1]][4,1],qz.tv[1],central.intv[[1]][4,2],qz.tv[2])
        tv.5<-t8[2,me]
        qz.tv.5<-central.intv[[2]][5,]
        tv.5<-c(tv.5,central.intv[[2]][1,1],central.intv[[2]][2,],central.intv[[2]][4,1],qz.tv.5[1],central.intv[[2]][4,2],qz.tv.5[2])
        qz.tv.95<-central.intv[[3]][5,]
        tv.95<-t8[3,me]
        tv.95<-c(tv.95,central.intv[[3]][1,1],central.intv[[3]][2,],central.intv[[3]][4,1],qz.tv.95[1],central.intv[[3]][4,2],qz.tv.95[2])
        
        names(tv.95)<-names(tv.5)<-names(tv)<-c("South/central England","Northumberland","Ireland","Poland","India","China","(average)","European segments (average)","African segments (average)","European segments if genome 100% African", "European segments if genome 100% European","African segments if genome ~100% African", "African segments if genome 100% European")
        tv<-tv[c("South/central England","Northumberland","Ireland","Poland","India","China","(average)","African segments (average)","European segments (average)","African segments if genome ~100% African","European segments if genome 100% African","African segments if genome 100% European","European segments if genome 100% European")]
        tv.5<-tv.5[c("South/central England","Northumberland","Ireland","Poland","India","China","(average)","African segments (average)","European segments (average)","African segments if genome ~100% African","European segments if genome 100% African","African segments if genome 100% European","European segments if genome 100% European")]
        tv.95<-tv.95[c("South/central England","Northumberland","Ireland","Poland","India","China","(average)","African segments (average)","European segments (average)","African segments if genome ~100% African","European segments if genome 100% African","African segments if genome 100% European","European segments if genome 100% European")]
		cols<-c('blue','blue','blue','blue','lightblue','green','red','red','blue','red','blue','red','blue')
		pdf(paste0(plot.dir,'/',pheno,'_coeff_with_confintv.pdf'))
        ba<-barplot(tv,col=cols,main=NA,ylab='Scaling of polygenic score prediction',las=2,ylim=range(c(0,0.55)),xaxt='n',border=NA)
        arrows(ba,tv.5,ba,tv.95,lwd=0.6,angle=90,code=3,length=0.02)
        #s.segments.v<-segments(x0=ba,y0=tv.5,x1=ba,y1=tv.95,lwd=0.6)
        #s.segments.vb<-segments(x0=ba-0.08,y0=tv.5,x1=ba+0.08,y1=tv.5,lwd=0.6)
        #s.segments.vt<-segments(x0=ba-0.08,y0=tv.95,x1=ba+0.08,y1=tv.95,lwd=0.6)
        text(ba,par('usr')[3],labels=names(tv),srt=35,adj=c(1.1,1.1),xpd=T,cex=0.45)
        dev.off()
    }
    mclapply(pheno.ids,bp,mc.cores=7)
}

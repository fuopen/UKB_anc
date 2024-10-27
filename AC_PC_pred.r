library(parallel)
library(colorRamps)

if(!exists('biobank.ac')){
    biobank.pc<-readRDS('v2_140PCs_487314.rds')
    biobank.ac<-readRDS('v2_127ACs_487314.rds')
    #biobank.pc[,1:2]
    #biobank.pc<-biobank.pc[,-(1:2)]
    colnames(biobank.pc)<-paste0('PC',1:140)
}

ac.pred.pc<-function(){
    app<<-mclapply(biobank.pc,function(x)lm(x ~ .,data=biobank.ac),mc.cores=12)
    app.pred<<-mclapply(app,function(x)predict(x,newdata=biobank.ac),mc.cores=12)
    app.cor<-mapply(cor,biobank.pc,app.pred)
    return(app.cor)
}

pc.pred.ac<-function(){
    ppa<<-mclapply(biobank.ac,function(x)lm(x ~ .,data=biobank.pc),mc.cores=12)
    ppa.pred<<-mclapply(ppa,function(x)predict(x,newdata=biobank.pc),mc.cores=12)
    ppa.cor<-mapply(cor,biobank.ac,ppa.pred)
    return(ppa.cor)
}

if(!exists('ppa')){
	ppa<-readRDS('predicted_ACs_by_140PCs.rds')
	ppa.pred<-ppa$ppa.pred
	ppa.cor<-ppa$ppa.cor
}

if(!exists('app')){
	app<-readRDS('predicted_PCs_by_127ACs.rds')
	app.pred<-app$app.pred
	app.cor<-app$app.cor
}

plot.ac.pred<-function(){
    png('AC_pred_17_140PCs2.png',width=18.5,height=18.5,units='in',res=1200,pointsize=4)
    par(mar=c(5,3,5,3)+0.2,oma=c(2,2,2,2),cex.main=5,cex.axis=3)
    cols<-rainbow(140)
    layout(matrix(1:144,ncol=12,byrow=T))
    for(i in 1:140){
		print(i)
        plot(biobank.pc[[i]],app.pred[[i]],main=paste0(names(biobank.pc)[i],':',round(app.cor[i]^2,3)),xlab=NA,ylab=NA,col=cols[i])
        abline(a=0,b=1,col='red',lwd=2.3)
    }
    dev.off()
}
plot.pc.pred<-function(){
    png('PCs_pred_127AC2.png',width=18.5,height=18.5,units='in',res=1200,pointsize=4)
    par(mar=c(5,3,5,3)+0.2,oma=c(2,2,2,2),cex.main=2.5,cex.axis=2)
    cols<-rainbow(127)
    layout(matrix(1:144,ncol=12,byrow=T))
	od<-order(ppa.cor,decreasing=T)
	biobank.ac<-biobank.ac[od]
	ppa.pred<-ppa.pred[od]
	ppa.cor<-ppa.cor[od]
    for(i in 1:127){
		print(i)
        plot(biobank.ac[[i]],ppa.pred[[i]],main=paste0(names(biobank.ac)[i],':',round(ppa.cor[i]^2,3)),xlab=NA,ylab=NA,col=cols[i])
        abline(a=0,b=1,col='red',lwd=2.3)
    }
    dev.off()
}

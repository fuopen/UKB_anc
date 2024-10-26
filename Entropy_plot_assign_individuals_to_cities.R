library(parallel)
library(sp)
library(fastcluster)
if(!exists('gb.spdf')){
	gb.spdf<-readRDS('v2_self_BI_487409.rds')
	bi.spdf<-spTransform(gb.spdf,CRS("+init=epsg:4326"))
}

if(!exists('England_Wales')){
	England_Wales<-readRDS('UK_main_cities_England_Wales.rds')
	ew<-spTransform(England_Wales,CRS("+init=epsg:4326"))
}

if(!exists('Glasgow')){
	Glasgow<-readRDS('Scotland/Glasgow.rds')
	sc.gl<-spTransform(Glasgow,CRS("+init=epsg:4326"))
}

if(!exists('Dundee')){
	Dundee<-readRDS('Scotland/Dundee.rds')
	sc.dd<-spTransform(Dundee,CRS("+init=epsg:4326"))
}

if(!exists('Edinburgh')){
	Edinburgh<-readRDS('Scotland/Edinburgh.rds')
	sc.ed<-spTransform(Edinburgh,CRS("+init=epsg:4326"))
}

if(!exists('Aberdeen')){
	Aberdeen<-readRDS('Scotland/Aberdeen.rds')
	sc.sb<-spTransform(Aberdeen,CRS("+init=epsg:4326"))
}

if(!exists('gb.counties')){
	gb.counties<-readRDS('GBR_adm2.rds')
	gb.counties<-spTransform(gb.counties,CRS('+init=epsg:4326'))
}

sampleInCity<-function(){
	ew.cities<-bi.spdf%over%ew
	gl<-bi.spdf%over%sc.gl
	dd<-bi.spdf%over%sc.dd
	ed<-bi.spdf%over%sc.ed
	sb<-bi.spdf%over%sc.sb
	city.in<-list(EnglandWales=ew.cities,Glasgow=gl,Dundee=dd,Edinburgh=ed,Aberdeen=sb)
}

sampleInCounty<-function(){
	county.in<-bi.spdf%over%gb.counties
}

get.city.ac<-function(){
	ac.city<-readRDS('EW_SC_4_cities.rds')
	ac.county<-readRDS('Counties.rds')
	city.names.plot<-list(
		DundeeOther=c('Fife','Angus','Perthshire and Kinross'),
		Aberdeen='',
		Edinburgh='',
		#Northumberland='',
		Hartlepool=c('Hartlepool','Stockon-on-Tees','Middlesbrough'),
		Yorkshire=c('East Riding of Yorkshire','North East Lincolnshire'),
		Birmingham='Birmingham',
		Norfolk='Norfolk',
		Suffolk='Suffolk',
		London='London',
		Oxfordshire='Oxfordshire',
		IW='Isle of Wight',
		Bristol='Bristol',
		Plymouth='Plymouth',
		Cardiff='Cardiff',
		Ceredigion=c('Ceredigion','Carmarthenshire'),
		Anglese=c('Anglese','Gwynedd'),
		Manchester='Manchester',
		Liverpool='Liverpool',
		WestYorkshire='West Yorkshire',
		#Cumbria='Cumbria',
		Glasgow=''
	)
	city.ids<-list(
		Bristol=which(!is.na(si.city$EnglandWales$tcity15nm) & as.character(si.city$EnglandWales$tcity15nm)=='Bristol'),
		Birmingham=which(!is.na(si.city$EnglandWales$tcity15nm) & as.character(si.city$EnglandWales$tcity15nm)=='Birmingham'),
		London=which(!is.na(si.city$EnglandWales$tcity15nm) & as.character(si.city$EnglandWales$tcity15nm)=='London'),
		Liverpool=which(!is.na(si.city$EnglandWales$tcity15nm) & as.character(si.city$EnglandWales$tcity15nm)=='Liverpool'),
		Aberdeen=which(!is.na(si.city$Aberdeen$Name)),
		Edinburgh=which(!is.na(si.city$Edinburgh$Name)),
		Glasgow=which(!is.na(si.city$Glasgow$Name)),
		Dundee=which(!is.na(si.city$Dundee$NAME))
	)	
	city.other<-city.names.plot[!names(city.names.plot)%in%names(city.ids)]	
	city.oids<-lapply(city.other,function(x)which(!is.na(ac.county$NAME_2) & ac.county$NAME_2%in%x))
	city.ids$Dundee<-unique(c(city.ids$Dundee,city.oids$DundeeOther))
	city.oids$DundeeOther<-NULL
	city.aids<-c(city.ids,city.oids)
	#city.aids<-list(c1=city.ids,c2=city.oids)
}

color=grDevices::colors()[grep('gr(a|e)y',grDevices::colors(),invert=T)]
if(!exists('new.colors')){
	gb.color<-readRDS('23_regions_color1.rds')
	g23.n<-sapply(strsplit("G_Anglia,G_SC_England,G_C_England,G_NC_England,G_N_Yorkshire,G_Northumberland,G_SE_England,G_Merseyside,G_S_England,G_SE_Wales,G_NE_Scot,G_SW_Wales,G_Devon,G_NW_Wales,G_S_Yorkshire,G_Cornwall_Tip,G_Orkney,G_Lincolnshire,G_Cumbria,G_NW_Scot,G_Noi,G_Ireland,G_Cornwall",','),function(x)x)
	names(gb.color)<-g23.n
	set.seed(20181223)
	other.colors<-color[!color%in% gb.color]
	other.colors<-sample(other.colors,ncol(bi.spdf)-length(gb.color),replace=T)
	names(other.colors)<-names(bi.spdf)[!names(bi.spdf)%in%names(gb.color)]
	new.colors<-c(gb.color,other.colors)
	new.colors<-new.colors[names(bi.spdf)]
	new.colors[['E_norwegian']]<-'chartreuse3'
	new.colors[['E_Switzerland_1']]<-'chartreuse4'
	new.colors[['E_Portugal']]<-'darkolivegreen2'
	new.colors[['E_TSI']]<-'gold'
	new.colors[['W_syrian']]<-'lightpink3'
	new.colors[['E_Denmark']]<-'springgreen'
	official.names<-c('Anglia','South Central England','Central England','North Central England','North Yorkshire','Northumberland','South East England','Merseyside','South England','South East Wales','North East Scotland','South West Wales','Devon','North West Wales','South Yorkshire','Cornwall Tip','Orkney','Lincolnshire','Cumbria','North West Scotland','North Ireland and South Scotland','Ireland','Cornwall')
	names(official.names)<-names(gb.color)
	official.others<-c('France','French Basque','Sardinian','North Italian','Italy1','Adygei','Russia','Ukraine','Bulgaria','Croatia','Cyprus','DutchlikeGerman','Germany', 'Greece', 'Lithuania', 'Mordovia', 'Norway', 'Poland', 'Sicily','Spain', 'Switzerland','Yugoslavia','Czech Republic','Romania', 'Hungary', 'Austria', 'BelgianlikeGerman', 'Portugal', 'PortlikeSpain', 'Belgium', 'Netherlands', 'Finland', 'Italy2', 'Tuscany', 'Biaka Pygmies', 'Mandenka', 'Mbuti Pygmies', 'LWK', 'Bedouin', 'Palestine', 'Hazara', 'Mozabite', 'Balochi', 'San', 'India', 'SriLanka', 'Pathan', 'Kalash', 'Burusho', 'Sindhi', 'Yizu', 'Miaozu', 'Daur', 'Hezhen', 'Xibo', 'Uygur', 'Dai', 'She', 'Tu', 'Yakut', 'Papuan', 'Melanesian', 'Lahu', 'Naxi', 'Cambodia', 'Pima', 'Colombia', 'Karitiana', 'Nama', 'GIH', 'MKK', 'Chechen', 'Chuvash', 'Syria', 'East Turkey', 'South Turkey', 'AMHARA', 'ANUAK', 'ARI', 'Somali', 'GUMUZ', 'amaXhosa', 'ColouredWellington', 'Karretjie', 'GuiGhanaKgal', 'Khomani', 'Khwe', 'Yoruba', 'Juhoan', 'Maya', 'Japan', 'China', 'Denmark', 'Sweden', 'Algeria', 'Egypt', 'Morocco', 'Libya', 'Nepal', 'Myanmar', 'Bangladesh', 'Thailand', 'Malaysia', 'Philippines')
	names(official.others)<-colnames(gb.spdf)[24:127]	
}

get.city.ancestry<-function(){
	city.ids<-readRDS('UKB_BI_ids_in_20_cities.rds')
	ac.new<-bi.spdf@data
	ac.new$E_italians<-ac.new$E_italians+ac.new$E_Italy_3+ac.new$E_North_Italian
	ac.new$E_Italy_3<-NULL
	ac.new$E_North_Italian<-NULL
	ac.new$E_French<-ac.new$E_French+ac.new$E_French_Basque
	ac.new$E_French_Basque<-NULL
	ac.new$E_german<-ac.new$E_german+ac.new$E_DutchlikeGerman+ac.new$E_BelgianlikeGerman
	ac.new$E_DutchlikeGerman<-NULL
	ac.new$E_BelgianlikeGerman<-NULL
	ac.new$E_Spain_3<-ac.new$E_Spain_3+ac.new$E_PortlikeSpain
	ac.new$E_PortlikeSpain<-NULL
	ac.meananc<-sapply(ac.new,function(x)mean(x>0.01))
	ac.new<-ac.new[,ac.meananc>0.01]
	city.anc<-lapply(city.ids,function(x)ac.new[x,])	
}
if(!exists('gcaa')){
	gcaa<-get.city.ancestry()
}

barplot.eachcity<-function(reg.name,pdir='GB_city_fig'){
	if(!dir.exists(pdir)){
		dir.create(pdir)
	}
	mat<-as.matrix(gcaa[[reg.name]])
	mat.which<-apply(mat,1,which.max)
	mat.reg<-colnames(mat)[mat.which]
	mro<-names(sort(table(mat.which),decreasing=T))
	mat.reg.ord<-colnames(mat)[as.integer(mro)]
	me<-colnames(mat)
	mes<-me[!me%in%mat.reg.ord]
	mat<-mat[,c(mat.reg.ord,mes)]
	mat.reg.list<-by(mat,mat.reg,function(x){
		x<-as.matrix(x)
		if(nrow(x)>2){
			x.id<-hclust.vector(x)
			x<-x[rev(x.id$order),]
		}
		return(x)
	})
	mat<-do.call(rbind,mat.reg.list[mat.reg.ord])
	png(paste0(pdir,'/',reg.name,'_strplot_',nrow(mat),'_order.png'),width=8.25,height=2.25,units='in',res=1200,pointsize=4)
	par(mar=c(0,0,0,0),oma=c(0,0,0,0))
	barplot(t(mat),col=new.colors[colnames(mat)],border=NA,space=0,ann=F,axes=F,xaxt='n',yaxt='n')
	dev.off()
}

make.figs<-function(){
	mu<-mclapply(names(gcaa),barplot.eachcity,mc.cores=4)	
}

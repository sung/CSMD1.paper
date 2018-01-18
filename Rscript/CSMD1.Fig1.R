
	my.RData=file.path("../RData",paste(my.tissue,my.cpg.type,"dt.meth.region.RData",sep="."))
	if(file.exists(my.RData)){
		cat("loading dt.meth.region.RData...\n")
		load(my.RData)
		cat("dt.meth.region loaded\n")
	}else{
		stop("RData not exist")
	}

	########################
	## Figure 1           ##
	## Beta Distribution  ##
	## via plot_grid      ##
	########################
	library("cowplot")
	dt.meth=dt.meth.region[Tissue==my.tissue & Region=="Sites"]
	dt.meth[,chr.type:=ifelse(chr=="X","X-chromosome","Autosomes")]
	p.site.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_line(stat="density", size=1) +
			ggtitle("CpG sites") +
			labs(x="",y="Density") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			facet_grid(.~chr.type) +
			theme_Publication() +
			theme(legend.position="none")

	dt.meth=dt.meth.region[Tissue==my.tissue & Region=="Genes"]
	dt.meth[,chr.type:=ifelse(chr=="X","X-chromosome","Autosomes")]
	p.gene.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_line(stat="density", size=1) +
			ggtitle("Gene-bodies") +
			labs(x="",y="") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			facet_grid(.~chr.type) +
			theme_Publication() +
			theme(legend.position="none")
	
	dt.meth=dt.meth.region[Tissue==my.tissue & Region=="Promo_15.05"]
	dt.meth[,chr.type:=ifelse(chr=="X","X-chromosome","Autosomes")]
	p.promoter.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_line(stat="density", size=1) +
			ggtitle("Promoters") +
			labs(x="% Methylation Level", y="Density") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			facet_grid(.~chr.type) +
			theme_Publication() +
			theme(legend.position="none")

	dt.meth=dt.meth.region[Tissue==my.tissue & Region=="CPGi"]
	dt.meth[,chr.type:=ifelse(chr=="X","X-chromosome","Autosomes")]
	p.cpgi.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_line(stat="density", size=1) +
			ggtitle("CpG islands") +
			labs(x="% Methylation Level", y="") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			facet_grid(.~chr.type) +
			theme_Publication() + 
			theme(legend.position=c(0.91,0.87)) 

	file.name<-file.path("../Figures/Beta.dist",paste("Fig1",my.tissue,my.cpg.type,time.stamp,"line.tiff",sep="."))
	tiff(filename=file.name, width=12, height=10, units="in",res=300, compression = 'lzw') #A4 size 
	cowplot::plot_grid(p.site.beta, p.gene.beta, p.promoter.beta, p.cpgi.beta, labels=c("A","B","C","D"), label_size=20, ncol = 2, nrow = 2)
	dev.off()



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
	p.site.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_line(stat="density", size=1) +
			labs(x="",y="Density") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			theme_Publication() +
			theme(legend.position="none")

	dt.meth=dt.meth.region[Tissue==my.tissue & Region=="Genes"]
	min.cpg.gene<-quantile(dt.meth[,num.sites], .25) # No. of CpG at 25%-percentile # min.cpg: 6
	p.gene.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_line(stat="density", size=1) +
			labs(x="",y="") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			theme_Publication() +
			theme(legend.position=c(0.91,0.91)) 
	
	dt.meth=dt.meth.region[Tissue==my.tissue & Region=="Promo_15.05"]
	min.cpg.promo<-quantile(dt.meth[,num.sites], .25) # No. of CpG at 25%-percentile # min.cpg: 4
	p.promoter.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_line(stat="density", size=1) +
			labs(x="% Methylation Level", y="Density") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			theme_Publication() +
			theme(legend.position="none")

	dt.meth=dt.meth.region[Tissue==my.tissue & Region=="CPGi"]
	min.cpg.cpgi<-quantile(dt.meth[,num.sites], .25) # No. of CpG at 25%-percentile # min.cpg: 28
	#p.cpgi.beta<-ggplot(dt.meth[num.sites>=min.cpg.cpgi], aes(meth*100, col=Gender)) + 
	p.cpgi.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_line(stat="density", size=1) +
			labs(x="% Methylation Level", y="") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			theme_Publication() + 
			theme(legend.position="none")

	file.name<-file.path("../Figures/Beta.dist",paste("Fig1",my.tissue,my.cpg.type,time.stamp,"line.tiff",sep="."))
	tiff(filename=file.name,width=10, height=8,units="in",res=300, compression = 'lzw') #A4 size 
	cowplot::plot_grid(p.site.beta, p.gene.beta, p.promoter.beta, p.cpgi.beta, labels=c("A","B","C","D"), label_size=15, ncol = 2, nrow = 2)
	dev.off()


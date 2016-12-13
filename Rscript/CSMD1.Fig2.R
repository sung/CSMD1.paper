
	my.RData=file.path("../RData",paste(my.tissue,my.cpg.type,"dt.meth.region.RData",sep="."))
	if(file.exists(my.RData)){
		cat("loading dt.meth.region.RData...\n")
		load(my.RData)
		cat("dt.meth.region loaded\n")
	}else{
		stop("RData not exist")
	}

	#######################
	## Figure 2          ##
	## Meth Diff Boxplot ##
	## via plot_grid     ##
	#######################
	library("cowplot")
	dt.meth=dt.meth.region[Tissue==my.tissue & Region=="Sites"]
	dt.meth.dcast.site=dcast(dt.meth, chr+strand+start+end~Gender, value.var="meth") # dcast: row(long) to column(wide)
	p.site.box<-ggplot(dt.meth.dcast.site, aes(chr, (Female-Male)*100)) + 
				geom_boxplot(outlier.shape=NA) + 
				labs(x="",y="% Methylation Difference\n   (more in male)     (more in female)") + 
				geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
				ylim(c(-50,50)) +
				scale_x_discrete(limits=my.chr.order) +
				theme_Publication()

	dt.meth=dt.meth.region[Tissue==my.tissue & Region=="Genes"]
	min.cpg.gene<-quantile(dt.meth[,num.sites], .25) # No. of CpG at 25%-percentile # min.cpg: 6
	dt.meth.dcast.gene=dcast(dt.meth, chr+strand+start+end~Gender, value.var="meth") # dcast: row(long) to column(wide)
	p.gene.box<-ggplot(dt.meth.dcast.gene, aes(chr, (Female-Male)*100)) + 
				geom_boxplot(outlier.shape=NA) + 
				labs(x="", y="") + 
				geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
				ylim(c(-50,50)) +
				scale_x_discrete(limits=my.chr.order) +
				theme_Publication()

	dt.meth=dt.meth.region[Tissue==my.tissue & Region=="Promo_15.05"]
	min.cpg.promo<-quantile(dt.meth[,num.sites], .25) # No. of CpG at 25%-percentile # min.cpg: 4
	dt.meth.dcast.promo=dcast(dt.meth, chr+strand+start+end~Gender, value.var="meth") # dcast: row(long) to column(wide)
	p.promoter.box<-ggplot(dt.meth.dcast.promo, aes(chr, (Female-Male)*100)) + 
					geom_boxplot(outlier.shape=NA) + 
					labs(x="Chromosome",y="% Methylation Difference\n   (more in male)     (more in female)") + 
					geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
					ylim(c(-50,50)) +
					scale_x_discrete(limits=my.chr.order) +
					theme_Publication()

	dt.meth=dt.meth.region[Tissue==my.tissue & Region=="CPGi"]
	min.cpg.cpgi<-quantile(dt.meth[,num.sites], .25) # No. of CpG at 25%-percentile # min.cpg: 28
	dt.meth.dcast.cpgi=dcast(dt.meth, chr+strand+start+end~Gender, value.var="meth") # dcast: row(long) to column(wide)
	p.cpgi.box<-ggplot(dt.meth.dcast.cpgi, aes(chr, (Female-Male)*100)) + 
					geom_boxplot(outlier.shape=NA) + 
					labs(x="Chromosome", y="") + 
					geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
					ylim(c(-50,50)) +
					scale_x_discrete(limits=my.chr.order) +
					theme_Publication()

	file.name<-file.path("../Figures/Meth.diff",paste("Fig2",my.tissue,my.cpg.type,time.stamp,"tiff",sep="."))
	tiff(filename=file.name,width=14, height=10,units="in",res=300, compression = 'lzw') #A4 size 
	cowplot::plot_grid(p.site.box, p.gene.box, p.promoter.box, p.cpgi.box, labels=c("A","B","C","D"), label_size=15, ncol = 2, nrow = 2)
	dev.off()

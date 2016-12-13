
	library(ggbio)
	#########################################################
	# Supplementary Figure
	# PLOT WGBS, SureSelect and RoadMap methylation profile #
	#########################################################
	load("../RData/CSMD1.meth.POPS.RoadMap.RData") # loading "rnb.diff.met"
	cat("wgbs.locus, sureselect.tech.replicate, sureselect.bio.replicate, and dt.tissue.locus loaded\n") 
	# wgbs.locus: WGoxBS for 2 males and 2 females
	# sureselect.tech.replicate: technical replicates
	# sureselect.bio.replicate: biological replicates
	# dt.tissue.locus: RoadMap methylation

	my.target<-"CSMD1.DMR2" # one of following
							# CSMD1.DMR:
							# CSMD1.DMR2:
   							# CSMD1.legacy.promo:
   							# CSMD1.putative.promo:
	my.contrast<-"Gender"
	mycol<-my.col[[my.contrast]] # colour defined from local.R
	my.gr.ensg <- reduce(rtracklayer::import.bed(target.region[[my.target]])) # target.region defined from local.R
	names(my.gr.ensg)<-my.target

	####################
	# transcript track #
	####################
	p.tr.reduced<-autoplot(hg19ensGene, which=my.gr.ensg, stat="reduce", color="brown", fill="brown") + theme_alignment() # hg19ensGene defined in 'Annotation.R'

	##########################
	# WGBS methylation GGBIO #
	##########################
	my.gr.wgbs<-as.data.frame(wgbs.locus[[my.target]])
	my.gr.wgbs$type='WGoxBS (N=4)' # isa 'data.frame'

	my.point.size=.8; my.alpha=0.4
	p.wgbs2<-ggplot(my.gr.wgbs, aes(start, score, col=Gender)) + 
				geom_point(size=my.point.size,alpha=my.alpha) +
				ylim(0,100) +
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
				labs(y="") +
				theme_Publication() + 
				theme(legend.position="top")
	####################
	## RoadMap        ##
	## Schultz et al. ##
	####################
	p.roadmap<-ggplot(rbindlist(dt.tissue.locus[[my.target]]), aes(start, methylation_level, col=Gender)) + geom_point(size=my.point.size,alpha=my.alpha) + ylim(0, 100) + 
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
				labs(y="% Methylation") +
				theme_Publication() + 
				theme(legend.position="none") + 
				facet_grid(Tissue ~ .)

	# wgbs, tech, bio, roadmap
	my.file.name<-"../Figures/Suppl/SupplFig.CSMD1.DMR+2exons.methylation.profiles.RoadMap.tiff"
	tiff(filename=my.file.name,width=10, height=15 ,units="in",res=300, compression = 'lzw') #A4 size 
	print(tracks(` `=p.tr.reduced, `Placenta`=p.wgbs2, `Schultz et al`=p.roadmap, heights = c(1,4,28), scale.height=2.5, label.text.cex = 1.5, xlim=my.gr.ensg))
	dev.off()

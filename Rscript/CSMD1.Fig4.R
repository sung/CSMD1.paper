
	library(ggbio)
	################################
	## Figure 4.                  ##
	## CSMD1 DMR + 2 exons        ## 
	## Methylation Profiles       ##
	################################
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

	################ 
	## Validation ##
	################ 
	# PT technical validation
	my.gr.tech<-as.data.frame(sureselect.tech.replicate[[my.target]])
	my.gr.tech$type="Technical (n=2)"
	# PT biological validation
	my.gr.bio<-as.data.frame(sureselect.bio.replicate[[my.target]])
	my.gr.bio$type="Biological (n=6)"
	# merger of two above
	df.validation<-rbind(my.gr.tech, my.gr.bio); 
	df.validation$type<-factor(df.validation$type,levels = c("Technical (n=2)", "Biological (n=6)")) # re-order this facet

	##############################################
	# PLOT WGBS & SureSelect methylation profile #
	##############################################
	my.point.size=1.7; my.alpha=0.5
	p.wgbs<-ggplot(my.gr.wgbs, aes(start, score, col=Gender)) + 
				geom_point(size=my.point.size,alpha=my.alpha) +
				ylim(0,100) +
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="", values=my.col[["Gender"]]) + 
				labs(y="% Methylation") +
				theme_Publication() + 
				theme(legend.position="top")
	p.validation<-ggplot(df.validation, aes(start, score, col=Gender)) + 
				geom_point(size=my.point.size,alpha=my.alpha) +
				ylim(0,100) +
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="", values=my.col[["Gender"]]) + 
				labs(y="% Methylation") +
				theme_Publication() + 
				#theme(legend.position=c(0.85,0.6)) +
				theme(legend.position="none") +
				facet_grid(type~ .)
	p.pt.meth<-tracks(`CSMD1`=p.tr.reduced, `WGoxBS (n=4)`=p.wgbs, `In-solution target capture`=p.validation, heights = c(1.2,4,7), scale.height=2.5, label.text.cex = 1.4, xlim=my.gr.ensg)
	my.file.name<-"../Figures/CSMD1.DMR.profile/Fig4.CSMD1.DMR+2exons.methylation.profiles.POPS.ver1.tiff"
	tiff(filename=my.file.name,width=15, height=8.5 ,units="in",res=300, compression = 'lzw') #A4 size 
	print(p.pt.meth)
	dev.off()


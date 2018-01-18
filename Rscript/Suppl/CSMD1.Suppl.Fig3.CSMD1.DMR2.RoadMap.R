
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
	p.tr.reduced<-autoplot(hg19ensGene, which=my.gr.ensg, stat="reduce", color="black", fill="black") + theme_alignment() # hg19ensGene defined in 'Annotation.R'

	##########################
	# WGBS methylation GGBIO #
	##########################
	my.gr.wgbs<-as.data.frame(wgbs.locus[[my.target]])
	my.gr.wgbs$Tissue='PT'

	my.point.size=.8; my.alpha=0.4
	p.wgbs2<-ggplot(my.gr.wgbs, aes(start, score, col=Gender)) + 
				geom_point(size=my.point.size,alpha=my.alpha) +
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
				scale_y_continuous(breaks=seq(from=0,to=100,by=50)) +
				labs(y="") +
				facet_grid(Tissue~ .) +
				theme_Publication() + 
				theme(legend.position="none")
				#theme(legend.position="top")
	####################
	## RoadMap        ##
	## Schultz et al. ##
	####################
	p.roadmap<-ggplot(rbindlist(dt.tissue.locus[[my.target]]), aes(start, methylation_level, col=Gender)) + geom_point(size=my.point.size,alpha=my.alpha) +
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
				scale_y_continuous(breaks=seq(from=0,to=100,by=50)) +
				labs(y="% Methylation") +
				facet_grid(Tissue ~ .) +
				theme_Publication() + 
				theme(legend.position="none") 

	#####################
	## RoadMap         ##
	## Broad Institute ##
	#####################
	dt.wgbs.BI<-rbind(
		`Hippocampus`=rbind(
						`Hippocampus.149`=as.data.table(as.data.frame(import.bed("../RData/GSM1112838_BI.Brain_Hippocampus_Middle.Bisulfite-Seq.149.GRCh37.CSMD1.bed.gz"))),
						`Hippocampus.150`=as.data.table(as.data.frame(import.bed("../RData/GSM916050_BI.Brain_Hippocampus_Middle.Bisulfite-Seq.150.GRCh37.CSMD1.bed.gz")))
						)[,.(`score`=round(mean(score)*100,2),`Gender`="Male",`Tissue`="HC"),"seqnames,end,strand"], # male
		`Fetal_Thymus`=as.data.table(as.data.frame(import.bed("../RData/GSM1172595_BI.Fetal_Thymus.Bisulfite-Seq.UW_H24943.GRCh37.CSMD1.bed.gz")))[,.(seqnames,end,strand,score=score*100,`Gender`="Female",`Tissue`="fTH")], # female
		`Fetal_Muscle`=as.data.table(as.data.frame(import.bed("../RData/GSM1172596_BI.Fetal_Muscle_Leg.Bisulfite-Seq.UW_H24996.GRCh37.CSMD1.bed.gz")))[,.(seqnames,end,strand,score=score*100,`Gender`="Female",`Tissue`="fMU")] # female
		#`LI`=as.data.table(as.data.frame(import.bed("../RData/GSM916049_BI.Adult_Liver.Bisulfite-Seq.3.GRCh37.CSMD1.bed.gz")))[,.(seqnames,end,strand,score=score*100,`Gender`="NA",`Tissue`="LI")] # unknown
	)
	p.meth.BI<-ggplot(dt.wgbs.BI[end>=start(ranges(my.gr.ensg)) & end<=end(ranges(my.gr.ensg))], aes(end, score, col=Gender)) + 
				geom_point(size=my.point.size,alpha=my.alpha) +
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
				scale_y_continuous(breaks=seq(from=0,to=100,by=50)) +
				labs(y="% Methylation") +
				facet_grid(Tissue~ .) +
				theme_Publication() + 
				theme(legend.position="none")

	##############
	## RoadMap  ##
	## UCSF-UBC ##
	##############
	dt.wgbs.UCSF<-rbind(
		`Brain_Germinal_Matrix`=as.data.table(as.data.frame(import.bed("../RData/GSM941747_UCSF-UBC.Brain_Germinal_Matrix.Bisulfite-Seq.HuFGM02.GRCh37.CSMD1.bed.gz")))[,.(seqnames,end,strand,score=score*100,`Gender`="Male",`Tissue`="GM")], # male
		`Breast_Luminal_Epithelial`=as.data.table(as.data.frame(import.bed("../RData/GSM1127125_UCSF-UBC.Breast_Luminal_Epithelial_Cells.Bisulfite-Seq.RM066.GRCh37.CSMD1.bed.gz")))[,.(seqnames,end,strand,score=score*100,`Gender`="Female",`Tissue`="LE")], # female
		`Testis_Spermatozoa`=as.data.table(as.data.frame(import.bed("../RData/GSM1127117_UCSF-UBC.Testis_Spermatozoa_Primary_Cells.Bisulfite-Seq.390ATA.GRCh37.CSMD1.bed.gz")))[,.(seqnames,end,strand,score=score*100,`Gender`="Male",`Tissue`="SP")] # male
	)
	p.meth.UCSF<-ggplot(dt.wgbs.UCSF[end>=start(ranges(my.gr.ensg)) & end<=end(ranges(my.gr.ensg))], aes(end, score, col=Gender)) + 
				geom_point(size=my.point.size,alpha=my.alpha) +
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
				scale_y_continuous(breaks=seq(from=0,to=100,by=50)) +
				labs(y="% Methylation") +
				facet_grid(Tissue~ .) +
				theme_Publication() + 
				theme(legend.position="none")

	#################################
	# Placenta(WGoxBS),            ##
	# Schultz et al (n=8 tissues), ## 
	# BI (n=4 tissues),            ##
	# UCSF-UBC (n=3 tissues)       ##
	my.file.name<-"../Figures/Suppl/Suppl.Fig.CSMD1.DMR+2exons.methylation.profiles.RoadMap.UCSD.BI.UCSF-UBC.tiff"
	tiff(filename=my.file.name,width=10, height=15 ,units="in",res=300, compression = 'lzw') #A4 size 
	print(tracks(` `=p.tr.reduced, `POPS`=p.wgbs2, `Schultz et al. (UCSD)`=p.roadmap, `Broad Institute`=p.meth.BI, `UCSF-UBC`=p.meth.UCSF, heights = c(1,3,24,9,9), scale.height=2.5, label.text.cex = 1.5, xlim=my.gr.ensg))
	dev.off()

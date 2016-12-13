
	my.exon.RData<-paste0("../RData/CSMD1.",TR_PREFIX,".POPS.RoadMap.Fetal.exon.RData")
	load(my.exon.RData)
	cat("POPS.ddsFpkm, RoadMap.ddsFpkm, RoadMap.colData, Fetal.ddsFpkm, and Fetal.colData loaded\n")
	# POPS.ddsFpkm: Fpkm of Boy.Girl.exon.GRCh37 AGA samples
	# RoadMap.ddsFpkms: Fpkm of RoadMap Consortium (placenta tissues) and Schultz data
	# RoadMap.colData: meta-info of RoadMap placenta & Schultz data 
	# Fetal.ddsFpkm: Fpkm of RoadMap fetal and Manchester embryo data 
	# Fetal.colData: meta-info of fetal and embryo data

	######################
	## Fetal and Embyro ##
	######################
	library(pheatmap)
	merged.fpkm<-cbind(Fetal.ddsFpkm, `PT`=rowMeans(POPS.ddsFpkm)) # isa 'matrix'

	ann_tissues<-data.frame(
					`Tissue`=factor(c(as.character(Fetal.colData$Origin), "Placenta-derived")),
					`Source`=factor(c(ifelse(Fetal.colData$Source=="Schultz", "Schultz et al.", ifelse(Fetal.colData$Source=="Consortium","RoadMap Consortium","Gerrad et al.")), "This Study")),
					row.names=colnames(merged.fpkm)
				)

	ann_colors<-list(`Tissue`=c(`Placenta-derived`='darkmagenta', `Embryo-Somatic`='yellow4', `Fetal-Somatic`='yellowgreen'),
					 `Source`=c(`RoadMap Consortium`=cbPalette2[2], `This Study`=cbPalette2[3], `Gerrad et al.`=cbPalette2[4]))

	# portrait (row; tissues, col: exons)
	my.filename<-"../Fgures/Suppl/Fig.pheatmap.exon.pops.roadmap.fetal.GRCh37.portrait.tiff"
	pheatmap(log(merged.fpkm[rev(rownames(merged.fpkm)),]/10^3+1), cluster_rows=F, cluster_cols=T, annotation_col=ann_tissues, annotation_colors=ann_colors, width=10, height=15, fontsize_col=10, fontsize_row=7, cellheight=9, filename=my.filename)

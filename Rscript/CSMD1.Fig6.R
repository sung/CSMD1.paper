	# part of Rscript/CSMD1.paper.figures.main.R

	my.exon.RData<-paste0("../RData/CSMD1.",TR_PREFIX,".POPS.RoadMap.Fetal.exon.RData")
	load(my.exon.RData)
	cat("POPS.ddsFpkm, RoadMap.ddsFpkm, RoadMap.colData, Fetal.ddsFpkm, and Fetal.colData loaded\n")
	# POPS.ddsFpkm: Fpkm of Boy.Girl.exon.GRCh37 AGA samples
	# RoadMap.ddsFpkms: Fpkm of RoadMap Consortium (placenta tissues) and Schultz data
	# RoadMap.colData: meta-info of RoadMap placenta & Schultz data 
	# Fetal.ddsFpkm: Fpkm of RoadMap fetal and Manchester embryo data 
	# Fetal.colData: meta-info of fetal and embryo data

	#########################
	# tissue-level pheatmap #
	#########################
	library(pheatmap)
	merged.fpkm<-cbind(RoadMap.ddsFpkm, `PT`=rowMeans(POPS.ddsFpkm)) # isa 'matrix'

	ann_tissues<-data.frame(
					`Tissue`=factor(c(ifelse(RoadMap.colData$Source=="Schultz", "Non-Placenta", "Placenta-derived"), "Placenta-derived")),
					`Source`=factor(c(ifelse(RoadMap.colData$Source=="Schultz", "Schultz et al.", ifelse(RoadMap.colData$Source=="Consortium","RoadMap Consortium","Gerrad et al.")), "This Study")),
					row.names=colnames(merged.fpkm)
				)
	ann_colors<-list(
					 `Tissue`=c(`Non-Placenta`='darkgreen',`Placenta-derived`='darkmagenta'),
					 `Source`=c(`Schultz et al.`=cbPalette2[1], `RoadMap Consortium`=cbPalette2[2], `This Study`=cbPalette2[3])
				)

	# portrait (row; tissues, col: exons)
	my.filename<-"../Figures/Pheatmap.exon/Fig.pheatmap.exon.pops.roadmap.GRCh37.portrait.tiff"
	pheatmap(log(merged.fpkm[rev(rownames(merged.fpkm)),]/10^3+1), cluster_rows=F, cluster_cols=T, annotation_col=ann_tissues, annotation_colors=ann_colors, width=10, height=15, fontsize_col=10, fontsize_row=7, cellheight=9, gaps_col=T, cutree_col=2, filename=my.filename)


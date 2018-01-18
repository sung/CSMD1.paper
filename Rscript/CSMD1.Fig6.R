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
	merged.fpkm<-cbind(RoadMap.ddsFpkm, `Placenta`=rowMeans(POPS.ddsFpkm)) # isa 'matrix'
	colnames(merged.fpkm)<-gsub("_"," ",colnames(merged.fpkm))

	ann_tissues<-data.frame(
					#`Tissue`=factor(c(as.character(RoadMap.colData$Origin), "Placenta-derived")),
					`Tissue`=factor(c(ifelse(grepl("Brain",rownames(RoadMap.colData)), "Brain", ifelse(RoadMap.colData$Origin=="Non-Placenta","Other-somatic",as.character(RoadMap.colData$Origin))),"Placenta-derived")),
					`Source`=factor(c(ifelse(RoadMap.colData$Source=="Schultz","Schultz et al.", ifelse(RoadMap.colData$Source=="BI","Broad Institute",as.character(RoadMap.colData$Source))), "This Study")),
					row.names=colnames(merged.fpkm)
				)
	ann_exons<-data.frame(
					`TSS`=ifelse(rownames(merged.fpkm)=="ENSE00002108991","Canonical TSS",ifelse(rownames(merged.fpkm)=="ENSE00002127078","Placental TSS",NA)),
					row.names=rownames(merged.fpkm)
						  )
	ann_colors<-list(
					 `Tissue`=c(`Other-somatic`=brewer.pal(3, "Dark2")[1], `Placenta-derived`=brewer.pal(3, "Dark2")[2], `Brain`=brewer.pal(3, "Dark2")[3]),
					 `Source`=c(`Schultz et al.`="grey0", `UCSF-UBC`=cbPalette2[2], `Broad Institute`=cbPalette2[5], `This Study`=cbPalette2[3]),
					 `TSS`=c(`Canonical TSS`="#377EB8",`Placental TSS`="#E41A1C")
				)

	# portrait (row; tissues, col: exons)
	my.filename<-"../Figures/Pheatmap.exon/Fig.pheatmap.exon.pops.roadmap.GRCh37.portrait.new.tiff"
	
	pheatmap(log(merged.fpkm[rev(rownames(merged.fpkm)),]/10^3+1), cluster_rows=F, cluster_cols=T, annotation_col=ann_tissues, annotation_row=ann_exons, annotation_colors=ann_colors, width=10, height=15, fontsize_col=10, fontsize_row=7, cellheight=9, gaps_col=T, cutree_col=4, filename=my.filename)

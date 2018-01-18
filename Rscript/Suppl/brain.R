	library(ggbio)

	dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type) # ida 'data.table' for Placenta 

	## CSMD1 Gene +- 1kb (Grch37) ##
	####################
	# transcript track #
	####################
	# hg19ensGene and gr.csmd.ext defined in 'lib/local.R'
	p.tr.reduced<-autoplot(hg19ensGene, which=gr.csmd.ext, stat="reduce", color="black", fill="black") + theme_alignment() 
	p.tss.exon<-autoplot(gr.csmd[c(my.pt.exon_id, my.canon.exon_id)], aes(col=exon_id)) + scale_color_manual(values=c("#377EB8","#E41A1C")) + theme_alignment() + theme(legend.position="none")

	###################
	# Placenta WGoxBS #
	###################
	my.point.size=1.5; my.alpha=0.4
	# dt.query is Ensembl (or NCBI) style coordinate (e.g. 8)
	# change this to UCSC-stype (e.g. chr8)
	dt.wgbs.pt<-rbind(
					`Female`=dt.query[V1==mapSeqlevels( seqlevels(gr.csmd),style="NCBI" ) & V2>=min(start(gr.csmd)) & V2<=max(end(gr.csmd))
						, .(`seqnames`=paste0("chr",V1), `end`=V2,`strand`=V3,`score`=round(V5.x/V6.x*100,2), `Gender`="Female",`Tissue`="Placenta")],
					`Male`=dt.query[V1==mapSeqlevels( seqlevels(gr.csmd),style="NCBI" ) & V2>=min(start(gr.csmd)) & V2<=max(end(gr.csmd))
						, .(`seqnames`=paste0("chr",V1), `end`=V2,`strand`=V3,`score`=round(V5.y/V6.y*100,2), `Gender`="Male", `Tissue`="Placenta")]
					)
	##################
	## RoadMap WGBS ##
	##################
		#`Hippocampus`=rbind(
		#				`Hippocampus.149`=as.data.table(as.data.frame(import.bed("../RData/GSM1112838_BI.Brain_Hippocampus_Middle.Bisulfite-Seq.149.GRCh37.CSMD1.bed.gz"))),
		#				`Hippocampus.150`=as.data.table(as.data.frame(import.bed("../RData/GSM916050_BI.Brain_Hippocampus_Middle.Bisulfite-Seq.150.GRCh37.CSMD1.bed.gz")))
		#				)[,.(`score`=round(mean(score)*100,2),`Gender`="Male",`Tissue`="HC"),"seqnames,end,strand"], # male
	dt.wgbs.brain<-rbind(
		`Hippocampus`=as.data.table(as.data.frame(import.bed("../RData/GSM916050_BI.Brain_Hippocampus_Middle.Bisulfite-Seq.150.GRCh37.CSMD1.bed.gz")))[,.(seqnames,end,strand,score=score*100,`Gender`="Male",`Tissue`="HC")],
		`Brain_Germinal_Matrix`=as.data.table(as.data.frame(import.bed("../RData/GSM941747_UCSF-UBC.Brain_Germinal_Matrix.Bisulfite-Seq.HuFGM02.GRCh37.CSMD1.bed.gz")))[,.(seqnames,end,strand,score=score*100,`Gender`="Male",`Tissue`="GM")] # male
	)

	p.meth<-ggplot(rbind(dt.wgbs.pt, dt.wgbs.brain), aes(end, score, col=Gender)) + 
				geom_point(size=my.point.size,alpha=my.alpha) +
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
				scale_y_continuous(breaks=seq(from=0,to=100,by=25)) +
				labs(y="% Methylation") +
				facet_grid(Tissue~ .) +
				theme_Publication() + 
				theme(legend.position="none")

	####################
	# RNA-Seq Coverage #
	####################
	# GRCh37.75.kjCSMD1.bedgraph.gz: UCSC-style (e.g. chr8) 
	dt.rna.cov.pt<-rbind(
		`Female`=as.data.table(as.data.frame(import.bedGraph("../RData/CTLF.GRCh37.75.CSMD1.bedgraph.gz"))),
		`Male`=as.data.table(as.data.frame(import.bedGraph("../RData/CTLM.GRCh37.75.CSMD1.bedgraph.gz")))
	)[,.(`score`=sum(score)),"seqnames,start,end,width,strand"]

	dt.rna.cov<-rbind(
		`Placenta.Female`=as.data.table(as.data.frame(import.bedGraph("../RData/CTLF.GRCh37.75.CSMD1.bedgraph.gz")))[,.(seqnames,start,end,strand,score,`Gender`="Female",`Tissue`="PT (Female)")],
		`Placenta.Male`=as.data.table(as.data.frame(import.bedGraph("../RData/CTLM.GRCh37.75.CSMD1.bedgraph.gz")))[,.(seqnames,start,end,strand,score,`Gender`="Male",`Tissue`="PT (Male)")],
		#`Hippocampus Middle 149`=as.data.table(as.data.frame(import.bedGraph("../RData/GSM1112836_BI.Brain_Hippocampus_Middle.mRNA-Seq.149.GRCh37.CSMD1.bedgraph.gz"))), # male
		`Hippocampus`=as.data.table(as.data.frame(import.bedGraph("../RData/GSM916094_BI.Brain_Hippocampus_Middle.mRNA-Seq.150.GRCh37.CSMD1.bedgraph.gz")))[,.(seqnames,start,end,strand,score,`Gender`="Male",`Tissue`="HC")], # male
		`Cerebellum`=as.data.table(as.data.frame(import.bedGraph("../RData/GSM1127104_UCSF-UBC.Brain_Cerebellum.mRNA-Seq.HuFGM01.GRCh37.CSMD1.bedgraph.gz")))[,.(seqnames,start,end,strand,score,`Gender`="Male",`Tissue`="CB")] # male
	)

	p.rna.cov<-ggplot(dt.rna.cov[score!=0], aes(xmin=start, xmax=end, ymin=0, ymax=score, col=Gender)) + 
				geom_rect(alpha=0.1) +
				scale_color_manual(name="Sex", values=my.col[["Gender"]]) + 
				labs(y="Average Coverage") +
				facet_grid(Tissue~ .) +
				theme_Publication() + 
				theme(legend.position="none")

	print(tracks(`CSMD1`=p.tr.reduced, `TSS`=p.tss.exon, `Bisulfite Sequencing`=p.meth, `RNA-Seq`=p.rna.cov, heights = c(1.3, 0.8, 9, 12), scale.height=2.5, label.text.cex = 1.3, xlim=gr.csmd.ext))

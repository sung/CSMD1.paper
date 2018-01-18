# part of Rscript/CSMD1.paper.figures.main.R
library(ggbio)

dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type) # ida 'data.table' for PT

## CSMD1 Gene +- 1kb (Grch37) ##
####################
# transcript track #
####################
# hg19ensGene and gr.csmd.ext defined in 'lib/local.R'
p.tr.reduced<-autoplot(hg19ensGene, which=gr.csmd.ext, stat="reduce", color="black", fill="black") + theme_alignment() 
p.tss.exon<-autoplot(gr.csmd[c(my.pt.exon_id, my.canon.exon_id)], aes(col=exon_id)) + scale_color_manual(values=c("#377EB8","#E41A1C")) + theme_alignment() + theme(legend.position="none")

####################
# WGBS methylation #
####################
my.point.size=1.5; my.alpha=0.4
# dt.query is Ensembl (or NCBI) style coordinate (e.g. 8)
# change this to UCSC-stype (e.g. chr8)
dt.csmd.meth<-dt.query[V1==mapSeqlevels( seqlevels(gr.csmd.ext),style="NCBI" ) & V2>=start(gr.csmd.ext) & V2<=end(gr.csmd.ext)
					   , .(`chr`=paste0("chr",V1), `start`=V2,`strand`=V3,`Female`=round(V5.x/V6.x*100,2), `Male`=round(V5.y/V6.y*100,2))]

p.pt.meth<-ggplot(melt(dt.csmd.meth[strand=="-"], id=1:3,measure=4:5,variable.name="Gender",value.name="Meth"), aes(start, Meth, col=Gender)) + 
			geom_point(size=my.point.size,alpha=my.alpha) +
			ylim(0,100) +
			geom_smooth(se=FALSE) + 
			scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
			labs(y="% Methylation") +
			theme_Publication() + 
			theme(legend.position="none")
####################
# RNA-Seq Coverage #
####################
# GRCh37.75.CSMD1.bedgraph.gz: UCSC-style (e.g. chr8) 
rna.cov.list<-list(
	`Female`=as.data.table(as.data.frame(import.bedGraph("../RData/CTLF.GRCh37.75.CSMD1.bedgraph.gz"))),
	`Male`=as.data.table(as.data.frame(import.bedGraph("../RData/CTLM.GRCh37.75.CSMD1.bedgraph.gz")))
)
dt.rna.cov<-rbindlist(lapply(names(rna.cov.list), function(i) rna.cov.list[[i]][,`Gender`:=i]))

my_label<-list(`Female`="Female (n=64)", `Male`="Male (n=67)")
my_labeller<-function(variable,value){
	return(my_label[value])
}   
p.rna.cov<-ggplot(dt.rna.cov[score!=0], aes(xmin=start, xmax=end, ymin=0, ymax=score, col=Gender, fill=Gender),alpha=0.1) + 
			geom_rect(alpha=0.1) +
			scale_fill_manual(name="Sex", values=my.col[["Gender"]]) +
			scale_color_manual(name="Sex", values=my.col[["Gender"]]) +
			labs(y="Average Coverage") +
			theme_Publication() + 
			theme(legend.position="none") + 
			facet_grid(Gender~ ., labeller=my_labeller)

my.file.name<-"../Figures/CSMD1.Meth.Exp.profile/Fig.CSMD1.methy.exp.profile.TSS.minus.strand.tiff"
tiff(filename=my.file.name,width=15, height=8.5 ,units="in",res=300, compression = 'lzw') #A4 size 
print(tracks(`CSMD1`=p.tr.reduced, `TSS`=p.tss.exon, `WGoxBS (n=4)`=p.pt.meth, `RNA-Seq (n=131)`=p.rna.cov, heights = c(1.3, 0.8, 3.5, 7), scale.height=2.5, label.text.cex = 1.3, xlim=gr.csmd.ext))
dev.off()

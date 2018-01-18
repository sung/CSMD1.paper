library(ggplot2)
library(GenomicFeatures)
library(rtracklayer)
library(ggbio)
library(data.table)

source("~/Pipelines/lib/theme_publish.R") # defines theme_Publication()

gr.csmd<-import("~/results/CSMD1.2016.paper/RData/Homo_sapiens.GRCh37.75.CSMD1.exon.union.sorted.ucsc.gtf") # ucsc style (e.g. chr8)
gr.csmd.ext<-promoters(range(reduce(gr.csmd)), upstream=1000, downstream=end(range(reduce(gr.csmd)))-start(range(reduce(gr.csmd)))+1+1000) # includes 1kb +TSS and 1kb +TES 

my.local.db<-"~/results/CSMD1.2016.paper/RData/hg19ensGene.sqlite" #annotation db
hg19ensGene<- loadDb(my.local.db) # from AnnotationDbi via GenomicFeatures
p.tr.reduced<-autoplot(hg19ensGene, which=gr.csmd.ext, stat="reduce", color="brown", fill="brown") + theme_alignment() 

####################
# WGBS methylation #
####################
my.point.size=1.5; my.alpha=0.5

dt.meth.brain.csmd1<-as.data.table(as.data.frame(import.bed("/home/ssg29/data/RoadMap/Consortium/Bisulfite-Seq/Brain/GSM1112838_BI.Brain_Hippocampus_Middle.Bisulfite-Seq.149.GRCh37.CSMD1.bed.gz")))
dt.exp.brain.csmd1<-as.data.table(as.data.frame(import.bedGraph("/home/ssg29/data/RoadMap/Consortium/RNA-Seq/Somatic/GSM1112836_BI.Brain_Hippocampus_Middle.mRNA-Seq.149.GRCh37.CSMD1.bedgraph.gz")))

p.brain.meth<-ggplot(dt.meth.brain.csmd1, aes(start, score*100)) + 
			geom_point(size=my.point.size,alpha=0.1) +
			ylim(0,100) +
			geom_smooth(se=FALSE) + 
			labs(y="% Methylation") +
			theme_Publication() + 
			theme(legend.position="none")

p.brain.exp<-ggplot(dt.exp.brain.csmd1, aes(xmin=start, xmax=end, ymin=0, ymax=score)) + 
			geom_rect(col="red",alpha=0.001) +
			labs(y="Average Coverage") +
			theme_Publication() + 
			theme(legend.position="none") 

print(tracks(`CSMD1`=p.tr.reduced, `Bisulfite-Seq`=p.brain.meth, `RNA-Seq`=p.brain.exp, heights = c(1.2,3.5,3.5), scale.height=2.5, label.text.cex = 1.3, xlim=gr.csmd.ext))

#!/usr/bin/Rscript --vanilla
# used by CSMD1.paper.figures.main.R
library(data.table)
library(GenomicFeatures)
library(biomaRt)
library(rtracklayer)
library(RColorBrewer) # for brewer.pal
library(ggplot2)
library(reshape2)

source("lib/theme_publish.R") # defines theme_Publication()
source("lib/multiplot.R") 

time.stamp <- format(Sys.time(), '%Y-%m-%d_%I_%M_%p') #2015-01-07_01_30PM
TR_PREFIX="GRCh37"
myGenome="hg19" # for goseq (run supportedGenomes())
ENS_VER=75 # 75 (for igenome), 82 (for manual downlaod)
my.local.db<-"../RData/hg19ensGene.sqlite" #annotation db
hg19ensGene<- loadDb(my.local.db) # from AnnotationDbi via GenomicFeatures
									# makeTxDbFromUCSC(genome = "hg19", tablename = "ensGene")
#biomart
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
myMart = grch37
gr.ensg<-genes(loadDb(my.local.db)) # from GenomicFeatures

# manual chr order (for geom_box)
my.chr.order=c(seq(1:22),"X")

# A colorblind-friendly palette
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#662506", "#997a00", "#800080", "#000000")
# The palette with black:
cbPalette2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

my.col=list(
	`Sex`=c(`M`="#D95F02", `F`="#1B9E77"), # orange/green (RnBeads default) for RNA-Seq data
	`Gender`=c(`Male`="#D95F02", `Female`="#1B9E77"), # orange/green (RnBeads default) for RoadMap via ggplot
)

my.cpg.type="CG" # CG, CHH, CHG
my.tissue="PT" # short for placenta
min.doc <- 10 # minimum depth of coverage
my.ensg='ENSG00000183117' # CSMD1
my.first.exon_id="ENSE00002127078" # the first exon at the putative TSS

gr.csmd<-import("../RData/Homo_sapiens.GRCh37.75.CSMD1.exon.union.sorted.ucsc.gtf") # ucsc style (e.g. chr8)
gr.csmd.ext<-promoters(range(reduce(gr.csmd)), upstream=1000, downstream=end(range(reduce(gr.csmd)))-start(range(reduce(gr.csmd)))+1+1000) # includes 1kb +TSS and 1kb +TES 

# UCSC-style (e.g. chr8)
target.region<-list(
	`CSMD1.DMR`="../RData/CSMD1.boy.girl.dmr.bed", # Original CSMD1 DMR
	`CSMD1.DMR2`="../RData/CSMD1.boy.girl.dmr+2exons.bed", # above + 5' 2 exons 
	`CSMD1.legacy.promo`="../RData/CSMD1.legacy.promo.4K.GRCh37.bed",
	`CSMD1.putative.promo`="../RData/CSMD1.placenta.putative.promoter.4K.GRCh37.bed"
)	

# Load PT DNA methylation data 
load.my.tissue.dt.merged<-function(my.tissue,my.cpg.type){
	## Load this RoadMap tissue
	## e.g. RData/PT.CG.dt.merged.RData
	my.new.RData<-file.path("../RData/", paste(my.tissue,my.cpg.type,"dt.merged.RData",sep="."))
	cat(paste0("Loading ", my.tissue, " ", my.cpg.type, "...\n"))
	print(system.time(load(my.new.RData))) # load list of data.table
	cat(paste(my.tissue,my.cpg.type,"dt.merged loaded\n",sep="."))
	this.tissue<-get(paste(my.tissue,my.cpg.type,"dt.merged",sep=".")) # isa list of data.table (STL2, STL3)
	eval(parse(text=paste0('rm(', my.tissue,'.',my.cpg.type,'.dt.merged)'))) # remove the original object from memory

	## Merge all chromosomes
	if(my.tissue!='PT'){
		this.tissue[["chrY"]] <- NULL # remove chrY
		this.tissue[["chrL"]] <- NULL # remove chrL
		cat("rbindlist for this.tissue...\n")
		print(system.time(dt.query<-rbindlist(this.tissue)))
	# Placenta
	}else{
		dt.query<-this.tissue
	}
	rm(this.tissue)
	#>dt.query
	#   V1    V2 V3  V4 V7 V5.x V6.x V5.y V6.y
	#1: 10 66269  + CGG  1   12   12   44   51
	#2: 10 66273  + CGG  1   14   14   36   53
	#3: 10 71723  + CGA  1    3   10    6   16
	dt.query[,End:=V2] # add this column
	#>dt.query
	#   V1    V2 V3  V4 V7 V5.x V6.x V5.y V6.y   End
	#1: 10 66269  + CGG  1   12   12   44   51 66269
	#2: 10 66273  + CGG  1   14   14   36   53 66273
	#3: 10 71723  + CGA  1    3   10    6   16 71723

	cat("some integer columns as numeric...\n")
	this.col=c("V5.x","V6.x","V5.y","V6.y")
	print(system.time(dt.query[, (this.col) := lapply(.SD,as.numeric), .SDcols=this.col]))
	return(dt.query)
}

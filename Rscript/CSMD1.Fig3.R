
	library(qqman)
	####################
	## Figure 3       ## 
	## Manhattan Plot ##
	#####################
	rnb.dm.rdata<-"../RData/rnb.meth.diff.PT.RData"
	cat(paste0("Loading rnb.dm.rdata: ",rnb.dm.rdata,"...\n"))
	load(rnb.dm.rdata) # loading "rnb.diff.met"
	cat("rnb.diff.met loaded\n") # isa list of data.frame

	# A. tiling 5K
	my.diff.met<-rnb.diff.met[["tiling"]]; select<-!is.na(my.diff.met$combinedRank) & my.diff.met$num.sites >= median(my.diff.met$num.sites)
	my.diff.met<-my.diff.met[select,]
	upper<-log(max(my.diff.met$combinedRank))
	t100.tiling<-upper-log(my.diff.met[100,"combinedRank"]) # combinedRank of top100
	top.diff.met<-apply(my.diff.met[my.diff.met$Chromosome=='chr8' & my.diff.met$Start>=2795000 & my.diff.met$End<=3020000,][1:20,c("Chromosome","Start","Strand")], 1, paste, collapse=".") # 20 CSMD1 5K regions
	p.tiling5k<-with(my.diff.met, 
				data.frame(
					SNP=paste(my.diff.met$Chromosome, my.diff.met$Start, my.diff.met$Strand, sep="."), # e.g. chr7.64540001.*
					CHR=sapply(substr(as.character(Chromosome),4,nchar(as.character(Chromosome))), # chromosome (e.g. chr1, chr2, chrX, chrM)
						function(i) if(i=="X"){i<-23}else if(i=="Y"){i<-24}else if(i=="M"){i<-25}else{i<-as.numeric(i)}),
					BP=Start,
					RANK=upper-log(combinedRank)
				))

	# B. genes 
	my.diff.met<-rnb.diff.met[["genes"]]; select<-!is.na(my.diff.met$combinedRank) & my.diff.met$num.sites >= median(my.diff.met$num.sites)
	my.diff.met<-my.diff.met[select,]
	upper<-log(max(my.diff.met$combinedRank))
	t100.genes<-upper-log(my.diff.met[100,"combinedRank"]) # combinedRank of top100
	p.genes<-with(my.diff.met, 
				data.frame(
					SNP=paste(my.diff.met$Chromosome, my.diff.met$Start, my.diff.met$Strand, sep="."), # e.g. chr7.64540001.*
					CHR=sapply(substr(as.character(Chromosome),4,nchar(as.character(Chromosome))), # chromosome (e.g. chr1, chr2, chrX, chrM)
						function(i) if(i=="X"){i<-23}else if(i=="Y"){i<-24}else if(i=="M"){i<-25}else{i<-as.numeric(i)}),
					BP=Start,
					RANK=upper-log(combinedRank)
				))

	# C. promoters 
	my.diff.met<-rnb.diff.met[["promoters"]]; select<-!is.na(my.diff.met$combinedRank) & my.diff.met$num.sites >= median(my.diff.met$num.sites)
	my.diff.met<-my.diff.met[select,]
	upper<-log(max(my.diff.met$combinedRank))
	t100.promoters<-upper-log(my.diff.met[100,"combinedRank"]) # combinedRank of top100
	p.promoters<-with(my.diff.met, 
				data.frame(
					SNP=paste(my.diff.met$Chromosome, my.diff.met$Start, my.diff.met$Strand, sep="."), # e.g. chr7.64540001.*
					CHR=sapply(substr(as.character(Chromosome),4,nchar(as.character(Chromosome))), # chromosome (e.g. chr1, chr2, chrX, chrM)
						function(i) if(i=="X"){i<-23}else if(i=="Y"){i<-24}else if(i=="M"){i<-25}else{i<-as.numeric(i)}),
					BP=Start,
					RANK=upper-log(combinedRank)
				))

	# D. cpgislands 
	my.diff.met<-rnb.diff.met[["cpgislands"]]; select<-!is.na(my.diff.met$combinedRank) & my.diff.met$num.sites >= median(my.diff.met$num.sites)
	my.diff.met<-my.diff.met[select,]
	upper<-log(max(my.diff.met$combinedRank))
	t100.cpgislands<-upper-log(my.diff.met[100,"combinedRank"]) # combinedRank of top100
	p.cpgislands<-with(my.diff.met, 
				data.frame(
					SNP=paste(my.diff.met$Chromosome, my.diff.met$Start, my.diff.met$Strand, sep="."), # e.g. chr7.64540001.*
					CHR=sapply(substr(as.character(Chromosome),4,nchar(as.character(Chromosome))), # chromosome (e.g. chr1, chr2, chrX, chrM)
						function(i) if(i=="X"){i<-23}else if(i=="Y"){i<-24}else if(i=="M"){i<-25}else{i<-as.numeric(i)}),
					BP=Start,
					RANK=upper-log(combinedRank)
				))

	# A-D in multi-pannel
	file.name<-file.path("../Figures/Manhattan",paste("Fig3.manhattan",time.stamp,"ALL.tiff",sep="."))
	tiff(filename=file.name,width=10, height=8,units="in",res=300, compression = 'lzw') #A4 size 

	par(mfrow=c(2,2))
	par(mar=c(4,4.5,4,1))
	manhattan(p.tiling5k, p="RANK", logp=FALSE, ylab="Rank", xlab="", genomewideline=FALSE, suggestiveline=FALSE, highlight=top.diff.met, chrlabs=c(1:22, "X"), main="5Kbp tiles", font.lab=2,  ylim=c(0,9), cex=1.5, cex.main=1.8, cex.lab=1.5, cex.axis=.8, yaxt='n') 
	abline(h=t100.tiling) 
	# https://gist.github.com/cboettig/709578
	mtext("A", side=3, adj=0, line=1.2, cex=1.8, font=2) # side = 3 means top margin; adj = 0 means left align; line= 1.2 moves the text up; font=2 to bold face
	axis(2, at=seq(0,9), las=1)

	par(mar=c(4,4.5,4,1))
	#mtext(adj=0)
	manhattan(p.genes, p="RANK", logp=FALSE, ylab="", xlab="", genomewideline=FALSE, suggestiveline=FALSE, chrlabs=c(1:22, "X"), main="Gene-bodies", font.lab=2, ylim=c(0,5), cex=1.5, cex.main=1.8, cex.lab=1.5, cex.axis=.8) 
	abline(h=t100.genes)
	mtext("B", side=3, adj=0, line=1.2, cex=1.8, font=2)

	par(mar=c(4,4.5,4,1))
	manhattan(p.promoters, p="RANK", logp=FALSE, ylab="Rank", genomewideline=FALSE, suggestiveline=FALSE, chrlabs=c(1:22, "X"), main="Promoters", font.lab=2, ylim=c(0,7), cex=1.5, cex.main=1.8, cex.lab=1.5, cex.axis=.8) 
	abline(h=t100.promoters)
	mtext("C", side=3, adj=0, line=1.2, cex=1.8, font=2)

	par(mar=c(4,4.5,4,1))
	manhattan(p.cpgislands, p="RANK", logp=FALSE, ylab="", genomewideline=FALSE, suggestiveline=FALSE, chrlabs=c(1:22, "X"), main="CpG islands", font.lab=2, ylim=c(0,8), cex=1.5, cex.main=1.8, cex.lab=1.5, cex.axis=.8, yaxt='n') 
	abline(h=t100.cpgislands)
	mtext("D", side=3, adj=0, line=1.2, cex=1.8, font=2)
	axis(2, at=seq(0,8), las=1) # las=1: text horizontal 
	dev.off()


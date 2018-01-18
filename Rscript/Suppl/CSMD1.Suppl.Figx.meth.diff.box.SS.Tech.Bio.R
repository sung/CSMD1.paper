
	dt.meth.list<-list()
	load(file.path("~/results/RoadMap/BS-Seq/RData/10X",paste("PT.SS.Tech",my.cpg.type,"dt.meth.region.RData",sep="."))) # load dt.meth.region
	dt.meth.list[[1]]<-dt.meth.region # same as dt.meth.region[cnt.C+cnt.T>=10]
	load(file.path("~/results/RoadMap/BS-Seq/RData/10X",paste("PT.SS.Bio",my.cpg.type,"dt.meth.region.RData",sep="."))) # load.dt.meth.region
	dt.meth.list[[2]]<-dt.meth.region # dt.meth.region[cnt.C+cnt.T>=30]

	dt.meth<-rbindlist(dt.meth.list)
	dt.meth[,`Chromosome Type`:=ifelse(chr=="X",'X-chromosome','Autosomes')]
	dt.meth.gender=data.table(dcast(dt.meth, Tissue+`Chromosome Type`+Region+chr+strand+start+end~Gender, value.var="meth")) 

	p.box<-ggplot(dt.meth.gender, aes(Region, (Female-Male)*100)) + 
				geom_boxplot(aes(fill=`Chromosome Type`), outlier.shape=NA, alpha=.7, size=.7, width=.5) + 
				labs(x="Genomic Regions",y="% Methylation Difference\n   (more in male)                     (more in female)") + 
				geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
				coord_cartesian(ylim=c(-50, 50)) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
				scale_y_continuous(breaks=seq(from=-50,to=50,by=10)) +
				scale_x_discrete(limits=c("Sites","Genes","Promo_15.05","CPGi"),labels=c("CpG sites","Gene-bodies","Promoters","CPG islands")) +
				scale_fill_manual(values=my.col[["Chromosome Types"]]) +
				facet_grid(.~Tissue, labeller=labeller(Tissue=c(`PT.SS.Bio`="SureSelect Bio-Rep (n=6)",`PT.SS.Tech`="SureSelect Tech-Rep (n=2)"))) +
				theme_Publication() +
				theme(legend.position="top", plot.title = element_text(hjust = 0, size = 15,face="bold"))

	file.name<-file.path("~/results/CSMD1.2016.paper/Figures/Suppl",paste("Suppl.Meth.diff.PT.SS.Tech.",time.stamp,"tiff",sep="."))
	tiff(filename=file.name,width=14, height=10,units="in",res=300, compression = 'lzw') #A4 size 
	print(p.box)
	dev.off()


	##
	dcast(dt.meth[,.(`Cs`=sum(cnt.C),`Ts`=sum(cnt.T),`Depth`=sum(cnt.C+cnt.T), `Meth`=round(sum(cnt.C)/sum(cnt.C+cnt.T)*100,2)), "Tissue,Gender,Region"], Tissue+Region~Gender,value.var="Meth")
	dcast(dt.meth[,.(`Cs`=sum(cnt.C),`Ts`=sum(cnt.T),`Depth`=sum(cnt.C+cnt.T), `Meth`=round(sum(cnt.C)/sum(cnt.C+cnt.T)*100,2)), "Chromosome Type,Tissue,Gender,Region"], `Chromosome Type`+Tissue+Region~Gender,value.var="Meth")

	##
	dt.meth.tissue.list<-list()
	load(file.path("~/results/RoadMap/BS-Seq/RData/10X",paste(my.cpg.type,"dt.meth.per.tissue.per.region.RData",sep="."))) # load.dt.meth.region
	dt.meth.tissue.list[[1]]<-rbindlist(dt.meth.per.tissue)
	load(file.path("~/results/RoadMap/BS-Seq/RData/10X",paste(my.cpg.type,"PT.SS.dt.meth.per.tissue.per.region.RData",sep="."))) # load.dt.meth.region
	dt.meth.tissue.list[[2]]<-rbindlist(dt.meth.per.tissue)

	dt.foo<-rbindlist(dt.meth.tissue.list)[context %in% c("Sites","Genes","Promo_15.05","CPGi")]
	dt.foo[,`chr.type`:=ifelse(V1=="X","X",ifelse(V1=="genome",NA,"autosomes"))]

	dt.bar<-dt.foo[,.(`No.CpG`=sum(num.sites),`Female.Cs`=round(sum(f.c)),`Female.Ts`=round(sum(f.t)),`Female.Meth`=round(sum(f.c)/sum(f.c+f.t)*100,2),`Male.Cs`=round(sum(m.c)),`Male.Ts`=round(sum(m.t)),`Male.Meth`=round(sum(m.c)/sum(m.c+m.t)*100,2), `Meth.diff`=round(sum(f.c)/sum(f.c+f.t)*100-sum(m.c)/sum(m.c+m.t)*100,2)), "chr.type,tissue,context"]
	dt.bar[is.na(chr.type),`chr.type`:="genome"]

	write.csv(merge(dt.bar, data.table(tissue=names(unlist(tissue.list)), `Tissue Name`=unlist(tissue.list)),all.x=TRUE), file="~/results/CSMD1.2016.paper/RData/meth.diff.by.tissues.csv", row.names=F, quote=F)

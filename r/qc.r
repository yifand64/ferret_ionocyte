# Given a QC file for a set of samples, produce plots that 
# check whether PCA or cluster membership is correlated
# strongly with QC (and is not biological variation)


# library(scales)
# library(reshape2)
# library(grid)
# library(gridExtra)
# library(ggplot2)

info("Loading quality control (QC) functions")

# for each cluster, calculate and visualise
# average
check.qc.in.clusters = function(qc, clusters)
{
	info("Checking the QC metrics for the given cluster partitioning: ")
	print(table(clusters))

	print("Cluster vector length:")
	print(length(clusters))

	print("QC table size:")
	print(dim(qc))



	# aggregate the QC data by taking the mean of every cluster:
	info("Aggregating data using the mean of each cluster..")
	aggdata <- aggregate(qc, by=list(clusters), FUN=mean, na.rm=TRUE)
	info("Done. ")
	print(head(aggdata))
	aggdata$Sample = aggdata$Group.1
	# print(aggdata)
	#info("Aggregated:")
	#print(aggdata)

	visualise.qc.summary(aggdata, 
						legend.key.title="Cluster",
						legend.font.size=12,
						output.file="qc_in_clusters.pdf")
	cat(sprintf("Cluster --  QC  check complete.\n"))
}

# check correlation of PCs with QC metrics
# if the correlation is above a threshold,
# a scatter plot is produced. input is a 
# pca object output by prcomp.
check.pcs = function(qc, 
					pca,
					groups,
					counts=NULL,
					add.count.string=F,
					n_components=4)
{
	
	initial_wd = getwd()
	qc_check_dir = "QC_Correlations"
	dir.create(qc_check_dir, showWarnings = FALSE)
	setwd(qc_check_dir)

	scores <- data.frame(groups, pca$x[,1:n_components])
	scores$sample_name <- rownames(scores) 
	

	if(!(is.null(qc)))
	{
		cat(sprintf("Checking correlation of PCs with QC metrics: \n"))
		# cat(sprintf("QC data: \n"))
		# print(head(qc))
		metrics = c("X.rRNA", "gene_count", "trans_mapped", "X.trans_mapped")
		metric_names = c("X.rRNA"="% Ribosomal RNA", "gene_count"="Genes Detected", "trans_mapped"="# Mapped Reads (Transcriptome)", "X.trans_mapped"="% Mapped Reads (Transcriptome")
		qc = qc[c("Sample", metrics)]
		if(add.count.string)
		{
			qc$Sample = paste("normalised_count", qc$Sample, sep="_")
		}
		colnames(qc) = c("sample_name", metrics)
		qc = merge(qc, scores, by="sample_name")
		
	}else
	{
		if(is.null(counts))
		{
			error("Must provide counts or QC!")
			return(FALSE)
		}
		cat(sprintf("No QC provided! \n"))
		cat(sprintf("Checking correlation of PCs with simple QC metrics: \n"))
		metrics = c("Total.Transcripts", "Genes.Detected")
		metric_names = c("Total.Transcripts" = "Total Transcripts", "Genes.Detected"="Genes Detected")
		qc = scores
		qc$Total.Transcripts = colSums(counts)
		qc$Genes.Detected = colSums(counts>0)
		print(head(qc))
	}


	for(i in 1:length(metrics))
	{
		qc_metric = metrics[[i]]
		for(j in 1:n_components)
		{
			pc = paste("PC",j, sep="")
			cat(sprintf("Checking correlation of %s with %s .. \n", qc_metric, pc))
			pearson = round(cor(qc[qc_metric], qc[pc], method="pearson"), 3)
			spearman = round(cor(qc[qc_metric], qc[pc], method="spearman"), 3)
			cat(sprintf("Spearman = %f, Pearson = %f \n", spearman, pearson))
			if(abs(spearman) > 0.5 | abs(pearson) > 0.5)
			{
				warn("Possible correlation. Printing scatter plot..")
				n = length(unique(qc$groups))
				correlation_label = paste("Pearson = ", pearson,"\n", "Spearman = ", spearman ,sep="")
				p = ggplot(qc, aes_string(x=qc_metric, y=pc, colour=groups)) + geom_point() + theme_bw() + 
					ggtitle(correlation_label) + xlab(metric_names[qc_metric]) +scale_colour_manual("", values=default.cols(n))
				ggsave(p, file=paste("Correlation", qc_metric,pc, ".pdf",sep="_"), height=6, width=6)
				# scale_x_continuous(xlab, breaks =  seq(xmin, xmax, 2), limits=c(xmin, xmax)) + 
				# scale_y_continuous(ylab, breaks =  seq(ymin, ymax, 2), limits=c(ymin, ymax)) 
				#p = p + annotate(geom="text", x=max(qc[qc_metric])/2, y=max(qc[pc])/2, label=correlation_label, size=3)
		
			
			}
		}
	}
	cat(sprintf("PC --  QC  check complete.\n"))
	setwd(initial_wd)
}


cluster.complexity <- function(seurat.obj, 
								clusters=NULL, 
								use.only=NULL, #specify the name of specific clusters to include
								use.cols=NULL, 
								pdf.output=F, 
								pdf.name="complexity.and.cell.type.pdf")
{
	info("Drawing cluster complexity violin plot..")
	
	if(pdf.output)
	{
		pdf(pdf.name)	
	}
	if(!is.null(use.only))
	{
		cells.use = which(clusters %in% use.only)
		clusters = clusters[cells.use]
		

	}
	x = seurat.obj@data.info
	if(!is.null(clusters))
	{
		x["clust.ident"] = clusters
	}else
	{
		x["clust.ident"] = x$DBclust.ident
	}
	print(colnames(x))
	colourCount = length(unique(x$clust.ident))

	getPalette = colorRampPalette(brewer.pal(9, "Set1"))
	if(is.null(use.cols))
	{
		use.cols = getPalette(colourCount)
	}
	grouped.by.mean = data.table(x)[,list(nGene=mean(nGene), 
		gene_count=mean(gene_count), gene_count=mean(gene_count), Est.Lib.Size=mean(Estimated.Library.Size)),by=clust.ident]

	use.cols = use.cols[order(grouped.by.mean$gene_count)]

	g = ggplot(x, aes(y=gene_count, x=reorder(as.factor(clust.ident), gene_count), fill=as.factor(clust.ident))) + 
		geom_violin() + theme_bw() + scale_fill_manual(values = use.cols, guide=FALSE)  + 
		ylab("Number of Genes detected") + xlab("Cell-type / Cluster") + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
	print(g)

	g = ggplot(x, aes(y=gene_count, x=reorder(as.factor(clust.ident), gene_count), fill=as.factor(clust.ident))) + 
		geom_boxplot() + theme_bw() + scale_fill_manual(values = use.cols, guide=FALSE)  + 
		ylab("Number of Genes detected") + xlab("Cell-type / Cluster") + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
	print(g)

	info("Drawing mean cluster complexity bar plot..")
	
	#print(grouped.by.mean)

	grouped.by.mean = data.frame(grouped.by.mean)
	g = ggplot(grouped.by.mean, aes(y=gene_count, x=reorder(as.factor(clust.ident), gene_count),fill=as.factor(clust.ident))) + 
		geom_bar(stat="identity") + theme_bw() + scale_fill_manual(values = use.cols, guide=FALSE)  + 
		ylab("Mean number of genes detected") + xlab("Cell-type / Cluster") + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
	print(g)

	if(pdf.output)
	{
		dev.off()
	}
	info("Done.")
}

# run clustering on the QC results to see which
# samples have the most similar/dissimilar QC
# profiles
cluster.qc = function(qc)
{
	# K-Means Clustering with 5 clusters
	fit <- kmeans(mydata, 5)

	# Cluster Plot against 1st 2 principal components

	# vary parameters for most readable graph
	clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE, 
	  	labels=2, lines=0)

	# Centroid Plot against 1st 2 discriminant functions
	library(fpc)
	plotcluster(mydata, fit$cluster)
}

# compare the clustering obtained by clustering
# the counts, and that by clustering the QC.
# If all the variation in the dataset is biological
# there should be no agreement. 
cluster.compare.qc.counts = function(qc, counts)
{
	# comparing 2 cluster solutions
	cluster.stats(d, fit1$cluster, fit2$cluster)
}

#possibly useful function:

grid_arrange_shared_legend <- function(...) {
    plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    grid.arrange(
        do.call(arrangeGrob, lapply(plots, function(x)
            x + theme(legend.position="none"))),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight))
}


reads.mapping.bar <- function(seurat.obj)
{
	if(length(unique(seurat.obj@data.info$orig.ident)) == 1)
	{
		warn("Only one batch, cannot draw mapping quality plot")
		return (FALSE)
	}
	source("marker.r")
	library(RColorBrewer)
	colnames(seurat.obj@data.info) = gsub("X.genome_mapped", "Genome", colnames(seurat.obj@data.info))
	colnames(seurat.obj@data.info) = gsub("X.rRNA", "Ribosomal.RNA", colnames(seurat.obj@data.info))
	colnames(seurat.obj@data.info) = gsub("X.trans_mapped", "Transcriptome", colnames(seurat.obj@data.info))
	pdf("BatchMappingQuality.pdf", width=14, height=8)
	bar(seurat.obj, group.by="orig.ident",
		features.plot=c("Genome", "Ribosomal.RNA", "Transcriptome"), stack.by.gene = F, 
		use.cols = brewer.pal(3, "Set3"), ylab=c("% of reads"), legend.title = "")
	dev.off()
}


# draw a bunch of barplots showing the quantities of QC
# metrics in a set of samples:
visualise.qc.summary <- function(data, 
									output.file="qc.summary.pdf",
									convert.to.numeric=T,
									single_page = FALSE,
									legend.key.title="Sample",
									legend.font.size = 3,
									title = "QC metrics")
{
	#ensure data is numeric:
	info("Generating QC plots..")
	if(convert.to.numeric)
	{
		info("Converting all necessary data to numeric")
		asNumeric <- function(x) as.numeric(levels(x))[x]#as.character(x))
		factorsNumeric <- function(d) modifyList(d, lapply(d[, sapply(d, is.factor)], asNumeric))
		cols = c("Intragenic.Rate", "Exonic.Rate","Intronic.Rate", "Intergenic.Rate","X.genome_mapped","X.trans_mapped","X.rRNA", "gene_count","fastq_frag_count","Estimated.Library.Size")
		formatted = c("% Intragenic", "% Exonic ","% Intronic", "% Intergenic","% Genome mapped","% Transcriptome mapped","% ribsomal RNA", "Gene Count", "Read Count", "Estimated Library Size")
		data[cols] <- factorsNumeric(data[cols])
	}
	data.keeping = na.omit(data[ , c("Sample", cols)])
	# print(data.keeping)
	info("Data dimensions: ")
	print(dim(data.keeping))

	
	# info("Melting data..")
	# long <- melt(data.keeping, id.vars=1)
	# info("Melt complete.")
	# info("Data dimensions: ")
	# print(dim(long))
	# print(long)
	# print(table(long$variable))

	pdf(output.file)
	# draw the bar plot
	#ggplot(long,aes(x = Sample, y = value)) + geom_bar(aes(fill = variable), stat = "identity",position = "dodge", width=0.5) + scale_y_discrete(breaks = c(0, 0.25,0.5,0.75,1))+
	#    scale_fill_brewer(palette="Set1", name="QC \n metric", labels=c(cols)) + theme_bw() + ggtitle(title) + scale_x_discrete("QC metrics"

	#print(data)
	colorss <- rev(colorRampPalette(brewer.pal(9,"YlGnBu"))(nrow(data)*1.2))
	ps = c()

	# some density histograms:
	q = data[, c("X.rRNA", "gene_count", "X.trans_mapped", "trans_mapped")]
	d = melt(q)
	levels(d$variable)= c("% Ribosomal RNA", "Genes Detected", "Transcriptome Mapping %", "Reads (transcriptome mapped)")
	g=ggplot(d, aes(x=value)) + geom_density() + facet_wrap(~variable, scale="free") + theme_bw()
	print(g)

	for (i in 1:length(cols))
	{
		metric_name = cols[i]
		# converting percentage data to fraction
		cat(sprintf("\"Printing bar graph for --> %s\"\n", metric_name))
		#print(data[metric_name])
		max = max(data[metric_name], na.rm=TRUE)
		print(max)
		if( max <= 1 )
		{
			print("Converting to percent..")
			data[metric_name] = data[metric_name] * 100
			#print("Converted.")
			#print (data[metric_name])
		}

		if (max > 100)
		{
			max_y=max
		}else
		{
			max_y=100
		}

		print("N- samples:")
		n_samples = length(data$Sample)
		print(n_samples)
		cat(sprintf("N- values [%s]:\n", metric_name))
		print(length(data[,metric_name]))

		if(n_samples > 200)
		{
			clean.name = gsub(".", " ", metric_name, fixed=T)
			clean.name = gsub("X", "%", clean.name, fixed=T)
			clean.name = gsub("_", " ", clean.name, fixed=T)

			warn(sprintf("WARN: Too many samples for bar plot. Drawing a histogram [%s]", clean.name))
			 p = ggplot(data, aes_string(x=metric_name))
			 p = p + geom_histogram() + theme_bw() + ylab("Number of samples/cells") + xlab(clean.name) +
				theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
				theme(axis.line = element_line(colour = "black"),
				    panel.grid.major = element_blank(),
				    panel.grid.minor = element_blank(),
				    panel.border = element_blank(),
				    panel.background = element_blank(), 
				    strip.background = element_blank(), 
				    strip.text.x = element_text(size = 20, face="bold"),
				    text = element_text(size=20, colour="gray22"))
			#p = p + scale_y_continuous(formatted[i], limits=c(0,max_y), breaks=seq(0,max_y,by=(max_y/20)))
			
		}else{
			cat(sprintf("Generating bar plot (ggplot)..\n"))
			p = ggplot(data, aes_string(x="Sample", y=metric_name, fill="Sample")) + geom_bar(stat="identity") + theme_bw() + scale_fill_manual(values=colorss, name=legend.key.title)
			p =p + theme(axis.line=element_blank(),
		      axis.text.x=element_blank(),
		      #axis.text.y=element_blank(),
		      axis.ticks.x=element_blank(),
		      axis.title.x=element_blank(),
		      #axis.title.y=element_blank(),
		      #legend.position="none",
		      #panel.background=element_blank(),
		      panel.border=element_blank())
		      #panel.grid.major=element_blank(),
		      #panel.grid.minor=element_blank(),
		      #plot.background=element_blank())
			p = p + scale_y_continuous(formatted[i], limits=c(0,max_y), breaks=seq(0,max_y,by=(max_y/20)))
			p = p + theme(axis.text.y = element_text(size=6))
			p = p + theme(legend.key.width=unit(0.1,"cm"), 
				legend.key.height=unit(0.1,"cm"), 
				legend.text = element_text(size = legend.font.size))
		}
 
		
		
		if(single_page){
			ps[[i]] = p
		}else{
			print(p)
		}
	}
	if (single_page){
		grid_arrange_shared_legend(ps)
	}
	#do.call(grid.arrange, ps)

	#look for correlation between transcriptome mapping rate and intergenic / intronic reads:
	cat(sprintf("Generating scatter between transcriptome mapping rate and intergenic..\n"))
	#
	print(ggplot(data, aes_string(x="Intergenic.Rate", y="X.trans_mapped"), label="Sample") + geom_point() + geom_text(label=data$Sample, vjust=3, size=1.5)+ theme_bw())

	cat(sprintf("Generating scatter between transcriptome mapped reads and number of genes [main filtering parameters]..\n"))
	#

	data["Log2.trans_mapped"] = log2(data$trans_mapped + 1)
	print(ggplot(data, aes_string(x="Log2.trans_mapped", y="gene_count"), label="Sample") + geom_point() + geom_text(label=data$Sample, vjust=3, size=1.5)+ theme_bw() + ylab("Gene count") +xlab("Transcriptome Mapped (Log2 + 1)")  )
	# 2d density plot
	# print(ggplot(data, aes_string(x="Intergenic.Rate", y="X.trans_mapped")) + stat_density2d(aes(fill = ..density..), geom="tile", contour=FALSE) + scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(9,"Spectral"))(50))))


	# for columns in quality_metric_cols, sort the samples by their values in that
	# column c, and then plot a scatter of their values

	quality_metric_cols = c("trans_mapped", "gene_count")
	# print(data$Sample)

	for (i in 1:length(quality_metric_cols)){
		c = quality_metric_cols[[i]]
		cat(sprintf("Constructing scatter and cumulative hist for %s .. \n", c))
		sorted = data[with(data, order(data[c])), ]
		print(head(sorted,n=2))
		sorted$Sample = factor(sorted$Sample, levels=sorted$Sample, ordered=TRUE)
		print(sorted$Sample)
		scatter = ggplot(sorted, aes_string(x="Sample", y=c)) + geom_point() + theme_bw() + ggtitle(paste(c, " for all samples", sep="")) + theme(text = element_text(size=6), axis.text.x = element_text(angle=90, vjust=1)) 
		print(scatter)


		cumulative_hist = ggplot(sorted, aes_string(x=c)) + geom_histogram(aes(y=cumsum(..count..)), binwidth=10)+ stat_bin(aes(y=cumsum(..count..)),geom="line",color="green") + theme_bw()
		print(cumulative_hist)
	}

	dev.off()
}


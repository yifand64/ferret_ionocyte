

#library(Heatplus)
#require(gplots)
#library(clusterSim)
library(fastcluster)
#require(devtools)
#require(plyr)
library(parallel)
library(cluster)
library(vegan)
library(gridExtra)
library(fpc)
homedir = path.expand("~")

possible.markers <- function(x, m=6, n=5)
{
	return (x[rowSums(x > m) > n,])
}


heatmap <- function(x, disp.min=-2.5, disp.max=2.5, scale=T, clip=T)
{
	if(scale)
	{
		x=t(scale(t(x)))  
	}
	if(clip)
	{
		x =minmax(x,min=disp.min,max=disp.max)
	}
	aheatmap(x)	
	return(x)
}

test.sparse.hclust <- function(x=NULL)
{
	if(is.null(x))
	{
		x = iris[, 1:4]
	}
	shc = sparse.hclust(x)
	# clusters.sparse <- cutreeDynamic(shc$hc, 
	# 			distM = shc$dists, 
	# 				deepSplit=2,
	# 				minClusterSize=5) 
	clusters.sparse = cutree(shc$hc, h=0.4)
	
	
	print(clusters.sparse)
	iris["sparse.hclust"] = clusters.sparse
	ggplot(iris, aes(x=Sepal.Length, y=Sepal.Width, fill=sparse.hclust)) + geom_point()
}

sparse.hclust <- function(x)
{
	library(sparcl)
	perm.out <- HierarchicalSparseCluster.permute(as.matrix(x), wbounds=c(1.5,2:6), nperms=5)
	sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists, wbound=perm.out$bestw, method="complete")
	plot(sparsehc)
	return(sparsehc)
}

drop.small = function(counts, 
					clusters, 
					seurat.obj = NULL,
					dropped.samples = -1, # samples in small clusters are marked with this symbol, and will be excluded from downstream anaylses
					min.cluster.size = 3)
{
	cat(sprintf("\n Dropping clusters with less than %i members .. \n", min.cluster.size))
	freq = table(clusters)
	cluster.names = unique(clusters)
	total.clusters = length(freq)
	for(i in 1: length(cluster.names))
	{
		#get the size of the ith cluster
		this.cluster = cluster.names[i]
		cat(sprintf("Checking \'%s\'..", this.cluster))

		cluster.indices = which(clusters==this.cluster)
		if(length(cluster.indices) < min.cluster.size)
		{
			cat(sprintf("\n 	WARN: Dropping cluster --> %s! [%i samples]\n", this.cluster, length(cluster.indices)))
			not.in.cluster = which(clusters != this.cluster)
			clusters[cluster.indices] = dropped.samples
			# print(colnames(counts)[cluster.indices])
			# cat(sprintf("Samples in cluster: %i, not: %i, total: %i \n", length(cluster.indices), length(not.in.cluster), ncol(counts)))
			# clusters = clusters[not.in.cluster]
			# counts = counts[, not.in.cluster]

			# print("Clusters:")
			# print(clusters)
		}else{
			cat(sprintf(" OK [%i cells]\n", length(cluster.indices)))
		}

	}
	clusters = factor(clusters) # remove 0 levels
	print("Returning clusters:")
	print(table(clusters))
	return (clusters)
}

# general purpose hierarchical clustering and heatmap generation
# function. 
cluster.hierarchical = function(counts,
						pdf.output          	= FALSE,
						cluster_rows 		    = TRUE,
						cluster_columns			= TRUE,
						column.label 		   	= FALSE,
						column.dendro 			= TRUE,
						cluster.get.pvalues 	= FALSE, # use pvClust to generate pvalues
						cluster.get.pvalues.parallel = TRUE, # run pvclust in parallel on maximum available cores
						heatmap.name        	= "heatmap.pdf",
						heatmap.draw        	= FALSE,
						heatmap.title       	= "hclust heatmap",
            			row.labels 			   	= TRUE,
            			row.label.size 			= 0.5,
            			row.dendro				= TRUE,
						z.score.norm.rows		= TRUE,
						z.score.before.cluster  = TRUE,     # if false scaling only helps visualisation and not clustering.
						clip.for.visual 		= TRUE,
						clip 					= 3,
						show.sample.type 		= TRUE, 	# annotates the heatmap showing sample type
						show.cluster.assignment = TRUE, 	# annotates the heatmap showing their cluster assignment
						show.gene.type			= TRUE, 	# if genes are annotated by a 'Type' column, coloured bars will label the rows
						use_dynamic_cut 		= TRUE, 	# cut dendrogram using dynamic tree cutting package
						dendrogram_cut_height 	= 0.8, #0.4
						use_manual_sample_names = FALSE,
            			genes 				    = NULL, #include only these genes
						genes.drop          	= NULL, #remove only these genes
						samples             	= NULL,
						sample.ident       		= NULL, # a named vector / matrix with sample group info (cell-types, batches, mice etc) that will be labeled if show.sample.type is on
						complete_data       	= NULL,
						raw.counts 				= NULL,
						min.cluster.size 		= 10,
						take.log 				= TRUE,
						legend.size          	= 0.75,
						legend.show 			= TRUE,	
						trim.names 				= FALSE,
						verbose   					= FALSE,
						count_string        		= "normalised_count",
						color.palette.heatmap 		= "RdBu",
						color.palette.gene 			= "Set2",
						color.palette.cluster 		= "Accent",
						color.palette.sample 		= c("Set1", "Greys", "Set3", "Paired", "Pastel1"), #if more attributes e.g batch are passed in, these palettes are used
						clustering_metric	 		= "pearson",
						use.manual.labels 			= NULL, # if a set of labels are provided they are used to refer to the clusters.
						clustering_method			= "average")
{

	#either a list of genes or a matrix of genes with annotations.
	# for the moment the annotations have to be in a 'Type' column.
	if(!is.null(genes))
	{
	    if("GENE_SYMBOL" %in% colnames(genes))
	    {
    		genes = subset(genes, !duplicated(genes$GENE_SYMBOL))
    		genes$GENE_SYMBOL = gsub(" ", "", genes$GENE_SYMBOL, fixed = TRUE) #spaces will mess up the match
    	}else{
    		genes = unique(genes)
    	}
	    
	    if("Type" %in% colnames(genes))
	    {
	    	info("Given gene list has type annotations, adding to the heatmap!")
	    	typed = TRUE
	    }else
	    {
	    	if(show.gene.type){
	    		info("Given gene list is not typed!")
	    	}
	    	typed = FALSE
	    }
	    if("GENE_SYMBOL" %in% colnames(genes))
	    {
    		genes.list = as.matrix(genes$GENE_SYMBOL)
		}else{
			genes.list = as.matrix(genes)
		}
	    
	   	if(verbose)
	   	{
	   		cat(sprintf("Using only the given %i genes.. \n", length(genes.list)))
	   	}
	    
		# drop all genes that don't appear in the list from the counts table
		counts = counts[rownames(counts) %in% genes.list, ]
		if(typed){
			# we need to drop any genes that don't actually appear in the frame
			# have to use the mapping from gene -> type defined by the initial
			# genes dataframe to label the genes that were actually found in the
			# counts table.
			remaining.genes = rownames(counts)
			matches = match(remaining.genes, genes$GENE_SYMBOL)
	
			gene_types = factor(genes$Type[matches])

			# keep this, good debug:
			# debug_df = data.frame(gene_types)
			# colnames(debug_df) = "Gene type"
			# debug_df["Gene name"] = factor(genes$GENE_SYMBOL[matches])
			# debug_df["Row names"] = rownames(counts)
			# print(debug_df)
		}
	}
	if(!is.null(genes.drop))
	{
	    cat(sprintf("Dropping the given %i genes.. \n", length(genes.drop)))
	    genes.drop = gsub(" ", "", genes.drop, fixed = TRUE)
	    counts = counts[!rownames(counts) %in% genes.drop, ]
	}
	if(!is.null(samples))
	{
		counts = counts[	, samples]
		sample.ident = sample.ident[samples]
	}
	if(verbose)
	{
		info("Dropping zero rows..")
	}
	counts = counts[which(rowSums(counts) > 0),]
	if(verbose)
	{
		info("Dropping zero columns..")
	}
	counts = counts[, colSums(abs(counts)) != 0]
  
    #print(head(counts))
    
  	info(sprintf("Number of samples --> %s", ncol(counts)))
	info(sprintf("Number of genes   --> %s ", nrow(counts)))
	#print(counts)
	# if(verbose)
	# {
		
	# 	cat(sprintf("Clustering, rows is %s ...\n", if(cluster_rows) "ON" else "OFF"))
	# }
    sample_names = colnames(counts)
	if(trim.names)
	{
		# trim the sample names for readability
		x = strsplit(sample_names, count_string) #  paste(count_string,"_"))
		trimmed_sample_names = list()
		if(column.label){
			for(i in 1:length(x))
			{
				trimmed_sample_names = c(trimmed_sample_names, x[[i]][2])
			}
			trimmed_sample_names = gsub("_", " ", trimmed_sample_names)
		}else{
			trimmed_sample_names = FALSE
		}
	}else
	{
		if(column.label){
			trimmed_sample_names = sample_names
		}else{
			trimmed_sample_names = FALSE
		}
	}
	

	# if(use_manual_sample_names)
	# {
	# 	trimmed_sample_names = c("Paneth", "Enteroendocrine", "Epcam +", "Stem Cells", "Trans-amplified")#c("Paneth [SI]", "Paneth [Org]", "Enteroendocrine [Org]", "Epcam + [Org]", "Stem cells [Org]",  "Trans-amplified [Org]", "Enteroendocrine [SI]", "Epcam + [SI]", "Stem cells [SI]", "Trans-amplified [SI]")
	# }

	# # # find sample classes:
	# # sample_class_title = "Cell type"
	

	if(is.null(sample.ident)){
		
		# info("Assigning samples.. This is deprecated and may fail")
		# # need to eventually deprecate this, and build a sample.ident vector/matrix
		# # by splitting the sample names on '_' and picking a field (like seurat does). 
		# # can then use mapvalues to map the found names to manual names if necessary.
		# sample.ident.labels = c("Paneth","Lgr5Hi","Lgr5lo","Entero", "Endo") #c("Organoid", "In vivo") #"CD103+ CD11b+", "CD103+ CD11b-", "CD103- CD11b+", "CD103- CD11b-", "Macrophage"
		# sample.ident.colors = colorRampPalette(brewer.pal(length(sample.ident.labels), color.palette.sample))(length(sample.ident.labels))
		# sample_class_strings = sample.ident.labels	# c("Org", "SI")  #"_CD103posCD11bpos_DC", "_CD103pos_","_CD11bpos_","_DN", "macro")
		# sample.ident = assign_samples(sample_names, 
  #                                       sample.ident.labels, 
  #                                       sample_class_strings, 
  #                                       "Unknown", 
  #                                       sample_colors=sample.ident.colors)
		sample.ident.labels = rep("1", length(sample_names))
		sample.ident = sample.ident.labels
	}else
	{
		info("Using provided sample ids..")
		sample.ident.labels = levels(sample.ident)
		sample_class_strings = sample.ident.labels
    
	}
	#print("Sample Ident:")
	#print(sample.ident)
	sample.ident= as.matrix(sample.ident)
	sample.color.matrix = matrix(0, nrow(as.matrix(sample.ident)), ncol(sample.ident))
	sample.ident.colors = list()
	#cat(sprintf("Given %i sample identifiers..\n", ncol(as.matrix(sample.ident))))
	for(j in 1:ncol(sample.ident))
	{
		n.types = length(unique(sample.ident[,j]))
		# cat(sprintf("ID %i has %i types..\n", j, n.types))
		colors = colorRampPalette(brewer.pal(n.types, color.palette.sample[j]))(n.types)
		#print(colors)
		sample.ident.colors[[j]] = colors
		# cat(sprintf("Mapping to %i \'%s\' colors..\n", n.types, color.palette.sample[j]))
		# sample.color.matrix[,i] = mapvalues(sample_assignments, from = sample.ident.labels, to = sample.ident.colors)
		sample.color.matrix[,j] = mapvalues(sample.ident[,j], from = unique(sample.ident[,j]), to = sample.ident.colors[[j]])
	}
	# print("Sample IDENT has colnames:")
	# print(colnames(sample.ident))
	colnames(sample.color.matrix) = colnames(sample.ident)
	# apply unique to the columns to a get a list of sample types
	# for each given attribute (batch, cell-type etc.). use split
	# to ensure its a list not a matrix

	# this step is necessary because apply(unique) returns a matrix
	# if all the columns have the same number of unique values, and 
	# list otherwise.
	x = apply(sample.ident, 2, unique)
	if(is.matrix(x))
	{
		sample.ident.labels = split(x, rep(1:ncol(x), each = nrow(x)))
	}else
	{
		sample.ident.labels = x
	}

	names(sample.ident.labels) = colnames(sample.ident)
	names(sample.ident.colors) = colnames(sample.ident)
	print(sample.ident.labels)
	
	# print("Sample coloring..")
	# print(sample.color.matrix)

	#find annotated gene classes. there has to be columns named annotation-x, annotation-y containing 'Yes' or a column named 'Type' containing the annotation name
	if(show.gene.type ){
		if(!is.null(complete_data)){
	      	n_gene_types = length(unique(complete_data[, "Type"]))
	  		gene_type_colors = colorRampPalette(brewer.pal(9, color.palette.gene))(n_gene_types)
	  		gene_assign_obj =  assign_genes(complete_data, gene_type_colors)
	  		gene.colors = gene_assign_obj[["assignments"]]
	  		gene_types = gene_assign_obj[["types"]]
		}else if("Type" %in% colnames(genes))
     	{
			#gene_types = genes[, "Type"]
			n_gene_types = length(unique(gene_types))  
			cat(sprintf("Found %i gene types.. \n", n_gene_types))
			gene_type_colors = colorRampPalette(brewer.pal(9, color.palette.gene))(n_gene_types)
			if(verbose)
			{
				info("Mapping gene types to colours..")
			}
			gene.colors = mapvalues(gene_types, from=unique(gene_types), to=gene_type_colors)
			gene.colors = as.matrix(as.character(gene.colors))
			rownames(gene.colors) = gene_types
			colnames(gene.colors) = "Gene Type"
			gene.colors = t(gene.colors)
          
	      }else{
	        print("Show gene type is on, but no gene types provided")
	        gene.colors = FALSE
	        show.gene.type   = FALSE
	        gene_types = NULL
	      }
	}else{
		if(verbose)
		{
			info("Show gene type is off. Not labeling categories of genes..")
		}
    	gene.colors 	= NULL
		show.gene.type  = FALSE
    	gene_types 		= NULL
	}
	# if(verbose)
	# {
	# 	print("Gene types:")
	# 	print(gene_types)
	# 	print("Gene assignments:")
	# 	print(gene.colors)
	# }
	# print("Gene types:")
	# print(t(gene_types))
	# print("Gene assignments:")
	# print(t(gene.colors))

	if(verbose)
	{
		cat(sprintf("\"# genes: %i, # samples: %i \"\n", nrow(counts), length(sample_names)))
	}

	# hclust clustering function
	run_clustering <- function(x, 
								clustering_metric = "pearson",
								cluster_method = "average") 
	{
		if(verbose)
		{
			cat(sprintf("Running clustering..\n"))
		}
		return (hclust(as.dist(1-cor(t(as.matrix(x)), method = clustering_metric)), method = cluster_method))
	}


	if(take.log){
		log.counts = log2(counts+1)
	}else
	{
		log.counts = counts #an earlier step in the pipeline must have done this already.
	}

	if(z.score.before.cluster )
	{
		info("Z-score normalising.. ")
		log.counts  <- t(scale(t(log.counts)))         # zscore -normalise
	}

	# if(center.data)
	# {
	# 	log.counts = scale(log.counts, scale = FALSE)
	# }
	
	# print(head(log.counts))

	if(cluster_columns){
		info("Clustering columns..")
		cl.col <- run_clustering(t(log.counts)) 
		column_dendrogram = as.dendrogram(cl.col)
	} else {
		column_dendrogram = FALSE
	}

	# cluster rows:
	if (cluster_rows){ 
		info("Clustering rows..")
		cl.row <- run_clustering(log.counts)
		if(use_dynamic_cut)
		{
			info("Dynamic tree cutting row dendrogram.. ")
			library(dynamicTreeCut)
			clusters.rows <- cutreeDynamic(cl.row, 
					distM = as.matrix(1-cor(as.matrix(t(log.counts)))), 
					deepSplit=2,
					minClusterSize=min.cluster.size) 
		}else
		{
			color_h_cut_r = max(cl.row$height)*dendrogram_cut_height
			clusters.rows <- cutree(cl.row, h=color_h_cut_r)
		}
		
		row_dendrogram  = as.dendrogram(cl.row)	
		
	}else{
		row_dendrogram <- FALSE
		
	}

	# sample and cluster colouring:
	# cluster output
	if(cluster_columns){
		color_h_cut_c = max(cl.col$height)*dendrogram_cut_height
		
		if(use_dynamic_cut){
			info("Dynamic tree cutting column dendrogram.. ")
			library(dynamicTreeCut)
			clusters.column <- cutreeDynamic(cl.col, 
					distM = as.matrix(1-cor(as.matrix(log.counts))), 
					deepSplit=2,
					minClusterSize=min.cluster.size) 
		}else{
			# cluster splitting (naive)!
			spl.r = get.splits(as.dendrogram(row.dendro), 0.36)
			clusters.column <- cutree(cl.col, h=color_h_cut_c) # k=n_colored_clusters_c
		}
		
		if(verbose)
		{
			print("Found clusters:")
    		print(table(clusters.column))
		}

		cluster_colors = colorRampPalette(brewer.pal(max(clusters.column), color.palette.cluster))(max(clusters.column))
		cluster_names = rep("Cluster ", max(clusters.column))
		if(!is.null(use.manual.labels))
		{
			cluster_names[1:length(use.manual.labels)] = use.manual.labels
			if(verbose)
			{
				print("Labeled: ")
				print(cluster_names)
			}
			n_clusters = length(table(clusters.column))
			if(length(use.manual.labels) < max(clusters.column))
			{
				cat(sprintf("WARN: %i clusters found and only %i manual names provided..\n", 
					n_clusters, length(use.manual.labels)))
				unlabeled = c(rep("", length(use.manual.labels)), 
					as.character(seq(length(use.manual.labels)+1, n_clusters)))
				cluster_names = paste(cluster_names, unlabeled, sep="")
				if(verbose)
				{
					print("Cluster names:")
					print(cluster_names)
				}
			}
		}else{
			# for(i in 1: max(clusters.column))
			# {
			# 	cluster_names[i] = paste("Cluster", i, sep=" ")
			# }
			cluster_names = unique(clusters.column)
		}

		# switch from integer labels to more meaningful names
		# (if provided, otherwise its just Cluster 1, etc.) 
		# int_labels = seq(1,max(clusters.column))
		# clusters.column = unlist(mapvalues(clusters.column, int_labels, cluster_names))
		# names(clusters.column) = sample_names

		for(i in 1: length(unique(clusters.column)))
		{
  			cluster.indices = which(clusters.column==cluster_names[i])
  			cluster.members = colnames(counts)[cluster.indices]
  			if(verbose)
  			{
  				cat(sprintf("%s: \n ", cluster_names[i]))
  				# print("Indices:")
  				# print(this.cluster.indices)
  				print("Names:")
  				print(cluster.members)
  			}
		}

		# label cells by their cluster assignment:
		if(heatmap.draw)
		{
			if(show.cluster.assignment & show.sample.type){
				cat(sprintf("Labeling cluster assignments.. \n"))
	      		cluster.color.labels = rep(0, length(clusters.column))
				cluster.color.labels = as.character(mapvalues(clusters.column, cluster_names, cluster_colors))
				if(verbose)
				{
					cat(sprintf("Building color matrix.. \n"))
				}
				# print(as.matrix(sample.color.matrix))
				column.color.matrix = cbind(as.matrix(sample.color.matrix), cluster.color.labels)
				colnames(column.color.matrix) = c(colnames(sample.color.matrix), "Cluster assignment")
			}else if(show.sample.type){
				column.color.matrix = as.matrix(sample.color.matrix)
				colnames(column.color.matrix) = colnames(sample.color.matrix)
			}else
			{
				column.color.matrix = NULL
			}
		}
		# ncolumn.color.matrixrows = nrow(column.color.matrix)
		# cat(sprintf("Color labeling matrix has %i rows. \n", ncolumn.color.matrixrows))
	}

	  # print("color matrix:")
	  # print(column.color.matrix)

	if(cluster_columns)
	{
		if(cluster_rows && row.dendro)
		{
			dendro = "both"
		}else{
			if(column.dendro)
			{
				dendro = "col"
			}else{
				dendro = "none"
			}
			
		}
	}else{
		if(cluster_rows && row.dendro)
		{
			dendro="row"
		}else
		{
			dendro="none"
		}
	}
	if(verbose)
	{
		cat(sprintf("Set dendrogram to %s ..\n", dendro))
	}
	
	clipped_data = log.counts

	if(z.score.norm.rows){
		# Scaling data (if not done before clustering)
		if(!z.score.before.cluster)
		{	
			info("Z-score normalising.. ")
			clipped_data <- t(scale(t(clipped_data)))         # zscore -normalise
		}
		key_name = "Z-Score"
	}else{
		key_name = "Log2 counts"
	}

	# enhances visual contrast:
	if(clip.for.visual)
	{
		cat(sprintf("Clipping z scores to +- %f ..\n", clip ))
		clipped_data <- pmin(pmax(clipped_data, -clip), clip)      ## Compressing data to max/min +/- clip
	}

	if(row.labels)
	{
		names = row.names(clipped_data)
	}else{
		names = ""
	}

	# key location
	# bottom = 2.5
	# left = 2.5
	# top = 1
	# right = 0

	# # key shape
	# key_a = 1.5
	# key_b = 0.5
	# key_c = 0
	if(verbose){
		print("Gene color matrix size:")
		print(dim(gene.colors))
		print("Data matrix size:")
		print(dim(clipped_data))
	}
	

	#---- Heatmap generation ---------#

	if(ncol(clipped_data) > 200 & heatmap.draw)
	{
		warn("Over 200 samples.. heatmap might not render properly..")
		plot.new()
		par(mar=c(7,4,4,2)+0.1) 
	}

	if(heatmap.draw)
	{
		if(pdf.output){
		
		if(ncol(clipped_data) > 200)
		{
			info("Rendering wide heatmap to PDF..")
			#png(file="large_heatmap.png", width=1600, height=1600)
			pdf(file=heatmap.name,  width=16, height=8.5)
		}else{
			info("Will render heatmap to PDF..")
			pdf(file=heatmap.name,  width=8.5, height=11)
		}
		}else
		{
			info("Will render heatmap to SCREEN..")
		}
	}
	

  if(heatmap.draw){

  	if(color.palette.heatmap != "SamR")
	{
	 	if(length(color.palette.heatmap) > 1)
	 	{
	 		info(sprintf("Building heatmap colors from gradient using %i colors: ", length(color.palette.heatmap)))
	 		info(sprintf(color.palette.heatmap))
	 		hm_cols = colorRampPalette(color.palette.heatmap)(512)
 		}else{
 			info(sprintf("Using brewer palette: %s", color.palette.heatmap))
 			hm_cols = rev(colorRampPalette(brewer.pal(11, color.palette.heatmap))(512))
 		}
	 	
	}else{
		info("Using SamR's heatmap colors..")
		hm_cols = get.hmap.col()
	}

  	hm = tryCatch(
  	{
      info("Building heatmap..")
      
    	hm <- heatmap.3(
      		clipped_data, 
      		main = heatmap.title, #"Subpopulations in single-cell data \n using markers defined by clusters", "SI and Organoid \n epithelial cell populations [HiSeq]",#
      		cex.main=0.75,
      		#col = colorRampPalette(c("green", "black", "red"))(n=256),
      		col  = hm_cols,
      		
      		# Dendrograms
      		Rowv = row_dendrogram,
      		Colv = column_dendrogram,
      		dendrogram = dendro,
      		
      		# #Color labels:
      		#RowSideColors = gene.colors,
      		ColSideColors = column.color.matrix,
      		#useRaster=TRUE,
      		margins		= c(6,16),#c(6,12),
      		#srtCol		= 90, #slant the column labels
      		trace		= "none", 
      		density.info= "none",#"colnhistogram", 
      		cexCol		= 1.5, 	# col label size
      		cexRow		= row.label.size,	# row label size
      		#labRow="",
      		labCol		= trimmed_sample_names,
      		labRow		= names,
      		scale 		= "none", #"row",
      		keysize 	= 1, 
      		KeyValueName=key_name,
    	)
    	# add legend for cluster / sample type / gene type labeling
    	if(legend.show){
			cat(sprintf("Legend settings: \n"))
			cat(sprintf("samples -- %s \n", show.sample.type))
			cat(sprintf("clusters -- %s \n", show.cluster.assignment))
			cat(sprintf("genes -- %s \n", show.gene.type))
			
			if(show.gene.type)
			{
				print(gene_types)
				gene.type.names = as.character(unique(gene_types))
				names(gene.type.names) = paste("GeneType",seq(1:length(gene.type.names)), sep="")
				print(unique(gene.type.names))
			}
			
			
			legend.names = NULL
	    	legend.colors = NULL

	    	if(show.cluster.assignment && show.sample.type && show.gene.type ){
	    		if(is.null(gene_types)){
	    			warn("Cannot show gene types, they are null..")
	    			sample.ident.labels$Cluster = cluster_names
		    	  	sample.ident.colors$Cluster = cluster_colors
		    	  	legend.names = unlist(sample.ident.labels)
		    	  	legend.colors = unlist(sample.ident.colors)
	    		}else{
	    			info("Drawing sample, cluster and gene type legend..")
	    			sample.ident.labels$Cluster = cluster_names
		    	  	sample.ident.colors$Cluster = cluster_colors

		    	  	legend.names = c(unlist(sample.ident.labels), gene.type.names)
		    	  	legend.colors = c(unlist(sample.ident.colors), gene_type_colors)	    
	    		}
	    	}else if (show.sample.type & show.cluster.assignment){
	    	  	info("Drawing cluster assignment and sample_type legend..")
	    	  	sample.ident.labels$Cluster = cluster_names
	    	  	sample.ident.colors$Cluster = cluster_colors
	    	  	legend.names = unlist(sample.ident.labels)
	    	  	legend.colors = unlist(sample.ident.colors)
	    	  	legend("topright",legend=legend.names, fill=legend.colors, border=FALSE, bty="n", y.intersp = 0.9* legend.size, cex=0.7* legend.size)
	    	  
	    	}else if (show.sample.type & show.gene.type){
	    		info("Drawing sample-type and gene-type legend..")
	    		print(as.character(unique(gene_types)))
	    		print(unlist(sample.ident.labels))
	    		legend.names = c(unlist(sample.ident.labels), gene.type.names)
	    		# print("Legend names")
	    		# print(legend.names)
	    	  	legend.colors = c(unlist(sample.ident.colors), gene_type_colors)	    
	    	}else if (show.cluster.assignment){
	    		info("Drawing cluster assignment legend..")
	    		legend("topright",legend=c(cluster_names), fill=c(cluster_colors), border=FALSE, bty="n", y.intersp = 0.8, cex=0.7)
	    	}else if (show.sample.type)
	    	{
	    		info("Drawing sample-type legend..")
	    	  	legend("topright",legend=c(sample.ident.labels), fill=c(sample.ident.colors), border=FALSE, bty="n", y.intersp = 0.8, cex=0.7)
	    
	    	}
	    	if(is.null(legend.names) | is.null(legend.colors))
	    	{
	    		warn("Show legend is on but nothing to show!")
    		}else{
    			# print("Legend names:")
    			# print(legend.names)

    			print(as.character(unique(gene_types)))
    			legend("topright",legend=legend.names, fill=legend.colors, border=FALSE, bty="n", y.intersp = 0.9* legend.size, cex=0.7* legend.size)
    		}

    	}else{
    		info("No legend")
    	}
      	
  },
  error=function(err){
	  	error("Heatmap generation failed!")
	  	cat(paste("    ", err))

   },
	warning=function(err) {
	    warn(err)
	},
	finally={
       
       if(pdf.output)
	  {
	      dev.off()
	  }
	})

}else{
	info("Not drawing heatmap..")
  
}
return(list("hc.row"=clusters.rows, "hc.col"=clusters.column))
}



# given a cluster boot object, draw a simple plot
# showing which clusters are stable.
plot.cluster.boot <- function(cboot_result, verbose=F)
{
	cluster.types = c("Dissolved", "Untrustworthy", "Meta-stable", "Stable", "Highly stable")
	colours = c("darkgrey", "darkorange4", "darkgoldenrod3",  "darkolivegreen3", "darkgreen")
	cutoffs = c(0, 0.5, 0.6, 0.75, 0.85, 1)
	x = cboot_result
	binned.stability = cut(x$bootmean, breaks=cutoffs)
	stability = x$bootmean
	cluster.size = table(x$partition)
	name = unique(x$partition)
	r = data.frame(cbind(stability, name, binned.stability, cluster.size))
	
	vals.in.use =sort(unique(binned.stability)) #only colour by the bins that are actually present
	g=ggplot(r, aes(x=name, y=stability)) + geom_bar(stat="identity", aes(fill = as.factor(binned.stability))) + 
		theme_bw() + scale_fill_manual("Cluster Stability", values=colours[vals.in.use], labels=cluster.types[vals.in.use]) + xlab("Cluster Index") + 
		ylab(sprintf("Cluster Stability [Bootstraps = %s]", x$B)) + scale_y_continuous(breaks=seq(0, 1, by=0.1))

	
	if(verbose)
	{
		r$binned.stability = NULL
		cat(" Cluster Bootstrap Info:\n")
		cat("=========================\n")
		print(r)
	}
	return(g)
}

f <- function(k, cluster.func=NULL, data=NULL, distM=NULL){
 	cb=clusterboot(data, clustermethod = cluster.func, distances=distM, k=k)
 	fdest = sprintf("cluster.stability.k_%s.txt", k)
 	stability = cb$bootmean
 	cluster = unique(cb$partition)
 	d = as.data.frame(cbind(cluster, stability))
 	d = d[with(d, order(cluster)), ]
 	write.table(d, file=fdest, quote=F, row.names=F)
 	#rval = list("plot"= plot.cluster.boot(cb), "cboot_result"=cb)
 	rval = plot.cluster.boot(cb)
 	return (rval)
}

# vary k over several values and plot the stability of clusters obtained (calculated
#	by bootstrapping)
cluster.stability.bootstrap <- function(data, cluster.func, distM=NULL, draw.plot=T, k.range=2:20, n.cores=1)
{
	
	if(is.null(distM))
	{
		distM = as.matrix(Dist(data, method="euclidean", nbproc=detectCores()))
	}
	plots = mclapply(k.range, f, cluster.func=cluster.func, data=data, distM=distM, mc.preschedule=FALSE, mc.set.seed=FALSE, mc.cores=n.cores)
	i = 1
	for(p in plots)
	{
		print(p)
		dest = sprintf("cluster.stability.k_%s.pdf", k.range[i])
		ggsave(p, file=dest, width=10, height=8)
		i = i + 1
	}
}

# run a (VERY SLOW) boot strapping over k and calculate all available
# cluster separation indices. 
cluster.how.many.bootstrap <- function(data, cluster.func, 
	distance.method="euclidean",k.range=2:10, B=20, draw.plot=T, parallel=F, pdf.output="cluster.metrics.pdf")
{

	n <- nrow(data)
	stat <- c("avg.silwidth", "g2", "g3", "pearsongamma", "dunn", "ch",
	  "within.cluster.ss", "wb.ratio")
	res <- array(dim=c(B, length(k.range), length(stat)),
	         dimnames=list(1:B, k.range, stat))
	for (b in 1:B) {
		cat("Iteration ", b, "\n")
		tmp <- data[sample(n, n, replace=TRUE), ]
		dd = as.matrix(Dist(tmp, method=distance.method, nbproc=detectCores()))
		if(parallel)
		{
			warn("not done yet.")
			return (FALSE)
		}else{
			for (i in k.range)
			{
				cat("	K=", i, "\n")
				if(as.character(substitute(cluster.func))=="graph.cluster")
				{
					cl = cluster.func(tmp, dm=dd, k=i)
				}else
				{
					cl = cluster.func(dd, k=i)
				}
				
				#print(table(cl$partition))
				res[b,i-1,] <- as.numeric(unlist(cluster.stats(dd, 
					cl$partition, G2=TRUE, G3=TRUE)[stat]))
			
			}
		}
		
	}
	df = as.data.frame.table(res) #merge 3d array into data.frame
	colnames(df) = c("bootstrap.Iteration", "k", "variable", "value")
	if(draw.plot)
	{
		info("Plotting..")
		# g = ggplot(df, aes(x=k, y=value)) + geom_violin() + theme_bw() + 
		# 	facet_wrap(~variable, scale="free")
		pdf(pdf.output, width=16, height=12)
		g = ggplot(df, aes(x=k, y=value)) + geom_boxplot() + theme_bw() + 
			facet_wrap(~variable, scale="free") + ylab("") 
		print(g)
		dev.off()
		#ggsave(g, file=pdf.output)
	}
	return(df)
}



#http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
# check/analyse the number of clusters in a dataset
cluster.how.many <- function(d, 
						run.affinity=F, 
						run.kmeans=T, 
						run.hclust=T,
						run.calinski=T, 
						run.mclust=F, # so slow!
						run.pamk=T,
						run.bclust=F,
						pdf.output=F,
						run.gap=T,
						pdf.name="clusters_how_many.pdf")
{
	if(pdf.output)
	{
		pdf(pdf.name)
	}


	
	cat("Running k-means clustering..\n")
	wss <- (nrow(d)-1)*sum(apply(d,2,var))
	  for (i in 2:15) wss[i] <- sum(kmeans(d,
	                                       centers=i)$withinss)
	plot(1:15, wss, type="b", log="y", xlab="Number of Clusters",
	     ylab="Within groups sum of squares")
	nclusts = data.frame(list("Method"="K-Means", "N.Clusters"="-"))
	


	if(run.pamk)
	{
		library(fpc)
		cat("Running PAM..\n")
		pamk.best <- pamk(d)
		cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
		
		#plot(pam(d, pamk.best$nc))
		
		cluster.vector = pamk.best$pamobject$clustering
		clusters = data.frame(cluster.vector)

		# run k-means with the optimal k (determined by PAM)
		cat("Running k-means with above optimal k .. \n")
		kmeans.clusters = kmeans(d, centers=pamk.best$nc)$cluster
		clusters = cbind(clusters, kmeans.clusters)
		colnames(clusters) = c("PAM.K", "K.Means")
		nclusts = rbind(nclusts,data.frame(list("Method"="Average silhouette", "N.Clusters"=as.character(pamk.best$nc))))
	}

	if(run.affinity)
	{
		library(apcluster)
		cat("Running affinity propagation..\n")
		d.apclus <- apcluster(negDistMat(r=2), d)
		cat("Affinity propagation optimal number of clusters:", length(d.apclus@clusters), "\n")
		nclusts = rbind(nclusts,data.frame(list("Method"="Affinity propagation", "N.Clusters"= as.character(length(d.apclus@clusters)))))
		#plot(d.apclus, d)
		clusters["Affinity"] = d.apclus@idx
	}

	if(run.mclust)
	{
		library(mclust)
		# Run the function to see how many clusters
		# it finds to be optimal, set it to search for
		# at least 1 model and up 20.
		cat("Running model based clustering (Mclust)..\n")
		d_clust <- Mclust(as.matrix(d), G=1:20)
		m.best <- dim(d_clust$z)[2]
		clusters["MClust"] = d_clust$classification
		cat("Model-based optimal number of clusters:", m.best, "\n")
		nclusts = rbind(nclusts,data.frame(list("Method"="Mclust", "N.Clusters"= as.character(m.best))))
		plot(d_clust, what = c("BIC"))
	}

	if(run.gap)
	{
		library(cluster)
		cat("Running gap statistic.. \n")
		cg = clusGap(d, kmeans, 10, B = 100, verbose = interactive())
		n.gap = with(cg, maxSE(Tab[,"gap"],Tab[,"SE.sim"]))
		plot(cg, main = "Gap statistic")
		nclusts = rbind(nclusts,data.frame(list("Method"="Gap-Statistic", "N.Clusters"= as.character(n.gap))))
	}

	if(run.bclust)
	{
		# a Bayesian clustering method, good for high-dimension data, more details:
		# http://vahid.probstat.ca/paper/2012-bclust.pdf
		cat("Running bayesian clustering..\n")
		library(bclust)
		x <- as.matrix(d)
		d.bclus <- bclust(x, transformed.par = c(0, -50, log(16), 0, 0, 0))
		viplot(imp(d.bclus)$var); 
		c = d.bclus$optim.alloc
		clusters["Bayesian"] = c
		plot(d.bclus); 
		ditplot(d.bclus)
		dptplot(d.bclus, scale = 20, horizbar.plot = TRUE,varimp = imp(d.bclus)$var, horizbar.distance = 0, dendrogram.lwd = 2)
		nclusts = rbind(nclusts, data.frame(list("Method"="Bayesian", "N.Clusters"= as.character(d.bclus$optim.clustno))))
		# I just include the dendrogram here
	}

	if(run.calinski)
	{
		require(vegan)
		fit <- cascadeKM(scale(d, center = TRUE,  scale = TRUE), 1, 10, iter = 1000)
		plot(fit, sortg = TRUE, grpmts.plot = TRUE)
		calinski.best <- as.numeric(which.max(fit$results[2,]))
		cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")
		nclusts = rbind(nclusts,data.frame(list("Method"="Calinski", "N.Clusters"= as.character(calinski.best))))
		# 5 clusters!
	}

	if(run.hclust)
	{
		hclusts = cluster.hierarchical(t(d))
		n_hclusts = length(unique(hclusts))
		clusters["hclust"] = hclusts
		cat("Hclust dynamic tree cut number of clusters:", n_hclusts, "\n")
		nclusts = rbind(nclusts, data.frame(list("Method"="hclust", "N.Clusters"= as.character(n_hclusts))))
	}

	# this dataframe has the optimal number of clusters, determined by each method.
	plot.new()
	grid.table(nclusts)


	if(pdf.output)
	{
		dev.off()
	}
	#eval = clusters.evaluate(d, clusters, verbose=TRUE)
	return (list("clusters"=clusters, "eval"=eval))
}

clusters.rename <- function(x, factor=FALSE)
{
	vals = unique(x)
	
	to = factor(seq(1:length(vals)))
	clusters = mapvalues(x, vals, to)
	if(!factor)
	{
		clusters = as.numeric(clusters)
	}
	return (clusters)
}

# For either a single clustering or multiple clusterings of 
# a given dataset, calculate three metrics to measure the 
# quality of the clustering. See https://cran.r-project.org/web/packages/clValid/vignettes/clValid.pdf
# for more details 

# The best clustering scheme essentially minimizes the Davies–Bouldin index
# nThe
#connectivity has a value between zero and ∞ and should be minimized.
# The Silhouette Width thus lies in the interval
# [−1, 1], and should be maximized.
clusters.evaluate <- function(data, clusters, verbose=F, distM=NULL)
{
	library(clValid)
	cat(sprintf("Evaluating clusterings produced by %i methods.. \n", ncol(clusters)))
	cat(sprintf(colnames(clusters) ))	
	cat(sprintf("\n \n"))

	if(is.null(distM))
	{
		if(verbose)
		{
			cat("Calculating distance matrix..\n")
		}
		distM = dist(data)
	}
	r = sapply(clusters, function(x){
		# rename clusters as integers
		cat("\n \n ")
		clusters = clusters.rename(x)
		
		if(verbose)
		{
			cat("Calculating connectivity..\n")
		}
		conn = connectivity(distance=distM, clusters)
		if(verbose)
		{
			cat(sprintf("%f\n", conn))
		}
		
		# gap = index.Gap()
		print(table(clusters))
		if(verbose)
		{
			cat("Calculating silhouette width..\n")
		}
		sil = index.S(distM,clusters) #silhouette(x, dist(data))
		if(verbose)
		{
			cat(sprintf("%f\n", mean(sil)))
		}

		if(verbose)
		{
			cat("Calculating Davies-Bouldin index..\n")
		}
		d=index.DB(data, clusters)
		db=d$DB
		if(verbose)
		{
			cat(sprintf("%f\n", db))
		}
		# dun = dunn(distance=distM, clusters)
		# if(verbose)
		# {
		# 	cat("Calculating Dunn index..\n")
		# 	cat(sprintf("%f\n", dunn))
		# }

		n = length(unique(clusters))
		rval = c(conn, mean(sil), db, n)	
		names(rval) = c("Connectivity - less is better", "Silhouette -- more is better", "Davies-Bouldin -- less is better", "N_Clusters")
		cat("Done.\n")
		return (rval)
		
	})
	r = data.frame(r)
	r = data.frame(t(r))
	r["Method"] = rownames(r)
	print(head(r))

	pca = prcomp(t(data))
	
	print("Size of clusters: ")
	print(dim(clusters))
	print("PCA:")
	print(dim(pca$rotation))
	cat("Adding PCs for visualising..\n")
	print(colnames(clusters))
	clusters = clusters
	#clusters.pca = melt(cbind(clusters, pca$rotation[, 1:2]), id=c("PC1", "PC2"))#rownames(r))
	clusters.pca = cbind(clusters, pca$rotation[, 1:2])
	#print(head(clusters.pca))
	
	print(clusters.pca)
	pdf("cluster_vis_pca.pdf")
	for(i in 1:ncol(clusters))
	{
		method = colnames(clusters)[i]
		cat(sprintf("PCA plot for %s\n" ,method))
		# print(clusters.pca[method])
		# print(as.character(clusters.pca[method]))
		# print(factor(clusters.pca[method]))
		clusters.pca[method] = factor(unlist(clusters.pca[method]))
		#print(clusters[method])
		g.pca=ggplot(clusters.pca, aes_string(x="PC1", y="PC2", colour=method)) + geom_point()  + ggtitle(paste(method, "clustering visualised via PCA")) + theme_bw()  #+ theme(plot.title = element_text(size=6))
		print(g.pca)
	}
	dev.off()
	

	


	cat("Melting eval data for plotting..\n")
	

	df = melt(r, id.vars = c("N_Clusters", "Method"))
	
	info.string = c("connectivity has a value between zero and ∞ and should be minimized.",
	 "Silhouette Width thus lies in the interval [−1, 1], and should be maximized")
	df["Info"] = mapvalues(df$variable, c("Connectivity", "Silhouette"), info.string)
	title.string = paste(info.string, sep="\n")
	print(title.string)
	#print(df)
	

	g=ggplot(df, aes(x=Method, y=value, size=N_Clusters)) + geom_point()  + facet_wrap(~variable, scale="free") + ggtitle(title.string) + theme(plot.title = element_text(size=6))

	ggsave(g, filename="cluster_evaluation.pdf")

	g=ggplot(df, aes(x=Method, y=value, size=N_Clusters)) + geom_point()  + facet_wrap(~variable, scale="free") + ggtitle(title.string) + theme(plot.title = element_text(size=6))

	cluster.stats(distM, clusters[ ,1], clusters[, 2])

	return (r)
}

# run hierarchical clustering with bootstrap to estimate p-values
# for clusters. note, will take a pretty long time for large datasets.
cluster.pvclust <- function(d, 
							parallel=T, 
							parallel.mpi=F, # mpi is unreliable on the LSF compute grid.
							nboot=1000, 
							pdf.output=F,
							pdf.name="pv_clust.pdf",
							min.p = 0.05, 
							n.cores=detectCores())
{
	library(pvclust)
	if(parallel)
	{
		if(parallel.mpi)
		{
			library(snow)
			cluster.type="MPI"
		}else
		{
			library(parallel)
			cluster.type="PSOCK"
		}
		
		cat(sprintf("Running parallel pvClust [%i cores]..\n", n.cores))
		cl <- makeCluster(n.cores, type=cluster.type)
		pv.clusters <- parPvclust(cl, d, nboot=nboot)
		stopCluster(cl)
	}else
	{
		pv.clusters = pvclust(d, nboot=nboot)
	}
	
	if(pdf.output)
	{
		pdf(pdf.name)
	}
	cat("Plotting..\n")
	plot(pv.clusters, cex.pv=0.25, cex=0.25)
	cat(sprintf("Drawing rectangles around clusters significant [p < %f]", min.p))
	pvrect(pv.clusters, alpha=1-min.p)
	
	if(pdf.output)
	{
		dev.off()
	}

	# pvpick(pv.clusters, alpha = 1-min.p, type = "geq")
	return(pv.clusters)
}
	

check.cluster.names <- function(clusters)
{
	for(i in 1:length(unique(clusters)))
    {
        cluster = unique(clusters)[[i]]
        if(suppressWarnings(!is.na(as.numeric(as.character(cluster)))))
        {
            warn(sprintf("Cluster %s is numeric. Adding 'Cluster' to the label..", cluster))
            old.name = cluster
            new.name = paste("Cluster", cluster)
            clusters = mapvalues(clusters, old.name, new.name)
        }
    }  
    return (clusters)
}




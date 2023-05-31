# Given a subset of a counts table, find the set of n most differentially
# expressed genes by comparing with a random sample of the rest of the 
# table. Can also show just the most expressed genes in that population.

#Population_info_class:
PopInfo <- setClass("PopInfo", slots = c(tested="DGEExact", 
 											#dge_list="DGEList", 
 											counts="data.frame"))



examine_subpopulation <- function(counts_table, 
									complete_table			= NULL,
									raw.counts 				= NULL,
									n_background_samples	= 0,
									sub.sample.groups.to    = 0,
									poi_cols 				= NULL, 
									poi_string				= NULL, 
									enriched_genes_only 	= TRUE,
									n_diff_genes			= 50, 
									min_p_val				= 0.05,
									name_of_population		= "POI",
									show_top				= TRUE, 
									single_cell_data 		= FALSE,
									seurat.obj 				= NULL,
									batch 					= NULL,
									verbose 				= TRUE,
									annotate.marker.genes   = F,
									expression_unit_string	= "normalised_count")
{
	

	if(is.null(poi_cols) & is.null(poi_string)){
		stop("ERROR: Must provide either a string to search for or a column of samples to examine!")
	}

	if(is.null(complete_table) & is.null(raw.counts) & single_cell_data){
		stop("Cannot run SCDE differential expression without raw counts values. Must provide either complete data frame or raw counts arguments!")
	}

	if(!is.null(poi_string))
	{
		cat(sprintf("Searching for sample names containing %s .. \n", poi_string))
		poi_cols = find_samples(counts_table, poi_string)
	}

	counts_table = data.frame(counts_table, check.names=F)
	if(!is.null(complete_table))
	{
		complete_table = data.frame(complete_table, check.names=F)
	}
	
	cat(sprintf("Examining %i samples.. \n", length(poi_cols)))
	#cat(sprintf("	%s \n", poi_cols))

	# sample 'size_of_random_sample' other samples and take their average as a background:
	if(is.null(poi_string)){
		sample_names = colnames(counts_table)
		# all_other_samples = sample_names[is.na(pmatch(sample_names, poi_cols))] # gets the *names* of the other samples,
		# but we want the indices:
		poi_names = colnames(counts_table)[poi_cols]
		all_other_samples = match(sample_names[is.na(pmatch(sample_names, poi_names))], sample_names) #http://stackoverflow.com/questions/2487922/how-can-i-get-the-complement-of-vector-y-in-vector-x
	}else{
		all_other_samples = find_samples(counts_table, poi_string, invert_search=TRUE, index=TRUE)
	}

	# print("All other samples:")
	# print(all_other_samples)

	# print("Cells in POI and Background [SHOULD BE NONE]")
	# print(intersect(all_other_samples, poi_cols))

	if(n_background_samples > 0 ){
		cat(sprintf("Selecting %i background samples.. \n"))
		background_cells = sample(all_other_samples, n_background_samples) #grep("Endo",colnames(counts_table))  
	}else
	{
		cat(sprintf("Comparing against all cells not in the population..\n"))
		background_cells = all_other_samples
	}
	background_cols = background_cells
	background_cells = colnames(counts_table)[background_cols]
	

	cat(sprintf("Using  %i samples as background.. \n", length(background_cells)))
	#cat(sprintf("		%s \n", background_cells))
	

	cat(sprintf("Single cell data is %s ..\n", if(single_cell_data){"ON"}else{"OFF"}))
	if(single_cell_data){
		cat("Loading SCDE.. \n")
		require(scde)
		
		# note, scde requires raw counts. can't use normalised counts.
		if(!is.null(complete_table))
		{
			cat(sprintf("Getting raw counts from complete data frame..\n"))
			raw_counts_table = get_raw_counts(counts_table, complete_table)

		}else if(!is.null(raw.counts)){
			raw_counts_table = data.frame(raw.counts, check.names=F)
		}else
		{
			error("Raw counts are required. Cannot run SCDE!")
			return (NULL)
		}
		

		# switch the background cell columns from normalised to raw:
		background_cells = gsub("normalised", "expected", background_cells, ignore.case=T)
		poi_names = colnames(raw_counts_table)[poi_cols]
		
		if(length(grep("expected_count", poi_names)) > 0)
		{
			poi_names = gsub("expected_count_", "", poi_names, ignore.case=T)
		}

		if(length(grep("expected_count", colnames(raw_counts_table))) > 0)
		{
			colnames(raw_counts_table) = gsub("expected_count_", "", colnames(raw_counts_table), ignore.case=T)
		}

		cat(sprintf("%i cells in population of interest (POI).. \n", length(poi_names)))
		cat(sprintf("%i cells in background group.. \n", length(background_cells)))

		cells_to_compare = c(background_cells, poi_names)
		cat(sprintf("Total cells to compare: %i \n", length(cells_to_compare)))
		# if(verbose)
		# {
		# 	print(cells_to_compare)
		# 	print(cells_to_compare %in% colnames(raw_counts_table))
		# }
		counts_for_scde = raw_counts_table[, cells_to_compare]
		

		# cleanup as per tutorial, i don't think it will do anything:
		# clean up the dataset
		cat(sprintf("Data dims before cleanup: %i x %i \n", nrow(counts_for_scde), ncol(counts_for_scde)))
		
		# omit genes that are never detected
		counts_for_scde <- counts_for_scde[rowSums(counts_for_scde)>0,];
		
		#warn("Not filtering very low coverage cells")
		min_counts = 1e4
		# omit cells with very poor coverage
		counts_for_scde <- counts_for_scde[, colSums(counts_for_scde)>min_counts]; 
		cat(sprintf("Data dims after: %i x %i \n", nrow(counts_for_scde), ncol(counts_for_scde)))

		# set up groups vector, that identifies all cells 
		# as either part of the POI or background.
		
		group.labels = setup.groups(poi_cols, 
								name_of_population, 
								pop2.indices=background_cols, 
								cell.names=colnames(counts_for_scde))
		# poi_count  = 0
		# bkgd_count = 0
		# for(i in 1:length(colnames(counts_for_scde)))
		# {
		# 	cell = colnames(counts_for_scde)[[i]]
		# 	if(cell %in% poi_names)
		# 	{
		# 		gps[[i]] = name_of_population
		# 		poi_count = poi_count + 1
		# 	}else{
		# 		gps[[i]] = "background"
		# 		bkgd_count = bkgd_count + 1
		# 	}
		# 	#cat(sprintf("Marking cell %i of %i [%s] as %s    [N_Background=%i, N_POI=%i].. \n", i, length(cells_to_compare), cell, gps[[i]], bkgd_count, poi_count))
			
		# }
		# gps = factor(gps, levels=c(name_of_population, "background"))
		# cat(sprintf("Groups: \n"))
		# print(gps)
		info("Comparing against background [SCDE]..")
		de.scde = find_differential_genes_scde(counts_for_scde, 
									groups=group.labels, 
									batch = batch,
									seurat.obj=seurat.obj,
									verbose=verbose,
									min_p_val=min_p_val,
									knn=use.knn.error.models,
									annotate.marker.genes=annotate.marker.genes,
									sub.sample.groups.to = sub.sample.groups.to, 
									name_of_population=name_of_population)

		return (de.scde)
	}else{
		library(edgeR)
		warn("edger compatibility is broken, code should arrive here")
		# NOTE. This takes an average within the population. Should probably use edgeR (or DeSeq?)
		# without taking averages.

		# get the expression average of the cells of interest:
		average_colname = paste(name_of_population,"average_counts",sep="_")
		#poi_col_names = colnames(counts_table)[poi_cols]
		#print(colnames(complete_table))
		compare = complete_table[, c("GENE_SYMBOL","ENSEMBL_ID", "transcript_id.s.")]
		cat(sprintf("\n \n \n"))
		cat(sprintf("Calculating average expression.. \n"))
		#print(head(counts_table[poi_cols]))
		compare[average_colname] = rowMeans(counts_table[poi_cols])
		#print(head(compare))
		cat(sprintf("Done. \n"))
		background_colname = paste("background", length(background_cells),"average_counts",sep="_")
		compare[background_colname] = rowMeans(counts_table[background_cells])
		print(head(compare))

		# list counts for 'n_diff_genes' highest expressed genes
		cat(sprintf("\n \n \n"))
		cat(sprintf("%i genes with highest expression in the population of interest.. \n", n_diff_genes))
		most_highly_expressed = head(compare[order(-compare[average_colname]),] , n=n_diff_genes)
		rownames(most_highly_expressed) = most_highly_expressed$GENE_SYMBOL
		most_highly_expressed = most_highly_expressed[c(average_colname, background_colname)]
		print(most_highly_expressed)
		make_heatmap(most_highly_expressed,
				title=paste(n_diff_genes, "most highly expressed genes in", name_of_population), #"Background signal: \n Most highly expressed in the background group", 
				filename="most_highly_expressed.pdf", cluster_rows=FALSE,
				sample_labels=c(name_of_population, "Random sample"))
		write.table(most_highly_expressed, file=paste(n_diff_genes,"_most_highly_expressed.txt", sep=""), sep="\t", row.names=TRUE)

		# call differentially expressed genes
		data_only = compare[c(average_colname, background_colname)]
		#print(data_only)

		#print(compare["GENE_SYMBOL"])

		rownames(data_only) = compare$GENE_SYMBOL
		edger_tested = find_differential_genes_edger(data_only, enriched_only=enriched_genes_only, name_of_population=name_of_population)
		rval = PopInfo(counts=compare,
						#dge_list=dgelist, 
						tested=edger_tested)
		return (rval)
	}
	
}


# find genes differentially expressed between 
# the populations defined by two clusters. Using SCDE or permutation test.
compare.clusters <- function(
							counts_table, 
							clusters,
							compare,
							filter.worst.percent=0,
							n.threads=4,
							use.knn.error.models=F,
							use.perm.test=F,
							o.ifm=NULL,		#precomputed SCDE error models
							o.prior=NULL, 	#precomputed SCDE expression prior
							complete_table			= NULL,
							raw.counts 				= NULL,
							n_background_samples	= 0,
							sub.sample.groups.to    = 0,
							n_diff_genes			= 0, 
							min_p_val 				= 0.05,
							n.randomizations 		= 100,
							show_top				= TRUE, 
							seurat.obj 				= NULL,
							batch 					= NULL,
							verbose 				= TRUE,
							annotate.marker.genes   = F
							)
{
	
	if(length(compare) != 2)
	{
		error("Compare must a vector of exactly two cluster IDs")
		return (FALSE)
	}
	cl1 = as.character(compare[1])
	cl2 = as.character(compare[2])
	info(sprintf("Finding differentially expressed genes between clusters %s and %s.. ", cl1, cl2))

	# note, scde requires raw counts. can't use normalised counts.
	if(!is.null(complete_table))
	{
		cat(sprintf("Getting raw counts from complete data frame..\n"))
		raw_counts_table = get_raw_counts(counts_table, complete_table)

	}else if(!is.null(raw.counts)){
		raw_counts_table = data.frame(raw.counts, check.names=F)
	}else
	{
		error("Raw counts are required. Cannot run SCDE!")
		return (NULL)
	}

	if(length(clusters) != ncol(raw_counts_table)){
		stop(sprintf("Cluster labels [%s] are not the same length as the counts table [%s]!", 
			length(clusters), ncol(raw_counts_table)))
	}
	# print(levels(clusters))
	# print(compare[1])
	# print(typeof(compare[1]))
	#poi -- population of interest
	poi_cols = which(clusters==cl1)
	background_cols = which(clusters == cl2)
	
	info("Counts values range:")
	print(range(raw_counts_table))

	# switch the background cell columns from normalised to raw:
	# background_cells = gsub("normalised", "expected", background_cells, ignore.case=T)
	
	background_cells = colnames(raw_counts_table)[background_cols]
	poi_names = colnames(raw_counts_table)[poi_cols]
	
	# if(length(grep("expected_count", poi_names)) > 0)
	# {
	# 	poi_names = gsub("expected_count_", "", poi_names, ignore.case=T)
	# }

	# if(length(grep("expected_count", colnames(raw_counts_table))) > 0)
	# {
	# 	colnames(raw_counts_table) = gsub("expected_count_", "", colnames(raw_counts_table), ignore.case=T)
	# }

	info(sprintf("%i cells in group 1  \n", length(poi_names)))
	#print(poi_names)

	info(sprintf("%i cells in group 2  \n", length(background_cells)))
	#print(background_cells)

	cells_to_compare = c(background_cells, poi_names)
	cat(sprintf("Total cells to compare: %i \n", length(cells_to_compare)))
	if(verbose)
	{
		print(cells_to_compare)
	}
	#counts_for_scde = raw_counts_table[, cells_to_compare]
	counts_for_scde = raw_counts_table

	# if(!is.null(o.ifm))
	# {
	# 	info("Selecting the appropriate error models from those provided")
	# 	groups <- rep(NA, ncol(raw_counts_table))
	# 	groups[c1] <- cluster1.name
	# 	groups[c2] <- cluster2.name

	# 	o.ifm = o.ifm[rownames(o.ifm) %in% cells_to_compare, ]
	# 	# print(rownames(o.ifm))
	# 	# print(cells_to_compare)

	# 	info("Calculating expression prior")
	# 	o.prior <- scde.expression.prior(models = o.ifm, counts = counts_for_scde, length.out = 400, show.plot = FALSE)
		
	# }

	# cleanup as per tutorial, i don't think it will do anything:
	# clean up the dataset
	cat(sprintf("Data dims before cleanup: %i x %i \n", nrow(counts_for_scde), ncol(counts_for_scde)))
	
	# omit genes that are never detected
	counts_for_scde <- counts_for_scde[rowSums(counts_for_scde) > 0,];
	min_counts = 1e4
	# omit cells with very poor coverage
	counts_for_scde <- counts_for_scde[, colSums(counts_for_scde) > min_counts]; 
	cat(sprintf("Data dims after: %i x %i \n", nrow(counts_for_scde), ncol(counts_for_scde)))

	# set up groups vector, that identifies all cells 
	# as either part of the POI or background or not included (NA)
	group.labels = setup.groups(poi_cols, 
								cl1, 
								cell.names=colnames(counts_for_scde),
								pop2.indices=background_cols, 
								pop2.name=cl2)
	
	# gps = rep(NA, length(colnames(counts_for_scde)))
	# names(gps) = colnames(counts_for_scde) 

	# poi_count  = 0
	# bkgd_count = 0
	# for(i in 1:length(colnames(counts_for_scde)))
	# {
	# 	cell = colnames(counts_for_scde)[[i]]
	# 	if(cell %in% poi_names)
	# 	{
	# 		gps[[i]] = compare[1]
	# 		poi_count = poi_count + 1
	# 	}else{
	# 		gps[[i]] = compare[2]
	# 		bkgd_count = bkgd_count + 1
	# 	}
	# 	#cat(sprintf("Marking cell %i of %i [%s] as %s    [N_Background=%i, N_POI=%i].. \n", i, length(cells_to_compare), cell, gps[[i]], bkgd_count, poi_count))
		
	# }
	# gps = factor(gps, levels=c(compare[1], compare[2]))
	if(use.perm.test)
	{
		de = compare.perm(counts)
	}else{
		de = find_differential_genes_scde(counts_for_scde, 
									n.cores=n.threads,
									o.ifm=o.ifm,
									o.prior=o.prior,
									n.randomizations=n.randomizations,
									knn=use.knn.error.models,
									groups=group.labels, 
									batch = batch,
									seurat.obj=seurat.obj,
									min_p_val=min_p_val,
									verbose=verbose,
									annotate.marker.genes = annotate.marker.genes,
									filter.worst.percent=filter.worst.percent,
									sub.sample.groups.to = sub.sample.groups.to, 
									name_of_population=cl1)
	}
	return(de)
}

#  extract the most different genes from a given cluster
examine_cluster <- function(counts_data, 
							clusters, 
							complete_data	 = NULL,
							verbose 		 = FALSE,
							raw.counts 		 = NULL,
							single_cell_data = TRUE,
							seurat.obj 	 	 = NULL,
							batch 			 = NULL, # if this is non null, and single_cell_data is on, then it is passed to SCDE to correct for batch effects.
							sample.n.cells   = 0,
							min_p_val 		 = 0.05,
							use.pairwise.corr = FALSE,
							compare.most.similar = FALSE, # more specific markers by searching for DE against most similar cluster, takes twice as long.
							most.similar.min.size = 10, # ignore small clusters when comparing to the most similar cluster.
							annotate.marker.genes = FALSE,
							cluster.id 	 = NULL)
{
	if(is.null(cluster.id))
	{
		cat(sprintf("ERROR: Must provide a cluster name to examine it!\n"))
		return (NULL)
	}
	cat(sprintf("Searching for markers in %s ..\n", cluster.id))
	if(compare.most.similar)
	{
		cat(sprintf("Will compare to all cells not in the cluster, and also specifically against most similar cluster of size greater than %s\n",most.similar.min.size))
	}else
	{
		cat(sprintf("Will compare to all cells not in the cluster as background..\n"))
	}
	
	indices_of_cluster = which(clusters == cluster.id)
	cluster_info_fname = paste(cluster.id, "members.txt", sep=" ")
	cluster.members = colnames(counts_data)[indices_of_cluster]
	cat(sprintf("There are %i samples in this cluster..\n", length(cluster.members)))
	cat(sprintf("Writing their names to %s .. \n", cluster_info_fname))
	write(cluster.members, file = cluster_info_fname)
	info("Comparing against background")
	scde.background = examine_subpopulation(counts_data,
									complete_table		 = complete_data,
									raw.counts 			 = raw.counts,
									poi_cols 			 = indices_of_cluster, 
									name_of_population	 = cluster.id,
									n_background_samples = 0, #if zero, take all others as background
									sub.sample.groups.to = sample.n.cells, #if zero, don't subsample (take all)
									n_diff_genes		 = 100, 
									enriched_genes_only	 = TRUE,
									min_p_val 			 = min_p_val,
									annotate.marker.genes = annotate.marker.genes,
									single_cell_data	 = single_cell_data,
									batch 				 = batch,
									seurat.obj			 = seurat.obj,
									verbose  			 = verbose,
									show_top			 = TRUE)
	
	sig.against.background = scde.background[["sig"]]
	sig.against.background.down = scde.background[["sig.down"]]
	info(sprintf("%s genes are significantly [p < %s] upregulated against the background", nrow(sig.against.background), min_p_val))
	rval = sig.against.background
	all = scde.background[["all"]]
	

	## If sig.against.background is NULL, write all genes and return NULL.
	if(is.null(sig.against.background))
	{
		warn("No signficantly upregulated genes against background. Writing complete list and exiting.")
		write.output(all, rval, cluster.id, de.sig.down=sig.against.background.down)
		return (NULL)
	}


	if(compare.most.similar)
	{
		compare.against = 1
		info(sprintf("Comparing against %i most similar cluster(s) [SCDE]", compare.against))
		
		if(!use.pairwise.corr & (nrow(sig.against.background) > 1))
		{
			sim.cluster = get.most.similar.cluster(counts_data, clusters, cluster.id,  marker.genes=sig.against.background$GENE_SYMBOL,
				method="marker.genes", min.cluster.size=most.similar.min.size)
		}else{
			sim.cluster = get.most.similar.cluster(counts_data, clusters, cluster.id,  method="average.pairwise.corr", 
				min.cluster.size=most.similar.min.size)
		}

		info(sprintf("Most similar cluster is: %s", sim.cluster))
		scde.clust.sim = compare.clusters(norm.counts, 
											clusters = clusters, 
											compare=c(cluster.id, sim.cluster), 
											raw.counts = raw.counts, 
											sub.sample.groups.to = sample.n.cells)
		sig.against.sim.cluster = scde.clust.sim[["sig"]]
		
		if(!is.null(sig.against.sim.cluster))
		{
			info(sprintf("%s genes are signficantly [p < %s] upregulated against the most similar cluster", nrow(sig.against.sim.cluster), min_p_val))
			info("Merging marker tables..")

			print(head(sig.against.background, n=2))

			print(head(sig.against.sim.cluster, n=2))

			sig.both = merge(sig.against.background, sig.against.sim.cluster, by="GENE_SYMBOL", suffixes=c(".background", ".sim.clust"))
			sig.both["Mean.Min.LogFC"] = rowMeans(sig.both[,c("Conservative.Estimate.background", "Conservative.Estimate.sim.clust")], na.rm=TRUE)
			sig.both = sig.both[order(sig.both$Mean.Min.LogFC, decreasing=T),]
			info(sprintf("%s genes are signficantly [p < %s] upregulated against both", nrow(sig.both), min_p_val))
			rownames(sig.both) = sig.both$GENE_SYMBOL

			info("Merging complete tables..")
			all.both = merge(scde.background[["all"]], scde.clust.sim[["all"]], by="GENE_SYMBOL", suffixes=c(".background", ".sim.clust"))
			all.both["Mean.Min.LogFC"] = rowMeans(all.both[,c("Conservative.Estimate.background", "Conservative.Estimate.sim.clust")], na.rm=TRUE)
			all.both = all.both[order(all.both$Mean.Min.LogFC, decreasing=T),]
			rownames(all.both) = all.both$GENE_SYMBOL

			#print(head(sig.both))
			rval = sig.both
			all = all.both
		}else
		{
			warn("No genes upregulated against similar cluster. returning only those up against background!")
			rval = sig.against.background
			all = scde.background[["all"]]
		}
	}

	# if annotate.marker.genes is true, this is already run during find_differential_genes_scde
	# if(annotate.marker.genes & !is.null(rval))
	# {
	# 	if(nrow(rval) > 1)
	# 	{
	# 		rval = annotate.markers(rval)
	# 	}
	# }	

	# WARN: if compare against most similar is on, then genes DOWNregulated against
	# the background are still written.
	write.output(all, rval, cluster.id, de.sig.down=sig.against.background.down)
 	return (rval)
}

find_differential_genes_edger <- function(counts_data, 
										sample_names=NULL, 
										n_genes=100, 
										min_p_val=0.05,
										name_of_population="POI",
										enriched_only=FALSE # look only for top enriched genes instead of just most different
									)
{

		library(ggplot2)
		library(edgeR)
		library(scales)
		
		cat(sprintf("Extracting %i differentially expressed genes.. \n", n_genes))
		# cat(sprintf("Data: \n"))
		# print(head(counts_data))
		# if there are replicates, use group to tell edgeR:
		# group <- c(rep("C", 4) , rep("T", 3))
		sample_names = colnames(counts_data)
		group <- c(1:length(sample_names))
		counts_matrix = as.matrix(counts_data)
		dge_list =  DGEList(counts=counts_matrix, group=group)
		#names(y)
		dge_list = calcNormFactors(dge_list)
		
		if(length(sample_names) > 2){
			# Output Multi-Dimensional Scaling Plot -- requires at least 3 samples..
			#An MDS plot measures the similarity of the samples and projects this measure into 2-dimensions
			pdf( "MDS_plot_1_ex1.pdf" , width = 7 , height = 7 ) # in inches
			plotMDS( dge_list , main = "MDS Plot for Count Data", labels = colnames( dge_list$counts ), cex=0.7 )
		}
		

		# Testing for DE genes:
		bcv <- 0.2 # see section 2.10 http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
		et <- exactTest(dge_list,dispersion=bcv^2)
		#print(et)

		#topTags(et, n=n_de, sort.by="p.value")
		cat(sprintf("Saving top DE genes.. \n"))
		if(enriched_only){
			#head(compare[order(-compare[average_colname]),] , n=n_diff_genes)
			# this gets us the genes with the most dramatic fold change, but they are often lowly expressed and not relevant.
			use_pvalue=TRUE
			if(use_pvalue){
				edger_top_tags = et$table[order(et$table["PValue"]),]
				edger_top_tags = edger_top_tags[edger_top_tags["logFC"] < 0, ]
			}else{
				edger_top_tags = et$table[order(et$table["logFC"]),]
				
			}
			top_DE_Genes = dge_list$counts[ rownames( edger_top_tags ) , ]
			


		}else{
			edger_top_tags = topTags( et , n = n_genes, sort.by="logFC" )
			top_DE_Genes = dge_list$counts[ rownames( edger_top_tags$table ) , ]
		}
		print("EdgeR Top Tags:")
		print(head(edger_top_tags, n=50))
		#write.table(edger_top_tags, file=paste(n_genes,"_edger_toptags.txt", sep=""), sep="\t")
		
		
		top_DE_Genes_for_output = merge(top_DE_Genes, edger_top_tags, by='row.names')

		top_DE_Genes_for_output = top_DE_Genes_for_output[top_DE_Genes_for_output$PValue < min_p_val, ]
		cat(sprintf("%i genes are DE (pval < %f) .. \n", nrow(top_DE_Genes_for_output), min_p_val))	
		print(head(top_DE_Genes_for_output, n=10))

		top_DE_Genes_for_output_sorted = top_DE_Genes_for_output[with(top_DE_Genes_for_output, order(PValue)), ]
		
		# cat(sprintf("After sorting, %i genes left.. \n", nrow(top_DE_Genes_for_output_sorted)))
		# print(head(top_DE_Genes_for_output_sorted, n=10))

		print("Top DE genes:")
		print(head(top_DE_Genes_for_output_sorted, n=50))

		colnames(top_DE_Genes_for_output_sorted)[1] <- "GENE_SYMBOL"
		write.table(top_DE_Genes_for_output_sorted, file=paste("top_DE_genes","_p_",min_p_val,".txt",sep=""), sep="\t", row.names=FALSE)
		
		make_heatmap(head(top_DE_Genes,n=n_genes), cluster_rows=FALSE,
			title=paste(n_genes, "top enriched genes in", name_of_population), #"Background signal: \n top 100 genes enriched in background average vs \n random sample of 100 cells", 
			filename="most_different_genes.pdf", sample_labels=c(name_of_population, "Background"))

		return (et)
		
}


add.empty.row <- function(data, name=NULL)
{
	temprow <- t(as.matrix(rep(NA,ncol(data))))
 
	# make it a data.frame and give cols the same names as data
	 
	newrow <- data.frame(temprow)
	colnames(newrow) <- colnames(data)
	 
	# rbind the empty row to data
	 
	 # print(ncol(data))
	 # print(ncol(newrow))
	data <- rbind(data, newrow)


	if(!is.null(name))
	{
		cat(sprintf("Naming new row ..%s \n", name))
		rownames(data)[nrow(data)] = name
	}

	return (data)
}




# to make the list of markers more informative, add the average expression of each gene
# in all clusters..
add.average.expression <- function(norm.counts, cluster.assignments)
{
	
	cat(sprintf("Adding average expression to %i marker lists..\n", length(unique(cluster.assignments))))
	cluster.names = unique(cluster.assignments)
	for(i in 1:length(cluster.names))
    {
        
        current.cluster = cluster.names[i]

        markers.file = paste(current.cluster, paste("markers_", current.cluster, ".txt", sep=""), sep="/")
        info(sprintf("Reading markers for cluster %i of %i from %s..", i, length(cluster.names), markers.file))
        
        if (isTRUE(file.exists(markers.file))) 
        {
        		markers = read.delim(markers.file)
		        
		        # drop some unimportant columns
				markers = markers[, c("GENE_SYMBOL", "Log2.Fold.Change.MLE", "p", "Important.Go.Terms", "Go.Terms", "Conservative.Estimate")]
				colnames(markers) = gsub("Log2.Fold.Change.MLE", "Estimated.Log2.Fold.Enrichment", colnames(markers))


		        cat(sprintf("Looking for expression of %i genes:\n", length(markers$GENE_SYMBOL)))
		        #print(markers$GENE_SYMBOL)
		        genes.indices =  rownames(norm.counts) %in% markers$GENE_SYMBOL

		        not.found = markers$GENE_SYMBOL[!(markers$GENE_SYMBOL %in% rownames(norm.counts)) ]
		        if(length(not.found) > 0)
		        {
		        	print("These genes not found!")
		        	print(not.found)
		        }

		        cat(sprintf("Found expression of %i genes.. \n", length(rownames(norm.counts)[genes.indices])))
		        #print(rownames(norm.counts)[genes.indices])

		        info("Average expression in other clusters..")
		        counts.markers = norm.counts[genes.indices, ]
		        counts.average.markers <- aggregate(t(counts.markers), by=list(cluster.assignments), FUN=mean, na.rm=TRUE)
				rownames(counts.average.markers) = counts.average.markers$Group.1
				counts.average.markers$Group.1 = NULL
				counts.average.markers = t(counts.average.markers)
				


				if(length(not.found) > 0)
				{
					warn("Some genes not found!!")
					for(i in 1:length(not.found))
					{
						missing.gene = not.found[i]
						info(sprintf("Adding empty row for %s ", missing.gene))


						counts.average.markers = add.empty.row(counts.average.markers, name=missing.gene)
						print(tail(counts.average.markers))
					}
				}

				#print(rownames(counts.average.markers))
				info("Ordering..")

				
				print(length(rownames(counts.average.markers)))
				print(length(markers$GENE_SYMBOL))

				# order to match the markers list (probably ranked by fold change or something)
				matches = match(rownames(counts.average.markers), markers$GENE_SYMBOL)
				#matches = match(markers$GENE_SYMBOL, rownames(counts.average.markers))
				#y[sort(order(y)[x])]

				# print("matches:")
				# print(matches)
				# print("max match:")
				# print(max(matches))

				# print(nrow(counts.average.markers))

				counts.average.markers = counts.average.markers[matches, ]
				counts.average.markers = data.frame(counts.average.markers)
				counts.average.markers["GENE_SYMBOL"] = rownames(counts.average.markers)
				with.average = merge(markers, counts.average.markers, by="GENE_SYMBOL")
				
				
				

				if("Mean.Min.LogFC" %in% colnames(with.average))
				{
					info("Sorting by Mean.Min.LogFC")
					with.average = with.average[order(with.average$Mean.Min.LogFC, decreasing=T),]
				}else
				{
					info("Sorting by Conservative.Estimate")
					with.average = with.average[order(with.average$Conservative.Estimate, decreasing=T),]
				}
				
				#print(head(with.average))
				# no need for this column in final output.
				with.average$Conservative.Estimate = NULL

				new.table.name = paste(current.cluster, paste("markers_", current.cluster, "_with_averages", ".txt", sep=""), sep="/")
				
				info("Writing new markers table to:")
				info(sprintf(new.table.name))
				write.table(with.average, file=new.table.name, sep="\t", row.names=F)
		}else
		{
			warn("Markers not found!")
		}

		cat("\n\n")
        
    }
}
    
setup.groups <- function(pop1.indices, pop1.name, cell.names,
			pop2.indices=NULL, 
			pop2.name="Background")
{

	g <- rep(NA, length(cell.names))
	names(g) = cell.names
	g[pop1.indices] <- pop1.name
	if(!is.null(pop2.indices))
	{
		g[pop2.indices] <- pop2.name
	}else
	{
		warn("No background cell indices specified, using all cells")
		all.other.cells = setdiff(1:length(cell.names), pop1.indices)
		g[all.other.cells] = pop2.name
	}
	g = factor(g)
	#info(sprintf("Length of groups vector is %s", length(g)))
	return (g)
}

make_heatmap <- function(data_matrix, 
						title="heatmap", 
						filename="heatmap.pdf",
						cluster_rows=TRUE, 
						cluster_cols=TRUE,
						scale_heatmap_rows=FALSE, 
						sample_labels=NULL)
{
	pdf(file=filename, width=8.5, height=11)
	cat(sprintf("Generating %s heatmap [Log2 counts] \n", title))
	library(gplots)
	data_matrix = data_matrix + 1
	data_matrix = log2(data_matrix)
	m = as.matrix(data_matrix)

	m.col <- cluster(t(m), info_string="column") # cluster colums
	if(cluster_rows){
		m.row <- cluster(m, info_string="row")	# and rows
		m_row_dend = as.dendrogram(m.row)
	}else{
		m_row_dend <- FALSE
	}

	if(scale_heatmap_rows){
		scale = "row"
		key_title = "Z-Score"
	}else{
		scale = "none"
		key_title = "Log2(Counts+1)"
	}

	if(!is.null(sample_labels))
	{
		labels = sample_labels
	}else
	{
		labels = colnames(data_matrix)
	}
	print(key_title)
	hm_diff_genes = heatmap.2(
		m,
		main=title, #paste(n_genes, " most differentially expressed genes"),
		col=rev(colorRampPalette(brewer.pal(11,"RdBu"))(512)),
		Rowv=m_row_dend,
		Colv=as.dendrogram(m.col),
		labRow=row.names(m),
		labCol=labels,
		trace="none",
	        density.info="none",#"histogram",
		scale=scale,
		cexRow=0.5,     # row label size
		key = TRUE,
		key.title = "HUH",
		key.xlab = "Log2(Counts+1)",
		key.ylab = "SUP",
		keysize=1,
	    #key.xlab==key_title,
		cexCol=0.8,       # col label size
		srtCol=45
	)
}

create.matrix <- function(rows, cols, use=NA) {
    
    matrix(rep(use, rows * cols), rows, cols)
   
    }



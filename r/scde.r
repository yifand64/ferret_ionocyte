
library(scde)

scde.test.one.gene <- function(gene, error.models, counts, prior, n.randomizations=100, groups=NULL)
{
	ediff = scde.test.gene.expression.difference(gene, error.models, counts, 
				groups=groups, prior=prior, n.randomizations=n.randomizations)
	ediff$p <- 2*pnorm(abs(ediff$Z), lower.tail=F) # 2-tailed p-value
	ediff$p.adj <- 2*pnorm(abs(ediff$cZ), lower.tail=F) # Adjusted to control for FDR
	return (ediff)

}

# use scde package to find genes differentially expressed in a group of 
# single cells. See http://pklab.med.harvard.edu/scde/Tutorials/
find_differential_genes_scde <- function(counts_data, 
										groups=NULL,
										batch=NULL,
										sample_names=NULL, 
										plot_n = 20,
										plot.top.DE.genes=F,
										o.ifm=NULL, 
										o.prior=NULL,
										n.cores=1,
										n_genes=100, 
										knn=F,
										filter.worst.percent=0,
										min_p_val=0.05,
										min.lower.bound = 0.5,
										n.randomizations=100,
										#min.average.fc = 1.5,
										name_of_population="POI",
										sub.sample.groups.to = 20, # only use to speed testing.
										seurat.obj = NULL,
										enriched_only=FALSE, # look only for top enriched genes instead of just most different
										verbose=F,
										annotate.marker.genes=T
									)
{

	library(scde)

	#subsampling to speed up calculations
	#---------------------------------------------------------------------------------
	cluster.indices = which(groups==name_of_population)#groups[groups==name_of_population]
	background.indices = which(groups!=name_of_population)#groups[groups != name_of_population]
	if(sub.sample.groups.to > 0)
	{
		cat(sprintf("Subsampling groups to %i cells!\n", sub.sample.groups.to))
		# print("Given groups:")
		# print(groups)
		
		if(length(cluster.indices) < sub.sample.groups.to)
		{
			cat(sprintf("No need to downsample cluster, cluster size [%i] is less than specified sub.sample [%i] \n",length(cluster.indices), sub.sample.groups.to))
			keep.in.cluster = cluster.indices
		}else{
			keep.in.cluster = sample(cluster.indices, sub.sample.groups.to)
		}
		if(length(background.indices) < sub.sample.groups.to)
		{
			cat(sprintf("No need to downsample background, background size [%i] is less than specified sub.sample [%i] \n",length(background.indices), sub.sample.groups.to))
			keep.in.background = background.indices
		}else{
			keep.in.background = sample(background.indices, sub.sample.groups.to)
		}

		if(verbose)
		{
			print("Cluster indices [keeping]:")
			print(keep.in.cluster)
			print("Background indices [keeping]:")
			print(keep.in.background)
			print("Keeping: ")
			print(keep)
		}
	}else
	{
		keep.in.cluster = which(groups==name_of_population)#groups[groups==name_of_population]
		keep.in.background = which(groups!=name_of_population)
	}

	if(filter.worst.percent>0)
	{
		warn(sprintf("Filtering out the worst %s%% of cells in each group", filter.worst.percent))
		clust.complexity = colSums(counts_data[,cluster.indices]>0)
		thresh = quantile(clust.complexity, filter.worst.percent/100)
		keep.in.cluster = cluster.indices[clust.complexity > thresh]
		info(sprintf("Dropped %s cells with less than %s genes detected in %s", length(cluster.indices)-length(keep.in.cluster), thresh,name_of_population))


		bkgd.complexity = colSums(counts_data[,background.indices]>0)
		thresh = quantile(bkgd.complexity, filter.worst.percent/100)
		keep.in.background = background.indices[bkgd.complexity > thresh]
		info(sprintf("Dropped %s cells with less than %s genes detected in background", length(background.indices)-length(keep.in.background), thresh))
		
	}

	# if(verbose){
	# 	print("Groups:")
	# 	print(groups)
	# 	print("names(groups):")
	# 	print(names(groups))
	# 	print("counts column names:")
	# 	print(colnames(counts_data))
	# }
	before = length(groups)
	keep = c(keep.in.cluster, keep.in.background)
	counts_data = na.omit(counts_data[, keep])
	groups = groups[keep]
	info(sprintf("Only considering %s of %s cells for this comparison", length(groups), before))


	# calculate models
	if(is.null(groups))
	{
		stop("No Groups!")
	}else{
		cat(sprintf("Using groups: \n"))
		print(table(groups))
		# print(head(groups, n=30))
	}
	

	if(is.null(o.ifm) | is.null(o.prior))
	{		
		scde.models = scde.compute.error.models(counts_data, 
												knn.models=knn,
												max.pairs=10000,
												verbose=verbose, 
												n.cores=n.cores, 
												groups=groups)
		info("Parsing SCDE results")
		o.ifm = scde.models[["error.models"]]
		counts_data = scde.models[["counts"]]
		o.prior = scde.models[["expression.prior"]]
		groups = groups[rownames(o.ifm)]
		counts_data = counts_data[, rownames(o.ifm)]
		
	}else
	{
		info("Using provided SCDE error models and expression prior")
	}

	if(!all(colnames(counts_data) %in% rownames(o.ifm)))
	{
		warn("Error models missing for some cells!")
		print(colnames(counts_data)[!colnames(counts_data) %in% rownames(o.ifm)])
	}else
	{
		if(!all(rownames(o.ifm) %in% colnames(counts_data)))
		{
			warn("Cells missing for some error models !")

			print(rownames(o.ifm)[!rownames(o.ifm) %in% colnames(counts_data)])

			print(head(rownames(o.ifm)))
			print(head(colnames(counts_data)))
		}else
		{
			info("Found error models for all cells")
		}
	}

	info(sprintf("Running differential expression tests using %s cores", n.cores))
	# run differential expression tests on all genes.
	
	if(!is.null(batch))
	{
		info("Using SCDE batch correction..")
		print("batch length:")
		print(length(batch))

		print("groups length")
		print(length(groups))

		ediff <- scde.expression.difference(o.ifm,
										counts_data,
										o.prior,
										groups=groups,
										batch = batch,
										n.randomizations=n.randomizations,
										n.cores=n.cores,
										verbose=1)
	
	}else{
		info(sprintf("Running SCDE expression difference test [n.randomizations=%s]", n.randomizations))
		ediff <- scde.expression.difference(o.ifm,
									counts_data,
									o.prior,
									groups=groups,
									n.randomizations=n.randomizations,
									n.cores=n.cores,
									verbose=1)

	}
	
	
	# top upregulated genes (tail would show top downregulated ones)

	#lb, mle, ub lower bound, maximum likelihood estimate, and upper bound of the 95

	#ce conservative estimate of expression-fold change (equals to the min(abs(c(lb,ub))), or 0 if the CI crosses the 0

	#Z uncorrected Z-score of expression difference

	#cZ expression difference Z-score corrected for multiple hypothesis testing using Holm procedure If 
	#batch correction has been performed (batch has been supplied), analogous data frames are returned 
	#in slots $batch.adjusted for batch-corrected results, and $batch.effect for the differences explained by batch effects alone.
	
	# print("Returned data frame columns: ")
	# print(head(ediff))
	# print(colnames(ediff))
	# print("lb  - lower bound")
	# print("mle - maximum likelihood estimate")
	# print("ub  - upper bound")
	# print("ce  - conservative estimate")
	# print("Z   - uncorrected Z score")
	# print("cZ  - corrected Z score")
	# print("p  - pvalue")
	# print("All quantities refer to log fold change..")

	cat(sprintf("Computing p-values..\n"))
	p.values <- 2*pnorm(abs(ediff$Z), lower.tail=F) # 2-tailed p-value
	p.values.adj <- 2*pnorm(abs(ediff$cZ), lower.tail=F) # Adjusted to control for FDR
	
	# keep only genes signicant below a given p-value
	significant.genes <- which(p.values.adj < min_p_val)
	# add the adjusted p-values to the table
	de.sig <- cbind(ediff[significant.genes, 1:6], p.values.adj[significant.genes])
	
	# rank by conservative estimate (same as lower/upper bound)
	de.sig <- de.sig[order(de.sig$ce, decreasing=T),]

	# add the p-values to the complete table, and rank
	de <- cbind(ediff[,1:6], p.values.adj)
	de <- de[order(de$ce, decreasing=T), ]	

	# debug
	# print("Top of de.sig")
	# print(head(de.sig))

	# print("Bottom:")
	# print(tail(de.sig))

	# grab significantly downregulated genes
	de.sig.down = de.sig[-1 * de.sig$ce > min.lower.bound, ]
	
	# debug
	# print("Top of de.sig.down")
	# print(head(de.sig.down))

	# print("Bottom of de.sig.down")
	# print(tail(de.sig.down))

	de.sig.down = de.sig.down[order(de.sig.down$ce, decreasing=F), ]

	# same for significantly upregulated genes.
	de.sig = de.sig[de.sig$ce > min.lower.bound,]
	
	info(sprintf("%s vs %s comparison:", names(table(groups))[1], names(table(groups))[2])) 
	info(sprintf("%i genes are DE, p<%f  ", nrow(de.sig), min_p_val))
	
	cat("\n")
	info(sprintf("%s vs %s comparison:", names(table(groups))[2], names(table(groups))[1]) ) 
	info(sprintf("%i genes are DE, p<%f ", nrow(de.sig.down), min_p_val))
	cat("\n")

	# give columns more intuitive names.
	cnames = c("Lower.bound","Log2.Fold.Change.MLE","Upper.Bound", "Conservative.Estimate","Z.Score", "Corrected.Z.score", "p")
	colnames(de) <- cnames
	colnames(de.sig) <- cnames
	colnames(de.sig.down) <- cnames

	# move the gene names to a column so the data looks better in e.g excel
	de["GENE_SYMBOL"] = rownames(de)
	de.sig["GENE_SYMBOL"] = rownames(de.sig)
	de.sig.down["GENE_SYMBOL"] = rownames(de.sig.down)
	
	if(nrow(de.sig) > 0){
		if(annotate.marker.genes)
		{
			de.sig = annotate.markers(de.sig)
		}
		if(plot.top.DE.genes){
			info("Plotting top markers..")
			#plot top marker gene:
			for(j in 1:max(plot_n, length(de.sig))){
				top_marker = rownames(de.sig)[j]
				initial_wd = getwd()
				dirname = "Expression prior plots [Top 10 markers]"
				dir.create(dirname, showWarnings = FALSE)
				setwd(dirname)
				if(!is.na(top_marker)){
					cat(sprintf("Generating plot of top marker. %i. %s.. \n", j, top_marker))
					pdf(paste(j, "_ranked_marker_", top_marker, "_", name_of_population, ".pdf", sep=""))
					scde.test.gene.expression.difference(top_marker, 
														models=o.ifm, 
														counts=counts_data, 
														prior=o.prior)
					dev.off()
				}
				setwd(initial_wd)
			}
		}
	info(sprintf("SCDE returning %s DE genes [%s up, %s down]", 
		nrow(de.sig) + nrow(de.sig.down), nrow(de.sig), nrow(de.sig.down)))
	}else{
		warn("No significantly enriched genes detected! (Not doing plots)")
	}
	
	compare_name = paste(names(table(groups))[1], names(table(groups))[2], sep="_vs_")
	rval <- new("SCDE.Results", 
					name=compare_name, 
					sig=de.sig, all=de, 
					sig.down=de.sig.down, min.p=min_p_val, 
					min.lower.bound=min.lower.bound, 
					error.models=o.ifm, 
					expr.prior=o.prior,
					groups=groups)
}


# compute error models
scde.compute.error.models <- function(raw.counts, 
										groups=NULL, 
										save.models.to="scde_error_models.Rdata", 
										save.prior.to="scde_expression_prior.Rdata",
										knn.models=F,
										min_counts = 1e4,
										max.pairs=5000,			
										k=0,
										n.cores=2,
										verbose=TRUE)
{
	cat(sprintf("\n"))
	info("Starting SCDE ")
	info(sprintf("Input counts is %s", paste(dim(raw.counts), collapse=", ")))
	info(sprintf("Using %i cores ", n.cores))
	
	# omit genes that are never detected
	keep.genes = which(rowSums(raw.counts) > 0)
	raw.counts <- raw.counts[keep.genes,];
	
	# omit cells with very poor coverage
	if(min_counts>0){
		info(sprintf("Filtering on coverage. Min Counts/UMIs [%s]", min_counts))
		
		keep.cells = which(colSums(raw.counts) > min_counts)
		raw.counts <- raw.counts[, keep.cells]

		info(sprintf("%s cells passed this filter", ncol(raw.counts)))
	}
	# groups = groups[keep.cells]

	# cat(sprintf("Groups vector has length: %s, and %s cells passed coverage filter\n", length(groups), ncol(raw.counts)))
	# scde expects integer counts
	info("Rounding all counts to integer values")
	rounded.counts = round(raw.counts)
	rounded.counts = apply(rounded.counts, 2, function(x) {storage.mode(x) <- 'integer'; x})
	# print("Input to SCDE:")
	if(verbose == "very"){
		verbose=2
	}else
	{
		if(verbose)
		{
			verbose=1
		}else
		{
			verbose=0
		}
	}
	
	# print(head(raw.counts[,1:5]))
	info(sprintf("Building error models for %s cells", ncol(rounded.counts)))
	info(sprintf("Max pairs (for cross-fitting cell models) set to --> %s", max.pairs))
	#o.ifm is a dataframe with error model coefficients for each cell (rows):
	

	if(knn.models)
	{
		if(k==0){k = round(sqrt(ncol(rounded.counts)))}
		info("Building knn error models")
		info(sprintf("	k = %s", k))
		o.ifm <- knn.error.models( counts=rounded.counts,
								groups=groups,
								k = k,
								n.cores=n.cores,
								save.model.plots=T,
								min.count.threshold = 2, 
								min.nonfailed = 5, 
								max.model.plots = 10,
								verbose=verbose);
	}else
	{
		info("Building error models")
		o.ifm <- scde.error.models( counts=rounded.counts,
								groups=groups, 
								n.cores=n.cores,
								threshold.segmentation=T,
								max.pairs=max.pairs,
								save.model.plots=F,
								save.crossfit.plots=F,
								verbose=verbose);
	}
	
	info("Error models:")
	print(head(o.ifm))

	# filter out cells that don't show positive correlation with
	# the expected expression magnitudes (very poor fits)
	info("Checking valid cells")
	valid.cells <- o.ifm$corr.a >0;
	print(table(valid.cells))
	o.ifm <- o.ifm[valid.cells,];

	info("Estimating gene expression prior..")
	o.prior <- scde.expression.prior(models=o.ifm,
										counts=rounded.counts,
										length.out=400,
										show.plot=F)
	info("Done!")
	rval = list("error.models"=o.ifm, "counts"=rounded.counts, "expression.prior"=o.prior)
	save(o.ifm, file=save.models.to)
	save(o.prior, file=save.prior.to)
	return (rval)
}


# visualise dropout probability:
scde.plot.dropout <- function(o.ifm, o.prior, groups)
{
	cat("Plotting dropout probability..\n")
	o.fail.curves <- scde.failure.probability(o.ifm, magnitudes=log((10^o.prior$x)-1))
	par(mfrow=c(1,1),mar = c(3.5,3.5,0.5,0.5), mgp = c(2.0,0.65,0), cex = 1);
	plot(c(),c(),xlim=range(o.prior$x),ylim=c(0,1),xlab="expression magnitude (log10)",ylab="drop-out probability")

	n.groups = length(unique(groups))
	cols = colorRampPalette(brewer.pal(n.groups, "Set1"))(n.groups)
	i = 1
	for (g in unique(groups))
	{
		color = cols[i]
		info(sprintf("Plotting dropout curve for %s in colour %s", g, color))
		invisible(apply(o.fail.curves[, which(groups==g)], 2, function(y) lines(x = o.prior$x, y = y,col = color)))
		i = i + 1
	}
	
}

# a wrapper to call the three distance types 
scde.dist <- function(raw.counts, o.ifm, o.prior, n.cores=detectCores(), dist.type="reciprocal")
{
	if(identical(dist.type, "reciprocal"))
	{
			rval = t(scde.dist.reciprocal(scde.models[["counts"]], scde.models[["error.models"]], scde.models[["expression.prior"]]))
	}else{
		if(identical(dist.type, "mode.rel"))
		{
			rval = t(scde.dist.mode.rel(scde.models[["counts"]], scde.models[["error.models"]], scde.models[["expression.prior"]]))
		}else
		{
			if(identical(dist.type, "direct"))
			{
				rval = t(scde.dist.direct.dropout(scde.models[["counts"]], scde.models[["error.models"]], scde.models[["expression.prior"]]))
			}
		}
	}
	return (rval)
}

scde.dist.direct.dropout <- function(raw.counts, o.ifm, o.prior, n.cores = detectCores(), n.simulations=1000)
{
	p.self.fail <- scde.failure.probability(models=o.ifm,counts=raw.counts)
	# simulate drop-outs
	# note: using 10 sampling rounds for illustration here. ~500 or more should be used.
	n.simulations <- 10; 
	k <- 0.9;
	require(boot)
	cell.names <- colnames(raw.counts); names(cell.names) <- cell.names;
	dl <- mclapply(1:n.simulations,function(i) {
	  scd1 <- do.call(cbind,lapply(cell.names,function(nam) {
	    x <- raw.counts[,nam];
	    # replace predicted drop outs with NAs
	    x[!as.logical(rbinom(length(x),1,1-p.self.fail[,nam]*k))] <- NA;
	    x;
	    }))
	  rownames(scd1) <- rownames(raw.counts); 
	  # calculate correlation on the complete observation pairs
	  cor(log10(scd1+1),use="pairwise.complete.obs");
	},mc.cores=n.cores)
	# calculate average distance across sampling rounds
	direct.dist <- as.matrix(1-Reduce("+",dl)/length(dl))
	return(direct.dist)
}

scde.dist.reciprocal <- function(raw.counts, o.ifm, o.prior, n.cores=detectCores())
{
	require(boot)
	cell.names <- colnames(raw.counts); names(cell.names) <- cell.names;
	o.fpm <- scde.expression.magnitude(o.ifm,counts=raw.counts);
	k <- 0.95;
	reciprocal.dist <- as.matrix(1-do.call(rbind,mclapply(cell.names,function(nam1) {
	  unlist(lapply(cell.names,function(nam2) {
	    # reciprocal probabilities
	    f1 <- scde.failure.probability(models=o.ifm[nam1,,drop=F],magnitudes=o.fpm[,nam2])
	    f2 <- scde.failure.probability(models=o.ifm[nam2,,drop=F],magnitudes=o.fpm[,nam1])
	    # weight factor
	    pnf <- sqrt((1-f1)*(1-f2))*k +(1-k); 
	    boot::corr(log10(cbind(raw.counts[,nam1],raw.counts[,nam2])+1),w=pnf)
	    }))
	},mc.cores=n.cores)),upper=F)
	return(reciprocal.dist)
}

# get the 'mode-relative' weighted distance between cells..
scde.dist.mode.rel <- function(raw.counts, o.ifm, o.prior, verbose=T)
{
	require(boot)
	if(verbose)
	{	
		print("Model:")
		print(head(o.ifm))

		print("Counts:")
		print(head(raw.counts))

		print("Expression prior:")
		print(head(o.prior))
	}
		# reclculate posteriors with the individual posterior modes 
	jp <- scde.posteriors(models=o.ifm,raw.counts,o.prior,return.individual.posterior.modes=T,n.cores=detectCores())
	# find joint posterior modes for each gene - a measure of MLE of group-average expression
	jp$jp.modes <- log(as.numeric(colnames(jp$jp)))[max.col(jp$jp)]
	p.mode.fail <- scde.failure.probability(models=o.ifm,magnitudes=jp$jp.modes)
	
	cat("Calculating dropout probabilities..\n")
	# get failure probabilities on the expresison range
	o.fail.curves <- scde.failure.probability(o.ifm,magnitudes=log((10^o.prior$x)-1))
	# get self-fail probabilities (at a given observed count)
	p.self.fail <- scde.failure.probability(models=o.ifm,counts=raw.counts)

	# weight matrix
	matw <- 1-sqrt(p.self.fail*sqrt(p.self.fail*p.mode.fail))
	# magnitude matrix (using individual posterior modes here)
	mat <- log10(exp(jp$modes)+1);

	cell.names <- colnames(raw.counts); names(cell.names) <- cell.names;
	# weighted distance
	cat("Generating weighted distance matrix..\n")
	#mode.fail.dist <- as.dist(1-do.call(rbind, mclapply(cell.names,function(nam1) {
	mode.fail.dist <- as.matrix(1-do.call(rbind, mclapply(cell.names,function(nam1) {
	  unlist(lapply(cell.names, function(nam2) {
	    corr(cbind(mat[,nam1],mat[,nam2]),w=sqrt(sqrt(matw[,nam1]*matw[,nam2])))
	  }))
	},mc.cores=detectCores())),upper=F);
	return (mode.fail.dist)
}

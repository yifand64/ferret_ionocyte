library(plyr)
library(crayon)
library(reshape)
library(RColorBrewer)

extract.field=function(string,field=1,delim="_", fixed=T) {
	return(strsplit(string,delim, fixed=fixed)[[1]][field])
}


# # source all the files in a directory. not used.
# source.dir <- function (path, pattern = "\\.[rR]$", env = NULL, chdir = TRUE) 
# {
#     files <- sort(dir(path, pattern, full.names = TRUE))
#     lapply(files, source, chdir = chdir)
# }

normalise.quantile <- function(counts, min.reads=5, within.lane=T)
{
	library(EDASeq)
	info("Using EDA quantile normalisation")
	#dds = DESeqDataSetFromMatrix(countData = counts, colData = NULL, design = NULL)
	
	data <- newSeqExpressionSet(counts=as.matrix(round(counts)))
	
	if(within.lane){data <- withinLaneNormalization(data,"gc", which="full")}
	dataNorm <- betweenLaneNormalization(data, which="full", offset=TRUE)
	normCounts(dataNorm)
}

winsorize <- function (x, fraction=0.05) {
   if(length(fraction) != 1 || fraction < 0 ||
         fraction > 0.5) {
      stop("bad value for 'fraction'")
   }
   lim <- quantile(x, probs=c(fraction, 1-fraction))
   x[ x < lim[1] ] <- lim[1]
   x[ x > lim[2] ] <- lim[2]
   x
}

RowVar_fast <- function(x) {
  library(Matrix)
  Matrix::rowSums((x - Matrix::rowMeans(x))^2)/(dim(x)[2] - 1)
}

aggregate_matrix_sum = function(x, group_by)
{
	rv = NULL
	for(g in unique(group_by)){
		if(is.null(rv)){rv = Matrix::colSums(x[group_by==g,])}else{
			rv = rbind(rv, Matrix::colSums(x[group_by==g,]))
		}
	}
	rv
}

get.variable.genes.umis <- function(umi.cts, residual.threshold=-0.25, UMIs.threshold=0, use.spline=F, batch=NULL, ret.plot=F, fit.spline=T, verbose=F)
{
	library(Matrix)
	library(mgcv)
	if(!is.null(batch))
	{
		v = as.vector(table(batch))
		total_transcripts = data.frame(as.matrix(t(Matrix.utils::aggregate.Matrix(t(umi.cts), groupings = batch, fun="sum"))))
		#total_transcripts = umi.cts %>% group_by(cyl) %>% summarize_all(.funs = sum) 
		detection_frac = Matrix.utils::aggregate.Matrix(t(umi.cts > 0), groupings = batch, fun="sum")
		detection_frac = data.frame(as.matrix(t(detection_frac / v)))
		test_genes = rownames(detection_frac)[Matrix::rowSums(detection_frac > 0) == length(unique(batch))]
		detection_frac = detection_frac[test_genes, ]
		total_transcripts = total_transcripts[test_genes, ]
		detection_frac$gene = rownames(detection_frac)
		total_transcripts$gene = rownames(total_transcripts)
		detection_frac = melt(detection_frac, id.vars="gene")
		colnames(detection_frac) = c("gene", "batch", "alpha")
		total_transcripts = melt(total_transcripts, id.vars="gene")
		colnames(total_transcripts) = c("gene", "batch", "UMIs")
		z = cbind(total_transcripts, detection_frac)[, c("gene", "batch", "alpha", "UMIs")]
		if(verbose) info("Fitting logistic GLM (controlling for batch covariate)")
		model.logit = glm(data = z, formula = alpha ~ log10(UMIs) + batch, family = binomial)
		#model.logit = robust::glmRob(data = z, formula = alpha ~ log10(UMIs), family = binomial)
		if(fit.spline){
			if(verbose) info("Fitting spline quantile regression (controlling for batch covariate)")
			model.gam = quantreg::rq(data = z, formula = alpha ~ splines::ns(log10(UMIs), df=15) + batch, tau=0.8)
			#model.gam = mgcv::gam(data = z, formula = alpha ~ s(log10(UMIs)), method="REML")
		}
	}else{
		if(verbose) info("Computing gene dection rates (alphas)..")
		z = data.frame(UMIs = Matrix::rowSums(umi.cts), alpha= Matrix::rowSums(umi.cts>0) / ncol(umi.cts))
		z = subset(z, UMIs > 0 | alpha > 0)
		if(verbose) info("Fitting GLMs..")
		model.logit = glm(data = z, formula = alpha ~ log10(UMIs), family = binomial)
		#model.logit = robust::glmRob(data = z, formula = alpha ~ log10(UMIs), family = binomial)
		if(fit.spline){
			model.gam = quantreg::rq(data = z, formula = alpha ~ splines::ns(log10(UMIs), df=15), tau=0.8)
			#model.gam = mgcv::gam(data = z, formula = alpha ~ s(log10(UMIs)), method="REML")
		}
	}
	
	
	if(use.spline & fit.spline){
		if(verbose) info("use.spline is ON. Using GAM fit (blue), logit in red")
		z$predicted = predict(object = model.gam, z, type="response")
		z$predicted.alternate = predict(object = model.logit, z, type="response")
		z$residual = model.gam$residuals
	}else{
		if(verbose)  info("use.spline is OFF. Using logit fit (blue), GAM in red")
		z$predicted = predict(model.logit, type="response") #predict(object = model.logit, z, type="response")
		z$residual = residuals(model.logit, type="response") #model.logit$residuals
		if(fit.spline){z$predicted.alternate = predict(object = model.gam, z, type="response")}
	}
	if(is.null(batch)) {z$gene = rownames(z)}
	outliers = subset(z, residual < residual.threshold & UMIs > UMIs.threshold)
	g = ggplot(z, aes(x=log10(UMIs), y=alpha, label=gene)) + geom_point(color="grey50", size=0.5, stroke=0) + 
		ylim(c(0,1)) + geom_line(aes(y=predicted), size=0.5, color="blue", linetype="dotted")  + 
		geom_text(data=outliers, color="black", size=1.5, vjust=2)
	if(fit.spline &! use.spline){geom_line(aes(y=predicted.alternate), size=0.5, color="red", linetype="dotted")}
	if(!is.null(batch)){g = g + facet_wrap(~batch)}
	if(!ret.plot){print(g)}
	rv = unique(unlist(lapply(rownames(outliers), extract.field, 1, delim="_")))
	if(ret.plot){return(list("var.genes"=rv, "plot"=g, "fit.data"=z, "logit"=model.logit))}
	rv
}

# from http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
get.variable.genes.brennecke <- function(ed, winsorize=F, show.n=0, fit.all=F, min.cv2=2, pdf=NULL, width=9, height=8, do.plot=T, p.thresh=0.05, verbose=T)
{
	library(statmod)
	library(Matrix)
	if(winsorize)
	{
		info("Winsorizing")
		ed <- t(apply(ed, 1, winsorize, fraction=2/ncol(ed)))
	}
	if(verbose) info("Calculating mean")
	means <- Matrix::rowMeans(ed)
	if(verbose) info("Calculating variance")
	vars <- RowVar_fast(ed)
	cv2 <- vars/means^2
	if(verbose) info("Fitting regression line")
	if(fit.all)
	{
		useForFit = names(means)
	}else
	{
		if(verbose) info("Selecting highly expressed")
		minMeanForFit <- unname( quantile( means[ which( cv2 > min.cv2 ) ], .95 ) )
		
		useForFit <- means >= minMeanForFit # & spikeins
		if(verbose) info(sprintf("Fitting only the %s genes with mean expression > %s", sum(useForFit), minMeanForFit))
	}
	fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ), cv2[useForFit] )
	a0 <- unname( fit$coefficients["a0"] )
	a1 <- unname( fit$coefficients["a1tilde"])
	#fit$coefficients
	if(do.plot){
		if(!is.null(pdf))
		{
			if(verbose) info(sprintf("Rendering variable genes fit to -->  %s", pdf))
			pdf(pdf, width=width, height=height)
		}
		par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2));
		
		if(verbose) info("Plotting")
	}
	# get fit line
	xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
	vfit <- a1/xg + a0
	# add fit line
	if(do.plot){lines( log(xg), log(vfit), col="black", lwd=3 )}
	
	df <- ncol(ed) - 1
	# add confidence interval
	if(do.plot){
		lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
		lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")
	}
	afit <- a1/means+a0
	varFitRatio <- vars/(afit*means^2)
	varorder <- order(varFitRatio, decreasing=T)
	oed <- ed[varorder,]
	pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
	adj.pval <- p.adjust(pval,"fdr")
	r = data.frame(rownames(ed), varFitRatio, pval, adj.pval)
	colnames(r) = c("Gene", "VarianceFitRatio", "p", "p.adj")
	v = r[!is.na(r$p.adj),]
	n.sig = sum(v$p.adj<p.thresh)
	if(verbose) info(sprintf("Found %s variable genes (p<0.05)", n.sig))

	# add top 100 genes
	if(do.plot){
		if(show.n>0){
			if(verbose) info(sprintf("Showing top %s genes", show.n))
			points(log(means[varorder[1:show.n]]),log(cv2[varorder[1:show.n]]),col=2)
		}else{
			if(verbose) info(sprintf("Showing all %s genes significant p<%s", n.sig, p.thresh))
			points(log(means[varorder[1:n.sig]]),log(cv2[varorder[1:n.sig]]),col=2)
		}
		
		if(!is.null(pdf)){
			dev.off()
		}
	}

	r = r[order(r$VarianceFitRatio, decreasing=T), ]
	r$Rank = 1:nrow(r)
	return(r)
	#return(oed)
}


# return the top 'n_genes' most variable/highly expressed genes from this
# subset of a counts table. 
# rank them by scoring the sum of their ranks. i.e the gene that is the 5th most 
# variable and 7th highest expressed will have score 12.
get.top <- function(counts, n_genes)
{
	vars = rowSums((counts - rowMeans(counts))^2)/(dim(counts)[2] - 1)
    means = rowMeans(counts)
    vars.rank = rank(-vars)
    means.rank = rank(-means)
    score = vars.rank + 0.5*means.rank
    df = data.frame(cbind(vars, means, vars.rank, means.rank, score))
    rownames(df) = rownames(counts)
    df = df[with(df, order(score)), ]
    #print(head(df))
    return (head(rownames(df), n_genes))
}

# grab a total of n_genes from a set of clusters
# where the clusters are all big enough to sample approximately equally,
# but exactly equal sampling will not give the correct number of genes.
take.equally <- function(counts, n_genes, clusters, rank=T)
{
	n_clusters = length(unique(clusters))
	n_each = n_genes / n_clusters
	q = n_genes - (floor(n_each) * n_clusters)
	p = n_clusters - q
	#print(table(clusters))
	info(sprintf("Will sample either %s from %s clusters and %s from the other %s", floor(n_each), p, ceiling(n_each), q))
	# found q and p by solving some simultaneous equations.
	sampled = list()
	for(i in 1:n_clusters)
    {
    	cluster.id = names(table(clusters))[i]
    	all.in.cluster = which(clusters == cluster.id)
    	if(i > p)
    	{
    		take.this.many = ceiling(n_each)
    	}else
    	{
    		take.this.many = floor(n_each)
    	}
    	if(rank)
    	{
    		sampled.from.cluster = match(get.top(counts[all.in.cluster, ], take.this.many), rownames(counts))
		}else{
			sampled.from.cluster = sample(all.in.cluster, take.this.many)
		}
    	
    	info(sprintf("	Sampled %s of the %s genes in cluster %s.. ", length(sampled.from.cluster), length(all.in.cluster), cluster.id))
    	# print(head(rownames(counts)[sampled.from.cluster], n=10))
    	sampled = unlist(c(sampled, sampled.from.cluster))
    }
    return(sampled)
}

# if sample.equally is ON, attempt to represent genes from each cluster as
# equally as possible. if off, cluster sizes are preserved.
# note. min.variance will take the fraction of most variable clusters.
# i.e 0.9 takes the 90% most variable, 0.1 the 10%
downsample.genes <- function(counts, 
								n_genes=0, 
								dfrac=0.5, 
								sample.equally=T, 
								min.cluster.size=60,
								rank.in.each.cluster=T, #if false, sample randomly. if true, pick the most variable/highly expressed in each.
								genes_of_interest=c("Lgr5", "Defa5", "Tff3", 
									"Anpep", "Smoc2", "Olfm4", "Muc2", "Mmp7", "Trpm5", "Scin", "Hck", "Car4", "A930009A15Rik", "Avil", "Vav1", "Defa4"),
								min.mean=0,
								min.variance=0)
{
	if(n_genes == 0)
	{
		n_genes = floor(nrow(counts) * dfrac)
	}
	info(sprintf("Downsampling from %s to %s rows", nrow(counts), n_genes))
	
	info("Clustering..")
	hc = cluster.hierarchical(counts, take.log=F, min.cluster.size=min.cluster.size)
	clusters = hc$hc.row
	names(clusters) = rownames(counts)
	#print(clusters)
	for(g in genes_of_interest)
	{
		info(sprintf("Gene: %s is in cluster %s", g, clusters[g]))
	}

    n_clusters = length(unique(clusters))
    size.of.smallest = min(table(clusters))

    mean.expr = aggregate(counts, by=list(clusters), FUN=mean, na.rm=TRUE)
    cluster.var = rowSums((mean.expr - rowMeans(mean.expr))^2)/(dim(mean.expr)[2] - 1)
    cluster.mean = rowMeans(mean.expr)
    cluster.size = table(clusters)

    cluster.info = data.frame(cbind(cluster.var, cluster.mean, cluster.size))
    rownames(cluster.info) = seq(1:n_clusters)
    info(sprintf("%s clusters detected. Smallest cluster has size: %s \n  ", n_clusters, size.of.smallest))

    print(cluster.info)
    cat("\n")

    if(min.variance > 0 | min.mean > 0)
    {
	    min.variance.level = quantile(cluster.info$cluster.var,  1-min.variance)
	    min.expr.level = quantile(cluster.info$cluster.mean,  1-min.mean)
		    # info(sprintf("Quantile = %f", q))
		    # nc = nc[nc$variance > q, ]

	    keep.clusters = rownames(cluster.info[cluster.var > min.variance.level & cluster.mean > min.expr.level, ])
	    

	    info(sprintf("%s of %s clusters (%s) are above the thresholds [variance=%s, mean=%s]..", length(unique(keep.clusters)),length(unique(clusters)), paste(keep.clusters, collapse=","), min.variance.level, min.expr.level))

	    var.counts = counts[clusters %in% keep.clusters, ] # throw away genes in low-variance clusters
	    if(nrow(var.counts) < n_genes)
		{
			error(sprintf("Variance threshold [top %s percent] is too stringent. Not enough genes in variable clusters [%s]!", min.variance*100, nrow(var.counts)))
			return (FALSE)
		}

	    clusters = clusters[match(rownames(var.counts), rownames(counts))]
	    counts = var.counts

	    for(g in genes_of_interest)
	    {
	    	info(sprintf("%s in the variable clusters -- %s", g, g %in% rownames(counts)))
	    }
	    info(sprintf("Sampling from the %s genes in the %s variable clusters..", nrow(counts), length(unique(clusters))))
    }

    

    n_clusters = length(unique(clusters))
    n_extra = n_genes %% n_clusters
    n_each = n_genes / n_clusters
    size.of.smallest = min(table(clusters))

    if(sample.equally)
    {
    	sampled = list()
    	if(size.of.smallest >= n_each)
	    {
	    	if(n_extra == 0)
	    	{
	    		info(sprintf("Will take exactly %s genes from each cluster..", n_each))
	    		
			    for(i in 1:n_clusters)
			    {
			    	all.in.cluster = which(clusters == i)
			    	if(rank.in.each.cluster)
			    	{
			    		sampled.from.cluster = match(get.top(counts[all.in.cluster, ], n_each), rownames(counts))
		    		}else
		    		{
		    			sampled.from.cluster = sample(all.in.cluster, n_each)
		    		}
			    	
			    	
			    	info(sprintf("Sampled %s from cluster %s.. ", length(sampled.from.cluster), i))
			    	sampled = unlist(c(sampled, sampled.from.cluster))
			    }
	    	}else
	    	{
		    	sampled = take.equally(counts, n_genes, clusters)
	    	}	    	
		}else
		{
			info("Cannot sample equally from clusters.")
			take.all = list()
			taken.already = 0
			cl = clusters
			for(i in 1:n_clusters)
	    	{
	    		n.cl = length(unique(cl))
	    		n_remaining = n_genes - taken.already
	    		needed.for.equal.sampling = floor(n_remaining / n.cl)
	    		smallest.size = min(table(cl))
	    		smallest.cluster = names(table(cl)) [which.min(table(cl))]
	    		#print(table(cl))
	    		#info(sprintf("Length of cluster vector: %s. Clusters: %s. %s genes taken from small clusters. %s genes needed for equal sampling. Smallest cluster: %s [%s genes] ", length(cl), n.cl, taken.already, needed.for.equal.sampling, smallest.cluster, smallest.size))
	    		
	    		#print(take.all)
	    		if(smallest.size <= needed.for.equal.sampling)
	    		{
	    			cl = cl[which(cl != smallest.cluster )] # remove genes in the smallest cluster.
	    			take.all = unlist(c(take.all, smallest.cluster)) # add the cluster to the list of clusters from whom to take all genes.
	    			taken.already = taken.already + length(which(clusters==smallest.cluster))
	    			#info(sprintf("taken already: %s", taken.already))
	    		}else
	    		{
	    			# take all genes from the take.all clusters, and take needed.for.equal.sampling from the rest of the clusters.
	    			sampled = list()
	    			#print("sup")
	    			take.equally =  unique(clusters)[which(!(unique(clusters) %in% take.all))]
	    			info(sprintf("Taking all genes from %s clusters [%s], and sampling equally from the %s others [%s]..", length(take.all), paste(take.all, collapse=","), length(take.equally), paste(take.equally, collapse=",")))
	    			if(length(take.all) > 0)
	    			{
	    			for(j in 1:length(take.all))
		    			{
		    				cluster.to.take = take.all[[j]]
		    				genes.to.take = which(clusters==take.all[[j]])
		    				#print(cluster.to.take)
		    				info(sprintf("	Taking all %s genes in cluster %s ..", length(genes.to.take), cluster.to.take))
		    				sampled = unlist(c(sampled, genes.to.take))
		    			}
	    			}
	    			sampled.equally = take.equally(counts, n_remaining, clusters[which(clusters %in% take.equally)])
	    			return (counts[c(sampled, sampled.equally), ])
	    		}
	    	}
    	}
    	return (counts[sampled, ])
	}else
	{
		warn("Sampling unequally not implemented yet!")
		return (NULL)
	}
}

drop.non.varying.genes <- function(norm.counts, rel.min.var=0.5, min.var=NULL)
{
	if(is.null(min.var))
	{
		info(sprintf("Keeping only the %f percent most variable genes..", 100*(1-min.var)))
		use.rel = T
	}else
	{
		info(sprintf("Keeping only genes with variance > %f", min.var))
		use.rel = F
	}
	
    variance = row.var(norm.counts)
    n.before = nrow(norm.counts)
    nc = data.frame(norm.counts)
    #info("Calculating variance..")
    nc["variance"] = variance
    
    if(use.rel)
    {
    	q = quantile(nc$variance,  min.var)
	    info(sprintf("Quantile = %f", q))
	    nc = nc[nc$variance > q, ]
    }else{
    	nc = nc[nc$variance > min.var, ]
    }
    
    nc["variance"] = NULL
    #print(head(norm.counts))
    rval = as.matrix(nc)

    n.dropped =  n.before - nrow(rval)
    info(sprintf("Dropped %i genes.. using %i most variable for further analysis..", n.dropped, nrow(rval)))
    return (rval)
}

check.paths <- function(files)
{
	ok = file.exists(files)
	if(!(all(ok)))
	{
		nope = files[which(ok == FALSE)] 
		error("Some files do not exist!")
		error(nope)
		return (FALSE)
	}
	if(any(duplicated(files)))
	{
		nope = files[which(duplicated(files))]
		error("Some files are duplicates!")
		error(nope)
		return(FALSE)
	}
	return (TRUE)
}


# utilities to load multiple rna-seq datasets and feed them (together)  into
# the rna-seq pipeline.
# samples file should have the format:
# Name	Counts.file QC.file
load.batch <- function(samples.file, 
						exclude=c("background", "population","No_cell", "Nocell", "No.cell", "Empty"),
						write.merged=T, 
						use.edger=F,
						use.tpm=F,
						n.samples.only=0,
						keep.all.genes=TRUE, #if false, genes that aren't detect in all samples are dropped
						plot.draw=T,
						plot.dir=".",
						min.trans.mapped=10,
						min.genes=2000,
						min.cells=5,
						is.expr=1,
						parallel.read=F,
						trim.samples.names=T,
						n.rows.only=0) #n.rows.only > 0 loads only the first n rows, for debug
{
	
	info("Loading samples file..")
	samples.info = read.delim(samples.file)
	n.samples = nrow(samples.info)
	counts.files = as.character(unlist(samples.info$Counts.file))
	qc.files = as.character(unlist(samples.info$QC.file))
	if(!(check.paths(counts.files) & check.paths(qc.files)))
	{
		error("Exiting. Fix samples file")
		quit()
	}else
	{
		info("Sample file ok!")
	}
	info(sprintf("Will load and merge %i datasets..\n", n.samples))
	info(sprintf("Loading counts tables.. [PARALLEL=%s]", parallel.read))
	counts = list()
	completes = list()
	if(use.tpm)
	{
		info("Use TPM is ON")
		count_string = "TPM"
	}else
	{
		info("Use TPM is OFF (using normalised count)")
		count_string = "expected_count" 
	}
	if(parallel.read)
	{
		print(counts.files )
		tables <- mclapply(counts.files, load_counts, 
											count_string=count_string,
											trim.sample.names=trim.samples.names, 
											exclude = exclude,
											n.rows.only=n.rows.only, 
											mc.cores = detectCores())
		for(i in 1:n.samples)
		{
			tb = tables[[i]]
			counts[[i]] = tb[["counts"]]
			print(head(rownames(tb[["counts"]])))
			completes[[i]] = tb[["data"]]
		}
	}else
	{
		for(i in 1:n.samples)
		{
			tb = load_counts(as.character(samples.info$Counts.file[i]), 
											trim.sample.names=trim.samples.names, 
											exclude = exclude,
											n.rows.only=n.rows.only)
			print("Loaded dataset:")
			
			counts[[i]] = tb[["counts"]]
			completes[[i]] = tb[["data"]]
		}
	}
	
	info("Loaded the following datasets:")
	all.cell.names = vector()
	for(i in 1:n.samples)
	{
		dataset.name = samples.info$Name[i]
		cat(sprintf("%s: ", dataset.name))
		counts[[i]] = data.frame(counts[[i]], check.names=F)
		
		print(dim(counts[[i]]))
		# cat("Range of Lgr5:\n")
		# print(range(counts[[i]]["Lgr5", ]))

		counts[[i]] ["GENE_SYMBOL"] = rownames(counts[[i]]) #  needed for the merge
		all.cell.names = c(all.cell.names, colnames(counts[[i]]))
		print(grep("Lgr5", counts[[i]]$GENE_SYMBOL, value=T, ignore.case=T))
	}

	info("Building groups vector..")
	n.samples.in.each.dataset = as.numeric(lapply(counts, ncol))
	n.samples.in.each.dataset = n.samples.in.each.dataset-1 #compensate for the gene_symbol column we had to add
	groups = factor(unlist(mapply(rep, samples.info$Name, n.samples.in.each.dataset))) 

	info("Merging counts tables..")

	# be careful not to convert characters (eg +, - in column names, as it can cause problems matching with QC file)
	merge.all <- function(x, y) {
    	merge(x, y, by="GENE_SYMBOL", all=keep.all.genes, check.names=F)
	}

	counts.merged <- Reduce(merge.all, counts)
	completes.merged <- Reduce(merge.all, completes)
	rownames(counts.merged) = counts.merged$GENE_SYMBOL
	counts.merged$GENE_SYMBOL = NULL
	completes.merged$ENSEMBL_ID = completes.merged$ENSEMBL_ID.x
	completes.merged$ENSEMBL_ID.x = NULL
	completes.merged$ENSEMBL_ID.y = NULL

	# probably unneccesary
	counts.merged = factorsNumeric(counts.merged)
	counts.merged = as.matrix(counts.merged)

	#print(head(colnames(counts.merged)))
	n.cells = length(all.cell.names)
	
	missing.cells = all.cell.names[!all.cell.names %in% colnames(counts.merged)]
	if(length(which(missing.cells != "GENE_SYMBOL"))>0)
	{
		error("Cells were lost in the merge!")
		error("Missing cells:")
		
		for(cn in missing.cells)
		{
			error(cn)
		}
	}else
	{
		info(sprintf("Successfully merged %s datasets. Merged dataset contains %s cells [%s genes]", n.samples, ncol(counts.merged), nrow(counts.merged)))
	}
	print(grep("A310", colnames(counts.merged), value=T))

	# for genes that aren't detected convert NA to 0
	if(keep.all.genes)
    {
            counts.merged[is.na(counts.merged)] <- 0
			completes.merged[is.na(completes.merged)] <- 0
    }

	print(dim(counts.merged))
	print(head(rownames(counts.merged)))

	info(sprintf("Loading %s QC tables..", n.samples))
	qc = list()
	for(i in 1:n.samples)
	{
		qc.each = read.delim(as.character(samples.info$QC.file[i]))
		# info("Using read'rs read_tsv (faster)")
		# warn("read_tsv will not convert % to X. so looking for e.g X.trans_mapped will fail")
		# qc.each = factorsNumeric(read_tsv(as.character(samples.info$QC.file[i])))
		print("Loaded QC data")
		# print(head(qc.each$Sample))
		# print(dim(qc.each))
		qc[[i]] = qc.each
	}

	info("Concatenating QC tables..")
	qc.merged =  do.call("rbind", qc)


	info(sprintf("Filtering. Keep cells with more than %s genes detected above %s counts", min.genes, is.expr))
	info(sprintf("		Keep genes expressed above %s counts in at least %s cells..", is.expr, min.cells))
	counts.merged.filtered = filter.counts(counts.merged, 
											qc=qc.merged,
											min.trans.mapped=min.trans.mapped,
											min.genes=min.genes, 
											min.cells=min.cells, 
											is.expr=is.expr, 
											draw.plot=plot.draw, 
											colour=groups)
	info("After filtering:")
	print(grep("A310", colnames(counts.merged.filtered), value=T))

	keep.samples = match(colnames(counts.merged.filtered), colnames(counts.merged))
	groups = groups[keep.samples]

	# if(write.merged)
	# {
	# 	info("Writing complete normalised counts (all genes) table..")
	# 	write.table(normalise.dseq(counts.merged), file="merged.norm.counts.all.genes.txt", sep="\t", quote=F)
	# }

	counts.merged = counts.merged.filtered

	info(sprintf("%s cells and %s genes passed filters.. ", ncol(counts.merged), nrow(counts.merged)))

	# print(head(counts.merged))

	if(write.merged)
	{
		write.table(counts.merged, file="merged.raw.counts.txt", sep="\t", quote=F)
	}


	cat(sprintf("Normalising across the %i datasets..\n", n.samples))
	if(use.edger)
	{
		warn("The edgeR pseudo counts inflate zeros a little bit which i'm uncomfortable with. Probably use DESeq2.")
		counts.merged.norm = normalise.edgr.tmm(counts.merged, groups=groups, write.norm.factors=T)
	}else
	{	
		if(use.tpm)
		{
			info("Use TPM is ON. Not normalising")
			counts.merged.norm = counts.merged
		}else{
			info("Use TPM is OFF. Normalising!")
			counts.merged.norm = normalise.dseq(counts.merged)
		}
		
	}

	if(write.merged)
	{
		write.table(counts.merged.norm, file="merged.norm.counts.txt", sep="\t", quote=F)
	}

	

	
	if(n.samples.only > 0)
	{
		keep.samples = sample(ncol(counts.merged), n.samples.only)
		counts.merged = counts.merged[, keep.samples]
		counts.merged.norm = counts.merged.norm[, keep.samples]
		groups = groups[keep.samples]
	}


	rval = list("norm.counts" = counts.merged.norm, "raw.counts" = counts.merged, "qc" = qc.merged, "groups"=groups, "data"=completes.merged) 
	return (rval)
}


load_counts <- function(expression_data_file, 
						count_string="expected_count", 
						exclude="nothing", 
						use.row.names=F,
						drop.genes.file=NULL,
						n.rows.only=0,
						counts.only=FALSE, #if this is true return only a matrix of raw RSEM counts (no FPKM, TPM etc.)
						trim.sample.names=FALSE)
{
	info(sprintf("Loading expression data from  --> %s", expression_data_file))

	if(n.rows.only > 0)
	{
		# data <- fread(expression_data_file, nrow=n.rows.only)
		# as.data.frame.matrix(data) 
		
		data <- read.delim(expression_data_file, nrow=n.rows.only, comment.char="", check.names=F)
		
	}else{
		# data <- fread(expression_data_file, nrow=n.rows.only)
		# data = as.data.frame.matrix(data) 
		#data <- read.delim(expression_data_file)
		if(use.row.names)
		{
			data <- read.delim(expression_data_file, check.names=F)
		}else
		{
			
			# this method is much faster, but I cant get it to stop throwing NAs into the data.
			# so stuck with read.delim
			#info("Using readr's read_tsv (faster)")
			#data <- data.frame(read_tsv(expression_data_file, na="NA"))
			data <- read.delim(expression_data_file, check.names=F)
			info(sprintf("%s NA values in table", sum(is.na(data))))
				
		}
		
	}

	if(!("GENE_SYMBOL" %in% colnames(data)))
	{
		warn("No GENE_SYMBOL column found. Assuming that the rownames are gene names, and the table is only counts..")
		data["GENE_SYMBOL"] = rownames(data)
		# warn("Added symbols: ")
		# print(head(data$GENE_SYMBOL))
		counts.df = TRUE
	}else
	{
		counts.df = FALSE 
		info(sprintf("Parsing counts columns from data frame using the count_string: %s", count_string))
	}

	
	n.before = nrow(data)
	data <- data[!duplicated(data[,"GENE_SYMBOL"]),] 
	n.after = nrow(data)
	if(n.after != n.before)
	{
		warn(sprintf("Dropped %s duplicates GENE_SYMBOL ", n.before-n.after))
	}

	if("transcript_id.s." %in% colnames(data))
	{
		info("Dropping duplicate transcript_id.s..")
		n.before = nrow(data)
		data <- data[!duplicated(data[,"transcript_id.s."]),] 
		n.after = nrow(data)
		if(n.after != n.before)
		{
			warn(sprintf("Dropped %s duplicates transcript_id ", n.before-n.after))
		}
	}
	info(sprintf("%i genes left.. ", nrow(data)))

	if(counts.df)
	{
		data = data.matrix(data)
		# is.nan.data.frame <- function(x)
  #               do.call(cbind, lapply(x, is.nan))
        data[is.na(data)] <- 0 

        print("Returning counts table, Nans?")
        print(any(is.na(data)))


		info(sprintf("Dropping zero rows \n"))
	  	non.zero.rows = which(rowSums(data) > 0)
		data = data[non.zero.rows,]

		info(sprintf("Dropping zero columns \n"))
	  	non.zero.cols = which(colSums(data) > 0)
		data = data[, non.zero.cols]

		info(sprintf("Dropping NA named columns \n"))
		data = data[, !is.na(colnames(data))]

		print(dim(data))
	  	print(corner(data))

		return (list("counts" = data, "data" = NULL)  )
	}
	
	cols = colnames(data)
	#print(head(cols))
	samples = cols[ grep(paste(count_string, collapse="|"), cols)]
	samples.excluded = samples[ grep(paste(exclude, collapse="|", ignore.case=TRUE), samples)]

	warn("Excluding these samples [matching the 'exclude' string]: ")
	print(samples.excluded)

	# drop the samples matching the 'exclude' string
	samples = samples[ grep(paste(exclude, collapse="|"), invert=TRUE, ignore.case=TRUE, samples)]
	
	cat(sprintf("Found %i samples: \n", length(samples)))
	print(samples)
	cat(sprintf("Excluding [DROPPING] %i samples: \n", length(samples.excluded)))
	print(samples.excluded)
	
	x = strsplit(samples, count_string)

	if(!is.null(drop.genes.file))
	{
		cat(sprintf("Loading list of genes to drop from f=%s ..\n", drop.genes.file))
		drop_info = read.delim(drop.genes.file)
		print(head(drop_info))
		print("Will drop.")
		drop_these = drop_info["GENE_SYMBOL"]
		print(head(drop_these))
		initial_genes = data["GENE_SYMBOL"]
		print("initial genes:")
		print(head(initial_genes))
		before = nrow(data)
		data <- data[ ! data$GENE_SYMBOL %in% drop_these, ]

		# data = data[ initial_genes[is.na(pmatch(initial_genes, drop_these))], ]
		cat(sprintf("Dropped %i genes found in the specified file.. ", before-nrow(data)))
		print(intersect(drop_these, initial_genes))
	}

	keeps <- c("GENE_SYMBOL", samples)

	#cat(sprintf("\"Selecting column:  %s\"\n", keeps))
	cat(sprintf("Excluding samples:  \n"))
	cat(sprintf("	 %s \n", exclude))
	counts_df <- data[keeps]
	data_cols <- samples
  	
	
	# # ensure data is numeric:
	# counts_df <- factorsNumeric(counts_df)
	

	# print("Looking for Lgr5")
	# counts_df = data.frame(counts_df, check.names=F)
	# print(grep("Lgr5", counts_df$GENE_SYMBOL, value=T))
	# print(sum(unlist(counts_df[counts_df$GENE_SYMBOL == "Lgr5", data_cols])))
	
	print("# NA values in table:")
	print(sum(is.na(counts_df)))	# cat(sprintf("Naming counts matrix rows with GENE_SYMBOL: \n"))
	# print(head(data$GENE_SYMBOL))

	# cat(sprintf("Before dropping rows, %i rows in complete, %i rows in counts .. \n", nrow(data), nrow(counts_df)))
	
	n.before = nrow(data)
	non.zero.rows = which(rowSums(counts_df[data_cols]) > 0)
	counts_df = counts_df[non.zero.rows,]
  	data = data[non.zero.rows,]
	n.after = nrow(data)
	if(n.after != n.before)
	{
		warn(sprintf("Dropped %s undetected genes (zero rows) ", n.before-n.after))
	}

	# print("Looking for Lgr5")
	# print(grep("Lgr5", counts_df$GENE_SYMBOL, value=T))
  	
  	n.before = ncol(data)
	counts_df = counts_df[data_cols]
  	non.zero.cols = colSums(abs(counts_df)) != 0
	counts_df = counts_df[, non.zero.cols]
  	data = data[, non.zero.cols]
	n.after = ncol(data)
	if(n.after != n.before)
	{
		warn(sprintf("Dropped %s samples with no genes (zero cols) ", n.before-n.after))
	}

	rownames(counts_df) = data$GENE_SYMBOL	
	counts_matrix = as.matrix(counts_df)
	print(corner(counts_matrix))
	# print(grep("Lgr5", rownames(counts_matrix), ignore.case=T))
	# cat("Range of Lgr5:\n")
	# print(range(counts_matrix["Lgr5", ]))


	if(trim.sample.names)
	{
		remove_start=function(string,start) return(strsplit(string,start)[[1]][2])
		c = paste(count_string, "_", sep="")
		colnames(counts_matrix) = factor(unlist(lapply(colnames(counts_matrix),remove_start,c)))

	}

	# info("Replacing strange charcacters with underscores..")
	# # clean up weird characters
	# #print(colnames(counts_matrix))
	# #print("============================ \n\n")
	# colnames(counts_matrix) = gsub(".", "", colnames(counts_matrix), fixed=TRUE)
	# colnames(counts_matrix) = gsub(" ", "", colnames(counts_matrix), fixed=TRUE)
	# #print(colnames(counts_matrix))

	
	if(counts.only)
	{
		return(counts_matrix)
	}
	rval  = list("counts" = counts_matrix, "data" = data)  

	return (rval)
}

get_raw_counts = function(counts, complete)
{
  	cat(sprintf("Retrieving raw counts.. \n"))
  	# print(head(colnames(counts)))
  	# print(head(colnames(complete)))
  	# print(tail(colnames(complete)))
  	# print("expected" %in% colnames(complete))
  	if("normalised" %in% colnames(counts)) # does this work?! careful
  	{
  		raw <- complete[ gsub("normalised", "expected", colnames(counts), ignore.case=T)]
	}else
	{
		if("expected" %in% colnames(counts)) # does this work?! careful
		{
			raw <- complete[ paste("expected_count_", colnames(counts), sep="")]
		}else
		{
			if(length(grep("expected", colnames(complete))) > 0)
			{
				raw <- complete[, grep("expected", colnames(complete))]
			}else
			{
				warn("Raw counts not provided! Returning NULL..")
				return (NULL)
			}
			
		}
		
	}
  	raw <- factorsNumeric(raw)
  	#cat(sprintf("L1=%i, L2=%i", nrow(raw), length(complete$GENE_SYMBOL)))
  	rownames(raw) = complete$GENE_SYMBOL
  
  	print("Dropping zero rows..")
  	raw = raw[which(rowSums(raw) > 0),]
  
  	print("Dropping zero columns..")
  	raw = raw[, colSums(abs(raw)) != 0]
  	
  	raw = as.matrix(raw)
  	return (raw)
}

load_qc <- function(qc_data_file)
{
  	cat(sprintf("Loading QC data from  --> %s\n", qc_data_file))
  	qc <- read.delim(qc_data_file)
  	#qc = factorsNumeric(qc)
  	cat(sprintf("Loaded QC data: \n"))
  	print(head(qc$Sample))
  	return (qc)
}

drop.genes.in.file= function(norm.counts, file)
{
    genes.table =  read.delim(file)
    genes.drop = genes.table$GENE_SYMBOL
    cat(sprintf("Dropping the given %i genes.. \n", length(genes.drop)))
    genes.drop = gsub(" ", "", genes.drop, fixed = TRUE)
    rval = norm.counts[!rownames(norm.counts) %in% genes.drop, ]
    return (rval)
}


normalise.edgr.tmm = function(counts, groups=NULL, write.norm.factors=F, use.pseudo.counts=F){
	library(edgeR)
	print("TMM normalizing [edgeR]..")
	# when there's no replicates, we need to hide this fact so edgeR can normalize:
	if(!is.null(groups))
	{
		y2 = DGEList(counts=counts, group=groups)
	}else
	{
		y2 = DGEList(counts=counts)
	}
	y2 = calcNormFactors(y2)
	print(y2$samples) # show normalization factors:
	if(write.norm.factors)
	{
		write.table(y2$samples, file="edger_norm_factors.txt", sep="\t")
	}
	y2 = estimateCommonDisp(y2, verbose=TRUE)
	if(use.pseudo.counts)
	{
		return (y2$pseudo.counts)
	}else
	{
		rval = 1e6*(counts / (y2$samples$lib.size * y2$samples$norm.factors))
	}
	
}

normalise.edgr.tmm = function(counts){
	print("TMM normalizing [edgeR]..")
	# when there's no replicates, we need to hide this fact so edgeR can normalize:
	y2 = DGEList(counts=counts)
	y2 = calcNormFactors(y2)
	
	print(y2$samples) # show normalization factors:
	y2 = estimateCommonDisp(y2, verbose=TRUE)
	return (y2$pseudo.counts)
}

# From Dillies et al (2012)
# http://bib.oxfordjournals.org/content/early/2012/09/15/bib.bbs046.full.pdf+html
# DESeq: This normalization method [14] is included in the DESeq Bioconductor package 
# (ver- sion 1.6.0) [14] and is based on the hypothesis that most genes are not DE. 
# A DESeq scaling factor for a given lane is computed as the median of the ratio, 
# for each gene, of its read count over its geometric mean across all lanes. 
# The underlying idea is that non-DE genes should have similar read counts across samples, 
# leading to a ratio of 1. Assuming most genes are not DE, the median of this ratio for 
# the lane provides an estimate of the correction factor that should be applied to all 
# read counts of this lane to fulfill the hypothesis.
normalise.dseq = function(counts, verbose=F)
{
	library(DESeq2)
	cat(sprintf("Normalising with DESeq2.. \n"))
	geoMeans <- apply(counts, 1, function(row) if (base::all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
	size.factors <- estimateSizeFactorsForMatrix(counts, geoMeans=geoMeans)
	if(verbose)
	{
		info("Size factors:")
		print(size.factors)
	}
	norm.counts = t( t(counts) / size.factors )
	return (norm.counts)
}

tpm <- function(counts, mult=10000, pseudocount=0)
{
	info("Running TPM normalisation")
	total.counts = colSums(counts)
	scaled.counts = t(t(counts+pseudocount) / total.counts) # this needs work! don't run with !=0 
	scaled.counts * mult
}

# From Dillies et al (2012)
# http://bib.oxfordjournals.org/content/early/2012/09/15/bib.bbs046.full.pdf+html
# Total count (TC): Gene counts are divided by the total number of mapped reads 
# (or library size) associated with their lane and multiplied by the mean total 
# count across all the samples of the dataset. Can use mean or median.
normalise.total.count <- function(counts, use="mean")
{
	total.counts = colSums(counts)
	scaled.counts = t(t(counts) / total.counts)
	if(use=="mean")
	{
		scale.factor = mean(total.counts)
	}else
	{
		scale.factor = median(total.counts)
	}
	info(sprintf("Normalised to %s total counts value [%s]", use, round(scale.factor, 2)))
	return ( t(t(scaled.counts) * scale.factor))
}

normalise.downsample.to.min <- function(counts)
{
	
}

rlog.dseq = function(counts)
{
	cat(sprintf("Running DESeq2's rlog transform.. \n"))
	SummarizedExperiment
	return (norm.counts)
}



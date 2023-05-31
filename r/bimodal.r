# utilities for analysing expression profiles of genes in single cell data,
# particularly searching for bimodally expressed genes.

# library(matrixStats)
# library(plyr)
# library(data.table)

row.var <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

fit.log.normal = function(data, 
				do.plot=T,
				pdf.output=T)
{
	library(fitdistrplus)

	f <- fitdist(data,"lnorm")
	if(do.plot)
	{
		plot(f)
	}
	return (f)
}


# for an expression matrix return the coefficient of variation of each gene
# across the samples (cells). 
# if the groups argument is provided, a matrix is returned (n_cols=n_groups)
#, where each column is the vector of CVs within a given group.
coeff.variation= function(m, groups=NULL, get.means=F)
{
	if(is.null(groups))
	{
		stds = rowSds(as.matrix(m), na.rm=TRUE)
		means = rowMeans(m)
		rval = stds/means
		return (rval)
	}

	agg.means <-aggregate(t(m), by=list(groups), FUN=mean, na.rm=TRUE)
	agg.means["Group.1"] = NULL
	agg.means = as.matrix(agg.means)

	agg.stds <- aggregate(t(m), by=list(groups), FUN=sd, na.rm=TRUE)
	agg.stds["Group.1"] = NULL
	agg.stds = as.matrix(agg.stds)


	print("Means:")
	print(head(agg.means[, 1:5]))

	print("Stds:")
	print(head(agg.stds[, 1:5]))
	
	cv = agg.stds / agg.means
	rownames(cv) = unique(groups)


	rval = cv
	if(get.means)
	{
		rownames(agg.stds) = unique(groups)
		rownames(agg.means) = unique(groups)
		rval = list("cv" = cv, "mean"=agg.means, "std"=agg.stds)
	}
	return (rval)
}

# use data.table to aggregate using the alpha value
# for each gene in each cluster
find.alpha  <-function(m, groups, threshold = 1)
{
	#m["GENE_SYMBOL"] = rownames(m)
	m = m > threshold
	agg.alpha <-aggregate(t(m), by=list(groups), FUN=sum)
	rownames(agg.alpha) = agg.alpha$Group.1
	agg.alpha["Group.1"] = NULL
	#agg.alpha["n.cells.in.cluster"] = n.in.groups
	return(agg.alpha)
}

# given a table of counts and a gene of interest, draw 
# a histogram of the gene expression over the samples. if
# groups and a group of interest is provided, then restrict
# to a 
make.hist <- function(counts, gene, groups=NULL, goi=NULL)
{

}


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

###
# eventually i'll implement rahul's method (described in supp. info)
# for the paper with alex shalek (2014). for now, calculate, for each
# gene, in each given group (cluster, probably) the CV and alpha, and
# rank them by the harmonic mean of 1/alpha and average expression
# 
###
find.bimodal.genes <- function(m, groups, output.dir="bimodal_genes", method="cv", min.alpha=0.1, max.alpha=0.8, write.tables=T)
{
	cat("Finding CVs..\n")
	n.in.groups = table(groups)
	var = coeff.variation(m, groups=groups, get.means=T)
	cv = var[["cv"]]
	mean = var[["mean"]]
	# print("mean:")
	# print(mean)
	std = var[["std"]]

	cv = cv[ order(row.names(cv)), ]
	mean = mean[ order(row.names(mean)), ]
	std = std[ order(row.names(std)), ]

	# print(head(cv))
	# calculate alpha -- fraction of the cells in a given group expressing a given
	# gene above a threshold.
	
	alpha = find.alpha(m, groups)
	
	alpha.scaled = alpha / n.in.groups
	if(identical(method, "cv"))
	{
		cat(sprintf("Generating tables of bimodal genes.."))
		for(i in 1:length(unique(groups)))
		{
			cluster.name = row.names(cv)[i]
			cat(sprintf("For %s..\n ", cluster.name))
			alpha.c = alpha[i, ]
			alpha.scaled.c = alpha.scaled[i, ]
			mean.c = mean[i, ]
			std.c = std[i, ]
			cv.c = cv[i, ]
			bimodal.c = data.frame(t(rbind(cv.c, alpha.c, alpha.scaled.c, mean.c, std.c)))
			colnames(bimodal.c) = c("CV", "N.cells", "Alpha", "Mean.expression.log2.count","Std.Dev")
			bimodal.c["Cv.expression"] = apply(bimodal.c[,c('CV', 'Mean.expression.log2.count')], 1, function(x) gm_mean(x) )
			#rowMeans(subset(bimodal.c, select = c(CV, Mean.expression.log2.count)), na.rm = TRUE)
			bimodal.c = bimodal.c[bimodal.c$Alpha > min.alpha, ]
			bimodal.c = bimodal.c[bimodal.c$Alpha < max.alpha, ]
			bimodal.c = bimodal.c[with(bimodal.c, order(-Cv.expression)), ]
			if(write.tables)
			{
				write.table(bimodal.c, file=sprintf("Cluster_%s_bimodal.txt",cluster.name), sep="\t")
			}
			print(head(bimodal.c))
		}
	}else{
		cat(sprintf("ERROR: Unknown method: %s \n", method))
	}



	# each column of cv contains the
}

vis.coeff.hists = function(m, 
							groups,
							output.file="cv_in_clusters.pdf",
							output.pdf=T,
							do.hists=T,
							do.violins=T)
{
	cat(sprintf("Visualising histograms of CVs..\n"))
	cv = coeff.variation(m, groups=groups)
	cv.m = melt(cv)
	if(output.pdf)
	{
		pdf(output.file, height=11, width=8.5)
	}
	if(do.hists)
	{
		hists = ggplot(cv.m, aes(x=value)) + geom_histogram() + 
			facet_wrap(~X1) + theme_bw()
		print(hists)
	}
	
	if(do.violins)
	{
		violins = ggplot(cv.m, aes(factor(X1), value)) + 
			geom_violin(scale="width", aes(fill = factor(X1))) +
			theme_bw() +  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
			coord_fixed()
		print(violins)
	}
	if(output.pdf)
	{
		dev.off()
	}
}




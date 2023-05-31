
require(Hmisc)
require(propagate)
get_covarying_genes <- function(seurat.obj, gene, method="pearson", min.p=0.05, cor.matrix=NULL)
{
	if(is.null(seurat.obj) & is.null(cor.matrix))
	{
		error("Must provide either a counts table or a correlation matrix.")
		return(FALSE)
	}	
	
	if(gene %in% rownames(seurat.obj@data))
	{
		gene.index = which(rownames(seurat.obj@data)==gene)
		counts = seurat.obj@data	
	}else
	{
		counts = rbind(seurat.obj@data, t(fetch.data(seurat.obj, gene)))
		gene.index = which(rownames(counts)==gene)
	}
	
	if(length(gene.index) == 0)
	{
		error(sprintf("Gene %s not found!", gene))
		return (FALSE)
	}

	if(is.null(cor.matrix))
	{
		info(sprintf("Using %s correlation to search for genes covarying with %s..", method, gene))
		info("Calculating covariance matrix..") 
		#cor = rcorr(as.matrix(t(counts)))
		cor.matrix = bigcor(as.matrix(t(counts)), method=method)
		colnames(cor.matrix) = rownames(counts)
		rownames(cor.matrix) = rownames(counts)
		info(sprintf("Built %s correlation matrix", dim(cor)))
		#print(head(cor[, 1:5]))
		#print(dim(cor))
	}
	

	#Pull out the vector of correlations:
	corr.data = cor.matrix[gene.index, ]
	names(corr.data) = rownames(counts)
	corr.data = data.frame(corr.data)
	colnames(corr.data) = method
	
	print(head(corr.data))
	
	

	# Calculate p values
	info("Calculating p values:")
	pvals = rep(0, nrow(counts))
	lower.bound = rep(0, nrow(counts))
	upper.bound = rep(0, nrow(counts))
	#gene.expr = as.numeric(counts[gene, ])
	gene.expr = as.numeric(unlist(fetch.data(seurat.obj, gene)))
	#print(gene.expr)
	for (i in 1:nrow(counts)) {
		gene.current = as.numeric(counts[i, ])
		r = corr.data[i, method]
		# if(r > 0.5)
		# {
			if(i %% 1000 == 0)
			{
				cat(sprintf("Calculating p for gene %i.. [r=%f]\n", i, r))
			}
			
			# print(dim(gene.expr))
			# print(dim(gene.current))
			# print(head(gene.expr))
			# print(head(gene.current))
			c = cor.test(gene.expr, gene.current, method=method)
			#print(c)
			pvals[i] = c$p.value
			if(method=="pearson"){
				lower.bound[i] = c$conf.int[1]
				upper.bound[i] = c$conf.int[2]
			}
			
		# }else
		# {
		# 	#cat(sprintf("Gene %i not correlated. [r=%f]\n", i, r))
		# 	pvals[i] = 0
		#	lower.bound[i] = 0
		#	upper.bound[i] = 0
		# }
	}
	corr.data["p"] = pvals
	if(method=="pearson")
	{
		corr.data["Lower.bound.95.CI]"] = lower.bound
		corr.data["Upper.bound.95.CI]"] = upper.bound
	}
	
	corr.data = data.frame(corr.data[order(-corr.data[,method]), ])
	print(head(corr.data))
	
	return (list("cor.mat"=cor.matrix, "corr.data"=corr.data))

}

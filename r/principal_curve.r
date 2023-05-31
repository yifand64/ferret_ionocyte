
# fit a curve to a tsne mapping (that has a main curve hopefully)
# and find genes that vary along it.
fit.principal.curve <- function(seurat.obj=NULL, data=NULL, method="spearman")
{
	library(princurve)
	pdf("principal_curve.pdf")
	if(is.null(data)){
		if(!is.null(seurat.obj))
		{	
			fit1 = principal.curve(as.matrix(seurat.obj@tsne.rot), plot.true = T)
			seurat.obj@data.info$Curve.Lambda = fit1$lambda
			curve.order = order(fit1$lambda, decreasing=T)
			ordered.counts = seurat.obj@data[,curve.order]
		}else{
			error("Provide a seurat object (with tSNE) or some data!")
		}
	}else{
		fit1 = principal.curve(as.matrix(data), plot.true = T)
		Curve.Lambda = fit1$lambda
		curve.order = order(fit1$lambda, decreasing=T)
		ordered.counts = data[,curve.order]
	}
	dev.off()
	
	coeff = 0.1
	linearly.increasing = coeff * 1:ncol(ordered.counts)
	

	pval = vector()
	corr = vector()
	lower.bound = vector()
	upper.bound = vector()
	for (i in 1:nrow(ordered.counts)) {
		
		x = as.numeric(ordered.counts [i, ])		 
		ct = cor.test(x, linearly.increasing, method=method)
		corr[i] = ct$estimate
		pval[i] = ct$p.value

		if(i %% 1000 == 0)
		{
			cat(sprintf("Calculating p for gene %i.. [r=%f]\n", i, corr[i]))
		}
	
		if(method=="pearson"){
			lower.bound[i] = ct$conf.int[1]
			upper.bound[i] = ct$conf.int[2]
		}
		
	}
	rval = data.frame(corr, pval)
	rownames(rval) = rownames(ordered.counts)
	
	if(method=="pearson")
	{
		rval[,"Lower.bound.95.CI]"] = lower.bound
		corr.data["Upper.bound.95.CI]"] = upper.bound
	}
	
	
	rval = list("corr.data"=rval[order(-rval[, "corr"]), ], "ordered.counts"=ordered.counts)

	if(!is.null(seurat.obj)){
		draw.varying.along.curve(rval, seurat.obj, curve.order)
	}
	
	return(rval)
}

# rval from above function (fit.principal.curve)
draw.varying.along.curve <- function(rval, seurat.obj, curve.order, cluster.cols.use=NULL, genes.for.heatmap=25, cols.use=NULL, clip=3)
{
	if(is.null(cluster.cols.use))
	{
		ncols = length(unique(seurat.obj@ident))
		cluster.cols.use = intense.cols(ncols) #colorRampPalette(brewer.pal(ncols, "Set1"))(ncols)
	}
	if(is.null(cols.use))
	{
		cols.use = rev(brewer.pal(9, "RdBu"))
	}
	

	# annotate the genes
	gene.type = c(rep("Upregulated", genes.for.heatmap), rep("Downregulated", genes.for.heatmap))
	gene.colors = data.frame(gene.type)
	colnames(gene.colors) = "Gene Type"
	gene.color.annotations = gene.colors
	gene.color.pal = colorRampPalette(brewer.pal(length(unique(gene.type)), "Set3"))(length(unique(gene.type)))

	# annotate the cells
	sample.colors = data.frame(factor(seurat.obj@ident[curve.order])) #, "Plate")])
	colnames(sample.colors) = c("Cell Type")#, "Plate")
	#sample.colors$Cluster.ID = factor(sample.colors$Cluster.ID)
	ann.colors = list("Cell Type"=cluster.cols.use, "Gene Type"=gene.color.pal)
	hmap.genes = c(head(rownames(rval$corr.data), n=genes.for.heatmap), tail(rownames(rval$corr.data), n=genes.for.heatmap))
	plot.data = t(scale(t(rval$ordered.counts[ hmap.genes, ])))
	if(clip != 0)
	{
		plot.data[plot.data > clip] = clip
		plot.data[plot.data < -clip] = -clip
	}
	par(mar=c(10.1, 10.1, 10.1, 10.1))
	nmf.options(grid.patch=TRUE) #avoid an extra blank page in the pdf.
	aheatmap(plot.data, Colv=NA, Rowv=NA, annColors = ann.colors, 
		annRow=gene.colors, annCol=sample.colors, cexCol=0, color=cols.use, filename="varying_genes.pdf", width=14, height=8)
}
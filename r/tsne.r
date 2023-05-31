# Generate visualisation of dimensionality reduction by tSNE and find
# clusters using DBSCAN / affinity propagation

library(tsne)
library(Rtsne)
library(fpc)
source("seurat.r")

# assumes the input data matrix has samples as columns and 
# measurements (genes) as rows  
convert_to_pca_scores <- function(data_matrix, pcs.use=1:20, center=T, scale=F, do.fast=T){
 	if(do.fast){
 		library(gmodels)
 		info("Running PCA [fast is ON] ")
		pca = fast.prcomp(data_matrix, center=center, scale=scale)
	}else{
		info("Running PCA [fast is OFF] ")
		pca = prcomp(data_matrix, center=center, scale=scale)
	}
 	#print("PCA object: ")
 	info("Using PCs:")
 	print(pcs.use)
 	rval = pca$x[, pcs.use]
 	#print(dim(rval))
 	return (rval)
}


fisheye.distort <- function(rotations)
{
	x = as.numeric(rotations[, 1])
	y = as.numeric(rotations[, 2])
	r = atan2(sqrt(x*x+y*y), 1);
	r = r/pi;  #/* -0.5 .. 0.5 */
	phi = atan2(y,x);

	u = r * cos(phi) + 0.5;	
	v = r * sin(phi) + 0.5;
	rval = cbind(u,v)
	rownames(rval) = rownames(rotations)
	return (rval)
}

# note, counts are usually either PCA rotations 
# or an adjusted distance matrix from SCDE
# for a discussion of possible inputs to tSNE see:
# https://lvdmaaten.github.io/tsne/#faq
do.tsne <- function(counts, 
					weights=NULL, # provide weights for the dimensions
					is.dist=F, #if true, convert to a distance matrix
					n.cores=1, 
					check_duplicates=T,
					plot.genes=c("Lgr5", "Cell.Cycle_mean_normalised", 
               							"MHC.II_mean_normalised", "Lyz1", "Defa3",
               							"Alpi", "Lct", "Fabp6",
               							"Chgb", "Neurog3",
               							"Fcgbp", "Muc2", "Dclk1", 
               							"Nrep", "Rac2"),
					#snow=F,
					n.runs=1, 
					#parallel=T, 
					whiten=F,
					max.iterations=200, 
					perplexity=15,  ### can be a vector.
					cluster.type="MPI",
					verbose=T,
					dbclust=F,
					seurat.obj=NULL, 
					barnes.hut=T)
{
	ptm <- proc.time()
	info(sprintf("About to run tSNE [%s times] on seurat object", n.runs))
	if(!is.null(weights))
	{
		if(length(weights) != ncol(counts))
		{
			error(sprintf("Length of weight vector [%s] must be the same as the number of dimensions in use [%s]", length(weights), ncol(counts)))
			return(F)
		}
		info("Using given weights")
		counts = sweep(counts,2,weights,`/`)
	}
	print(seurat.obj)
	info(sprintf("Barnes-Hut approximation is: %s ", barnes.hut))

	# make the necessary directories.
	for(j in 1:length(perplexity))
    {
        perp = perplexity[j]
        info(sprintf("Perplexity--%i  [%i of %i]\n", perp, j, length(perplexity)))
        dirname = sprintf("tSNE_perp_%i", perp)
        dir.create(dirname, showWarnings = FALSE)
    }

    run.id = 1:n.runs
    run.perplexity = perplexity
    runs = expand.grid(run.id, run.perplexity)
    colnames(runs) = c("Run.ID", "Perplexity")
    runs$Working.Dir = paste0("tSNE_perp_", runs$Perplexity)

	if(n.cores > 1)
	{
	# 	if(!snow)
	# 	{
			library(foreach)
			library(doParallel)
			cl<-makePSOCKcluster(n.cores, outfile="")
			registerDoParallel(cl, n.cores)
			
#			doParallel(cores=n.cores, outfile="")
			info(sprintf("Running in parallel [doParallel] using %s cores \n", n.cores))
			
			foreach(run.index=1:nrow(runs), .export=c("tsne.iteration", "info", "set.tsne", 
				"feature.plot.scale", "default.cols", "material.heat", "kelly.cols", "kelly"), .packages=c("apcluster",
				"tsne","Rtsne","fpc","ggplot2","RColorBrewer","plyr", "futile.logger", "shape"))  %dopar%
			{
				run.info = runs[run.index,]
				tsne.iteration(counts, run.info[["Run.ID"]], tsne_max_iters=max.iterations, verbose=verbose, is.dist=is.dist,
					tsne_perplexity=run.info[["Perplexity"]], seurat.obj=seurat.obj, barnes.hut=barnes.hut, plot.genes=plot.genes,
					working.dir=run.info[["Working.Dir"]])
			}
			stopCluster(cl)
			registerDoSEQ() # unregister the parallel backend	
	}
	else
	{
		info("Running sequentially")
		for(run.index in 1:nrow(runs)) {

			run.info = runs[run.index,]
			tsne.iteration(counts, run.info[["Run.ID"]], check_duplicates=check_duplicates, tsne_max_iters=max.iterations, verbose=verbose, is.dist=is.dist,
				tsne_perplexity=run.info[["Perplexity"]], seurat.obj=seurat.obj, barnes.hut=barnes.hut, plot.genes=plot.genes,
				working.dir=run.info[["Working.Dir"]])
		}
	}
	y = proc.time() - ptm
	info(sprintf("tSNE run(s) completed. Runtime: %s s\n ", signif(y[["elapsed"]], 3)))

}

tsne.scanpy <- function()
{

}

tsne.iteration <-function(counts, 
							tsne_run=0, 
							is.dist=F,
							check_duplicates=T,
							working.dir=NULL,
							save.tsne.rotation=T, 
							plot.genes 					=c("Lgr5", "Cell.Cycle_mean_normalised", 
               							"MHC.II_mean_normalised", "Lyz1", "Defa3",
               							"Alpi", "Lct", "Fabp6",
               							"Chgb", "Neurog3",
               							"Fcgbp", "Muc2", "Dclk1", 
               							"Nrep", "Rac2"),
							tsne_perplexity				= 15,
							tsne_max_iters				= 500,
							tsne_whiten					= FALSE,
							clusterer 					= "dbscan", # hierarchical, affinity, dbscan
							cluster  					= FALSE,
							load=F,
							barnes.hut 					= FALSE,
							verbose						= verbose,
							dbscan_eps 					= 6, #12
							dbscan_min_pts 				= 5,
							k 							= 6,
							show_sample_name 			= TRUE,
							border_points				= FALSE,
							label_samples 				= TRUE,
							show_legend 				= FALSE,
							seurat.obj 					= NULL,
							use_rainbow_colours			= TRUE, # otherwise, brewer palette accent is used.
							tsne_min_cost				= 0.05)
{
		#MPI seems to require this. different R workspace on each slave process i guess
		if(load)
		{	
			source("util.r")	
			source("seurat.r")
		}
		init.dir = getwd()
		if(!is.null(working.dir)){
			setwd(working.dir)
		}else{
			working.dir="."
		}
		
		cat(sprintf("\n\n\n"))
		info(sprintf("Starting tSNE run %i -- [Perp=%f, Iters=%f, MinCost=%f, Barnes-Hut=%s, WorkingDir=%s]\n",
			tsne_run, tsne_perplexity, tsne_max_iters, tsne_min_cost, barnes.hut, working.dir))
		#print(head(counts))
		#dist(tsne_input, method="minkowski")
		
		names=colnames(counts);
		if(is.dist){
			info("Converting to distance matrix")
			counts=as.dist(counts)
		}

		if(barnes.hut)
		{
			info("Using Barnes-Hut approximation for faster tSNE")
			#print(corner(counts))
			barnes_hut_tsne = Rtsne(counts, 
								check_duplicates=check_duplicates,
								pca=FALSE, #dont run PCA again
								initial_dims = nrow(counts), 
								perplexity = tsne_perplexity, 
								max_iter = tsne_max_iters, 
								verbose=verbose,
								whiten=tsne_whiten, 
								min_cost=tsne_min_cost)
			tsne.rot = barnes_hut_tsne$Y
		}else
		{
			tsne.rot = tsne(counts, 
							initial_dims = nrow(counts), 
							perplexity = tsne_perplexity, 
							max_iter = tsne_max_iters, 
							whiten=tsne_whiten, 
							min_cost=tsne_min_cost)
		}
		
		print(dim(tsne.rot))
		print(dim(counts))
		if(is.dist)
		{
			rownames(tsne.rot) = names
		}else{
			if(nrow(tsne.rot)==nrow(counts)){rownames(tsne.rot) = rownames(counts)}else
			{
				rownames(tsne.rot) = colnames(counts)
			}
		}
		
		if(save.tsne.rotation)
		{
			rotations.file = paste("tsne_run", tsne_run, "rotations.txt", sep="_")
			#print(rotations.file)
			info(sprintf("Writing tSNE rotations to %s", rotations.file))
			write.table(tsne.rot, file=rotations.file, row.names=T, col.names=T, sep="\t", quote=F)
		}
		
		if(cluster){
			info(sprintf("Clustering method: %s ", clusterer))
			info(sprintf("Clustering .. \n"))

			if(clusterer == "affinity"){
				library(apcluster)
				apresult = apcluster(negDistMat(r=2), tsne.rot)
				#print(apresult)
				# probably a better way of getting the cluster indices out of the apresult object but i 
				# cant see it rn

				examplars = unique(apresult@exemplars)
				sample_to_examplar_map = apresult@idx
				sample_to_cluster_index_map = list()
				for (i in 1:length(sample_to_examplar_map)) {
					sample_to_cluster_index_map[[i]] = which(examplars==sample_to_examplar_map[[i]])
				}

				clusters = unlist(sample_to_cluster_index_map)

			}else if(clusterer=="hierarchical"){
				

				#hierarchical clustering:
				d <- dist(tsne.rot)
				hc <- hclust(d, method="complete")
				clusters <- cutree(hc, k=k)
			}else if(clusterer == "dbscan"){

				#set.seed(665544)
				dbscan_output = dbscan(tsne.rot, eps=dbscan_eps, MinPts=dbscan_min_pts)
				clusters = dbscan_output$cluster
				
			}else{
				cat(sprintf("ERROR:  Unknown clustering method: %s \n", clusterer))	
			}

			cat(sprintf("Found %i clusters:" , length(unique(clusters))))
			print(table(clusters))

			
			output_df = data.frame(tsne.rot)
			colnames(output_df) = c("tsne_1", "tsne_2")

			# trim 'normalised_count' off the sample names
			if(length(grep("expected_count", rownames(output_df))) > 0)
			{
				count_string = "expected_count"
			}else{
				if(length(grep("expected_count", rownames(output_df))) > 0)
				{
					count_string = "normalised_count"
				}else
				{
					cat("WARN: Could not detect a known count_string, cannot trim sample names!\n")
					count_string= NULL
				}
			}

			if(!is.null(count_string))
			{
				x = strsplit(rownames(output_df), count_string) #  paste(count_string,"_"))
				trimmed_sample_names = list()
				for(i in 1:length(x))
				{
					trimmed_sample_names = c(trimmed_sample_names, x[[i]][2])
				}
				trimmed_sample_names = gsub("_", " ", trimmed_sample_names)
			}else
			{
				trimmed_sample_names = rownames(output_df)
			}
			

			output_df <- cbind(sample_name = rownames(output_df), output_df)
			output_df["trimmed_sample_names"] = unlist(trimmed_sample_names)
			rownames(output_df) <- NULL

			output_df["cluster"] = as.factor(clusters)
			# cat(sprintf("Built output data frame. First 6 rows: \n"))
			# print(head(output_df))

			#colour points according to cluster
			if(use_rainbow_colours){
				info(sprintf("Using rainbow colors.. \n"))
				myColors  = rainbow(length(unique(output_df$cluster)))
			}else{
				myColors <- colorRampPalette(brewer.pal(9,"Accent"))(length(levels(output_df$cluster))) 
			}
			names(myColors) <- levels(output_df$cluster)
			

			info(sprintf("Plotting.. \n"))
			cairo_pdf(paste("tsne_iteration_", tsne_run, ".pdf", sep=""), width=8, height=7)
			plot = ggplot(output_df, aes(x=tsne_1, y=tsne_2, colour=cluster, label=trimmed_sample_names))
			if(border_points){
				colScale <- scale_fill_manual(name = paste("cluster-",clusterer,sep=""),values = myColors)
				plot = plot + geom_point(aes(fill = factor(cluster), shape=21, size=3), size=2)  +  scale_shape_identity() + scale_colour_manual(values=rep("black",length(unique(output_df$cluster))))
			}else{
				colScale <- scale_colour_manual(name = paste("cluster-",clusterer,sep=""),values = myColors)
				plot = plot + geom_point(aes(fill=factor(cluster)))
			}
			library(cowplot)
			plot = plot + coord_fixed() + theme_cowplot() + theme(aspect.ratio=1) + colScale
			
			if(label_samples){
					plot = plot + geom_text(data = output_df, vjust=5, size=0.5, colour="black")
			}
			if(!show_legend){
				plot = plot + theme(legend.position="none")
			}
			#print(plot(tsne, col=clusters))
			print(plot)
	        if(!is.null(seurat.obj))
	        {
	               info(sprintf("Drawing feature plot of genes %s", paste(plot.genes, collapse=", ")))
	               s=set.tsne(seurat.obj, data.frame(tsne.rot))
	               rownames(s@tsne.rot) = seurat.obj@cell.names
	               s@data.info$DBclust.ident = clusters
	               feature.plot.scale(s, plot.genes)
	       }else
	       {
	               print("Seurat obj is NULL, cannot draw feature plot")
	       }
	       dev.off() 
		}
		 
       setwd(init.dir)      
       return (TRUE)
}

#### deprecated. use do.tsne with a vector argument for perplexity. [June 2016]
many.tsne <- function(	seurat.obj, 
						weights=NULL,
						data=NULL,
						genes=NULL,
						perp.values=c(5, 6, 7, 8, 9, 10, 12, 15), 
						pc.values=seq(7, 16), 
						iters=10000,
						n.var.genes=0,
						var.genes.pval=0.05,
						use.var.genes.only = T,
						n.cores = 4,
						n.runs=12,
						barnes.hut=T)
{
	
    info("Running multiple tSNE iterations [takes ages obv]..")
    
    #perplexity = 25
    initial_wd = getwd()
    
    # info("Norm counts contains Nans?")
    # print(any(is.na(norm.counts )))

    if(is.null(data)){
    	norm.counts = seurat.obj@data
    	if(is.null(genes))
	    {
	    	if(use.var.genes.only)
		    {
		         info("Using most variable genes only")
		         #var.genes = head(rownames(get.variable.genes( seurat.obj@raw.data)), n.var.genes)#get.var.genes(norm.counts)
		         vg = get.variable.genes(2^norm.counts -1, winsorize = F, fit.all = F)
		         if(n.var.genes >0)
		         {
		         	var.genes = head(rownames(vg), n.var.genes)
		         	info(sprintf("Using %s most variable genes", length(var.genes)))
		         }else
		         {
		         	svg = vg[vg$p.adj<var.genes.pval,] # significantly variable
		         	var.genes = rownames(svg)
		         	info(sprintf("Using %s genes significantly variable (p < %s)", length(var.genes), var.genes.pval))
		         }
		         
		         
		         norm.counts.var.genes = norm.counts[var.genes,]
		         use.counts = norm.counts.var.genes
		    }else
		    {
		         use.counts = norm.counts
		    }
	    }else{
	    	use.counts = na.omit(norm.counts[genes,])
	    }
    }else{
    	use.counts=data
    }
   

    #print(var.genes)
    
    
    info("Size of counts df:")
    print(dim(use.counts))    
    # print(any(is.na(use.counts )))
    # print(head(use.counts [, 1:20]))
    # print(range(use.counts ))

    info("About to start. Seurat obj:")
    print(seurat.obj)       

    if(is.null(pc.values))
    {
    	pc.runs = 1
    }else
    {
    	pc.runs = length(pc.values)
    }

    for(k in 1:pc.runs)
    {
        
        if(!is.null(pc.values)){
        	n.pcs = pc.values[k]
        	info(sprintf("PCs--%i  [%i of %i] \n", n.pcs, k, length(pc.values)))
        	use.data = convert_to_pca_scores(t(use.counts), pcs=1:n.pcs)
    	}else
    	{
    		use.data = use.counts
    	}
        
        for(j in 1:length(perp.values))
        {
            perplexity = perp.values[j]
            info(sprintf("Perplexity--%i  [%i of %i]\n", perplexity, j, length(perp.values)))
            if(!is.null(pc.values))
            {
            	dirname = sprintf("tSNE_PCs_%i_perp_%i", n.pcs, perplexity)
        	}else
        	{
        		dirname = sprintf("tSNE_perp_%i", perplexity)
        	}
            
            dir.create(dirname, showWarnings = FALSE)
            setwd(dirname)
            do.tsne(use.data, 
            			weights=weights,
						n.runs=n.runs, 
            			n.cores=n.cores,
						max.iterations=iters, 
            			perplexity=perplexity, 
            			seurat.obj=seurat.obj, 
						parallel=(n.cores>1),
            			barnes.hut=barnes.hut)
            setwd(initial_wd)
        }
    }
}

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




info("Loading seurat dependencies ")
library(Seurat)
library(shape)
library(ggvis)
library(plyr)
library(NMF) # nmf.options(grid.patch=TRUE
info("Loading seurat functions ")

humpCt=function(x, min=0) {
  return(length(x[x>min]))
}
extract.field=function(string,field=1,delim="_PLATE", fixed=F) return(strsplit(string,delim, fixed=fixed)[[1]][field])

findNGene=function(data,is.expr=1) {
  toRet=unlist(lapply(1:ncol(data),function(x)gtCut(data[,x],cutoff=is.expr)))  
  names(toRet)=colnames(data)
  return(toRet)
}

cleanSeurat <- function(obj, clean.dim.red=T)
{
	cn = colnames(obj@raw.data)
	obj@cell.names = cn
	colnames(obj@data) = cn
	if(clean.dim.red){
		rownames(obj@tsne.rot) = cn
		rownames(obj@pca.rot) = cn
	}
	obj
}

### added this function to convert from Seurat1 objects (now defunct, will not run with Seurat 3 loaded)
convert_to_seurat2 <- function(obj, UMI.data=T, run.dr=T){
  rv<-CreateSeuratObject(counts=obj@raw.data,min.cells = 0, min.genes = 0)
  if(UMI.data){
    ### log2(TPM+1) normalization
    rv=NormalizeData(object = rv,normalization.method = "LogNormalize", scale.factor = 10000,display.progress = TRUE)
  }else{
    rv@data <- log2(rv@raw.data + 1)
  }
  ### scale data
  rv <- ScaleData(object = rv, display.progress = TRUE)
  rv <- FindVariableGenes(object = rv, do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 4, y.cutoff = 0.8,y.high.cutoff = Inf)
  if(run.dr){
  	rv <- RunPCA(object = rv,pcs.compute=50, pc.genes = rv@var.genes, do.print = F)
  	rv <- RunTSNE(object = rv, dims.use = 1:5, do.fast = TRUE,perplexity=20, check_duplicates=F)
  	old_tsne = as.matrix(obj@tsne.rot)
  	rownames(old_tsne) <- rv@cell.names
  	rv@dr$tsne@cell.embeddings <- old_tsne
  }
  rv@meta.data<-obj@data.info
  return(rv)
}

convert_to_seurat3 <- function(obj, UMI.data=T, run.dr=F){
	  rv<-CreateSeuratObject(counts=obj@raw.data,min.cells = 0, min.features = 0)
	  if(UMI.data){
	    ### log2(TPM+1) normalization
	    rv=NormalizeData(object = rv,normalization.method = "LogNormalize", scale.factor = 10000,display.progress = TRUE)
	  }else{
	    rv@data <- log2(rv@raw.data + 1)
	  }
	  	### scale data
	  	rv <- FindVariableFeatures(object = rv, do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 4, y.cutoff = 0.8,y.high.cutoff = Inf)
	    rv <- ScaleData(rv, features=VariableFeatures(rv), display.progress=T)
	  if(run.dr){
		  	rv <- RunPCA(object = rv,pcs.compute=50, pc.genes = rv@var.genes, do.print = F)
		  	rv <- RunTSNE(object = rv, dims.use = 1:5, do.fast = TRUE,perplexity=20, check_duplicates=F)
		  	old_tsne = as.matrix(obj@tsne.rot)
		  	rownames(old_tsne) <- colnames(rv)
		  	rv@dr$tsne@cell.embeddings <- old_tsne
	  }else{
	  		#rv@reductions = obj@dr # fails. can't work out how to extract the dimensionality reductions from Seurat V1 objects.
	  }
	  rv@meta.data<-obj@data.info
	  return(rv)
}

dens.plot <- function(object, feature, reduction.use="tsne",
	interpolate=T,
	pt.size=3,
	pt.cols=NULL,
	dens.cols = jim_special_16,
	cex.all=15, grid.col="grey95", background="white", 
	use.image=T, n=200, alpha=1)
{

	# -- Select the point coordinates
	if (reduction.use=="pca") {
	  dims = c("PC1", "PC2")
	}
	if (reduction.use=="tsne") {
	  dims=c("tSNE_1", "tSNE_2")
	}
	if (reduction.use=="ica") {
	  dims = c("IC1", "IC2")
	}

	library(ggtern)
	x = fetch.data(object, c(dims, feature))
	print(head(x))
	print(dim(x))
	dens <- kde2d.weighted(x[, dims[1]], x[,dims[2]], w=x[,feature], n=n)
	dfdens <- data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z))
	
	r = range(unlist(x[,feature]))
	if(is.null(pt.cols)){
		pt.cols=hmap.cols(r[1], 0.4*r[2], r[2])
	}

	g = ggplot(x, aes_string(x=dims[1], y=dims[2])) 
	if(use.image){
		#image(dens, col=dens.cols)
		info(sprintf("Alpha=%s", alpha))
		g = g + geom_raster(data=dfdens, aes(x=x, y=y, fill=z), alpha=alpha, interpolate = interpolate) 
	}else{
		g = g+ stat_contour(aes(x=x, y=y, z=z, fill=..level.., alpha=..level..), 
		data= dfdens, geom="polygon")
	}
	
	g = g + geom_point(aes_string(colour=feature), size=pt.size) + 
			scale_fill_gradientn(sprintf("%s\n weighted density", feature), colors=dens.cols, guide="colorbar", labels=NULL) +
			guides(alpha=FALSE) + theme_bw() +
			scale_colour_gradientn(sprintf("%s\n expression", feature), colors=pt.cols) + 
			theme(axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"),
			 	axis.line.y = element_line(color="gray15", size = 0.75, lineend="round")) + 
			theme(
				# strip.background = element_blank(), 
		#     	# strip.text.x = element_text(size = cex.strip, face="bold"),
				panel.grid.major = element_line(colour = grid.col),
				panel.grid.minor = element_line(colour = grid.col),
				panel.border= element_blank(),
				panel.background =  element_rect(fill = background),
				panel.border= element_blank(),
				text = element_text(size=cex.all, colour="gray22", face="bold.italic")
			) 
	g
}


### look for a set of features in a seurat object
### if one is missing, print a warning and return a clean set. if not are present,
### return false.
check.data <- function(object, features)
{
    var.options=c("data.info","pca.rot","ica.rot","tsne.rot")    
    clean = vector()

    for (my.var in features) {
      data.use=data.frame()
      if (my.var %in% rownames(object@data)) {
        clean = c(clean, my.var)
      } else {
        for(i in var.options) {
          found = F
          eval(parse(text=paste("data.use = object@",i,sep="")))
          if (my.var %in% colnames(data.use)) {
          	found = T #ok
          	clean = c(clean, my.var)
            break;
          }
        }
        if(!found){warn(sprintf("%s not found!", my.var))}
      }    
 	}
 	if(length(clean)==0)
 	{
 		error("No features found!")
 		return(FALSE)
 	}
    return(clean)
}


subsetSeuratObj <- function(object, cells.use, 
								run.tsne=F, 
								run.tpm=T,
								get.variable.genes=F, 
								remove.nonexpressed=T,
								max.p=0.05)
{
	
	if(!all(rownames(object@data) %in% rownames(object@raw.data)))
	{
		stop("Cannot subset seurat object, @data and @raw.data have different gene names!")
	}
	colnames(object@raw.data) = object@cell.names
	object@raw.data = object@raw.data[, cells.use]
	object@data.info = object@data.info[cells.use,]

	if(remove.nonexpressed){
		object@raw.data = object@raw.data[rowSums(object@raw.data) > 0,]
	}

	if(run.tpm) {
		object@data = data.frame(log2(1+tpm(object@raw.data)), check.names=F)
	}else{
		info("run.tpm is OFF, @data is log2(@raw.data + 1)")
		object@data = data.frame(log2(object@raw.data+1), check.names=F)
	}

	if(!all(dim(object@raw.data) == dim(object@data))) stop ("@raw.data and @data are different size! How? ")
#	object@data = object@data[rownames(object@raw.data), colnames(object@raw.data)]
	info("Row z-scored data in @scale.data")
	object@scale.data = t(scale(t(object@data)))
	object@cell.names = colnames(object@raw.data)
	colnames(object@data) = colnames(object@raw.data)

	info(sprintf("Built new seurat object with %s cells, %s genes" , length(object@cell.names), nrow(object@data)))
	if(get.variable.genes)
	{
		var.genes = get.variable.genes(object@raw.data)
		s@var.genes = var.genes$Gene[var.genes$p.adj < max.p]
		s@var.genes = object@var.genes[!is.na(object@var.genes)]
		info(sprintf("Found %s variable genes (p<%s) for new dataset", length(object@var.genes), max.p))
	}

	if(run.tsne)
	{
		object = Seurat::pca(object, pc.genes=object@var.genes)
		object = run_tsne(object)
	}
	return (object)
}

#override seurat's find_markers, but allow specification of any ID.
get.markers.seurat <- function(object, group.of.interest, group.by=NULL,genes.use=NULL, # group.2=NULL,
								thresh.use=log(2), test.use="bimod")
{
	genes.use=set.ifnull(genes.use,rownames(object@data))
     
    if(group.by %in% colnames(object@data.info))
    {
    	cells.1=object@cell.names[object@data.info[, group.by]==group.of.interest]       
    }else
    {
    	error(sprintf("Given object has no \'%s\' id", group.by))
    	return(FALSE)
    }
    # in case the user passed in cells instead of identity classes
    if (length(group.of.interest>1)&&any(group.of.interest%in%object@cell.names)) {
      cells.1=ainb(group.of.interest,object@cell.names)
    }
    
    cells.2 = object@cell.names
    # if NULL for group.2, use all other cells
    # if (is.null(group.2)) {
    #   cells.2=object@cell.names
    # }
    # else {
    #   cells.2=which.cells(object,group.1)
    # }
    # if (length(group.2>1)&&any(group.2%in%object@cell.names)) {
    #   cells.2=ainb(group.2,object@cell.names)
    # }
    cells.2=anotinb(cells.2,cells.1)
    
    info(sprintf("%s cells in group1. %s cells in group2", length(cells.1), length(cells.2)))
    #error checking
    if (length(cells.1)==0) {
      print(paste("Cell group 1 is empty - no cells with identity class", group.of.interest))
      return(NULL)
    }
    if (length(cells.2)==0) {
      print(paste("Cell group 2 is empty - no cells with identity class", group.2))
      return(NULL)
    }
    
    if (test.use=="bimod") to.return=diffExp.test(object,cells.1,cells.2,genes.use,thresh.use) 
    if (test.use=="roc") to.return=marker.test(object,cells.1,cells.2,genes.use,thresh.use) 
    if (test.use=="t") to.return=diff.t.test(object,cells.1,cells.2,genes.use,thresh.use) 
    if (test.use=="tobit") to.return=tobit.test(object,cells.1,cells.2,genes.use,thresh.use) 
    
    to.return$gene = rownames(to.return)
    return(to.return)
}

#get.markers.all.seurat <- function(object)

display.info <- function(v, dec.places=2)
{
	info(sprintf("Mean: %s, Median: %s, Max: %s, Min: %s, Std.Dev: %s", round(mean(v), dec.places), round(median(v), dec.places), 
		round(max(v), dec.places), round(min(v), dec.places), round(sd(v), dec.places)))
}
# copied seurats filtering. keep cells that have at least min.genes genes expressed
# above threshold is.expr, and keep genes expressed above is.expr in at least min.cells
# cells.
filter.counts <- function(counts, 
						min.genes=2000, 
						min.cells=3, 
						min.trans.mapped=10, #fraction.
						qc=NULL,
						show.n.dropped=10,
						draw.plot=F,
						colour.by=NULL,
						is.expr=1)
{
	
	info("Filter.counts starting")
	info(sprintf("Provided counts table --> %s", paste(dim(counts), collapse=",")))
	info(sprintf("Finding number of expressed genes per cell [expression threshold -- %s count(s)]", is.expr))
	num.genes = colSums(counts>is.expr)
	
	#display the cutoffs
	x = data.frame(sort(num.genes))
	colnames(x) = "nGene"
	x$rank = 1:nrow(x)
	x$cut = x$nGene > min.genes
	n.pass = sum(x$cut)
	ggsave(ggplot(x, aes(x=rank, y=nGene, color=cut)) + geom_point() +scale_colour_manual("Passed filter", values=brewer.pal(3, "Set1")[c(1,3)]) + 
		geom_hline(yintercept=min.genes) + ggtitle(sprintf("Complexity filtering [Min. Genes=%s, %s cells passed]", min.genes,n.pass)) + theme_bw(), filename="Genes_cutoff.pdf", width=9, height=9)

	nas = sum(is.na(names(num.genes)))
	if(nas > 0 )
	{
		warn(sprintf("There are %s NAs in num.genes", nas))
	}
	display.info(num.genes)
	info(sprintf("Number of reads per cell"))
	num.reads = colSums(counts)
	display.info(num.reads)

	nas = sum(is.na(names(num.reads)))
	if(nas > 0 )
	{
		warn(sprintf("There are %s NAs in num.reads", nas))
	}
    vis = data.frame(num.genes, num.reads)
    vis["cell.name"] = colnames(counts)
    # info("Filter.counts using cell names: " )
    # print(head(colnames(counts)))
    if(!is.null(colour.by))
    {
    	vis["batch"] = colour.by
    }
    
    initial.genes = nrow(counts)
    if(!is.null(qc) & min.trans.mapped > 0)
    {
    	before = ncol(counts)
    	info("Mapping rates:")
    	display.info(qc$X.trans_mapped)
    	info(sprintf("Keeping only cells with transcriptome mapping fraction > %s", min.trans.mapped))
    	
    	# info("Using QC sample names:")
	    # print(head(qc$Sample))
	    # print(head(colnames(counts)))
	   	#qc = qc[qc$Sample %in% colnames(counts), ]
	    qc = qc[match(colnames(counts), qc$Sample), ]

	    check = data.frame(as.character(unlist(qc$Sample)), 
	    	colnames(counts), 
	    	as.character(unlist(qc$Sample))==colnames(counts), 
	    	check.names=F)
	    colnames(check) = c("QC.Sample.Name", "Counts.Table.Column", "Match")
	    
	   
	    # return(FALSE)

	    counts.not.qc = colnames(counts)[!(colnames(counts) %in% qc$Sample)]


	    if(length(counts.not.qc) > 0)
	    {
	    	info("QC dataframe size:")
		    info(paste(dim(qc), collapse=", "))
		    info("Counts dataframe size:")
		    info(paste(dim(counts), collapse=", "))
		    
		    matched.ok = colnames(counts)[colnames(counts) %in% qc$Sample]
			error("These samples were found successfully in QC data frame!")
		    info(paste(matched.ok, collapse=", "))

		    error("These samples aren't in the QC data frame!")
		    info(paste(counts.not.qc, collapse=", "))

		    info("Check frame:")
		    print(head(check, n=10))

		    warn("Filtering failed! Not all samples found in QC")
		    #return(FALSE)
	    }
	    

    	cells.dropped = which(qc$X.trans_mapped < min.trans.mapped)
    	counts = counts[, qc$X.trans_mapped > min.trans.mapped]
    	#info(sprintf("After filtering counts table size [%s cells dropped]", cells.dropped))
    	#print(dim(counts))

    	info("Recalculating genes detected and num reads")
    	num.genes = colSums(counts>is.expr)
		num.reads = colSums(counts)
    	
    	vis["trans.mapped.percent"] = qc$X.trans_mapped
    	vis["rRNA"] = qc$X.rRNA
    	info(sprintf("Dropped %s cells with low mapping rates, showing %s of them ", before-ncol(counts), show.n.dropped))
    	
    	show = qc[, c("X.trans_mapped", "Sample")]
    	show$FACS = vis$batch
    	show = show[cells.dropped,]
    	#print(head(show, show.n.dropped))

    	drop.log.file = sprintf("dropped_cells_trans_mapped_under_%s.txt", min.trans.mapped)
    	info(sprintf("Wrote the names of the other dropped cells [mapping rate] to %s", drop.log.file))
    	write.table(show, file=drop.log.file, sep="\t", row.names=F, quote=F)
    	
    	info("Batch composition of cells dropped for low mapping rate:")
    	print(table(as.character(unlist(show$FACS))))

    	info(sprintf("%s cells passed transcriptome mapping rate filtering", ncol(counts)))
    	
    }else
    {
    	cells.dropped = list()
    }


    cells.use = names(num.genes[which(num.genes>min.genes)]) 
    if(!is.null(qc) & min.trans.mapped > 0)
    {
    	info("Transcriptome Mapping rate stats for passing cells:")
		display.info(qc$X.trans_mapped[qc$Sample %in% cells.use])
    }
    

    cells.dropped.complexity = which(num.genes<min.genes)
	cells.dropped = c(cells.dropped, cells.dropped.complexity)

	drop.log.file = sprintf("dropped_cells_ngenes_under_%s.txt", min.genes)
	# print(cells.dropped.complexity)
	# print(head(vis))
	if(length(cells.dropped.complexity) > 0)
	{
		info(sprintf("Wrote the names of %s the dropped cells [complexity] to %s", 
			length(cells.dropped.complexity), drop.log.file))
		write.table(vis[cells.dropped.complexity, ], file=drop.log.file, sep="\t", row.names=F, quote=F)
	}

	vis["dropped"] = 1:nrow(vis) %in% cells.dropped
    counts=counts[,cells.use]
    num.cells=apply(counts,1, humpCt, min=is.expr)
    genes.use=names(num.cells[which(num.cells>min.cells)])
    counts=counts[genes.use,]

    info(sprintf("%s cells passed complexity filtering", length(cells.use)))
	info(sprintf("Complexity of remaining cells: "))
	display.info(findNGene(counts, is.expr))

	info(sprintf("Dropped %s genes not expressed above %s in at least %s cells", 
		initial.genes-nrow(counts), is.expr, min.cells))

	
	info(sprintf("Final counts matrix dimension: %s genes, %s cells", nrow(counts), ncol(counts)))
	
	

    if(draw.plot)
    {
    	info("Plotting")	
    	plot.dir = "QC_filtering"
    	dir.create(plot.dir, showWarnings=FALSE)
    	
    	filter.scatter = paste("filter.mingene.", min.genes, ".min.cells", min.cells,".pdf", sep="")
    	filter.scatter = paste(plot.dir, filter.scatter, sep="/")

    	filter.contour = paste("filter.mingene.", min.genes, ".min.cells", min.cells,"_contour.pdf", sep="")
    	filter.contour = paste(plot.dir, filter.contour, sep="/")

    	filter.contour.rRna = paste("filter.mingene.", min.genes, ".min.cells", min.cells,"_contour_rRNA.pdf", sep="")
    	filter.contour.rRna = paste(plot.dir, filter.contour.rRna, sep="/")
    	
    	filter.genes.hist = paste("filter.gene.hist.", min.genes,".pdf", sep="")
    	filter.genes.hist = paste(plot.dir, filter.genes.hist, sep="/")

    	filter.tmapped.hist = paste("filter.tmapped.hist.", min.trans.mapped,".pdf", sep="")
    	filter.tmapped.hist = paste(plot.dir, filter.tmapped.hist, sep="/")

    	print(filter.scatter)

    	if("trans.mapped.percent" %in% colnames(vis))
    	{
    		if("batch" %in% colnames(vis))
    		{
    			n.cols = length(unique(vis$batch))
    			cols = colorRampPalette(brewer.pal(n.cols, "Set1"))(n.cols) #default.cols(n.cols)

    			# regular scatter plot (labeled with cell name)
    			g = ggplot(vis, aes(x=num.genes, y=trans.mapped.percent, colour=batch, label=cell.name)) + 
	    			geom_hline(yintercept=min.trans.mapped) + geom_vline(xintercept=min.genes) + 
	    			geom_point() + theme_bw() + ggtitle("Single cell QC filtering") + 
	    			geom_text(aes(label=cell.name),hjust=0.5, vjust=2, size=1, colour="black") + scale_colour_manual(values=cols)

	    		# scatter plot showing cell density, no name labels	
    			g.contour = ggplot(vis, aes(x=num.genes, y=trans.mapped.percent)) + geom_vline(xintercept=min.genes) + 
					geom_hline(yintercept = min.trans.mapped) + xlab("Number of genes detected") + ylab("% of reads mapping to transcriptome") + 
					geom_polygon(aes(fill = ..level..), stat = "density2d", alpha = .2) + theme_bw() + 
					scale_fill_continuous("Density of Cells", low = "#56B1F7", high = "#132B43", guide=FALSE ) + 
					geom_point(data=vis, aes(x=num.genes, y=trans.mapped.percent, colour=batch)) + 
					scale_colour_manual("Batch", values=cols) + theme(
						axis.line = element_line(colour = "black"),
					    panel.border = element_blank(),
					    panel.background = element_blank(), 
					    strip.background = element_blank(), 
					    strip.text.x = element_text(size = 16, face="bold"),
					    text = element_text(size=16, colour="gray22"))

				# scatter plot showing cell density, no name labels	, size points by rRNA
    			g.contour.rRna = ggplot(vis, aes(x=num.genes, y=trans.mapped.percent, size=rRNA)) + geom_vline(xintercept=min.genes) + 
					geom_hline(yintercept = min.trans.mapped) + xlab("Number of genes detected") + ylab("% of reads mapping to transcriptome") + 
					geom_polygon(aes(fill = ..level..), stat = "density2d", alpha = .2) + theme_bw() + 
					scale_fill_continuous("Density of Cells", low = "#56B1F7", high = "#132B43", guide=FALSE ) + 
					geom_point(data=vis, aes(x=num.genes, y=trans.mapped.percent, colour=batch)) + 
					scale_colour_manual("Batch", values=cols) + theme(
						axis.line = element_line(colour = "black"),
					    panel.border = element_blank(),
					    panel.background = element_blank(), 
					    strip.background = element_blank(), 
					    strip.text.x = element_text(size = 16, face="bold"),
					    text = element_text(size=16, colour="gray22")) 


    			# note, need seperate density geoms for solid line and transparaent fill, see https://github.com/hadley/ggplot2/issues/1523
    			g.density = ggplot(vis, aes(x=num.genes, fill=batch, colour=batch)) + geom_density(aes(x=num.genes, fill=batch), color = "NA", alpha = 0.3) +  
    				geom_density(aes(x=num.genes, colour=batch), fill = NA) + ggtitle("Filtering by cell complexity") + 
    				theme_bw() + geom_vline(xintercept=min.genes) + xlab("Number of genes detected") + 
    				ylab("Proportion of cells (density estimate)") + scale_fill_manual("Batch", values=cols) +
    				scale_colour_manual("Batch", values=cols, guide=FALSE ) + theme(
    					  axis.line=element_blank(),
					      #axis.text.x=element_blank(),
					      axis.text.y=element_blank(),
					      axis.ticks.y=element_blank(),
					      #axis.title.x=element_blank(),
					      #axis.title.y=element_blank(),
					      #legend.position="none",
					      #panel.background=element_blank(),
					      panel.border=element_blank(),
					      text = element_text(size=16, colour="gray22"))

    			info(sprintf("Drawing vertical line at %s", min.trans.mapped))
				g.density.tmap = ggplot(vis, aes(x=trans.mapped.percent, fill=batch, colour=batch)) + geom_vline(xintercept = min.trans.mapped) + 
					geom_density(aes(x=trans.mapped.percent, fill=batch), color = "NA", alpha = 0.3) +  
    				geom_density(aes(x=trans.mapped.percent, colour=batch), fill = NA) + ggtitle("Filtering by transcriptome mapping rate") + 
    				xlab("% reads mapping to transcriptome") + ylab("Proportion of cells (density estimate)") + theme_bw() + 
    				scale_fill_manual("Batch", values=cols) + scale_colour_manual("Batch", values=cols, guide=FALSE ) + 
    				theme(axis.line=element_blank(),
					      #axis.text.x=element_blank(),
					      axis.text.y=element_blank(),
					      axis.ticks.y=element_blank(),
					      #axis.title.x=element_blank(),
					      #axis.title.y=element_blank(),
					      #legend.position="none",
					      #panel.background=element_blank(),
					      panel.border=element_blank(),
					      text = element_text(size=16, colour="gray22"))
    		}else{
    			g = ggplot(vis, aes(x=num.genes, y=trans.mapped.percent, colour=dropped, label=cell.name)) + 
    			geom_hline(yintercept=min.trans.mapped) + geom_vline(xintercept=min.genes) + 
	    		geom_point() + theme_bw() + ggtitle("Single cell QC filtering") + 
	    		geom_text(aes(label=cell.name),hjust=0.5, vjust=2, size=1, colour="black")
	    		g.density = ggplot(vis, aes(x=num.genes, fill=dropped, colour=dropped)) + geom_density(alpha=0.05) +
    				ggtitle("Filtering by Genes Detected") + theme_bw() + geom_vline(xintercept=min.genes)  
    			g.density.tmap = NULL

    			g.contour = ggplot(vis, aes(x=num.genes, y=trans.mapped.percent)) + geom_point() + geom_vline(xintercept=2000) + 
					geom_hline(yintercept = 35) + xlab("Number of Genes Detected") + ylab("% of reads mapping to transcriptome") + 
					geom_polygon(aes(fill = ..level..), stat = "density2d", alpha = .2) + theme_bw() + 
					scale_fill_continuous("Density Level", low = "#56B1F7", high = "#132B43", guide=FALSE ) 

				g.contour.rRna  = NULL
    		}
		}else{
			g = ggplot(vis, aes(x=num.genes, y=num.reads, colour=dropped)) + geom_point()
			g.density = ggplot(vis, aes(x=num.genes, fill=dropped, colour=dropped)) + geom_density(alpha=0.05) +
    				ggtitle("Filtering by Genes Detected") + theme_bw() + geom_vline(xintercept=min.genes)  
    		g.density.tmap = NULL
    		g.contour=NULL
    		g.contour.rRna = NULL
		}
		ggsave(g, file=filter.scatter, width=12, height=8)
		ggsave(g.density, file=filter.genes.hist, width=12, height=8)
		
		if(!is.null(g.contour))
		{
			ggsave(g.contour, file=filter.contour, width=12, height=8)
		}

		if(!is.null(g.contour.rRna ))
		{
			ggsave(g.contour.rRna, file=filter.contour.rRna, width=12, height=8)
		}
		
		
		if(!is.null(g.density.tmap))
		{
			ggsave(g.density.tmap, file=filter.tmapped.hist, width=12, height=8)
		}
    }

    return (counts)
}


feature.plot.scale <- function(object,
								features.plot,
								raw.data=F, # if set to true, plot the data in @raw.data, usually UMI counts (10X) or RSEM TPM (SS2)
								pdf.output=NULL,
								width=12,
								height=7,
								pc.1=1,
								pc.2=2,
								clip = 0, #if a number between 0-1, large values are clipped to increase colour contrast.
								legend.cex=1,
								show.cluster.title=F,
								borders=F,
								caps=TRUE,
								fix.scale=FALSE, #use same colour scale for all plots
								cells.use=NULL,
								pt.size=1,
								cols.use=NULL,
								cube.helix.cols=FALSE,
								pch.use=16,
								reduction.use="tsne",
								nCol=NULL, 
								round.breaks.to.nearest=100, #only used for values > 20, i.e QC metrics, not log counts
								heat=F,
								palettes=NULL, #optionally colour each Cluster ID with a different brewer palette
								cluster.ids=NULL,
								use.celltype.cols=F, # standardise colour allocation of known celltypes e.g. Goblet, Paneth
								cols.cluster=NULL,
								show.clusters=TRUE,
								color.skew=0.25,
								zero.min=F, 
								show.labels=T) 
{
		if(is.null(cluster.ids)){features.plot=c(features.plot); cluster.ids="ident"}
		object@data.info$ident = object@ident
		features.plot = check.data(object, features.plot)
		
		if(!is.null(palettes))
		{
			if(length(palettes) != length(cluster.ids))
			{
				error(sprintf("Provide the same number of palettes [%s] as cluster IDs [%s]", length(palettes), length(cluster.ids)))
				return (FALSE)
			}
		}
		if(is.null(cells.use))
		{
			cells.use = colnames(object@data)
		}
        dim.code="PC"
 		
 		if(show.clusters){features.plot = c(cluster.ids, features.plot)}
        #if(show.clusters){features.plot = setdiff(c(cluster.ids, features.plot), "ident")}

        if (is.null(nCol)) {
          nCol=2
          if (length(features.plot)==1) nCol=1
          if (length(features.plot)>6) nCol=3
          if (length(features.plot)>9) nCol=4
        }        
        info(sprintf("Using %s columns", nCol)) 
        num.row=floor(length(features.plot)/nCol-1e-5)+1
        par(mfrow=c(num.row, nCol))
        
        # -- Select the point coordinates
        if (reduction.use=="pca") {
          data.plot=object@pca.rot[cells.use,]
          dim.code="PC"
        }
        if (reduction.use=="tsne") {
          data.plot=object@tsne.rot[cells.use,]
          dim.code="tSNE_"
        }
        if (reduction.use=="ica") {
          data.plot=object@ica.rot[cells.use,]
          dim.code="IC"
        }
        if (reduction.use=="diffusionmap"){
        	dc.cols = grep("DC", colnames(object@data.info), value=T)
        	#print(dc.cols)
        	data.plot=object@data.info[cells.use,dc.cols]
          	dim.code="DC"
        }
        
        # -- Get data to plot, axis labels, etc.
        ident.use=as.factor(object@ident[cells.use])
        data.plot$ident=ident.use
        x1=paste(dim.code,pc.1,sep=""); x2=paste(dim.code,pc.2,sep="")
        data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
        data.plot$pt.size=pt.size
        data.use=data.frame(t(fetch.data.updated(object, features.plot, cells.use=cells.use, use.raw=raw.data)))
        feature.count = 1
        info(sprintf("Features to plot: %s", paste(features.plot, collapse=", ")))
        info(sprintf("Cluster IDs: %s", paste(cluster.ids, collapse=", ")))
        use.cols.cluster=cols.cluster
        if(is.null(cols.use)){cols.use = material.heat(20)}
        for(i in features.plot) {
			data.gene=na.omit(data.frame(data.use[i,]))
			if(show.clusters & feature.count<=length(cluster.ids))
			{
				# print("sup")
				cell.ident = cluster.ids[feature.count]
				#print(cell.ident)
				if(cell.ident=="ident")
				{
					clustering = object@ident
					data.gene  = object@ident
				}else{
					clustering = object@data.info[,cluster.ids[feature.count]]
				}
				ncol = length(unique(clustering))
				#info(sprintf("Setting up %s colors to show \'%s\' IDs..", ncol, cell.ident))
				if(is.null(use.cols.cluster))
				{
					if(!is.null(palettes))
					{
						cols.cluster = palettes[[feature.count]]#brewer.pal(ncol, palettes[feature.count])[1:ncol] #distinct.cols(ncol)
					}else
					{
						cols.cluster = default.cols(ncol) #distinct.cols(ncol)
					}
				}else
				{
					# if more than one cluster id is specified, we need to reset the colour
					# palette for each plot, because they might have different numbers of clusters
					if(ncol > length(use.cols.cluster))
					{
						info(sprintf("Not enough colours given for id \'%s\', (%s needed, %s given)", 
							cluster.ids[feature.count], ncol, length(use.cols.cluster)))
						cols.cluster = default.cols(ncol)
					}else
					{
						cols.cluster = use.cols.cluster[1:ncol]
					}
				}
				## If need to match colours with other plots (overrides specified colours)
				if(use.celltype.cols)
				{
					if(!is.null(cols.cluster))
					{
						warn("Overriding provided colours and matching with known cell types")
					}
					info("Using standardised cell-type colours")
					use.cols.cluster = get.celltype.cols(clustering)
					cols.cluster = use.cols.cluster
				}
				order.cells = order(data.gene)
				dp = data.plot[ order.cells, ]
				data.gene = as.factor(unlist(data.gene[order.cells]))
				data.gene = factor(data.gene)
				data.col = as.character(mapvalues(data.gene, unique(data.gene), cols.cluster))
			}else
			{
				info(sprintf("Plotting feature: %s", i)) 
				data.gene = as.numeric(as.matrix(data.gene))
				if(zero.min){min.display=0}else{min.display=min(data.gene)}
				if(clip > 0)
				{
					data.gene = minmax(data.gene, min(data.gene), max(data.gene) * clip)
				}
				if(sum(data.gene != 0) == 0){
					data.cut = rep(1, length(data.gene))
				}else{
					data.cut = as.numeric(as.factor(cut(data.gene, breaks = length(cols.use))))
				}
				if(zero.min){data.range=c(0, max(data.gene))}else{data.range = range(data.gene)}
				if(max(data.gene > 20))
				{
					# not a gene, probably qc metric
					b = (max(data.gene) - min.display) / 10
				}else{
					if(max(data.gene > 5))
					{
						b=2
					}else{
						b=1
					}
				}
				breaks = seq(floor(data.range[1]), floor(max(data.gene)), by=b)
				if(max(data.gene > 20)){
					breaks = round.breaks(breaks, round.breaks.to.nearest)
					is.log2 = FALSE
				}else
				{
					is.log2 = TRUE # assume small values are log2 count
				}
				data.col=cols.use[data.cut]
			}
			# bottom, left, top, right margins in lines of text.
			par("mar"=c(5, 8, 4, 8))
			if(show.clusters & feature.count <= length(cluster.ids))
			{
				if(show.cluster.title)
				{
					title = paste(cluster.ids[feature.count],"IDs")
				}else{
					title = ""
				}

				if(borders)
				{
					plot(dp$x, dp$y,col=data.col,cex=pt.size,pch=pch.use,main=title,xlab=x1,ylab=x2)
				}else{
					plot(dp$x, dp$y,col=data.col,cex=pt.size,pch=pch.use,main=title,xlab=x1,ylab=x2, bty="l")
				}
				par(xpd = TRUE)# allow legend outside plot area
				x=max(dp$x) + 10
				y=max(dp$y) 
				legend(x=x, y=y,
					#inset=c(-0.8,-0.02), 
					cex=legend.cex, 
					legend = unique(data.gene) ,
					border=F, # no border for color cells
					bty = "n" , # no border for legend
					fill = cols.cluster)
				par(xpd=FALSE)
				
			}else
			{
				if(caps)
				{
					title = toupper(i)
				}else
				{
					title = i
				}
				if(borders)
				{
					plot(data.plot$x,data.plot$y, col=data.col, cex=pt.size, pch=pch.use, main=title, xlab=x1, ylab=x2	)
				}else{
					plot(data.plot$x,data.plot$y, col=data.col, cex=pt.size, pch=pch.use, main=title, xlab=x1, ylab=x2, bty="l")
				}
				if(show.labels)
				{
					if(!raw.data)
					{
						title = "Log2(TPM+1)"
					}else
					{
						title = "Log2 UMI count"
					}
					if(is.log2)
					{
						shape::colorlegend(col=cols.use, zlim=data.range, zval=breaks,  left=F, posy=c(0.05,0.85), posx=c(0.85, 0.88), 
							main=title, main.cex=legend.cex, cex=legend.cex*1.5)
					}else
					{
						shape::colorlegend(col=cols.use, zlim=data.range, zval=breaks,  left=F, posy=c(0.05,0.85), posx=c(0.85, 0.88),  cex=legend.cex*1.5)
					}
				}else{
					shape::colorlegend(col=cols.use, zlim=NULL, zval=NULL, posy=c(0,0.8),  cex=legend.cex*1.5) #main="Log2(count)",
				}
			}
			feature.count = feature.count+1
        }
        #par(mfrow=c(1,1))
      
      # if(!is.null(pdf.output))
      # {
      # 		dev.off()
      # 		info("Writing pdf")
      # } 
  }



gene.histogram <- function(seurat.obj, gene, density=F, log=T, by.cluster=F, binwidth=0.25)
{
	d = data.frame(seurat.obj@data[gene,])
	d = t(d)
	d = data.frame(d)
	d["cluster"] = seurat.obj@ident
	if(density)
	{
		g = ggplot(d, aes_string(x=gene)) + geom_density() + xlab("Log2(Counts+1)") 
	}else{
		g = ggplot(d, aes_string(x=gene)) + geom_histogram(binwidth=binwidth, drop=T) + xlab("Log2(Counts+1)") 
	}
	
	if(log &! density)
	{
		g = g +  scale_y_log10() #scale_y_sqrt()
	}

	g = g + theme_bw()

	if(by.cluster)
	{
		g = g + facet_wrap(~cluster)
	}

	return (g)
}

cell.plot <- function(object, cells1, cells2, gene.ids=NULL,col.use="black",nrpoints.use=Inf,pch.use=16,cex.use=0.5,do.ident=FALSE,...) {
            gene.ids=set.ifnull(gene.ids,rownames(object@data))
            if(length(cells1)==1){
            	c1=as.numeric(object@data[gene.ids,cells1])
            	xlab=c1
            }else{
            	c1 = as.numeric(rowMeans(object@data[gene.ids,cells1]))
            	xlab=""
            }

            if(length(cells1)==1){
            	c2=as.numeric(object@data[gene.ids,cell2])
            	ylab=c2
            }else{
            	c2 = as.numeric(rowMeans(object@data[gene.ids,cells2]))
            	ylab=""
            }
            
            gene.cor=round(cor(c1,c2),2)
            smoothScatter(c1,c2,xlab=xlab,ylab=ylab,col=col.use,nrpoints=nrpoints.use,pch=pch.use,cex=cex.use,main=gene.cor)
            if (do.ident) {
              identify(c1,c2,labels = gene.ids)
            }
          }

get.var.genes <- function(counts, take.log=T,
							min.cells = 1, 
							seurat.names.field = 3,
							var.genes.cutoff.x = 1.5,
                            var.genes.cutoff.y = 1.5,
                            min.genes = 1)
{
	info("Setting up Seurat objects.. ")
	if(take.log){
		seurat.obj = new("seurat", raw.data=data.frame(log2(counts+1))) #colcount.data=counts, ident.fxn=getStat3)
	}else{
		seurat.obj = new("seurat", raw.data=data.frame(counts)) 
	}
	info(sprintf("Setting up seurat object. Min cells=%i, min genes=%i ..\n", min.cells, min.genes))
	seurat.obj = Seurat::setup(seurat.obj, 
						names.field = seurat.names.field,
						names.delim = "_",
						project="finding.variable.genes", 
						min.cells = min.cells, 
						min.genes = min.genes, 
						#min.counts=30,
						calc.noise=FALSE, 
						is.expr=0)
	info("Finding variable genes..")
	seurat.obj=mean.var.plot(seurat.obj,y.cutoff = var.genes.cutoff.y, x.low.cutoff = var.genes.cutoff.x,fxn.x = expMean,fxn.y = logVarDivMean)
	cat(sprintf("Found %i variable genes.. \n", length(seurat.obj@var.genes)))
	return (seurat.obj@var.genes)
}

get.pc.genes <- function( seurat.obj, pc.use=1, num.genes=10)
{
	data.pc=seurat.obj@pca.x;
    data.pc=data.pc[order(data.pc[,pc.use]),]
    genes.1=head(rownames(data.pc),num.genes); genes.2=(tail(rownames(data.pc),num.genes))
    genes.use=unique(c(genes.1,genes.2))
    return (genes.use)
}

# see the tutorial at:
# http://www.satijalab.org/clustertutorial1.html
seurat.setup  = function(  counts, 
                            genes.of.interest=c("Lgr5","Lyz1","Scin","Enpep","Hck"), 
                            brewer.palette="OrRd",
                            project.name="Gut.Circuits",
                            run.tsne=FALSE,
                            run.pca = FALSE,
                            do.plot=FALSE,
                            reduce.dims.for.tsne=6,
                            use.ica.for.tsne=FALSE,
                            file.output=FALSE,
                            pca.plots = TRUE,
                            take.log = TRUE,
                            pca.components=50,
                            pca.randomized=T,
                            is.expr = 0.1,
                            min.cells = 5,
                            min.genes = 2000,
                            tsne.rotations=NULL,
                            tsne.iterations=2000,
                            seurat.names.field = 3, # search for groups in the names by splitting on '_', i.e if 1, Lgr5_AS_A77 will get group Lgr5.
                            plot.n.components = 5,
                            get.variable.genes=T,
                            var.genes.min.dispersion = 2,
                            var.genes.max.dispersion = 100,
                            var.genes.min.expression = 0.25,
                            var.genes.max.expression = 100,
                            run.jackstraw=TRUE,
                            draw.genes.of.interest=FALSE,
                            cell_type_label=NULL)
{

	info("Setting up Seurat objects.. ")
	
	# use check.names=F to make sure that the cell names are not altered.
	# if this argument is not there, e.g '+' is replaced by '.' 
	if(take.log){
		info("Taking log2 transform..")
		seurat.obj = new("seurat", raw.data=data.frame(log2(counts+1), check.names=F))#, raw.data=data.frame(counts)) 
	}else
	{
		seurat.obj = new("seurat", raw.data=data.frame(counts, check.names=F)) #, raw.data=data.frame(counts)) 
	}
		
	# print("Raw data column names:")
	# print(head(colnames(seurat.obj@raw.data)))

	cat(sprintf("Seurat filtering. Min cells=%i, min genes=%i ..\n", min.cells, min.genes))
	seurat.obj = Seurat::setup(seurat.obj, 
						names.field = seurat.names.field,
						names.delim = "_",
						project=project.name, 
						min.cells = min.cells, 
						min.genes = min.genes, 
						#min.counts=30,
						calc.noise=FALSE, 
						is.expr=is.expr)
	# print("Cell names:")
	# print(head(seurat.obj@cell.names))
	# print(seurat.obj)

	# info("Plotting two cells.. ")
	# if(do.plot)
	# {
	#  	pdf("cell_plot1.pdf",width=11, height=8.5)
	# }
	# cellPlot(seurat.obj,seurat.obj@cell.names[1],seurat.obj@cell.names[2],do.ident = FALSE)
	# cellPlot(seurat.obj,seurat.obj@cell.names[3],seurat.obj@cell.names[4],do.ident = FALSE)
	# if(do.plot)
	# {
	# 	dev.off()
	# }

	info("Checking for duplicate cell-names")
	if(any(duplicated(seurat.obj@cell.names)))
	{
		warn(sprintf("%s duplicate cell names found. Removing them!", sum(duplicated(seurat.obj@cell.names))))
		seurat.obj = subsetSeuratObj(seurat.obj, cells.use=seurat.obj@cell.names[!duplicated(seurat.obj@cell.names)], run.tsne=F, get.variable.genes=F)
	}

	if(get.variable.genes){
		info("Finding variable genes..")
		if(do.plot)
		{
		    pdf("variable_genes.pdf", width=11, height=8.5)
		}
		seurat.obj=mean.var.plot(seurat.obj,
								y.cutoff = var.genes.min.dispersion,
								y.high.cutoff = var.genes.max.dispersion, 
								x.low.cutoff = var.genes.min.expression,
								x.high.cutoff = var.genes.max.expression,
								fxn.x = expMean,
								fxn.y = logVarDivMean)
		info(sprintf("Found %i variable genes.. ", length(seurat.obj@var.genes)))
	}
	

	if(do.plot)
	{
	    dev.off()
	}
	if(run.pca){
		if(pca.randomized){
			library(rsvd)
			info(sprintf("Doing PCA (Randomized. k=%s)", pca.components)) 
			pca = rpca(t(seurat.obj@data), center=T, scale=T, retx=T, k=pca.components)
			info("PCA Done.")
			seurat.obj@pca.obj = list(sdev=pca$sdev, rotation=pca$rotation, x=pca$x, center=pca$center, scale=pca$scale)
		}else{
			info("Doing PCA (prcomp), center=T, scale=T")
			pca = prcomp(t(seurat.obj@data),center=T,scale=T)
			info("PCA Done.")
			seurat.obj@pca.obj = list(sdev=pca$sdev, rotation=pca$rotation, x=pca$x, center=pca$center, scale=pca$scale)
		}
		# Rahul's code has the labels wrong
		seurat.obj@pca.rot = data.frame(pca$x)
		seurat.obj@pca.x = data.frame(pca$rotation)
		pdf(paste("pc.loaded.genes.pdf"),width=11, height=8.5)
		viz.pca(seurat.obj,1:6)
		dev.off()
	}
  	

	
	# info("Generating pc heatmap..")
	# pdf(paste("pc_heatmap.pdf"),width=11, height=8.5)
	# pcHeatmap(seurat.obj,pc.use = 1,do.balanced = FALSE)
	# dev.off()
	# info("Heatmap done.")

	# if(pca.plots){
	# 	initial_wd = getwd()
	# 	info("Drawing PCA plots")
	# 	dirname = "PCA"
	# 	dir.create(dirname, showWarnings = FALSE)
	# 	setwd(dirname)
	
		
	# 	x = c(1:plot.n.components)
	# 	pairs = combn(x, 2)
	# 	# print(pairs)  
	# 	n_pairs = ncol(pairs)
	# 	cat(sprintf("Generating %i plots for PCs 1-%i.. \n", n_pairs, plot.n.components))
	# 	for(p in 1:n_pairs){
	# 		pair = pairs[, p]
	# 		i = pair[1]
	# 		j = pair[2] 
	# 		if(do.plot)
	# 		{
	# 		    pdf(paste("pc",i,j,".pdf"),width=11, height=8.5)
	# 		}
	# 		pca.plot(seurat.obj,i,j, pt.size = 4)
	# 		if(do.plot)
	# 		{
	# 		    dev.off()
	# 		}
			
	# 	}
	# 	setwd(initial_wd)
	# 	#visualise pcs:
	# 	)
	# 	}
	# }else{
	# 	info("PCA plots is OFF..")
	# }
  
	if(run.jackstraw){
	    info("Running jackstraw..")
		# Do 200 random samplings to find significant genes, each time randomly permute 1% of genes
		# This returns a 'p-value' for each gene in each PC, based on how likely the gene/PC score woud have been observed by chance
		# Note that in this case we get the same result with 200 or 1000 samplings, so we do 200 here for expediency
		seurat.obj=jackStraw(seurat.obj,num.replicate = 200,do.print = TRUE) 

		# The jackStraw plot compares the distribution of P-values for each PC with a uniform distribution (dashed line)
		# 'Significant' PCs will have a strong enrichment of genes with low p-values (solid curve above dashed line)
		# In this case PC1-9 are strongly significant
		if(do.plot){
			pdf("jackStraw.pdf",width=20, height=16)
			print(jackStrawPlot(seurat.obj,PCs = 1:40))
			dev.off()
		}
	}

	# --------- set or calculate tSNE rotations   # ------
	if(!is.null(tsne.rotations))
	{
		info("Reading tSNE rotations from given file:")
		info(sprintf(tsne.rotations))
		rotations = read.delim(tsne.rotations)
		asNumeric <- function(x) as.numeric(levels(x))[x]#as.character(x))
		factorsNumeric <- function(d) modifyList(d, lapply(d[, sapply(d, is.factor)], asNumeric))
		rotations = factorsNumeric(rotations)
		info("Read tSNE rotations. Dimensions: ")
		print(dim(rotations))
		seurat.obj = set.tsne(seurat.obj, rotations)
		info("tSNE rotations set.")
		print(head(seurat.obj@tsne.rot))

	}else{
		if(run.tsne){
			info("Running tSNE..")
		
			if(use.ica.for.tsne)
			{
				reduction = "ica"
				info("Running ICA..")
				seurat.obj = ica(seurat.obj)
			}else
			{
				reduction = "pca"
			}
			info(sprintf("	Using %i %s components..\n", reduce.dims.for.tsne, reduction))
			script("tsne")

			write.table(pca$x[,1:reduce.dims.for.tsne], file="running_tsne_on.txt",  sep="\t", quote=F)

			perp=15
			do.tsne(pca$x[,1:reduce.dims.for.tsne], perplexity=perp, max.iterations=500)
			seurat.obj = set.tsne(seurat.obj, read.delim(sprintf("tSNE_perp_%s/tsne_run_1_rotations.txt", perp)))
			rownames(seurat.obj@tsne.rot) = seurat.obj@cell.names
			colnames(seurat.obj@data) = seurat.obj@cell.names

			# seurat.obj = run_tsne(seurat.obj, 
			# 					reduction.use=reduction, 
			# 					perplexity=15,
			# 					do.fast=T,
			# 					dims.use = 1:reduce.dims.for.tsne, 
			# 					max_iter=500)
			if(do.plot)
			{
			    initial_wd = getwd()
			    dirname = "tSNE"
			    dir.create(dirname, showWarnings = FALSE)
			    setwd(dirname)
				pdf("tsne.pdf",width=11, height=8.5)
			}
			info("Drawing tSNE plot..")
			tsne.plot(seurat.obj, label.pt.size = 2, pt.size = 2)
			if(do.plot)
			{
				dev.off()
			}
		}
		
	}
	#-----------------------------------------------------------------


  	

	if(draw.genes.of.interest){
		info("Drawing violin plots..")
		if(do.plot)
		{
		    pdf(paste("violins.pdf"),width=11, height=8.5)
		}
		vlnPlot(seurat.obj, genes.of.interest)
		if(do.plot)
		{
		    dev.off()
		}

		info("Drawing feature plots..")
		if(do.plot)
		{
		pdf(paste("feature_plots.pdf"),width=11, height=8.5)
		}
		feature.plot.scale(seurat.obj,
						genes.of.interest, 
						pt.size = 1, 
						cols.use=rev(colorRampPalette(brewer.pal(9,brewer.palette))(20)))
		if(do.plot)
		{
		    dev.off()
		}
		setwd(initial_wd)
	}

	info("Seurat setup complete!")
	return(seurat.obj)
}


seurat.cluster = function(seurat.obj, pcs.for.tsne=6, do.plot=T, dbclust.eps=6, find.markers=T)
{
	info(sprintf("Density based clustering based on tSNE map.. [eps=%f]", dbclust.eps))
	#print(seurat.obj)
	# Density cluster the tSNE map - note that the G.use parameter is the density parameter for the clustering - lower G.use to get finer settings
	# Cells which are 'unassigned' are put in cluster 1 - though in this case there are none
	# Assigned cluster will be placed in the 'DBClust.ident' field of seurat.obj@data.info. Putting set.ident=TRUE means that the assigned clusters will also be stored in seurat.obj@ident
	seurat.obj= DBclust_dimension(seurat.obj,1,2, reduction.use = "tsne",G.use = dbclust.eps, set.ident = TRUE)
	if(do.plot)
	{
	  	pdf(paste("tsne_db_clust.pdf"),width=11, height=8.5)
	}
	n.clusters = length(unique(seurat.obj@data.info$DBclust_dimension))
	tsne.plot(seurat.obj, label.pt.size = 2, pt.size = 2)
	if(do.plot)
	{
		  dev.off()
	}
	

	# Build a phylogenetic tree, and rename/reorder cluster names according to their position on the tree
	# See help for details on tree building strategy
	# This gives closely related clusters similar cluster IDs, which is occasionally useful for visualization later on
	# Assigned cluster will be placed in the 'tree.ident' field of seurat.obj@data.info, and also stored in seurat.obj@ident
	
	# info("Building tree based on DB clusters..")
	# pdf("seurat_cluster_tree.pdf")
 #  	seurat.obj = buildClusterTree(seurat.obj, do.reorder = TRUE, reorder.numeric = TRUE,pcs.use = 1:pcs.for.tsne)
 #  	dev.off()
  
	# View the t-SNE plot with the new labels, in a slightly different format
	if(do.plot)
	{
	    pdf(paste("tsne_db_clust_labeled.pdf"), width=11, height=8.5)
	}
	tsne.plot(seurat.obj, do.label = TRUE, label.pt.size = 1, pt.size=1)
	initial_wd = getwd()
	if(do.plot)
	{
    	dev.off()
      	setwd(initial_wd)
	}

	if(find.markers && n.clusters > 1)
	{
		info("Looking for markers..")
		markers.all=find_all_markers(seurat.obj, thresh.test = 3, test.use = "roc", do.print = TRUE)
		print(head(markers.all))
		write.table(markers.all, "seurat_markers.txt", sep="\t")
		markers.use=subset(markers.all, avg_diff>0&power>0.5)$gene
		if(do.plot)
		{
		  	pdf(paste("seurat_markers_heatmap.pdf"),width=11, height=8.5)
		}
		info("Building marker gene heatmap..")
		doHeatMap(seurat.obj,
			genes.use = markers.use,
			slim.col.label = TRUE,
			remove.key = TRUE,
			cexRow=0.2)
		if(do.plot)
		{
			dev.off()
		}

		# draw the feature-heatmap thing.
		# Final visualization! Splits a 'feature plot' into clusters, very useful for seeing lots of info across many clusters
		# Applied here to PC scores. 
		
		info("Drawing feature heatmap..")
		pdf("seurat_feature_heatmap.pdf")
		pcs.plot=paste("PC",1:10,sep="")
		feature.heatmap(seurat.obj,pcs.plot,cols.use = heat.colors(10), pt.size = 2)
		dev.off()
	}
  	
	return (seurat.obj)
}

# wrote this function by looking at the Seurat.R code, and seeing
# how the function 'mclust_dimension' sets the @ident attribute.
set.cluster.labels <- function(seurat.obj, cluster.labels)
{
	to.set=unlist(cluster.labels) #as.numeric(unlist(cluster.labels))
    data.names=names(seurat.obj@ident)
    seurat.obj@data.info[data.names,"m"]=to.set
    seurat.obj@ident=factor(to.set)
    names(seurat.obj@ident)=data.names              
    return (seurat.obj)
}

# replace the tSNE rotations in a seurat object
set.tsne <- function(seurat.obj, rotations)
{
	colnames(rotations) = c("tSNE_1","tSNE_2")
	if(nrow(rotations) == length(seurat.obj@cell.names))
	{
		seurat.obj@tsne.rot = rotations
	}else
	{
		warn("Differing numbers of rows in given tSNE rotations and current object!")
		warn("Trying to match up rownames..")
		print("tSNE rot:")
		print(head(rownames(rotations)))
		print("Cell names:")
		print(head(seurat.obj@cell.names))
		cn = seurat.obj@cell.names
		# a hack. fix later

		if( (length(grep("PLATE", toupper(cn))) > 0) & 
			(length(grep("PLATE", toupper(rownames(rotations)))) == 0)
		   )
		{
			warn("Cells have PLATE info but tSNE does not. Trying to match anyway..")
			
			cn = as.character(unlist(lapply(cn, extract.field)))
			# print(head(cn))
			# print(head(rownames(rotations)))
		}

		indices = match(toupper(cn), toupper(rownames(rotations)))
		indices = na.omit(indices) # ignore tSNE rows that refer to cells we don't have anymore
		
		# ignore cells not mentioned in the tSNE
		keep.cells = seurat.obj@cell.names[seurat.obj@cell.names %in% rownames(rotations)]
		seurat.obj = subsetData(seurat.obj, cells.use= keep.cells)

		# #DEBUG:
		# print(indices)
		# print(length(indices))
		check = data.frame(cbind(rownames(rotations)[indices], seurat.obj@cell.names, cn))
		colnames(check) = c("tSNE-row","seurat-cell", "seurat-cell-corrected")
		print(head(check))
		rots = rotations[indices,]
		rownames(rots) = seurat.obj@cell.names
		seurat.obj@tsne.rot= rots
	}

	return(seurat.obj)
}

# split a cluster (into two) based on the expression 
# of a given gene
split.cluster <- function(seurat.obj, id, gene, frac=0.5)
{
	cat(sprintf("Splitting cluster %s on %s expression..\n", id, gene))
	cat(sprintf("Expression range of %s: ", gene))
	r = range(seurat.obj@data[gene, ])
	cat(r)
	cutoff = as.numeric(quantile(seurat.obj@data[gene, ],prob=1-frac))  # r[1] + frac * (r[2] - r[1])
	cat(sprintf("\n Cutoff-value: %f [Threshold = %f] \n", cutoff, frac))
	
	cells.in.cluster = seurat.obj@data[gene, which(seurat.obj@ident==id)]
	


	#cat(sprintf("Expression of cells: \n"))

	cells.low = which(seurat.obj@ident==id & seurat.obj@data[gene, ] < cutoff)
	
	#print(cells.low)
	cells.high = which(seurat.obj@ident==id & seurat.obj@data[gene, ] >= cutoff)
	cat(sprintf("Split cluster %s [%i cells] into %i %s-low and %i %s-high.. \n", id, length(cells.in.cluster), length(cells.low), gene, length(cells.high), gene))

	label.low = paste(id, gene, "low", sep="-")
	label.high = paste(id, gene, "high", sep="-")
	
	levels(seurat.obj@ident) = c(levels(seurat.obj@ident), label.low, label.high)

	seurat.obj@ident[cells.low] = label.low
	seurat.obj@ident[cells.high] = label.high



	return (seurat.obj)
}

#(re) draw feature plots of top_n markers
draw.feature.plots <- function(seurat.obj, cluster.assignments, top_n=10)
{
	
	cat(sprintf("Adding average expression to %i marker lists..\n", length(unique(cluster.assignments))))
	cluster.names = unique(cluster.assignments)
	n.clusters = length(cluster.names)
	n.types = length(unique(seurat.obj@data.info$orig.ident))

	cols.clusters = sample(colorRampPalette(brewer.pal(9, "Set1"))(n.clusters))
	cols.orig =  colorRampPalette(brewer.pal(9, "Accent"))(n.types)
	for(i in 1:length(cluster.names))
    {
        
        current.cluster = cluster.names[i]
        markers.file = paste(current.cluster, paste("markers_", current.cluster, ".txt", sep=""), sep="/")
        info(sprintf("Reading markers for cluster %i of %i from %s..", i, length(cluster.names), markers.file))
        
        if (isTRUE(file.exists(markers.file))) 
        {
        	markers = read.delim(markers.file)
        	plot_n = top_n
	        
	       	for(k in 1:top_n){
				top_marker = markers$GENE_SYMBOL[k]
				if(!is.na(top_marker)){
					
				}else
				{
					warn(sprintf("Marker %i is Nan. Can only plot the top %i! \n", k, k-1))
					plot_n = k-1
					break
				}
			}	 	
	    
			# draw the plot:
	    	out.file = paste(current.cluster, paste("feature_plot_", current.cluster, "_big", ".pdf", sep=""), sep="/")
			pdf(out.file, width=20, height=20)
			info(sprintf("Writing feature plot [%i markers] to %s", plot_n, out.file))
			feature.plot.scale(seurat.obj, markers$GENE_SYMBOL[1:plot_n], cols.use=rev(get.hmap.col()))
			tsne.plot(seurat.obj, cols.use = cols.clusters)
			s = set.all.ident(seurat.obj, "orig.ident")
			tsne.plot(s, cols.use=cols.orig)
			dev.off()

		}else{
			warn(sprintf("No markers in %s!", current.cluster))
		}

	}
}
#http://stackoverflow.com/questions/25018598/add-a-plot-title-to-ggvis
add_title <- function(vis, ..., x_lab = "X units", title = "Plot Title") 
{
  add_axis(vis, "x", title = x_lab) %>% 
    add_axis("x", orient = "top", ticks = 0, title = title,
             properties = axis_props(
               axis = list(stroke = "white"),
               labels = list(fontSize = 0)
             ), ...)
}






merge.clusters <- function(clustering, clusters.to.merge, new.name=NULL)
{
	if(length(clustering) < 2)
	{
		cat("ERROR: Must provide 2 or more cluster ID's to merge!")
		return (clusterings)
	}

	i = 1
	if(!is.null(new.name)){
		use.id = new.name
		levels(clustering) = c(levels(clustering), use.id)
		clustering[which(clustering == clusters.to.merge[1])] = use.id
	}else
	{
		use.id = clusters.to.merge[1]
	}
	for(id in clusters.to.merge)
	{
		if(i > 1)
		{
			cat(sprintf("Merging cluster %s into %s ..\n", id, use.id))
			clustering[which(clustering == id)] = use.id
		}
		i = i + 1 
		
	} 
	return (factor(clustering))
}

rename.cluster <- function(clustering, clusters.to.rename)
{
	if(length(clustering) < 2)
	{
		cat("ERROR: Must provide 2 or more cluster ID's to rename!")
		return (clusterings)
	}
	id = clusters.to.rename[1]
	new_id = clusters.to.rename[2]
	cat(sprintf("Renaming cluster %s to %s ..\n", id, new_id))
	levels(clustering) = c(levels(clustering), new_id)
	clustering[which(clustering==id)] <- new_id
	return(factor(clustering))
}

# # provide a list of cluster IDs to merge
# # all cells in all clusters will be given
# # the ID of the first cluster in the list.
# merge.clusters <- function(seurat.obj, clusters)
# {
# 	if(length(clusters) < 2)
# 	{
# 		cat("ERROR: Must provide 2 or more cluster ID's to merge!")
# 		return (seurat.obj)
# 	}

# 	i = 1
# 	use.id = clusters[1]
# 	for(id in clusters)
# 	{
# 		if(i > 1)
# 		{
# 			cat(sprintf("Merging cluster %s into %s ..\n", id, use.id))
# 			cat(levels(seurat.obj@ident))
# 			seurat.obj@ident[which(seurat.obj@ident == id)] = use.id
# 		}
# 		i = i + 1 
		
# 	} 

# 	return (seurat.obj)

# }



# given a *NAMED* vector/dataframe of values (must be the same length as
# the number of cells) store it in the seurat object. useful for
# feature plots etc.
set.attribute <- function(seurat.obj, attrib.data)
{
    attrib.name=colnames(attrib.data)
    if(ncol(attrib.data) > 1)
    {
    	 if(all(colnames(attrib.data) %in% colnames(seurat.obj@data.info)))
    	 {
    	 	 seurat.obj@data.info[colnames(attrib.data)] = attrib.data
	 	}else{
	 		 seurat.obj@data.info = cbind(seurat.obj@data.info, attrib.data)
	 	}
    	
	}else{
		 seurat.obj@data.info[attrib.name]=attrib.data     
	}
   	print(colnames(seurat.obj@data.info))
    return (seurat.obj)
}

set.ifnull=function(x,y) {
  if(is.null(x)) x=y
  return(x)
}

# combine two seurat objects
combineData <- function(s1, s2, 
	keep.all=T)   ##if true, genes in one but not the other are included. if false, only genes in both are kept.
	#renormalise=F)  
{

	#script("samples")
	s = s1
   	
   	info(sprintf("Merging S1 [%s] with S2 [%s]", paste(dim(s1@raw.data), collapse=", "), paste(dim(s2@raw.data), collapse=", ")))
    # s1.data = s1@raw.data
    # s1.data$gene = rownames(s1.data)

    # s2.data = s2@raw.data
    # s2.data$gene = rownames(s2.data)

    # set up raw counts by merging
    
    if(!keep.all)
    {
    	#combined.raw=merge(s1.data, s2.data, by="gene", all=keep.all)
    	common.genes = intersect(rownames(s1@raw.data), rownames(s2@raw.data))
    	combined.raw = cbind.data.frame(s1@raw.data[common.genes,], s2@raw.data[common.genes,])
    }else{
    	one.only = setdiff(rownames(s1@raw.data), rownames(s2@raw.data))
    	two.only = setdiff(rownames(s2@raw.data), rownames(s1@raw.data))
    	s1@raw.data[two.only,] <- 0
    	s2@raw.data[one.only,] <- 0
    	common.genes = intersect(rownames(s1@raw.data), rownames(s2@raw.data))
    	if(!(length(common.genes) == nrow(s1@raw.data) & length(common.genes) == nrow(s1@raw.data))) stop("Rows don't match! Error")
    	combined.raw = cbind.data.frame(s1@raw.data[common.genes,], s2@raw.data[common.genes,])
    	

    }
    
    combined.raw[is.na(combined.raw)] <- 0
    info(sprintf("Merged raw counts table has dimensions %s ", paste(dim(combined.raw), collapse=", ")))

    #print(corner(combined.raw))
    #rownames(combined.raw) = combined.raw$gene
    #combined.raw$gene = NULL
    
    n = nrow(combined.raw)
    combined.raw = combined.raw[rowSums(combined.raw) > 0, ]
    info(sprintf("%s genes that are all zero were removed", n - nrow(combined.raw)))

    info(sprintf("There are %s NA values in the new table", sum(is.na(combined.raw))))

    s@raw.data = combined.raw
    s@data = data.frame(log2(tpm(s@raw.data)+1))
 #    # normalise and take log
 #    if(renormalise)
 #    {	
 #    	info("Recomputing Log2TPM")
 #    	s@data = data.frame(log2(normalise.dseq(combined) + 1))
	# }else
	# {
	# 	info("Renormalise is OFF. Merging previously normalised @data")
	# 	# set up raw counts by merging
	#     s1.data = s1@data
	#     s1.data$gene = rownames(s1.data)

	#     s2.data = s2@data
	#     s2.data$gene = rownames(s2.data)
	#     combined=merge(s1.data, s2.data, by="gene", all=keep.all)
	#     combined[is.na(combined)] <- 0
	#     info(sprintf("Merged counts table has dimensions %s ", paste(dim(combined), collapse=", ")))

	#     print(corner(combined))
	#     rownames(combined) = combined$gene
	#     combined$gene = NULL
	#     s@data= combined
	# }
    

    s@cell.names = colnames(s@raw.data)
    colnames(s@data) = s@cell.names
    info("Scaling data")
    s@scale.data=t(scale(t(s@data),center=T,scale=T))
    s@scale.data=s@scale.data[complete.cases(s@scale.data),]
    s@ident=c(s1@ident, s2@ident)

    s@tsne.rot=data.frame(matrix(NA, ncol=2, nrow=ncol(s@data)))
    s@pca.rot=data.frame(matrix(NA, ncol=2, nrow=ncol(s@data)))
    
    s@data.info = rbind.fill(s1@data.info, s2@data.info)
    rownames(s@data.info) = s@cell.names
    # object@gene.scores=data.frame(object@gene.scores[cells.use,]); colnames(object@gene.scores)[1]="nGene"; rownames(object@gene.scores)=colnames(object@data)
    # object@data.info=data.frame(object@data.info[cells.use,])
    # object@mix.probs=data.frame(object@mix.probs[cells.use,]); colnames(object@mix.probs)[1]="nGene"; rownames(object@mix.probs)=colnames(object@data)
    
    return(s)
                 
}


# adapated version of the cor.matrix function in seurat.
 cor.matrix <- function(object, cor.genes=NULL, title="Correlation matrix",
 						cell.inds=NULL, do.k=FALSE,
 						pdf.output=TRUE, pdf.name="cor.matrix.pdf", clusters=NULL,
 						k.seed=1,k.num=4,vis.low=(-1),vis.high=1,vis.one=0.8,pcs.use=1:3,
 						col.use=rev(colorRampPalette(brewer.pal(9,"RdBu"))(64))) {

info("Doing cor.matrix") 	
cor.genes=set.ifnull(cor.genes,object@var.genes)
cell.inds=set.ifnull(cell.inds,colnames(object@data))
cor.genes=cor.genes[cor.genes%in%rownames(object@data)]
data.cor=object@data[cor.genes,cell.inds]
cor.matrix=cor((data.cor))


if(is.null(clusters))
{
	set.seed(k.seed)
	kmeans.cor=kmeans(cor.matrix,k.num)
    cor.matrix=cor.matrix[order(kmeans.cor$cluster),order(kmeans.cor$cluster)]
    kmeans.names=rownames(cor.matrix)
    clusters.use = kmeans.cor$cluster[kmeans.names]
    row.annot=data.frame(cbind(as.factor(clusters.use) ,object@pca.rot[kmeans.names,pcs.use]))
    cluster.col.name = "KMeans.Cluster.ID"
    object@kmeans.cell=list(kmeans.cor)
    ann.colors = list(KMeans.Cluster.ID=default.cols(length(unique(clusters.use))))
}else{
	do.k=FALSE
	cor.matrix=cor.matrix[order(clusters),order(clusters)]
	cluster.names = rownames(cor.matrix)
	clusters.use = clusters[cluster.names]
	row.annot=data.frame(cbind(clusters.use, object@pca.rot[cluster.names,pcs.use]))
	cluster.col.name = "Cluster.ID"
	ann.colors = list(Cluster.ID=default.cols(length(unique(clusters.use))))
}
colnames(row.annot)=c(cluster.col.name,paste("PC",pcs.use,sep=""))


cor.matrix[cor.matrix==1]=vis.one
cor.matrix=minmax(cor.matrix,min = vis.low,max=vis.high)

if(pdf.output)
{
	pdf(pdf.name, width=16, height=12)
}
nmf.options(grid.patch=TRUE) #avoid an extra blank page in the pdf.
if (do.k) 
{
	aheatmap(cor.matrix, col=col.use,Rowv=NA,Colv=NA,annRow=row.annot, main=title, annColors=ann.colors)
}else{
	aheatmap(cor.matrix, col=col.use,annRow=row.annot, main=title, annColors=ann.colors)
}

	if(pdf.output)
	{
		dev.off()
	}

	return(object)
}    



fetch.data.updated <- function(object, vars.all=NULL,cells.use=NULL,use.imputed=FALSE, use.scaled=FALSE, use.raw=FALSE) {
    #print(vars.all)
    cells.use=set.ifnull(cells.use,object@cell.names)
    data.return=data.frame(row.names = cells.use)
    data.expression=object@data; 
    if (use.imputed) data.expression=object@imputed; 
    if (use.scaled) data.expression=object@scale.data;
    if (use.raw) data.expression=log2(object@raw.data+1);
    var.options=c("data.info","pca.rot","ica.rot","tsne.rot")
    data.expression=t(data.expression)
    for (my.var in vars.all) {
		data.use=data.frame()
		if (my.var %in% colnames(data.expression)) {
			data.use=data.expression
		} else {
			for(i in var.options) {
			  	eval(parse(text=paste("data.use = object@",i,sep="")))
			  	if (my.var %in% colnames(data.use)) {
			    	break;
			  	}
			}
		}              
		if (ncol(data.use)==0) {
			print(paste("Error : ", my.var, " not found", sep=""))
			return(0);
		}
		cells.use=ainb(cells.use,rownames(data.use))
		#print(length(cells.use))
		#print(my.var)
		data.add=data.use[cells.use,my.var]
		#print(dim(data.add))
		#print(dim(data.return))
		data.return=cbind(data.return,data.add)  
    }
    colnames(data.return)=vars.all
    rownames(data.return)=cells.use
    return(data.return)
}




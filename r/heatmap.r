library(gplots)

# read as many genes from genes.to.read as necessary 
# to reach n unique genes
read.unique <- function(genes.to.read, existing.genes, n=10)
{
	rval = vector()
 	for(g in genes.to.read)
 	{
 		if(!(g %in% existing.genes))
 		{
 			rval = c(rval, g)
 		}
 		if(length(rval) > n-1 )
 		{
 			return (factor(rval))
 		}
 	}
 	warn(sprintf("Could not find %s unique genes, returning %s", n, length(rval)))
 	length(rval) = n
 	return (factor(rval))
}

MinMax <- function(data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}

### internal function to use heatmap.2 to draw heatmaps with gaps between the groups (colsep)
hmap <- function(seurat.obj, genes, scale=T, disp.min=-2.5, disp.max=2.5, 
	column.labels.angle=45,
	cexRow=1.3, 
	groups.only=NULL, 
	show.mean.expr=T,
	RowSideColors=NULL,
	cols.mean.expr=colorRampPalette(rev(brewer.pal(9, "Greens")))(50), 
	cexCol=NULL, 
	cols=NULL, 
	cols.gene=NULL,
	shuffle.after.trim=F,
	group.by, 
	colsep.col="white", 
	colsep.width=c(1,1), 
	Rowv=NA)
{

	info("Starting hmap function [colsep heatmap.2]")
	info(sprintf("disp.min = %s, data dims=%s, genes to use= %s", disp.min, paste(dim(GetAssayData(seurat.obj)), collapse=","), nrow(genes)))
	info(sprintf("Groups.only -> %s", paste(groups.only, collapse=", ")))
	s = seurat.obj
	#print("GENES:")
	#print(g)
	

	if(show.mean.expr & !is.null(RowSideColors)){
		stop("Heatmap.2 can only draw either 1 set of RowSideColors. So show.mean.expr and show.gene.label cannot both be TRUE if heatmap.2 is used")
	}

	if(scale)
	{
		#rowSide colors showing mean expression:
		if(show.mean.expr){
			cols.use = cols.mean.expr 
			d=Matrix::rowMeans(GetAssayData(s)[g[g%in%rownames(s)],])

			data.cut = as.numeric(as.factor(cut(d, breaks = length(cols.use))))
			breaks = seq(0, floor(max(d)), by=1)
			
			if(show.mean.expr){RowSideColors=rev(cols.use)[data.cut]}else{
				RowSideColors=rev(c("white", "white", "white"))[data.cut]
			}
		}

		
		if(is.null(cols)){cols = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)}
		doHeatMap.seurat(s, genes.use = genes,  scale=scale, disp.min=disp.min, disp.max=disp.max, cexRow=cexRow, cexCol=cexCol,
	          col.use = cols, group.by = group.by, RowSideColors=as.character(unlist(RowSideColors)), 
	          Rowv=Rowv, groups.only=groups.only, shuffle.after.trim=shuffle.after.trim,
	          sepwidth=colsep.width, srtCol=column.labels.angle, remove.key=F, gene.label.cols=cols.gene,
	          sepcolor=colsep.col, density.info="none", margins=c(9,10),
	          keysize=0.8, key.xlab="Normalised Log2-expression", key.title="")
		if(show.mean.expr){
			colorlegend(col=rev(cols.use), zlim=range(d), zval=breaks,  left=T, posy=c(0.45,0.7), posx=c(0.09, 0.11), 
            main="Mean\nLog2-expression", main.cex=0.8, cex=.8)
		}
		
	}else{
		if(disp.min < 0)
		{
			warn("Disp.min should not be negative for unscaled plot")
		}
		if(show.mean.expr){RowSideColors = rep("white", length(g))}
		if(is.null(cols)){cols = colorRampPalette(brewer.pal(9, "GnBu"))(50)}
		doHeatMap.seurat(s, genes.use = genes,  scale=scale, disp.min=disp.min, disp.max=disp.max,
	          col.use = cols, group.by = group.by, cexRow=cexRow, cexCol=cexCol, Rowv=Rowv, 
	          RowSideColors=as.character(unlist(RowSideColors)), shuffle.after.trim=shuffle.after.trim,
	          sepwidth=colsep.width, srtCol=column.labels.angle, remove.key=F, groups.only=groups.only, gene.label.cols=cols.gene,
	          sepcolor="grey20", density.info="none", margins=c(9,10),
	          keysize=0.8, key.xlab="Log2-expression", key.title="")
	}
}

doHeatMap.seurat <-function(object,cells.use=NULL,
							genes.use=NULL, 
							scale=T,
							disp.min=-2.5,
							disp.max=2.5, 
							Rowv=NA,
							groups.only=NULL, 
							shuffle.after.trim=F,
							RowSideColors=NULL, 
							gene.label.cols=NULL,
							do.return=FALSE,
							col.use=pyCols,
							group.by=NULL,
							remove.key=FALSE,...) {
			d = GetAssayData(object)
			
			info(sprintf("Starting heatmap.seurat. Data dims=%s, Disp.min=%s, disp.max=%s", 
				paste(dim(d), collapse=","), disp.min, disp.max))
			if(is.na(Rowv))
			{
				dgram = "none"
			}else{
				info("Using received row dendrogram")
				dgram = "row"
				Rowv = as.dendrogram(Rowv)
				print(Rowv)
			}
			

			
            cells.use=colnames(object)
            genes.found=unname(intersect(genes.use$GENE_SYMBOL,rownames(d)))
            genes.not.found = unname(setdiff(genes.use$GENE_SYMBOL, rownames(d)))
            cells.use=intersect(cells.use,colnames(object))
            cells.ident=Idents(object)[cells.use]
            if (!is.null(group.by)) cells.ident=factor(FetchData(object,group.by)[,1])
            if(is.null(groups.only)){groups.only=levels(cells.ident)}


            cells.ident=factor(cells.ident, levels=groups.only) 
            if(shuffle.after.trim){
        		cells.use=cells.use[order(cells.ident, sample(1:length(cells.ident)))]
        	}else{
    			cells.use=cells.use[order(cells.ident)]
    		}
            
            genes.ident=factor(genes.found) 
            	
            print("B")
            #print(genes.found.indices)
            #print(RowSideColors)
            print(table(RowSideColors))
            print("C")
            if(!is.null(RowSideColors)){
            	print("Recieved row labels")
            	#print(RowSideColors)
            	gene.labels = RowSideColors
            	#print(match(genes.found, genes.use$GENE_SYMBOL))
            	gene.labels = gene.labels[match(genes.found, genes.use$GENE_SYMBOL)] #maintain the order
            	## map the gene labels to their colors
            	RowSideColors = mapvalues(gene.labels, unique(gene.labels), gene.label.cols[1:length(unique(gene.labels))])
            	print("After mapping to colors")
            	print(table(RowSideColors))
            }else{
            	info("No RowSideColors found, generating white labels!")
            	RowSideColors = rep("white", length(genes.found))
            	gene.labels = RowSideColors
            }

            print("D")
            if(is.null(gene.label.cols)){gene.label.cols = distinct.phrogz(length(unique(gene.labels)))}
            
            
            cid = as.character(cells.ident)[order(cells.ident)]
            #print(cid)
            ColSideColors = mapvalues(cid, unique(groups.only), gene.label.cols[1:length(unique(groups.only))])
            # info("Genes found:")
            # print(genes.found)
            # info("Genes NOT found:")
            # print(genes.not.found)

            if(any(!genes.found %in% rownames(d)))
            {
            	error("Somehow these genes not found in table:")
            	print(genes.found[!genes.found %in% rownames(d)])
            	stop("Stopping")
            }

            if(any(!cells.use %in% colnames(d)))
            {
            	error("Somehow these cells not found in table:")
            	print(cells.use[!cells.use %in% colnames(d)])

            	print("colnames are:")
            	print(colnames(d))
            	stop("Stopping")
            }

            data.use=as.matrix(d[genes.found,cells.use])
            info("Data dimensions:")
            print(dim(data.use))
            if(scale){data.use = t(scale(t(data.use)))}
            info("Z-score range:")
            print(range(data.use))

            data.use=MinMax(data.use,min=disp.min,max=disp.max)
            data.use=as.matrix(data.use)
            #print(corner(data.use))

            vline.use=NULL;
            colsep.use=NULL
            rowsep.use=NULL
            hmFunction=heatmap.2
            if (remove.key) hmFunction=heatmap2NoKey

            colsep.use=cumsum(table(cells.ident))
            gene.labels = factor(gene.labels, levels=unique(gene.labels)) ## make sure the order of the levels respects the reordering based on which genes are present.
	        rowsep.use=unname(cumsum(table(gene.labels)))

	        info("Gene labels:")
	        print(table(gene.labels))

	        #print("cumsum")
	        #print(gene.labels)

	        #print(rowsep.use)
			col.lab=rep("",length(cells.use))
			col.lab[round(cumsum(table(cells.ident))-table(cells.ident)/2)+1]=levels(cells.ident)


			info(sprintf("Rendering heatmap.."))
			hmFunction(data.use,Rowv=Rowv,Colv=NA,trace = "none", dendrogram=dgram, col=col.use, colsep = colsep.use, labCol=col.lab, rowsep=rowsep.use, RowSideColors=RowSideColors, ColSideColors=ColSideColors, ...)

			if(!is.null(RowSideColors)){
	        	legend("left",legend=unique(gene.labels), fill=unique(RowSideColors), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
	        }
	        info("Done")
			if (do.return) {
				return(data.use)
			}
  }



average.heatmap <- function(counts, 
							groups, 
							genes, 
							ret.data=F,
							show.text=F,
							pdf.output=NULL, 
							aggregate.by=mean, # median, 
							hclust.method="ward.D", 
							hclust.distance="pearson",
							txt=NULL,
							width=6,
							height=8,
							border=F,
							transpose=F,
							annot.cols=NULL,
							annot.cols.column=NULL,
							annot.column=NULL,
							annot.row=NULL,
							draw.heatmap=T,
							nmf.heatmap=T,
							groups.only=NULL, 
							order.data=NULL, # if null, genes will be either clustered or displayed alphabetically.
							order.by=NULL, # order by expression in a given group (descending)
							cols=NULL,
							color.skew=0.25, # move the midpoint of colour scale for contrast
							disp.min=0, 
							disp.max=10,
							Rowv=NULL, 
							Colv=NULL,
							cexRow=0.7, 
							cexCol=0.3, 
							fontsize=12,
							show.background=F,
							scale = F, 
							balance.colors=T,
							heatmap.title="")
{
	groups = unlist(groups)
	if(!is.null(genes))
	{
		if(!is.null(annot.row)){
			annot.row = data.frame(annot.row)
			cn = colnames(annot.row)
			annot = cbind.data.frame(genes, annot.row)
			annot = annot[!duplicated(genes),]
			genes = as.character(annot$genes)
		}else{
			genes = unique(genes)
		}
		
		keep.genes.in.cts = which(rownames(counts) %in% genes)
		keep.genes.in.list = which(genes %in% rownames(counts))
		
		counts = counts[keep.genes.in.cts, ]

		if(!is.null(annot.row)){
			annot = annot[keep.genes.in.list,]
			genes = as.character(annot$genes)
			print(annot)
			annot.row = data.frame(annot[, cn])
			colnames(annot.row) = cn
			print(annot.row)
		}else{
			genes = genes[keep.genes.in.list]
		}
		
		# reorder 
		counts = counts[genes, ]
		info(sprintf("Shrunk counts table to %i genes ", nrow(counts)))
		info(sprintf("Row annotation has length %s ", length(annot.row)))
	}else
	{
		counts = counts
	}
	if(!is.null(groups.only))
	{
		info("Only keeping samples in these groups: ")
		info(paste(groups.only, collapse=", "))
		keep = which(groups %in% groups.only)
		if(show.background)
		{
			groups[!(groups %in% groups.only)] = "Background"
		}else{
			counts = counts[, keep]
			groups=groups[keep]
		}
	}
		
	# inflate sparse matrices:
	counts <- data.frame(as.matrix(counts))	
	info(sprintf("NCol=%i, Length of groups=%i", ncol(counts), length(groups)))
	cts <- aggregate(t(counts), by=list(groups), FUN=aggregate.by)
	rownames(cts) = cts$Group.1
	cts$Group.1 = NULL
	cts = t(cts)

	rows.non.nan = which(rowSums(is.na(cts)) == 0)
	cts = cts[rows.non.nan,]
	if(!is.null(annot.row)){
		annot.row[rows.non.nan,]	
	}

	if(!draw.heatmap)
	{
		return(cts)
	}

	#draw a heatmap of average cluster expression of markers:
	if(!is.null(pdf.output))
	{
		info("Rendering markers average heatmap to PDF..")
		pdf(pdf.output, width=width, height=height)
	}

	c.a = cts

	if(!is.null(groups.only))
	{
		c.a = c.a[, groups.only]
	}

	if(scale)
	{
		info("Z-score normalising..")
		c.a = t(scale(t(c.a)))
		legend="Z-score"
	}else{
		legend="Log2(Counts+1)"
	}

	# if(is.null(cols))
	# {
	# 	maxv = max(c.a) / color.clip.ratio
	# 	range= range(c.a) / color.clip.ratio
	# 	info("Setting up colours.. ")
	# 	info(sprintf("Colour range: %s", range))
	# 	cols = get.hmap.col(range.val=range, 
	# 						mid.val=color.skew*maxv, 
	# 						steps=c("navyblue", "blue", "cyan", "yellow", "red"), 
	# 						n.steps.between=c(3,6,10,6))
	# }

	if(is.null(cols))
	{
		cols = material.heat(20)
	}
	if(tolower(hclust.distance) %in% c("pearson", "spearman", "kendall")){
		dm = as.dist(1-cor(t(c.a), method=hclust.distance)) 
	}else{
		dm = dist(c.a, method=hclust.distance)
	}

	if(is.null(Rowv))
	{
		Rowv = hclust(dm, method=hclust.method)
	}

	if(is.null(Rowv))
	{
		Rowv = hclust(dm, method=hclust.method)
	}

	c.a[c.a > disp.max] = disp.max
	c.a[c.a < disp.min] = disp.min

	print(range(c.a))

	if(scale & balance.colors){
		info("Balancing colors")
		#print(min(c.a))
		#print(max(c.a))
		breaks = unique(c(seq(min(c.a), 0, length=100), seq(0, max(c.a), length=100)))
		#print(breaks)
	}else{
		breaks = NA
	}

	# if(nmf.heatmap)
	# {
	info(sprintf("Drawing %s heatmap", paste(dim(c.a), collapse="x")))
	#par(mar=c(10.1, 10.1, 10.1, 10.1))
	#nmf.options(grid.patch=TRUE) #avoid an extra blank page in the pdf. (not needed in later NMFversions (>0.23))
	if(transpose)
	{
		c.a = t(c.a)
		#sc = "col"
		tmp = Colv
		Colv=Rowv 
		Rowv=tmp
	}else
	{
		#sc = "row"
	}
	
	###DEBUG:
	# m = data.frame(c.a)
	# m$label = annot.row
	# print(m)

	if(border){
		b = list(cell = TRUE)
	}else{
		b = NA
	}
	if(show.text)
	{
		txt=round(c.a, 2)
	}else{
		txt=NULL
	}
	if(scale)
	{
		aheatmap(c.a, scale="none", txt=txt, color=cols, cexRow=cexRow, 
			cexCol=cexCol, Rowv=Rowv, Colv=Colv, fontsize=fontsize, 
			annCol=annot.column, annRow=annot.row, annColors =annot.cols, 
			border=b, breaks=breaks, main=heatmap.title)
	}else
	{
		aheatmap(c.a, scale="none", txt=txt, color=cols, cexRow=cexRow, 
			cexCol=cexCol, Rowv=Rowv, Colv=Colv, fontsize=fontsize, 
			annCol=annot.column, annRow=annot.row, annColors =annot.cols, 
			border=b, breaks=breaks, main=heatmap.title)
	}
	if(!is.null(pdf.output))
	{
		dev.off()
	}
	if(ret.data){return (c.a)}
}

### smooth data matrix (row-wise) with loess
smooth.data <- function(df, fn="loess", loess.span=10, ma.window=10)
{
	if(!toupper(fn) %in% c("LOESS", "MA", "SPLINE"))
	{
		stop("Smoothing function (fn) must be one of loess, ma (moving average) or spline (not working)")
	}
	sm = function(x){
		dmmy = 1: length(x)
		if(fn=="loess"){
			ls = loess(x ~ dmmy, span=loess.span)
			predict(ls)
		}else{
			if(fn=="spline")
			{
				smooth.spline(x, y=dmmy,control = list(trace = TRUE, tol = 1e-6, low = -1.5))
			}else{
				ma <- function(x,n=ma.window){stats::filter(x,rep(1/n,n), sides=2)}
				unlist(ma(x)[1:length(x)])
			}
		}
		
	}
	r = data.frame(t(apply(df, 1, sm)))
	colnames(r) = colnames(df)
	return(r)
}

#soi - signatures.of.interest
heatmap.signature <- function(seurat.obj,
	genes=NULL,
	sigs=NULL, 
	soi=NULL, 
	ret.data=F,
	order.by=NULL,
	group.by = NULL, 
	groups.only=NULL, # only draw cells from these groups
	second.group.by=NULL,
	cols.group2=NULL,
	cols.mean.expr=NULL,
	show.gene.label=T,
	show.clustering=T,
	show.mean.expr=T,
	show.cell.quality=T,
	show.batch=T,
	show.reads=T,
	k=0,
	h=0,
	colsep=F, # seperate groups
	colsep.width=c(0.1,0.1), 
	colsep.col="white",
	hclust.dynamic.deepSplit=2,
	hclust.dynamic.minclustersize=10,
	hclust.distance="pearson",
	hclust.method="ward.D", #hclust.method = ward.D@2
	cluster.names=NULL,
	var.only=F,
	use.raw=F,
	disp.min=-2.5, 
	disp.max=2.5, 
	disp.max.nGene = 7000,
	disp.min.nGene = 3000,
	do.pam.with.k=NULL, #use PAM clustering with given k
	trim.lists=25, 
	shuffle.after.trim=F,
	filter.genes.min.expr=0, 
	filter.genes.min.cells=0, ## if > 0, will remove genes expressed in less than this many cells
	cexRow=0, 
	cexCol=0, 
	fontsize=12,
	cols.cluster=NULL,
	cols.gene=NULL,
	cols.groups=NULL,
	cols.ngene=brewer.pal(9, "Greys"),
	cols=NULL,
	scale=T,
	center=T,
	pdf.output=NULL,
	width=10,
	height=7,
	Rowv=NULL, #provide dendrogram, F for cluster no dendro, NA for no cluster no dendro
	Colv=NULL) #provide dendrogram, F for cluster no dendro ,NA for no cluster no dendro
{
	if(!is.null(groups.only))
	{
		if(is.null(group.by))
		{
			error("If groups.only is specified, you have to specify a variable to group.by")
			return(FALSE)
		}
		cells.in.clusters = colnames(seurat.obj)[seurat.obj@meta.data[, group.by] %in% groups.only]
        if(length(cells.in.clusters) < length(colnames(seurat.obj)))
        {
        	info("Extracting cells of interest")
	        seurat.obj = subset(seurat.obj, cells=cells.in.clusters)
	        info(sprintf("Reduced to %s cells in %s groups", length(colnames(seurat.obj)), paste(groups.only, collapse=", ")))
	        #print(table(seurat.obj@meta.data[, group.by]))
        }
	}else
	{
		if(is.null(group.by))
		{
			groups.only=unique(unlist(Idents(seurat.obj)))
		}else
		{
			groups.only=unique(unlist(FetchData(seurat.obj, group.by)))
		}
	}

	if(is.null(cols.groups))
	{
		info(sprintf("Setting up %s groups colours", length(groups.only)))
		cols.groups=default.cols(length(groups.only))
	}
	if(use.raw)
	{
		counts = log2(GetAssayData(seurat.obj, slot="counts")+1)
	}else
	{
		counts = GetAssayData(seurat.obj)
	}
	
	if(is.null(soi))
	{
		soi = colnames(sigs)
	}
	
	# print("initial:")
	# print(sigs)

	if(!is.null(sigs) && !is.null(soi))
	{
		if(filter.genes.min.cells>0)
		{
			genes = as.character(melt(sigs, id.vars=NULL)$value) ### this step should put the genes in the same order as the groups.
			genes = na.omit(genes)
			genes = genes[genes != ""]
			info(sprintf("Filtering out genes expressed below %s (log2TPM) in %s cells", filter.genes.min.expr, filter.genes.min.cells))
			gf = rownames(counts)[toupper(rownames(counts)) %in% toupper(genes)]
			if(!is.null(groups.only))
			{
				use.cells = colnames(counts)[unlist(FetchData(seurat.obj, group.by)) %in% groups.only]
			}else{
				use.cells = colnames(counts)
			}
			kg = Matrix::rowSums(counts[gf,use.cells] > filter.genes.min.expr)>filter.genes.min.cells
			expressed.genes = rownames(counts[kg,])
			keep.genes = genes[toupper(genes) %in% toupper(rownames(counts[kg,]))]
			print(keep.genes)
			sigs = prune.df.with.list(sigs, keep.genes)
		}
		if(trim.lists>0)
		{
			sigs = head(sigs, n=trim.lists)
			if(shuffle.after.trim){
				sigs <- sigs[sample(nrow(sigs)),] #### so the best genes aren't at the top of each block.
			}
		}
		#print(head(sigs))

		#print(soi)
		if(length(soi)>1){ ## if soi is a single signature, this reordering destroys the dataframe.
			sigs.reordered = sigs[,soi]
		}else{
			sigs.reordered = data.frame(sigs[,soi])
			colnames(sigs.reordered) = soi
		}
		

		# print(head(sigs.reordered))
		genes = melt(sigs.reordered, id.vars=NULL) ### this step should put the genes in the same order as the groups.
		genes = na.omit(genes)
		colnames(genes) = c("Type", "GENE_SYMBOL")
		#print(head(genes))
		#stop("done")
		#print("after selecting sigs of interest:")
		#print(genes)
	}else
	{
		if(var.only)
		{
			genes = data.frame(seurat.obj@var.genes)
			colnames(genes) = "GENE_SYMBOL"
			genes["Type"] = rep("Variable.Seurat", nrow(genes))
		}else
		{
			genes = data.frame(genes)
			if(ncol(genes)==1)
			{
				colnames(genes) = "GENE_SYMBOL"
			}

		}
	}

	# select the subset of genes.
	genes$GENE_SYMBOL = gsub(" ", "", unlist(as.character(genes$GENE_SYMBOL)), fixed=T)	
	genes = subset(genes, !duplicated(genes$GENE_SYMBOL))
	info(sprintf("Using %s selected genes", nrow(genes)))
	genes$GENE_SYMBOL = gsub(" ", "", genes$GENE_SYMBOL, fixed = TRUE)
	
	### this variable is only used to keep track of which genes are not found.
	original.genes.list = unlist(genes$GENE_SYMBOL)
	
	# keep gene indices, in the given gene list (orig), and in the counts table (cts)
	use.genes.orig = which(toupper(unlist(genes$GENE_SYMBOL)) %in% toupper(rownames(counts)))
	use.genes.cts = which(toupper(rownames(counts)) %in% toupper(unlist(genes$GENE_SYMBOL)))
	genes.list.orig = unlist(genes$GENE_SYMBOL)[use.genes.orig]
	

	# Make sure to use the gene names from the counts table, not the given list, to prevent case or other disagreements.
	genes = genes[use.genes.orig,]
	#genes = genes[match(toupper(genes.list.orig), toupper(rownames(counts)[use.genes.cts])),] ### preserve the original gene order
	if("Type" %in% colnames(genes)) genes = genes[order(genes$Type), ]


	# genes = data.frame(genes.list, stringsAsFactors=F)
	# colnames(genes) = "GENE_SYMBOL"

	d = as.matrix(counts[genes$GENE_SYMBOL, ])
	info(sprintf("Found %s of them in counts table", nrow(d)))
	#print(genes)
	not.found = original.genes.list[!toupper(original.genes.list) %in% toupper(rownames(d))]
	#print(not.found)
	#print(rownames(d))
	#print(genes.list)
	if(length(not.found)>0){
		warn("Could not find these genes:")
		warn(paste(not.found, collapse=", "))
	}

	#print(range(colSums(d)))

	nc.before = ncol(d)
	nr.before = nrow(d)

	keep.genes = which(Matrix::rowSums(d)>0)

	#keep.cells = 1:ncol(d)
	keep.cells = which(Matrix::colSums(d)>0)
	d = d[keep.genes, keep.cells]
	#print(keep.cells)
	if(nc.before!=ncol(d))
	{
		warn(sprintf("Dropped %s samples, they don't express the selected genes", nc.before-ncol(d)))
	}

	if(nr.before!=nrow(d))
	{
		warn(sprintf("Dropped %s genes, they are not present in the selected samples", nr.before-nrow(d)))
	}
	remaining.genes = toupper(rownames(d))
	genes = data.frame(genes)

	matches = match(remaining.genes, toupper(genes$GENE_SYMBOL))
	gene_types = factor(genes$Type[matches])

	# set up gene annotation colours
	gene.color.annotations = data.frame(gene_types)
	colnames(gene.color.annotations) = "Signature Type"
	
	if(show.mean.expr)
	{
		gene.color.annotations$Mean.Expr = rowMeans(d)[keep.genes]
		colnames(gene.color.annotations) = gsub("Mean.Expr", "Mean Expr", colnames(gene.color.annotations))
		if(is.null(cols.mean.expr))
		{
			cols.mean.expr = brewer.pal(9, "Greens")
		}
	}

	if(is.null(cols.gene))
	{
		cols.gene = colorRampPalette(brewer.pal(length(unique(genes$Type)), "Set3"))(length(unique(genes$Type)))
	}
	# set up sample annotation colours
	if(!is.null(second.group.by))
	{
		col.annot = data.frame(seurat.obj@meta.data[,c(second.group.by, group.by)])[ keep.cells,]
		colnames(col.annot) = c(second.group.by, group.by)
		col.annot[,second.group.by]= factor(col.annot[,second.group.by])
		n.groups.2 = length(unique(col.annot[,second.group.by]))
	}else
	{
		keep.cells = unname(unlist(keep.cells))
		if(!is.null(group.by))
		{
			info(sprintf("Grouping cells by %s", group.by))
			col.annot = FetchData(seurat.obj, c(group.by, "nFeature_RNA")) [keep.cells,]
			col.annot[, group.by] = factor(col.annot[, group.by], levels=groups.only)
		}else{
			col.annot = FetchData(seurat.obj, c("orig.ident", "nFeature_RNA"))[keep.cells,]
		}
	}

	if(is.null(cols.group2) & !is.null(second.group.by))
	{
		info(sprintf("Setting up %s colours for second group %s", n.groups.2, second.group.by))
		library(colorspace)
		cols.group2 = rainbow_hcl(n.groups.2, start = 60, end = 240)
	}
	info(sprintf("Using %s gene types", length(unique(gene.color.annotations))))
	
	# do clustering
	if(!is.null(order.by))
	{
		ordby = FetchData(seurat.obj, order.by)[ keep.cells,]
		info(sprintf("Ordering columns by %s", order.by))
		cell.order = order(ordby, decreasing=T)
		d = d[, cell.order]
		Colv = NA
	}else
	{	

		if(is.null(Colv)){
			info("ere")
			if(sum(is.na(d)) > 0)
			{
				error(sprintf("There are %s NAs in the data! Exiting", sum(is.na(d))))
				return(FALSE)
			}

			if(hclust.distance=="pearson")
			{
				dm = as.dist(1-cor(d))
			}else
			{
				if(hclust.distance=="spearman"){
					dm = as.dist(1-cor(d, method="spearman"))
				}else{
					dm = dist(t(d), method=hclust.distance)
				}
				
			}
			#print(range(colSums(d)))
			#print(range(rowSums(d)))
			info("Clustering columns")
			Colv = hclust(dm, method=hclust.method)
		}else{
			if(is.na(Colv)){
				info("Colv set to NA. Not generating column clustering")
			}else{
				stop("Why is Colv neither NULL nor NA")
			}
		}
		
	}

	if(is.null(Rowv))
	{
		if(hclust.distance=="pearson")
		{
			dm.rows = as.dist(1-cor(t(d)))
		}else
		{

			if(hclust.distance=="spearman"){
				dm.rows = as.dist(1-cor(t(d), method="spearman"))
			}else{
				dm.rows = dist(d, method=hclust.distance)
			}
		}
		info("Clustering rows")
		Rowv = hclust(dm.rows, method=hclust.method)
		#print(Rowv)
	}

	if(!is.na(Colv) & is.null(order.by)) {    #& is.null(Colv)){
		info("Now here")
		if(k!=0){
			cl = cutree(Colv, k=k)
		}else
		{
			if(h!=0){
				cl = cutree(Colv, h=h)
				k = length(unique(cl))
			}else{

				info("No k or h provided, using dynamicTreeCut")
				library(dynamicTreeCut)
				cl <- dynamicTreeCut::cutreeDynamic(Colv, 
					distM=as.matrix(dm),
					deepSplit=hclust.dynamic.deepSplit,
					minClusterSize=hclust.dynamic.minclustersize)
				k = length(unique(cl))
			}
		}
		if(!is.null(cluster.names))
		{
			info("Using given cluster names")
			cl = mapvalues(cl, 1:k, cluster.names)
			names(cl) =colnames(seurat.obj)
			#print(table(cl))
		}

		#facs.colors=c("orange","brown","brown2","#4DAF4A","#377EB8","purple")
		if(show.clustering)
		{
			info("Showing the hclust labels")
			col.annot$Cluster = factor(cl)
		}
	}else{
		cl=NULL
	}
	

	b = as.factor(seurat.obj@meta.data$Batch)[keep.cells]
	n.batches = length(levels(b))
	

	if(show.cell.quality)
	{
		col.annot$nFeature_RNA = seurat.obj@meta.data$nFeature_RNA[keep.cells]
		col.annot$nFeature_RNA[col.annot$nFeature_RNA > disp.max.nGene] = disp.max.nGene
		col.annot$nFeature_RNA[col.annot$nFeature_RNA < disp.min.nGene] = disp.min.nGene
		if(show.reads){col.annot$nUMIs.log10 = log10(1+colSums(seurat.obj@raw.data))[keep.cells]}
		if(show.batch){col.annot$Batch = b}
	}

	if(is.null(cols.cluster)){cols.cluster=default.cols(k)}
	ann.colors = list("Signature Type" = cols.gene, 
						"Mean Expr" = cols.mean.expr,
						"Plate"  = brewer.pal(9, "Greys"),
						"Cluster"= cols.cluster,
						"nFeature_RNA"= cols.ngene, 
						"nUMIs.log10"= brewer.pal(9, "RdPu"),
						"Batch"= distinct.cols(n.batches))
	print(ann.colors)
	#print(table(cols.groups))

	if(!is.null(group.by)){
		ann.colors[[group.by]] = cols.groups
	}else{
		#remove a column we added just so the dataframe didn't reduce to a vector
		col.annot$DBclust.ident = NULL 
	}
	#col.annot$nFeature_RNA = NULL
	
	if(!is.null(second.group.by))
	{
		ann.colors[[second.group.by]]=cols.group2
	}

	if(!show.gene.label & !show.mean.expr)
	{
		gene.color.annotations=NULL
	}else{
		if(!show.gene.label)
		{
			gene.color.annotations[, "Signature Type"] <- NULL
		}else{
			info("Show gene label is on!")
			print(table(gene.color.annotations))
		}
	}
	#print(ann.colors)
	if(scale | center)
	{
		info(sprintf("Scaling [center=%s, scale=%s]", center, scale))
		draw.d = t(scale(t(d), center=center, scale=scale))  
		draw.d = MinMax(draw.d,min=disp.min,max=disp.max)
		print(range(draw.d))
		print(cexCol)
	}else{
		draw.d = d
		draw.d = MinMax(draw.d,min=disp.min,max=disp.max)
		# library(rje)
		# cols= cubeHelix(15) #, start=.5, r=-.75)) #
	}
	
	
	if(is.null(cols))
	{
		if(scale)
		{
			cols = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
		}else
		{
			library(colorspace)
			#cols= rev(heat_hcl(12, c = c(80, 30), l = c(30, 90), power = c(1/5, 2))) 
			cols = material.heat(20) #matplotlib.viridis #colorRampPalette(brewer.pal(9, "YlGnBu"))(50)
		}
	}
	info(sprintf("Drawing %s x %s heatmap. Range is (%s, %s)", nrow(draw.d), ncol(draw.d), min(draw.d), max(draw.d)))

	
	if(!is.null(order.by))
	{
		col.annot = col.annot[cell.order,]
	}

	# print(genes)
	# genes = as.character(genes$GENE_SYMBOL) ## updated recently from genes$GENE_SYMBOL
	
	if(is.null(pdf.output))
	{
		if(colsep){
			#print(Rowv)
        
			#heatmap.3(draw.d, trace = "none",col=cols,colsep = colsep.use,labCol=col.lab,cexCol=0.2+1/log10(length(unique(cells.ident))))	
			genes = rownames(draw.d)
			hmap(seurat.obj, genes=genes, scale=scale, group.by=group.by, cols=cols, disp.min=disp.min, disp.max=disp.max, RowSideColors=gene.color.annotations,
					cexRow=cexRow, cexCol=cexCol, Rowv=Rowv, shuffle.after.trim=shuffle.after.trim, groups.only=groups.only, show.mean.expr=show.mean.expr, cols.gene=cols.gene,
					colsep.width=colsep.width, colsep.col=colsep.col)

		}else{
			aheatmap(draw.d, annRow = gene.color.annotations, annColors = ann.colors, annCol=col.annot,
				Rowv=Rowv, Colv=Colv, color=cols, legend = T, cexRow=cexRow, cexCol=cexCol, fontsize=fontsize)
		}
		
	}else
	{
		# nmf.options(grid.patch=TRUE) #avoid an extra blank page in the pdf.  (not needed in later NMFversions (>0.23))
		info(sprintf("Rendering to %s", pdf.output))
		#print(Rowv)
		pdf(pdf.output, width=width, height=height)
		if(colsep){
			hmap(seurat.obj, genes=genes, scale=scale, group.by=group.by, cols=cols, disp.min=disp.min, disp.max=disp.max, RowSideColors=gene.color.annotations,
				Rowv=Rowv, show.mean.expr=show.mean.expr,  cols.gene=cols.gene, shuffle.after.trim=shuffle.after.trim,
				cexRow=cexRow, cexCol=cexCol, colsep.width=colsep.width, colsep.col=colsep.col, groups.only=groups.only)
		}else{
			aheatmap(draw.d, annRow = gene.color.annotations, annColors = ann.colors, annCol=col.annot,
				 Rowv=Rowv, Colv=Colv, color=cols, legend = T, cexRow=cexRow, cexCol=cexCol, fontsize=fontsize)
		}
		
		dev.off()
	}
	if(ret.data){return(list("data"=d, "cluster"=cl, "dendro"=Colv, "draw.data"=draw.d))}
	
}

heatmap.gaps <- function(data, groups)
{
	groups = factor(groups)
	colsep.use=cumsum(table(groups))
      
	col.lab=rep("",nrow(data))
	col.lab[round(cumsum(table(groups))-table(groups)/2)+1]=levels(groups)
	print(colsep.use)
	print(col.lab)
	#heatmap.2(as.matrix(data),Rowv=NA,Colv=NA,trace = "none",colsep = colsep.use,labCol=col.lab,cexCol=0.2+1/log10(length(unique(groups))))

	heatmap.2(as.matrix(data), trace="none", 
		Rowv=NA,
		Colv=NA,
		#hclustfun = function(x)hclust(x, method="ward.D"), 
		#distfun = function(x) as.dist(1-cor(t(x))), 
		col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
		colsep = colsep.use,
		labCol = col.lab,
		cexCol=0.2+1/log10(length(unique(groups))))
}



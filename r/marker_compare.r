
library(Rmisc)
concise.cols.scde = c("GENE_SYMBOL", "FC.Z.Mean", "Z.Lower95", "Z.Upper95", "Z.Stouffer", 
		"p.Max", "p.Fisher", "FC.Lower95", "FC.Upper95", "FC.Mean", "Mean.Log2.Ct", "Mean.Log2.Ct.Bkgrnd",
		"Fraction.Expressing", "Adjusted.Expression")

concise.cols = c("GENE_SYMBOL", "p.Max", "Q.Max",  "p.Fisher", "Q.Fisher", "FC.Min", "FC.Lower95", "FC.Mean", 
			"Mean.Log2.Ct", "Mean.Log2.Ct.Bkgrnd",
		"Fraction.Expressing")


weighted.z <- function(zscores){
	# combine the Z-scores and weight lower ones to penalise genes that
	# are less specific:
	w <- function(x){if(x<10){10-x}else{0.1}}
	if(sum(zscores)==0 | length(zscores)==0 | is.infinite(sum(zscores)))
	{
		return(0)
	}
	weights = unlist(lapply(zscores, w))
	#print(weights)
	(weights %*% zscores) / sqrt(sum(weights^2))
}

# Stouffer's Z-score method
unweighted.z <- function(zscores)
{
	if(sum(zscores)==0 | length(zscores)==0 | is.infinite(sum(zscores)))
	{
		return(0)
	}
	sum(zscores) / sqrt(length(zscores))
}

# Fisher's combined probability test
combine.p.fisher <- function(pvals)
{
	keep <- (pvals > 0) & (pvals <= 1)
    lnp <- log(pvals[keep])
    chisq <- (-2) * sum(lnp)
    df <- 2 * length(lnp)

    if (length(lnp) < 2) {
        #warn("Must have at least two valid p values")
        return(1)
    }
    	
  	pchisq(chisq, df, lower.tail = FALSE)
}

bayesCI <- function(x, R=100)
{
	b = bayesboot(x, weighted.mean, R=R)
	s=summary(b)
	r = c(s[[3]][4], s[[3]][1], s[[3]][3])
	names(r) = c("upper", "mean", "lower")
	return(r)
}

# detailed comparison between cells in a set of n groups,
# calculate pairwise DE genes, and search for the most 
# discriminative up and down-regulated genes between the
# groups. 
# Uses SCDE or permutation/t.test and does pairwise comparisons
# between all groups.
markers.compare <- function(cluster.labels,
							compare.groups=NULL,
							use.scde=F,
							test="wilcox",
							use.knn=F,
							use.k=10,
							use.max.pairs=5000,
							n.threads=4,
							o.ifm=NULL,
							o.prior=NULL,
							use.genes=NULL,
							min.cluster.size = 0,
							precalculate.error.models = T,
							complete_table			= NULL,
							norm.counts 			= NULL,
							raw.counts 				= NULL,
							n_background_samples	= 0,
							sub.sample.groups.to    = 0,
							n_diff_genes			= 0, 
							min.p.val 				= 0.05,
							show_top				= TRUE, 
							seurat.obj 				= NULL,
							batch 					= NULL,
							verbose 				= T,
							annotate.marker.genes   = F)
{

	if(is.null(norm.counts) & is.null(raw.counts))
	{
		stop("Provide some counts!")
	}
	info(sprintf("Starting markers-compare [Test=%s]", test))
	if(is.null(use.genes))
	{
		warn("use.genes is NULL. Testing ALL genes (will impact multiple hypothesis testing by FDR)!")
		if(!check_yes_no()){
			stop()
		}
	}
	
	
	if(is.null(compare.groups))
	{
		compare.groups = factor(unique(cluster.labels))
	}
	info(sprintf("Comparing DE in these groups:"))
	#print(compare.groups)
	for (g in compare.groups){
		if(g %in% unique(cluster.labels))
		{
			cat(sprintf("	%s\n", g))
		}else
		{
			error(sprintf("%s is not one of the given clusters!", g))
			return(FALSE)
		}
	}
	
	if(use.scde){
		names(cluster.labels) = colnames(raw.counts)
		info("Size of input counts matrix:")
		print(dim(raw.counts))
		info("Running SCDE pre-processing")
		if(is.null(o.ifm) | is.null(o.prior))
		{
			if(precalculate.error.models)
			{	
				use.cells = which(cluster.labels %in% compare.groups)
				info(sprintf("Will fit SCDE models for %s cells in %s clusters", length(use.cells), length(compare.groups)))
				
				scde.models = scde.compute.error.models(raw.counts[, use.cells],
						groups=cluster.labels, 
						verbose=verbose, 
						n.cores=n.threads, 
						max.pairs=use.max.pairs,
						knn.models=use.knn, 
						k=use.k)
				
				o.ifm = scde.models[["error.models"]]
				counts_data = scde.models[["counts"]]
				o.prior = scde.models[["expression.prior"]]
				info(sprintf("Built error models for %s cells", ncol(counts_data)))
			}else
			{
				# recalculate the error models for each pairwise comparison. much slower
				o.ifm = NULL
				o.prior = NULL
			}	
		}else
		{
			info("Using provided error models")

		}
		if(!is.null(o.ifm) & !is.null(o.prior))
		{
			before = ncol(raw.counts)
			cluster.labels = cluster.labels[rownames(o.ifm)]
			raw.counts = raw.counts[, rownames(o.ifm)]
			norm.counts = norm.counts[, rownames(o.ifm)]
			info(sprintf("Removing %s cells that we don't have error models for", before-ncol(raw.counts)))
		}
	}else
	{
		names(cluster.labels) = colnames(norm.counts)
	}
	
	if(min.cluster.size > 0){
		info(sprintf("Minimum cluster size set to -> %s ", min.cluster.size))
		keep = table(cluster.labels) > min.cluster.size
		compare.groups = names(table(cluster.labels))[keep]
		dropped.groups = names(table(cluster.labels))[!keep]
		info(sprintf("%s groups dropped: %s", length(dropped.groups), paste(dropped.groups, collapse=", ")))
	}

	comps = combn(compare.groups, 2)
	
	n_compare = ncol(comps)
	info(sprintf("Running %s pairwise comparisons", n_compare))
	rval = list()
	for(i in 1:n_compare)
	{
		pair = unlist(comps[,i])
		info(sprintf("Comparing %s with %s [%s of %s] [using %s threads]", pair[1], pair[2], i, n_compare, n.threads))
		if(use.scde)
		{
			#compare.clusters is the old function, meant for use with SCDE
			de.results = compare.clusters(norm.counts, cluster.labels, pair, o.ifm=o.ifm, o.prior=o.prior,
				raw.counts = raw.counts, sub.sample.groups.to = sub.sample.groups.to, 
				verbose=verbose, min_p_val=min.p.val, n.threads=n.threads)
		}else{
			#compare.groups is newer, intended for permutation, t/wilcoxon or MAST
			de.results = compare.groups(norm.counts, groups=cluster.labels, compare=pair, n.cores=n.threads, test=test, np=1000, subsample=sub.sample.groups.to)
		}
		
		rval[[i]] = de.results
	}
	return (list("pairs"=comps, "diffs"=rval, "error.models"=o.ifm, "expression.prior"=o.prior))
}

# calculate the harmonic mean score (same as F1) of a two 
f1 <- function(x)
{
	1/rowMeans(1/x)
}

gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

gm <- function(x)
{
	apply(x, 1, gm_mean)
}



### make a pretty excel spreadsheet with the signatures, that can be used as a supplementary table.
### provide a summary (all the signatures for all celltypes) and a named list of marker tables with more
### details for each cell type
make.excel.table <- function(summary,
							list.of.marker.tables, 
							working.dir=".", 
							max.sheet.length=1000,
							plain_fonts = T, 
							analysis.name="signatures")
{
	
	library(xlsx)
	list.of.marker.tables = list.of.marker.tables[order(names(list.of.marker.tables))] ## put sheets in alphebetical order based on sheet names (celltypes)
	info(sprintf("Starting 'make.excel.table'. Given list of %s tables", length(list.of.marker.tables)))
	
	excel.dest = paste0(working.dir, "/", sprintf("%s.xlsx", analysis.name))
	warn("Requesting 16GB of RAM to write the excel sheet (xlsx)")
	options(java.parameters = "-Xmx16000m")
	if (file.exists(excel.dest)){file.remove(excel.dest)}

	wb <- createWorkbook(type="xlsx")
	if(plain_fonts){
		TITLE_STYLE     <- CellStyle(wb) 
		SUB_TITLE_STYLE <- CellStyle(wb)
		TABLE_ROWNAMES_STYLE <- CellStyle(wb)
		TABLE_COLNAMES_STYLE <- CellStyle(wb) 
    }else{
    	TITLE_STYLE     <- CellStyle(wb) + Font(wb, heightInPoints=16, isBold=TRUE, underline=1, color = "blue")
		SUB_TITLE_STYLE <- CellStyle(wb) + Font(wb, heightInPoints=14, isItalic=TRUE, isBold=FALSE, color = "blue")
		TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) + Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") + 
								Border(color="black", position=c("RIGHT"), 
								           pen=c("BORDER_THIN")) 
		TABLE_COLNAMES_STYLE <- CellStyle(wb) + #+ Font(wb, isBold=TRUE,heightInPoints=14) +
		    Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
		    Border(color="black", position=c("TOP", "BOTTOM"), 
		           pen=c("BORDER_THIN", "BORDER_THICK")) 
    }
    
    CELL_STYLE <-  CellStyle(wb) + Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER")

	xlsx.addTitle<-function(sheet, rowIndex, title, titleStyle){
		rows <-createRow(sheet,rowIndex=rowIndex)
		sheetTitle <-createCell(rows, colIndex=1)
		setCellValue(sheetTitle[[1,1]], title)
		setCellStyle(sheetTitle[[1,1]], titleStyle)
	}

	### start with a summary of all the markers.
	sheetname = "Summary"
	table = summary
	sheet <- createSheet(wb, sheetName = sheetname)
	xlsx.addTitle(sheet, rowIndex=2, 
	              title=sprintf("Summary of marker genes [%s]", analysis.name),
	              titleStyle = TITLE_STYLE)
	df  = table
	rownames(df) <- NULL
	df = df[, order(colnames(df))]
	colnames(df) = gsub(".", " ", colnames(df), fixed=T)
	n = ncol(df)
	addDataFrame(df, sheet, startRow=6, startColumn=1, 
			 colStyle = CELL_STYLE, row.names=FALSE,
             colnamesStyle = TABLE_COLNAMES_STYLE,
             rownamesStyle = TABLE_ROWNAMES_STYLE)
	setColumnWidth(sheet, colIndex=1:n, colWidth=15)

	jgc <- function() ##java garbage collection
	{
		info("Garbage collecting..")
		gc()
		.jcall("java/lang/System", method = "gc")
	} 

	info(sprintf("Max sheet length set to: %s", max.sheet.length))
	### add all the individual marker lists
	for(table.name in names(list.of.marker.tables))
	{
		
		#jgc() ## try to be memory-efficient
		sheetname = gsub(".", " ", table.name, fixed=T)
		table = list.of.marker.tables[[table.name]]
		sheet <- createSheet(wb, sheetName = sheetname)
		xlsx.addTitle(sheet, rowIndex=2, 
		              title=analysis.name,
		              titleStyle = TITLE_STYLE)
		xlsx.addTitle(sheet, rowIndex=4, 
			      title=sprintf("Marker genes -- %s", table.name),
			      titleStyle = SUB_TITLE_STYLE)
		df  = table[1:max.sheet.length,]
		info(sprintf("Adding data frame (%s) to \'%s\' sheet", paste(dim(df), collapse=" x "), sheetname))
		rownames(df) <- NULL
		
		addDataFrame(df, sheet, startRow=6, startColumn=1, 
				 colStyle = CELL_STYLE, row.names=FALSE,
	             colnamesStyle = TABLE_COLNAMES_STYLE,
	             rownamesStyle = TABLE_ROWNAMES_STYLE)
		setColumnWidth(sheet, colIndex=1:8, colWidth=15)
		setColumnWidth(sheet, colIndex=9:10, colWidth=22)


		#write.xlsx2(, file=, sheetName=sheetname, append=TRUE,row.names=F)
	}
	info(sprintf("Saving to --> %s", excel.dest))
	saveWorkbook(wb, excel.dest)	
}


#test
# use the results of pairwise calculations
# to rank genes, finding the most significantly
# up and down regulated genes for each group
# if need, pass the name of another column name to 
# produce an extra table ranked on that column 
# e.g 'p.Fisher'
markers.rank <- function(compare.obj, counts, clusters, 
							scde=F, 
							max.Q=0.1,
							min.fc=0.25, 
							top.n=500,
							fc.metric="FC.Min", # can also be FC.Lower95
							Q.metric="Q.Max", # can also be Q.Fisher
							fc.col="log2fc",
							min.frac.express=0)
							# extra.rank.on=NULL, 
							# extra.rank.thresh=NULL) #"Log2.Fold.Change.MLE") 
{
	
	if(!fc.metric %in% c("FC.Min", "FC.Lower95"))
	{
		stop("fc.metric must be 'FC.Min' or 'FC.Lower95'")
	}
	if(!Q.metric %in% c("Q.Max", "Q.Fisher"))
	{
		stop("Q.metric must be 'Q.Max' or 'Q.Fisher'")
	}

	info("Starting markers.rank..")
	pairs = compare.obj$pairs
	diffs = compare.obj$diffs
	groups = unique(c(as.character(pairs)))
	n_groups = length(groups)
	summary = data.frame(matrix(NA, nrow = top.n, ncol = n_groups))
	colnames(summary) = groups

	summary.fc.mean = summary
	summary.fc.lower = summary
	summary.fc.min = summary
	#summary.fc.min.filtered.pmax = summary
	summary.filtered = summary
	summary.lb.fc.filtered.Q = summary

	### named list of tables to put into the excel spreadsheet 
	tabs = list()

	for(i in 1:n_groups)
	{
		# pair = unlist(pairs[,i])
		group = groups[i]
		cat("\n\n")
		info(sprintf("Ranking markers for %s [%s of %s]", group, i, n_groups))
		cat("--------------------------------------------------------------\n")
		pairs.up = as.character(unlist(pairs[1, ]))
		pairs.down = as.character(unlist(pairs[2, ]))
		# find the list of comparisons done for this group, keeping track
		# of the direction of fold change
		info(sprintf("Looking for comparisons containing %s", group))
		comparisons.up = which(pairs.up == group) #grep(sprintf("^%s$", as.character(group)), pairs.up) # use the anchors ^$ for an exact grep match
		comparisons.down = which(pairs.down == group) #grep(sprintf("^%s$", as.character(group)), pairs.down)
		comparisons = unlist(c(comparisons.up, -comparisons.down))

		# merge the pairwise comparisons to calculate average ranks 
		merged = NULL
		info(sprintf("Merging %s pairwise comparisons", length(comparisons)))
		print(comparisons)
		start.dir = getwd()
		group.clean = make.names(group)
		dir.create(group.clean)
		info(sprintf("Moving to %s directory", group.clean))
		setwd(group.clean)
		compare.count = 0
		for(j in comparisons)
		{	
			compare.fail = FALSE
			compare.count = compare.count + 1
			if(is.null(diffs[[abs(j)]]))
			{
				compare.fail = T
			}else{
				d.all =  diffs[[abs(j)]]@all
	 			compare.name = diffs[[abs(j)]]@name
	 			info(sprintf("Comparison |%s of %s|: %s", compare.count, length(comparisons), compare.name))
	 			compare.first = strsplit(compare.name, "_vs_")[[1]][1]
	 			compare.second = strsplit(compare.name, "_vs_")[[1]][2]
	 			if(compare.first != group)
	 			{
	 				compare.name = paste(compare.second, compare.first, sep="_vs_")
	 				warn(sprintf("Flipping the fold change, renamed to %s ", compare.name))
	 				d.all <- d.all[order(d.all[, fc.col], decreasing=T), ]
	 				#directional.cols = which(!colnames(d.all) %in% c("p", "p.adj", "GENE_SYMBOL"))
	 				directional.cols = grep("fc", colnames(d.all), ignore.case=T)
		 			d.all[, directional.cols] = d.all[, directional.cols] * -1
	 			}
	 			# # rename the potentially ambiguously named 'p' column
	 			# pval.col = which(colnames(d.all)=="p")
	 			# colnames(d.all)[pval.col] = "p.adj"
	 		}

			if(!compare.fail)
			{
				# keep only relevant columns and add a rank column
				if(scde){
					d.all = d.all[ , c("GENE_SYMBOL", fc.col, "p.adj", "Corrected.Z.score")]
				}else{
				#	info(sprintf("fc.col = %s", fc.col))
				#	d.all = d.all[ , c("GENE_SYMBOL", "p", "p.adj", fc.col, "Lower.bound", "Upper.bound")]
				#}
					use.cols = c("GENE_SYMBOL", "p", "p.adj", fc.col)
					if("Lower.bound" %in% colnames(d.all) & "Upper.bound" %in% colnames(d.all))
					{
						use.cols = c(use.cols, "Lower.bound", "Upper.bound")
					}
					d.all = d.all[ , use.cols]
				}				
				#deal with '/s in names'
				compare.name = gsub("/", ".", compare.name, fixed=T)
				write.table(d.all, file=paste("all_", compare.name, ".txt", sep=""), sep="\t", row.names=F, quote=F)
				if(is.null(merged))
				{
					merged = d.all
					first.compare.name = compare.name 
				}else
				{
					#info(sprintf("Merging %s comparisons", compare.name))
					merged = merge(merged, d.all, by="GENE_SYMBOL", all=T, suffixes=c("",paste("_",compare.name, sep="")))
				}
			}else{
				warn("Comparison failed. Skipping")
			}
		}

		merged = do.rank(merged, group, clusters, counts, fc.col=fc.col, do.fast=T)


		# the first comparison wasn't merged, so it doesn't have a suffix. add one
		info("Adding suffix to first comparison")
		unlabeled = c(fc.col, "p.adj", "Corrected.Z.score")
		unlabeled.indices = which(colnames(merged) == fc.col | 
		colnames(merged) == "p.adj" | colnames(merged) == "Corrected.Z.score")
		colnames(merged)[unlabeled.indices] = paste(unlabeled, first.compare.name, sep="_")

		markers.table = paste("markers_", group.clean, "_full.txt", sep="")
		markers.table.small = paste("markers_", group.clean, ".txt", sep="")
		info("Writing tables")
		# merged[unranked, ] = merged[unranked[order(merged[unranked, "Average.Rank"])], ]
		write.table(merged, file=markers.table, sep="\t", quote=F, row.names=F)
		if(scde){
			write.table(merged[, concise.cols.scde], file=markers.table.small, sep="\t", quote=F, row.names=F)
		}else{
			write.table(merged[, concise.cols], file=markers.table.small, sep="\t", quote=F, row.names=F)
		}

		
		
		if(scde){
			merged = merged[order(merged$Z.Stouffer, decreasing=T), ]
		}
		# print(dim(summary))
		# print(dim(merged))
		if(nrow(summary) > nrow(merged))
		{
			summary = summary[1:nrow(merged),]
			summary.fc.mean = summary.fc.mean[1:nrow(merged),]
			summary.fc.lower = summary.fc.mean[1:nrow(merged),]

			top.n = nrow(merged)
		}

		info(sprintf("Adding %s column to summary [FC.Mean]", group))
		summary[, group] = head(merged$GENE_SYMBOL, n=top.n)
		merged = merged[order(merged$FC.Mean, decreasing=T), ]
		
		info(sprintf("Adding %s column to summary [FC.Lower95]", group))
		summary.fc.mean[, group] = head(merged$GENE_SYMBOL, n=top.n)
		merged = merged[order(merged$FC.Lower95, decreasing=T), ]
		summary.fc.lower[, group] = head(merged$GENE_SYMBOL, n=top.n)

		info(sprintf("Adding %s column to summary [FC.Min]", group))
		###SORT 
		info(sprintf("SORTING ON --->  %s", fc.metric))
		if(!fc.metric %in% colnames(merged))
		{
			print(colnames(merged))
			stop(sprintf("%s not a valid column to sort on!", fc.metric))
		}
		merged = merged[order(merged[, fc.metric], decreasing=T), ]
		
		ranked.genes = head(merged$GENE_SYMBOL, n=top.n)
		length(ranked.genes) = nrow(summary.fc.min)
		summary.fc.min[, group] = ranked.genes

		## store the ranked (by minFC) in a list so it can be put in the excel sheet
		tabs[[group]] = merged[, concise.cols]

		#========================================================================
		### Filter genes on maximum p-value, and must have positive FC, to 
		### define signature genes for a given cell-type
		# info(sprintf("Filtering [max.p.max=%s, min.fc.min=%s]", max.p.max, min.fc.min))
		# merged.filtered = merged[merged$p.Max < max.p.max, ]
		# merged.filtered = merged.filtered[merged.filtered$FC.Min > min.fc.min, ]
		
		# info(sprintf("%s filtered [p.max < %s, min.fc <%s] marker genes for %s-cells:", nrow(merged.filtered), max.p.max, min.fc.min, group))
		# write.table(merged.filtered, file=sprintf("markers_%s_filtered_pmax_%s_minfc_%s.txt", group.clean, max.p.max, min.fc.min), sep="\t", quote=F, row.names=F)

		# filtered.genes = head(merged.filtered$GENE_SYMBOL, n=write.top.summary)
		# length(filtered.genes) = nrow(summary.fc.min.filtered.pmax)
		# summary.fc.min.filtered.pmax[, group] = filtered.genes
		#========================================================================


		#========================================================================
		### Filter genes on maximum (FDR) Q-value, with positive FC, to define
		### more stringent signatures. This should replace pmax.
	
		info(sprintf("FC metric: %s", fc.metric))
		info(sprintf("Q metric: %s", Q.metric))

		info(sprintf("Filtering on %s < %s", Q.metric, max.Q))
		if(!Q.metric %in% colnames(merged))
		{
			print(colnames(merged))
			stop(sprintf("%s not a valid column to sort on!", Q.metric))
		}
		info(sprintf("Filtering on %s > %s", fc.metric, min.fc))
		merged.filtered = merged[merged[, Q.metric] < max.Q, ]
		info(sprintf("%s genes passed (%s removed)", nrow(merged.filtered), nrow(merged)-nrow(merged.filtered)))
		merged.filtered = merged.filtered[merged.filtered[, fc.metric] > min.fc, ]
		info(sprintf("%s filtered [%s < %s, %s > %s] marker genes for %s-cells:", nrow(merged.filtered), Q.metric, max.Q, fc.metric, min.fc, group))
		
		merged.filtered = merged.filtered[merged.filtered$Fraction.Expressing > min.frac.express, ]
		write.table(merged.filtered, file=sprintf("markers_%s_filtered_%s_%s_%s_%s_%s_frac.exprs.txt", 
			group.clean, Q.metric, max.Q, fc.metric, min.fc, min.frac.express), sep="\t", quote=F, row.names=F)
		
		
		filtered.genes = head(merged.filtered$GENE_SYMBOL, n=top.n)
		length(filtered.genes) = nrow(summary.filtered)
		summary.filtered[, group] = filtered.genes
		
		
		cat("\n\n")
		cat("\n\n")

		setwd(start.dir)
	}
	if(scde){
		write.table(summary, file=sprintf("top_%s_markers_Z_Stouffer.txt", top.n), sep="\t", row.names=F, quote=F)
		write.table(summary.fc.mean, file=sprintf("top_%s_Mean_FC.txt", top.n), sep="\t", row.names=F, quote=F)
	}else{
		write.table(summary.fc.min, file=sprintf("top_%s_%s.txt", top.n, fc.metric), sep="\t", row.names=F, quote=F)
		write.table(summary.filtered, file=sprintf("top_%s_filtered_%s_%s_%s_%s_frac.expr_%s.txt", top.n, fc.metric, min.fc, Q.metric, max.Q, min.frac.express), sep="\t", row.names=F, quote=F)
	}
	
	write.table(summary.fc.lower, file=sprintf("top_%s_Lower95_FC.txt", top.n), sep="\t", row.names=F, quote=F)

	# expl = data.frame(concise.cols)
	# expl$Explanation = c("-", 
	# 	"Mean of Z StN and FC StN", 
	# 	"Fold change (FC) signal-to-noise (StN), mean/std of Conservative Estimate for FC",
	# 	"Z-score signal-to-noise (StN), mean/std of corrected Z-score",
	# 	"Aggregate corrected Z-score (Stouffer's Z-score)", 
	# 	"Minimum corrected Z-score", 
	# 	"Maximum p value",
	# 	"Mean log2 count (expression level) in the group",
	# 	"Mean log2 count (expression level) in the background",
	# 	"Fraction of cells in the group expressing the gene", 
	# 	"Mean log2 count (expression level) in the group divided by fraction expressing")
	# colnames(expl) = c("Column/Metric", "Explanation")
	# write.table(expl, file="explanation.txt", quote=F, row.names=F)

	
	info(sprintf("Completed ranking markers for %s groups", n_groups))
	return(list("summary.filtered"=summary.filtered, "tabs"=tabs, "Q.metric"=Q.metric, 
		"max.Q"=max.Q, "fc.metric"=fc.metric, "min.fc"=min.fc, "min.frac.express"=min.frac.express))
}


do.marker.plots <- function(seurat.obj, clusters, n.markers = 20)
{	
	init.dir = getwd()
	for(i in 1: length(unique(clusters)))
    {
        cluster.id = as.character(unique(clusters)[i])
        
        cat("\n \n")
        cat(sprintf("Examining \'%s\', cluster %i of %i..\n", cluster.id, i, length(unique(clusters))))
        setwd(cluster.id)
        info(sprintf("Moving to %s", getwd()))
        
        markers.file = sprintf("markers_%s.txt", cluster.id)
        if(file.exists(markers.file))
		{
			m = read.delim(markers.file, check.names=F)
			rownames(m) <- m$GENE_SYMBOL
			plot.markers(seurat.obj=seurat.obj, de.sig=m, clusters=clusters, secondary.id="Batch")         
            if(nrow(m) < n.markers)
            {
                n.markers = nrow(m)
            }
            av = average.heatmap(seurat.obj@data, 
                          clusters, genes=unique(head(m$GENE_SYMBOL, n=n.markers)), scale=T,  
                          pdf.output=TRUE, pdf.name=paste(cluster.id,".Top", n.markers,"Markers.Z.Score.pdf", sep=""))
            av = average.heatmap(seurat.obj@data, 
                          clusters, genes=unique(head(m$GENE_SYMBOL, n=n.markers)), scale=F,  
                          pdf.output=TRUE, pdf.name=paste(cluster.id,".Top", n.markers,"Markers.pdf", sep=""))           
		}else{
          	warn(sprintf("No markers detected for %s", cluster.id))
        }
        setwd(init.dir)
    }
}


do.rank <- function(merged, group, clusters, counts, fc.col, expression.thresh=0, scde=F, do.fast=F)
{
	
	# add info on expression level 
	# -----------------------------------
	
	info(sprintf("Ranking %s genes to find markers for %s", nrow(counts), group))
	cells = which(clusters==group)
	
	if(length(cells)==0){
		stop("No cells in this group !? Exiting")
	}
	not.cells = which(clusters!=group)
	group.counts = counts[,cells]
	bkgd.counts = counts[, not.cells]
	is.expr = 0
	#print(cells)
	info("Computing mean expression in background")
	Mean.Log2.Ct.Bkgrnd = Matrix::rowMeans(bkgd.counts)
	info(sprintf("Computing mean expression in %s", group))
	Mean.Log2.Ct = Matrix::rowMeans(group.counts)
	Fraction.Expressing = Matrix::rowSums(group.counts > is.expr) / length(cells)
	expression.info = data.frame(Mean.Log2.Ct, Mean.Log2.Ct.Bkgrnd, Fraction.Expressing)
	expression.info$GENE_SYMBOL = rownames(expression.info)
	rownames(expression.info) = NULL
	merged = merge(merged, expression.info, by="GENE_SYMBOL")
	before = nrow(merged)
	merged = merged[merged$Mean.Log2.Ct>expression.thresh,]
	info(sprintf("Removed %s genes that are not expressed above %s log2(counts) in %s", before-nrow(merged), expression.thresh, group))

	info("Ranking genes based on merged comparisons")
	#rank.cols = grep("Rank", colnames(merged))
	fc.cols = grep(fc.col, colnames(merged))
	pval.cols = c(which(colnames(merged)=="p"), grep("p_", colnames(merged)))
	Qval.cols = grep("p.adj", colnames(merged))
	
	if(scde){
		info("Found zscore columns: ")
		print(colnames(merged)[zscore.cols])
		zscore.cols = grep("Corrected.Z.score", colnames(merged))
		info("Calculating bounds on Z-score (95-CI")
		Z.CI.95 = apply(merged[, zscore.cols], 1, bayesCI)
		merged$Z.Lower95 = Z.CI.95[3,]
		merged$Z.Upper95 = Z.CI.95[1,]
		merged$Z.Mean = Z.CI.95[2, ]
		merged$Z.Stouffer = apply(merged[, zscore.cols], 1, unweighted.z)
	}
	

	# info("Found FC columns: ")
	# print(colnames(merged)[fc.cols])

	# info("Found p-val columns:")
	# print(colnames(merged)[pval.cols])

	# info("Found Q-val columns:")
	# print(colnames(merged)[Qval.cols])
	
	merged[is.na(merged)] <- 0
	
	# Calculate metrics
	# ----------------------

	#merged$F1.Like = f1(merged[ c("FC.StN", "Z.StN")])

	# p.CI.95 = apply(merged[, pval.cols], 1, bayesCI)
	# merged$p.Lower95 = p.CI.95[3,]
	# merged$p.Upper95 = p.CI.95[1,]
	merged$p.Max = apply(merged[, pval.cols], 1, max, na.rm=T)
	merged$Q.Max = apply(merged[, Qval.cols], 1, max, na.rm=T)
	merged$p.Fisher = apply(merged[, Qval.cols], 1, combine.p.fisher)
	merged$Q.Fisher = p.adjust(merged$p.Fisher)

	if(do.fast){
		info("Calculating FC confidence intervals (CI-95)")
		FC.CI.95 = apply(merged[, fc.cols], 1, CI)
	}else{
		info("Calculating FC confidence intervals (bayes bootstrap)")
		FC.CI.95 = apply(merged[, fc.cols], 1, bayesCI)
	}
	
	merged$FC.Lower95 = FC.CI.95[3, ]
	merged$FC.Upper95 = FC.CI.95[1, ]
	merged$FC.Mean = FC.CI.95[2, ]
	merged$FC.Min = apply(merged[, fc.cols], 1, min, na.rm=T)
	if(scde){
		merged$FC.Z.Mean = rowMeans(merged[,c("Z.Stouffer", "FC.Mean")])
	}
	
	
	# rank markers
	# ----------------------------------
	if(scde){
		info("SCDE is ON. Ranking markers by Z.Lower95 ")
		merged = merged[order(merged$Z.Lower95 , decreasing=T), ]
	}else{
		info("SCDE is OFF. Ranking markers by FC.Min ")
		fcmin.rank = order(merged$FC.Min , decreasing=T)
		#print(length(fcmin.rank))
		#print(nrow(merged))
		merged = merged[fcmin.rank, ]
	}
	info("Reordering columns")
	#info("Adding adjusted expression column")
	#merged = merged[order(merged$Z, decreasing=T), ]
	#merged$Adjusted.Expression = merged$Mean.Log2.Ct / merged$Fraction.Expressing
	merged = merged[ , order(colnames(merged))]
	col_idx <- grep("GENE_SYMBOL", colnames(merged))
	others <- grep("GENE_SYMBOL", colnames(merged), invert=T)
	merged = merged[, c(col_idx, others)]
	info(sprintf("Returning %s ranked genes", nrow(merged)))
	return (merged)
}



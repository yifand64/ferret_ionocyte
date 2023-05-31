
# #test
# # use the results of pairwise SCDE calculations
# # to rank genes, finding the most significantly
# # up and down regulated genes for each group
# markers.rank.debug <- function(compare.obj)
# {
# 	pairs = compare.obj$pairs
# 	diffs = compare.obj$diffs
# 	groups = unique(c(pairs))
# 	n_groups = length(groups)
# 	#info(sprintf("Ranking most specific markers for %s ", groups))
# 	# print(groups)
# 	# print(pairs)
# 	for(i in 1:n_groups)
# 	{
# 		# pair = unlist(pairs[,i])
# 		group = groups[i]
# 		cat("\n\n")
# 		info(sprintf("Ranking markers for %s [%s of %s]", group, i, n_groups))
# 		cat("--------------------------------------------------------------\n")
# 		pairs.up = as.character(unlist(pairs[1, ]))
# 		pairs.down = as.character(unlist(pairs[2, ]))
# 		# find the list of comparisons done for this group, keeping track
# 		# of the direction of fold change
# 		info(sprintf("Looking for comparisons containing %s", group))

# 		comparisons.up = which(pairs.up == group) #grep(sprintf("^%s$", as.character(group)), pairs.up) # use the anchors ^$ for an exact grep match
# 		comparisons.down = which(pairs.down == group) #grep(sprintf("^%s$", as.character(group)), pairs.down)

# 		print("Comparisons UP:")
# 		print(comparisons.up)
# 		print("Comparisons DOWN:")
# 		print(comparisons.down)
# 		comparisons = unlist(c(comparisons.up, -comparisons.down))

# 		# merge the pairwise comparisons to calculate average ranks 
# 		merged = NULL

# 		info(sprintf("Merging %s pairwise comparisons", length(comparisons)))
# 		print(comparisons)

# 		start.dir = getwd()
# 		dir.create(group)

# 		info(sprintf("Moving to %s directory", group))
# 		setwd(group)
		
# 		# #-----
# 		# # TODO
# 		# #	this whole section unnecessarily distinguishes between up and 
# 		# #	downregulated genes. should simply merge the all_ genes and
# 		# #	then use the sign of the fold change to select top UP and DOWN
# 		# #	don't need @sig or @sig.down at all

# 		for(i in comparisons)
# 		{	
# 			compare.fail = FALSE
# 			if(i < 0)
# 			{
# 				if(!is.null(diffs[[i*-1]]))
# 				{
# 		 			d.all =  diffs[[i*-1]]@all
# 		 			directional.cols = which(colnames(d.all)!= "p" & colnames(d.all) != "GENE_SYMBOL")
# 		 			d.all[, directional.cols] = d.all[, directional.cols] * -1
# 		 			# rename the potentially ambiguously named 'p' column
# 		 			pval.col = which(colnames(d.all)=="p")
# 		 			colnames(d.all)[pval.col] = "p.adj"
# 		 			compare.name = diffs[[i*-1]]@name

# 		 			# reorder the fold change label
# 		 			compare.name = paste(strsplit(compare.name, "_vs_")[[1]][2], 
# 		 				strsplit(compare.name, "_vs_")[[1]][1], sep="_vs_")
# 		 			#compare.name = sprintf("%s_flipped", compare.name) #this is only for debug
# 		 			info(sprintf("For %s", compare.name))
# 		 			warn("Flipping the fold change")
# 		 			d.all <- d.all[order(d.all$Lower.bound, decreasing=T), ]	 
# 				}else
# 				{
# 					compare.fail = T
# 				}	
# 			}else
# 			{
# 				if(is.null(diffs[[i]]))
# 				{
# 					compare.fail = T
# 				}else{
# 					d.all = diffs[[i]]@all
# 					# rename the potentially ambiguously named 'p' column
# 		 			pval.col = which(colnames(d.all)=="p")
# 		 			colnames(d.all)[pval.col] = "p.adj"
# 					compare.name = diffs[[i]]@name
# 					d.all <- d.all[order(d.all$Lower.bound, decreasing=T), ]	 
# 					info(sprintf("For %s", compare.name))
# 				}
				
# 			}
# 			if(!compare.fail)
# 			{
# 				# write.output(d.all, d.up, compare.name, de.sig.down=d.down)
# 				write.table(d.all, file=paste("all_", compare.name, ".txt", sep=""), sep="\t", row.names=F, quote=F)
# 				# keep only relevant columns and add a rank column
# 				d.all = d.all[ , c("GENE_SYMBOL", "Log2.Fold.Change.MLE", "p.adj")]
# 				d.all[, paste("Rank", compare.name, sep="_")] = 1:nrow(d.all)	

# 				if(is.null(merged))
# 				{
# 					merged = d.all
# 					first.compare.name = compare.name 
# 				}else
# 				{
# 					info(sprintf("Merging %s comparisons", compare.name))
# 					merged = merge(merged, d.all, by="GENE_SYMBOL", all=T, suffixes=c("",paste("_",compare.name, sep="")))
# 				}
# 			}else{
# 				warn("Comparison failed. Skipping")
# 			}
# 		}

		
# 		info("Ranking genes based on merged comparisons")
# 		rank.cols = grep("Rank", colnames(merged))
# 		fc.cols = grep("Fold.Change", colnames(merged))
# 		pval.cols = grep("p.adj", colnames(merged))
# 		# info("Found pval columns: ")
# 		# print(colnames(merged)[pval.cols])

# 		# info("Found rank columns: ")
# 		# print(colnames(merged)[rank.cols])

# 		#merged$Average.Rank = rowMeans(merged[, rank.cols])
# 		info("Calculating average and median ranks")
# 		merged$Average.Rank = rowMeans(merged[, rank.cols], na.rm=T)
# 		merged$Median.Rank = apply(merged[, rank.cols], 1, median, na.rm=T)

# 		# replace NAs by Infs so we can sort the genes that appear 
# 		# in only one pairwise comparisons
# 		merged[is.na(merged)] = Inf
# 		merged$Min.Rank = apply(merged[, rank.cols],1, min)
		
# 		merged$max_p = apply(merged[, pval.cols], 1, max, na.rm=T)
# 		merged$min_fc = apply(merged[, fc.cols], 1, min, na.rm=T)
		
# 		# the best markers have high min_fc and low max_p, so score them
# 		# by their distance (in the max_p/min_fc plane) from the point (0, max(min_fc))
# 		maxfc=max(merged$min_fc)
# 		merged$Score = 10 - sqrt((maxfc*merged$max_p)^2 + (merged$min_fc-maxfc)^2)
# 		merged$Score[merged$max_p==1] <- 0
# 		merged$Score[merged$Score < 0] <- 0
# 		merged$Score[merged$min_fc < 0] <- 0

# 		# rank by score
# 		merged = merged[order(merged$Score, decreasing=T),]


# 		# visualise the best markers
# 		g=ggplot(merged, aes(x=max_p, y=min_fc, label=GENE_SYMBOL, colour=Score)) + 
# 			geom_point() + theme_bw() + 
# 			geom_text(data=subset(merged, max_p < 1 & min_fc > 0), vjust=2, cex=1.5, colour="black") + 
# 			scale_colour_gradientn(colours=colorRampPalette(brewer.pal(9,"OrRd"))(50)) + 
# 			ylab("Minimum Estimated Log2 FC") + 
# 			xlab("Maximum p value")
# 		ggsave(g, filename="markers.pdf", width=7, height=6)		
# 		# replace the Infs by NAs again
# 		is.na(merged) <- do.call(cbind, lapply(merged, is.infinite))
# 		rank.cols = c("GENE_SYMBOL", grep("Rank", colnames(merged), value=T))

# 		# above we set any gene that has a maximum pvalue = 1 (non-significant),
# 		# to have a zero score. This is going to be most genes, so we can't just
# 		# rank by the score. But the score does rank the *best* markers.
# 		best.markers = which(merged$Score > 0)
# 		unranked = which(merged$Score == 0)

# 		upreg.against.all = length(best.markers)
# 		info(sprintf("%s genes UP regulated in all pairwise comparisons.", upreg.against.all))

# 		# the first comparison wasn't merged, so it doesn't have a suffix. add one
# 		unlabeled = c("Log2.Fold.Change.MLE", "p.adj")
# 		unlabeled.indices = which(colnames(merged) == "Log2.Fold.Change.MLE" | colnames(merged) == "p.adj")
# 		colnames(merged)[unlabeled.indices] = paste(unlabeled, first.compare.name, sep="_")

# 		markers.table = paste("markers_", group, ".txt", sep="")
# 		markers.table.best = paste("markers_", group, "_BEST.txt", sep="")

# 		info(sprintf("Writing BEST marker genes to %s", markers.table.best))
# 		info(sprintf("Writing all ranked marker genes to %s", markers.table))

# 		merged = merged[ , order(colnames(merged))]
# 		col_idx <- grep("GENE_SYMBOL", colnames(merged))
# 		others <- grep("GENE_SYMBOL", colnames(merged), invert=T)
# 		merged = merged[, c(col_idx, others)]

# 		merged[unranked, ] = merged[unranked[order(merged[unranked, "Average.Rank"])], ]
# 		write.table(merged, file=markers.table, sep="\t", quote=F, row.names=F)
# 		write.table(merged[best.markers, ], file=markers.table.best, sep="\t", quote=F, row.names=F)

# 		cat("\n\n")
# 		setwd(start.dir)
# 	}
# 	info(sprintf("Completed ranking markers for %s groups", n_groups))
# }



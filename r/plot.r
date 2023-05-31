

library(shiny)
library(ggvis)
library(ggplot2)
library(reshape2)
library(cowplot)
library(dplyr)


# input should be a categorical vector (factor or character). See 'waffle chart' here:
#https://r-statistics.co/Top50-Ggplot2-Visualizations-MasterList-R-Code.html
## NOTE: nrows is fixed at 10. i haven't worked out why this breaks for other values of nrows
waffle_plot = function(x, subtitle="Frequency of category levels", caption="", colors=NULL)
{
	nrows=10
	## Prep data (nothing to change here)
	df <- expand.grid(y = 1:nrows, x = 1:nrows)
	y = table(x) * ((nrows*nrows)/(length(x)))
	categ_table <- floor(y)
	if(is.null(colors)){colors=default.cols(length(categ_table))}
	# simple algorithm to make sure rounded values add to 100. 
	# Adapted from https://revs.runtime-revolution.com/getting-100-with-rounded-percentages-273ffa70252b
	diff = 100 - sum(categ_table)
	decimal_part = y - floor(y)
	for(i in 1:diff){
		categ_table[order(decimal_part, decreasing=T)[i]] = categ_table[order(decimal_part, decreasing=T)[i]] + 1
	}

	df$category <- factor(rep(names(categ_table), categ_table))  
	## Plot
	ggplot(df, aes(x = x, y = y, fill = category)) + 
	        geom_tile(color = "black", size = 0.5) +
	        scale_x_continuous(expand = c(0, 0)) +
	        scale_y_continuous(expand = c(0, 0), trans = 'reverse') +
	        scale_fill_manual("", values=colors) +
	        labs(title="Waffle Chart", subtitle=subtitle,
	             caption=caption) + 
	        theme(panel.border = element_rect(size = 2),
	              plot.title = element_text(size = rel(1.2)),
	              axis.text = element_blank(),
	              axis.title = element_blank(),
	              axis.ticks = element_blank(),
	              legend.title = element_blank(),
	              legend.position = "right")
}


### internal, called by FeatureP
single_plot = function(df, verbose=F, pt.size=1.25, light.theme=T, max.value=0, density=F, pt.stroke=0.1, cols=NULL, max.value.rel=0.75, reduc.dims=c("tSNE_1", "tSNE_2")) {
	if(light.theme){
		plot_theme = theme_cowplot()
	}else{
		plot_theme = theme(panel.grid.major = element_line(colour = "grey10"), panel.grid.minor = element_line(colour = "grey11"), 
			panel.border= element_blank(),panel.background =  element_rect(fill = "grey12"))
	}
	g = as.character(unique(df$variable))
	formatting = theme(plot.title = element_text(size = 14, face = "bold.italic")) + plot_theme
	if(verbose){val.num = as.numeric(as.character(df$value)); info(sprintf("Plotting %s", g))
		}else{val.num = suppressWarnings(as.numeric(as.character(df$value)))}
	if(verbose){dtype = check.class(df$value);
		info(sprintf("Feature %s has datatype: %s", g, dtype))
		info(sprintf("NAs?=%s, All integers?=%s, All positive?=%s, Less than 100?=%s", any(is.na(val.num)), all(val.num %% 1 == 0), min(val.num) >= 0, max(val.num) < 100))}
	if(any(is.na(val.num)) | (all(val.num %% 1 == 0) & max(val.num) < 100)) ## if there's NAs, the value has strings. if they are all integers and in the 0-100 range, probably cluster labels
	{
		n = length(unique(df$value))
		if(verbose){info(sprintf("Feature %s is discrete [%s values]", g, n))}
		df$value <- as.character(df$value)
		if(is.null(cols)){cols=default.cols(n)}
		color_theme = scale_colour_manual("", values=cols, guide=guide_legend(override.aes = list(size=3))) 
	}else{
		if(light.theme){if(is.null(cols)){cols=weather_heat(25)}}else{if(is.null(cols)){cols=image_cols}}
		if(verbose){info(sprintf("Feature %s is continuous", g))}
		if(max.value.rel > 0){
			if(max.value > 0){stop("Set max.value or max.value.rel, not both!")}
			max.value = max.value.rel * max(val.num) #quantile(val.num, max.value.rel) 
		}
		if(max.value > 0){
			if(verbose){info(sprintf("Clipping values over %s", max.value))}
			val.num[val.num > max.value] <- max.value
		}
		df$value <- val.num
		color_theme = scale_color_gradientn("",colors = cols)
		if(all(val.num == 0))
		{
			info(sprintf("Feature %s is not detected in any cell!", g))
			color_theme = scale_fill_manual(values = c("red"), label = "zero")
		}
	}
	rv = ggplot(df, aes_string(x=reduc.dims[1], y=reduc.dims[2], color="value")) + geom_point(size=pt.size, stroke=pt.stroke) + 
		color_theme + ggtitle(g) + formatting  + xlab(reduc.dims[1]) + ylab(reduc.dims[2]) 
	if(density){rv = rv + stat_density2d()}
	return(rv)
}

### Convenience method to switch between UMAP and tSNE.
set_default_reduction = function(obj, reduc)
{
	obj@misc = reduc
	obj
}

### Seurat-2 compatible plotting function. replaces feature.plot.scale
FeatureP <- function(obj, features, cols=NULL, pt.size=0.5, verbose=F, light.theme=F, density=F, max.value=0, use.raw=F, pt.stroke=0.1, max.value.rel=0.75, ncol=NULL, reduction.use=NULL)
{	
   if(is.null(reduction.use)){reduction.use = names(obj@reductions)[length(obj@reductions)]}
   identifier = toupper(grep(reduction.use, names(obj@reductions), ignore.case=T, value=T))
   if(identifier == "PCA"){identifier = "PC"} 
   if(identifier == "TSNE"){identifier = "tSNE"} 
   reduc.dims=paste0(identifier, sep="_", 1:2)
   pList = list()
   # print(features)
   for(i in 1:length(features))
   {
	   		f = features[i]
	   		if(verbose){info(sprintf("Iteration %s. Feature: %s", i, f))}
	   		slot = "data"; if(use.raw){slot = "counts"}
	   		d = FetchData(obj,c(reduc.dims, f) , slot=slot)
	   		#print(head(d))
	   		pList[[i]] <- single_plot(melt(d,  id.vars=reduc.dims), max.value.rel=max.value.rel, reduc.dims=reduc.dims,
	   			max.value=max.value, verbose=verbose, light.theme=light.theme, pt.size=pt.size, pt.stroke=pt.stroke, cols=cols, density=density)
   }
   cowplot::plot_grid(plotlist = pList, ncol = ncol)
}


### watch out, can't do pt.stroke=0 with cairo_pdf
ViolinP <-function(obj, features, group.by="ident", cols=NULL, pt.size=0.5, verbose=F, significance=T, significance.compare.all=F, 
	significance.map_level=F, ylab="log2(TPM+1)", pt.stroke=0, significance.test=wilcox.test, significance.textsize=2.5, legend=F,
	pt.col="grey25", summary=T, sort=F, light.theme=T, x.axis.angle=45, show.median=F, nrow=NULL, groups.only=NULL)
{
   	if(legend){
   		formatting = theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
					  strip.text.x = element_text(size = 14, face = "bold.italic"), 
					  strip.background = element_blank()) 
   	}else{
   		formatting = theme(axis.text.x = element_text(angle = x.axis.angle, hjust = 1), 
					  strip.text.x = element_text(size = 14, face = "bold.italic"), 
					  strip.background = element_blank())
   	}
	if(significance){library(ggsignif)}
	if(pt.size>0){library(ggbeeswarm)}
	d = FetchData(obj, c(features, group.by))
	if(!is.null(groups.only)){d = d[d[, group.by] %in% groups.only,]}
	d = melt(d, id.vars=group.by)
	if(!is.null(groups.only)){
		info("Resetting levels"); 
		d[, group.by] = factor(d[, group.by], levels=groups.only)
	}else{
		if(is.numeric(d[, group.by])) {# preserve numeric character if its there:
					d[, group.by] = factor(as.character(d[, group.by]), levels=sort(unique(d[, group.by]))) ## necessary to ensure that grouping by numeric cluster ids doesn't give an error
		}else{
			orig_levels = levels(factor(d[, group.by]))
			d[, group.by] = factor(as.character(d[, group.by]), levels=orig_levels) ## necessary to ensure that grouping by numeric cluster ids doesn't give an error
		}
	}
	if(sort){
		y = aggregate(d$value, by=list(d[, group.by]), FUN=mean)
		d[, group.by] = factor(d[, group.by], levels=levels(d[, group.by])[y$Group.1[order(y$x, decreasing = T)]])
	}
	if(is.null(cols)){cols=weather_heat(length(unique(d[,group.by])))} #material.heat(length(unique(d[,group.by])))}
	if(legend){legend.use="legend"}else{legend.use=F}
	g = ggplot(d, aes_string("x"=group.by, "y"="value", "fill"=group.by)) + 
		geom_violin(scale="width",adjust=1,trim=TRUE, alpha=0.85, color=NA) + 
		scale_fill_manual(values=cols, guide=legend.use) + xlab("") + ylab("") + theme_cowplot()
	if(is.null(nrow)){g = g + facet_wrap(~variable, scales="free")}else{g = g + facet_wrap(~variable, scales="free", nrow=nrow)}
	if(pt.size>0){g = g + geom_quasirandom(varwidth = TRUE, size=pt.size, stroke=pt.stroke, colour=pt.col, show.legend=FALSE, method = "pseudorandom")}
	if(summary){g = g + stat_summary(color="black", fun.data="mean_cl_boot", fun.args=list(conf.int=.99), size=0.8, geom="errorbar", width=0.1) + 
						stat_summary(color="white", fun.data="mean_cl_boot", fun.args=list(conf.int=.99), size=0.5, geom="errorbar", width=0.1) + 
						stat_summary(color="black", geom="point", fun.y=mean, size=3, stroke=pt.stroke) + 
						stat_summary(color="white", geom="point", fun.y=mean, size=2.5, stroke=pt.stroke)}
	if(significance){
		group.levels = levels(d[,group.by])
		if(significance.compare.all){compar = combn(group.levels, 2)}else{
			#print(group.levels)
			compar = t(cbind(group.levels[-length(group.levels)], group.levels[-1])) # just consecutive pairs
		}
		compar = as.list(as.data.frame(compar, stringsAsFactors=F)) # convert to list
		#print(compar)
		g = g + geom_signif(comparisons = compar, map_signif_level = significance.map_level, test=significance.test, textsize=significance.textsize)
	}
	if(show.median){g = g + stat_summary(color="blue", geom="point", fun.y=median, shape=5, size=3)}
	g + formatting + ylab(ylab) #+ scale_y_continuous(expand = c(0,0))
}


interact_plot = function(obj, gene, cell, celltype.id="res.1_merged", condition.id="condition", condition.levels=NULL,
								condition.id.2=NULL, batch.id="batch", pt.size=0, pt.stroke=0, pdf=F)
{
    d = FetchData(obj, c(condition.id, batch.id, celltype.id, gene, condition.id.2))
    d$celltype = as.character(d[, celltype.id])
    d$celltype[d$celltype != cell] <- "All cells"
    if(is.null(condition.levels)){condition.levels=unique(d$condition)}
    d$condition = factor(d$condition, levels=condition.levels)
    colnames(d)[4] = "gene"
    poly = geom_violin(scale="width",adjust=1,trim=TRUE) #geom_boxplot(outlier.size = 0.1, width=0.2) 
    g = ggplot(d, aes(x=celltype, y=gene)) + poly + 
    	#geom_point(size=pt.size, stroke=pt.stroke) + 
    	geom_quasirandom(varwidth = TRUE, size=pt.size, stroke=pt.stroke, show.legend=FALSE, method = "pseudorandom", width=0.9) +
    	geom_signif(comparisons = list(c("All cells", cell))) + ylab(gene) + xlab("") + 
    	stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) + 
    	stat_summary(fun.y=mean, geom = "point", color="red")
    if(is.null(condition.id.2)){g = g + facet_wrap(as.formula(paste("~", condition.id)))}else{
    	g = g + facet_grid(as.formula(paste(condition.id, "~", condition.id.2)))
    }
    if(!pdf){return(g)}else{
    	cairo_pdf(sprintf("%s_%s.pdf", cell, gene)); print(g); dev.off()
    }
}

dot.plot <- function(data, groups, 
	genes=NULL, 
	groups.only=NULL,
	cols=material.heat(50), 
	border.col="black",
	mu.max=0,
	transpose=F,
	alpha.max=0,
	cex.all=16, 
	size.range=c(1,6),
	cluster=T,
	ret.data=F,
	alpha.per.group=T,
	dist.method="euclidean",
	x.axis.angle=90)
{

	print(table(groups))

	
	if(!is.null(genes)){
		info(sprintf("%s of %s genes found in data", length(intersect(genes, rownames(data))), length(genes)))
		genes = intersect(genes, rownames(data))
		data = data[genes,]
	}

	#print(length(groups))
	#print(dim(data))
	#print(rownames(data))
	info("Calculating mu")
	mean.non.zero = function(x){ if(sum(x>0) > 0) mean(x[x>0]) else 0}
	mu = group.means(data, groups, fn=mean.non.zero)
	#print(mu)
    if(sum(is.na(mu)) > 0)
    {
    	warn("Removing NAs from mu!")
    	mu = na.omit(mu)
    }
    mu[!is.finite(mu)] <- 0
    #print(mu)

    info("Calculating alpha")
    if(alpha.per.group){
		alpha = group.means(counts = data, groups=groups, fn = function(x){sum(x>1)/length(x)})
		if(sum(is.na(alpha)) > 0)
	    {
	    	warn("Removing NAs from alpha!")
	    	alpha = na.omit(alpha)
	    }
    }else{
    	### TODO: calculate alpha across the whole data, not in each group.	
    }
    mmu = melt(mu)
    malpha = melt(alpha)
    #print(head(mmu))
    #print(head(malpha))
    #print(malpha)
    merge.by = colnames(mmu)[1:2]
    if(!all(merge.by == colnames(malpha)[1:2])){stop("Something went wrong during the merge!")}

    df = merge(mmu, malpha, suffixes=c("_mu", "_alpha"), by=merge.by) 
    if(ret.data){return(df)}
    if(cluster){
    	ord = hclust( dist(df[,3:4], method = dist.method), method = "ward.D" )$order
		
	}else{
		ord = 1:nrow(df)
		df[,merge.by[1]] = factor(df[,merge.by[1]], levels=rev(genes))
	}
    df[,merge.by[1]] = factor(df[,merge.by[1]], levels=unique(df[,merge.by[1]])[ord])

    
    if(!is.null(groups.only)){

    	df = df[df[,merge.by[2]] %in% groups.only, ]
    	df[merge.by[2]] = factor(df[,merge.by[2]], levels=groups.only)

    }
    #print(df)
    if(mu.max > 0){
    	vals = as.numeric(df$value_mu)
    	vals[vals>mu.max] <- mu.max
    	df$value_mu = vals
    }
    if(alpha.max > 0){
    	vals = as.numeric(df$value_alpha)
    	vals[vals>alpha.max] <- alpha.max
    	df$value_alpha = vals
    }
    if(transpose){g = ggplot(df, aes_string(y=merge.by[2], x=merge.by[1]))}else{g = ggplot(df, aes_string(y=merge.by[1], x=merge.by[2]))}
	g + theme_bw() + 
				geom_tile(color = "black", fill="white") +
				geom_point(aes(fill = value_mu, size =value_alpha), shape=21, color=border.col) + 
				scale_fill_gradientn("Mean expression (mu)", colors=cols) + 
				theme(axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"),
						axis.text.x = element_text(angle = x.axis.angle, hjust = 1),
						axis.line.y = element_line(color="gray15", size = 0.75, lineend="round")) + 
				   
				theme(panel.grid.major = element_blank(),
			    	panel.grid.minor = element_blank(),
					panel.border= element_blank(),
			 		panel.background = element_blank(),
			    	strip.background = element_blank(), 
			    	# strip.text.x = element_text(size = cex.strip, face="bold"),
			    	# axis.title.x = element_text(vjust=-0.5),
			    	legend.key = element_blank(),
			    	text = element_text(size=cex.all, color="gray15")) +
				scale_size_continuous("Fraction of cells\nexpressing (alpha)", range=size.range, guide=guide_legend(override.aes = list(color="black"))) + xlab("") + ylab("") 
}

# ## Another possible palette for Atlas paper.
# celltype.cols.material = list("Goblet"=material.700[1], 
# 	"Paneth"=material.700[17], 
# 	"Stem"=material.700[3], 
# 	"Endocrine"=material.700[19],
# 	"Tuft"=material.700[5],
# 	"TA.Early"=material.700[15], 
# 	"TA.G1"=material.700[7],
# 	"TA.G2"=material.700[8],
# 	"TA"=material.700[11],
# 	"Enterocyte.Progenitor.Early"=material.700[9], 
# 	"Enterocyte.Progenitor.Late"= material.700[10], 
# 	"Enterocyte.Immature.Proximal"=material.700[13], 
# 	"Enterocyte.Immature.Distal"=material.700[12], 
# 	"Enterocyte.Mature.Proximal"=material.700[4], 
# 	"Enterocyte.Mature.Distal"=material.700[2], 
# 	"Enterocyte.Immature"=material.700[6],
# 	"Enterocyte.Mature"=material.700[16],
# 	"M-Like Tuft"=material.700[18],
# 	"Microfold"=material.700[18]
# 	)


### generate rounded breaks for plotting axes
round.breaks <- function(x, nearest=100)
{
	unlist(lapply(x, function(y){(y %/% nearest) * nearest}))
}


get.nice.breaks <- function(vals=NULL, min.value=NULL, max.value=NULL, breaks.every=100)
{

	if(!is.null(vals)){
		min.value = min(vals)
		max.value = max(vals)
	}else{
		if(is.null(min.value) | is.null(max.value))
		{
			error("Provide a vector of values or a min.value and max.value")
			return(FALSE)
		}
	}
	print(breaks.every)
	#info(sprintf("Getting breaks every %s", breaks.every))
	m = min.value %/% breaks.every
	n = 1 + (max.value %/% breaks.every)
	m:n * breaks.every
	
}


## Match cell-type names to colours to standardise plots.
## Assumes the input is a factor.
get.celltype.cols <- function(groups)
{
	if(!is.factor(groups))
	{
		warn("Coercing to a factor!")
		groups = factor(groups, levels=unique(groups))
	}
	known = names(celltype.cols)
	celltypes = levels(groups)
	cols = unlist(unique(celltype.cols))

	if(!any(celltypes %in% known)){
		error("Could not match any levels in the given groups factor to known celltypes!")
		error("Given: ")
		error(paste(celltypes, collapse=", "))
		error("Known:")
		error(paste(known, collapse=", "))
		return (NULL)
	}
	rval = vector()
	i = 1
	for (ct in celltypes)
	{
		if(ct %in% known){rval[i]=celltype.cols[[ct]]}else{
			rval[i] = NA
		}
		i = i + 1
	}
	no.match = celltypes[is.na(rval)]
	if(length(no.match) > 0){
		warn(sprintf("Could not find a known cell type match for celltypes %s", paste(no.match, collapse=", ")))
	}
	
	already.used = rval[!is.na(rval)]
	warn(sprintf("Already used these colours %s", paste(already.used, collapse=", ")))

	available =  cols[!cols %in% already.used]
	warn(sprintf("Still have these %s", paste(available, collapse=", ")))

	no.match.ids = which(is.na(rval))
	for(i in 1:length(no.match))
	{
		j = no.match.ids[i]
		rval[j] = available[i]
	}
	return(rval)
}




info("Loading plotting functions")
set1.cols = function(n){colorRampPalette(brewer.pal(n, "Set1"))(n)}


violin <- function(seurat.obj, features.plot, 
	xlab="",
	ylab=NULL,
	group.by=NULL,
	group.by.second=NULL, # use this to facet the plot. only one feature at at time though.
	point.size=0, # not actually 0, v small
	scale.free.facet=F,
	groups.only=NULL,
	groups.only.second=NULL,
	groups.rename=NULL,
	color.plot = T,
	use.cols=NULL,
	use.celltype.cols=F, # standardise colour allocation of known celltypes e.g. Goblet, Paneth
	background="white",
	grid.col="white",
	legend=F,
	legend.title="Group",
	sort=F,
	sort.fn=mean,
	pdf=NULL,
	axis.labels.x=T,
	axis.labels.y=T,
	cex.strip=13,
	cex.all=12,
	width=8,
	height=7,
	title=NULL,
	show.points.col="grey20",
	show.points.beeswarm.method="pseudorandom", #can also be 'tukeyDense' or 'tukey'
	show.points="density", #can also be none, or jitter
	ret.plot=T, #return the plot object
	y.breaks.every=NULL,
	nrow=NULL, ## if needed, specify the number of rows for multiple plots (e.g 2x2 grid -> nrow=2, 1x4 grid -> nrow=1)
	y.max=NULL,
	y.min=NULL,
	raw.data=F, ### if true, pull data from @raw.data
	violin.show.summary=T,
	violin.summary.fn=median,
	summary.bar.width=0.5,
	violin.alpha=1,
	boxplot.minimal=T,
	boxplot.jitter.width=2.5,
	boxplot.width=0.9,
	boxplot.notch=F,
	boxplot=F,
	test=F, #if true, use geom_signif to run and display significance stars
	text.col="black")
{
	
	if(!	tolower(show.points) %in% c("density", "none", "jitter")){stop("show.points must be one of 'density', 'jitter', or 'none'")}else{
		show.points=tolower(show.points)}

	features.plot = check.data(seurat.obj, features.plot)

	if(is.null(group.by))
	{
		groups = unlist(seurat.obj@ident)
	}else{
		groups = unlist(fetch.data(seurat.obj, group.by))
	}
	
	if(!is.null(group.by.second))
	{
		data = fetch.data.updated(seurat.obj, c(group.by.second,features.plot), use.raw=raw.data)
	}else{
		info(sprintf("FETCHING. raw.data=%s", raw.data))
		data = fetch.data.updated(seurat.obj, features.plot, use.raw=raw.data)
	}
	data$group = factor(groups)
	
	#print(head(data))
	if(!is.null(groups.only))
	{
		data = data[data$group %in% groups.only,]
		
		if(!is.null(groups.rename))
		{
			data$group = mapvalues(data$group, groups.only, groups.rename)
			groups.only = groups.rename
		}
		#print(unique(dm$group))
		#print(groups.only)
		data$group = factor(data$group, levels=groups.only)
		# print(levels(dm$group))
		# levels(dm$group) = 
		groups = data$group
	}

	if(!is.null(groups.only.second))
	{	
		keep = which(unlist(data[, group.by.second]) %in% groups.only.second)
		data = data[keep,]
	}

	## If no colours are provided
	n.groups = length(unique(data$group))
	if(is.null(use.cols))
	{
		use.cols = default.cols(n.groups)
	}


	## If need to order by value
	if(sort)
	{
		f = features.plot[1]
		info(sprintf("Sorting by %s of %s", as.character(substitute(sort.fn)), f))
		x = unlist(data[, f])
		y = aggregate(x~groups, FUN=sort.fn)
		colnames(y) = c("Group", "Mean")
		
		y = y[order(y$Mean, decreasing=T),]
		use.cols = use.cols[order(y$Mean, decreasing=T)]

		new.levels = unname(unlist(y$Group))
		print("New levels:")
		print(new.levels)
		data$group = factor(data$group, levels=new.levels)
	}

	## If need to match colours with other plots (overrides specified colours)
	if(use.celltype.cols)
	{
		if(!is.null(use.cols))
		{
			warn("Overriding provided colours and matching with known cell types")
		}
		info("Using standardised cell-type colours")
		use.cols = get.celltype.cols(data$group)
	}

	if(!is.null(group.by.second))
	{
		dm = melt(data, id.vars=c("group", group.by.second))
	}else{
		dm = melt(data, id.vars="group")
	}

	if(!is.null(title))
	{
		levels(dm$variable) = title
	}

	## Clean up variable names:
	# dm$variable = gsub("_", " ", dm$variable, fixed=T)
	# dm$variable = gsub(".", " ", dm$variable, fixed=T)
	# features.plot = gsub("_", " ", features.plot, fixed=T)
	# features.plot = gsub(".", " ", features.plot, fixed=T)

	## Keep features.plot in original order
	dm$variable = factor(dm$variable, levels=features.plot)
	# choose tick frequency
	# y.range =  max(dm$value) - min(dm$value)
	# y.breaks = seq(floor(min(dm$value)), floor(max(dm$value)), by=y.range/15)
	# if(max(dm$value)>20) y.breaks = round.breaks(y.breaks, round.y.breaks.to.nearest)
	if(is.null(y.breaks.every)){
		if(max(dm$value)>20){ ### this is unlikely to be an expression level as nothing has log2 > 20 (a million counts/tpm)
			y.breaks.every=100  
		}else{
			y.breaks.every=1
		}
	}

	if(max(dm$value)>20){### this is unlikely to be an expression level as nothing has log2 > 20 (a million counts/tpm)
		if(is.null(ylab))
		{
			ylab = features.plot 
		}
	}else{
		if(is.null(ylab))
		{
			if(raw.data){
				ylab = paste(features.plot, "Raw expression Level (Log2 UMI count)")
			}else{
				ylab = paste(features.plot, "Expression Level (Log2 TPM)")
			}
		}
	}
	

	y.breaks = get.nice.breaks(vals=dm$value, breaks.every=y.breaks.every)

	if(color.plot)
	{
		g = ggplot(dm, aes(y=value, x=group, fill=group, colour=group))
		median.bar.col = "black"
	}else{
		g = ggplot(dm, aes(y=value, x=group), color="black")
		median.bar.col = "black"
	}

	## Start building the plot
	#library(cowplot)
	g = g + ylab(ylab) + xlab(xlab) + theme_cowplot() + 
		theme(#axis.line.x = element_line(color="black", size = 0.75, lineend="round"),
				#axis.line.y = element_line(color="black", size = 0.75, lineend="round"),
				#panel.grid.major = element_line(colour = grid.col),
				#panel.grid.minor = element_line(colour = grid.col),
				#panel.border= element_blank(),
				#panel.background =  element_rect(fill = background),
				strip.background = element_blank(), 
				legend.key = element_blank(),
				strip.text.x = element_text(size = cex.strip, face="bold", colour=text.col),
				axis.title.x = element_text(vjust=-0.5, color=text.col),
				axis.text.y = element_text(color=text.col),
				axis.text.x = element_text(angle = 45, hjust = 1, color=text.col),
				text = element_text(size=cex.all, color=text.col)) 

	if(test){ # use ggsignif to run some wilcoxon rank-sum tests and show significance stars
		library(ggsignif)
		x = unique(dm$group)
		y = y = combn(x, 2)
		g = g + geom_signif(color="black", 
							test=wilcox.test, 
							comparisons = mapply(c, y[1,], y[2,], SIMPLIFY = FALSE), 
							map_signif_level=TRUE)
	}

	if(length(features.plot) == 1){ # cant manually scale the breaks for a facetted plot
		g = g + scale_y_continuous(breaks = y.breaks) #tick.freq.y
	}
		
	if(length(n.groups) > length(use.cols)){
		cols = colorRampPalette(use.cols)(n.groups) # interpolate if we dont have enough colours
	}else{
		cols = use.cols
	}
	

	if(legend)
	{
		info(sprintf("Adding %s legend", legend.title))
		g = g + scale_fill_manual(legend.title, values = cols, guide = "legend") + 
			scale_colour_manual(legend.title, values = cols, guide="legend") #+
			#guides(colour=guide_legend(title=legend.title)) + theme(legend.key = element_blank())
	}else{
		g = g + scale_fill_manual(legend.title, values = cols, guide=FALSE) +
			scale_colour_manual(legend.title, values = cols, guide=FALSE) 
	}


	if(!axis.labels.x){
		g = g + theme(axis.ticks.x=element_blank(), 
  						axis.text.x=element_blank())
	}
	if(!axis.labels.y){
		g = g + theme(axis.ticks.y=element_blank(), 
  						axis.text.y=element_blank())
	}


	if(boxplot)
	{
		if(boxplot.minimal)
		{
			g = g + geom_boxplot(outlier.colour = NULL, outlier.shape = NA, width=boxplot.width, notch=boxplot.notch, aes(colour=group)) +
					stat_summary(geom = "crossbar", show.legend=FALSE, width=boxplot.width*summary.bar.width, fatten=1.5, color=median.bar.col, 
						fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) })
		}else{
			g = g + geom_boxplot(width=boxplot.width)
		}
		
	}else
	{
		if(violin.show.summary)
		{
			g = g + geom_violin(alpha=violin.alpha, scale="width",adjust=1,trim=TRUE)
		}else{
			g = g + geom_violin()
		}
	}

	if(show.points=="jitter")
	{
		#g = g + geom_point(position=position_jitterdodge(dodge.width=0.6, jitter.width = boxplot.jitter.width), size=point.size) #alpha=0.25, 
		g = g + geom_jitter(height=0,size=point.size, colour=show.points.col, show.legend=FALSE)
	}else{
		if(show.points == "density"){
			library(ggbeeswarm)
			g = g + geom_quasirandom(varwidth = TRUE, size=point.size, colour=show.points.col, show.legend=FALSE, method = show.points.beeswarm.method)
		}else{
			if(show.points != "none"){stop(sprintf("Weird error! show.points is an impossible value [%s]", show.points))}
		}
	}

	if(violin.show.summary &! boxplot)
	{
		g = g + stat_summary(geom = "crossbar", show.legend=FALSE, width=boxplot.width*summary.bar.width, fatten=2.5, color=median.bar.col, 
						fun.data = function(x){ return(c(y=violin.summary.fn(x), ymin=violin.summary.fn(x), ymax=violin.summary.fn(x))) })
	}

	if(!is.null(group.by.second))
	{
		if(length(features.plot) > 1)
		{
			error("Cannot plot multiple features with two group.by variables. Remove 'group.by.second' argument or reduce 'features.plot' to 1.")
		}else{
			if(!scale.free.facet)
			{
				facet.scales = "fixed"
			}else{
				facet.scales = "free"
			}
			info(sprintf("Set facet.scale=%s", facet.scales))
			g = g + facet_wrap(as.formula(paste("~", group.by.second)), scales=facet.scales, nrow=nrow)
		}
		
	}else{
		if(length(features.plot) > 1){
			if(!scale.free.facet)
			{
				facet.scales = "fixed"
			}else{
				facet.scales = "free"
			}
			info(sprintf("Set facet.scale=%s", facet.scales))
			g = g + facet_wrap(~variable, scales=facet.scales, nrow=nrow)
		}
		
	}


	

	if(!is.null(title))
	{
		g = g + ggtitle(title)
	}

	if(!is.null(y.max) | !is.null(y.min))
	{
		if(is.null(y.max)){y.max=max(dm$value)}
		if(is.null(y.min)){y.min=min(dm$value)}
		g= g + coord_cartesian(ylim=c(y.min, y.max))
	}

	if(!is.null(pdf))
	{
		cairo_pdf(filename=pdf, width=width, height=height); print(g); dev.off()
	}else{
		if(ret.plot){
			return(g)
		}else{
			print(g)
		}
	}
	

	# if(!is.null(groups.only))
	# {
	# 	keep.cells = seurat.obj@cell.names[which(groups %in% groups.only)]
	# 	#print(head(keep.cells))
	# 	seurat.obj = subsetData(seurat.obj, cells.use=keep.cells)
	# 	#print(seurat.obj)
	# }
	
	# groups = unlist(fetch.data(seurat.obj, group.by))
	# n.groups = length(unique(groups))
	# x = vlnPlot(seurat.obj, features.plot, 
	# 						group.by = group.by,
	# 						size.x.use=8, 
	# 						size.y.use=8, 
	# 						size.title.use=10, 
	# 						cols.use=
	# 						do.ret = T)[[1]] + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
	
	# print(x)
	
}
check.class <- function(x)
{
	#info("Checking variable type")
	x <- as.character(unlist(x))
	if(any(is.na(as.numeric(x)))){
		rval="non-numeric"
	}else{
		if(sum(as.numeric(x) %% 1) == 0){
			rval="integer"
		}else{
			rval="double"
		}
	}
	#info(sprintf("Found %s", rval))
	rval
}


###unfinished! sep 14, '16
violin.double <- function(s.obj, gene, compare, celltype.id, max.dens=2, groups.only=NULL, group.by="Condition",width=9, height=7, pdf=NULL)
{
	#make a dataframe with three columns: 
	# y - value
	# x - celltype
	# m - condition

	my_data = fetch.data(s.obj, c(gene, celltype.id, group.by))
	
	if(!is.null(groups.only))
	{
		my_data = my_data[grep(paste(groups.only, collapse="|"), my_data[, celltype.id]),]
	}
	

	labs = unique(my_data[, celltype.id])
	n.cells = length(labs)
	colnames(my_data) = c("y", "x", "m")

	#library(dplyr)
	pdat <- my_data %>%
	  group_by(x, m) %>%
	  do(data.frame(loc = density(.$y)$x,
	                dens = density(.$y)$y))

	bks = 0:(n.cells-1)
	info("Breaks are: ")
	print(bks)
  	#Flip and offset densities for the groups
	pdat$dens <- ifelse(pdat$m == compare[1], pdat$dens * -1, pdat$dens)
	
	for(i in 1:n.cells){
		celltype = labs[i]
		info(sprintf("Adding %s to %s", (i-1), celltype))
		pdat$dens <- ifelse(pdat$x == celltype, pdat$dens + (i-1), pdat$dens)
	}

	# print(range(pdat$dens))
	# print(pdat$dens [abs(pdat$dens)>max.dens])

	# # remove singular density estimates:
	#pdat$dens [abs(pdat$dens)>max.dens] <- 0

	# print(pdat)

	
	#print(labs)
	info("plotting")
	g=ggplot(pdat, aes(dens, loc, fill = m, group = interaction(m, x))) + 
	  geom_polygon(colour="grey20", size=1) +
	  scale_x_continuous(breaks = bks, labels = labs) +
	  scale_fill_manual("", values=(c("cornflowerblue", "lightcoral"))) + 
	  ylab('Log2 expression level') + xlab("") + 
	  theme_bw() + ggtitle(gene) + 
	  theme(axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"),
	  axis.line.y = element_line(color="gray15", size = 0.75, lineend="round")) + 
		theme(
		    panel.border= element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
		    strip.background = element_blank(), 
		    plot.margin = unit(c(1.2,1.2,1.2,1.2), "cm"),
		    text = element_text(size=24, colour="gray22"))
	if(is.null(pdf))
	{
		print(g)
	}else{
		ggsave(g, filename=pdf, width=width,height=height)
	}
}

scatter <- function(seurat.obj, x.var, y.var, 
	xlab=NULL,
	ylab=NULL,
	zlab=NULL,
	raw.data=F,
	black.rim.pts=F,
	background=NULL,
	scale.alpha=F,
	highlight=NULL,
	highlight.cols=NULL, 
	split.by=NULL,
	groups.only = NULL,
	show.marginals=F,
	show.legend=T,
	show.cor=F,
	label=NULL,
	# NEED to consolidate these variables
	# -------------------------------------
	cluster.id = "DBclust.ident",
	colour.by=NULL, # can be a vector, will cause the plot to be facetted.
	colour.by.group=NULL,
	# -------------------------------------
	disp.max=0,
	disp.min=0,
	cex.strip=20,
	cex.all=20,
	grid.col=NULL,
	size.by=NULL,
	shape.by=NULL,
	legend.title.shape=NULL,
	legend.title=NULL,
	legend.point.size=4,
	legend.cex=10,
	legend.vjust=0,
	min.zero=F, # force the colour scale to start at zero, rather than the minumum value
	night=F,
	cols=NULL,
	use.celltype.cols=F, # standardise colour allocation of known celltypes e.g. Goblet, Paneth
	point.size=3,
	width=9, 
	height=7,
	ellipse=F,
	title="",
	hull=F,
	ret.plot=T, #return the plot object
	ret.data=F, #return the mappings for colour, size, shape etc used to generate the plot
	density=F, 
	density.cols=rev(brewer.pal(11,"Spectral")),
	density.col.pts=F,
	contour=T,
	interactive=F,
	pdf=NULL,
	fit.method="lm",
	add.curve=NULL, #plot (probably principal) curve. provide the output of principal.curve (its called fit1 in the example)
	fit.line=F)
{
	if (length(grep("PC", x.var))>0) {
		x.var = gsub("_", "", x.var)
		x.var = gsub(" ", "", x.var)
		pc = as.numeric(strsplit(x.var, "PC")[[1]][2])
		info(sprintf("Fetching PC %s", pc))
      	x = unlist(seurat.obj@pca.rot[,pc])
    }else{
    	x = fetch.data.updated(seurat.obj, x.var, use.raw=raw.data)
    }
	
	if (length(grep("PC", y.var))>0) {
		y.var = gsub("_", "", y.var)
		y.var = gsub(" ", "", y.var)
		pc = as.numeric(strsplit(y.var, "PC")[[1]][2])
		info(sprintf("Fetching PC %s", pc))
      	y=unlist(seurat.obj@pca.rot[,pc])
    }else{
    	y = fetch.data.updated(seurat.obj, y.var, use.raw=raw.data)
    }

    if(is.null(colour.by.group))
    {
    	colour.by.group = colour.by
    }
	df = data.frame(x,y)
	colnames(df) = c(x.var, y.var)
	if(is.null(size.by))
	{
		size.by= "point.size"
		#info(sprintf("Point size set to: %s", point.size))
		df[ , size.by] = point.size
	}else{
		info(sprintf("Sizing pts by %s ", size.by))
		df[ , size.by] = fetch.data(seurat.obj, size.by)
	}

	if(night)
	{
		if(is.null(background)){background="grey20"}
		if(is.null(grid.col)){grid.col="grey26"}
		if(is.null(cols)){cols=image_cols}
	}else
	{
		if(is.null(background)){background="white"}
		if(is.null(grid.col)){grid.col="white"}
	}

	pearson = cor.test(unlist(x), unlist(y), method="pearson")
	pearson.cor = round(pearson$estimate, 3)
	pearson.p = signif(pearson$p.value, 3)
	
	spearman = cor.test(unlist(x), unlist(y), method="spearman")
	spearman.cor = round(spearman$estimate, 3)
	spearman.p = signif(spearman$p.value,3)

	if(is.null(xlab))
	{
		xlab=x.var
	}
	if(is.null(ylab))
	{
		ylab=y.var
	}
	# cat(sprintf("Spearman = %f, Pearson = %f \n", spearman, pearson))
	correlation_label = paste("Pearson = ", pearson.cor," (p=", pearson.p,") \n", "Spearman = ", spearman.cor, " (p=", spearman.p, ")",sep="")
	
	if(!is.null(colour.by))
	{
		info(sprintf("Colouring by %s", colour.by))
		
		if(colour.by=="ident"){df[, colour.by]=seurat.obj@ident}else{
			if(length(colour.by)==1){
				df[, colour.by] = fetch.data(seurat.obj, colour.by)[,1]
			}else{
				df[, colour.by] = fetch.data(seurat.obj, colour.by)
			}
		}
		
		if(disp.max > 0)
		{
			above.max = which(df[, colour.by] > disp.max)
			df[above.max, colour.by] <- disp.max
		}

		if(disp.min < 0)
		{
			below.min = which(df[, colour.by] < disp.min)
			df[below.min, colour.by] <- disp.min
		}
		
		colnames(df) = c(x.var, y.var, size.by, colour.by)
		datatype.to.colour = check.class(df[, colour.by])
		numeric.data = datatype.to.colour == "double" 
		if(datatype.to.colour == "integer"){
			if(max(as.numeric(df[, colour.by])) > 50) #small integers are likely factors like cluster labels
			{
				info(sprintf("But, there are large values [max=%s]. Its probably numeric", max(df[, colour.by])))
				numeric.data = T # large integers are likely measurements like genes detected
			}else{
				info("Values are small. Assuming they reflect cluster labels.")
				df[,colour.by] = as.character(df[, colour.by])
				print(head(df))
				print(table(df[, colour.by]))
			}
		}

		info(sprintf("Numeric colourby data is %s", numeric.data))
	}else{
		numeric.data = F
	}

	if(!is.null(split.by))
	{
		df[ , split.by] = fetch.data(seurat.obj, split.by)
	}

	rownames(df) = seurat.obj@cell.names
	
	if(!is.null(label))
	{
		seurat.obj@data.info$Name = seurat.obj@cell.names
		df[, label] = fetch.data(seurat.obj, label)
	}else{
		
		df[, "Name"] = seurat.obj@cell.names
	}

	if(!is.null(shape.by))
	{
		#allowed.shapes = 15:18
		shape.ids = unlist(fetch.data(seurat.obj, shape.by))
		#n.shapes = length(unique(shape.ids))
		
		#shape.vector = mapvalues(shape.ids, unique(shape.ids), allowed.shapes[1:n.shapes])
		#print(table(shape.vector))
		df[,shape.by]=factor(shape.ids)

	}else{
		shape.by = "shape.by"
		df[,shape.by]=rep(1, nrow(df))
		shape.by = NULL
	}

	if(!is.null(groups.only))
	{
		#cluster.id = colour.by
		#print(cluster.id)
		keep.cells = which(unlist(fetch.data(seurat.obj, colour.by)) %in% groups.only)
		#print(keep.cells)
		df = df[keep.cells, ]
		info(sprintf("Only plotting %s the cells in clusters %s", length(keep.cells), paste(groups.only, collapse=" ,")))
		
		df[, colour.by] = factor(df[, colour.by], levels=groups.only)
		# if(colour.by %in% colnames(df))
		# {
		# 	print(table(df[, colour.by]))
		# }
	} 
	if(is.null(legend.title))
	{
		legend.title = colour.by
	}
	if(is.null(legend.title.shape))
	{
		legend.title.shape = shape.by
	}

	if(interactive)
	{
		if(is.null(colour.by)){colour.by="ident"}
		info("Calling gene plot")
		gene.plot(seurat.obj, x.var, y.var, id.use=colour.by, fit.line=fit.line)
	}else
	{

		# print(head(df))
		# print(sapply(df, typeof))
		# print(sapply(df, is.factor))

		if(!is.null(highlight))
		{
			if(!is.null(highlight.cols))
			{
				if(length(highlight.cols) != length(highlight))
				{
					error("Highlight cols must be the same length as the number of groups to highlight")
					return (FALSE)
				}
			}else
			{
				ncols = length(highlight)
				highlight.cols= default.cols(ncols)
			}

			info(sprintf("Highlighting %s", highlight))
			if(colour.by=="ident"){groups=seurat.obj@ident}else{
				groups = fetch.data(seurat.obj, colour.by)
			}
			
			not.in = which(!unlist(groups) %in% highlight)
			#groups = unlist(groups) %in% highlight
			
			# print(table(groups))
			#colour.by = cluster.id
			#df[, colour.by] = factor(mapvalues(groups, c(rep(T, length(highlight),F), c(highlight, "All cells")))
			groups.vec = unlist(groups)
			levels(groups.vec) = c(levels(groups.vec), "All Cells")
			groups.vec[not.in] = "All Cells"
			df[, colour.by] = factor(groups.vec)

			#print(table(df[, colour.by]))
			colnames(df) = c(x.var, y.var, size.by, colour.by)
		}

		# sort into alphabetical order so legend is ordered
		if(!numeric.data){
			df[, colour.by] = factor(df[, colour.by], levels=sort(unique(df[, colour.by])))
		}
		## If need to match colours with other plots (overrides specified colours)
		if(use.celltype.cols)
		{
			
			if(!is.null(cols))
			{
				warn("Overriding provided colours and matching with known cell types")
			}
			info("Using standardised cell-type colours")
			cols = get.celltype.cols(df[, colour.by])
		}

		# Clean up variable names (removing '.'s and '_'s)
		if(!numeric.data & !is.null(colour.by)){
			groups =  df[,colour.by]
			
			#print(table(groups))
			# print("Order")
			# print(unique(groups))
			# print(levels(groups))
			cleaned.groups = gsub(".", " ", groups, fixed=T)
			cleaned.groups = factor(cleaned.groups, levels=gsub(".", " ", levels(factor(groups)), fixed=T))
			
			df[,colour.by] = cleaned.groups
			# print(unique(cleaned.groups))
			#print(levels(cleaned.groups))
		}
		library(cowplot)
		p = ggplot(df, aes_string(x=x.var, y=y.var, label=label)) + 
				xlab(xlab) + ylab(ylab) + #+ theme_bw() 
				theme(axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"),
       				 axis.line.y = element_line(color="gray15", size = 0.75, lineend="round")) + 
				#  + 
				theme(panel.grid.major = element_line(colour = grid.col), 
			     	panel.grid.minor = element_line(colour = grid.col),
				 	panel.border= element_blank(),
			  		panel.background =  element_rect(fill = background),
			     	strip.background = element_blank(), 
			     	strip.text.x = element_text(size = cex.strip, face="bold"),
			     	legend.key = element_blank(),
			     	text = element_text(size=cex.all, colour="gray22")) 

		if(show.cor)
		{
			p = p + ggtitle(correlation_label)
		}else{
			p = p + ggtitle(title)
		}

		if(ellipse | hull)
		{
			# print(head(df))
			# print(table(df[, group.colour.col]))
			if(ellipse)
			{
				p = p + stat_ellipse(aes_string(colour=colour.by), size=1.5, alpha=0.5) 
			}else{
				#data=df[which(df$group.colour.col!=NA),], inherit.aes=F,
				#print(table(df[ ,group.colour.col]))
				p = p + stat_chull(aes_string(colour=colour.by), fill=NA, size=2) 
			}
			#p = p + scale_colour_manual("Cell Type", values=default.cols(length(unique(df[,group.colour.col]))))
		}

		if(!is.null(colour.by) | !(is.null(highlight)))
		{
				if((!is.null(highlight) | numeric.data) & scale.alpha) # if we need to vary the alpha,
				{			
					if(size.by=="point.size")
					{
						if(black.rim.pts)
						{
							p = p + geom_point(aes_string(fill = colour.by, alpha=colour.by, shape=21), size=point.size) #, 
						}else{
							p = p + geom_point(aes_string(colour = colour.by, alpha=colour.by, shape=shape.by), size=point.size) #, 
						}
						p = p + scale_size(guide='none')
					}else{
						if(black.rim.pts)
						{
							p = p + geom_point(aes_string(fill = colour.by, alpha=colour.by, shape=21, size=size.by)) 
						}else{
							p = p + geom_point(aes_string(colour = colour.by, alpha=colour.by, size=size.by, shape=shape.by))  
						}
						p = p  + scale_size_continuous(size.by)
					}
					
				}else{
					if(size.by=="point.size")
					{
						if(black.rim.pts)
						{
							p = p + geom_point(aes_string(fill = colour.by, shape=21), size=point.size)
						}else{
							p = p + geom_point(aes_string(colour = colour.by, shape=shape.by), size=point.size)
						}
						p = p + scale_size(guide='none')
					}else{
						if(black.rim.pts)
						{
							p = p + geom_point(aes_string(fill = colour.by, shape=21, size=size.by))
						}else{
							p = p + geom_point(aes_string(colour = colour.by, size=size.by, shape=shape.by))
						}
						p = p  + scale_size_continuous(size.by)
					}
				}
				if(black.rim.pts){p = p + scale_shape_identity()}
				#p = p + scale_shape_identity()
				
				if(is.null(shape.by))
				{
					info("Not scaling shape")
					#p = p + scale_shape_manual("Cluster", values=c(16))
					#p = p + scale_shape_manual("Cluster", values=c(15, 17, 16, 18))
				}else{
					info("Scaling shape manual")
					p = p + scale_shape_manual(legend.title.shape, values=c(15, 19, 17, 18))
				}

				if(numeric.data)
				{
					info(sprintf("Colouring by %s (numeric data)", colour.by))
					if(is.null(cols)){cols = material.heat(20)}#colorRamps::matlab.like(20)}
					if(min.zero)
					{
						lims = c(0, max(unlist(df[, colour.by])))
					}else
					{
						lims = range(unlist(df[, colour.by]))
					}
					info("Drawing colour bar")
					
					if(black.rim.pts){
						p = p + scale_fill_gradientn(legend.title, colours=cols, guide=guide_colorbar(nbin=100, raster=T, title.vjust=legend.vjust), limits=lims)
					}else{
						#p = p + scale_colour_gradientn(legend.title, colours=cols, guide="colorbar", limits=lims)
						p = p + scale_colour_gradientn(legend.title, colours=cols, guide=guide_colorbar(nbin=100, raster=T, title.vjust=legend.vjust), limits=lims) 
					}
					
					if(scale.alpha)
					{
						info("Scaling alpha")

						p = p + scale_alpha_continuous(range=c(0.2, 1))
					}
					p = p + guides(alpha=FALSE)
					
				}else
				{	
					if(is.null(highlight))
					{
						n.values = length(unique(unlist(df[,colour.by])))
						info(sprintf("Colouring by %s (non-numeric data. %s discrete values)", colour.by,n.values))
						if(is.null(cols) | length(cols) < n.values)
						{
							#warn(sprintf("Only have %s colours, but there are %s values, switching to intense palette", length(cols), n.values))
							cols = default.cols(n.values)
						}
					}else{
						info("Setting manual alpha scaling for highlight")
						cols = c(highlight.cols, "grey83")
						#
						print(cols)
						print(levels(df[, colour.by]))
						print(table(df[, colour.by]))
						# print(cols)
						# print(alphas)
						if(scale.alpha){
							alphas = c(rep(1, length(highlight)), 0.6)
							p = p + scale_alpha_manual(values=alphas)
						}
						
					}
					
					if(black.rim.pts){
						p = p + scale_fill_manual(legend.title, values=cols) + guides(alpha=FALSE) 
					}else{
						p = p + scale_colour_manual(legend.title, values=cols) + guides(alpha=FALSE) 
					}
				}
			
				#p = p + geom_point(aes_string(colour=colour.by.colname)) + scale_colour_gradientn(colours=cols) 
				if(density)
				{
					warn(sprintf("Cannot draw density and colour by %s!", colour.by))
					density=F
				}
					
				if(show.legend)
				{
					if(black.rim.pts)
					{
						if(!numeric.data){
						p = p + guides(fill = guide_legend(override.aes = list(shape = 21, size=point.size), title=legend.title)) + 
							guides(colour = guide_legend(override.aes = list(size=3)))
						}
					}else{
						p = p + guides(title=legend.title)
						if(!numeric.data){
							p = p + guides(colour = guide_legend(override.aes = list(size=legend.point.size)))
						}
					}
					
				}else{
					p = p + guides(colour=FALSE) + theme(legend.position = "none") 
				}
		}else
		{
			if(size.by=="point.size")
			{
				p = p + geom_point(size=point.size) + scale_size_continuous(size.by)
			}else{
				p = p + geom_point(aes_string(size=size.by)) + scale_size_continuous(size.by)
			}
			
		}
		
		if(density)
		{
			# k <- kde2d(x, y, n=200)
			# image(k, col=cols, alpha=0.2)
			if(density.col.pts)
			{
				use.kde2d=FALSE
				x = unlist(x)
				y = unlist(y)
				if(use.kde2d)
				{
					library(MASS)
					info("Calculating 2D-density estimate using kde2d")
					
					# Calculate 2d density over a grid
					#print(head(x))
					#print(head(y))
					
					dens <- kde2d(x, y, n = 200) 
					# create a new data frame of that 2d density grid
					# (needs checking that I haven't stuffed up the order here of z?)
					gr <- data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z))
					names(gr) <- c("xgr", "ygr", "zgr")
					# Fit a model
					mod <- loess(zgr~xgr*ygr, data=gr)
					# Apply the model to the original data to estimate density at that point
					info("Predicting density value for each point")
					df[, "density"] <- predict(mod, newdata=data.frame(xgr=x, ygr=y))
					p = p + scale_colour_gradientn(colours=density.cols, guide=FALSE)
				}else{
					info("Calculating 2D-density estimate using densCols")
					df[, "density"] <- densCols(x, y, colramp=colorRampPalette(density.cols))
					p = p + geom_point(data=df, aes(color=density)) + scale_color_identity() 
				}
				
				#p = p + stat_density2d(aes(colour = ..density..), contour = FALSE, geom="point") 
			}else{
				if(contour)
				{
					p = p + stat_density2d(aes(fill=..level.., alpha=..level..), geom='polygon', colour='black', show.legend = FALSE) + 
					scale_fill_gradientn(colours=density.cols, guide=FALSE)
				}else{
					# fade with alpha..
					# p = p + stat_density2d(aes(fill = ..density..,alpha=..density..), contour=FALSE, show.legend = FALSE, geom="tile",) + 
					# scale_fill_gradientn(colours=density.cols, guide=FALSE)
					p = p + stat_density2d(aes(fill = ..density..), contour=FALSE, show.legend = FALSE, geom="tile",) + 
					scale_fill_gradientn(colours=density.cols, guide=FALSE)
				}
			}
		}

		if(fit.line)
		{
			p = p + stat_smooth(method = fit.method)
		}
		if(!is.null(split.by))
		{
			p = p + facet_wrap(as.formula(paste("~", split.by)))
		}

		if(!is.null(label))
		{
			info(sprintf("Labeling points by %s", label))
			p = p + geom_text(size=1, vjust=4)
		}

		if(!is.null(add.curve))
		{
			x = data.frame(add.curve$s, add.curve$lambda)
			colnames(x) = c("tSNE_1", "tSNE_2", "lambda")
			x = x[order(x$lambda),]
			#print(head(x))
			p = p + geom_path(data = x, aes(x=tSNE_1, y=tSNE_2, colour=lambda), size=2, lineend="round")
		}

		if(show.marginals){

			info("Drawing marginal densities")
			warn("This isn't finished.")
			#### marginal density plots [see 
			#### http://www.sthda.com/english/wiki/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page-r-software-and-data-visualization	]
			# Marginal density plot of x (top panel)
			xdensity <- ggplot(df, aes_string(x.var, colour=colour.by)) + 
				geom_density(alpha=.5) + 
			  	theme(legend.position = "none") + scale_colour_manual("", values=cols) 

			# Marginal density plot of y (right panel)
			ydensity <- ggplot(df, aes_string(y.var, colour=colour.by)) + coord_flip() + 
				geom_density(alpha=.5) + 
			  	theme(legend.position = "none")  + scale_colour_manual("", values=cols) 

			blankPlot <- ggplot(df, aes_string(colour=colour.by))+geom_blank(aes(1,1)) + scale_colour_manual(legend.title, values=cols) + guides(alpha=FALSE) +
				theme(
					plot.background = element_blank(), 
					panel.grid.major = element_blank(),
					panel.grid.minor = element_blank(), 
					panel.border = element_blank(),
					panel.background = element_blank(),
					axis.title.x = element_blank(),
					axis.title.y = element_blank(),
					axis.text.x = element_blank(), 
					axis.text.y = element_blank(),
					axis.ticks = element_blank(),
					axis.line = element_blank()
				)

			library("cowplot")
			# grid.arrange(xdensity, blankPlot, p, ydensity, 
			# 	ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))

			# ggdraw() +
			#   draw_plot(xdensity, 0, .5, 1, .5) +
			#   draw_plot(p, 0, 0, .5, .5) +
			#   draw_plot(ydensity, .5, 0, .5, .5) +
			#   draw_plot_label(c("A", "B", "C"), c(0, 0, 0.5), c(1, 0.5, 0.5), size = 15)

			p = plot_grid(xdensity, NULL, p, ydensity, labels = c("A", "B", "C", "D"), ncol = 2)
			
			#return(p)
		}

		if(ret.data){
			return(df)
		}

		if(is.null(pdf))
		{
			if(ret.plot){
				return(p)
			}else{
				print(p)
			}
		}else{
			ggsave(p, filename=pdf, width=width, height=height)
		}
		
	}
	
}

# a more flexible and interactive version of genePlot
# in Seurat using ggvis
gene.plot <- function(seurat.obj, gene1, gene2, gene3=NULL, fit.line=F, id.use="orig.ident")
{
	if(is.null(gene3))
	{
		#g = data.frame(t(seurat.obj@data[c(gene1, gene2),]))
		g = fetch.data(seurat.obj, c(gene1, gene2, id.use))
	}else{
		g = fetch.data(seurat.obj, c(gene1, gene2, gene3, id.use))
	}
	g["Name"] = rownames(g)
	
	g = cbind(g, seurat.obj@data.info)
	
	if(!is.null(gene3))
	{
		library("lattice")
		#print(head(g))
		cloud(gene3~gene1*gene2,data=g)
	}else
	{
		g$id <- 1:nrow(g)  # Add an id column to use ask the key

		all_values <- function(x) {
		  if(is.null(x)) return(NULL)
		  row <- g[g$id == x$id, ]
		  paste0(names(row), ": ", format(row), collapse = "<br />")
		}

		rownames(g) = NULL
		#print(head(g))
		title = paste("Plotting ", gene1, "vs", gene2, "  (Log2 Counts, R =", round(cor(g[gene1], g[gene2]), digits=2), ")")
		print(title)
		if(fit.line)
		{
			g %>% ggvis(prop("x", as.name(gene1)), prop("y", as.name(gene2))) %>% 
			 layer_smooths(method = "lm", se = TRUE) %>% layer_points(key := ~id) %>% add_tooltip(all_values, "hover") %>% add_title(title = title)  
		}else{
			#fill=~Type
			g %>% ggvis(prop("x", as.name(gene1)), prop("y", as.name(gene2)), prop("fill", as.name(id.use)), key := ~id) %>%  
			 layer_points() %>% add_tooltip(all_values, "hover") %>% add_title(title = title, x_lab=gene1) # %>% add_axis("x", title=paste(gene1, "(Log2 count)"))   %>%
			 #add_axis("y",title=paste(gene2, "(Log2 count)"))
		}
	}	
}


violin_vertical <- function(seurat.obj, features.plot, 
	stretch.first=2, 
	cex.groups=12, 
	cex.genes=18, 
	group.by=NULL, 
	groups.only=NULL,
	ylab.max=12, ymax.use=NULL, 
	do.ret=FALSE, 
	sort.by=NULL, 
	sort.by.decreasing=T, 
	sort.fn=median, 
	use.imputed=FALSE,
	adjust.use=1, 
	boxplot=F, 
	show.axis.scale=T,
	border=F, 
	border.col="grey22", 
	border.width=1, 
	summary.show=T,
	summary.fn=median, 
	boxplot.width=0.9, 
	summary.bar.thick.inner=1, 
	summary.bar.thick.outer=2.5, 
	summary.bar.width.inner=0.35, 
	summary.bar.width.outer=0.5, 
	summary.bar.colour.inner="white", 
	summary.bar.colour.outer="grey10", 
	ret.plot=F,
	point.size=1,cols.use=NULL,do.jitter=FALSE,use.raw=FALSE,...) {
    
    if(!is.null(group.by)){
    	info(sprintf("Grouping by %s", group.by))
    	seurat.obj@ident = unlist(fetch.data(seurat.obj, group.by))
    	#print(table(seurat.obj@ident))
    	
    }  
    groups = seurat.obj@ident

    if(!is.null(groups.only))
    {
		keep.cells = which(seurat.obj@ident %in% groups.only)	
		if(length(keep.cells)<length(seurat.obj@cell.names)){
			seurat.obj = subsetSeuratObj(seurat.obj, cells.use=seurat.obj@cell.names[keep.cells], 
    		run.tsne=F, get.variable.genes=F)
		}
    	
    	seurat.obj@ident = groups[keep.cells]
    	seurat.obj@ident = factor(seurat.obj@ident, levels=rev(groups.only))

    	info(sprintf("Restricting plot to the %s cells in %s", length(keep.cells), paste(groups.only, collapse=", ")))
    	#print(table(seurat.obj@ident))
    }

    n.found = sum(features.plot %in% rownames(seurat.obj@data))
    info(sprintf("Found %s of the %s features.plot in rownames of count data", n.found, length(features.plot)))
    # print(features.plot)
    # print(corner(seurat.obj@data))


    if(!is.null(sort.by))
    {
    	info(sprintf("Sorting by %s", sort.by))
    	sort.df = t(fetch.data(seurat.obj, c(sort.by)))
    	 #seurat.obj@data[c(sort.by),]
    	sort.order = order(group.means(sort.df, groups=seurat.obj@ident, fn=sort.fn)[sort.by,], decreasing=sort.by.decreasing)
    	seurat.obj@ident = factor(seurat.obj@ident, levels=levels(seurat.obj@ident)[sort.order])
    }

    nf = length(features.plot)
    #size.use = point.size
    print("Drawing first plot")
    plots = VlnPlot(seurat.obj, features.plot, do.ret=TRUE, size.title.use=cex.genes)

    if(border){b=element_rect(size=border.width, colour=border.col)}else{b=element_blank()}
    plots[[1]]=plots[[1]]+coord_flip()+ scale_y_continuous(expand = c(0,0.1))+ labs(x=NULL, y=NULL)+
      theme(axis.text.y = element_text(angle=0, size=cex.groups),
            legend.position="none", panel.border=b, panel.margin=unit(c(0,0,0,0),"cm"), axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"))
    if(!show.axis.scale){
    	plots[[1]]=plots[[1]]+ theme(line = element_blank(), axis.text.x = element_blank())
    }
    if(!is.null(cols.use))
    {
    	plots[[1]] = plots[[1]] + scale_fill_manual(values=rev(cols.use))
    }
    if(summary.show){
    	plots[[1]] = plots[[1]] + stat_summary(geom = "crossbar", show.legend=FALSE, width=summary.bar.width.outer, fatten=summary.bar.thick.outer, color=summary.bar.colour.outer, 
						fun.data = function(x){ return(c(y=summary.fn(x), ymin=summary.fn(x), ymax=summary.fn(x))) }) +
						stat_summary(geom = "crossbar", show.legend=FALSE, width=summary.bar.width.inner, fatten=summary.bar.thick.inner, color=summary.bar.colour.inner, 
						fun.data = function(x){ return(c(y=summary.fn(x), ymin=summary.fn(x), ymax=summary.fn(x))) })
    }
    for(i in 2:nf){

		if(border){b=element_rect(size=border.width, colour=border.col)}else{b=element_blank()}
		plots[[i]]=plots[[i]]+coord_flip()+ scale_y_continuous(expand = c(0,0.1))+ labs(x=NULL, y=NULL)+
		theme(legend.position="none",axis.text.y = element_blank(), #plot.margin = unit(c(0,0,0,0),"cm"),
		      panel.border=b, panel.margin=unit(c(0,0,0,0),"cm"), axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"))
		if(!show.axis.scale){
    		plots[[i]]=plots[[i]]+ theme(line = element_blank(), axis.text.x = element_blank())
    	}
	    if(!is.null(cols.use))
	    {
	    	plots[[i]] = plots[[i]] + scale_fill_manual(values=rev(cols.use))
	    }
	    if(summary.show){
    	plots[[i]] = plots[[i]] + stat_summary(geom = "crossbar", show.legend=FALSE, width=summary.bar.width.outer, fatten=summary.bar.thick.outer, color=summary.bar.colour.outer, 
						fun.data = function(x){ return(c(y=summary.fn(x), ymin=summary.fn(x), ymax=summary.fn(x))) }) +
						stat_summary(geom = "crossbar", show.legend=FALSE, width=summary.bar.width.inner, fatten=summary.bar.thick.inner, color=summary.bar.colour.inner, 
						fun.data = function(x){ return(c(y=summary.fn(x), ymin=summary.fn(x), ymax=summary.fn(x))) })
    	}
    }

    library(cowplot)
    print("Arranging plots (cowplot)..")
    rval = plot_grid(plotlist=plots, align = 'h', nrow=1, rel_widths = c(stretch.first, rep(1, nf)))
    if(ret.plot){
    	return(rval)
    }else{
    	print(rval)
    }
    
}

many_bars <- function(seurat.obj, features.plot)
{
    plots = list()
    for(i in 1:length(features.plot)){
      f = features.plot[i]
      info(sprintf("Plotting %s ", f))
      p = bar(seurat.obj, f, dodgeplot = F, do.ret=T)
      plots[[i]]=p #+ scale_y_continuous(expand = c(0,0.1))+ labs(x=NULL, y=NULL)+
        # theme(line = element_blank(),axis.text = element_blank(),legend.position="none",plot.margin = unit(c(0,0,0,0),"cm"),
        #       panel.border=element_rect(size=2), panel.margin=unit(c(0,0,0,0),"cm"))
    }
    do.call("grid.arrange",c(plots, nrow=length(plots)))
}


histogram <- function(seurat.obj, features.plot, 
						group.by=NULL, 
						facet=T,
						ret.plot=F,
						xlab=NULL, 
						ylab=NULL,
						density=T, 
						cex.all=16,
						cex.strip=12,
						cumulative=F,
						scaled=F,
						use.cols=NULL, 
						binwidth=NULL, 
						interactive=F,
						width=9, 
						height=7, 
						pdf=NULL, 
						fill=T)
{
	info(sprintf("Generating histogram for %s", paste(features.plot, collapse=", ")))
	df = Seurat::fetch.data(seurat.obj, features.plot)
	colnames(df) = gsub(".", "_", colnames(df), fixed=T) # for some reason periods fuck things up
	genes = gsub(".", "_", features.plot, fixed=T)
	if(is.null(xlab)){xlab=features.plot}
	if(!is.null(group.by))
	{
		df$group = as.factor(unlist(fetch.data(seurat.obj, group.by)))
	}else{
		df$group = "All Cells"
	}
	ncols = length(unique(df$group))
	if(is.null(use.cols))
	{
		use.cols = default.cols(ncols) 
	}
	dm = melt(df, id.vars="group")
	print(head(dm))
	p = ggplot(dm, aes(x=value, colour=group, fill=group))  + theme_bw() + 
			theme(axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"),
			 axis.line.y = element_line(color="gray15", size = 0.75, lineend="round")) + 
			theme(
				panel.border= element_blank(),
		    	strip.background = element_blank(), 
		    	strip.text.x = element_text(size = cex.strip, face="bold"),
		    	axis.title.x = element_text(vjust=-0.5),
		    	#plot.margin = unit(c(1.2,1.2,1.2,1.2), "cm"),
		    	axis.text.x = element_text(angle = 45, hjust = 1),
		    	text = element_text(size=cex.all, colour="gray22")) + xlab(xlab) + ylab(ylab)
	if(facet){
		p = p + facet_wrap(variable~group, scale="free") 
	}
			

			#scale_fill_manual(group.by, values=use.cols) + 
	p = p + guides(alpha=FALSE) 
	if(fill)
	{
		if(density)
		{
			p = p + 
				geom_density(aes(x=value, fill=group), color = "NA", alpha = 0.3) +  
    			geom_density(aes(x=value, colour=group), fill = NA) +
			 	scale_fill_manual(group.by, values=use.cols) + 
				scale_colour_manual(group.by, values=use.cols) #values=rep("black", ncols)) 
			if(cumulative){
				warn("Density is ON, ignoring cumulative argument.")
			}
		}else{
			if(is.null(binwidth))
			{
				binwidth = (max(dm$value) - min(dm$value)) / 30
			}
			if(cumulative)
			{
				if(!scaled)
				{
					p = p + geom_histogram(binwidth=binwidth, aes(y=cumsum(..count..)))
				}else{
					p = p + stat_ecdf()
				}
				p = p + facet_wrap(variable~group, scale="free") +
						scale_fill_manual(group.by, values=use.cols) + 
						scale_colour_manual(group.by, values=rep("black", ncols)) 		
  				# stat_bin(aes(y=cumsum(..count..)),geom="line",color="green")
			}else{
				p = p + geom_histogram(binwidth=binwidth) + facet_wrap(variable~group, scale="free") + 
						scale_fill_manual(group.by, values=use.cols) + 
						scale_colour_manual(group.by, values=rep("black", ncols)) 	
			}
		}
	}else{
		warn("if fill argument is FALSE, code is unfinished.")
		if(density)
		{
			p = p + geom_line(lineend="round", stat="density",size=3) + 
				scale_colour_manual(group.by, values=use.cols) + 
				scale_fill_manual(group.by, values=rep(NA, ncols)) 
		}else{
			if(is.null(binwidth))
			{
				binwidth = (max(dm$value) - min(dm$value)) / 30
			}
			p = p + geom_histogram(binwidth=binwidth) 
		}
		
	}
	
	if(!is.null(pdf))
	{
		ggsave(p, file=pdf, width=width, height=height)
	}else
	{
		if(ret.plot){
			return(p)
		}else{
			print(p)
		}
		
		if(interactive)
		{
			library(plotly)
			ggplotly()
		}
	}

}

# draw a stacked bar plot with standard error bars (95% CI)
# assumes there's a groups column, and we want to compare the
# averages of 'genes' in 'groups'
# provide counts data frame and groups vector (cluster assignments probably)
bar <- function(seurat.obj=NULL, features.plot, 
					df = NULL,
					aggregate.by=mean,
					scale.to.group=NULL, #scale values to this group 
					group.by=NULL,
					facet.wrap.by=NULL,
					stack.by.gene=T,
					show.value=T,
					sort=F,
					sort.fn=mean,
					cex.label=3,
					cex.all=12,
					cex.strip=12,
					error.bar.width=0.25,
					line.thickness=1,
					use.cols=NULL,
					use.celltype.cols=F, # standardise colour allocation of known celltypes e.g. Goblet, Paneth
					legend.title=NULL,
					pdf=NULL,
					width=12,
					height=7,
					xlab=NULL, 
					ylab=NULL,
					borders=T,
					dodgeplot=T,
					do.ret=F,
					title=NULL,
					error="CI.95",  # can also be "SEM", which is 95CI/1.96
					#space.between.stacks=1,
					show.more.than.max.factor=1.4,
					groups.only=NULL)#
{
	#df = t(counts[rownames(counts) %in% genes, ])
	if(!is.null(seurat.obj))
	{
		df = Seurat::fetch.data(seurat.obj, vars.all=features.plot)
		if(is.null(group.by))
		{
			groups = unlist(fetch.data(seurat.obj, group.by))
		}else{
			groups = unlist(seurat.obj@ident)
		}
		
	}else{
		if(!is.null(df))
		{
			info(sprintf("No seurat object provided. Using the given dataframe (dims %s)", paste(dim(df), collapse="x")))
			#print(head(df))
			groups = unlist(as.character(df[,group.by]))
			info(sprintf("Groups [%s]:", group.by))
			print(table(groups))
			df = data.frame(df[,features.plot])
			colnames(df) = features.plot
			#print(head(df))
		}else
		{
			stop("Must provide non-null seurat object or data frame!")
		}
		
	}
	if(error=="SEM")
	{
		info("Displaying standard error of the mean (SEM)")
	}
	
	
	if(!is.null(scale.to.group))
	{
		scale.to.group = gsub("_", " ", scale.to.group, fixed=T)
	}

	if(is.null(groups.only))
	{
		#print(unique(groups))
		groups.only = sort(as.character(unique(groups)))
	}else{

		info("Adjusting groups")
		groups = groups[groups %in% groups.only]
	}

	if(sort)
	{
		f = features.plot[1]
		info(sprintf("Sorting by %s of %s", as.character(substitute(sort.fn)), f))
		x = unlist(df[, f])
		y = aggregate(x~groups, FUN=sort.fn)
		colnames(y) = c("Group", "Mean")
		
		y = y[order(y$Mean, decreasing=T),]
		print(y)
		new.levels = unname(unlist(y$Group))
		print("New levels:")
		print(new.levels)
		groups = factor(groups, levels=new.levels)
	}else{
		new.levels = groups.only
	}

	if(is.null(use.cols))
	{
		if(stack.by.gene)
		{
			ncols = length(unique(groups.only))
			use.cols = default.cols(ncols) #
		}else
		{
			ncols = length(unique(genes))

		}
		use.cols = default.cols(ncols)
	}	

	## If need to match colours with other plots (overrides specified colours)
	if(use.celltype.cols)
	{
		if(!is.null(use.cols))
		{
			warn("Overriding provided colours and matching with known cell types")
		}
		info("Using standardised cell-type colours")
		use.cols = get.celltype.cols(groups)
	}

	# clean up some strings
	info("Cleaning up column names")
	colnames(df) = gsub(".", "_", colnames(df), fixed=T) # for some reason periods fuck things up
	print(colnames(df))
	genes = gsub(".", "_", features.plot, fixed=T)
	groups = gsub("_", " ", groups, fixed=T)



	agg.mean <- aggregate(df, by=list(groups), FUN=aggregate.by, na.rm=TRUE)
	agg.mean = agg.mean[agg.mean$Group.1 %in% groups.only,]
	rownames(agg.mean) = agg.mean$Group.1
	agg.mean$groups=NULL
	agg.mean$Group.1 = NULL
	colnames(agg.mean) = paste(colnames(agg.mean), "mean", sep=".")


	agg.sd <- aggregate(df, by=list(groups), FUN=sd, na.rm=TRUE)
	agg.sd = agg.sd[agg.sd$Group.1 %in% groups.only,]
	rownames(agg.sd) = agg.sd$Group.1
	agg.sd$Group.1 = NULL
	agg.sd$groups=NULL
	colnames(agg.sd) = paste(colnames(agg.sd), "SD", sep=".")

	data = merge(agg.sd, agg.mean, by="row.names")
	#print(data)

	n = rep(0, length(groups.only))
	CI.95 = rep(0, length(groups.only))
	for(i in 1:nrow(data))
	{
		n[i] = length(which(groups==data$Row.names[i]))
	}


	for(i in 1:length(genes))
	{
		s = data[,paste(genes[i], "SD", sep=".")]
		if(error=="CI.95")
		{
			ci = qnorm(0.975)*s/sqrt(n) # qnorm(0.975) is ~1.96
		}else{
			if(error=="SEM")
			{
				ci = s/sqrt(n)
			}else{
				stop("Error argument must be 'SEM' or 'CI.95'")
			}
		}
		
		# cat(sprintf("ci's over groups for %s\n", genes[i]))
		# print(ci)
		data[paste(genes[i],"CI_95", sep=".")] = ci
	}
	
	# get the data in the right format!
	data = data[, -grep(".SD", colnames(data))]	
	rownames(data) = NULL
	print(data)
	data = rename(data, c("Row.names"="group"))
	data$group = factor(data$group, levels=new.levels)
	#print(levels(groups))
	data = melt(data, id.vars="group")
	
	#split the variable column (have to escape the period, its a regex wildcard):
	data = data.frame(data)
	data$measure <- lapply(strsplit(as.character(data$variable), "\\."), "[", 2)
	data$gene <- lapply(strsplit(as.character(data$variable), "\\."), "[", 1)
	data$variable = NULL

	#dcast will fail if columns aren't factors
	data$measure = as.factor(unlist(data$measure))
	data$gene = as.factor(unlist(data$gene))
	data$group = as.factor(unlist(data$group))
	data <- dcast(data, group + gene  ~ measure)#Mean + CI_95  ~ measure)
	
	# preserve plotting order
	data$gene = factor(data$gene, levels=features.plot)

	#print(data)
	if(!is.null(scale.to.group))
	{
		
		# if(is.null(scale.to.group))
		# {
		# 	scale.to.group = groups.only[i]
		# }
		info("Calculating relative quantification")
		info(sprintf("Scaling to %s", scale.to.group))
		if(!scale.to.group %in% data$group)
		{
			stop(sprintf("%s not present in given groups!", scale.to.group))
		}
		d.baseline = data[data$group==scale.to.group,]
		
		#print(d.baseline)
		#scale = d.baseline$mean # vector of the baseline values for each gene.
		for(i in 1:length(genes))
		{
			s = d.baseline[d.baseline$gene == genes[i], ]$mean # careful, s is a baseline, not a scale factor
			# print("Baseline:")
			# print(d.baseline)
			# print("Before scaling:")
			# print(data[data$gene==genes[i], ])

			
			if(s == 0)
			{
				warn(sprintf("Could not scale to %s, %s is not expressed!", scale.to.group, genes[i]))
			}else{
				info(sprintf("Scaling %s [BaseLine (%s)=%s]", genes[i], scale.to.group, s))
				data[data$gene==genes[i], ]$mean = data[data$gene==genes[i], ]$mean / s
				data[data$gene==genes[i], ]$CI_95 = data[data$gene==genes[i], ]$CI_95 / s
			}
			

			#print("After scaling:")
			# print(data[data$gene==genes[i], ])
		}
		y.lab = "Relative Quantification\n"
	}else
	{
		y.lab = "Log2 TPM + 1\n" 
	}
	if(!is.null(ylab))
	{
		y.lab = ylab
	}
	#print(data)	
	
	info("Building plot..")

	#plot params
	if(max(data$mean) < 10){breaks=0.5}else{breaks=1}
	dodge <- position_dodge(width=0.9)
	limits <- aes(ymax = mean + CI_95, ymin= ymin)
	# make sure error bars don't go through the floor
	data$ymin = data$mean - data$CI_95
	data$ymin[data$ymin < 0] <- 0

	#setup ggplot object
	if(stack.by.gene)
	{
		p <- ggplot(data, aes(fill=group, y=mean, x=gene))# width param provides more spacing between stacked bars for clarity  
	}else
	{
		p <- ggplot(data, aes(fill=gene, y=mean, x=group))
	}	

	#plot the error bars first, to get UP half only
	if(dodgeplot)
	{
		p=p + geom_errorbar(limits, position=dodge, width=error.bar.width, size=line.thickness) + theme_bw() +
			ylab(y.lab) + xlab(xlab) + scale_fill_manual(values=use.cols, name=legend.title) 
	}else
	{	
		p=p + geom_errorbar(limits, width=error.bar.width, size=line.thickness) + theme_bw() +
			ylab(y.lab) + xlab(xlab) + scale_fill_manual(values=use.cols, name=legend.title) 
	}

	if(dodgeplot)
	{
		p = p + geom_bar(position="dodge", stat="identity", colour="black", size=line.thickness) + 
				scale_y_continuous(breaks=seq(0, max(data$mean) * show.more.than.max.factor, breaks), limits=c(0, max(data$mean) * show.more.than.max.factor))
	}else{
		p = p + geom_bar(stat="identity", colour="black", size=line.thickness) + 
				scale_y_continuous(breaks=seq(0, max(data$mean) * show.more.than.max.factor, breaks), limits=c(0, max(data$mean) * show.more.than.max.factor))
	}
	
	if(is.null(legend.title))
	{
		if(stack.by.gene)
		{
			legend.title = "Cell Type"
		}else
		{
			legend.title = "Gene"
		}	
	}
	
	# info(sprintf("Using %s colours..", length(unique(use.cols))))
	#p=p + geom_bar(position=dodge) + geom_errorbar(limits, position=dodge, width=0.25) + theme_bw() +
	
	

	if(!borders)
	{
		p = p + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
			panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks=element_line(colour="black"))
	}
	
	if(show.value & dodgeplot)
	{
		p = p + geom_text(aes(label = round(mean, 2), y = mean + (1.5*CI_95)), position=dodge, size = cex.label, vjust=-1)
	}
	# p = p + theme(axis.title.y = element_text(vjust=2.5), # move the y title away from the plot a bit
	# 		axis.text.x = element_text(angle = 45, hjust = 1, size=12), 
	# 		axis.text.y = element_text(size=12), 
	# 		plot.margin = unit(c(1.2,1.2,1.2,1.2), "cm"), text = element_text(size=text.size)) + 
	# 		theme(axis.line = element_line(colour = "black"),
	# 	    panel.grid.major = element_blank(),
	# 	    panel.grid.minor = element_blank(),
	# 	    panel.border = element_blank(),
	# 	    panel.background = element_blank(), 
	# 	    text = element_text(size=20, colour="gray22", face="bold"))  


    p = p + theme_bw() + 
		theme(axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"),
			 axis.line.y = element_line(color="gray15", size = 0.75, lineend="round")) + 
		theme(
			panel.border= element_blank(),
	    	strip.background = element_blank(), 
	    	strip.text.x = element_text(size = cex.strip, face="bold"),
	    	axis.title.x = element_text(vjust=-0.5),
	    	plot.margin = unit(c(1.2,1.2,1.2,1.2), "cm"),
	    	axis.text.x = element_text(angle = 45, hjust = 1),
	    	text = element_text(size=cex.all, colour="gray22")) 

    if(!is.null(title))
    {
    	p = p + ggtitle(title)
    }

	if(do.ret)
	{
		return(p)
	}
	if(is.null(pdf))
	{
		print(p)
	}else{
		ggsave(p, filename=pdf, width=width, height=height)
	}
	return(p)

}

#data.frame(spline.poly(y[chull(y),], 1000))

#convex hulls for highlighting points
#-----------------------------------

StatChull <- ggproto("StatChull", Stat,
  compute_group = function(data, scales) {
    # cat("\n\n\n\n")
    # info("Computing convex hull")
    # print(data)
    # info("Data dims:")
    # print(dim(data))
    # if(all(df$colour) == NA)
    # {
    # 	return (NULL)
    # }
    hull = data[chull(data$x, data$y), , drop = FALSE]
	# print("hull")
	# print(head(hull))    
    
    smooth.hull = data.frame(spline.poly(hull[, c("x", "y")], 1000))
    colnames(smooth.hull) = c("x", "y")
    cval = unlist(hull$colour)[1]
    # print("smooth hull")
    # print(head(smooth.hull))
    # info(sprintf("Setting colour to %s ", cval))
    smooth.hull$colour = rep(cval, nrow(smooth.hull))
    smooth.hull$PANEL = rep(unlist(hull$PANEL)[1], nrow(smooth.hull))
    smooth.hull$group = rep(unlist(hull$group)[1], nrow(smooth.hull))
    return (smooth.hull)
  },
  required_aes = c("x", "y")
)

stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


# Splining a polygon.
#
#   The rows of 'xy' give coordinates of the boundary vertices, in order.
#   'vertices' is the number of spline vertices to create.
#              (Not all are used: some are clipped from the ends.)
#   'k' is the number of points to wrap around the ends to obtain
#       a smooth periodic spline.
#
#   Returns an array of points. 
# 
spline.poly <- function(xy, vertices, k=5, ...) {
    # Assert: xy is an n by 2 matrix with n >= k.
    # Wrap k vertices around each end.
    n <- dim(xy)[1]
    if (k >= 1) {
        data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
    } else {
        data <- xy
    }
    # Spline the x and y coordinates.
    data.spline <- spline(1:(n+2*k), data[,1], n=vertices, ...)
    x <- data.spline$x
    x1 <- data.spline$y
    x2 <- spline(1:(n+2*k), data[,2], n=vertices, ...)$y
    # Retain only the middle part.
    cbind(x1, x2)[k < x & x <= n+k, ]
}

### plot the output of a linear regression.
### downloaded from 
### https://susanejohnston.wordpress.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/
ggplotRegression <- function (fit) {
	require(ggplot2)

	ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
	  geom_point() +
	  stat_smooth(method = "lm", col = "red") +
	  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
	                     "Intercept =",signif(fit$coef[[1]],5 ),
	                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))
}

# draw plots to help filter low-quality cells in scRNA-seq data
filter_plots = function(counts, gene_min = 1000, gene_max = 0, mito_max = 1)
{
		mito.genes <- grep(pattern = "^mt-", x = rownames(counts), value = TRUE, ignore.case = T)
		info("Computing mitochondrial read content")
		percent.mito <- Matrix::colSums(counts[mito.genes, ])/Matrix::colSums(counts)
		info("Computing genes detected per cell")
		nGene = Matrix::colSums(counts > 0)
		if(gene_max == 0){gene_max = max(nGene) + 100}
		x = data.frame(percent.mito, nGene)

		info("Filtering")
		x$pass_filter = x$nGene > gene_min & x$nGene < gene_max & x$percent.mito < mito_max
		info("Drawing plots")
		ggsave(ggplot(x, aes(x=nGene)) + geom_histogram(fill="dodgerblue") + geom_vline(xintercept = gene_min, color="red") + theme_cowplot() + 
			geom_vline(xintercept = gene_max, color="red"), filename="filter_nGene.pdf")

		ggsave(ggplot(x, aes(x=percent.mito)) + geom_histogram(fill="dodgerblue") +  theme_cowplot() + scale_x_continuous(labels=scales::percent) + 
			geom_vline(xintercept = mito_max, color="red"), filename="filter_PercentMito.pdf")

		ggsave(ggplot(x, aes(x=percent.mito, y=nGene, color=pass_filter)) + geom_point(stroke=0, size=0.5) +  theme_cowplot() + 
			geom_vline(xintercept = mito_max, color="red") + geom_hline(yintercept = gene_min, color="red") + 
			geom_hline(yintercept = gene_max, color="red") + scale_color_manual("Pass filtering", values=c("grey75", "green")) +  scale_x_continuous(labels=scales::percent) + 
			guides(colour = guide_legend(override.aes = list(size=6))) + stat_density2d(color="black", size=0.5), filename="filter_both.pdf")

}




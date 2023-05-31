library(crayon)
#library(MASS)
library(RColorBrewer)
library(futile.logger)
library(plyr)
#library(data.table)
library(ggsignif)

LOGGING_TO <<- "console"

# script <- function(name)
# {
# 	homedir = path.expand("~")
# 	path = sprintf("../../../dev/adam/rna_seq/r/%s.r", homedir, name)
# 	# info(sprintf("Sourcing script: %s", path))
# 	source(path)
# }

discretise <- function(x, num.bins=4)
{
	library(Hmisc)
	info("Discretising")
	rval = t(apply(x, 1, function(a){as.integer(cut2(a, g=num.bins))-1}))
	colnames(rval) = colnames(x)
	rval
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))} # scale a range to 0-1

# draw a bar plot of the cluster labels of unsorted (bulk) cells
bulk.cell.ids <- function(seurat.obj, clusters)
{
	if(length(grep("BULK", seurat.obj@data.info$orig.ident)) > 0)
	{
		cat("Drawing bar plot of cluster labels of bulk cells\n")
		bulk.cells = whic
		h(seurat.obj@data.info$orig.ident == "BULK")
		x = data.frame(table(clusters[bulk.cells]))
		x$Freq.P = x$Freq / length(clusters[bulk.cells])
		g=ggplot(x, aes(x = reorder(Var1, -Freq.P), y=100*Freq.P, fill = factor(Var1))) + 
			geom_bar(stat="identity") + 
			scale_fill_manual(values=intense.cols(length(unique(x$Var1)))) + 
			theme_bw() + ggtitle("Putative Cell Type IDs of 124 Bulk IECs") + 
			scale_y_continuous(breaks=seq(0, 35, by=5)) + xlab("") + ylab("Frequency")  + 
			theme(axis.text.x = element_text(angle = 45, hjust = 1)); 
		ggsave(g, file="Bulk.cell.type.frequencies.pdf", width=10, height=10)
	}else
	{
		cat("No BULK cells in object, cannot draw bar plot of cluster labels\n")
	}
	
}

log.to.file <- function(filename="default.log")
{
	if(is.na(filename) | is.null(filename) | filename=="")
	{
		log.to.console()
	}else{
		flog.warn("Switching to log file %s", filename)
		flog.appender(appender.file(filename))
		LOGGING_TO <<- filename
	}

}

log.to.console <- function()
{
	if(LOGGING_TO != "console"){flog.warn("Switching logging to console")} # write in the log file that we switched
	flog.appender(appender.console())
	flog.warn("Switching logging to console")
	LOGGING_TO <<- "console"
}

# reset logging variables when this script is sourced.
log.to.console()


log.to.console.and.file <- function(filename="default.log")
{
	flog.warn("Logging [TEE] to console and log file %s", filename)
	m = flog.appender(appender.tee(filename))
	LOGGING_TO <<- c("console", filename)
}

asNumeric <- function(x) as.numeric(levels(x))[x]#as.character(x))
factorsNumeric <- function(d) modifyList(d, lapply(d[, sapply(d, is.factor)], asNumeric))

am.logging.to.file <- function()
{
	if(length(LOGGING_TO)==1){
		return(LOGGING_TO != "console")
	}else{
		return(TRUE)
	}
}

# some logging utils.
info <- function(text, ...)
{
	#cat(sprintf(paste(Sys.time(),"INFO:", text,"\n")))
	if(is.data.frame(text))
	{
		if(am.logging.to.file()){
			#info("Logging data frame:") #debug
			if(length(LOGGING_TO)==1){
				log.file = LOGGING_TO
			}else{
				print(text)  ## TEE logging is on, write the dataframe to the console.
				log.file = LOGGING_TO[2]
			}
			max.print <- getOption('max.print')
			print(dim(text))
			print(log.file)
			options(max.print=nrow(text) * ncol(text))
			sink(file=log.file, append=T)
			print(text)
			sink()
			options(max.print=max.print)
		}else{
			print(text)
		}
	}else{
		flog.info(text)
	}
}

warn <- function(text, ...)
{
	flog.warn(crayon::yellow(text))
}

error <-function(text, ...)
{
	flog.error(crayon::red(text))
}

check_yes_no <- function()
{
    message("Continue? Y/N")
    x <- readLines(n = 1)
    return(toupper(x) %in% c("Y", "YES"))
}

alpha <- function(x, thresh=0){sum(x > thresh) / length(x)}

group.means <- function(counts, groups, fn=mean, use.data.table=F)
{
	if(use.data.table)
	{
		#info("Using data.table")
		dt <- data.table(t(counts), check.names=F)
		dt$group = groups
		r = data.frame(dt[, lapply(.SD, mean), by=group])
		rownames(r) = r$group
		r$group = NULL
		r = t(r)
	}else{
		#info("Using aggregate")
		counts <- aggregate(t(counts), by=list(groups), FUN=fn)
		rownames(counts) = counts$Group.1
		counts$Group.1 = NULL
		r = t(counts)
	}
	#info("Done")
	return(r)
}

source("color.r", chdir = T)
source("plot.r", chdir = T)

### given a dataframe full of strings (Gene ids), and a list of strings (gene ids),
### return a new dataframe preserving the column-wise order, but removing
### any ids that are not in the list.
prune.df.with.list <- function(df, l, remove.dups=F)
{
	cols = list()
	for(i in 1:ncol(df))
	{
		old.col = df[,i]
		new.col = as.character(old.col[toupper(old.col) %in% toupper(l)])
		cols[[i]] = new.col
	}
	names(cols) = colnames(df)

	### this block removes duplicates by placing each gene in the list in which it is ranked the highest.
	if(remove.dups){
		clean.cols = list()
		for(i in 1:length(cols))
		{
			old.col = cols[[i]]
			new.col = NULL
			for(j in 1:length(old.col[!is.na(old.col)]))
			{
				id = old.col[j]
				best.in.other = min(get.ranks(id, cols))
				rmv = best.in.other < j
				if(!rmv){
					if(is.null(new.col)){
						new.col = as.character(id)
					}else{
						new.col = c(new.col, as.character(id))
					}
				}
				warn(sprintf("%s is at rank %s in current list and highest %s in another list. Remove=%s", id, j, best.in.other, rmv))
			}
			print(new.col)
			clean.cols[[i]] = new.col
		}
		cols = clean.cols
	}
	n.obs <- sapply(cols, length)
	seq.max <- seq_len(max(n.obs))
	if(max(seq.max) < 2){
		error("Found:")
		print(cols)
		error("Not enough hits to make a dataframe.")
		return(NULL)
	}
	d = data.frame(sapply(cols, "[", i = seq.max))
	colnames(d) = colnames(df)
	return(d)
}

### used by prune.df.with.list. Return the ranks of an id in a list of lists.
### list1 = c("A", "B", "C", "D")
### list2 = c("B", "A", "D", "C")
### list3 = c("E", "C", "B", "F")
### m = list(list1,list2,list3)
### e.g get.ranks("A", ) --> 1, 2, integer(0)

get.ranks <- function(id, lists)
{	
	f = function(x){grep(id, x, ignore.case = T)}
	rval = unlist(lapply(lists, f))
	info(sprintf("get.ranks for %s --> %s", id, paste(rval, collapse=", ")))
	return(rval)
}


#wget http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt -O human_mouse_orthologs_JAX.txt
convert_mouse_to_human <- function(mouse_genes)
{
	info("Reading JAX ortholog table")
	table_file = "~/ref/tables/HOM_MouseHumanSequence.rpt"
	if(file.exists(table_file))
	    ortho = read.delim(table_file)
	else{
	    error("Couldn't load JAX table")
	    error("Looked here:")
	    error(table_file)
	    stop("Make sure broad server is mounted")
	}
	info("Matching IDs")
	ortho_mm = subset(ortho, Common.Organism.Name == "mouse, laboratory")
	ortho_hs = subset(ortho, Common.Organism.Name == "human")
	homol_ids = subset(ortho_mm, Symbol %in% mouse_genes)$HomoloGene.ID
	as.character(subset(ortho_hs, HomoloGene.ID %in% homol_ids)$Symbol)
}

#wget http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt -O human_mouse_orthologs_JAX.txt
convert_human_to_mouse <- function(human_genes)
{
	info("Reading JAX ortholog table")
	ortho = read.delim("~/ref/tables/HOM_MouseHumanSequence.rpt")
	info("Matching IDs")
	ortho_mm = subset(ortho, Common.Organism.Name == "mouse, laboratory")
	ortho_hs = subset(ortho, Common.Organism.Name == "human")
	homol_ids = ortho_hs[match(human_genes, ortho_hs$Symbol),]$HomoloGene.ID
	as.character(ortho_mm[match(homol_ids, ortho_mm$HomoloGene.ID),]$Symbol)
}

# equalize a grouping variable by randomly sampling
resample <- function(groups, max.group.size=100, return.index=F, verbose=T)
{
    keep = c()
    levs = levels(groups)
    if(is.null(levs)){levs = unique(groups)} # if groups is a character vector it may not have levels
    for(g in levs) # levels is preferable
    {
        in.group = which(groups==g)
        n.in.group = sum(groups==g)
        if(verbose) info(sprintf("%s in %s group", n.in.group, g))
        if(n.in.group > max.group.size)
        {
            if(verbose) info(sprintf("Resizing to %s", max.group.size))
            keep = c(keep, sample(in.group, max.group.size))
        }else{
            keep = c(keep, in.group)
        }
    }
    if(return.index){
      return(keep)
    }else{
      return(factor(groups[keep]))
    }
}

# function for creating breaks for a diverging color map
# that center on zero
diverging_breaks = function(x = c(rnorm(100, 2, 2), rnorm(100, -0.5, 0.5)), n=25, eps=0.01, center=0, center_region=NULL)
{
	#if((n %% 2) == 0){stop("n should be odd for symmetric color breaks")}
	if(min(x) >= 0){stop("diverging breaks should only be used with data containing positive and negative values")}
	if(max(x) <= 0){stop("diverging breaks should only be used with data containing positive and negative values")}
	k =  (n-1)/2
	pos = 1:k * max(x) / k #quantile(x[x < 0], k, probs=1:k / k)
	neg = k:1 * min(x) / k #quantile(x[x > 0], k, probs=1:k / k)
	#if(is.null(center_region)){center_region = mean()}
	pos[k] = pos[k] + eps #taking the exact max as the break limit means values are outside the limits
	neg[1] = neg[1] - eps
	c(neg, 0, pos)
}

# use the above function to generate the colors to match the breaks to colors 
# n - total colors required, must be even
# eps - extends the breaks this far past the maximum
# brewer_col_range - uses this range of brewer colors, 9 is most contrast, 1 is least (will). 
diverging_colors = function(x = c(rnorm(100, 2, 2), rnorm(100, -0.5, 0.5)), n=25, eps=0.01, center=0, center_col="white", palette = "RdBu", brewer_col_range = 9)
{
	if((n %% 2) != 0){stop("n should be even for symmetric color breaks")}
	if(min(x) >= 0){stop("diverging breaks should only be used with data containing positive and negative values")}
	if(max(x) <= 0){stop("diverging breaks should only be used with data containing positive and negative values")}
	brks = diverging_breaks(x, n=n+1)
	lookup = list("Spectral"=c("YlOrRd", "YlGnBu"),
					"RdYlGn"=c("Reds", "Greens"),
					"RdYlBu"=c("YlOrRd", "Blues"),
					"RdGy"=c("Reds", "Greys"),
					"RdBu"=c("Reds", "Blues"),
					"PuOr"=c("Oranges", "Purples"),
					"PRGn"=c("Purples", "Greens"),
					"PiYG"=c("RdPu", "Greens"),
					"BrBG"=c("YlOrBr", "PuBuGn"))
	if(!palette %in% names(lookup)){stop("Palette must be one of the diverging brewer palettes: Spectral, RdYlGn, RdYlBu, RdGy, RdBu, PuOr, PRGn, PiYG, BrBG")}
	cols = lookup[[palette]]
	k = floor(n/2) -1 
	list("colors"=c(rev(colorRampPalette(brewer.pal(9, cols[2])[1:brewer_col_range])(k)), rep(center_col, 2) , colorRampPalette(brewer.pal(9, cols[1])[1:brewer_col_range])(k)),
		"breaks"=brks)
}

# test whether the genes that are upregulated are enriched [hypergeomtric] for
# genes from the given signatures
# @param 'geneset' is any list of genes, e.g from a DE test.
# @param 'marker.genesets' is a named dataframe where the columns are known gene sets to test for enrichment, e.g. cell-cycle, apoptosis
marker.enrichment.test <- function(geneset, marker.genesets, n.total.genes, do.plot=T, verbose=T)
{
	library(cowplot)
	geneset = unlist(as.character(geneset))
	geneset = geneset[!is.na(geneset)]
	m = length(geneset)
	n = n.total.genes #sum(rowSums(counts)>10) # don't use all genes, not a realistic null

	nl = ncol(marker.genesets)
	
	if(verbose){info(sprintf("Testing %s DE genes against %s marker genesets", m, nl))}
	
	p = vector(mode = "numeric", length = nl)
	num.olaps = vector(mode = "integer", length = nl)
	genes = vector(mode = "character", length = nl)
	genes.in.sig = vector(mode = "character", length = nl)
	num.genes.in.sig = vector(mode = "integer", length = nl)
	cells = vector(mode = "character", length = nl)
	i = 1
	for(celltype in colnames(marker.genesets))
	{
		#https://www.biostars.org/p/15548/
		#phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
		#where q=size of overlap-1; m=number of upregulated genes in experiment #1; 
		#n=(total number of genes on platform-m); k=number of upregulated genes in experiment #2.
		mg = unique(unlist(unname(as.character(marker.genesets[, celltype]))))
		mg = mg[!is.na(mg)]
		olap = toupper(as.character(geneset)) %in% toupper(mg)
		olap.genes = geneset[olap]
		olap.str = paste(olap.genes, collapse=", ")
		mg.str = paste(mg, collapse=", ")
		no = length(olap.genes)	
		if(no>0)
		{
				if(verbose){info(sprintf("%s enriched genes were found in the \'%s\' marker set. [%s]", no, celltype, olap.str))}
		}
		
		q = sum(olap) -1
		k = length(mg)
		#info(sprintf("Running hypergeometric test. N=%s, overlap-1=%s, m=%s, k=%s", n, q,m, k))
		p[i]   = phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)

		#info(sprintf("Setting element %s of olap to %s", i, no))
		num.olaps[i]  = no
		num.genes.in.sig[i] = k
		genes.in.sig[i] = mg.str
		genes[i] = olap.str
		cells[i] = celltype
		i = i + 1
	}

	res = data.frame(cells, p, num.olaps, genes, num.genes.in.sig, genes.in.sig, m)
	colnames(res) = c("MarkerGeneSet", "p", "N. genes enriched", "Genes", "N. genes in sig.", "Genes in sig", "N.genes DE")

	if(do.plot)
	{
		g=ggplot(res, aes(x=MarkerGeneSet, y=-log10(p))) + 
			geom_bar(stat="identity") + xlab("") + geom_hline(yintercept=-log10(0.05), 
				linetype="dotted", color="firebrick2") + 
			theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
			ylab("Probability of enrichment [-log10(p)]")

	}
	return (list("results"=res, "plot"=g))
}

#https://stackoverflow.com/questions/15201305/how-to-convert-a-list-consisting-of-vector-of-different-lengths-to-a-usable-data
convert_list_to_df = function(x)
{
		n.obs <- sapply(x, length)
		seq.max <- seq_len(max(n.obs))
       data.frame(sapply(x, "[", i = seq.max))
}

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

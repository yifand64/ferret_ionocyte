
# Hack to fix problems setting up clusters in RStudio
if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" ) {
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
}

source("mast_single_cell.r")
# Store the results of SCDE comparisons
library(pbapply)
DE.Results <- setClass("DE.Results", slots = c(
 											sig="data.frame", 
 											sig.down="data.frame", 
 											sig.tfs="data.frame",
 											all="data.frame", 
 											min.p="numeric", 
 											error.models="data.frame", 
 											expr.prior="data.frame", 
 											mast.glm="ZlmFit",
 											min.lower.bound="numeric", 
 											name="character", 
 											groups="character", 
 											batch="character"))

setMethod("show", "DE.Results",
          function(object) {
            cat("An object of class ", class(object), "\n", sep = "")
            cat(" ", nrow(object@sig), " upregulated genes, and ",
                nrow(object@sig.down), " downregulated p<", object@min.p, "\n", sep = "")
            invisible(NULL)
          }
)

### compare two ECDFs (http://stackoverflow.com/questions/27233738/finding-location-of-maximum-d-statistic-from-ks-test)
compare <- function(x, y) {
  n <- length(x); m <- length(y)
  w <- c(x, y)
  o <- order(w)
  z <- cumsum(ifelse(o <= n, m, -n))
  i <- which.max(abs(z))
  w[o[i]]
}

# m is the output of get.markers
count.de <- function(m, fdr=0.05)
{
	n =ncol(m$pairs)
	info(sprintf("For %s comparisons", n))
	for(i in 1:n)
	{
		name = paste(m$pairs[,i], collapse="_vs_")
		d = d= m$diffs[[i]]@all
		num.de = nrow(subset(d, p.adj < fdr))
		if(num.de < 10) {warn(sprintf("Name: %s, num-DE: %s", name, num.de))}
		else{
			info(sprintf("Name: %s, num-DE: %s", name, num.de))
		}
	}
}


### tests if two bimodal distributions have siginificantly different proportional densities in the two modes (fisher's exact)
### note, that if the distributions are unimodal, the hartigan's dip-test p-value will be high and the test should be ignored.
do.bimodal.test <- function(x, y)
{
    library(diptest)
    library(mclust)
    dt.x = diptest::dip.test(x)
    dt.y = diptest::dip.test(y)
    
    if(dt.x$p.value > 0.8 | dt.y$p.value > 0.8){
        warn(sprintf("p-value for bimodality is %s for X and %s for Y [hartigan's dip-test]!", dt.x$p.value, dt.y$p.value))
    }else{
        info(sprintf("p-value for bimodality is %s for X and %s for Y [hartigan's dip-test]", dt.x$p.value, dt.y$p.value))
    }
    
    if(length(x)<10 | length(y)<10){
    	error(sprintf("Mclust sometimes fails on small samples. Don't use this test with very short vectors! Given lengths: %s, %s", length(x), length(y)))
    	return(NULL)
    }
    
    info("Fitting GMM for x")
    x.gmm = Mclust(x, G = 2)
    
    info("Fitting GMM for y")
    y.gmm = Mclust(y, G = 2)
    
    if(length(unique(y.gmm$classification))==1 | length(unique(x.gmm$classification))==1){
    	error("Mclust put everything in the same cluster! Test failed!")
    	return(NULL)
    }

    if(mean(y[y.gmm$classification==2]) > mean(y[y.gmm$classification==1]))
    {
        high.mode.y = 2
        
    }else{
        high.mode.y = 1
    }
    
    if(mean(x[x.gmm$classification==2]) > mean(x[x.gmm$classification==1]))
    {
        high.mode.x = 2
    }else{
        high.mode.x = 1
    }
    
    n.high.x = sum(x.gmm$classification==high.mode.x)
    n.low.x = sum(x.gmm$classification!=high.mode.x)
    
    n.high.y = sum(y.gmm$classification==high.mode.y)
    n.low.y = sum(y.gmm$classification!=high.mode.y)
    
    contig.table = data.frame(c(n.high.x, n.low.x), c(n.high.y, n.low.y))
    
    line.col = "grey50"
    line.thick = 1
    lab.loc.x = 1
    lab.loc.y = 0.2
    cex.all = 12
    x.mean = mean(x)
    y.mean = mean(y)
    tt = t.test(x,y)
    rst = wilcox.test(x,y)
    #print(head(d))
    
    
    x.high.mean = mean(x[x.gmm$classification==high.mode.x])
    x.low.mean = mean(x[x.gmm$classification!=high.mode.x])
    y.high.mean = mean(y[y.gmm$classification==high.mode.y])
    y.low.mean = mean(y[y.gmm$classification!=high.mode.y])

    info("------------")
    info(sprintf("X-high: %s", signif(x.high.mean, 3)))
    info(sprintf("X-low: %s", signif(x.low.mean, 3)))
    info(sprintf("Y-high: %s", signif(y.high.mean, 3)))
    info(sprintf("Y-low: %s", signif(y.low.mean, 3)))
    info("------------")
    
    length(x) = length(y)
    d = melt(cbind.data.frame(x,y))
    g.hist = ggplot(d, aes(x=value,fill=variable, color=variable)) + theme_bw() + geom_density(alpha=0.2) + # ggtitle(sprintf("%s (Density hist)", plot.title)) + 
      theme(axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"),
            axis.line.y = element_line(color="gray15", size = 0.75, lineend="round")) + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border= element_blank(),
            panel.background =  element_blank(),
            legend.key = element_blank(),
            text = element_text(size=cex.all, colour="gray22")) +
      geom_vline(xintercept=x.high.mean, lty=2, colour=brewer.pal(9, "Paired")[5]) +
      geom_vline(xintercept=x.low.mean, lty=2, colour=brewer.pal(9, "Paired")[6]) +
      
      geom_vline(xintercept=y.high.mean, lty=2, colour=brewer.pal(9, "Paired")[3]) +
      geom_vline(xintercept=y.low.mean, lty=2, colour=brewer.pal(9, "Paired")[4]) +
      
      scale_fill_manual("", values=brewer.pal(9, "Set1")) + 
      scale_color_manual("", values=brewer.pal(9, "Set1"), guide=F) + ggtitle("Density hist for bimodal test")
    #xlab(xlab) + ylab(ylab)
      #geom_segment(aes(x = x.mean, y = lab.loc.y, xend = y.mean, yend = lab.loc.y), color = line.col, size = line.thick, arrow = arrow(ends="both", length = unit(0.01, "npc"))) +
      #annotate("text", x = lab.loc.x, y = lab.loc.y, label = sprintf("t-test p = %s\n rank-sum p = %s", signif(tt$p.value, 3), signif(rst$p.value, 3)), hjust=0) + 
      
    p.bimod.fisher = fisher.test(contig.table)$p.value
    #print(contig.table)
    return(list("p.bimod.fisher"=p.bimod.fisher, "p.wilcox"=rst$p.value, 
    	"p.ttest"=tt$p.value, "p.dip"=max(c(dt.x$p.value, dt.y$p.value)), 
    	"plot"=g.hist))
}


do.ks.test <- function(d, do.plot=T, cex.all=12, line.col="grey10", line.thick=0.5, xlab="X", 
	ylab="Cumulative Probability", pdf=NULL, w=8, h=7, plot.title=NULL, run.bimod=F)
{
	if(ncol(d)!=2)
	{
		stop("Input data frame must have two columns")
	}

	if(is.null(plot.title)){
		plot.title = paste(colnames(d), collapse=" vs ")
	}

	x.test = unlist(d[,1]) [!is.na(unlist(d[,1]))]
	y.test = unlist(d[,2]) [!is.na(unlist(d[,2]))]

	x.mean = mean(x.test)
	y.mean = mean(y.test)
	direction = sign(x.mean - y.mean)

	e.x <- ecdf(x.test)
	e.y <- ecdf(y.test)
	
	u = compare(x.test, y.test)
	max.d = abs(e.x(u) - e.y(u))

	d = melt(d)

	r = ks.test(x.test,y.test)
	groups = c(rep("X", length(x.test)), rep("Y", length(y.test)))
	tt = t.test(c(x.test,y.test)~ groups)
	rst = wilcox.test(c(x.test,y.test)~ groups)


	if(run.bimod)
	{
		bimod = do.bimodal.test(x.test, y.test)
	}else{
		bimod = NA
	}

	if(do.plot)
	{
		# find and show point of maximal difference
		#
		x.lab.loc = 0.9*max(d$value, na.rm=T)

		info(sprintf("Maximum difference at %s [D=%s]. Xlab at %s", u, max.d, x.lab.loc))
		
		g = ggplot(d, aes(x=value,color=variable)) + stat_ecdf() + ggtitle(sprintf("%s (CDF)", plot.title)) + 
			theme(axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"),
       				 axis.line.y = element_line(color="gray15", size = 0.75, lineend="round")) + 
				theme(panel.grid.major = element_blank(),
			    	panel.grid.minor = element_blank(),
					panel.border= element_blank(),
			 		panel.background =  element_blank(),
			    	legend.key = element_blank(),
			    	text = element_text(size=cex.all, colour="gray22")) +
			geom_vline(xintercept=u, lty=2)  + 
			geom_segment(aes(x = u, y = e.y(u), xend = u, yend = e.x(u)), color = line.col, size = line.thick, arrow = arrow(ends="both", length = unit(0.01, "npc"))) +
			annotate("text", x = x.lab.loc, y = 0.1, label = sprintf("D = %s\n p = %s", signif(r$statistic, 3), signif(r$p.value, 3))) + 
			scale_color_manual("", values=brewer.pal(9, "Set1")) + 
			xlab(xlab) + ylab(ylab)

		lab.loc.x = 0.2 + max(c(x.mean, y.mean), na.rm=T)

		max.dens.x = max(density(x.test)$x)
		max.dens.y = max(density(y.test)$x)
		lab.loc.y = 0.1 * max(c(max.dens.x, max.dens.y))

		info(sprintf("max densities. x=%s, y=%s", max.dens.x, max.dens.y))

		g.hist = ggplot(d, aes(x=value,fill=variable, color=variable)) + theme_bw() + geom_density(alpha=0.2) + ggtitle(sprintf("%s (Density hist)", plot.title)) + 
			theme(axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"),
       				 axis.line.y = element_line(color="gray15", size = 0.75, lineend="round")) + 
				theme(panel.grid.major = element_blank(),
			    	panel.grid.minor = element_blank(),
					panel.border= element_blank(),
			 		panel.background =  element_blank(),
			    	legend.key = element_blank(),
			    	text = element_text(size=cex.all, colour="gray22")) +
			geom_vline(xintercept=x.mean, lty=2, colour=brewer.pal(9, "Set1")[1]) +
			geom_vline(xintercept=y.mean, lty=2, colour=brewer.pal(9, "Set1")[2]) +
			geom_segment(aes(x = x.mean, y = lab.loc.y, xend = y.mean, yend = lab.loc.y), color = line.col, size = line.thick, arrow = arrow(ends="both", length = unit(0.01, "npc"))) +
			annotate("text", x = lab.loc.x, y = lab.loc.y, label = sprintf("t-test p = %s\n rank-sum p = %s", signif(tt$p.value, 3), signif(rst$p.value, 3)), hjust=0) + 
			scale_fill_manual("", values=brewer.pal(9, "Set1")) + 
			scale_color_manual("", values=brewer.pal(9, "Set1"), guide=F) + 
			xlab(xlab) + ylab(ylab)
			
		if(is.null(pdf)){
			print(g)
		}else{
			info(sprintf("Saving plot to %s", pdf))
			# ggsave(g, filename=pdf, width=w, height=h)
			# ggsave(g.hist, filename=pdf, width=w, height=h)
			if(run.bimod){
				l = list(g, g.hist, bimod$plot)
			}else{
				l = list(g, g.hist)
			}
			
			pdf(pdf, width=w, height=h)
			invisible(lapply(l, print))
			dev.off()
		}
	}


	rval = list(direction=direction,
				p.ks=r$p.value, 
				d.statistic=r$statistic,
				p.ttest=tt$p.value, 
				t.statistic=tt$statistic, 
				p.wilcox=rst$p.value, 
				w.statistic=rst$statistic)
	if(run.bimod){
		if(!is.null(bimod)){
			rval$p.bimod=bimod$p.bimod.fisher
			rval$p.dip = bimod$p.dip
		}else{
			rval$p.bimod=NA
			rval$p.dip = NA
		}
	}else{
		rval$p.bimod=NA
		rval$p.dip = NA
	}
	return(rval)
}

# plot
plot.DE.volcano <- function(DE.table,title, top.genes.up=20, top.genes.down=0){
  library(ggrepel)
  df <- DE.table
  
  if(top.genes.down==0){top.genes.down=top.genes.up}

  df.up <- df[df$log2fc>0,]
  pval.cutoff.up <- sort(df.up$p,decreasing = F)[top.genes.up]
  label.genes.up <- subset(df.up,subset= p < pval.cutoff.up)

  df.down <- df[df$log2fc<0,]
  pval.cutoff.down <- sort(df.down$p,decreasing = F)[top.genes.down]
  label.genes.down <- subset(df.down,subset= p < pval.cutoff.down)
  


  ggplot(data=df,aes(log2fc,-log10(p))) + 
    geom_point(size=0.5) + 
    geom_point(data=label.genes.up,aes(log2fc,-log10(p)),size=1,colour='red') +
    geom_point(data=label.genes.down,aes(log2fc,-log10(p)),size=1,colour='green') +

    geom_text_repel(data=label.genes.up,aes(log2fc,-log10(p),label=rownames(label.genes.up))) + 
    geom_text_repel(data=label.genes.down,aes(log2fc,-log10(p),label=rownames(label.genes.down))) + 

    labs(title=title,x='log2 fold change',y='-log10(pvalue)') +
    theme(panel.grid.major = element_blank(),
  		    panel.grid.minor = element_blank(),
  				panel.border= element_blank(),
  		 		panel.background =  element_rect(fill = 'white'),
  	    	strip.background = element_blank(), 
  	    	strip.text.x = element_text(size = 12, face="bold"),
  	    	axis.title.x = element_text(vjust=-0.5),
  	    	legend.key = element_blank(),
  	    	text = element_text(size=12, colour="gray22"),
  				axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"),
  			  axis.line.y = element_line(color="gray15", size = 0.75, lineend="round")) 
}


bayes.mean.diff <- function(x, groups, R = 100, do.plot=F)
{
	library(bayesboot)
	groups= factor(groups)
	g1 = levels(groups)[1]
	g2 = levels(groups)[2]

	b1 <- bayesboot(x[groups==g1], weighted.mean, use.weights = TRUE, R=R)
	b2  <- bayesboot(x[groups==g2], weighted.mean, use.weights = TRUE, R=R)

	# Calculating the posterior difference
	# (here converting back to a bayesboot object for pretty plotting)
	b_diff <- as.bayesboot(b1 - b2)
	if(do.plot){
		plot(b_diff)
	}
	s=summary(b_diff)
	## mean, hdilow, hdihigh (95% highest density interval)
	m = s[[3]][1]
	l = s[[3]][3]
	u = s[[3]][4]

	##quantiles
	q2.5 = s[[3]][5]
	q25 = s[[3]][6]
	median = s[[3]][7]
	q75 = s[[3]][8]
	q97.5 = s[[3]][9]

	return(list("difference.posterior"=b_diff, "estimated_diff"=m, "lower.bound"=l, "upper.bound"=u, "q2.5"=q2.5, "q25"=q25, "median"=median, "q75"=q75,"q97.5"=q97.5))
}

do.anova.test <- function(data, groups, gene.index, gene=NULL, seurat.obj=NULL, groups.second=NULL, ret.fit=T)
{
	if(is.null(gene))
	{
		gene = rownames(data)[gene.index]
	}
	if(is.null(seurat.obj)){
		x = unlist(data[gene,])
	}else{
		x = unlist(fetch.data(seurat.obj, gene))
	}
	
	if(!is.null(groups.second))
	{
		info("Running Two-way ANOVA")
		fit <- aov(x ~ groups * groups.second)
		pval = summary(fit)[[1]][["Pr(>F)"]][1]
	}else{
		fit <- aov(x ~ groups)
		pval = summary(fit)[[1]][["Pr(>F)"]][1]
	}
	
	
	if(ret.fit){
		return(list(gene=gene, p=pval, fit=fit))
	}else{
		return(list(gene=gene, p=pval))
	}
	
}

find.de.genes.anova <- function(data, groups, verbose=T, n.cores=1, debug.do.only.first=0)
{
	
	genes = vector()
	pvals = vector()
	d = data
	if(debug.do.only.first>0){
		d = d[1:debug.do.only.first,]
	}
	
	if(n.cores>1){
		library(parallel)
		library(foreach)
		library(doParallel)
		cl<-makePSOCKcluster(n.cores, outfile="")
		registerDoParallel(cl, n.cores)
		rval=tryCatch({		
			chunksize = nrow(d)/n.cores
			vals = split(1:nrow(d), ceiling(seq_along(1:nrow(d))/chunksize))
			p <- foreach(run.id=1:n.cores,.export=c("bayes.mean.diff", "alpha", "get.non.zero", 
				"info", "compare.groups", "do.anova.test", "flog.info"), .combine=c) %dopar% {
				# home = path.expand("~")
				# source(paste(home, "dev/adam/rna_seq/r/permutation.r", sep="/"))
				# source(paste(home, "dev/adam/rna_seq/r/util.r", sep="/"))
				info(sprintf("Worker thread #%s starting", run.id))
				v = vals[[run.id]]
				init.val = v[1]
				lapply(v, function(i) {
					x = unlist(d[i,])
					gene.num = i - init.val
					if(gene.num %% 200==0 && verbose && run.id==1)
					{
						info(sprintf("%s%% complete. Tested %s [%s of %s]", 
							100*round(gene.num/length(v), 3), rownames(d)[i], i, length(v)))
					}
					do.anova.test(d, groups, i, ret.fit=F)
				})
			}
			rv = data.frame(matrix(unlist(p), nrow=nrow(d), byrow=T),stringsAsFactors=FALSE)
			colnames(rv) = names(p[[1]])	
			list("results"=rv, "cluster"=cl)
		}
		, interrupt = function(ex) {
			cat("An interrupt was detected.\n");
			print(ex);
			list("cluster"=cl)
		}, error = function(ex) {
			cat("An error was detected.\n");
			print(ex);
			list("cluster"=cl)
		}, finally = {
			list("cluster"=cl)
		})
		cat("Unregistering parallel backend.. ");
		stopCluster(rval$cluster)
		registerDoSEQ()
		cat("Done.\n")
		return(rval)

	}else{
		i = 1
		for(gene in rownames(d))
		{
			y = do.anova.test(d, groups, gene=gene)
			genes[i] = y$gene
			pvals[i] = y$p

			i = i + 1
			if(i %% 200==0 && verbose)
			{
				info(sprintf("%s%% complete. Tested %s [%s of %s] --  pval=%s", 
					100*round(i/nrow(d), 3), rownames(d)[i], i, nrow(d), y$p))
			}
		}
		return(data.frame(genes, pvals))
	}
}


# calculate summary statistics and run the relevant statistical test.
do.test <- function(x, groups, i, d, perm=F, wilcox=T, verbose=T, np=100, estimate.diff=F, estimate.diff.boostraps=100)
{
	gp1.lab = levels(groups)[1] 
	gp2.lab = levels(groups)[2]
	y = list()
	a1 = alpha(x[groups==gp1.lab])
	a2 = alpha(x[groups==gp2.lab])
	y[[sprintf("%s.alpha",gp1.lab)]] <- a1
	y[[sprintf("%s.alpha",gp2.lab)]] <- a2
	y$alpha.diff = a1 - a2
	y[[sprintf("%s.mean", gp1.lab)]] <- mean(x[groups==gp1.lab])
	y[[sprintf("%s.mean", gp2.lab)]] <- mean(x[groups==gp2.lab])
	y[[sprintf("%s.mu", gp1.lab)]] <- mean(get.non.zero(x, groups, gp1.lab))
	y[[sprintf("%s.mu", gp2.lab)]] <- mean(get.non.zero(x, groups, gp2.lab))
	y$log2fc = mean(x[groups==gp1.lab]) - mean(x[groups==gp2.lab])
	if(estimate.diff)
	{
		b = bayes.mean.diff(x, groups, R=estimate.diff.boostraps)
		y$Log2FC.mean.bb = b$estimated_diff
		y$Lower.bound = b$lower.bound
		y$Upper.bound = b$upper.bound
		y$q2.5 = b$q2.5
		y$q25 = b$q25 
		y$q75 = b$q75
		y$q97.5= b$q97.5
	}
	if(!perm)
	{
		if(wilcox){
			y$p = wilcox.test(x~groups)$p.value
			if(is.na(y$p)){y$p=1}
		}else{
			y$p = t.test(x~groups)$p.value
		}
	}else
	{
		z=perm.test(x, groups,n=np)
		y$p = z$p.fc
		y$p.alpha = z$p.alpha
	}
	return(y)
}


# NOTE doMC does not work with RStudio.
# for a two factor group variable, test each gene for DE
compare.groups <- function(d, groups, 
	compare=NULL, 
	return_mast_glm=FALSE,
	np=50, 
	n.cores=1,  
	subsample=0, ## take only a few samples from each
	subsample_batch_equal = F, # if true, subsample equally from the batches
	batch=NULL, # correct for other factors. should eventually be a dataframe (to cope with multiple), but for now a factor (for a single covariate).
	first.n.genes=0, verbose=T,
	use.genes=NULL, 
	test="wilcox", # can also be perm, or MAST.
	mast.component="H", #D (discrete), C (continuous) or H (hurdle) for combined pvalues from MAST
	estimate.diff=F, #use the bayesian bootstrap to get bounds on difference between groups
	estimate.diff.boostraps=100,
	filter.min.Q=0.05, 
	filter.min.alpha=0.25, 
	filter.min.mu=1)
{

	use.mast = F
	use.wilcox = F
	use.perm = F
	switch(test,
		wilcox={
			info("Using Wilcoxon Rank-sum")
			if(!is.null(batch)){stop("Can only do batch-correction with test=mast")}
			use.wilcox = T
			},
		permutation={
			if(!is.null(batch)){stop("Can only do batch-correction with test=mast")}
			info(sprintf("Using permutation test [n=%s]", np))
			use.perm = T
		},
		mast={
			info("Using MAST test (hurdle model)")
			use.mast = T
			info("Done")
		},
		{
			if(!is.null(batch)){stop("Can only do batch-correction with test=mast")}
			info("Using two-sample t-test")
		}
	)
	if(!is.null(use.genes))
	{	
		matches = match(toupper(use.genes), toupper(rownames(d)))
		d = d[matches[!is.na(matches)], ]
		info(sprintf("Shrunk counts table to %i genes ", nrow(d)))
		nas = sum(is.na(d))
		if(nas > 0)
		{	
			warn(sprintf("%s NAs in table!", nas))
		}
	}
	ptm <- proc.time()
	
	if(is.null(compare))
	{
		if(length(unique(groups))!=2)
		{
			stop("If compare is NULL, groups vector must have exactly two labels!")
		}
		compare = unique(groups)
	}else
	{
		if(length(compare)==1)
		{
			info(sprintf("Comparing %s against a background of all other groups", compare[1]))
			g = as.character(groups)
			g[groups==compare[1]] = compare[1]
			g[groups!=compare[1]] = "Background"
			groups = g
			compare = c(compare[1], "Background")
		}else
		{
			if(!is.null(batch)){batch = batch[groups %in% compare]}
			d = d[, groups %in% compare]
			groups = groups[groups %in% compare]
			info(sprintf("Comparing only the %s samples in %s", ncol(d), paste(compare, collapse=", ")))
		}
	}
	groups = factor(groups, levels=compare)
	g1 = which(groups==unique(groups)[1])
	g2 = which(groups==unique(groups)[2])
	if(length(g1) == 0 | length(g2)==0){stop("One or more groups is empty!")}
	if(subsample>0)
	{
		info(sprintf("Subsampling to %s", subsample))
		if(subsample_batch_equal){
			if(is.null(batch)){stop("Must provide batch labels if you want to subsample batches equally")}
			warn("Equal batch sampling isn't finished, its only approximate and probably buggy")
			b1 = batch[g1]
			b2 = batch[g2]
			if(length(g1) > subsample) b1i = resample(b1, max.group.size=floor(subsample / length(unique(b1))), verbose=F, return.index=T)else{b1i=1:length(g1)}
			if(length(g2) > subsample) b2i = resample(b2, max.group.size=floor(subsample / length(unique(b2))), verbose=F, return.index=T)else{b2i=1:length(g2)}
			keep = c(g1[b1i], g2[b2i])
		}else{
			if(length(g1) > subsample) g1 = sample(g1, subsample)
			if(length(g2) > subsample) g2 = sample(g2, subsample)
			keep = c(g1,g2)
		}
		d = d[, keep]
		groups = groups[keep]
		if(!is.null(batch)) batch = batch[keep]
		info(sprintf("Table is now %s", paste(dim(d), collapse=" x ")))
		
	}
	print(table(groups))
	fc = vector()
	p = vector()
	if(first.n.genes>0)
	{
		d = d[1:first.n.genes,]
	}
	if(use.mast){
		###MAST single-cell
		info("Loading MAST")
		library(MAST)
		### multi-core
		options(mc.cores=n.cores)
		### build 'colData' to store covariates
 		numGen<-dim(d)[1] 
		fdat<-matrix(rep(1,numGen))
		cngeneson = scale(Matrix::colSums(d > 0)) ### normalize for difference in cell complexity
		colData = cbind.data.frame(groups, cngeneson)
		if(!is.null(batch))
		{
			colData$batch = batch
		}
		rownames(colData) = colnames(d)

		#print(head(colData))
		### set up MAST single-cell assay object:
		sca = MAST::FromMatrix(as.matrix(d), colData, data.frame(fdat))
		info(sprintf("Fitting MAST 'hurdle' model [n.cores=%s]", n.cores))
		if(!is.null(batch)){
			info(sprintf("Batch covariate has %s levels:", length(unique(colData$batch))))
			print(table(colData$batch))
			info("Distributed over the groups:")
			print(table(groups, colData$batch))
			zlm.output <- zlm(~groups + cngeneson + batch, sca)
		}else{
			zlm.output <- zlm(~groups + cngeneson, sca)
		}
		info("Running MAST likelihood-ratio test (LRT)")
		grp_string = paste0('groups', compare[2])
		summaryCond <- summary(zlm.output, doLRT=grp_string)
		res <- summaryCond$datatable
		res <- data.frame(res)
		rv <- makeNice(res, val=grp_string, comp=mast.component)
		#return(list("rv"=rv, "res"=res))
		### keep MASTs FC estimate, and add the fold change of means
		rv$log2fc_MAST = -1 * rv$log2fc
		mg1 = Matrix::rowMeans(d[rownames(rv), groups==levels(groups)[1]])
		mg2 = Matrix::rowMeans(d[rownames(rv), groups==levels(groups)[2]])
		rv$log2fc = mg1 - mg2	
		### remove a couple superfluous columns
		rv$Gene <- NULL
		rv$padj <- NULL
		rv = rv[rownames(d), ]
		options(mc.cores=1)
		#return(res)
	}else{
		zlm.output = new(Class = "ZlmFit") #empty
		export.obj.names = c("bayes.mean.diff", "alpha", "get.non.zero", "info", "compare.groups", "do.test", "flog.info")
		
		### Here wrote code for pblapply, but its really slow. back to the original [May, 2018]
		# if(n.cores==1)
		# {
		
		# cl <- makeCluster(as.integer(n.cores))
		# clusterExport(cl, export.obj.names)
		# info(sprintf("Running Wilcoxon test for each gene [using %s cores]", n.cores))
		# p = pblapply(1:nrow(d), function(i) {
		# 	x = unlist(d[i,])
		# 	y= do.test(x, groups, i, d, perm=use.perm, wilcox=use.wilcox, verbose=verbose, np=np, estimate.diff=estimate.diff, estimate.diff.boostraps=estimate.diff.boostraps)
		# }, cl=cl)
		# info("Done!")

		# rv = data.frame(matrix(unlist(p), nrow=nrow(d), byrow=T),stringsAsFactors=FALSE)
		# colnames(rv) = names(p[[1]])
		# info("Cleaning up cluster resources")
		# cl <- NULL
		# gc()
		# info("Done!")

		## the below code works, but i switched to pblapply for simplicity. [March, 2018]
		#}else{
		
			library(parallel)
			library(foreach)
			library(doParallel)
			cl<-makePSOCKcluster(n.cores, outfile="")
			registerDoParallel(cl, n.cores)
			rval=tryCatch({		
					chunksize = nrow(d)/n.cores ### break the list of genes into chunks, one for each thread.
					vals = split(1:nrow(d), ceiling(seq_along(1:nrow(d))/chunksize))
					p <- foreach(run.id=1:n.cores,.export=export.obj.names, .combine=c) %dopar% {
						library(Matrix)
						v = vals[[run.id]]
						init.val = v[1]
						lapply(v, function(i) {
							gene.num = i - init.val
							x = unlist(d[i,])
							if(gene.num %% 200==0 && verbose && run.id==1)
							{
								info(sprintf("%s%% complete. Tested %s [%s of %s]", 
									100*round(gene.num/length(v), 3), rownames(d)[i], i, length(v)))
							}
							do.test(x, groups, i, d, perm=use.perm, wilcox=use.wilcox, verbose=run.id==1, np=np, estimate.diff=estimate.diff, estimate.diff.boostraps=estimate.diff.boostraps)
						})
					}
					rv = data.frame(matrix(unlist(p), nrow=nrow(d), byrow=T),stringsAsFactors=FALSE)
					colnames(rv) = names(p[[1]])	
					list("results"=rv, "cluster"=cl)
			}
			, interrupt = function(ex) {
				cat("An interrupt was detected.\n");
				print(ex);
				list("cluster"=cl)
			}, error = function(ex) {
				cat("An error was detected.\n");
				print(ex);
				list("cluster"=cl)
			}, finally = {
				list("cluster"=cl)
			})
			rv = rval$results
			cat("Unregistering parallel backend.. ");
			stopCluster(rval$cluster)
			registerDoSEQ()
			cat("Done.\n")
	}

	info("Processing DE results")
	info(sprintf("FDR correcting [%s p-values]", length(rv$p)))
	rv$p.adj = p.adjust(rv$p, method = "fdr")
	if("p.alpha" %in% colnames(rv)){
		rv$p.alpha.adj = p.adjust(rv$p.alpha, method = "fdr")
	}
	rv$GENE_SYMBOL = rownames(d)
	rv=rv[, c(ncol(rv), 1:ncol(rv)-1)]

	info("Filtering results")
	initial = nrow(rv)
	if("p.alpha" %in% colnames(rv)){
		rvf = rv[rv$p.alpha.adj<filter.min.Q | rv$p.adj<filter.min.Q, ]
	}else{
		rvf = rv[rv$p.adj<filter.min.Q, ]
	}
	
	info(sprintf("Dropped %s non-significant genes [Q<%s]", initial-nrow(rvf), filter.min.Q))
	initial = nrow(rvf)
	if(length(grep("mu", colnames(rvf))) > 0)
	{
		up.and.high = which(rvf$log2fc > 0 & (rvf[, sprintf("%s.mu", compare[1])] > filter.min.mu))
		down.and.high = which(rvf$log2fc < 0 & (rvf[, sprintf("%s.mu",compare[2])] > filter.min.mu))
		rvf = rvf[c(up.and.high, down.and.high), ]
		info(sprintf("Dropped %s lowly expressed genes [mu<%s]", initial-nrow(rvf), filter.min.mu))

		initial = nrow(rvf)
		up.and.frequent = which(rvf$log2fc > 0 & (rvf[, sprintf("%s.alpha", compare[1])] > filter.min.alpha))
		down.and.frequent = which(rvf$log2fc < 0 & (rvf[, sprintf("%s.alpha",compare[2])] > filter.min.alpha))
		info(sprintf("Dropped %s rarely detected genes [alpha<%s]", initial-nrow(rvf), filter.min.alpha))
		rvf = rvf[c(up.and.frequent, down.and.frequent), ]
	}
	
	info(sprintf("Sorting %s filtered genes", nrow(rvf)))
	rvf = rvf[order(abs(rvf$log2fc), decreasing=T),]
	
	su = rvf[rvf$log2fc>0,]
	sd = rvf[rvf$log2fc<0,]
	
	rownames(rv) = rv$GENE_SYMBOL

	tf_file = "~/ref/gene_lists/gut_circuits/List_of_TFs_MusMusculus_AnimalTFDB.txt"
	if(file.exists(tf_file)){
	tfs = read.delim(tf_file)
	rv_tfs = rvf[ rvf$GENE_SYMBOL %in% tfs$Gene_symbol,]
	}else{rv_tfs=data.frame()}
	info(sprintf("%s genes up, %s genes down. %s TFs DE.", nrow(su), nrow(sd), nrow(rv_tfs)))

	print("Runtime:")
	print(proc.time() - ptm)

	compare_name = paste(compare[1], compare[2], sep="_vs_")
	if(!return_mast_glm){zlm.output=new("ZlmFit")}
	if(is.null(batch)){batch=factor("NULL")}
	rval <- new("DE.Results", 
					mast.glm=zlm.output,
					name=as.character(compare_name), 
					all=rv, 
					sig=su, 
					sig.down=sd, 
					sig.tfs=rv_tfs,
					min.p=filter.min.Q, 
 					groups=as.character(groups),
 					batch=as.character(batch))
	info("Cleaning up resources")
	sca <- NULL
	d <- NULL
	rvf <- NULL
	gc()
	return(rval)
}


# return the fraction of non-zero values
alpha <- function(x, na.rm=T)
{
	sum(x>0)/length(x)
}

# wrapper to avoid returning NAs
get.non.zero <- function(x, groups, label)
{
	cond = groups==label & x>0
	if(any(cond)){
		return(x[cond])
	}else{
		return(0)
	}
}

# permutation test for differences in the median values.
# if test.non.zero.values is true, the test is performed for both:
#   (a) the fraction, alpha, of non-zero values
#	(b) the median of the non-zero values
perm.test <- function(x, groups, n=100, statistic=mean, test.non.zero.values=F, alpha.thresh=0)
{
	y = length(x)
	#cat(sprintf("Running permutation test [n=%s]\n", n))
	if(length(unique(groups))!=2)
	{
		cat("ERROR: Groups vector must have exactly two labels!\n")
		return(FALSE)
	}

	gp1.lab = levels(groups)[1] 
	gp2.lab = levels(groups)[2]
	
	mu = c(statistic(get.non.zero(x, groups, gp1.lab)), statistic(get.non.zero(x, groups, gp2.lab)))
	mean.expr = c(mean(x[groups==gp1.lab]), mean(x[groups==gp2.lab]))

	if(test.non.zero.values)
	{
		fc = statistic(get.non.zero(x, groups, gp1.lab)) - statistic(get.non.zero(x, groups, gp2.lab)) # in log space, ratio is a difference
	}else{
		fc = statistic(x[groups==gp1.lab]) - statistic(x[groups==gp2.lab])
	}

	a1 = alpha(x[groups==gp1.lab])
	a2 = alpha(x[groups==gp2.lab])
	ad = a1 - a2
	
	if(fc<1){flip.fc=T}else{flip.fc=F}

	if(ad<0){flip.alpha=T}else{flip.alpha=F}
	
	success=0
	success.alpha=0
	for(i in 1:n)
	{
		shuf.gps = sample(groups, y) #permuted group labels
		
		if(test.non.zero.values)
		{
			fc.test =  statistic(get.non.zero(x, shuf.gps, gp1.lab)) - statistic(get.non.zero(x, shuf.gps, gp2.lab))
		}else{
			fc.test =  statistic(x[shuf.gps==gp1.lab]) - statistic(x[shuf.gps==gp2.lab])
		}

		# print("group1 nonzero vals:")
		# print(get.non.zero(x, shuf.gps, gp1.lab))
		# print("Group2 nonzero vals:")
		# print(get.non.zero(x, shuf.gps, gp2.lab))
		# print(median(get.non.zero(x, shuf.gps, gp1.lab)))
		# print(fc.test)

		a.test = alpha(x[shuf.gps==gp1.lab]) - alpha(x[shuf.gps==gp2.lab])
		
		if(flip.fc){
			if(fc.test <= fc){success=success+1}
		}else{
			if(fc.test >= fc){success=success+1}
		}

		if(flip.alpha)
		{
			if(a.test <= ad){success.alpha=success.alpha+1}
		}else{
			if(a.test >= ad){success.alpha=success.alpha+1}
		}
		
		#cat(sprintf("Permutation %s: Stat=%s, Successes=%s\n", i, fc.test, success))
	}
	return(list("p.fc"=success/n, "p.alpha"=success.alpha/n, "log2fc"==fc))
}

#plot distance between means against and n samples
test.permutation.test <- function()
{
	d = seq(0,5, by=0.2)
	n = 10:20
	vp = matrix(nrow=length(d), ncol=length(n))
	colnames(vp) = n
	rownames(vp) = d
	vtt = matrix(nrow=length(d), ncol=length(n))

	colnames(vtt) = n
	rownames(vtt) = d

	for(j in 1:length(n))
	{
		npg = n[j]	 #num per group
		g = rep(c("A", "B"), each=npg)
		cat(sprintf("%s per group. %s%% done\n", npg, 100*(j/length(n))))
		for(i in 1:length(d))
		{
			
			x1 = mvrnorm(n=npg, mu = 100, Sigma = 2)
			x2 = mvrnorm(n=npg, mu = 100+d[i], Sigma = 2)
			x = c(x1,x2)
			#cat(sprintf("d=%s\n", d[i]))
			vp[i,j] = perm.test(x, groups=g, n=10000)$p.value
			vtt[i,j] = t.test(x~g)$p.value
			vmw[i,j] = wilcox.test(x~g)$p.value
		}
	}
	pdf('test.pdf', width=16, height=13)
	par(mfrow=c(1,2))
	d = -log10(vp)
	max.disp=4
	d[d>max.disp] <- max.disp
	aheatmap(d, main = "Permutation test", Rowv=NA, Colv=NA)

	d = -log10(vtt)
	d[d>max.disp] <- max.disp
	aheatmap(d, main = "t-Test", Rowv=NA, Colv=NA)

	dev.off()
	list("perm"=vp, "t.test"=vtt)
}




### Different way of finding markers. Fit an RF to each gene and rank by OOB prediction error when predicting the cell type of interest.
### Empirically seems to perform better than ranking by hold change. Adam, May 29, 2021
### Note can't get 'future' to run fast - use only pblapply
### @param obj is a seurat object
### @param group is the factor level of interest
### @param group_by is a clustering name, or other factor in the metadata table
get_markers_rf = function(obj, group, group_by, min_frac_expr=0.05, ncores=parallel::detectCores(logical=F), max.group.size=0) #, parallel_engine = "pblapply")
{
	#if(!parallel_engine %in% c("future", "pblapply")){stop("parallel_engine must be one of 'future' or 'pblapply'")}
	# options(future.globals.maxSize= 2000*1024^2)
	tictoc::tic()
	umis = GetAssayData(obj, slot = "counts")
	index = obj@meta.data[, group_by] == group
	y = umis[, index]
	x = Matrix::rowSums(y > 0) / ncol(y)
	x = x[x > min_frac_expr]
	genes = names(x)
	if(max.group.size > 0){
		use_cells = resample(obj@meta.data[, group_by], max.group.size=max.group.size, verbose=F, return.index=T)
		obj = subset(obj, cells=colnames(obj)[use_cells])
	}
	obj@meta.data[, "is_cell_of_interest"] = obj@meta.data[, group_by] == group
	cat(sprintf("Fitting random forest to %s genes in %s cells [using %s cores]\n", length(genes), ncol(obj), ncores))
	#if(parallel_engine=="pblapply"){
		#future::plan(future::sequential)
		#cl = parallel::makeCluster(ncores, setup_timeout = 0.5)
		#parallel::clusterExport(cl=cl, list("obj"), envir=environment())
		rv = pbapply::pblapply(genes, function(g) {
			d=Seurat::FetchData(obj, c(g, "is_cell_of_interest")); 
			colnames(d)[1]="gene"; 
			ranger::ranger(formula=is_cell_of_interest ~ gene, data=d)$prediction.error}, cl = ncores)
	# Could compute other prediction metrics like F1 score
	# pred = predict(rf, FetchData(nasal, c("Ltc4s", "is_tuft_gland")))
	# cm = caret::confusionMatrix(factor(pred$predictions), factor(as.numeric(nasal$is_tuft_gland)))
	
	# Future Parallel processing didn't seem to work well:	
	# }else{
	# 	future::plan(future::multisession, workers = ncores)
	# 	rv = furrr::future_pmap(list(genes), function(g) {
	# 		d=FetchData(obj, c(g, "is_cell_of_interest")); 
	# 		colnames(d)[1]="gene"; 
	# 		ranger::ranger(formula=is_cell_of_interest ~ gene, data=d)$prediction.error}, .progress=T)
	# }
	rv = data.frame("rf_oob_err"=unlist(rv))
	rv$gene = genes
	rv$frac_expr = unname(unlist(x))
	rv = rv %>% arrange(rf_oob_err)
	tictoc::toc()
	rv
}

### Different way of finding markers. Fit an RF to each gene and rank by OOB prediction error when predicting the cell type of interest.
### Empirically seems to perform better than ranking by hold change. Adam, May 29, 2021
### Note can't get 'future' to run fast - use only pblapply
### @param obj is a seurat object
### @param group is the factor level of interest
### @param group_by is a clustering name, or other factor in the metadata table
get_markers_xgb = function(obj, group, group_by, min_frac_expr=0.05, ncores=parallel::detectCores(logical=F), max.group.size=0) #, parallel_engine = "pblapply")
{
	#if(!parallel_engine %in% c("future", "pblapply")){stop("parallel_engine must be one of 'future' or 'pblapply'")}
	# options(future.globals.maxSize= 2000*1024^2)
	tictoc::tic()
	umis = GetAssayData(obj, slot = "counts")
	index = obj@meta.data[, group_by] == group
	y = umis[, index]
	x = Matrix::rowSums(y > 0) / ncol(y)
	x = x[x > min_frac_expr]
	genes = names(x)
	if(max.group.size > 0){
		use_cells = resample(obj@meta.data[, group_by], max.group.size=max.group.size, verbose=F, return.index=T)
		obj = subset(obj, cells=colnames(obj)[use_cells])
	}
	obj$is_cell_of_interest = obj@meta.data[, group_by] == group
	cat(sprintf("Fitting random forest to %s genes in %s cells [using %s cores]\n", length(genes), ncol(obj), ncores))
	#if(parallel_engine=="pblapply"){
		#future::plan(future::sequential)
		#cl = parallel::makeCluster(ncores, setup_timeout = 0.5)
		#parallel::clusterExport(cl=cl, list("obj"), envir=environment())
		rv = pbapply::pblapply(genes, function(g) {
			d=Seurat::FetchData(obj, c(g, "is_cell_of_interest")); 
			colnames(d)[1]="gene"; 
			ranger::ranger(formula=is_cell_of_interest ~ gene, data=d)$prediction.error}, cl = ncores)
	# Could compute other prediction metrics like F1 score
	# pred = predict(rf, FetchData(nasal, c("Ltc4s", "is_tuft_gland")))
	# cm = caret::confusionMatrix(factor(pred$predictions), factor(as.numeric(nasal$is_tuft_gland)))
	
	# Future Parallel processing didn't seem to work well:	
	# }else{
	# 	future::plan(future::multisession, workers = ncores)
	# 	rv = furrr::future_pmap(list(genes), function(g) {
	# 		d=FetchData(obj, c(g, "is_cell_of_interest")); 
	# 		colnames(d)[1]="gene"; 
	# 		ranger::ranger(formula=is_cell_of_interest ~ gene, data=d)$prediction.error}, .progress=T)
	# }
	rv = data.frame("rf_oob_err"=unlist(rv))
	rv$gene = genes
	rv$frac_expr = unname(unlist(x))
	rv = rv %>% arrange(rf_oob_err)
	tictoc::toc()
	rv
}








#### [TODO, implement find markers with XGBoost]

#### [TODO, implement multi-class classification with best discrimination of a specified class]
# https://rpubs.com/mharris/multiclass_xgboost

#xgm <- xgboost(data = as.matrix(GetAssayData(nasal)["Nrgn",]), 
                   # label = as.numeric(factor(nasal$cell_labels)), 
                  #  max.depth = 2, eta = 1, nthread = 8, nrounds = 10, objective = "multi:softmax", num_class = 1+length(unique(nasal$cell_labels)))


#numberOfClasses <- length(unique(dat$Site))
#xgb_params <- list("objective" = "multi:softprob",
#                   "eval_metric" = "mlogloss",
#                   "num_class" = numberOfClasses)
#nround    <- 50 # number of XGBoost rounds
#cv.nfold  <- 5

# Fit cv.nfold * cv.nround XGB models and save OOF predictions
#cv_model <- xgb.cv(params = xgb_params,
#                  data = train_matrix, 
#                  nrounds = nround,
#                  nfold = cv.nfold, 
#                  verbose = FALSE,
#                  prediction = TRUE)


### Different way of finding markers. Fit an RF to each gene and rank by OOB prediction error when predicting the cell type of interest.
### Empirically seems to perform better than ranking by hold change. Adam, May 29, 2021
### Note can't get 'future' to run fast - use only pblapply
### @param obj is a seurat object
### @param group is the factor level of interest
### @param group_by is a clustering name, or other factor in the metadata table
# get_markers_xgb = function(obj, group, group_by, min_frac_expr=0.05, ncores=parallel::detectCores(logical=F), max.group.size=0) #, parallel_engine = "pblapply")
# {
# 	#if(!parallel_engine %in% c("future", "pblapply")){stop("parallel_engine must be one of 'future' or 'pblapply'")}
# 	# options(future.globals.maxSize= 2000*1024^2)
# 	tictoc::tic()
# 	umis = GetAssayData(obj, slot = "counts")
# 	index = obj@meta.data[, group_by] == group
# 	y = umis[, index]
# 	x = Matrix::rowSums(y > 0) / ncol(y)
# 	x = x[x > min_frac_expr]
# 	genes = names(x)
# 	if(max.group.size > 0){
# 		use_cells = resample(obj@meta.data[, group_by], max.group.size=max.group.size, verbose=F, return.index=T)
# 		obj = subset(obj, cells=colnames(obj)[use_cells])
# 	}
# 	obj$is_cell_of_interest = obj@meta.data[, group_by] == group
# 	cat(sprintf("Fitting random forest to %s genes in %s cells [using %s cores]\n", length(genes), ncol(obj), ncores))
# 		# bstSparse <- xgboost(data = train$data, label = train$label, max.depth = 2, eta = 1, nthread = 2, nrounds = 2, objective = "binary:logistic")

# 		rv = pbapply::pblapply(genes, function(g) {
# 			d=Seurat::FetchData(obj, c(g, "is_cell_of_interest")); 
# 			colnames(d)[1]="gene"; 
# 			ranger::ranger(formula=is_cell_of_interest ~ gene, data=d)$prediction.error}, cl = ncores)

# 	rv = data.frame("rf_oob_err"=unlist(rv))
# 	rv$gene = genes
# 	rv$frac_expr = unname(unlist(x))
# 	rv = rv %>% arrange(rf_oob_err)
# 	tictoc::toc()
# 	rv
# }


# Test xgboost
# xgm <- xgboost(data = as.matrix(GetAssayData(nasal)["mach",]), 
#                    label = as.numeric(nasal$cell_labels=="Tuft-MVC"), 
#                    max.depth = 2, eta = 1, nthread = 2, nrounds = 4, objective = "binary:logistic")





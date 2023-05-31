

binned.histogram <- function(counts, num.bins=20)
{
	dc = apply(t(counts), 1, function(a){as.integer(cut2(a, g=num.bins))-1})
	dcm = melt(dc, id.vars=NULL)
}


cdf.plot <- function(counts, colour.by=NULL)
{
	# draw cumulative distribution factor to visualise sample-sample differences
	
	info("Drawing CDF")
	if(!is.null(colour.by))
	{
		counts$id = colour.by
		ctm = melt(counts, id.vars="id")
	}else
	{
		ctm = melt(counts)
	}
	
	print(colnames(ctm))
	n.col = length(unique(ctm$variable))
	ggplot(ctm,aes(log2(value+1),colour=variable)) + stat_ecdf(size=1) +  
		scale_colour_manual("Sample", values=default.cols(n.col)) + 
		ggtitle("Cumulative distribution") + theme_bw() +
		xlab("Log2 RSEM expected count") + ylab("") + 
		theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) + 
		theme(plot.title=element_text(vjust=-3)) + 
		theme(axis.line = element_line(colour = "black"),
		    panel.grid.major = element_blank(),
		    panel.grid.minor = element_blank(),
		    panel.border = element_blank(),
		    panel.background = element_blank(), 
		    text = element_text(size=20, colour="gray22"))
}
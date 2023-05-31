#functions for generating a tree-like figure to represent cell-lineage
library(Hmisc)
library(vegan)
library(igraph)
# for a given cell, find the closest cluster centres, and its location
# on the line between them.
process.cell <- function(cell.name, centres=NULL, counts=NULL, distance.method="euclidean")
{
	if(is.null(centres) | is.null(counts)){stop("Must provide cluster centres and counts data")}
	info(sprintf("Processing cell %s", cell.name))
	cell = counts[,cell.name]
	x = t(cbind(cell,centres))
	if(distance.method %in% c("pearson","spearman"))
	{
		d = 1-as.matrix(cor(t(x), method=distance.method), method=distance.method)[1, ]
	}else
	{
		d = as.matrix(dist(x, method=distance.method))[1, ]
	}
	d = sort(d)
	c1n = names(d)[2]
	c2n = names(d)[3]
	c1 = centres[,c1n]
	c2 = centres[,c2n]
	f = abs(drop(norm(cell-c1) %*% norm(c2-c1)))
	return(list("closest"=c1n, "second_closest"=c2n, "fraction"=f))
}

# pass the seurat object to generate of the spanning tree
build.map <- function(s, group.by, use.pca=T)
{
	#n.cells = 10
	info("Finding cluster centres")

	if(use.pca){
		cluster.centres = group.means(t(s@pca.rot[,1:10	]), groups=s@data.info[, group.by])
		st = make.spanning.tree(s, cluster.centres, discretise=F)
	}else{
		cluster.centres = group.means(s@data, groups=s@data.info[, group.by])
		st = make.spanning.tree(s, cluster.centres)
	}
	

	
	graph.data = to.igraph(st$tree, colnames(cluster.centres))
	
	info("Draw using ggnet2")
	library(sna)
	library(GGally)
	library(intergraph)
	gg = ggnet2(graph.data$igraph, alpha = 0.75, mode = graph.data$layout, label = T)

	info("Processing cells (allocating them to a centre)")
	cell.info = sapply(s@cell.names, process.cell, centres=cluster.centres, counts=s@data)

	l = graph.data$layout
	colnames(l) = c("x", "y")
	rownames(l) = V(graph.data$igraph)$name

	print("Graph layout:")
	print(l)

	info("Calculating coordinates")
	d = get.cell.locations(cell.info, l) #, graph.data$igraph)

	spanning.tree.plot(d, l, graph.data$edges)

	return(list("ci"=cell.info, "layout"=l, 
		"centres"=cluster.centres, "igraph"=graph.data$igraph, 
		"edges"=graph.data$edges))
}

norm <- function(v)
{
	v / sqrt(sum(v^2))
}

make.spanning.tree <- function(s, cluster.centres, num.bins=4, distance.method="euclidean", do.plot=T, discretise=T)
{
	#clusters = s@data.info$gcl.manual
	if(discretise){
		info("Discretising")
		v = t(apply(t(cluster.centres), 1, function(a){as.integer(cut2(a, g=num.bins))-1}))
		rownames(v) = colnames(cluster.centres)
	}else{v=cluster.centres}
	info("Calculating distance matrix")
	d = dist(v, method=distance.method)
	info("Finding spanning tree")
	s = spantree(d)
	if(do.plot){info("Plotting");plot(s, type='t', labels=colnames(cluster.centres))}

	return(list("tree"=s,"centres"=cluster.centres))
}

# convert a minimum spanning tree to an igraph graph.
to.igraph <- function(s_tree, node.names)
{
	info("Converting vegan MST to igraph")
	l = s_tree$kid
	#kid
	#From vegan docs: The child node of the parent, starting from parent number two. 
	#If there is no link from the parent, value will be NA and tree 
	#is disconnected at the node.
	edgelist = NULL
	for(i in 1:length(l)+1)
	{
		
		info(sprintf("Node %s is the parent of child %s", i, l[i-1])) # -1 because vegan docs says 
		edge = c(node.names[i], node.names[l[i-1]])
		info(sprintf("Edge: %s <--> %s", edge[1], edge[2]))
		
		if(is.null(edgelist))
		{
			edgelist =  as.matrix(t(edge))
		}else
		{
			edgelist = rbind(edgelist, as.matrix(t(edge)))
		}
		#print(edgelist)
	}
	g = graph_from_edgelist(edgelist)
	layout = plot.graph(g, node.size = 10, layout.type = "kamada")
	return (list("igraph"=g, "layout"=layout, "edges"=edgelist))
}

# using the info from cell info, and the locations of the
# nodes of the spanning tree from the igraph layout, work
# out where each cell should go, and return a dataframe
# 			x	y	cluster
# CELL1 	1.2	-3.2	Quiescent.Stem
# CELL2		0.5	1.2		Cycling.Stem
get.cell.locations <- function(ci, layout, do.only=0)
{
	ci = t(ci)
	if(do.only > 0)
	{
		ci = head(ci, n=do.only)
	}
	info(sprintf("Finding locations (x-y) for %s cells", nrow(ci)))
	locs = t(apply(ci, 1, loc, layout=layout))
	cbind.data.frame(ci, locs)
}

loc <- function(ci, layout)
{
	#print(ci)
	cl1 = layout[ci[["closest"]],]
	cl2 = layout[ci[["second_closest"]],]
	cl1 + ci[["fraction"]]*(cl2-cl1)
}

#return the end points of all the edges of an edgelist and a layout
get.edges <- function(e, l)
{
	info("Edgelist")
	print(e)
	print("Layout")
	print(l)
	d = NULL
	for(i in 1:nrow(e))
	{
		n1 = e[i,][1]
		n2 = e[i,][2]
		x = l[n1,][1]
		y = l[n1,][2]
		x1 = l[n2,][1]
		y1 = l[n2,][2]

		#info(sprintf("Edge #%s, %s <--> %s. %s located at (%s,%s), %s at (%s,%s)", i, n1, n2, n1, x, y, n2, x1, y1))
		
		if(is.null(d))
		{
			d = data.frame(x=x, y=y, x1=x1, y1=y1)
		}else
		{
			d = rbind(d, c(x,y,x1,y1))
		}
	}
	print(d)
	return(d)
}

spanning.tree.plot <- function(ci, layout, edgelist, plot.title="Spanning Tree", show.cells=T)
{
	 l = data.frame(layout)
	 l$cell = rownames(l)
	 ci$dens = densCols(ci$x,ci$y, colramp=colorRampPalette(rev(brewer.pal(11, "Spectral"))))
	 
	 edge.df = get.edges(edgelist, layout)
	 info("Drawing spanning tree plot")
	 print(head(ci))
	 g =ggplot(ci, aes(x=x, y=y)) + theme_bw() + scale_colour_identity() + 
	 	geom_point(data=l, size=4, alpha=0.4, colour="red") + 
	 	geom_text(data=l, aes(label=cell), vjust=3, size=2) + 
	 	xlab("") + ylab("") + ggtitle(plot.title) +
	 	geom_segment(data=edge.df, aes(x=x, y=y, xend=x1, yend=y1), color="grey", size=.3, linetype="solid", na.rm=TRUE) + 
	 	theme(#axis.line = element_line(colour = "black"),
		    panel.grid.major = element_blank(),
		    axis.text.x=element_blank(),
          	axis.text.y=element_blank(),axis.ticks=element_blank(),
		    panel.grid.minor = element_blank(),
		    panel.border = element_blank(),
		    panel.background = element_blank(), 
		    strip.background = element_blank(), 
		    text = element_text(size=12, colour="gray22"))

 	if(show.cells)
 	{
 		g = g + geom_point(size=0.3, aes(col=dens))
 	}
 	g
}


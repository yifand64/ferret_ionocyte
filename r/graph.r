


### updated infomap graph-clustering. much more scalable than the old version
# 1. does not compute the full distance matrix so can run on large data
# 2. uses approximate knn from largeVis so should be faster (haven't tested)
infomap <- function(d, k, 	n.cores=1, 
							neighbors=NULL,
							dm = "Euclidean", 
							similarity = NULL,
							perplexity = 100000,
							trials=10, 
							plot_nng=F)
{
    library(largeVis)
    library(igraph)
    source("util.r")
    if(!dm %in% c("Euclidean", "Cosine"))
    {
    	stop("dm must be one of 'Euclidean' or 'Cosine'. largeVis will segfault if it isn't!")
    }
    # soon should parallelize this function over K.
 #    cl <- makeCluster(as.integer(n.cores))
	# clusterExport(cl, "igraph")
	# p = pblapply(1:length(k), function(i) {
	# 		k.use = k[]
	# 		im = 
	# 	}, cl=cl)
	# 	info("Done!")

    #if(is.null(tree_threshold)){ tree_threshold = ncol(d); info(sprintf("Setting treeThreshold to %s", tree_threshold))}
    if(is.null(neighbors))
    {
	    info(sprintf("Building approximate kNN graph using largeVis. [nsamples=%s, ncores=%s, k=%s, trials=%s, distance=%s]", nrow(d), n.cores, k, trials, dm))
	    neighbors <- randomProjectionTreeSearch(t(d), K = k, max_iter = 5, threads = n.cores, distance_method = dm)
    }else{
    	info(sprintf("Using provided kNN graph (k=%s)", nrow(neighbors)))
    }
    info("Converting to adjacency matrix")
    edges <- buildEdgeMatrix(data = t(d), neighbors = neighbors, distance_method = dm, verbose=T)
	# rm(neighbors)
	# gc()
	am <- buildWijMatrix(edges, perplexity=perplexity)
    g = igraph::graph_from_adjacency_matrix(am, mode="undirected", weighted=T)
    if(!is.null(similarity)){
    	info(sprintf("Computing %s graph", similarity))
    	g = similarity.graph(g, type=similarity)
    }
    
    info("Running infomap")
    r = igraph::infomap.community(g, modularity = F, nb.trials = trials)
    n_cl = length(unique(r$membership))
    size = mean(table(r$membership))
    info(sprintf("Done! Found %s clusters, average size: %s", n_cl, size))
    igraph::V(g)$name <- rownames(d)
    if(plot_nng){cairo_pdf("nng.pdf"); plot(g, vertex.color=r$membership); dev.off()}
    list("cluster"=r$membership, "igraph_obj"=r, "nng"=g)
}



# # graph based algorithm and visualisation stuff.
# library(igraph)
# library(matrixcalc)

# graphk <- function(counts, dm, k.range=2:50)
# {
# 	vals = list()
# 	i = 1
# 	for(k in k.range)
# 	{
# 		g = graph.cluster(counts, dm=dm, k=k)
# 		s = silhouette(g$partition, dist=dm)
# 		val = mean(s[,"sil_width"])
# 		vals[[i]] = val
# 		i = i + 1
# 		cat(sprintf("k=%s\n", k))
# 	}
	
# 	return (unlist(vals))
# }



# test.k.vs.n.clusters <- function(counts, k.range=2:10, method="euclidean", sim.type=NULL, plot=T, pdf.output=NULL, use.markov=TRUE)
# {
# 	louvain = list()
# 	louvain.Q = list()
# 	markov = list()
# 	markov.Q = list()
# 	infomap = list()
# 	infomap.Q = list()

# 	info(sprintf("Calcuting %s distance matrix.. ", method))
# 	if(method != "cosine")
# 	{
# 		dm = as.matrix(amap::Dist(counts, method=method, nbproc=detectCores()))
# 	}else
# 	{	
# 		dm = cos.dist(counts)
# 	}	
# 	count=1
# 	for(k in k.range)
# 	{
# 		info(sprintf("K-setting %s of %s..", count, length(k.range)))
# 		knn = knn.graph(counts, k=k, distM=dm)
# 		if(!is.null(sim.type))
# 		{
# 			to.cluster = similarity.graph(knn, type=sim.type)
# 		}else{
# 			to.cluster = knn
# 		}

# 		if(use.markov)
# 		{
# 			m = cluster.markov(to.cluster)
# 			clusters.m = m$Cluster
# 			m.q = igraph::modularity(to.cluster, clusters.m)
# 		}
		
# 		clusters.l = multilevel.community(as.undirected(to.cluster))$membership
# 		clusters.i = infomap.community(to.cluster, modularity=TRUE)$membership

# 		#modularity scores help evaluate clustering performance.
# 		i.q = igraph::modularity(to.cluster, clusters.i)
# 		l.q = igraph::modularity(to.cluster, clusters.l)

		
# 		louvain.Q = unlist(c(louvain.Q, l.q))
# 		infomap.Q = unlist(c(infomap.Q, i.q))
# 		if(use.markov)
# 		{
# 			markov.Q = unlist(c(markov.Q, m.q))
# 			markov = unlist(c(markov, length(unique(clusters.m))))
# 		}
		
# 		louvain = unlist(c(louvain, length(unique(clusters.l))))
# 		infomap = unlist(c(infomap, length(unique(clusters.i))))
# 		cat("			---------------\n\n")
# 		#print(louvain.Q)
# 		count = count + 1
# 	}
# 	if(use.markov)
# 	{
# 		x = data.frame(cbind(k.range, louvain, louvain.Q, markov, markov.Q, infomap, infomap.Q))
# 	}else
# 	{
# 		x = data.frame(cbind(k.range, louvain, louvain.Q, infomap, infomap.Q))
# 	}
	
# 	# print(x)
# 	if(plot)
# 	{
# 		sep="_"
# 		extract.field=function(string,field=1,delim="_") {return(strsplit(string,delim)[[1]][field])}
# 		y = x
# 		info("Plotting clusters vs k..")
# 		y$louvain = paste(y$louvain, y$louvain.Q, sep=sep)
# 		if(use.markov)
# 		{
# 			y$markov = paste(y$markov, y$markov.Q, sep=sep)
# 			y$markov.Q = NULL
# 		}
		
# 		y$infomap = paste(y$infomap, y$infomap.Q, sep=sep)
# 		y$louvain.Q = NULL
# 		y$infomap.Q = NULL

# 		df = melt(data.frame(y), id.vars = "k.range")
# 		df["N.Cluster"] = as.numeric(unlist(lapply(as.character(df$value),extract.field,1)))
# 		df["Modularity"] = as.numeric(unlist(lapply(as.character(df$value),extract.field,2)))
# 		df$value = NULL

# 		# print(df)
		
# 		if(!is.null(pdf.output))
# 		{
# 			pdf(pdf.output)
# 		}
# 		log.breaks = seq(0, max(df$N.Cluster), by = 5)
# 		log.breaks[1] = 1 # watch out for inf otherwise
# 		g = ggplot(df, aes(x=k.range, y=N.Cluster, color=variable, size=Modularity))  + geom_line(size=2) + geom_point() + 
# 			xlab("K (nearest neighbors)") + scale_shape_identity() +
# 			ylab("N Clusters") + theme_bw() + scale_color_brewer("Clustering Alg", palette="Set1") + scale_y_log10(breaks = log.breaks) + 
# 			ggtitle("Clusters detected [Log scale] as a function of  k (nearest neighbors). \n Parition modularity is shown as sizes.")
# 		print(g)
# 		# tried geom_point(aes(fill = factor(variable), shape=21, size=Modularity)) to get black border, but didn't work.
# 		g = ggplot(df, aes(x=k.range, y=N.Cluster, color=variable, size=Modularity))  + geom_line(size=2) + geom_point() + xlab("K (nearest neighbors)") + 
# 			ylab("N Clusters") + theme_bw() + scale_color_brewer("Clustering Alg", palette="Set1") + scale_shape_identity() + 
# 			scale_y_continuous(breaks = seq(0, max(df$N.Cluster), by = 5))+ 
# 			ggtitle("Clusters detected as a function of  k (nearest neighbors). \n Parition modularity is shown as sizes.")
# 		print(g)

# 		if(!is.null(pdf.output))
# 		{
# 			dev.off()
# 		}
# 	}
# 	return (x)
# }

# # run a graph clustering pipeline on a data matrix. kNN -> similarity -> community detection. 
# # this is basically phenograph.
# test.graph.basic <- function(counts, k=0, knn.type="euclidean", sim.type="invlogweighted", plot=F)
# {
# 	knn = knn.graph(counts, k=k, dist.type=knn.type)
# 	sim = similarity.graph(knn, type=sim.type)
# 	#sim = delete.edges(sim, which(E(sim)$weight <=.1))
	
	
# 	clusters.m = cluster.markov(sim)$Cluster
# 	clusters.l = multilevel.community(as.undirected(sim))$membership
# 	clusters.i = infomap.community(sim)$membership

# 	info(sprintf("Louvain clustering found %s clusters..", length(unique(clusters.l))))
# 	info(sprintf("Infomap clustering found %s clusters..", length(unique(clusters.i))))
# 	pdf(paste(sim.type,".edge.weight.hist.pdf",sep="")); hist(E(sim)$weight); dev.off()
	
# 	if(plot){
# 		before = ecount(sim)
# 		sim.for.plot = delete.edges(sim, which(E(sim)$weight <=.2))
# 		after = ecount(sim.for.plot)
# 		info("Generating graph layouts..")
# 		k = layout.fruchterman.reingold.grid(knn)
# 		l = layout.fruchterman.reingold.grid(sim.for.plot)
# 		info(sprintf("Removed %s edges for plotting [%s remain]..", before-after, after))
# 		pdf("test.knn.markov.pdf", width=12, height=12); plot.graph(knn, layout=k, node.clustering = clusters.m, node.size=250, node.label.size=0.4); dev.off()
# 		pdf("test.knn.louvain.pdf", width=12, height=12); plot.graph(knn, layout=k, node.clustering = clusters.l, node.size=250, node.label.size=0.4); dev.off()
		
# 		pdf("test.sim.infomap.pdf", width=12, height=12); plot.graph(sim.for.plot, layout=l, node.clustering = clusters.i, node.size=500, node.label.size=0.2); dev.off()
# 		pdf("test.sim.markov.pdf", width=12, height=12); plot.graph(sim.for.plot, layout=l, node.clustering = clusters.m, node.size=500, node.label.size=0.2); dev.off()
# 		pdf("test.sim.louvain.pdf", width=12, height=12); plot.graph(sim.for.plot, layout=l, node.clustering = clusters.l, node.size=500, node.label.size=0.2); dev.off()
		
# 		pdf("test.community.sim.louvain.pdf", width=12, height=12); plot.graph(sim.for.plot, node.clustering = clusters.l, community=T, node.size=5); dev.off()
# 	}
	

# }

# plot.graph <- function(g, 
# 	node.names=NULL, 
# 	node.clustering=NULL, 
# 	node.label.size=2, 
# 	use.cols=NULL,
# 	layout=NULL,
# 	layout.type=NULL,
# 	grid=F, 
# 	simplify=F, 
# 	community=F,
# 	node.size=500,
# 	max.edges.for.plot=100000,
# 	edge.width=1)
# {
# 	if(is.null(node.names))
# 	{
# 		#node.names = c(rep("", vcount(g)))
# 		node.names = V(g)$name
# 	}
# 	V(g)$label.cex = node.label.size
	
# 	if(!is.null(node.clustering))
# 	{
# 		if(length(node.clustering) != vcount(g))
# 		{
# 			error("Length of cluster vector must match number of nodes in graph!")
# 			return (FALSE)
# 		}
# 		n.colors = length(unique(node.clustering))
# 		if(is.null(use.cols))
# 		{
# 			use.cols = intense.cols(n.colors)	
# 		}
# 		if(community) # simplified graph showing community structure only
# 		{
# 			g <- contract.vertices(g, node.clustering)
# 			g <- igraph::simplify(g, remove.loops=TRUE)
# 			V(g)$color = use.cols
# 		}else
# 		{
# 			V(g)$color = mapvalues(node.clustering, unique(node.clustering), use.cols)
# 		}
		
# 		#E(g)$weight <- 1
		
# 	}else
# 	{
# 		if(community)
# 		{
# 			error("Cannot draw a community graph without being passed a graph clustering!")
# 			return (FALSE)
# 		}
# 	}

# 	if(simplify)
# 	{
# 		info("Simplifying graph..")
# 		g = igraph::simplify(g)
# 	}

# 	if(ecount(g) > max.edges.for.plot)
# 	{
# 		#TODO: should trim graph to the maximum edge num.
# 		warn(sprintf("Graph contains more than %s edges. Trimming some..",max.edges.for.plot))
# 		before = ecount(g)
# 		g = delete.edges(g, which(E(g)$weight <=.35))
# 		after = ecount(g)
# 		info(sprintf("Removed %s edges for plotting [%s remain]..", before-after, after))
# 	}

# 	if(is.null(layout))
# 	{
# 		info(sprintf("Laying out graph [%s nodes]..", vcount(g)))
# 		if(is.null(layout.type))
# 		{
# 			info("Using autolayout. Try: fruchterman.reingold, reingoild.tilford, kamada, or spring for better results")
# 			l = layout.auto(g)
# 		}else{
# 			if(layout.type=="fruchterman.reingold"){l = layout.fruchterman.reingold(g)}
# 			else{if(layout.type=="reingold.tilford"){
# 					g = simplify(g)
# 					l <- layout.reingold.tilford(g, circular=T)}else{
# 					if(layout.type=="kamada"){l = layout.kamada.kawai(g)}else{
# 						if(layout.type == "spring"){l = layout.spring(g)}else{
# 							error("unknown layout type")
# 						}
# 					}
# 				}
# 			}
# 		}
# 	}else{
# 		l = layout
# 	}
	
# 	info("Plotting..")
# 	plot(g, layout=l, vertex.size=node.size, edge.width=edge.width, edge.arrow.size=0,
# 			vertex.label=node.names, rescale=FALSE, xlim=range(l[,1]), ylim=range(l[,2]),
# 		 	vertex.label.dist=1)
# 	return(l)
# }


# pc.dist <- function(data, pcs.use=1:10, distance.method="euclidean", scale=F, center=T)
# {
# 	library(amap)
# 	library(parallel)
# 	info("Running PCA..")
# 	pca = prcomp(data, scale=scale, center=center)
# 	pca.d = pca$x[, pcs.use]
# 	info(sprintf("Calculating %s distance matrix..", distance.method))
# 	dm = as.matrix(amap::Dist(pca.d, method=distance.method, nbproc=detectCores()))
# 	return (dm)
# }

# graph.type can be jaccard, invlogweighted or dice, community detect
# can be louvain, infomap or markov. 
graph.cluster <- function(	data, 
							graph.type="knn", # can be threshold (binarise the distance matrix), jaccard or knn.
							graph.threshold=0.5,
							dm=NULL,
							community.detect="infomap", 
							distance.method="euclidean",
							k=0)
{
	
	
	library(amap)
	if(k==0)
	{
		k=floor(log2(nrow(data)))
	}
	if(is.null(dm))
	{
		info(sprintf("No distance matrix provided. Calculating %s..", distance.method))
		if(distance.method!="cosine")
		{
			dm = as.matrix(amap::Dist(data, method=distance.method, nbproc=detectCores()))
		}else
		{
			dm = as.matrix(cos.dist(data))
		}
	}
	if(graph.type=="threshold"){	
		print(range(dm))	
		am =  (as.matrix(dm) < graph.threshold)  + 0 #adding 0 converts logical to numeric
		#to.cluster = igraph::simplify(graph.adjacency(am))
		to.cluster = graph.adjacency(am)
		#plot.graph(to.cluster)
	}else{
		knn = knn.graph(data, k=k, distM=dm)
		if(!identical(toupper("knn"), toupper(graph.type)))
		{
			to.cluster = similarity.graph(knn, type=graph.type)
		}else{
			to.cluster = knn
		}
	}
	

	if(identical(toupper(community.detect), toupper("markov")))
	{
		r = cluster.markov(to.cluster)
		clusters = r$Cluster
		
	}else{
		if(identical(toupper(community.detect), toupper("louvain")))
		{
			r = multilevel.community(as.undirected(to.cluster))
			clusters = r$membership
		}else{
			if(identical(toupper(community.detect), toupper("infomap")))
			{
				r = infomap.community(to.cluster, modularity=TRUE)
				clusters = r$membership
			}else{
				error(sprintf("Unknown community detection method: %s", community.detect))
				return (FALSE)
			}
		}
	}
	n.clusters =length(unique(clusters))
	
	f = function(i){as.vector(clusters==i)}
	clist= lapply(1:n.clusters, f)
	m = igraph::modularity(to.cluster, clusters)
	return (list("result"=r,
		"clustermethod"=paste(graph.type, "-graph clustering [", community.detect,"]", sep=""), 
		"nc"=n.clusters, 
		"modularity"=m, 
		"clusterlist"=clist,		
		"partition"=clusters))
}

# cluster.graph <- function(data, clusters, 
# 								genes=NULL,
# 								edge.threshold=0.7, 
# 								cluster.sim=NULL, 
# 								use.pca = T,
# 								pcs.use=1:25)
# {
# 	if(!is.null(genes)){data=data[genes,]}
# 	if(is.null(cluster.sim))
# 	{
# 		if(use.pca)
# 		{
# 			cluster.sim = cluster.similarity.matrix(data, clusters, dm=pc.dist(t(data), pcs.use=pcs.use))# requires ~/dev/adam/rna_seq/r/marker.r
# 		}else
# 		{
# 			cluster.sim = cluster.similarity.matrix(data, clusters)# requires ~/dev/adam/rna_seq/r/marker.r
# 		}
# 	}
# 	print(cluster.sim)
# 	if(use.pca)
# 	{	
# 		#if usa.pca is true, the distances are distances, so closely related clusters have small values in the matrix.
# 		am =  (cluster.sim < edge.threshold)  + 0 #adding 0 converts logical to numeric
# 	}else
# 	{
# 		# otherwise they are correlations, close values are large
# 		am =  (cluster.sim > edge.threshold)  + 0 #adding 0 converts logical to numeric
# 	}
	
# 	g =  simplify(graph.adjacency(am))
# 	n = length(unique(clusters))
# 	cols = colorRampPalette(brewer.pal(n, "Set1"))(n)
# 	plot.graph(g, 
# 		node.clustering = 1:n,
# 		use.cols = cols,
# 		grid=T, 
# 		node.size= 3 * sqrt(as.matrix(table(clusters))), 
# 		node.label.size = 1)
# 	return(g)
# }

# following phenograph, smooth the knn graph using the
# jaccard coefficient. basically construct a weighted
# graph where the weight of an edge between cells i,j
# is the number of neighbor cells they have in common.
# types: jaccard, dice, invlogweighted. 
# see: http://www.inside-r.org/packages/cran/igraph/docs/similarity
similarity.graph <- function(knn.g, type="jaccard", simplify=T)
{
	types = c("invlogweighted", "dice", "jaccard")
	if(!type %in% types){stop(sprintf("Type must be one of %s", paste(types, collapse=", ")))}	
	info(sprintf("Building %s-similarity graph..", type))
	if(type=="jaccard")
	{
		s = igraph::similarity.jaccard(knn.g)
	}else
	{
		if(type=="dice")
		{
			s = igraph::similarity.dice(knn.g)
		}else
		{
			if(type=="invlogweighted")
			{
				s = igraph::similarity.dice(knn.g)
			}
		}
	}
	g = graph.adjacency(s, weighted=TRUE)
	#V(g)$name = V(knn.g)$name
	info(sprintf("%s-similarity graph computed. Average degree: %s", type, mean(degree(g))))
	if(simplify)
	{
		info("Simplifying graph..")
		g = igraph::simplify(g)
	}
	info(sprintf("Average degree after self- and multi-edges removed: %s", mean(degree(g))))
	return (g)
}

knn.graph <- function(data, k=0, dist.type="euclidean", distM=NULL, verbose=F)
{
	library(cccd) 
	if(k==0)
	{
		k = floor(sqrt(nrow(data))/2)
	}
	if(verbose)
	{
		info(sprintf("Building %s-nearest [%s] neighbor graph..", k, dist.type))
	}
	
	if(is.null(distM))
	{
		if(dist.type == "pearson" | dist.type == "spearman" | dist.type=="euclidean")
		{
			n.cores = detectCores()
			if(verbose)
			{
				info(sprintf("Calculating %s distances.. [Using %s cores]", dist.type, n.cores))
			}
			dmat = as.matrix(amap::Dist(data, method=dist.type, nbproc=n.cores))
			#dmat = 1-cor(as.matrix(data), method=dist.type)
			g <- nng(dx=dmat,k=k)
		}else{
			if(dist.type == "cosine")
			{
				dmat = as.matrix(cos.dist(data))
				g <- nng(dx=dmat,k=k)
			}else
			{
				g <- nng(x=as.matrix(data),k=k, method=dist.type)
			}
		}
	}else
	{
		g <- nng(dx=distM,k=k)
	}
	
	V(g)$name = rownames(data)
	if(verbose)
	{
		info(sprintf("%s %s-NN computed. Average degree: %s", dist.type, k, mean(degree(g))))
	}
	return(g)
}

cluster.markov <- function(g)
{
	library(MCL)
	info("Building adjacency matrix..")
	num.nodes = vcount(g)
	a = get.adjacency(g, sparse = F)
	info("Running Markov clustering..")

	m = tryCatch({
	    m = mcl(a, addLoops=TRUE)
	    info(sprintf("Found %s clusters: ", length(table(m$Cluster))))
	    return(m)
	}, warning = function(w) {
	    warn(w)
	    return (NULL)
	}, error = function(e) {
		error("Markov clustering failed!")
	    error(e)
	    m = list("Cluster" = rep(0, num.nodes))
	    info(sprintf("Returning all nodes in %s cluster ..", length(table(m$Cluster))))
	    return(m)
	}, finally = {
	    
	})
	
	return(m)
}



# cos.dist <- function(X, corr=FALSE){
# 	X = as.matrix(X)
# 	if(corr){ X = apply(X, 2, function(x){ x-mean(x) }) }
# 	tosolve = diag(sqrt(diag(t(X)%*%X)))
# 	if(is.singular.matrix(tosolve))
# 	{
# 		return (0)
# 	}else
# 	{
# 		denom = solve(tosolve)	
# 	}
# 	rval = denom%*%(t(X)%*%X)%*%denom 
# 	rownames(rval) = colnames(X)
# 	colnames(rval) = colnames(X)
# 	return( rval)
# } 

# methods for 'topological data analysis' or TDA.

require(MASS)
require(amap)
require(hexbin)
require(RColorBrewer)

# A Reeb graph (named after Georges Reeb) is a mathematical object reflecting the evolution of the level sets of a real-valued function on a manifold.[1] 
# Using 'filter functions' which are real-valued, i.e the stemness score, l-infinity centrality, whatever, and then binning the data, we can generate
# a Reeb graph that tells something about the shape of a dataset. Inspired by:
# http://www.nature.com/articles/srep01236

reeb.graph <- function(filtered.data, resolution=10, percent.overlap=5, cluster=F)
{
	d = as.data.frame(filtered.data)
	dim = ncol(filtered.data)
	info(sprintf("Building Reeb graph from % dimensional filtered data..", dim))
	info(sprintf("	Resolution=%s, Percent.Overlap=%s", resolution, percent.overlap))

	cols = rev(colorRampPalette(brewer.pal(11, "Spectral"))(25))

	if(dim==2)
	{
		g = ggplot(d, aes_string(x=colnames(d)[1], y=colnames(d)[2])) + stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) + 
			scale_fill_gradientn(colours=cols)
		print(g)
		#bin <- hexbin(d[,1], d[,2])


		nbins <- resolution
		x.bin <- seq(floor(min(x)), ceiling(max(x)), length=nbins)
		y.bin <- seq(floor(min(y)), ceiling(max(y)), length=nbins)
	}
	
}

# given data, n observations (columns) of m variables (rows)
# return a vector of length n, the distance to its most distant neighbor 
# this is l.infinity centrality as defined in:
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3566620/
l.infinity.centrality <- function(data, distM=NULL, distance="euclidean")
{
	if(is.null(distM))
	{
		info(sprintf("Computing %s distance matrix..", distance))
		distM = as.matrix(Dist(t(data), method=distance, nbproc=detectCores()))
	}
	print(dim(distM))
	return (unlist(apply(distM, 2, max)))
}